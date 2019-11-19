//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include "RBFInterface.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>

#include <tetgen.h>

#include <ctime>

//const double RBFInterface::EPSILON = 1.0e-3;
const double RBFInterface::SMALL_EPSILON = 1.0e-6;

RBFInterface::RBFInterface(std::vector<vec3> myData,
                           const vec3& myOrigin, const vec3& mySize, const vec3& mySpacing,
                           const double myOffset, AxisList myAxis,
                           const bool compute2DConvexHull,
                           const bool invertSeedOrder, Kernel kernel) :
  thresholdValue_(0),
  origin_(myOrigin),
  size_(mySize),
  spacing_(mySpacing),
  offset_(myOffset),
  axisList_(myAxis),
  compute2DConvexHull_(compute2DConvexHull),
  invertSeedOrder_(invertSeedOrder),
  kernel_(kernel)
{
  time_t tstart, tend;
  tstart = time(0);
  if ( this->invertSeedOrder_ )
  {
    // inplace
    std::reverse( myData.begin(), myData.end() );
  }

  for (int i = 0; i < myData.size(); i++)
  {
//    std::cerr << "point: " << myData[i][0] << ", " << myData[i][1] << ", " << myData[i][2] << std::endl;
    // TODO: would be better to skip this step...
    // Input point components:
    this->points_x_.push_back(myData[i][0]); // X component
    this->points_y_.push_back(myData[i][1]); // Y component
    this->points_z_.push_back(myData[i][2]); // Z component
    this->threshold_.push_back(this->thresholdValue_);
  }

  // TODO: error reporting???
  // TODO: would it be better to break this out, or have constructor throw exception???
  if ( ! this->compute2DConvexHull_ )
  {
    create3DSurface();
  }
  else
  {
    create2DSurface();
  }
  tend = time(0);
  cout << "RBF Interface constructor took " << difftime(tend, tstart) << " second(s)." << endl;
}

//RBFInterface::~RBFInterface()
//{
//std::cerr << "RBFInterface::~RBFInterface()" << std::endl;
//  delete this->surfaceData_;
//}

void RBFInterface::create3DSurface()
{
  // TODO: debug print
  std::cerr << "Calling Tetgen..." << std::endl;
  tetgenio in, out;
  in.numberofpoints = static_cast<int>( this->points_x_.size() );
  in.pointlist = new REAL[in.numberofpoints * 3];
  for ( size_t i = 0; i < in.numberofpoints; ++i )
  {
    in.pointlist[i*3] = this->points_x_[i];
    in.pointlist[i*3+1] = this->points_y_[i];
    in.pointlist[i*3+2] = this->points_z_[i];
  }

  try
  {
    std::string args("BE");
    tetrahedralize(const_cast< char* >( args.c_str() ), &in, &out);
  }
  catch (std::exception& e)
  {
    std::cerr << "TetGen failed to generate a mesh: " << e.what() << std::endl;
    throw;
  }
  catch (...)
  {
    std::cerr << "TetGen failed to generate a mesh" << std::endl;
    throw;
  }

  const size_t NUMBER_POINTS = out.numberofpoints;
  const size_t NUMBER_TRI_FACES = out.numberoftrifaces;
  const size_t NUMBER_TRI_POINTS = 3;

#ifndef NDEBUG
  // trifaces == convex hull
  std::cerr << "# tri faces=" << NUMBER_TRI_FACES << std::endl;
  for (size_t i = 0; i < out.numberoftrifaces; ++i)
  {
    std::cerr << "trifacelist " << i << ": " << out.trifacelist[i*NUMBER_TRI_POINTS] << " "
                                             << out.trifacelist[i*NUMBER_TRI_POINTS+1] << " "
                                             << out.trifacelist[i*NUMBER_TRI_POINTS+2] << std::endl;
  }

  std::cerr << "# points=" << NUMBER_POINTS << std::endl;
  for (size_t i = 0; i < NUMBER_POINTS; ++i)
  {
    std::cerr << "pointlist " << i << ": " << out.pointlist[i*NUMBER_TRI_POINTS] << " "
                                           << out.pointlist[i*NUMBER_TRI_POINTS+1] << " "
                                           << out.pointlist[i*NUMBER_TRI_POINTS+2] << std::endl;
  }
#endif

  IndexList* listOfIntsPerVertex = new IndexList[NUMBER_POINTS];
  for (size_t i = 0; i < NUMBER_TRI_FACES; ++i)
  {
    for (size_t j = 0; j < NUMBER_TRI_POINTS; ++j)
    {
      size_t index = out.trifacelist[i*NUMBER_TRI_POINTS+j];
      listOfIntsPerVertex[index].push_back(i);
    }
  }

#ifndef NDEBUG
  for (int i = 0; i < NUMBER_POINTS; ++i)
  {
    std::cerr << i << ": ";
    for (int j = 0; j < listOfIntsPerVertex[i].size(); ++j)
    {
      std::cerr << listOfIntsPerVertex[i][j] << " ";
    }
    std::cerr << std::endl;
  }
#endif

  // normal calculation on face
  // For triangle p1, p2, p3 and vectors U = p2 - p1, V = p3 - p1, then normal N = UxV:
  //  Nx = UyVz - UzVy
  //  Ny = UzVx - UxVz
  //  Nz = UxVy - UyVx

  vec3* normalsPerFace = new vec3[NUMBER_TRI_FACES];
  for (size_t i = 0; i < NUMBER_TRI_FACES; ++i)
  {
    //const size_t i1 = out.trifacelist[i*3+2];
    const size_t i1 = out.trifacelist[i*3];
    const size_t i2 = out.trifacelist[i*3+1];
    //const size_t i3 = out.trifacelist[i*3];
    const size_t i3 = out.trifacelist[i*3+2];

    vec3 p1( out.pointlist[i1*NUMBER_TRI_POINTS], out.pointlist[i1*NUMBER_TRI_POINTS+1], out.pointlist[i1*NUMBER_TRI_POINTS+2] );
    vec3 p2( out.pointlist[i2*NUMBER_TRI_POINTS], out.pointlist[i2*NUMBER_TRI_POINTS+1], out.pointlist[i2*NUMBER_TRI_POINTS+2] );
    vec3 p3( out.pointlist[i3*NUMBER_TRI_POINTS], out.pointlist[i3*NUMBER_TRI_POINTS+1], out.pointlist[i3*NUMBER_TRI_POINTS+2] );

    vec3 u = p2 - p1;
    vec3 v = p3 - p1;

    vec3 n(
           ( u.y() * v.z() ) - ( u.z() * v.y() ),
           ( u.z() * v.x() ) - ( u.x() * v.z() ),
           ( u.x() * v.y() ) - ( u.y() * v.x() )
          );
    normalsPerFace[i] = normalize(n, SMALL_EPSILON);
std::cerr << "normalsPerFace[" << i << "]=" << normalsPerFace[i] << ", len=" << length(normalsPerFace[i]) << std::endl;
  }

  // TODO: initialize in constructor?
  this->surfaceData_ = new ScatteredData(this->points_x_, this->points_y_, this->points_z_, this->threshold_, this->axisList_);

  // TODO: this code mirrors the 2D convex hull method...refactor?
  // TODO: difficult to keep the ScatteredData point and function vectors in sync...

  this->surfaceData_->surfacePoints_[0].clear();
  this->surfaceData_->surfacePoints_[1].clear();
  this->surfaceData_->surfacePoints_[2].clear();

  this->surfaceData_->leftovers_[0].clear();
  this->surfaceData_->leftovers_[1].clear();
  this->surfaceData_->leftovers_[2].clear();

  this->surfaceData_->fnc_.clear();

  // convex hull...

  for (size_t i = 0; i < NUMBER_POINTS; ++i)
  {
    // indices of points not in convex hull
    if ( listOfIntsPerVertex[i].size() == 0 )
    {
std::cerr << "Point " << i << " in leftovers." << std::endl;
      this->surfaceData_->leftovers_[0].push_back( out.pointlist[i*NUMBER_TRI_POINTS] );
      this->surfaceData_->leftovers_[1].push_back( out.pointlist[i*NUMBER_TRI_POINTS+1] );
      this->surfaceData_->leftovers_[2].push_back( out.pointlist[i*NUMBER_TRI_POINTS+2] );
    }
  }

  const size_t M = this->surfaceData_->leftovers_[0].size();

  std::vector< vec3 > normalsPerVertex;
  for (size_t i = 0; i < NUMBER_POINTS; ++i)
  {
    if ( listOfIntsPerVertex[i].size() == 0 ) continue;

std::cerr << "Point " << i << " in surface points." << std::endl;
    this->surfaceData_->surfacePoints_[0].push_back( out.pointlist[i*NUMBER_TRI_POINTS] );
    this->surfaceData_->surfacePoints_[1].push_back( out.pointlist[i*NUMBER_TRI_POINTS+1] );
    this->surfaceData_->surfacePoints_[2].push_back( out.pointlist[i*NUMBER_TRI_POINTS+2] );
    this->surfaceData_->fnc_.push_back( this->thresholdValue_ );

    vec3 tmpVec;
    for (size_t j = 0; j < listOfIntsPerVertex[i].size(); ++j)
    {
      tmpVec += normalsPerFace[ listOfIntsPerVertex[i][j] ];
    }

//    normalsPerVertex.push_back( normalize(tmpVec/3, SMALL_EPSILON) );
    normalsPerVertex.push_back( normalize(tmpVec, SMALL_EPSILON) );
  }

  this->surfaceData_->origSize_ = this->surfaceData_->surfacePoints_[0].size();
std::cerr << "#points=" << this->surfaceData_->surfacePoints_[0].size() << ", "
          << "#leftovers=" << this->surfaceData_->leftovers_[0].size() << std::endl;

for (int i = 0; i < normalsPerVertex.size(); ++i)
{
  std::cerr << "normalsPerVertex[" << i << "]=" << normalsPerVertex[i] << ", len=" << length(normalsPerVertex[i]) << std::endl;
}

  // iterate through list of points not on hull, add to list as zero points
  for (size_t i = 0; i < M; ++i)
  {
    this->surfaceData_->surfacePoints_[0].push_back(this->surfaceData_->leftovers_[0][i]);
    this->surfaceData_->surfacePoints_[1].push_back(this->surfaceData_->leftovers_[1][i]);
    this->surfaceData_->surfacePoints_[2].push_back(this->surfaceData_->leftovers_[2][i]);
  }
  this->surfaceData_->fnc_.insert(this->surfaceData_->fnc_.end(), M, this->thresholdValue_);

  const size_t N = this->surfaceData_->origSize_;
  const size_t DIM_3D = 3;
  const int NORMAL_IN = 10, NORMAL_OUT = -10;

  for (size_t i = 0; i < N; i++)
  {
//    vec3 point;
    for (size_t j = 0; j < DIM_3D; j++)
    {
      // TODO: check endpoint of normal from this->surfaceData_->surfacePoints_[j][i] + this->offset_ * normalsPerVertex[i][j] for inside c hull
      // generated printed warning
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] + this->offset_ * normalsPerVertex[i][j]);

//std::cerr << j << ": " << this->surfaceData_->surfacePoints_[j][i] << ", " << this->surfaceData_->surfacePoints_[j][i] + this->offset_ * normalsPerVertex[i][j] << std::endl;
//      point[j] = this->surfaceData_->surfacePoints_[j][i] + this->offset_ * normalsPerVertex[i][j];
//      this->surfaceData_->surfacePoints_[j].push_back( point[j] );
    }

    // check endpoint of normal from this->surfaceData_->surfacePoints_[j][i] + this->offset_ * normalsPerVertex[i][j] for inside c hull
    // printed warning, member variable for testing...
    //pointInsideConvexHull( point );

    // normals point inward
    this->surfaceData_->fnc_.push_back(NORMAL_IN);

    for (int j = 0; j < DIM_3D; j++)
    {
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] - this->offset_ * normalsPerVertex[i][j]);
    }

    // normals pointing outward
    this->surfaceData_->fnc_.push_back(NORMAL_OUT);
  }

  createRasterizedSurface();

  delete [] listOfIntsPerVertex;
  delete [] normalsPerFace;
}

//bool RBFInterface::pointInsideConvexHull( const vec3& point )
//{
////  std::cerr << point << std::endl;
//  return false;
//}

// driver
void RBFInterface::create2DSurface()
{
  // TODO: initialize in constructor?
  this->surfaceData_ = new ScatteredData(this->points_x_, this->points_y_, this->points_z_, this->threshold_, this->axisList_);

  //this->surfaceData_->axisInformation_ = myAxis;
  this->surfaceData_->compute2DHull();
  // TODO: this is bad - ScatteredData should set this to maintain correct internal state!!!
  this->surfaceData_->origSize_ = this->surfaceData_->surfacePoints_[0].size();

  // TODO: does not include points outside convex hull (if computed)...
  augmentNormalData();

  createRasterizedSurface();
}

void RBFInterface::createRasterizedSurface()
{
  // TODO: make local?
  this->rbf_ = new RBF(this->surfaceData_, kernel_);
  this->rbf_->setDataReduction(All);

  // Construct RBFs
  this->rbf_->computeFunction();  // TODO: throws exception if internal code used...

  // Fill the values into the vector.
  // In the first loop, we initialize the matrix with all values set to -100.
  // In the second loop, we change the values from -100 to the correct rasterData_ if the point in the domain described above.

  this->rasterData_.resize(static_cast<int>(this->size_[0]));
  for (int i = 0; i < this->size_[0]; i++)
  {
    this->rasterData_[i].resize(static_cast<int>(this->size_[1]));
    for (int j = 0; j < this->size_[1]; j++)
    {
      this->rasterData_[i][j].resize(static_cast<int>(this->size_[2]), -100);
    }
  }

  for (int i = 0; i < this->size_[0]; i++)
  {
    //vec3 location = this->origin_ + this->spacing_[0] * i * vec3::unitX;
    //if (location[0]<myMin[0]||location[0]>myMax[0])
    //  continue;
    for (int j = 0; j < this->size_[1]; j++)
    {
      //location = this->origin_ + this->spacing_[1]*j*vec3::unitY;
      //if (location[1]<myMin[1]||location[1]>myMax[1])
      //  continue;
      for (int k = 0; k < this->size_[2]; k++)
      {
        // TODO: shadowing variables in outer scopes!!!
        //location = this->origin_ + this->spacing_[0]*i*vec3::unitX + this->spacing_[1]*j*vec3::unitY + this->spacing_[2]*k*vec3::unitZ;
        vec3 location = this->origin_ + this->spacing_[0] * i * vec3::unitX + this->spacing_[1] * j * vec3::unitY + this->spacing_[2] * k * vec3::unitZ;

        //if (location[2]<myMin[2]||location[2]>myMax[2])
        //  continue;
        //std::cout<<"Computing Val ... "<<std::endl;
        //double myVal = mySurface->computeValue(location);

        double myVal = this->rbf_->computeValue(location);

        //printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
        this->rasterData_[i][j][k] = myVal;
      }
    }
  }

  //delete this->rbf_;
}

// TODO: move this and findNormalAxis to new class?
void RBFInterface::augmentNormalData()
{
  const size_t DIM_3D = 3;
  const int NORMAL_IN = 10, NORMAL_OUT = -10;
  const size_t N = this->surfaceData_->origSize_;
  const size_t M = this->surfaceData_->leftovers_[0].size();

  // iterate through list of points not on hull, add to list as zero points
  for (int i = 0; i < M; ++i)
  {
    this->surfaceData_->surfacePoints_[0].push_back(this->surfaceData_->leftovers_[0][i]);
    this->surfaceData_->surfacePoints_[1].push_back(this->surfaceData_->leftovers_[1][i]);
    this->surfaceData_->surfacePoints_[2].push_back(this->surfaceData_->leftovers_[2][i]);
  }
  this->surfaceData_->fnc_.insert(this->surfaceData_->fnc_.end(), M, this->thresholdValue_);

  // size of surfaceData_->surfacePoints_ entries tripled with this code
  for (int i = 0; i < N; i++)
  {
    vec3 myNormal = findNormalAxis(i);
    myNormal = normalize(myNormal, SMALL_EPSILON);
//std::cerr << "normal: [ " << myNormal << " ]" << std::endl;
    //vec3 inNorm;
    for (int j = 0; j < DIM_3D; j++)
    {
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] + this->offset_ * myNormal[j]);
      //inNorm[j] = this->surfaceData_->surfacePoints_[j][i] + this->offset_ * myNormal[j];
      //inNormals.push_back(inNorm);
    }

    // normals point inward
    this->surfaceData_->fnc_.push_back(NORMAL_IN);

    //vec3 outNorm;
    for (int j = 0; j < DIM_3D; j++)
    {
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] - this->offset_ * myNormal[j]);
      //outNorm[j] = this->surfaceData_->surfacePoints_[j][i] - this->offset_ * myNormal[j];
      //outNormals.push_back(outNorm);
    }

    // normals pointing outward
    this->surfaceData_->fnc_.push_back(NORMAL_OUT);
  }
}

vec3 RBFInterface::findNormalAxis(const int n)
{
  //printf("here\n");
  const int TOT = this->surfaceData_->origSize_;
  int prev = (n-1) >= 0 ? n-1 : TOT-1; // wrap
  int next = (n+1) < TOT ? n+1 : 0; // wrap
  //axis_t myAxis = this->surfaceData_->axisInformation_[n];
  axis_t myAxis = this->surfaceData_->updatedAxisInformation_[n];

  // TODO: why is this needed? AK 09/13/2016
  while (fabs(this->surfaceData_->surfacePoints_[myAxis][prev] - this->surfaceData_->surfacePoints_[myAxis][n]) > SMALL_EPSILON)
  {
    prev = (prev-1) >= 0 ? prev-1 : TOT-1; // wrap
  }
  
  while(fabs(this->surfaceData_->surfacePoints_[myAxis][next] - this->surfaceData_->surfacePoints_[myAxis][n]) > SMALL_EPSILON)
  {
    next = (next+1) < TOT ? next+1 : 0; // wrap
  }
//  printf("Computing normals (2): prev=%d, n=%d, next=%d, TOT=%d\n", prev,n,next,TOT); fflush(stdout);

  // normals from points either on convex hull boundary or in complete dataset if convex hull not used...

  // X == 0, Y == 1, Z == 2
  vec3 a(this->surfaceData_->surfacePoints_[0][n],
         this->surfaceData_->surfacePoints_[1][n],
         this->surfaceData_->surfacePoints_[2][n]);
  vec3 b(this->surfaceData_->surfacePoints_[0][prev],
         this->surfaceData_->surfacePoints_[1][prev],
         this->surfaceData_->surfacePoints_[2][prev]);
  vec3 c(this->surfaceData_->surfacePoints_[0][next],
         this->surfaceData_->surfacePoints_[1][next],
         this->surfaceData_->surfacePoints_[2][next]);

  //vec3 tangent = b-c;
  vec3 tan1 = normalize(b-a, SMALL_EPSILON);
  vec3 tan2 = normalize(a-c, SMALL_EPSILON);
  vec3 tangent = (tan1 + tan2) / 2;

  double ret_x, ret_y, ret_z;
  vec3 ret(tangent);
  //rotate by 90 degrees on the x-y plane
  switch(myAxis)
  {
    case 0:
      ret[1] = ret_y = -tangent[2];
      ret[2] = ret_z = tangent[1];
      break;
    case 1:
      ret[2] = ret_z = -tangent[0];
      ret[0] = ret_x = tangent[2];
      break;
    case 2:
      ret[0] = ret_x = -tangent[1];
      ret[1] = ret_y = tangent[0];
      break;
  }
  return ret;
}