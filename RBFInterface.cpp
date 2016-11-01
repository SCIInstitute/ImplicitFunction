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
#include <algorithm>

//#define TETLIBRARY   // Required definition for use of tetgen library
#include <tetgen.h>

const double RBFInterface::EPSILON = 1.0e-3;
const double RBFInterface::SMALL_EPSILON = 1.0e-6;

RBFInterface::RBFInterface(std::vector<vec3> myData,
                           vec3 myOrigin, vec3 mySize, vec3 mySpacing,
                           double myOffset, std::vector<axis_t> myAxis,
                           const bool useConvexHull, const bool compute2DConvexHull,
                           const bool invertSeedOrder, Kernel kernel) :
  thresholdValue_(0),
  useConvexHull_(useConvexHull),
  compute2DConvexHull_(compute2DConvexHull),
  invertSeedOrder_(invertSeedOrder),
  kernel_(kernel)
{
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

  if ( this->useConvexHull_ && ! this->compute2DConvexHull_ )
  {
    Create3DSurface();
  }
  else
  {
    CreateSurface(myData, myOrigin, mySize, mySpacing, myOffset, myAxis);
  }
}

void RBFInterface::Create3DSurface()
{
  // TODO: debug print
  std::cerr << "Calling Tetgen..." << std::endl;
  tetgenio in, out;
  in.numberofpoints = this->points_x_.size();
  in.pointlist = new REAL[in.numberofpoints * 3];
  for ( size_t i = 0; i < in.numberofpoints; ++i )
  {
    in.pointlist[i*3] = this->points_x_[i];
    in.pointlist[i*3+1] = this->points_y_[i];
    in.pointlist[i*3+2] = this->points_z_[i];
  }

  // no faces...
  try
  {
    tetrahedralize("BE", &in, &out);
  }
  catch (std::exception& e)
  {
    std::cerr << "TetGen failed to generate a mesh: " << e.what() << std::endl;
    return;
  }
  catch (...)
  {
    std::cerr << "TetGen failed to generate a mesh" << std::endl;
    return;
  }

  // trifaces == convex hull
  std::cerr << "# tri faces=" << out.numberoftrifaces << std::endl;
  for (int i = 0; i < out.numberoftrifaces; ++i)
  {
    std::cerr << "trifacelist " << i << ": " << out.trifacelist[i*3] << " "
                                             << out.trifacelist[i*3+1] << " "
                                             << out.trifacelist[i*3+2] << std::endl;
  }

  std::cerr << "# points=" << out.numberofpoints << std::endl;
  for (int i = 0; i < out.numberofpoints; ++i)
  {
    std::cerr << "pointlist " << i << ": " << out.pointlist[i*3] << " "
                                           << out.pointlist[i*3+1] << " "
                                           << out.pointlist[i*3+2] << std::endl;
  }

  std::vector< std::vector< int > > listOfIntsPerVertex(out.numberofpoints);
  for (int i = 0; i < out.numberoftrifaces; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      int index = out.trifacelist[i*3+j];
      listOfIntsPerVertex[index].push_back(i);
    }
  }

  for (int i = 0; i < listOfIntsPerVertex.size(); ++i)
  {
    std::cerr << i << ": ";
    for (int j = 0; j < listOfIntsPerVertex[i].size(); ++j)
    {
      std::cerr << listOfIntsPerVertex[i][j] << " ";
    }
    std::cerr << std::endl;
  }

  // normal calculation on face
  // For triangle p1, p2, p3 and vectors U = p2 - p1, V = p3 - p1, then normal N = UxV:
  //  Nx = UyVz - UzVy
  //  Ny = UzVx - UxVz
  //  Nz = UxVy - UyVx

  std::vector< vec3 > normalsPerFace;
  for (int i = 0; i < out.numberoftrifaces; ++i)
  {
    const int i1 = out.trifacelist[i*3];
    const int i2 = out.trifacelist[i*3+1];
    const int i3 = out.trifacelist[i*3+2];

    vec3 p1( out.pointlist[i1], out.pointlist[i1+1], out.pointlist[i1+2] );
    vec3 p2( out.pointlist[i2], out.pointlist[i2+1], out.pointlist[i2+2] );
    vec3 p3( out.pointlist[i3], out.pointlist[i3+1], out.pointlist[i3+2] );

    vec3 u = p2 - p1;
    vec3 v = p3 - p1;

    vec3 n(
           ( u.y() * v.z() ) - ( u.z() * v.y() ),
           ( u.z() * v.x() ) - ( u.x() * v.z() ),
           ( u.x() * v.y() ) - ( u.y() * v.x() )
          );
    normalsPerFace.push_back(n);
  }

  std::vector< vec3 > normalsPerVertex(out.numberofpoints);
  for (int i = 0; i < listOfIntsPerVertex.size(); ++i)
  {
    const int tmpSize = listOfIntsPerVertex[i].size();
    vec3 tmpVec;
    std::cerr << "i=" << i << " tmpSize=" << tmpSize << " ";
    for (int j = 0; j < tmpSize; ++j)
    {
      tmpVec += normalsPerFace[ listOfIntsPerVertex[i][j] ];
      std::cerr << "j=" << j << " normalsPerFace=" << normalsPerFace[ listOfIntsPerVertex[i][j] ] << " ";
    }
    std::cerr << std::endl;
    normalsPerVertex[i] = normalize(tmpVec, 1.0e-6);
  }

  for (int i = 0; i < normalsPerVertex.size(); ++i)
  {
    std::cerr << normalsPerVertex[i] << std::endl;
  }

}


// TODO: why is myData arg???
// driver
void RBFInterface::CreateSurface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing, double myOffset, std::vector<axis_t> myAxis)
{
  // This code figures out the bounds for the data selected
  std::vector<double>::iterator minx = std::min_element(this->points_x_.begin(), this->points_x_.end());
  std::vector<double>::iterator miny = std::min_element(this->points_y_.begin(), this->points_y_.end());
  std::vector<double>::iterator minz = std::min_element(this->points_z_.begin(), this->points_z_.end());

  std::vector<double>::iterator maxx = std::max_element(this->points_x_.begin(), this->points_x_.end());
  std::vector<double>::iterator maxy = std::max_element(this->points_y_.begin(), this->points_y_.end());
  std::vector<double>::iterator maxz = std::max_element(this->points_z_.begin(), this->points_z_.end());

  vec3 myMin(*minx, *miny, *minz), myMax(*maxx, *maxy, *maxz);

  // TODO: magic number
  myMin = myMin - 0.05 * mySize;
  myMax = myMax + 0.05 * mySize;

  this->surfaceData_ = new ScatteredData(this->points_x_, this->points_y_, this->points_z_, this->threshold_, myAxis);
  //this->surfaceData_->axisInformation_ = myAxis;
  if ( this->useConvexHull_ && this->compute2DConvexHull_ ) // TODO: not sure this is necessary...
  {
    this->surfaceData_->compute2DHull();
    // TODO: this is bad - ScatteredData should set this to maintain correct internal state!!!
    this->surfaceData_->origSize_ = this->surfaceData_->surfacePoints_[0].size();
  }
  else
  {
    this->surfaceData_->origSize_ = this->points_x_.size();
  }

  // TODO: does not include points outside convex hull (if computed)...
  augmentNormalData(myOffset);
  this->rbf_ = new RBF(this->surfaceData_, kernel_);
  this->rbf_->setDataReduction(All);

  // Construct RBFs
  //mySurface->computeRBF();
  this->rbf_->computeFunction();
  
  // sanity check
  for (int i = 0; i < this->surfaceData_->fnc_.size(); i++)
  {
    // X == 0, Y == 1, Z == 2
    // TODO: also how to convert to vec3 structure...
    vec3 location(this->surfaceData_->surfacePoints_[0][i],
                  this->surfaceData_->surfacePoints_[1][i],
                  this->surfaceData_->surfacePoints_[2][i]);
    //double myVal = mySurface->computeValue(location);

    double myVal = this->rbf_->computeValue(location);
//    printf("myVal: %lf, fnc: %lf\n", myVal, this->surfaceData_->fnc_[i]);
    double error  = fabs(myVal - this->surfaceData_->fnc_[i]);
    if (error > EPSILON)
    {
      printf("ERROR (numerics): %lf\n", error);
      fflush(stdout);
    }
  }
  
  // Fill the values into the vector.
  // In the first loop, we initialize the matrix with all values set to -100.
  // In the second loop, we change the values from -100 to the correct rasterData_ if the point in the domain described above.
  
  this->rasterData_.resize(static_cast<int>(mySize[0]));
  for (int i = 0; i < mySize[0]; i++)
  {
    this->rasterData_[i].resize(static_cast<int>(mySize[1]));
    for (int j = 0; j < mySize[1]; j++)
    {
      this->rasterData_[i][j].resize(static_cast<int>(mySize[2]), -100);
    }
  }
  
  for (int i = 0; i < mySize[0]; i++)
  {
    //vec3 location = myOrigin + mySpacing[0] * i * vec3::unitX;
    //if (location[0]<myMin[0]||location[0]>myMax[0])
    //  continue;
    for (int j = 0; j < mySize[1]; j++)
    {
      //location = myOrigin + mySpacing[1]*j*vec3::unitY;
      //if (location[1]<myMin[1]||location[1]>myMax[1])
      //  continue;
      for (int k = 0; k < mySize[2]; k++)
      {
        // TODO: shadowing variables in outer scopes!!!
        //location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
        vec3 location = myOrigin + mySpacing[0] * i * vec3::unitX + mySpacing[1] * j * vec3::unitY + mySpacing[2] * k * vec3::unitZ;

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
}

//vec3 RBFInterface::findNormal(ScatteredData *data, int n)
//{
//  int tot = data->origSize;
//  int prev = (n-1)>=0?n-1:tot-1;
//  int next = (n+1)<tot?n+1:0;
//  
//  while(data->surfacePoints_[2][prev]!=data->surfacePoints_[2][n])
//  {
//    prev = (prev-1)>=0?prev-1:tot-1;
//  }
//  
//  while(data->surfacePoints_[2][next]!=data->surfacePoints_[2][n])
//  {
//    next = (next+1)<tot?next+1:0;
//  }
//  printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);
//  
//  vec3 a(data->surfacePoints_[0][n], data->surfacePoints_[1][n], data->surfacePoints_[2][n]);
//  vec3 b(data->surfacePoints_[0][prev], data->surfacePoints_[1][prev], data->surfacePoints_[2][prev]);
//  vec3 c(data->surfacePoints_[0][next], data->surfacePoints_[1][next], data->surfacePoints_[2][next]);
//  
//  vec3 tangent = b-c;
//  //rotate by 90 degrees on the x-y plane
//  double ret_x = -tangent[1];
//  double ret_y = tangent[0];
//  vec3 ret(ret_x, ret_y, tangent[2]);
//  
//  return ret;
//}
//
//vec3 RBFInterface::findSphericalNormal(ScatteredData *data, int n)
//{
//  vec3 ret(0,0,0);
//  for(int j=0; j<3; j++)
//    // TODO: hardcoded normal scaling?
//    ret[j] = (data->surfacePoints_[j][n] - data->centroid[j])/10;
//  
//  return ret;
//}

// TODO: move this and findNormalAxis to new class?
void RBFInterface::augmentNormalData(const double myOffset)
{
  const size_t DIM_3D = 3;
  const int NORMAL_IN = 10, NORMAL_OUT = -10;
  const size_t N = this->surfaceData_->origSize_;
  const size_t M = this->surfaceData_->leftovers_[0].size();

  if ( this->useConvexHull_ )
  {
    // iterate through list of points not on hull, add to list as zero points
    for (int i = 0; i < M; ++i)
    {
      this->surfaceData_->surfacePoints_[0].push_back(this->surfaceData_->leftovers_[0][i]);
      this->surfaceData_->surfacePoints_[1].push_back(this->surfaceData_->leftovers_[1][i]);
      this->surfaceData_->surfacePoints_[2].push_back(this->surfaceData_->leftovers_[2][i]);
    }
    this->surfaceData_->fnc_.insert(this->surfaceData_->fnc_.end(), M, this->thresholdValue_);
  }

  // size of surfaceData_->surfacePoints_ entries tripled with this code
  for (int i = 0; i < N; i++)
  {
    vec3 myNormal = findNormalAxis(i);
    myNormal = normalize(myNormal, SMALL_EPSILON);
//std::cerr << "normal: [ " << myNormal << " ]" << std::endl;
    vec3 inNorm;
    for (int j = 0; j < DIM_3D; j++)
    {
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] + myOffset * myNormal[j]);
      inNorm[j] = this->surfaceData_->surfacePoints_[j][i] + myOffset * myNormal[j];
      inNormals.push_back(inNorm);
    }

    // normals point inward
    this->surfaceData_->fnc_.push_back(NORMAL_IN);

    vec3 outNorm;
    for (int j = 0; j < DIM_3D; j++)
    {
      this->surfaceData_->surfacePoints_[j].push_back(this->surfaceData_->surfacePoints_[j][i] - myOffset * myNormal[j]);
      outNorm[j] = this->surfaceData_->surfacePoints_[j][i] - myOffset * myNormal[j];
      outNormals.push_back(outNorm);
    }

    // normals pointing outward
    this->surfaceData_->fnc_.push_back(NORMAL_OUT);
  }

//  std::cerr << "points x component: ";
//  for (int i = 0; i < this->surfaceData_->surfacePoints_[0].size(); ++i)
//  {
//    std::cerr << this->surfaceData_->surfacePoints_[0][i] << " ";
//  }
//  std::cerr << std::endl;
//
//  std::cerr << "points y component: ";
//  for (int i = 0; i < this->surfaceData_->surfacePoints_[1].size(); ++i)
//  {
//    std::cerr << this->surfaceData_->surfacePoints_[1][i] << " ";
//  }
//  std::cerr << std::endl;
//
//  std::cerr << "points z component: ";
//  for (int i = 0; i < this->surfaceData_->surfacePoints_[2].size(); ++i)
//  {
//    std::cerr << this->surfaceData_->surfacePoints_[2][i] << " ";
//  }
//  std::cerr << std::endl;
}

vec3 RBFInterface::findNormalAxis(const int n)
{
  //printf("here\n");
  const int TOT = this->surfaceData_->origSize_;
  int prev = (n-1) >= 0 ? n-1 : TOT-1; // wrap
  int next = (n+1) < TOT ? n+1 : 0; // wrap
  //axis_t myAxis = this->surfaceData_->axisInformation_[n];
  axis_t myAxis = this->surfaceData_->updatedAxisInformation_[n];

//printf("Computing normals (1): prev=%d, n=%d, next=%d, TOT=%d\n", prev,n,next,TOT); fflush(stdout);
//std::cerr << fabs(this->surfaceData_->surfacePoints_[myAxis][prev] - this->surfaceData_->surfacePoints_[myAxis][n]) << ", " << fabs(this->surfaceData_->surfacePoints_[myAxis][next] - this->surfaceData_->surfacePoints_[myAxis][n]) << std::endl;

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