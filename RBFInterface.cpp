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

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/RboxPoints.h>
#include <libqhullcpp/QhullFacetList.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

const double RBFInterface::EPSILON = 1.0e-3;
const double RBFInterface::SMALL_EPSILON = 1.0e-6;

RBFInterface::RBFInterface(std::vector<vec3> myData,
                           vec3 myOrigin, vec3 mySize, vec3 mySpacing,
                           double myOffset, std::vector<axis_t> myAxis,
                           const bool useConvexHull, const bool compute2DConvexHull,
                           Kernel kernel) :
  thresholdValue_(0),
  useConvexHull_(useConvexHull),
  compute2DConvexHull_(compute2DConvexHull),
  kernel_(kernel)
{
//for (int i = 0; i < myAxis.size(); ++i) {
//std::cerr << "axis " << i << ": " << myAxis[i] << std::endl;
//}
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
  const size_t POINT_LIST_SIZE = this->points_x_.size() * 3;
  double *pointList = new double[POINT_LIST_SIZE];
  for (size_t i = 0; i < this->points_x_.size(); ++i)
  {
    size_t j = i*3;
    pointList[j] = this->points_x_[i];
    pointList[j+1] = this->points_y_[i];
    pointList[j+2] = this->points_z_[i];
  }

  // TODO: debug print
  std::cerr << "Calling Qhull..." << std::endl;
  //orgQhull::RboxPoints rbox("D3");
  orgQhull::Qhull qhull("input", 3, POINT_LIST_SIZE, pointList, "");
  //qhull.runQhull(rbox, "");
  qhull.outputQhull();

  if ( qhull.hasQhullMessage() )
  {
    std::cerr << "\nResults of qhull\n" << qhull.qhullMessage();
    qhull.clearQhullMessage();
  }
  orgQhull::QhullFacetList facets = qhull.facetList();
  std::cout << "\nFacets created by Qhull::runQhull()\n" << facets;

  delete [] pointList;
  // TODO: same normal calc?
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