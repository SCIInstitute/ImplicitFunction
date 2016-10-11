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

#include <vector>
#include <cstdio>
#include <algorithm>
#include <cmath>

#include "ScatteredData.h"
#include "vec3.h"
#include "ETSP.h"
#include "2DConvexHull.h"

using std::vector;

//int ScatteredData::myAxis = 2;
ScatteredData::ScatteredData()
{
}

ScatteredData::ScatteredData(const std::vector<double>& points_x,
                             const std::vector<double>& points_y,
                             const std::vector<double>& points_z,
                             const std::vector<double>& func,
                             const std::vector<axis_t>& axisInfo):
  axisInformation_(axisInfo),
  updatedAxisInformation_(axisInfo)
{
  setData(points_x, points_y, points_z, func);
}

void ScatteredData::setData(const std::vector<double>& points_x,
                            const std::vector<double>& points_y,
                            const std::vector<double>& points_z,
                            const std::vector<double>& func)
{
  this->surfacePoints_[X] = points_x; this->surfacePoints_[Y] = points_y; this->surfacePoints_[Z] = points_z; this->fnc_ = func;

  //for(int i=0; i<a.size(); i++)
  //{
  //	printf("%d %d\n", this->surfacePoints_[0][i], a[i]); fflush(stdout);
  //}
  this->origSize_ = points_x.size();

// centroid only used for finding spherical normal (code not used)
//  vec3 temp(0,0,0);
//  for (int i = 0; i < this->origSize_; i++)
//  {
//    for(int j = 0; j < 3; j++)
//    {
//      temp[j] += this->surfacePoints_[j][i];
//    }
//  }
//  for (int j = 0; j < 3; j++)
//    this->centroid_[j] = temp[j] / this->origSize_;
}

//void ScatteredData::computeOrdering()
//{
//  SDmultisort();
//  //printf("Sorted\n");
//  //for(int i=0; i<inputData_.size(); i++)
//  //	printf("%d %lf %lf %lf\n", i, inputData_[i][0], inputData_[i][0], inputData_[i][2]);
//
//  int count = 0;
//  //printf("%d %d\n", inputData_.size(), this->surfacePoints_[0].size()); fflush(stdout);
//  while(count != this->surfacePoints_[0].size())
//  {
//    int start = count;
//    //printf("%lf %d\n", inputData_[count][myAxis], this->surfacePoints_[0].size());
//    while(count != this->surfacePoints_[0].size())
//    {
//      //printf("%lf %lf %lf\n", inputData_[count][myAxis], inputData_[start][myAxis],fabs(inputData_[count][myAxis]-inputData_[start][myAxis]));
//      count++;
//      axis_t sortingAxis = this->axisInformation_[start];
//      if(this->axisInformation_[count]!=this->axisInformation_[start])
//        break;
//      if(fabs(inputData_[count][sortingAxis]-inputData_[start][sortingAxis]) > 1e-6)
//        break;
//    }
//    int end = count-1;
//    printf("BCD Start: %d End %d\n", start, end);
//    int myAxis;
//    switch(this->axisInformation_[start])
//    {
//      case X: myAxis = 0; break;
//      case Y: myAxis = 1; break;
//      case Z: myAxis = 2; break;
//        
//    }
//    ETSP reorder(inputData_, start, end, myAxis);
//    for(int i=start; i<=end; i++)
//    {
//      inputData_[i][0] = this->surfacePoints_[0][i]=reorder.data[reorder.order[i-start]][0];
//      inputData_[i][1] = this->surfacePoints_[1][i]=reorder.data[reorder.order[i-start]][1];
//      inputData_[i][2] = this->surfacePoints_[2][i]=reorder.data[reorder.order[i-start]][2];
//    }
//  }
//  //for(int i=0; i<inputData_.size(); i++)
//  //	printf("%lf %lf %lf\n", inputData_[i][0], inputData_[i][0], inputData_[i][2]);
//}

void ScatteredData::SDsort()
{
  vec3Sorter sortObject;
  sortObject.axisToSort=this->axisInformation_[0];

  for (int i = 0; i < this->surfacePoints_[0].size(); i++)
  {
    vec3 a(this->surfacePoints_[X][i], this->surfacePoints_[Y][i], this->surfacePoints_[Z][i]);
    this->inputData_.push_back(a);
  }
  std::sort(this->inputData_.begin(), this->inputData_.end(), sortObject);
}

void ScatteredData::SDmultisort()
{
  vec3Sorter sortObject;
  std::vector<vec3> dataX, dataY, dataZ;
  for (int i = 0; i < this->surfacePoints_[0].size(); i++)
  {
    vec3 a(this->surfacePoints_[X][i], this->surfacePoints_[Y][i], this->surfacePoints_[Z][i]);
    switch(this->axisInformation_[i])
    {
      case X:
        dataX.push_back(a);
        break;
      case Y:
        dataY.push_back(a);
        break;
      case Z:
        dataZ.push_back(a);
        break;
    }
  }

  this->inputData_.clear();
  sortObject.axisToSort = X;
  std::sort(dataX.begin(), dataX.end(), sortObject);

  sortObject.axisToSort = Y;
  std::sort(dataY.begin(), dataY.end(), sortObject);

  sortObject.axisToSort = Z;
  std::sort(dataZ.begin(), dataZ.end(), sortObject);

  this->axisInformation_.clear();
  for(int i = 0; i < dataX.size(); i++)
  {
    this->inputData_.push_back(dataX[i]);
    this->axisInformation_.push_back(X);
  }

  for(int i = 0; i < dataY.size(); i++)
  {
    this->inputData_.push_back(dataY[i]);
    this->axisInformation_.push_back(Y);
  }

  for(int i = 0; i < dataZ.size(); i++)
  {
    this->inputData_.push_back(dataZ[i]);
    this->axisInformation_.push_back(Z);
  }
}


void ScatteredData::compute2DHull()
{
  std::vector<int> reorderedData;

  // TODO: debug print
  for(int i = 0; i < this->surfacePoints_[0].size(); i++)
    printf("%d %lf %lf %lf\n", i, this->surfacePoints_[0][i], this->surfacePoints_[1][i], this->surfacePoints_[2][i]);

  SDmultisort();
  printf("Sorted\n");

  // TODO: debug print
  for (int i = 0; i < this->inputData_.size(); i++)
    printf("%d %lf %lf %lf\n", i, this->inputData_[i][0], this->inputData_[i][1], this->inputData_[i][2]);
  
  int count = 0;
  //printf("%d %d\n", this->inputData_.size(), this->surfacePoints_[0].size()); fflush(stdout);
  while ( count != this->surfacePoints_[0].size() )
  {
    int start = count;
    //printf("%lf %d\n", this->inputData_[count][myAxis], this->surfacePoints_[0].size());
    while ( count != this->surfacePoints_[0].size() )
    {
      //printf("%lf %lf %lf\n", this->inputData_[count][myAxis], this->inputData_[start][myAxis],fabs(this->inputData_[count][myAxis]-this->inputData_[start][myAxis]));
      count++;
      axis_t sortingAxis = this->axisInformation_[start];
      if (this->axisInformation_[count] != this->axisInformation_[start])
        break;

      if (fabs(this->inputData_[count][sortingAxis] - this->inputData_[start][sortingAxis]) > 1e-6)
        break;
    }
    int end = count-1;

    // TODO: debug print
    printf("ABC Start: %d End %d\n", start, end);

    int myAxis = -1;
    switch(this->axisInformation_[start])
    {
      case X: myAxis = 0; break;
      case Y: myAxis = 1; break;
      case Z: myAxis = 2; break;
    }

    vector<vec3> inPoints;
    inPoints.resize(end-start+1);
    for (int i = start; i <= end; i++)
    {
      inPoints[i-start][0] = this->inputData_[i][0];
      inPoints[i-start][1] = this->inputData_[i][1];
      inPoints[i-start][2] = this->inputData_[i][2];
    }

    //printf("inPoints written\n");
    vector<int> reorder = getConvexHull(inPoints, myAxis);

    // TODO: simplify data transfer

    printf("Convex Hull found\n");
    // indices of points in convex hull
    for(int i = 0; i < reorder.size(); i++)
    {
      printf("%d ", reorder[i]+start);
      reorderedData.push_back(reorder[i]+start);
    }
    printf("\n");
  }

  std::sort( reorderedData.begin(), reorderedData.end() );

  //rewrite everything
  std::vector<double> newx[3], tmp[3], newfnc;

  //std::vector<vec3> myNewData;
  //std::vector<axis_t> newAxisInformation;
  //printf("Rewriting everything\n");

  this->updatedAxisInformation_.clear();

  for(int i = 0; i < reorderedData.size(); i++)
  {
    const int j = reorderedData[i];
    newx[X].push_back(this->inputData_[j][0]);
    newx[Y].push_back(this->inputData_[j][1]);
    newx[Z].push_back(this->inputData_[j][2]);

    newfnc.push_back(this->fnc_[j]);

    this->convexHullData_.push_back( this->inputData_[j] ); // TODO: never actually used except for debug prints...
    this->updatedAxisInformation_.push_back( this->axisInformation_[j] );
  }

  this->surfacePoints_[X].clear();
  this->surfacePoints_[Y].clear();
  this->surfacePoints_[Z].clear();

  this->surfacePoints_[X] = newx[X];
  this->surfacePoints_[Y] = newx[Y];
  this->surfacePoints_[Z] = newx[Z];

  this->fnc_.clear();
  this->fnc_ = newfnc;

  for (int i = 0; i < this->inputData_.size(); ++i)
  {
    auto it = std::find( reorderedData.begin(), reorderedData.end(), i);
    if ( it == reorderedData.end() )
    {
      tmp[X].push_back(this->inputData_[i][0]);
      tmp[Y].push_back(this->inputData_[i][1]);
      tmp[Z].push_back(this->inputData_[i][2]);
    }
  }

  this->leftovers_[X].clear();
  this->leftovers_[Y].clear();
  this->leftovers_[Z].clear();

  this->leftovers_[X] = tmp[X];
  this->leftovers_[Y] = tmp[Y];
  this->leftovers_[Z] = tmp[Z];

  // TODO: this appears to replace the original input data with convex hull
  // points
  // Points are being ignored
  //this->inputData_.clear();
  //this->inputData_=myNewData;
  //
  // TODO: overwrites original axis information...
  // Just reordered?
  //this->axisInformation_.clear();
  //this->axisInformation_ = newAxisInformation;

  printf("Convex hull points\n");
  for(int i = 0; i < this->surfacePoints_[0].size(); i++)
  	printf("%d %lf %lf %lf\n", i, this->surfacePoints_[0][i], this->surfacePoints_[1][i], this->surfacePoints_[2][i]);

//  for (int i = 0; i < this->convexHullData_.size(); i++)
//    printf("%d %lf %lf %lf\n", i, this->convexHullData_[i][0], this->convexHullData_[i][1], this->convexHullData_[i][2]);

  printf("Points inside convex hull\n");
  for (int i = 0; i < this->leftovers_[0].size(); i++)
    printf("%d %lf %lf %lf\n", i, this->leftovers_[0][i], this->leftovers_[1][i], this->leftovers_[2][i]);

//  this->origSize_ = this->inputData_.size();

  // TODO: this gets overwritten in RBFInterface???
  this->origSize_ = this->convexHullData_.size();
}
