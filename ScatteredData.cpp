#include <vector>
#include <cstdio>
#include <algorithm>
#include "ScatteredData.h"
#include "vec3.h"
#include "ETSP.h"
#include "2DConvexHull.h"
#include <cmath>

using std::vector;

//int ScatteredData::myAxis = 2;
ScatteredData::ScatteredData()
{
}

ScatteredData::ScatteredData(const vector<double>& a,
                             const vector<double>& b,
                             const vector<double>& c,
                             const vector<double>& d)
{
  setData(a, b, c, d);
}

void ScatteredData::setData(const vector<double>& a,
                            const vector<double>& b,
                            const vector<double>& c,
                            const vector<double>& d)
{
  this->x_[0] = a; this->x_[1] = b; this->x_[2] = c; fnc = d;
  //for(int i=0; i<a.size(); i++)
  //{
  //	printf("%d %d\n", this->x_[0][i], a[i]); fflush(stdout);
  //}
  origSize = a.size();
  vec3 temp(0,0,0);
  for (int i = 0; i < origSize; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      temp[j] += this->x_[j][i];
    }
  }
  for (int j = 0; j < 3; j++)
    centroid[j] = temp[j] / origSize;
}

//void ScatteredData::computeOrdering()
//{
//  SDmultisort();
//  //printf("Sorted\n");
//  //for(int i=0; i<inputData_.size(); i++)
//  //	printf("%d %lf %lf %lf\n", i, inputData_[i][0], inputData_[i][0], inputData_[i][2]);
//
//  int count = 0;
//  //printf("%d %d\n", inputData_.size(), this->x_[0].size()); fflush(stdout);
//  while(count != this->x_[0].size())
//  {
//    int start = count;
//    //printf("%lf %d\n", inputData_[count][myAxis], this->x_[0].size());
//    while(count != this->x_[0].size())
//    {
//      //printf("%lf %lf %lf\n", inputData_[count][myAxis], inputData_[start][myAxis],fabs(inputData_[count][myAxis]-inputData_[start][myAxis]));
//      count++;
//      axis_t sortingAxis = axisInformation[start];
//      if(axisInformation[count]!=axisInformation[start])
//        break;
//      if(fabs(inputData_[count][sortingAxis]-inputData_[start][sortingAxis]) > 1e-6)
//        break;
//    }
//    int end = count-1;
//    printf("BCD Start: %d End %d\n", start, end);
//    int myAxis;
//    switch(axisInformation[start])
//    {
//      case X: myAxis = 0; break;
//      case Y: myAxis = 1; break;
//      case Z: myAxis = 2; break;
//        
//    }
//    ETSP reorder(inputData_, start, end, myAxis);
//    for(int i=start; i<=end; i++)
//    {
//      inputData_[i][0] = this->x_[0][i]=reorder.data[reorder.order[i-start]][0];
//      inputData_[i][1] = this->x_[1][i]=reorder.data[reorder.order[i-start]][1];
//      inputData_[i][2] = this->x_[2][i]=reorder.data[reorder.order[i-start]][2];
//    }
//  }
//  //for(int i=0; i<inputData_.size(); i++)
//  //	printf("%lf %lf %lf\n", inputData_[i][0], inputData_[i][0], inputData_[i][2]);
//}

void ScatteredData::SDsort()
{
  vec3Sorter sortObject;
  sortObject.axisToSort=axisInformation[0];
  for(int i=0; i<this->x_[0].size(); i++)
  {
    vec3 a(this->x_[0][i],this->x_[1][i],this->x_[2][i]);
    inputData_.push_back(a);
  }
  std::sort(inputData_.begin(), inputData_.end(), sortObject);
}

void ScatteredData::SDmultisort()
{
  vec3Sorter sortObject;
  std::vector<vec3> dataX, dataY, dataZ;
  for(int i=0; i<this->x_[0].size(); i++)
  {
    vec3 a(this->x_[0][i],this->x_[1][i],this->x_[2][i]);
    switch(axisInformation[i])
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
  inputData_.clear();
  sortObject.axisToSort=X;
  std::sort(dataX.begin(), dataX.end(), sortObject);
  sortObject.axisToSort=Y;
  std::sort(dataY.begin(), dataY.end(), sortObject);
  sortObject.axisToSort=Z;
  std::sort(dataZ.begin(), dataZ.end(), sortObject);
  axisInformation.clear();
  for(int i=0; i<dataX.size(); i++)
  {
    inputData_.push_back(dataX[i]);
    axisInformation.push_back(X);
  }
  for(int i=0; i<dataY.size(); i++)
  {
    inputData_.push_back(dataY[i]);
    axisInformation.push_back(Y);
  }
  for(int i=0; i<dataZ.size(); i++)
  {
    inputData_.push_back(dataZ[i]);
    axisInformation.push_back(Z);
  }
}


void ScatteredData::compute2DHull()
{
std::cerr << "ScatteredData::compute2DHull()" << std::endl;

  std::vector<int> reorderedData;
  for(int i = 0; i < this->x_[0].size(); i++)
    printf("%d %lf %lf %lf\n", i, this->x_[0][i], this->x_[1][i], this->x_[2][i]);
  SDmultisort();
  printf("Sorted\n");
  for(int i = 0; i < inputData_.size(); i++)
    printf("%d %lf %lf %lf\n", i, inputData_[i][0], inputData_[i][1], inputData_[i][2]);
  
  int count = 0;
  //printf("%d %d\n", inputData_.size(), this->x_[0].size()); fflush(stdout);
  while ( count != this->x_[0].size() )
  {
    int start = count;
    //printf("%lf %d\n", inputData_[count][myAxis], this->x_[0].size());
    while ( count != this->x_[0].size() )
    {
      //printf("%lf %lf %lf\n", inputData_[count][myAxis], inputData_[start][myAxis],fabs(inputData_[count][myAxis]-inputData_[start][myAxis]));
      count++;
      axis_t sortingAxis = axisInformation[start];
      if (axisInformation[count] != axisInformation[start])
        break;
      if (fabs(inputData_[count][sortingAxis] - inputData_[start][sortingAxis]) > 1e-6)
        break;
    }
    int end = count-1;
    printf("ABC Start: %d End %d\n", start, end);
    int myAxis;
    switch(axisInformation[start])
    {
      case X: myAxis = 0; break;
      case Y: myAxis = 1; break;
      case Z: myAxis = 2; break;
        
    }
    vector<vec3> inPoints;
    inPoints.resize(end-start+1);
    for(int i=start; i<=end; i++)
    {
      inPoints[i-start][0] = inputData_[i][0];
      inPoints[i-start][1] = inputData_[i][1];
      inPoints[i-start][2] = inputData_[i][2];
    }
    //printf("inPoints written\n");
    vector<int> reorder = getConvexHull(inPoints, myAxis);
    printf("Convex Hull found\n");
    for(int i = 0; i < reorder.size(); i++)
    {
      printf("%d ", i, reorder[i]+start);
      reorderedData.push_back(reorder[i]+start);
    }
    printf("\n");
  }
  //rewrite everything
  std::vector<double> newx[3], newfnc;
  //std::vector<vec3> myNewData;
  std::vector<axis_t> newAxisInformation;
  //printf("Rewriting everything\n");

  for(int i = 0; i < reorderedData.size(); i++)
  {
    int j = reorderedData[i];
    newx[0].push_back(inputData_[j][0]);
    newx[1].push_back(inputData_[j][1]);
    newx[2].push_back(inputData_[j][2]);
    newfnc.push_back(fnc[j]);
    convexHullData_.push_back(inputData_[j]);
    newAxisInformation.push_back(axisInformation[j]);
  }

  this->x_[0].clear();
  this->x_[1].clear();
  this->x_[2].clear();
  this->x_[0] = newx[0];
  this->x_[1] = newx[1];
  this->x_[2] = newx[2];
  fnc.clear();
  fnc = newfnc;

  // TODO: this appears to replace the original input data with convex hull
  // points
  // Points are being ignored
  //inputData_.clear();
  //inputData_=myNewData;

  axisInformation.clear();
  axisInformation = newAxisInformation;
  //for(int i=0; i<this->x_[0].size(); i++)
  //	printf("%d %lf %lf %lf\n", i, this->x_[0][i], this->x_[1][i], this->x_[2][i]);
  printf("Convex hull points\n");
  for(int i=0; i<inputData_.size(); i++)
    printf("%d %lf %lf %lf\n", i, inputData_[i][0], inputData_[i][1], inputData_[i][2]);
  origSize = inputData_.size();
}
