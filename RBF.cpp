#include "ScatteredData.h"
#include "RBF.h"
#include "SparseMatrix.h"
#include "vec3.h"
#include "LinearSolver.h"
#include "FMM.h"

#include <cmath>
#include <cstdio>

// STL Includes
#include <vector>
#include <utility>
#include <algorithm>

using std::vector;
using std::pair;


RBF::RBF(ScatteredData *myData, Kernel myKernel) :
  kernel_(myKernel),
  completeData_(myData)
{
}

void RBF::setKernel(Kernel myKernel)
{
  kernel_ = myKernel;
}

void RBF::setData(ScatteredData *myData)
{
  completeData_ = myData;
}

void RBF::setAcceleration(Acceleration myAcceleration)
{
  acceleration = myAcceleration;
  switch(acceleration)
  {
    case FastMultipole:
      fmm_ = new FMM();
      break;
    default:
      break;
  }
}

void RBF::setDataReduction(DataReduction myDataReduction)
{
  dataReduction = myDataReduction;
}

void RBF::computeFunction()
{
  // TODO: memory leak?
  data_ = new ScatteredData;
  switch(dataReduction)
  {
    case All:
      data_->setData(completeData_->x[0], completeData_->x[1], completeData_->x[2], completeData_->fnc);
      computeFunctionForData();
      break;
      
    case Random:
      bool *added;
      int n = completeData_->fnc.size();
      added = new bool[n];
      printf("%d\n", n);
      
      for(int i=0; i<n; i++)
        added[i]=false;
      
      for(int i=0; i<25; i++)
      {
        int j = rand()%n;
        if(added[j])
        {
          i--;
          continue;
        }
        added[j]=true;
        data_->x[0].push_back(completeData_->x[0][j]);
        data_->x[1].push_back(completeData_->x[1][j]);
        data_->x[2].push_back(completeData_->x[2][j]);
        data_->fnc.push_back(completeData_->fnc[j]);
        printf("%d %lf %lf %lf %lf\n", j, completeData_->x[0][j],completeData_->x[1][j],completeData_->x[2][j],completeData_->fnc[j]);
      }
      vector<pair<double, int> > error;
      bool smallError=false;
      while(!smallError)
      {
        computeFunctionForData();
        computeErrorForData(error);
        std::sort(error.begin(), error.end());
        std::reverse(error.begin(), error.end());
        
        printf("Largest error: %lf\n", error[0].first);
        if(error[0].first < 1e-2)
        {
          smallError=true;
          continue;
        }
        for(int k=0; k<5; k++)
        {
          //printf("Error %d: %lf\n", k, error[k].first);
          if(error[k].first>1e-2 || error[k].first!=error[k].first)
          {
            int j = error[k].second;
            printf("Adding data_ point %d\n", j);
            added[j]=true;
            data_->x[0].push_back(completeData_->x[0][j]);
            data_->x[1].push_back(completeData_->x[1][j]);
            data_->x[2].push_back(completeData_->x[2][j]);
            data_->fnc.push_back(completeData_->fnc[j]);
          }
        }
        
      }
      printf("Total no. of data_ point: %d\n",  data_->fnc.size()); fflush(stdout);
  }
}

void RBF::computeFunctionForData()
{
  switch(acceleration)
  {
    case FastMultipole:
      fmmComputeFunction();
    case None:
    default:
      // TODO: move to function
      int n = data_->fnc.size();
      printf("Solving linear equations: \n"); fflush(stdout);
      LinearSolver rbfSolver;
      SparseMatrix rbfMatrix(n);
      printf("Constructing matrix ... "); fflush(stdout);
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<n; j++)
        {
          //printf("%d %d ", i,j); fflush(stdout);
          double val = computeKernel(i,j);
          //printf("%lf\n", val); fflush(stdout);
          rbfMatrix.push_back(i,j,val);
        }
      }
      printf("Done\n"); fflush(stdout);
      rbfSolver.setMatrix(&rbfMatrix);
      printf("Running BiCGSTAB Iterations ... "); fflush(stdout);
      rbfSolver.biCGStab(data_->fnc, coeff);
      printf("Done\n"); fflush(stdout);
  }
}

double RBF::computeValue(vec3 x)
{
  switch(acceleration)
  {
    case FastMultipole:
      return fmmComputeValue(x);
    case None:
    default:
      // TODO: move to function
      double sum=0;
      for(int i=0; i<coeff.size(); i++)
        sum+=coeff[i]*computeKernel(i, x);
      return sum;
  }
}

void RBF::computeErrorForData(vector<pair<double, int> > &error)
{
  int n = completeData_->fnc.size();
  error.clear();
  for (int i=0; i<n; i++)
  {
    vec3 x(completeData_->x[0][i],completeData_->x[1][i],completeData_->x[2][i]);
    double err = completeData_->fnc[i]-computeValue(x);
    error.push_back(std::make_pair(err,i));
  }
}

double RBF::computeKernel(int i, int j)
{
  double r = sqrt( (data_->x[0][i] - data_->x[0][j])*(data_->x[0][i] - data_->x[0][j]) +
                  (data_->x[1][i] - data_->x[1][j])*(data_->x[1][i] - data_->x[1][j]) +
                  (data_->x[2][i] - data_->x[2][j])*(data_->x[2][i] - data_->x[2][j]) );
  
  return computeRadialFunction(r);
  
}

double RBF::computeKernel(int i, vec3 b)
{
  double r = sqrt( (data_->x[0][i] - b[0])*(data_->x[0][i] - b[0]) +
                  (data_->x[1][i] - b[1])*(data_->x[1][i] - b[1]) +
                  (data_->x[2][i] - b[2])*(data_->x[2][i] - b[2]) );
  
  return computeRadialFunction(r);
}

double RBF::computeRadialFunction(double r)
{
  double c = 0.1;
  switch(kernel_)
  {
    case Gaussian:
      r = r*0.01;
      return 1.0/sqrt(r*r + c*c);
      break;
    case ThinPlate:
      return r*r*log(r+c);
      break;
    case MultiQuadratic:
      return sqrt(r*r + c*c);
    default:
      return r;
      break;
  }
  return 0;
}







//FMM Codes
void RBF::fmmComputeFunction()
{
  fmmBuildTree();
}

void RBF::fmmBuildTree()
{
  vector<int> myIndices;
  int n = data_->x[0].size();
  for(int i=0; i<n; i++)
    myIndices.push_back(i);
  fmm_->tree = new BHNode();
  fmm_->tree->box.min = vec3::zero;
  vec3 corner(100,100,100);
  fmm_->tree->box.max = corner;
  fmmBuildTree(myIndices, fmm_->tree);
  //printf("Tree Built\n");
  //fmmPrintTree(fmm_->tree, 0);
}

void RBF::fmmPrintTree(BHNode *myNode, int stack)
{
  if (stack>4)
    return;
  for(int j=0; j<stack; j++)
  {
    printf(" ");
  }
  printf("%d %d %d\n", myNode, myNode->index, myNode->pts.size());
  for(int i=0; i<8; i++)
  {
    if(myNode->nodes[i]!=NULL)
      fmmPrintTree(myNode->nodes[i], stack+1);
  }
}

void RBF::fmmBuildTree(vector<int> &myPoints, BHNode *myNode)
{
  
  //printf("[%lf %lf %lf] [%lf %lf %lf] %d\n", myNode->box.min[0], myNode->box.min[1], myNode->box.min[2], myNode->box.max[0], myNode->box.max[1], myNode->box.max[2], myPoints.size());
  vector<int> children[8];
  int n = myPoints.size();
  
  myNode->index = fmm_->numOfNodes;
  fmm_->numOfNodes += 1;
  fmm_->nodePointer.push_back(myNode);
  
  myNode->mass = n;
  myNode->leaf = (n<=1)?true:false;
  myNode->center = vec3::zero;
  
  myNode->coeff = 1; //REPLACE
  //add all the coefficients
  
  
  for(int i=0; i<myPoints.size(); i++)
    myNode->pts.push_back(myPoints[i]);
  
  for(int i=0; i<n; i++)
  {
    vec3 location(data_->x[0][myPoints[i]], data_->x[1][myPoints[i]],data_->x[2][myPoints[i]]);
    myNode->center = myNode->center + (location/n);
  }
  
  if (n==1)
  {
    for(int i=0; i<8; i++)
      myNode->nodes[i]=NULL;
    return;
  }
  if(length(myNode->box.max - myNode->box.min) < 1e-6)
  {
    for(int i=0; i<8; i++)
      myNode->nodes[i]=NULL;
    return;
  }
  
  vec3 mid((myNode->box.getMin()+myNode->box.getMax())/2);
  for(int i=0; i<n; i++)
  {
    int octant=0;
    //FIND OCTANTS
    if(data_->x[0][myPoints[i]] > mid[0])
      octant+=1;
    if(data_->x[1][myPoints[i]] > mid[1])
      octant+=2;
    if(data_->x[2][myPoints[i]] > mid[2])
      octant+=4;
    
    //printf("%d %d %d %lf %lf %lf\n", i,octant, myPoints[i], data_->x[0][myPoints[i]],data_->x[1][myPoints[i]], data_->x[2][myPoints[i]]);
    
    children[octant].push_back(myPoints[i]);
  }
  for(int i=0; i<8; i++)
  {
    if (children[i].size() == 0)
    {
      myNode->nodes[i] = NULL;
      continue;
    }
    
    myNode->nodes[i] = new BHNode();
    
    int hash = i;
    for(int j=0; j<3; j++)
    {
      if(hash%2==0)
      {
        myNode->nodes[i]->box.min[j]=myNode->box.min[j];
        myNode->nodes[i]->box.max[j]=mid[j];
      }
      else
      {
        myNode->nodes[i]->box.min[j]=mid[j];
        myNode->nodes[i]->box.max[j]=myNode->box.max[j];
      }
      hash /= 2;
    }
    
    
    fmmBuildTree(children[i], myNode->nodes[i]);
  }
}

double RBF::fmmComputeValue(vec3 x)
{
  double val = fmmComputeValueRecurse(x, fmm_->tree);
  return val;
}


double RBF::fmmComputeValueRecurse(vec3 x, BHNode *myNode)
{
  double val=0;
  
  if(length(myNode->center - x)/length(myNode->box.min-myNode->box.max) < 0)  //far away
  {
    val = myNode->coeff*fmmComputeKernel(x, myNode);
  }
  else if (!myNode->leaf) // leaf
  {
    for(int i=0; i<8; i++)
    {
      if(myNode->nodes[i] != nullptr)
      {
        val+=fmmComputeValueRecurse(x, myNode->nodes[i]);
      }
    }
  }
  else //close, but not a leaf
  {
    int n=myNode->pts.size();
    for(int i=0; i<n; i++)
    {
      val+=coeff[myNode->pts[i]]*computeKernel(myNode->pts[i], x);
    }
  }
  return val;
}

double RBF::fmmComputeKernel(vec3 b, BHNode *myNode)
{
  double r = length(myNode->center - b);
  
  return computeRadialFunction(r);
}


double RBF::fmmComputeKernel(BHNode *myNode, vec3 b)
{
  return fmmComputeKernel(b, myNode);
}
