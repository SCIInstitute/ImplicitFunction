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

const double RBF::EPSILON = 1.0e-2;

RBF::RBF(ScatteredData *myData, Kernel myKernel) :
  completeData_(myData),
  kernel_(myKernel),
  acceleration_(None),
  dataReduction_(All)
{
}

RBF::~RBF()
{
  delete this->data_;
  delete this->fmm_;
}

void RBF::setAcceleration(Acceleration myAcceleration)
{
  this->acceleration_ = myAcceleration;
  switch(this->acceleration_)
  {
    case FastMultipole:
      // TODO: memory leak?
      fmm_ = new FMM();
      break;
    default:
      break;
  }
}

void RBF::setDataReduction(DataReduction myDataReduction)
{
  this->dataReduction_ = myDataReduction;
}

void RBF::computeFunction()
{
  // TODO: memory leak?
  data_ = new ScatteredData;
  switch(this->dataReduction_)
  {
    case All:
      data_->setData(this->completeData_->surfacePoints_[0],
                     this->completeData_->surfacePoints_[1],
                     this->completeData_->surfacePoints_[2],
                     this->completeData_->fnc_);
      computeFunctionForData();
      break;
      
    case Random:
      bool *added;
      const int N = this->completeData_->fnc_.size();
      added = new bool[N];
      printf("%d\n", N);
      
      for (int i = 0; i < N; i++)
        added[i] = false;

      // TODO: magic number
      for (int i = 0; i < 25; i++)
      {
        int j = rand() % N;
        if( added[j] )
        {
          i--;
          continue;
        }
        added[j] = true;
        this->data_->surfacePoints_[0].push_back(this->completeData_->surfacePoints_[0][j]);
        this->data_->surfacePoints_[1].push_back(this->completeData_->surfacePoints_[1][j]);
        this->data_->surfacePoints_[2].push_back(this->completeData_->surfacePoints_[2][j]);
        this->data_->fnc_.push_back(this->completeData_->fnc_[j]);
        //printf("%d %lf %lf %lf %lf\n", j, this->completeData_->surfacePoints_[0][j],this->completeData_->surfacePoints_[1][j],this->completeData_->surfacePoints_[2][j],this->completeData_->fnc_[j]);
      }

      vector<pair<double, int> > error;
      bool smallError = false;
      while (! smallError)
      {
        computeFunctionForData();
        computeErrorForData(error);
        std::sort(error.begin(), error.end());
        std::reverse(error.begin(), error.end());
        
        printf("Largest error: %lf\n", error[0].first);
        if (error[0].first < EPSILON)
        {
          smallError = true;
          continue;
        }
        for (int k = 0; k < 5; k++)
        {
          //printf("Error %d: %lf\n", k, error[k].first);
          if (error[k].first > EPSILON || error[k].first != error[k].first)
          {
            int j = error[k].second;
            //printf("Adding data point %d\n", j);
            added[j] = true;
            this->data_->surfacePoints_[0].push_back(this->completeData_->surfacePoints_[0][j]);
            this->data_->surfacePoints_[1].push_back(this->completeData_->surfacePoints_[1][j]);
            this->data_->surfacePoints_[2].push_back(this->completeData_->surfacePoints_[2][j]);
            this->data_->fnc_.push_back(this->completeData_->fnc_[j]);
          }
        }
      }
      //printf("Total no. of data_ point: %d\n",  this->data_->fnc_.size()); fflush(stdout);
      delete [] added;
      break;
  }
}

void RBF::computeFunctionForData()
{
  switch(this->acceleration_)
  {
    case FastMultipole:
      fmmBuildTree();
      break;
    case None:
    default:
      // TODO: move to function
      const int N = this->data_->fnc_.size();
      printf("Solving linear equations: \n"); fflush(stdout);
      LinearSolver rbfSolver;
      SparseMatrix rbfMatrix(N);
      //printf("Constructing matrix...\n"); fflush(stdout);
      for (int i = 0; i < N; i++)
      {
        for (int j = 0; j < N; j++)
        {
          //printf("%d %d ", i,j); fflush(stdout);
          double val = computeKernel(i, j);
          //printf("%lf\n", val); fflush(stdout);
          rbfMatrix.push_back(i, j, val);
        }
      }
      //printf("Done\n"); fflush(stdout);
      rbfSolver.setMatrix(&rbfMatrix);
      printf("Running BiCGSTAB Iterations ... "); fflush(stdout);
      rbfSolver.biCGStab(this->data_->fnc_, this->coeff_);

      break;
  }
}

double RBF::computeValue(const vec3& x)
{
  switch(this->acceleration_)
  {
    case FastMultipole:
      return fmmComputeValue(x);
    case None:
    default:
      // TODO: move to function
      double sum = 0;
      for(int i = 0; i < this->coeff_.size(); i++)
        sum += this->coeff_[i]*computeKernel(i, x);
      return sum;
  }
}

void RBF::computeErrorForData(vector<pair<double, int> > &error)
{
  const int N = this->completeData_->fnc_.size();
  error.clear();
  for (int i = 0; i < N; i++)
  {
    vec3 x(this->completeData_->surfacePoints_[0][i],
           this->completeData_->surfacePoints_[1][i],
           this->completeData_->surfacePoints_[2][i]);
    double err = this->completeData_->fnc_[i]-computeValue(x);
    error.push_back( std::make_pair(err, i) );
  }
}

double RBF::computeKernel(int i, int j)
{
  double r = sqrt( (this->data_->surfacePoints_[0][i] - this->data_->surfacePoints_[0][j]) *
                   (this->data_->surfacePoints_[0][i] - this->data_->surfacePoints_[0][j]) +  // x
                   (this->data_->surfacePoints_[1][i] - this->data_->surfacePoints_[1][j]) *
                   (this->data_->surfacePoints_[1][i] - this->data_->surfacePoints_[1][j]) +  // y
                   (this->data_->surfacePoints_[2][i] - this->data_->surfacePoints_[2][j]) *
                   (this->data_->surfacePoints_[2][i] - this->data_->surfacePoints_[2][j]) ); // z

  return computeRadialFunction(r);
  
}

double RBF::computeKernel(int i, const vec3& b)
{
  double r = sqrt( (this->data_->surfacePoints_[0][i] - b[0])*(this->data_->surfacePoints_[0][i] - b[0]) +  // x
                   (this->data_->surfacePoints_[1][i] - b[1])*(this->data_->surfacePoints_[1][i] - b[1]) +  // y
                   (this->data_->surfacePoints_[2][i] - b[2])*(this->data_->surfacePoints_[2][i] - b[2]) ); // z

  return computeRadialFunction(r);
}

double RBF::computeRadialFunction(double r)
{
  const double C = 0.1;
  const double SCALE = 0.01;

  switch(kernel_)
  {
    case Gaussian:
      r = r * SCALE;
      return 1.0/sqrt(r*r + C*C);
      break;
    case ThinPlate:
      return r*r*log(r+C);
      break;
    case MultiQuadratic:
      return sqrt(r*r + C*C);
    default:
      return r;
      break;
  }
  return 0;
}

//FMM Codes

void RBF::fmmBuildTree()
{
  vector<int> myIndices;
  const int N = this->data_->surfacePoints_[0].size();

  for (int i = 0; i < N; i++)
    myIndices.push_back(i);

  fmm_->tree = new BHNode();
  fmm_->tree->box_.min_ = vec3::zero;

  vec3 corner(100, 100, 100);
  fmm_->tree->box_.max_ = corner;
  fmmBuildTree(myIndices, fmm_->tree);
  //printf("Tree Built\n");
  //fmmPrintTree(fmm_->tree, 0);
}

void RBF::fmmPrintTree(BHNode *myNode, int stack)
{
  if (stack > 4)
    return;

  for (int j = 0; j < stack; j++)
  {
    printf(" ");
  }

  printf("%d %d %d\n", myNode, myNode->index_, myNode->pts_.size());

  for (int i = 0; i < 8; i++)
  {
    if (myNode->nodes[i] != nullptr)
      fmmPrintTree(myNode->nodes[i], stack+1);
  }
}

void RBF::fmmBuildTree(vector<int> &myPoints, BHNode *myNode)
{
  const double SMALL_EPSILON = 1.0e-6;
  //printf("[%lf %lf %lf] [%lf %lf %lf] %d\n", myNode->box_.min_[0], myNode->box_.min_[1], myNode->box_.min_[2], myNode->box_.max_[0], myNode->box_.max_[1], myNode->box_.max_[2], myPoints.size());
  vector<int> children[8];
  const int N = myPoints.size();
  
  myNode->index_ = fmm_->numOfNodes;
  fmm_->numOfNodes += 1;
  fmm_->nodePointer.push_back(myNode);
  
  myNode->mass_ = N;
  myNode->leaf_ = (N <= 1) ? true : false;
  myNode->center_ = vec3::zero;
  
  myNode->coeff_ = 1; //REPLACE
  //add all the coefficients
  
  for (int i = 0; i < N; i++)
    myNode->pts_.push_back( myPoints[i] );
  
  for (int i = 0; i < N; i++)
  {
    vec3 location(this->data_->surfacePoints_[0][myPoints[i]], this->data_->surfacePoints_[1][myPoints[i]], this->data_->surfacePoints_[2][myPoints[i]]);
    myNode->center_ = myNode->center_ + (location/N);
  }
  
  if (N == 1)
  {
    for (int i = 0; i < 8; i++)
      myNode->nodes[i] = nullptr;

    return;
  }

  if (length(myNode->box_.max_ - myNode->box_.min_) < SMALL_EPSILON)
  {
    for (int i = 0; i < 8; i++)
      myNode->nodes[i] = nullptr;

    return;
  }

  vec3 mid( (myNode->box_.getMin() + myNode->box_.getMax())/2 );
  for (int i = 0; i < N; i++)
  {
    int octant = 0;
    //FIND OCTANTS
    if (this->data_->surfacePoints_[0][myPoints[i]] > mid[0])
      octant += 1;
    if (this->data_->surfacePoints_[1][myPoints[i]] > mid[1])
      octant += 2;
    if (this->data_->surfacePoints_[2][myPoints[i]] > mid[2])
      octant += 4;
    
    //printf("%d %d %d %lf %lf %lf\n", i,octant, myPoints[i], this->data_->surfacePoints_[0][myPoints[i]],this->data_->surfacePoints_[1][myPoints[i]], this->data_->surfacePoints_[2][myPoints[i]]);
    
    children[octant].push_back(myPoints[i]);
  }

  for (int i = 0; i < 8; i++)
  {
    if (children[i].size() == 0)
    {
      myNode->nodes[i] = nullptr;
      continue;
    }

    myNode->nodes[i] = new BHNode();

    int hash = i;
    for (int j = 0; j < 3; j++)
    {
      if (hash % 2 == 0)
      {
        myNode->nodes[i]->box_.min_[j] = myNode->box_.min_[j];
        myNode->nodes[i]->box_.max_[j] = mid[j];
      }
      else
      {
        myNode->nodes[i]->box_.min_[j] = mid[j];
        myNode->nodes[i]->box_.max_[j] = myNode->box_.max_[j];
      }
      hash /= 2;
    }

    fmmBuildTree(children[i], myNode->nodes[i]);
  }
}

double RBF::fmmComputeValue(const vec3& x)
{
  double val = fmmComputeValueRecurse(x, fmm_->tree);
  return val;
}


double RBF::fmmComputeValueRecurse(const vec3& x, BHNode *myNode)
{
  double val = 0;

  // TODO: div by zero?
  if (length(myNode->center_ - x)/length(myNode->box_.min_ - myNode->box_.max_) < 0)  //far away
  {
    val = myNode->coeff_ * fmmComputeKernel(x, myNode);
  }
  else if (! myNode->leaf_) // leaf
  {
    for (int i = 0; i < 8; i++)
    {
      if (myNode->nodes[i] != nullptr)
      {
        val += fmmComputeValueRecurse(x, myNode->nodes[i]);
      }
    }
  }
  else //close, but not a leaf
  {
    const int N = myNode->pts_.size();
    for (int i = 0; i < N; i++)
    {
      val += this->coeff_[myNode->pts_[i]] * computeKernel(myNode->pts_[i], x);
    }
  }
  return val;
}

double RBF::fmmComputeKernel(const vec3& b, BHNode *myNode)
{
  double r = length(myNode->center_ - b);
  
  return computeRadialFunction(r);
}