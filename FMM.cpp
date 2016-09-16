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
using std::vector;

FMM::FMM()
{
  theta=0.75;
}

FMM::~FMM()
{
  freeTheTree(tree);
}

void FMM::freeTheTree(BHNode *myNode)
{
  if ( myNode == nullptr )
    return;

  // TODO: magic number
  for(int i = 0; i < 8; i++)
    freeTheTree(myNode->nodes[i]);
}

void FMM::setTheta(double myTheta)
{
  theta = myTheta;
}
