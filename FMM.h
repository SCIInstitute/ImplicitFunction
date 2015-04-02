#ifndef _FMM_H_
#define _FMM_H_

#include "ScatteredData.h"
#include "SparseMatrix.h"
#include "vec3.h"
#include "BHNode.h"

// STL Includes
#include <vector>
using std::vector;

class FMM
{
	public:
	FMM();
	~FMM();

	//variables
	BHNode *tree;
	double theta;
	int numOfNodes;

	vector<BHNode*> nodePointer;

	//functions to set the variable
	void setTheta(double myTheta);



	void freeTheTree(BHNode *myNode);
};

#endif //_FMM_H_

