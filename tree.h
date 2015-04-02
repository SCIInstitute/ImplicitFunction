#ifndef _STree_H_
#define _STree_H_

#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"
#include "horizon.h"
#include "fault.h"
#include "tear.h"

#include <vector>

using std::vector;

class StructureTree
{
	public:
	StructureTree *left, *right;
	bool isFault;
	Fault *fault;
	Horizon *horizon;
	Tear *tear;
	int material;
	int index;

	double computeValue(vec3 x);
};


#endif //_STree_H_
