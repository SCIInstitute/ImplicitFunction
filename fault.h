#ifndef _FAULT_H_
#define _FAULT_H_

#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"

#include <vector>

using std::vector;

class Fault
{
	public:
	ScatteredData *data;
	RBF *rbf;
	Fault(ScatteredData *myData, RBF *myrbf);

	double computeSide(vec3 x);
	double computeValue(vec3 x);
	void move(vec3 x);
	void computeRBF();
	
};


#endif //_FAULT_H_
