#ifndef _TEAR_H_
#define _TEAR_H_

#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"
#include "Surface.h"
#include "fault.h"

#include <vector>

using std::vector;

class Tear
{
	public:
	Surface *primary, *secondary;
	Fault *clip, *fault;

	Tear();
	Tear(Surface *myPrimary, Surface *mySecondary, Fault *myClip);
	double computeSide(vec3 x);
	double computeValue(vec3 x);
	void move(vec3 x);
};

#endif //_TEAR_H_ 
