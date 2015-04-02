#ifndef _SAMPLEDATA_H_
#define _SAMPLEDATA_H_

// STL Includes
#include <vector>

#include "ScatteredData.h"
#include "vec3.h"

using std::vector;

void readScatteredData(char file[], ScatteredData *myData);
void readFunctionData(ScatteredData *myData);
void sphereData(ScatteredData *myData, vec3 center, double r, double theta1, double theta2, double phi1, double phi2, int number);


#endif //_SAMPLEDATA
