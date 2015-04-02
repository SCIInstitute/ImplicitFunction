#include<cstdio>
#include"vec3.h"
#include<vector>
#include"fault.h"


Fault::Fault(ScatteredData *myData, RBF *myRBF)
{
	data = myData;
	rbf = myRBF;
}

double Fault::computeSide(vec3 x)
{
	return rbf->computeValue(x);
}
double Fault::computeValue(vec3 x)
{
	return computeSide(x);
}

void Fault::move(vec3 x)
{
}

void Fault::computeRBF()
{
	rbf->computeFunction();
}
