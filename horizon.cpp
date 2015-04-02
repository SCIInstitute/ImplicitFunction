#include<cstdio>
#include"vec3.h"
#include<vector>
#include"horizon.h"

Horizon::Horizon()
{
}

Horizon::Horizon(ScatteredData *myData, RBF *myRBF)
{
	data = myData;
	rbf = myRBF;
}

double Horizon::computeSide(vec3 x)
{
	return rbf->computeValue(x);
}

double Horizon::computeValue(vec3 x)
{
	return computeSide(x);
}

void Horizon::move(vec3 x)
{
}

void Horizon::computeRBF()
{
	rbf->computeFunction();
}
