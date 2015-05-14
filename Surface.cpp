#include<cstdio>
#include"vec3.h"
#include<vector>
#include"Surface.h"

Surface::Surface()
{
}

Surface::Surface(ScatteredData *myData, RBF *myRBF)
{
	data = myData;
	rbf = myRBF;
}

double Surface::computeSide(vec3 x)
{
	return rbf->computeValue(x);
}

double Surface::computeValue(vec3 x)
{
	return computeSide(x);
}

void Surface::move(vec3 x)
{
}

void Surface::computeRBF()
{
	rbf->computeFunction();
}
