#include<cstdio>
#include"vec3.h"
#include<vector>
#include"Surface.h"
#include"tear.h"
#include"fault.h"

Tear::Tear()
{
}

Tear::Tear(Surface *myPrimary, Surface *mySecondary, Fault *myClip)
{
	primary = myPrimary;
	secondary = mySecondary;
	clip = myClip;
}

double Tear::computeSide(vec3 x)
{
	double a = clip->computeValue(x);
	double b = fault->computeValue(x);

	double weight;
	if(a == 0 || (b/a) > 1)
		weight = 1;
	else 
		weight = 0.5*(b/a) + 0.5;
	
	double ret = weight * primary->computeValue(x) + (1-weight)*secondary->computeValue(x);
	return ret;
}

double Tear::computeValue(vec3 x)
{
	return computeSide(x);
}

void Tear::move(vec3 x)
{
}

