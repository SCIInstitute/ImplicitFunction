// STL Includes
#include <vector>

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"

#include <cstdlib>
#include <cmath>
#include <cstdio>

void readScatteredData(char file[], ScatteredData *myData)
{
	myData = NULL;
}

void readFunctionData(ScatteredData *myData)
{

	vec3 center(50, 50, 50);
	printf("Collecting Data ... "); fflush(stdout);
	sphereData(myData, center, 40, -0.25, 0.25, -0.25, +0.25, 100);
	printf("Done\n"); fflush(stdout);
}

void sphereData(ScatteredData *myData, vec3 center, double r, double theta1, double theta2, double phi1, double phi2, int number)
{
	std::vector<double> a,b,c,d;
	
	for(int i=0; i<3*number; i++)
	{
		int myRand[2];
		myRand[0] = rand()%101;
		myRand[1] = rand()%101;

		double theta = theta1 + (theta2 - theta1)*(myRand[0]/100.0);
		double phi   =   phi1 + (  phi2 -   phi1)*(myRand[1]/100.0);
		double depth;

		switch(i%3)
		{
			case 0:
				depth=-1;				
				break;
			case 1:
				depth=0;				
				break;
			default:
				depth=1;				
		}

		a.push_back(center[0] + (r + depth)*sin(theta)*cos(phi));
		b.push_back(center[1] + (r + depth)*sin(theta)*sin(phi));
		c.push_back(center[2] + (r + depth)*cos(theta));
		d.push_back(depth*1);
		printf("%d %lf %lf %lf %lf\n", i, a[i], b[i], c[i], d[i]);
	}
	
	myData->setData(a, b, c, d);
}

