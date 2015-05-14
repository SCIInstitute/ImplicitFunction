#include <string>
#include <algorithm>
#include <vector>
#include <cmath>

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"
#include "RBF.h"
#include "Surface.h"
#include "fileIO.h"
#include "API.h"

using std::string;

void API::CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	mySurfaceData = new ScatteredData();
	readSurfaceDataFile(filename, mySurfaceData);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	//mySurfaceRBF->setDataReduction(Random);
	myKernel = ThinPlate;
	mySurface = new Surface(mySurfaceData, mySurfaceRBF);

	//Construct RBFs
	mySurface->computeRBF();

	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);
	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);
	value.resize((int)(mySampling[0]));
	for(int i=0; i<mySampling[0]; i++)
	{
		//if(i%10==0)
			//printf("%d/100 done\n", i); fflush(stdout);
		value[i].resize((int)(mySampling[1]));
		for(int j=0; j<mySampling[1]; j++)
		{
			//if(j%10==0)
			//	printf("\t%d/100 done\n", j); fflush(stdout);
			value[i][j].resize((int)(mySampling[2]));
			for(int k=0; k<mySampling[2]; k++)
			{
				//if(k%10==0)
				//	printf("\t\t%d/100 done\n", k); fflush(stdout);
				vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
				//std::cout<<"Computing Val ... "<<std::endl;
				double myVal = mySurface->computeValue(location);
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
				value[i][j][k]=myVal;
			}
		}
	}
}

API::API(string filename, string dimensions)
{
	//read from the dimentsion file here	TODO
	vec3 myOrigin(-30, -50, 80);
	vec3 mySize(60, 50, 10);
	vec3 mySampling(100, 100, 100);
	CreateSurface(filename, myOrigin, mySize, mySampling);
}

API::API()
{
}

API::API(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	CreateSurface(myData, myOrigin, mySize, mySampling);
}

vector<vector<vector<double> > >API::CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	vector<double> a,b,c,d;
	for(int i=0; i<myData.size(); i++)
	{
		a.push_back(myData[i][0]);
		b.push_back(myData[i][1]);
		c.push_back(myData[i][2]);
		d.push_back(0);
	}

	vector<double>::iterator minx = std::min_element(a.begin(), a.end());
	vector<double>::iterator miny = std::min_element(b.begin(), b.end());
	vector<double>::iterator minz = std::min_element(c.begin(), c.end());
	vector<double>::iterator maxx = std::max_element(a.begin(), a.end());
	vector<double>::iterator maxy = std::max_element(b.begin(), b.end());
	vector<double>::iterator maxz = std::max_element(c.begin(), c.end());
	vec3 myMin(*minx, *miny, *minz), myMax(*maxx, *maxy, *maxz);
	myMin = myMin - 0.05*mySize;
	myMax = myMax + 0.05*mySize;

	mySurfaceData = new ScatteredData(a,b,c,d);
	augmentNormalData(mySurfaceData);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	//mySurfaceRBF->setDataReduction(Random);
	myKernel = ThinPlate;
	mySurface = new Surface(mySurfaceData, mySurfaceRBF);

	//Construct RBFs
	mySurface->computeRBF();

	//sanity check
	for(int i=0; i<mySurfaceData->fnc.size(); i++)
	{
		vec3 myLocation(mySurfaceData->x[0][i], mySurfaceData->x[1][i], mySurfaceData->x[2][i]);
		double myVal = mySurface->computeValue(myLocation);
		double error  = fabs(myVal - mySurfaceData->fnc[i]);
		if (error>1e-3)
		{
			printf("%lf\n", error);
			fflush(stdout);
		}
	}

	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);
	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);

	value.resize((int)(mySampling[0]));
	for(int i=0; i<mySampling[0]; i++)
	{
		value[i].resize((int)(mySampling[1]));
		for(int j=0; j<mySampling[1]; j++)
		{
			value[i][j].resize((int)(mySampling[2]), -100);
		}
	}

	for(int i=0; i<mySampling[0]; i++)
	{
		vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX;
		if (location[0]<myMin[0]||location[0]>myMax[0])
			continue;
		for(int j=0; j<mySampling[1]; j++)
		{
			location = myOrigin + mySpacing[1]*j*vec3::unitY;
			if (location[1]<myMin[1]||location[1]>myMax[1])
				continue;
			for(int k=0; k<mySampling[2]; k++)
			{
				location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
				if (location[2]<myMin[2]||location[2]>myMax[2])
					continue;
				//std::cout<<"Computing Val ... "<<std::endl;
				double myVal = mySurface->computeValue(location);
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
				value[i][j][k]=myVal;
			}
		}
	}
	
	return value;
}

vec3 API::findNormal(ScatteredData *data, int n)
{
	int tot = data->x[0].size();
	int prev = (n-1)>=0?n-1:tot-1;
	int next = (n+1)<tot?n+1:0;

	while(data->x[2][prev]!=data->x[2][n])
	{
		prev = (prev-1)>=0?prev-1:tot-1;
	}
	
	while(data->x[2][next]!=data->x[2][n])
	{
		next = (next+1)<tot?next+1:0;
	}
	//printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);

	vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
	vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
	vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
	vec3 one = b-a;
	vec3 two = c-a;
	vec3 ret = one+two;
	return ret;
}

void API::augmentNormalData(ScatteredData *data)
{
	int n = data->x[0].size();
	for(int i=0; i<n; i++)
	{
		vec3 myNormal = findNormal(data, i);
		myNormal = normalize(myNormal);
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] + myNormal[j]);
		}
		data->fnc.push_back(10);

		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] - myNormal[j]);
		}
		data->fnc.push_back(-10);
	}
}
