#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"
#include "SparseMatrix.h"
#include "LinearSolver.h"
#include "RBF.h"
#include "fault.h"
#include "horizon.h"
#include "structure.h"
#include "API.h"


using std::vector;

void constructRBFforSampleData()
{
	ScatteredData *myData;
	vector<vector<vector<double> > > value;
	value.resize(100);	

	myData = new ScatteredData();

	readFunctionData(myData);
	Kernel myKernel = ThinPlate;
	RBF myRBF(myData, myKernel);
	myRBF.setAcceleration(FastMultipole);
	myRBF.setDataReduction(Random);
	myRBF.computeFunction();
	for(int i=0; i<100; i++)
	{
		if(i%10==0)
			printf("%d/100 done\n", i); fflush(stdout);
		value[i].resize(100);
		for(int j=0; j<100; j++)
		{
			value[i][j].resize(100);
			for(int k=0; k<100; k++)
			{
				vec3 location = i*vec3::unitX + j*vec3::unitY + k*vec3::unitZ;
				double val = myRBF.computeValue(location);
				value[i][j][k]=val;
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], val); fflush(stdout);
			}
		}
	}

	vec3 dim(100,100,100);
	vec3 spacing(1,1,1);
	string filename = "sphere.nrrd";
	//writeNrrdFile(filename, value, dim, spacing);
}

void aSimpleCase()
{
	//Take the input
	ScatteredData *mySurfaceData;
	mySurfaceData = new ScatteredData();

	//readSurfaceDataFile("./segmentation/sample_points_DEMRI.txt", mySurfaceData);
	vector<Horizon*> mySurface;

	RBF *mySurfaceRBF;
	Kernel myKernel = ThinPlate;

	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	mySurfaceRBF->setDataReduction(Random);

	mySurface.push_back(new Horizon(mySurfaceData, mySurfaceRBF));


	//Construct RBFs
	mySurface[0]->computeRBF();

	
	//Compute material
	vector<vector<vector<double> > > value;
	value.resize(100);
	vec3 myOrigin(-30, -50, 80);
	for(int i=0; i<100; i++)
	{
		if(i%10==0)
			printf("%d/100 done\n", i); fflush(stdout);
		value[i].resize(100);
		for(int j=0; j<100; j++)
		{
			//if(j%10==0)
			//	printf("\t%d/100 done\n", j); fflush(stdout);
			value[i][j].resize(100);
			for(int k=0; k<100; k++)
			{
				//if(k%10==0)
				//	printf("\t\t%d/100 done\n", k); fflush(stdout);
				vec3 location = myOrigin + 0.6*i*vec3::unitX + 0.5*j*vec3::unitY + 0.1*k*vec3::unitZ;
				//std::cout<<"Computing Val ... "<<std::endl;
				double myVal = mySurface[0]->computeValue(location);
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
				value[i][j][k]=myVal;
			}
		}
	}

	vec3 dim(100,100,100);
	vec3 spacing(0.6,0.5,0.1);
	string filename = "surface.nrrd";
	//writeNrrdFile(filename, value, dim, spacing);
}

void getAPI()
{
	API *myAPI;
	myAPI = new API("./segmentation/sample_points_DEMRI.txt", "dummy");
}
int main()
{
	//constructRBFforSampleData();
	//aSimpleCase();
	getAPI();
	//aComplicatedCase();
	//aSimpleCase2();
	//aComplicatedCase2();
	return 0;
}

