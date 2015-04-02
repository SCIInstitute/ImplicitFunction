#ifndef _API_H_
#define _API_H_

#include <vector>
#include <string>

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"
#include "RBF.h"
#include "horizon.h"

using std::string; 

class API
{
	public:
	ScatteredData *mySurfaceData;
	Horizon *mySurface;
	RBF *mySurfaceRBF;
	Kernel myKernel;
	vector<vector<vector<double> > > value;

	void CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	vector<vector<vector<double> > >CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	API(string filename, string dimensions);
	API(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	API();

	private:

	void augmentNormalData(ScatteredData *data);
	vec3 findNormal(ScatteredData *data, int n);
};

#endif //_API_H_ 
