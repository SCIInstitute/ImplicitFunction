#ifndef _RBF_H_
#define _RBF_H_

#include "ScatteredData.h"
#include "SparseMatrix.h"
#include "vec3.h"
#include "FMM.h"

// STL Includes
#include <vector>
#include <utility>
using std::vector;
using std::pair;

enum Kernel {Guassian, ThinPlate, MultiQuadratic};
enum Acceleration {None, FastMultipole};
enum DataReduction {All, Random};
class RBF
{
	public:
	RBF();
	RBF(ScatteredData *myData, Kernel myKernel);

	void setKernel(Kernel myKernel);
	void setData(ScatteredData *myData);
	void setAcceleration(Acceleration myAcceleration);
	void setDataReduction(DataReduction myDataReduction);

	void computeFunction();
	double computeValue(vec3 x);
	private:
	ScatteredData *data, *completeData;
	Kernel kernel;
	vector<double> coeff;
	Acceleration acceleration;
	DataReduction dataReduction;

	FMM *fmm;

	void computeFunctionForData();
	void computeErrorForData(vector<pair<double, int> > &error);

	double computeKernel(int i, int j);
	double computeKernel(int i, vec3 b);
	double computeRadialFunction(double r);


	void fmmComputeFunction();
	void fmmBuildTree();
	void fmmPrintTree(BHNode* myNode, int stack);
	void fmmBuildTree(vector<int> &myPoints, BHNode *myNode);
	double fmmComputeValue(vec3 x);
	double fmmComputeValueRecurse(vec3 x, BHNode *myNode);
	double fmmComputeKernel(vec3 b, BHNode *myNode);
	double fmmComputeKernel(BHNode *myNode, vec3 b);
	
};

#endif //_RBF_H_
