#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

// STL Includes
#include <vector>
#include "SparseMatrix.h"

using std::vector;

class LinearSolver
{
public:
	LinearSolver();
	LinearSolver(SparseMatrix *myMat);
	void setMatrix(SparseMatrix *myMat);

	//Sparse Linear Solvers
	vector<double> biCGStab(vector<double> &b);
	void biCGStab(vector<double> &b, vector<double> &x);
private:
	SparseMatrix *A;
	double norm(vector<double> &a);
	double norm(vector<double> &a, vector<double> &b);
	void SpMV(vector<double> &a, vector<double> &b);
};

#endif //_LINEARSOLVER_H_

