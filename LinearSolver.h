//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

// STL Includes
#include <vector>
#include "SparseMatrix.h"

class LinearSolver
{
public:
	LinearSolver();
	LinearSolver(SparseMatrix *myMat);
	void setMatrix(SparseMatrix *myMat);

	//Sparse Linear Solvers
	std::vector<double> biCGStab(std::vector<double> &b);
	void biCGStab(std::vector<double> &b, std::vector<double> &x);
private:
	SparseMatrix *A;
	double norm(std::vector<double> &a);
	double norm(std::vector<double> &a, std::vector<double> &b);
	void SpMV(std::vector<double> &a, std::vector<double> &b);
};

#endif //_LINEARSOLVER_H_

