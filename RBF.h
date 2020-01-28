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

#ifndef _RBF_H_
#define _RBF_H_

#include "ScatteredData.h"
#include "SparseMatrix.h"
#include "vec3.h"
#include "FMM.h"

// STL Includes
#include <vector>
#include <utility>
#include <chrono>
#include <memory>

enum Kernel { Gaussian, ThinPlate, MultiQuadratic };
enum Acceleration { None, FastMultipole };
enum DataReduction { All, Random };

struct ScopedTimer
{
	explicit ScopedTimer(const char* message) : message_(message)
	{
		start = std::chrono::high_resolution_clock::now();
	}
	~ScopedTimer()
	{
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << message_ << " took " << std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count() << " ms." << std::endl;
	}
private:
	std::string message_;
	std::chrono::high_resolution_clock::time_point start;
};

class RBF
{
public:
	RBF(const ScatteredData *myData, Kernel myKernel);

	void setAcceleration(Acceleration myAcceleration);
	void setDataReduction(DataReduction myDataReduction);

	void computeFunction();
	double computeValue(const vec3& x) const;

private:
	const ScatteredData* completeData_;
	ScatteredData data_;
	Kernel kernel_;
	std::vector<double> coeff_;
	Acceleration acceleration_;
	DataReduction dataReduction_;

	std::unique_ptr<FMM> fmm_;

	void computeFunctionForData(); // throws std::runtime_error
	void computeErrorForData(std::vector<std::pair<double, int> > &error);

	double computeKernel(int i, int j) const;
	double computeKernel(int i, const vec3& b) const;
	double computeRadialFunctionOnSquaredDistance(double r2) const;
	double computeSumOfAllKernels(const vec3& b) const;

	mutable double* coeffPtr_;
	mutable size_t kernelBufferSize_;

	void fmmBuildTree();
	void fmmPrintTree(BHNode* myNode, int stack);
	void fmmBuildTree(std::vector<int> &myPoints, BHNode *myNode);
	double fmmComputeValue(const vec3& x) const;
	double fmmComputeValueRecurse(const vec3& x, BHNode *myNode) const;
	double fmmComputeKernel(const vec3& b, BHNode *myNode) const;

  static const double EPSILON;
};

#endif //_RBF_H_
