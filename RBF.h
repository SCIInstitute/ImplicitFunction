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
#include <memory>

enum Kernel { Gaussian, ThinPlate, MultiQuadratic };
enum Acceleration { None, FastMultipole };
enum DataReduction { All, Random };

class RBF
{
public:
	RBF(const ScatteredData *myData, Kernel myKernel);

	void setAcceleration(Acceleration myAcceleration);
	void setDataReduction(DataReduction myDataReduction);

	void computeFunction();
	double computeValue(const vec3& x);

private:
	const ScatteredData* completeData_;
	std::unique_ptr<ScatteredData> data_;
	Kernel kernel_;
	std::vector<double> coeff_;
	Acceleration acceleration_;
	DataReduction dataReduction_;

	std::unique_ptr<FMM> fmm_;

	void computeFunctionForData(); // throws std::runtime_error
	void computeErrorForData(std::vector<std::pair<double, int> > &error);

	double computeKernel(int i, int j);
	double computeKernel(int i, const vec3& b);
	double computeRadialFunction(double r);

	void fmmBuildTree();
	void fmmPrintTree(BHNode* myNode, int stack);
	void fmmBuildTree(std::vector<int> &myPoints, BHNode *myNode);
	double fmmComputeValue(const vec3& x);
	double fmmComputeValueRecurse(const vec3& x, BHNode *myNode);
	double fmmComputeKernel(const vec3& b, BHNode *myNode);

  static const double EPSILON;
};

#endif //_RBF_H_
