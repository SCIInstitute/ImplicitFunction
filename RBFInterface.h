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

#ifndef _RBFInterface_H_
#define _RBFInterface_H_

#include <vector>
#include <string>

#include "ScatteredData.h"
#include "vec3.h"
#include "RBF.h"
#include "Surface.h"

//using std::string; 
typedef vector< vector< vector<double> > > DataStructure;

class RBFInterface
{
public:
	ScatteredData *mySurfaceData;
	Surface *mySurface;
	RBF *mySurfaceRBF;
	Kernel myKernel;
	DataStructure value;

	void CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing, double myOffset);
	RBFInterface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing, double myOffset);

	double getThresholdValue() const { return thresholdValue; }
	void setOffset(double myOffset);

private:
	void augmentNormalData(ScatteredData *data, double myOffset);
	vec3 findNormal(ScatteredData *data, int n);

  double thresholdValue;
  double offset;

  static const double EPSILON;
};

#endif //_RBFInterface_H_ 
