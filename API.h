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

#ifndef _API_H_
#define _API_H_

#include <vector>
#include <string>
using std::vector;
using std::string;

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"
#include "RBF.h"
#include "Surface.h"

class API
{
	public:
	ScatteredData *mySurfaceData;
	Surface *mySurface;
	RBF *mySurfaceRBF;
	Kernel myKernel;
	std::vector<std::vector<std::vector<double> > > value;

	void CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	std::vector<std::vector<std::vector<double> > >CreateSurface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	API(string filename, string dimensions);
	API(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	API();

	private:

	void augmentNormalData(ScatteredData *data);
	vec3 findNormal(ScatteredData *data, int n);
	vec3 findSphericalNormal(ScatteredData *data, int n);
	vec3 findNormalAxis(ScatteredData *data, int n);
};

#endif //_API_H_ 
