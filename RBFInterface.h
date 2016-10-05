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

#include "ScatteredData.h"
#include "vec3.h"
#include "RBF.h"

typedef std::vector< std::vector< std::vector<double> > > DataStorage;

class RBFInterface
{
public:
	RBFInterface(std::vector<vec3> myData,
               vec3 myOrigin, vec3 mySize, vec3 mySpacing,
               double myOffset, std::vector<axis_t> myAxis, const bool use2DConvexHull=false, Kernel kernel=ThinPlate);

//  RBFInterface(std::vector<vec3> myData, 
//               vec3 myOrigin, vec3 mySize, vec3 mySpacing,
//               double myOffset, std::vector<axis_t> myAxis, bool use2DConvexHull=false, Kernel kernel=ThinPlate);

	double getThresholdValue() const { return thresholdValue_; }
  const DataStorage getRasterData() const { return rasterData_; }

private:
  void CreateSurface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing, double myOffset, std::vector<axis_t> myAxis);

  void augmentNormalData(const double myOffset);
	vec3 findNormalAxis(const int n);

  double distanceToLineSegment(const vec3& a, const vec3& b, const vec3& p);

  ScatteredData *surfaceData_;
  RBF *rbf_;
  DataStorage rasterData_;

  const double thresholdValue_;
//  double offset_;
  const bool use2DConvexHull_;
	Kernel kernel_;
  std::vector<double> points_x_, points_y_ , points_z_, threshold_;

  // change to inside or outside
  std::vector<vec3> inNormals, outNormals;

  static const double EPSILON;
  static const double SMALL_EPSILON;
};

#endif //_RBFInterface_H_ 
