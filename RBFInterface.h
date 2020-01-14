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
typedef std::vector<size_t> IndexList;
typedef std::vector<axis_t> AxisList;

class RBFInterface
{
public:
	RBFInterface(const std::vector<vec3>& myData,
               const vec3& myOrigin, const vec3& mySize, const vec3& mySpacing,
               const double myOffset, AxisList myAxis,
               const bool compute2DConvexHull=true,
               const bool invertSeedOrder=false, Kernel kernel=ThinPlate);

	double getThresholdValue() const { return thresholdValue_; }
  const DataStorage getRasterData() const { return rasterData_; }
  const ScatteredData* getSurfaceData() const { return surfaceData_.get(); }

private:
  // 2D calculation
  void create2DSurface(); // propagates exceptions
  void create3DSurface();
  void augmentNormalData();
  vec3 findNormalAxis(const int n);
  void createRasterizedSurface();

  std::unique_ptr<ScatteredData> surfaceData_;
  DataStorage rasterData_;

  const double thresholdValue_;
  const vec3 origin_;
  const vec3 size_;
  const vec3 spacing_;
  const double offset_;
  AxisList axisList_;
  const bool compute2DConvexHull_;
  const bool invertSeedOrder_;
  const Kernel kernel_;
  std::vector<double> points_x_, points_y_ , points_z_, threshold_;

  // change to inside or outside
  std::vector<vec3> inNormals, outNormals;

  static const double SMALL_EPSILON;
};

#endif //_RBFInterface_H_
