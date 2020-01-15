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

#ifndef _SCATTEREDDATA_H_
#define _SCATTEREDDATA_H_

// STL Includes
#include <vector>
#include "vec3.h"

enum axis_t { X = 0, Y, Z };

class ScatteredData
{
public:
  // TODO: what happens when vector sizes don't match???
  ScatteredData();
  ScatteredData(const std::vector<double>& points_x,
                const std::vector<double>& points_y,
                const std::vector<double>& points_z,
                const std::vector<double>& func,
                const std::vector<axis_t>& axisInfo);

  void setData(const std::vector<double>& points_x,
               const std::vector<double>& points_y,
               const std::vector<double>& points_z,
               const std::vector<double>& func);

  // TODO: make private
  // leftovers_ is set of points inside convex hull
  //
  // surfacePoints_ is array of vectors of point x components, point y components and point z components
  // If convex hull is computed, then surfacePoints_ is points on convex hull boundary
  //
  // fnc_ is vector of single threshold value and constant values (positive, negative) at indices
  // where normal data has been pushed to surfacePoints_ vector
  //
  // TODO: make surfacePoints_ and leftovers_ vector of vec3??
  // cohesive data structure would be better...
  std::vector<double> leftovers_[3], fnc_;
  std::vector<axis_t> axisInformation_, updatedAxisInformation_;

	void compute2DHull();
	int origSize_;

  vec3 surfacePoint(size_t i) const;

private:
  void SDsort();
	void SDmultisort();
public: //TODO: private
  std::vector<double> surfacePoints_[3];
private:
  std::vector<vec3> inputData_, convexHullData_;
};

struct vec3Sorter
{
  axis_t axisToSort;
  int myAxis;
  bool operator() (const vec3& a, const vec3& b)
  {
    if (axisToSort == X) myAxis = 0;
    if (axisToSort == Y) myAxis = 1;
    if (axisToSort == Z) myAxis = 2;
    return (a[myAxis] < b[myAxis]);
  }
};

#endif //_SCATTEREDDATA_H_
