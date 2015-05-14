#ifndef _SCATTEREDDATA_H_
#define _SCATTEREDDATA_H_

// STL Includes
#include <vector>
#include "vec3.h"

using std::vector;

class ScatteredData
{
	public:
	ScatteredData();
	ScatteredData(vector<double> a, vector<double> b, vector<double> c, vector<double> d);
	void setData(vector<double> a, vector<double> b, vector<double> c, vector<double> d);
	vector<double> x[3], fnc;
	vector<vec3> myData;
        void computeOrdering();
        static int myAxis;

	private:
        void SDsort();
};

struct vec3Sorter {
  bool operator() (vec3 a,vec3 b) 
  { 
    return (a[ScatteredData::myAxis]<b[ScatteredData::myAxis]);
  }
};

#endif //_SCATTEREDDATA_H_
