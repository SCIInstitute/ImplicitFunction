#ifndef _SCATTEREDDATA_H_
#define _SCATTEREDDATA_H_

// STL Includes
#include <vector>
using std::vector;

class ScatteredData
{
	public:
	ScatteredData();
	ScatteredData(vector<double> a, vector<double> b, vector<double> c, vector<double> d);
	void setData(vector<double> a, vector<double> b, vector<double> c, vector<double> d);
	
	vector<double> x[3], fnc;
};

#endif //_SCATTEREDDATA_H_
