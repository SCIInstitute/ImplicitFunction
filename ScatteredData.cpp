#include <vector>
#include <cstdio>
#include "ScatteredData.h"

using std::vector;

ScatteredData::ScatteredData()
{
}

ScatteredData::ScatteredData(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
	setData(a,b,c,d);
}

void ScatteredData::setData(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
	x[0] = a; x[1] = b; x[2] = c; fnc = d;
	/*for(int i=0; i<a.size(); i++)
	{
		printf("%d %d\n", x[0][i], a[i]); fflush(stdout);
	}*/
}

