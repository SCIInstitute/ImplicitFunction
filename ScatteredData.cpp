#include <vector>
#include <cstdio>
#include <algorithm>
#include "ScatteredData.h"
#include "vec3.h"
#include "ETSP.h"
#include <cmath>

using std::vector;

//int ScatteredData::myAxis = 2;
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
	//for(int i=0; i<a.size(); i++)
	//{
	//	printf("%d %d\n", x[0][i], a[i]); fflush(stdout);
	//}
	origSize = a.size();
}

void ScatteredData::computeOrdering()
{
	SDmultisort();
	//printf("Sorted\n");
	//for(int i=0; i<myData.size(); i++)
	//	printf("%lf %lf %lf\n", myData[i][0], myData[i][0], myData[i][2]);

	int count = 0;
	//printf("%d %d\n", myData.size(), x[0].size()); fflush(stdout);
	while(count != x[0].size())
	{
		int start = count;
		//printf("%lf %d\n", myData[count][myAxis], x[0].size());
		while(count != x[0].size())
		{
			//printf("%lf %lf %lf\n", myData[count][myAxis], myData[start][myAxis],fabs(myData[count][myAxis]-myData[start][myAxis]));
			count++;
			int sortingAxis = axisInformation[start];
			if(axisInformation[count]!=axisInformation[start])
				break;
			if(fabs(myData[count][sortingAxis]-myData[start][sortingAxis]) > 1e-6)
				break;
		}
		int end = count-1;
		//printf("Start: %d End %d\n", start, end);
		ETSP reorder(myData, start, end, axisInformation[start]);
		for(int i=start; i<=end; i++)
		{
			myData[i][0] = x[0][i]=reorder.data[reorder.order[i-start]][0];
			myData[i][1] = x[1][i]=reorder.data[reorder.order[i-start]][1];
			myData[i][2] = x[2][i]=reorder.data[reorder.order[i-start]][2];
		}
	}
	//for(int i=0; i<myData.size(); i++)
	//	printf("%lf %lf %lf\n", myData[i][0], myData[i][0], myData[i][2]);
}

void ScatteredData::SDsort()
{
	vec3Sorter sortObject;
	sortObject.axisToSort=axisInformation[0];
	for(int i=0; i<x[0].size(); i++)
	{
		vec3 a(x[0][i],x[1][i],x[2][i]);
		myData.push_back(a);
	}
	std::sort(myData.begin(), myData.end(), sortObject);
}

void ScatteredData::SDmultisort()
{
	vec3Sorter sortObject;
	std::vector<vec3> dataX, dataY, dataZ;
	for(int i=0; i<x[0].size(); i++)
	{
		vec3 a(x[0][i],x[1][i],x[2][i]);
		switch(axisInformation[i])
		{
			case 0:
				dataX.push_back(a);
				break;
			case 1:
				dataY.push_back(a);
				break;
			case 2:
				dataZ.push_back(a);
				break;
		}
		
	}
	sortObject.axisToSort=0;
	std::sort(dataX.begin(), dataX.end(), sortObject);
	sortObject.axisToSort=1;
	std::sort(dataY.begin(), dataY.end(), sortObject);
	sortObject.axisToSort=2;
	std::sort(dataZ.begin(), dataZ.end(), sortObject);
	axisInformation.clear();
	for(int i=0; i<dataX.size(); i++)
	{
		myData.push_back(dataX[i]);
		axisInformation.push_back(0);
	}
	for(int i=0; i<dataY.size(); i++)
	{
		myData.push_back(dataY[i]);
		axisInformation.push_back(1);
	}
	for(int i=0; i<dataZ.size(); i++)
	{
		myData.push_back(dataZ[i]);
		axisInformation.push_back(2);
	}
}

