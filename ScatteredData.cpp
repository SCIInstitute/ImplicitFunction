#include <vector>
#include <cstdio>
#include <algorithm>
#include "ScatteredData.h"
#include "vec3.h"
#include "ETSP.h"
#include "2DConvexHull.h"
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
	vec3 temp(0,0,0);
	for(int i=0; i<origSize; i++)
	{
		for(int j=0; j<3; j++)
		{
			temp[j]+=x[j][i];
		}
	}
	for(int j=0; j<3; j++)
		centroid[j]=temp[j]/origSize;
}

void ScatteredData::computeOrdering()
{
	SDmultisort();
	//printf("Sorted\n");
	//for(int i=0; i<myData.size(); i++)
	//	printf("%d %lf %lf %lf\n", i, myData[i][0], myData[i][0], myData[i][2]);

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
			axis_t sortingAxis = axisInformation[start];
			if(axisInformation[count]!=axisInformation[start])
				break;
			if(fabs(myData[count][sortingAxis]-myData[start][sortingAxis]) > 1e-6)
				break;
		}
		int end = count-1;
		printf("BCD Start: %d End %d\n", start, end);
		int myAxis;
		switch(axisInformation[start])
		{
			case X: myAxis = 0; break;
			case Y: myAxis = 1; break;
			case Z: myAxis = 2; break;

		}
		ETSP reorder(myData, start, end, myAxis);
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
			case X:
				dataX.push_back(a);
				break;
			case Y:
				dataY.push_back(a);
				break;
			case Z:
				dataZ.push_back(a);
				break;
		}
		
	}
        myData.clear();
	sortObject.axisToSort=X;
	std::sort(dataX.begin(), dataX.end(), sortObject);
	sortObject.axisToSort=Y;
	std::sort(dataY.begin(), dataY.end(), sortObject);
	sortObject.axisToSort=Z;
	std::sort(dataZ.begin(), dataZ.end(), sortObject);
	axisInformation.clear();
	for(int i=0; i<dataX.size(); i++)
	{
		myData.push_back(dataX[i]);
		axisInformation.push_back(X);
	}
	for(int i=0; i<dataY.size(); i++)
	{
		myData.push_back(dataY[i]);
		axisInformation.push_back(Y);
	}
	for(int i=0; i<dataZ.size(); i++)
	{
		myData.push_back(dataZ[i]);
		axisInformation.push_back(Z);
	}
}


void ScatteredData::compute2DHull()
{
   
        std::vector<int> reorderedData;
	for(int i=0; i<x[0].size(); i++)
		printf("%d %lf %lf %lf\n", i, x[0][i], x[1][i], x[2][i]);
	SDmultisort();
	printf("Sorted\n");
	for(int i=0; i<myData.size(); i++)
		printf("%d %lf %lf %lf\n", i, myData[i][0], myData[i][1], myData[i][2]);

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
			axis_t sortingAxis = axisInformation[start];
			if(axisInformation[count]!=axisInformation[start])
				break;
			if(fabs(myData[count][sortingAxis]-myData[start][sortingAxis]) > 1e-6)
				break;
		}
		int end = count-1;
		printf("ABC Start: %d End %d\n", start, end);
		int myAxis;
		switch(axisInformation[start])
		{
			case X: myAxis = 0; break;
			case Y: myAxis = 1; break;
			case Z: myAxis = 2; break;

		}
                vector<vec3> inPoints;
                inPoints.resize(end-start+1);
                for(int i=start; i<=end; i++)
                {
                  inPoints[i-start][0] = myData[i][0];
                  inPoints[i-start][1] = myData[i][1];
                  inPoints[i-start][2] = myData[i][2];
                }
                //printf("inPoints written\n");
		vector<int> reorder = getConvexHull(inPoints, myAxis);
                printf("Convex Hull found\n");
		for(int i=0; i<reorder.size(); i++)
		{
                        printf("%d ", i, reorder[i]+start);
			reorderedData.push_back(reorder[i]+start);
		}
                printf("\n");
	}
	//rewrite everything
	std::vector<double> newx[3], newfnc;
	std::vector<vec3> myNewData;
	std::vector<axis_t> newAxisInformation;
        //printf("Rewriting everything\n");
        for(int i=0; i<reorderedData.size(); i++)
        {
          int j = reorderedData[i];
          newx[0].push_back(myData[j][0]);
          newx[1].push_back(myData[j][1]);
          newx[2].push_back(myData[j][2]);
          newfnc.push_back(fnc[j]);
          myNewData.push_back(myData[j]);
          newAxisInformation.push_back(axisInformation[j]);
        }
        x[0].clear();
        x[1].clear();
        x[2].clear();
        x[0] = newx[0];
        x[1] = newx[1];
        x[2] = newx[2];
        fnc.clear();
        fnc = newfnc;
        myData.clear();
        myData=myNewData;
        axisInformation.clear();
        axisInformation = newAxisInformation;
	//for(int i=0; i<x[0].size(); i++)
	//	printf("%d %lf %lf %lf\n", i, x[0][i], x[1][i], x[2][i]);
        printf("Convex hull points\n");
	for(int i=0; i<myData.size(); i++)
		printf("%d %lf %lf %lf\n", i, myData[i][0], myData[i][1], myData[i][2]);
        origSize = myData.size();
}
