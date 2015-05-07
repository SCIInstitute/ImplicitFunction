#ifndef _ETSP_H_
#define _ETSP_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include "vec3.h"

using std::vector;

class ETSP
{
	public:
	ETSP(vector<vec3> &myData, int start, int end, int myplane);
	int plane;
	vector<vec3> data;
        vector<int> order;

	private:
	vector<vector<int> > graph;  //adjeceny list of the MST
        vector<int> inverse_mapped_order;
	void MST();
	vector<int> match();
	void orderFromMST();
	void traverse(int curr, vector<bool> &visited);
	void orderFromMatch();
	double areaTri(int i, int j);
	double computeArea();
};

#endif //_ETSP_H_
