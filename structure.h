#ifndef _STRUCTURE_H_
#define _STRCUTURE_H_

#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"
#include "horizon.h"
#include "fault.h"
#include "tree.h"

#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>

using std::vector;
using std::string;

class Structure
{
	
	vector<Fault*> fault;
	vector<Horizon*> horizon;

	public:
	Structure(vector<Fault*> &myFault, vector<Horizon*> &myHorizon);

	StructureTree *tree;
	vector<StructureTree*> nodes;
	vector<vector<int> > neighbors;


	void constructTree(string filename);    //The input has to come from the user
	void readStructureFile(string filename, vector<int> &nodeToFeature, vector<vector<int> > &children, vector<int> &whichFeature);
	void computeRBFs();
	std::pair<int,double> returnHorizon(vec3 x);
};


#endif //_STRUCTURE_H_
