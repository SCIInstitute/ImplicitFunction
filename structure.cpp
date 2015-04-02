#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"
#include "horizon.h"
#include "structure.h"
#include "tree.h"

#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <utility>

using std::vector;
using std::string;
using std::cout;
using std::endl;

Structure::Structure(vector<Fault*> &myFault, vector<Horizon*> &myHorizon)
{
	for(int i=0; i<myFault.size(); i++)
		fault.push_back(myFault[i]);
	for(int i=0; i<myHorizon.size(); i++)
		horizon.push_back(myHorizon[i]);
}

void Structure::constructTree(string filename)
{

	vector<int> nodeToFeature;
	vector<vector<int> > children;
	vector<int> whichFeature;

	readStructureFile(filename, nodeToFeature, children, whichFeature);	
	int n = children.size();
	
	nodes.resize(n);
	for(int i=0; i<n; i++)
	{
		nodes[i] = new StructureTree();
		nodes[i]->index = i;
		//cout<<nodeToFeature[i]<<" "<<whichFeature[i]<<endl;
	}
	tree = nodes[0];

	for(int i=0; i<n; i++)
	{
		//cout<<i<<" "<<whichFeature[i]<<" ";
		//if(children[i].size()>0)
		//	cout<<children[i][0]<<" "<<children[i][1]<<" ";
		//cout<<nodeToFeature[i]<<endl;
		//cout<<endl;
		switch(whichFeature[i])
		{
			case 0:
				nodes[i]->fault = fault[nodeToFeature[i]];
				nodes[i]->horizon = NULL;
				break;
			case 1:
				nodes[i]->horizon = horizon[nodeToFeature[i]];
				nodes[i]->fault = NULL;
				break;
			case 2:
				nodes[i]->material = nodeToFeature[i];
		}
		if(children[i].size()>0)
		{
			nodes[i]->left  = nodes[children[i][0]];
			nodes[i]->right = nodes[children[i][1]];
		}
		else
		{
			nodes[i]->left  = NULL;
			nodes[i]->right = NULL;
		}
	}
}

std::pair<int,double> Structure::returnHorizon(vec3 x)
{
	//traverse tree and return horizon
	double myVal;
	StructureTree *currNode, *prevNode;
	currNode = tree;

	while(currNode!=NULL)
	{
		prevNode = currNode;
		double myVal = currNode->computeValue(x);
		if(myVal>=0)
			currNode = currNode->left;
		else
			currNode = currNode->right;
	}
	int myMat = prevNode->material;

	double myDist=1e10;
	int index = prevNode->index;
	int n = neighbors[index].size();
	for(int i=0; i<n; i++)
	{
		double tempVal = nodes[neighbors[index][i]]->computeValue(x);
		myDist = std::min(myDist,fabs(tempVal));
	}
	return(std::make_pair(myMat,myDist));
	
}

void Structure::computeRBFs()
{
	for (int i=0; i<fault.size(); i++)
		fault[i]->computeRBF();
	for (int i=0; i<horizon.size(); i++)
		horizon[i]->computeRBF();
}

void Structure::readStructureFile(string filename, vector<int> &nodeToFeature, vector<vector<int> > &children, vector<int> &whichFeature)
{
	std::cout <<"Reading file '" << filename <<"'"<<std::endl;
	std::ifstream datafile(filename.c_str());
	if(datafile.is_open())
	{
		std::cout<<"Opened"<<std::endl;
		//number of nodes in the tree
		int n;
		datafile>>n;
		nodeToFeature.resize(n);
		children.resize(n);
		for(int i=0; i<n; i++)
			children[i].clear();
		whichFeature.resize(n);

		//node -> feature correspondace as well as informing us 
		//which among those is a fault(0), horizon(1), or a surface(2). 
		int isHorizon;
		for(int i=0; i<n; i++)
		{
			datafile>>nodeToFeature[i];
		}
		for(int i=0; i<n; i++)
		{
			datafile>>whichFeature[i];
		}

		//parent-child relationship
		int m;
		datafile>>m;
		for(int i=0; i<m; i++)
		{
			int j;
			datafile>>j;
			children[j].resize(2);
			datafile>>children[j][0]>>children[j][1];
		}

		//neighborinformation
		neighbors.resize(n);
		for(int i=0; i<n; i++)
		{
			int m;
			datafile>>m;
			//std::cout<<i<<" ";
			for(int j=0; j<m; j++)
			{
				int p;
				datafile>>p;
				neighbors[i].push_back(p);
				//std::cout<<p<<" ";
			}
			//std::cout<<std::endl;
		}
		cout<<"Done"<<endl;

		
	}
	else
	{
		std::cout <<"Could not open '" << filename <<"'"<<std::endl;
	}
}

