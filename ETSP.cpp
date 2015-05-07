#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include "vec3.h"
#include "ETSP.h"
#include <utility>
#include <cstdio>

using std::vector;
using std::priority_queue;
using std::pair;
using std::make_pair;


ETSP::ETSP(vector<vec3> &myData, int start, int end, int myplane)
{
  for(int i=start; i<=end; i++)
  {
    vec3 temp(myData[i]);
    data.push_back(temp);
  }
  plane = myplane;
  //printf("Data copied\n");
  MST();
  //printf("MST Done\n");
  orderFromMST();
  //printf("Order\n");
  orderFromMatch();
  //printf("Reorder\n");
  //printf("Area: %lf\n", computeArea());
  if (computeArea()<0)
    std::reverse(order.begin(), order.end());
  //for(int i=0; i<order.size(); i++)
  //  printf("%d ",order[i]);
  //printf("\n");
}

double ETSP::computeArea()
{
  double area=0;
  for(int i=0; i<order.size(); i++)
  {
    area+=areaTri(i, (i+1)%order.size());
  }
  return area;
}

double ETSP::areaTri(int i, int j)
{
  vec3 a(data[order[i]]), b(data[order[j]]), n(0,0,0);
  n[plane] = 1;
  return dot(cross(a,b), n);
}

void ETSP::MST()
{
  priority_queue<pair<double, pair<int,int> > > myQueue;

  vector<bool> visited;
  visited.resize(data.size(),false);
  graph.resize(data.size());

  int visitedCount=1;
  visited[0]=true;

  for(int i=1; i<data.size(); i++)
  {
    double dist = length(data[0] - data[i]);
    myQueue.push(make_pair(-dist, make_pair(0,i)));
    //printf("Adding %d: %lf %lf %lf\n", i, data[i][0], data[i][1], dist);
  }

  while(visitedCount<data.size())
  {
    pair<double, pair<int,int> > edge = myQueue.top();
    int i = edge.second.first;
    int j = edge.second.second;
    myQueue.pop();

    if (visited[i] && visited[j])
      continue;

    int toAdd;
    if(visited[i])
      toAdd = j;
    else
      toAdd = i;
    visited[toAdd]=true;
    visitedCount++;
    //printf("Popped %d\n", toAdd);

    for(int i=0; i<data.size(); i++)
    {
      if(!visited[i])
      {
        double dist = length(data[toAdd] - data[i]);
        myQueue.push(make_pair(-dist, make_pair(toAdd,i)));
        //printf("Adding %d: %lf %lf %lf\n", i, data[i][0], data[i][1], dist);
      }
    }

    graph[toAdd].push_back(i+j-toAdd);
    graph[i+j-toAdd].push_back(toAdd);
	
  }

}

vector<int> ETSP::match()
{
  vector<int> oddNodes;
  vector<int> bestMatch;
  for(int i=0; i<data.size(); i++)
  {
    if(graph[i].size()%2==1)
    {
      oddNodes.push_back(i);
      bestMatch.push_back(i);
    }
  }
  
  //for(int i=0; i<oddNodes.size(); i++)
  //  printf("%d ",oddNodes[i]);
  //printf("\n");
  if(oddNodes.size()>0)
  {
    double best=1e10;
    do 
    {
       double curr = 0;
       for(int i=0; i<oddNodes.size(); i+=2)
       {
         curr+=length(data[oddNodes[i]] - data[oddNodes[i+1]]);
       }
       //printf("\t%lf %lf\n", curr,best);
       if(curr<best)
       {
         best=curr;
         for(int j=0; j<oddNodes.size(); j++)
         {
           bestMatch[j]=oddNodes[j];
         }
         //printf("\tnow: %lf %lf\n", curr,best);
       }
    } while (std::next_permutation(oddNodes.begin(), oddNodes.end()) );
  }
  //for(int i=0; i<bestMatch.size(); i++)
  //  printf("%d ",bestMatch[i]);
  return bestMatch;
}

void ETSP::orderFromMST()
{
  vector<bool> visited;
  visited.resize(data.size(),false);
  int curr=0;
  inverse_mapped_order.resize(data.size());
  traverse(0, visited);
  /*for(int i=0; i<order.size(); i++)
    printf("%d ",order[i]);
  printf("\n");*/
}

void ETSP::traverse(int curr, vector<bool> &visited)
{
  visited[curr]=true;
  order.push_back(curr);
  inverse_mapped_order[curr] = order.size()-1;
  for(int i=0; i<graph[curr].size(); i++)
  {
    if(visited[graph[curr][i]])
      continue;
    traverse(graph[curr][i], visited);
  }
}

void ETSP::orderFromMatch()
{
  vector<int> myMatch = match();
  //for(int i=0; i<myMatch.size(); i++)
  //  printf("%d ",myMatch[i]);
  //printf("\n");
  for(int i=0; i<myMatch.size(); i+=2)
  {
    int start = inverse_mapped_order[myMatch[i]];
    int end = inverse_mapped_order[myMatch[i+1]];
    if(start>end)
    {
      int temp=start;
      start=end;
      end=temp;
    }
    //printf("reordering: %d %d\n", start, end);
    for(int j=start+1; j<end; j++, end--)
    {
      int temp = order[j];
      order[j] = order[end];
      order[end] = temp;
      inverse_mapped_order[order[j]] = j;
      inverse_mapped_order[order[end]] = end;
    }
  }
  //for(int i=0; i<order.size(); i++)
  //  printf("%d ",order[i]);
  //printf("\n");
}

