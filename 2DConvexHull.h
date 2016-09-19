#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>
#include "vec3.h"

using std::vector;
using std::stack;
 
struct Point
{
    double x, y;
    int index;
};
 
// A utility function to find next to top in a stack
Point nextToTop(stack<Point> &S)
{
    Point p = S.top();
    S.pop();
    Point res = S.top();
    S.push(p);
    return res;
}
 
void swap(Point &p1, Point &p2)
{
    Point temp = p1;
    p1 = p2;
    p2 = temp;
}
 
double distSq(Point p1, Point p2)
{
    return (p1.x - p2.x)*(p1.x - p2.x) +
          (p1.y - p2.y)*(p1.y - p2.y);
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 -> p, q and r are colinear
// 1 -> Clockwise
// 2 -> Counterclockwise
int orientation(Point p, Point q, Point r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
    return (val > 0)? 1: 2; // clock or counterclock wise
}
 
// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void *vp1, const void *vp2)
{
   Point *p1 = (Point *)vp1;
   Point *p2 = (Point *)vp2;

   Point p0; p0.x=0; p0.y=0;
 
   // Find orientation
   int o = orientation(p0, *p1, *p2);
   if (o == 0)
     return (distSq(p0, *p2) >= distSq(p0, *p1))? -1 : 1;
 
   return (o == 2)? -1: 1;
}
 
// Prints convex hull of a set of n points.
stack<Point> convexHull(Point points[], int n)
{
   // Find the bottommost point
   int ymin = points[0].y, min = 0;
   for (int i = 1; i < n; i++)
   {
     int y = points[i].y;
 
     // Pick the bottom-most or chose the left
     // most point in case of tie
     if ((y < ymin) || (ymin == y &&
         points[i].x < points[min].x))
        ymin = points[i].y, min = i;
   }
 
   // Place the bottom-most point at first position
   swap(points[0], points[min]);
 
   // Sort n-1 points with respect to the first point.
   Point p0;
   for(int i=1; i<n; i++)
   {
     points[i].x -= points[0].x;
     points[i].y -= points[0].y;
   }
   p0.x = points[0].x = 0;
   p0.y = points[0].y = 0;

   qsort(&points[1], n-1, sizeof(Point), compare);
 
   // If two or more points make same angle with p0,
   // Remove all but the one that is farthest from p0
   int m = 1; 
   for (int i=1; i<n; i++)
   {
       while (i < n-1 && orientation(p0, points[i], points[i+1]) == 0)
          i++;
       points[m] = points[i];
       m++;  // Update size of modified array
   }
 
   stack<Point> S;
   if (m < 3) return S;
 
   S.push(points[0]);
   S.push(points[1]);
   S.push(points[2]);
 
   // Process remaining n-3 points
   for (int i = 3; i < m; i++)
   {
      // Keep removing top while the angle formed by
      // points next-to-top, top, and points[i] makes
      // a non-left turn
      while (orientation(nextToTop(S), S.top(), points[i]) != 2)
         S.pop();
      S.push(points[i]);
   }
   return S;
}

vector<int> getConvexHull(vector<vec3> &inPoints, int dim)  
//dim refers to the dimention that needs to be ignored
{
  Point *myPoints;
  myPoints = new Point[inPoints.size()];
  for(int i = 0; i < inPoints.size(); i++)
  {
    int index = 0;
    for(int j = 0; j < 3; j++)
    {
      if (j == dim) continue;
      if (index == 0) myPoints[i].x = inPoints[i][j];
      if (index == 1) myPoints[i].y = inPoints[i][j];
      index++;
    }
    myPoints[i].index = i;
  }
  stack<Point> myStack = convexHull(myPoints, inPoints.size());

  vector<int> ret;
  while (!myStack.empty())
  {
    Point p = myStack.top();
    ret.push_back(p.index);
    myStack.pop();
  }
  return ret;
}

