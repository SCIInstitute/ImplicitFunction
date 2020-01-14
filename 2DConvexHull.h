//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#ifndef TWODCONVEXHULL
#define TWODCONVEXHULL 1

#include <iostream>
#include <stack>
#include <vector>
#include <stdlib.h>

#include "vec3.h"

struct Point
{
  double x, y;
  int index;
};

// A utility function to find next to top in a stack
Point nextToTop(std::stack<Point> &S)
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

double distSq(const Point& p1, const Point& p2)
{
  return (p1.x - p2.x) * (p1.x - p2.x) +
         (p1.y - p2.y) * (p1.y - p2.y);
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 -> p, q and r are colinear
// 1 -> Clockwise
// 2 -> Counterclockwise
int orientation(const Point& p, const Point& q, const Point& r)
{
  int val = (q.y - p.y) * (r.x - q.x) -
            (q.x - p.x) * (r.y - q.y);

  if (val == 0) return 0;  // colinear
  // TODO: make enum or consts...
  return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void *vp1, const void *vp2)
{
  // TODO: C-style cast...
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
// TODO: points gets overwritten...
std::stack<Point> convexHull(std::vector<Point>& points)
{
  const int n = points.size();
  // Find the bottommost point
  int ymin = points[0].y, min = 0;
  for (int i = 1; i < n; i++)
  {
    int y = points[i].y;

    // Pick the bottom-most or chose the left
    // most point in case of tie
    if ( (y < ymin) || (ymin == y && points[i].x < points[min].x) )
      ymin = points[i].y, min = i;
  }

  // Place the bottom-most point at first position
  std::swap(points[0], points[min]);

  // Sort n-1 points with respect to the first point.
  Point p0;
  for(int i = 1; i < n; i++)
  {
    points[i].x -= points[0].x;
    points[i].y -= points[0].y;
  }
  p0.x = points[0].x = 0;
  p0.y = points[0].y = 0;

  qsort(&points[1], n-1, sizeof(Point), compare);

  // Process first 3 points

  // If two or more points make same angle with p0,
  // Remove all but the one that is farthest from p0
  int m = 1;
  for (int i = 1; i < n; i++)
  {
    while ( (i < n-1) && orientation(p0, points[i], points[i+1]) == 0 )
      i++;
    points[m] = points[i];
    m++;  // Update size of modified array
  }

  std::stack<Point> S;
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

// dim refers to the dimension that needs to be ignored
std::vector<int> getConvexHull(const std::vector<vec3> &inPoints, const int dim)
{
  std::vector<Point> myPoints(inPoints.size());
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

  // TODO: myPoints gets overwritten...
  auto myStack = convexHull(myPoints);

  std::vector<int> ret;
  while (! myStack.empty() )
  {
    Point p = myStack.top();
    ret.push_back(p.index);
    myStack.pop();
  }

  return ret;
}

#endif
