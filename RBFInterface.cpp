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

#include "RBFInterface.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

const double RBFInterface::EPSILON = 1.0e-3;

RBFInterface::RBFInterface(std::vector<vec3> myData,
                           vec3 myOrigin, vec3 mySize, vec3 mySpacing,
                           double myOffset, std::vector<axis_t> myAxis, Kernel kernel) :
  thresholdValue(0),
  myKernel(kernel)
{
	CreateSurface(myData, myOrigin, mySize, mySpacing, myOffset, myAxis);
}

void RBFInterface::setOffset(double myOffset)
{
	offset = myOffset;
}

void RBFInterface::CreateSurface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing, double myOffset, std::vector<axis_t> myAxis)
{
	std::vector<double> a,b,c,d;
	for(int i=0; i<myData.size(); i++)
	{
    std::cerr << "point : " << myData[i][0] << ", " << myData[i][1] << ", " << myData[i][2] << std::endl;
		a.push_back(myData[i][0]);
		b.push_back(myData[i][1]);
		c.push_back(myData[i][2]);
		d.push_back(thresholdValue);
	}

  // This code figures out the bounds for the data selected

  std::vector<double>::iterator minx = std::min_element(a.begin(), a.end());
  std::vector<double>::iterator miny = std::min_element(b.begin(), b.end());
  std::vector<double>::iterator minz = std::min_element(c.begin(), c.end());
  std::vector<double>::iterator maxx = std::max_element(a.begin(), a.end());
  std::vector<double>::iterator maxy = std::max_element(b.begin(), b.end());
  std::vector<double>::iterator maxz = std::max_element(c.begin(), c.end());
  vec3 myMin(*minx, *miny, *minz), myMax(*maxx, *maxy, *maxz);
  myMin = myMin - 0.05*mySize;
  myMax = myMax + 0.05*mySize;

	mySurfaceData = new ScatteredData(a,b,c,d);
	mySurfaceData->axisInformation = myAxis;
        mySurfaceData->compute2DHull();
	mySurfaceData->origSize = mySurfaceData->x[0].size();

	augmentNormalData(mySurfaceData, myOffset);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	mySurfaceRBF->setDataReduction(All);

  // TODO: let caller pick the kernel
	myKernel = ThinPlate;
	mySurface = new Surface(mySurfaceData, mySurfaceRBF);

	// Construct RBFs
	mySurface->computeRBF();

  // sanity check
  for(int i=0; i<mySurfaceData->fnc.size(); i++)
  {
    vec3 myLocation(mySurfaceData->x[0][i], mySurfaceData->x[1][i], mySurfaceData->x[2][i]);
    double myVal = mySurface->computeValue(myLocation);
    double error  = fabs(myVal - mySurfaceData->fnc[i]);
    if (error > EPSILON)
    {
      printf("ERROR (numerics): %lf\n", error);
      fflush(stdout);
    }
  }

  // Fill the values into the vector. In the first loop, we initialize the matrix with all values set to -100 (or - inf). In the second loop, we change the values from -100 to the correct value if the point in the domain described above.

  value.resize(static_cast<int>(mySize[0]));
  for(int i=0; i<mySize[0]; i++)
  {
    value[i].resize(static_cast<int>(mySize[1]));
    for(int j=0; j<mySize[1]; j++)
    {
      value[i][j].resize(static_cast<int>(mySize[2]), -100);
    }
  }

  for(int i=0; i<mySize[0]; i++)
  {
    vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX;
    if (location[0]<myMin[0]||location[0]>myMax[0])
      continue;
    for(int j=0; j<mySize[1]; j++)
    {
      location = myOrigin + mySpacing[1]*j*vec3::unitY;
      if (location[1]<myMin[1]||location[1]>myMax[1])
        continue;
      for(int k=0; k<mySize[2]; k++)
      {
        location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
        if (location[2]<myMin[2]||location[2]>myMax[2])
          continue;
        //std::cout<<"Computing Val ... "<<std::endl;
        double myVal = mySurface->computeValue(location);
        //printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
        value[i][j][k]=myVal;
      }
    }
  }
}


vec3 RBFInterface::findNormal(ScatteredData *data, int n)
{
	int tot = data->origSize;
	int prev = (n-1)>=0?n-1:tot-1;
	int next = (n+1)<tot?n+1:0;

	while(data->x[2][prev]!=data->x[2][n])
	{
		prev = (prev-1)>=0?prev-1:tot-1;
	}

	while(data->x[2][next]!=data->x[2][n])
	{
		next = (next+1)<tot?next+1:0;
	}
	printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);

	vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
	vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
	vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
	
	vec3 tangent = b-c;
        //rotate by 90 degrees on the x-y plane
        double ret_x = -tangent[1];
        double ret_y = tangent[0];
        vec3 ret(ret_x, ret_y, tangent[2]);

	return ret;
}

vec3 RBFInterface::findSphericalNormal(ScatteredData *data, int n)
{
	vec3 ret(0,0,0);
	for(int j=0; j<3; j++)
		ret[j] = (data->x[j][n] - data->centroid[j])/10;

	return ret;
}

void RBFInterface::augmentNormalData(ScatteredData *data, double myOffset)
{
	printf("here\n");
	int n = data->origSize;
	for(int i=0; i<n; i++)
	{
		vec3 myNormal = findNormalAxis(data, i);
		myNormal = normalize(myNormal);
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] + myOffset*myNormal[j]);
		}
		data->fnc.push_back(10);

		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] - myOffset*myNormal[j]);
		}
		data->fnc.push_back(-10);
	}
}


vec3 RBFInterface::findNormalAxis(ScatteredData *data, int n)
{
	//printf("here\n");
	int tot = data->origSize;
	int prev = (n-1)>=0?n-1:tot-1;
	int next = (n+1)<tot?n+1:0;
	axis_t myAxis = data->axisInformation[n];

	while(data->x[myAxis][prev]!=data->x[myAxis][n])
	{
		prev = (prev-1)>=0?prev-1:tot-1;
	}

	while(data->x[myAxis][next]!=data->x[myAxis][n])
	{
		next = (next+1)<tot?next+1:0;
	}
	//printf(" see: %d %d %d %d\n", prev,n,next,tot); fflush(stdout);

	vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
	vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
	vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
	
	vec3 tangent = b-c;
	double ret_x, ret_y, ret_z;
	vec3 ret(tangent);
        //rotate by 90 degrees on the x-y plane
	switch(myAxis)
	{
		case 0:
			ret[1] = ret_y = -tangent[2];
			ret[2] = ret_z = tangent[1];
			break;
		case 1:
			ret[2] = ret_z = -tangent[0];
			ret[0] = ret_x = tangent[2];
			break;
		case 2:
			ret[0] = ret_x = -tangent[1];
			ret[1] = ret_y = tangent[0];
			break;
	}
	return ret;
}
