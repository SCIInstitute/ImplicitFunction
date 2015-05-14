#ifndef _FILEIO_H_
#define _FILEIO_H_


#include <cstdio>
#include <vector>
#include <string>
#include <fstream>

#include "vec3.h"

using std::vector;
using std::string;

void writeNrrdFile(string filename, vector<vector<vector<double> > > &myArray, vec3 myDim, vec3 mySpacing)
{
        std::cout << "Writing file '" << filename << "'" << std::endl;
	std::ofstream nrrd_file(filename.c_str(), std::ofstream::binary);

        if(nrrd_file.is_open())
        {
            nrrd_file << "NRRD0001" << std::endl;
            nrrd_file << "# Complete NRRD file format specification at:" << std::endl;
            nrrd_file << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
            nrrd_file << "type: float" << std::endl;
            nrrd_file << "dimension: 3" << std::endl;
            nrrd_file << "sizes: " << myDim[0] << " " << myDim[1] << " " << myDim[2] << std::endl;
            nrrd_file << "axis mins: " << 0 << ", " << 0 << ", " << 0 << std::endl;
            nrrd_file << "spacings: " << mySpacing[0] << " " << mySpacing[1] << " " << mySpacing[2] << std::endl;
            nrrd_file << "centerings: cell cell cell" << std::endl;
            nrrd_file << "endian: little" << std::endl;
            nrrd_file << "encoding: raw" << std::endl;
            nrrd_file << std::endl;

            // write data portion
            for(int k=0; k < myDim[2]; k++)
            {
                for(int j=0; j < myDim[1]; j++)
                {
                    for(int i=0; i < myDim[0]; i++)
                    {
                        float val = myArray[i][j][k];
                        nrrd_file.write((char*)&val, sizeof(float));
                    }
                }
            }

            nrrd_file.close();
        }
}

void augmentData(ScatteredData *data)
{
	int n = data->x[0].size();
	int increment[] = {0,0,100};
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] + increment[j]);
		}
		data->fnc.push_back(10);
	}
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] - increment[j]);
		}
		data->fnc.push_back(-10);
	}
}

vec3 findNormal(ScatteredData *data, int n)
{
	int tot = data->x[0].size();
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
	//printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);

	vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
	vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
	vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
	/*vec3 one = b-a;
	vec3 two = c-a;
	vec3 ret = one+two;*/

	vec3 tangent = b-c;
	//rotate by 90 degrees on the x-y plane
	double ret_x = -tangent[1];
	double ret_y = tangent[0];
	vec3 ret(ret_x, ret_y, tangent[2]);
	return ret;
}

void augmentNormalData(ScatteredData *data)
{
	int n = data->x[0].size();
	for(int i=0; i<n; i++)
	{
		vec3 myNormal = findNormal(data, i);
		myNormal = normalize(myNormal);
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] + myNormal[j]);
		}
		data->fnc.push_back(10);

		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] - myNormal[j]);
		}
		data->fnc.push_back(-10);
	}
}

void readDataFile(string filename, ScatteredData *data)
{
	std::cout <<"Reading file '" << filename <<"'"<<std::endl;
	std::ifstream datafile(filename.c_str());
	if(datafile.is_open())
	{
		std::cout<<"Opened"<<std::endl;
		double myData[3];
		while(!datafile.eof())
		{
			datafile>>myData[0]>>myData[1]>>myData[2];
			for(int i=0; i<3; i++)
				data->x[i].push_back(myData[i]);
			data->fnc.push_back(0);
		}
		std::cout<<"Done"<<std::endl;
		std::cout<<"Augmenting data"<<std::endl;
		augmentData(data);
		std::cout<<"Done"<<std::endl;
	}
	else
	{
		std::cout <<"Could not open '" << filename <<"'"<<std::endl;
	}
}


void readSurfaceDataFile(string filename, ScatteredData *data)
{
	std::cout <<"Reading file '" << filename <<"'"<<std::endl;
	std::ifstream datafile(filename.c_str());
	if(datafile.is_open())
	{
		std::cout<<"Opened"<<std::endl;
		double myData[3];
		while(!datafile.eof())
		{
			datafile>>myData[0]>>myData[1]>>myData[2];
			for(int i=0; i<3; i++)
				data->x[i].push_back(myData[i]);
			data->fnc.push_back(0);
		}
		for(int i=0; i<3; i++)
			data->x[i].pop_back();
		data->fnc.pop_back();
		std::cout<<"Done"<<std::endl;
		data->computeOrdering();
		std::cout<<"Augmenting data"<<std::endl;
		augmentNormalData(data);
		std::cout<<"Done"<<std::endl;
	}
	else
	{
		std::cout <<"Could not open '" << filename <<"'"<<std::endl;
	}
}



#endif //_FILEIO_H_
