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

#ifndef _STRUCTURE_H_
#define _STRCUTURE_H_

#include "ScatteredData.h"
#include "RBF.h"
#include "vec3.h"
#include "Surface.h"
#include "fault.h"
#include "tree.h"

#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>

class Structure
{
	
	std::vector<Fault*> fault;
	std::vector<Surface*> horizon;

	public:
	Structure(std::vector<Fault*> &myFault, std::vector<Surface*> &myHorizon);

	StructureTree *tree;
	std::vector<StructureTree*> nodes;
	std::vector<std::vector<int> > neighbors;


	void constructTree(std::string filename);    //The input has to come from the user
	void readStructureFile(std::string filename, std::vector<int> &nodeToFeature, std::vector<std::vector<int> > &children, std::vector<int> &whichFeature);
	void computeRBFs();
	std::pair<int,double> returnHorizon(vec3 x);
};


#endif //_STRUCTURE_H_
