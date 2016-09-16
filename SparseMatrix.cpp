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

// STL Includes
#include <vector>
#include "SparseMatrix.h"

using std::vector;

SparseMatrix::SparseMatrix()
{
}

SparseMatrix::SparseMatrix(int n)
{
  col.resize(n);
  val.resize(n);
}

void SparseMatrix::resize(int n)
{
  col.resize(n);
  val.resize(n);
}

void SparseMatrix::push_back(int myRow, int myCol, double myVal)
{
  col[myRow].push_back(myCol);
  val[myRow].push_back(myVal);
}

void SparseMatrix::multiply(vector<double> &b, vector<double> &c)
{
  int n = b.size(), i, j;
  c.clear();
  c.resize(n);

  for(i = 0; i<this->col.size(); i++)
  {
    int row, col;
    row = i;
    for(j=0; j<this->col[i].size(); j++)
    {
      col = this->col[i][j];
      //cout<<" Row, Col = "<<row<<", "<<col<<endl;
      c[row] += this->val[i][j]*b[col];
    }
  }
}

std::vector<double> SparseMatrix::multiply(vector<double> &b)
{
  vector<double> c;
  multiply(b,c);
  return c;
}

