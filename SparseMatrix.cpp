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

