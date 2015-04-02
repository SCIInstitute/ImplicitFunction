#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

// STL Includes
#include <vector>
using std::vector;

class SparseMatrix
{
public:
    SparseMatrix();
    SparseMatrix(int n);
    void push_back(int myRow, int myCol, double myVal);
    void resize(int n);
    void multiply(vector<double> &b, vector<double> &c);
    vector<double> multiply(vector<double> &b);


private:
    vector< vector<int> > col;
    vector< vector<double> > val;
};

#endif //_SPARSEMATRIX_H_


