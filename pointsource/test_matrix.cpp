/*
 * Test Matrix of matrix.h
 *
 * 2018-01-07
 */

#include <iostream>
#include <vector>

#include "matrix.h"

using namespace std;


int main(void)
{
    int rows = 3, cols = 4;
    Matrix<double> mat(rows, cols);
    mat(1,1) = 1.7;
    mat(2,2) = 5.4;

    cout << "Matrix[" << rows << "," << cols << "]:" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << mat(i, j) << ", ";
        }
        cout << endl;
    }

    // vector of matrices
    cout << "-----------------------------------------------" << endl;;
    cout << "Vector of matrices ..." << endl;
    vector<Matrix<double> > mats;
    mats.push_back(mat);
    mats.push_back(Matrix<double>(3,3));
    mats.push_back(Matrix<double>(4,4));
    mats.push_back(Matrix<double>(5,5));
    mats.push_back(Matrix<double>(6,3));

    cout << "Number of matrices: " << mats.size() << endl;
    cout << "Iterator:" << endl;
    for (vector<Matrix<double> >::iterator it = mats.begin();
         it != mats.end();
         ++it) {
        cout << "Matrix: [" << it->rows() << ","
             << it->cols() << "]" << endl;
    }

    cout << "Indexer:" << endl;
    for (vector<Matrix<double> >::size_type i = 0;
         i < mats.size();
         ++i) {
        cout << "Matrix: [" << mats[i].rows() << ","
             << mats[i].cols() << "]" << endl;
    }

    return 0;
}
