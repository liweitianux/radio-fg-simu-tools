/*
 * Matrix based on valarray
 *
 * Credit:
 * * https://stackoverflow.com/a/2188001
 * * https://github.com/CaptGreg/SenecaOOP345-attic/blob/master/valarray-matrix.cpp
 *
 * References: http://www.stroustrup.com/matrix.c
 */

#ifndef MATRIX_H_
#define MATRIX_H_


#include <valarray>

template <typename T>
class Matrix {
public:
    // create an empty matrix
    Matrix(std::size_t rows, std::size_t cols):
        rows_(rows),
        cols_(cols),
        data_(rows*cols)
    { }

    // get the number of rows
    std::size_t rows() const { return rows_; }
    // get the number of columns
    std::size_t cols() const { return cols_; }
    // get a copy of the data
    std::valarray<T> array() const { return data_; }

    // retrieve the data from row `r`
    std::valarray<T> row(std::size_t r) const {
        return data_[std::slice(r * cols(), cols(), 1)];
    }

    // retrieve the data from column `c`
    std::valarray<T> col(std::size_t c) const {
        return data_[std::slice(c, rows(), cols())];
    }

    // retrieve reference to the data from row `r`
    std::slice_array<T> row(std::size_t r) {
        return data_[std::slice(r * cols(), cols(), 1)];
    }

    // retrieve reference to the data from column `c`
    std::slice_array<T> col(std::size_t c) {
        return data_[std::slice(c, rows(), cols())];
    }

    // basic item reference
    T& operator() (std::size_t r, std::size_t c) {
        return data_[r * cols() + c];
    }

    // basic item retrieval
    T operator() (std::size_t r, std::size_t c) const {
        return row(r)[c];
    }

private:
    size_t rows_;
    size_t cols_;
    std::valarray<T> data_;
};


#endif  /* MATRIX_H_ */
