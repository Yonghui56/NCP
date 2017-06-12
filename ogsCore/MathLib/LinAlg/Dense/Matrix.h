/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Matrix.h
 *
 * Created on 2011-05-24 by Thomas Fischer
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <new>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <cassert>

namespace MathLib {

/**
 * Matrix represents a dense matrix for a numeric data type.
 */
template <class T> class Matrix
{
public:
   Matrix();
   Matrix (size_t rows, size_t cols);
   Matrix (size_t rows, size_t cols, const T& val);
   Matrix (const Matrix<T> &src);

   ~Matrix ();

   void resize(size_t rows, size_t cols);
   size_t getNRows () const { return nrows; }
   size_t getNCols () const { return ncols; }
   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy ( T alpha, const T* x, T beta, T* y) const;

   /**
    * Assignment operation for all matrix entities
    * @param a
    * @return
    */
   Matrix<T>& operator= (const Matrix<T> &src);

   /**
    * Set operation for all matrix entities
    * @param a
    * @return
    */
   Matrix<T>& operator= (T a);
   /**
    * Matrix scalar multiplication
    * @param a
    * @return
    */
   void operator*= (T a);
   /**
    * Matrix scalar divide
    * @param a
    * @return
    */
   void operator/= (T a);
   
   /**
    * Matrix vector multiplication
    * @param x
    * @return
    */
   T* operator* (const T *x) const;
   /**
    * Matrix matrix addition.
    * @param mat
    * @return
    */
   Matrix<T>* operator+ (const Matrix<T>& mat) const throw (std::range_error);
   /**
    * Matrix matrix addition.
    * @param mat
    * @return
    */
   Matrix<T>& operator+= (const Matrix<T>& mat) throw (std::range_error);
   /**
    * Matrix matrix subtraction
    * @param mat
    * @return
    */
   Matrix<T>* operator- (const Matrix<T>& mat) const throw (std::range_error);

   /**
    * Matrix matrix multiplication \f$ C = A \cdot B\f$
    * @param mat the matrix \f$ B \f$
    * @return the matrix \f$ C \f$
    */
   Matrix<T>* operator* (const Matrix<T>& mat) const throw (std::range_error);

   /**
    * matrix transpose
    * @return the transpose of the matrix
    */
   Matrix<T>* transpose() const; // HB & ZC

   /**
    * matrix transpose
    * @param mat a matrix where transpose result is stored
    */
   void transpose(Matrix<T> &mat) const;

   Matrix<T>* getSubMatrix (size_t b_row, size_t b_col, size_t e_row, size_t e_col) const throw (std::range_error);

   /**
    * overwrites values of the matrix with the given sub matrix
    * @param b_row the first row
    * @param b_col the first column
    * @param sub_mat the sub matrix
    */
   void setSubMatrix (size_t b_row, size_t b_col, const Matrix<T>& sub_mat) throw (std::range_error);

   /**
    * add values of the matrix with the given symmetric matrix
    * @param pos row and column positions in the matrix
    * @param mat the given matrix
    */
   void add(std::vector<size_t> &pos, const Matrix<T>& mat);

   /**
    * add values of the matrix with the given matrix
    * @param row_pos row positions in the matrix
    * @param col_pos column positions in the matrix
    * @param mat the given matrix
    */
   void add(std::vector<size_t> &row_pos, std::vector<size_t> &col_pos, const Matrix<T>& mat);

   /**
    * multiply values of this matrix A with the given matrix B and add results into a matrix C, i.e. C += A*B
    * @param matB the matrix B
    * @param matC the matrix C
    */
   void multiply(const Matrix<T> &matB, Matrix<T> &matC) const;

   /**
    * multiply values of this transposed matrix A with the given matrix B and add results into a matrix C, i.e. C += a*A^T*B
    * @param matB the transposed matrix B
    * @param matC the matrix C
    * @param a the scalar factor. default = 1.0
    */
   void transposeAndMultiply(const Matrix<T> &matB, Matrix<T> &matC, double a=1.0) const;
   /**
    * C += a*A^T*v*B
    * @param matB the transposed matrix B
    * @param matC the matrix C
    * @param a the scalar factor. default = 1.0
    */
   void transposeAndMultiply(const Matrix<T> &matB, T*row_vec, Matrix<T> &matC, double a=1.0) const;

   /**
    * return a determinant of the matrix
    * @return determinant
    */
   T determinant() const;

   /**
    * return an inverse matrix
    * @return inverse matrix
    */
   Matrix<T>* inverse() const;
   /**
    * return an inverse matrix
    * @param inverse matrix
    */
   void inverse(Matrix<T>*) const;

   inline T & operator() (size_t row, size_t col) throw (std::range_error);
   inline T & operator() (size_t row, size_t col) const throw (std::range_error);

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

   T const* getData () { return data; }

private:
   // zero based addressing, but Fortran storage layout
   //inline size_t address(size_t i, size_t j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline size_t address(size_t i, size_t j) const { return i*ncols+j; };

   size_t nrows;
   size_t ncols;
   T *data;
};

template<class T> Matrix<T>::Matrix ()
    : nrows (0), ncols (0), data (0)
{}

template<class T> Matrix<T>::Matrix (size_t rows, size_t cols)
      : nrows (rows), ncols (cols), data (new T[nrows*ncols])
{}

template<class T> Matrix<T>::Matrix (size_t rows, size_t cols, const T &val)
    : nrows (rows), ncols (cols), data (new T[nrows*ncols])
{
    (*this) = val;
}

template<class T> Matrix<T>::Matrix (const Matrix& src) :
    nrows (src.getNRows ()), ncols (src.getNCols ()), data (new T[nrows * ncols])
{
   for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
         data[address(i,j)] = src (i, j);
}

template<class T> Matrix<T>& Matrix<T>::operator= (const Matrix<T> &src)
{
    assert(nrows == src.getNRows() && ncols==src.getNCols());
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
            data[address(i,j)] = src (i, j);
    return *this;
}

template <class T> Matrix<T>::~Matrix ()
{
   delete [] data;
}

template<class T> void Matrix<T>::resize (size_t rows, size_t cols) 
{
    if (rows!=nrows || cols!=ncols) {
        delete [] data;
        data = new T[rows*cols];
        nrows = rows;
        ncols = cols;
    }
}

template<class T> Matrix<T>& Matrix<T>::operator= (T a)
{
   for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
         data[address(i,j)] = a;

    return *this;
}

template<class T> void Matrix<T>::operator*= (T a)
{
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
            data[address(i,j)] *= a;

}

template<class T> void Matrix<T>::operator/= (T a)
{
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
            data[address(i,j)] /= a;

}

template<class T> void Matrix<T>::axpy ( T alpha, const T* x, T beta, T* y) const
{
   for (size_t i(0); i<nrows; i++) {
      y[i] += beta * y[i];
      for (size_t j(0); j<ncols; j++) {
         y[i] += alpha * data[address(i,j)] * x[j];
      }
   }
}

template<class T> T* Matrix<T>::operator* (const T *x) const
{
    T *y (new T[nrows]);
    for (size_t i(0); i < nrows; i++) {
        y[i] = 0.0;
        for (size_t j(0); j < ncols; j++) {
            y[i] += data[address(i, j)] * x[j];
        }
    }

    return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator+ (const Matrix<T>& mat) const throw (std::range_error)
{
    // make sure the two matrices have the same dimension.
    if (nrows != mat.getNRows() || ncols != mat.getNCols())
        throw std::range_error("Matrix::operator+, illegal matrix size!");

    Matrix<T>* y(new Matrix<T> (nrows, ncols));
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            (*y)(i, j) = data[address(i, j)] + mat(i, j);
        }
    }

    return y;
}

template<class T> Matrix<T>& Matrix<T>::operator+= (const Matrix<T>& mat) throw (std::range_error)
{
    // make sure the two matrices have the same dimension.
    if (nrows != mat.getNRows() || ncols != mat.getNCols())
        throw std::range_error("Matrix::operator+, illegal matrix size!");

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            data[address(i, j)] += mat(i, j);
        }
    }

    return *this;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator- (const Matrix<T>& mat) const throw (std::range_error)
{
    // make sure the two matrices have the same dimension.
    if (nrows != mat.getNRows() || ncols != mat.getNCols())
        throw std::range_error("Matrix::operator-, illegal matrix size!");

    Matrix<T>* y(new Matrix<T> (nrows, ncols));
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            (*y)(i, j) = data[address(i, j)] - mat(i, j);
        }
    }

    return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator* (const Matrix<T>& mat) const throw (std::range_error)
{
    // make sure the two matrices have the same dimension.
    if (ncols != mat.getNRows())
        throw std::range_error(
                "Matrix::operator*, number of rows and cols should be the same!");

    size_t y_cols(mat.getNCols());
    Matrix<T>* y(new Matrix<T> (nrows, y_cols, T(0)));

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < y_cols; j++) {
            for (size_t k = 0; k < ncols; k++)
                (*y)(i, j) += data[address(i, k)] * mat(k, j);
        }
    }

    return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::transpose() const
{
    Matrix<T>* y(new Matrix<T> (ncols, nrows));

    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
//            y->data[y->address(j, i)] = data[address(i, j)];
            (*y)(j,i) = data[address(i, j)];
        }
    }
    return y;
}

template<class T> void Matrix<T>::transpose(Matrix<T> &mat) const
{
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            mat(j,i) = data[address(i, j)];
        }
    }
}

template<class T> Matrix<T>* Matrix<T>::getSubMatrix(
        size_t b_row, size_t b_col,
        size_t e_row, size_t e_col) const throw (std::range_error)
{
    if (b_row >= e_row | b_col >= e_col)
        throw std::range_error ("Matrix::getSubMatrix() illegal sub matrix");
    if (e_row > nrows | e_col > ncols)
        throw std::range_error ("Matrix::getSubMatrix() illegal sub matrix");

    Matrix<T>* y(new Matrix<T> (e_row-b_row, e_col-b_col));
    for (size_t i=b_row; i<e_row; i++) {
        for (size_t j=b_col; j<e_col; j++) {
            (*y)(i-b_row, j-b_col) = data[address(i, j)];
        }
    }
    return y;
}

template<class T> void Matrix<T>::setSubMatrix(
        size_t b_row, size_t b_col, const Matrix<T>& sub_mat) throw (std::range_error)
{
    if (b_row + sub_mat.getNRows() > nrows | b_col + sub_mat.getNCols() > ncols)
        throw std::range_error ("Matrix::setSubMatrix() sub matrix to big");

    for (size_t i=0; i<sub_mat.getNRows(); i++) {
        for (size_t j=0; j<sub_mat.getNCols(); j++) {
            data[address(i+b_row, j+b_col)] = sub_mat(i,j);
        }
    }
}

template<class T> void Matrix<T>::add(std::vector<size_t> &pos, const Matrix<T>& mat)
{
    add(pos, pos, mat);
}

template<class T> void Matrix<T>::add(std::vector<size_t> &row_pos, std::vector<size_t> &col_pos, const Matrix<T>& mat)
{
    for (size_t i=0; i<row_pos.size(); i++) {
        for (size_t j=0; j<col_pos.size(); j++) {
            data[address(row_pos[i], col_pos[j])] += mat(i,j);
        }
    }
}

template<class T> void Matrix<T>::multiply(const Matrix<T> &matB, Matrix<T> &matC) const
{
    for (size_t rowA=0; rowA<this->nrows; rowA++) {
        for (size_t colB=0; colB<matB.getNCols(); colB++) {
            for (size_t colArowB=0; colArowB<this->ncols; colArowB++) {
                matC(rowA, colB) += (*this)(rowA, colArowB) * matB(colArowB, colB);
            }
        }
    }
}

template<class T> void Matrix<T>::transposeAndMultiply(const Matrix<T> &matB, Matrix<T> &matC, double a) const
{
    assert(this->nrows==matB.getNRows());
    assert(this->ncols==matB.getNCols());
    assert(this->ncols==matC.getNCols());
    assert(this->ncols==matC.getNRows());

    // C += A^T * B
    for (size_t colA=0; colA<this->ncols; colA++) {
        for (size_t colB=0; colB<matB.getNCols(); colB++) {
            for (size_t rowArowB=0; rowArowB<this->nrows; rowArowB++) {
                matC(colA, colB) += a * (*this)(rowArowB, colA) * matB(rowArowB, colB);
            }
        }
    }
}

template<class T> void Matrix<T>::transposeAndMultiply(const Matrix<T> &matB, T*row_vec, Matrix<T> &matC, double a) const
{
    Matrix<T> m(1, matB.getNCols());
    m = .0;
    // M = v * B
    for (size_t colB=0; colB<matB.getNCols(); colB++) {
        for (size_t rowB=0; rowB<matB.getNRows(); rowB++) {
            m(0, colB) += row_vec[rowB]*matB(rowB, colB);
        }
    }
    // C += a* A^T * M
    transposeAndMultiply(m, matC, a);
    //for (size_t colA=0; colA<this->ncols; colA++) {
    //    for (size_t colB=0; colB<m.getNCols(); colB++) {
    //        for (size_t colAcolB=0; colAcolB<this->nrows; colAcolB++) {
    //            matC(colA, colB) += a * (*this)(colA, colAcolB) * m(colAcolB, colB);
    //        }
    //    }
    //}
}

template<class T> T Matrix<T>::determinant() const
{
    assert(nrows==ncols);
    assert(nrows<4);

    T det;
    if (nrows==1) {
        det = (*this)(0,0);
    } else if (nrows==2) {
        det = (*this)(0,0)*(*this)(1,1);
        det -= (*this)(0,1)*(*this)(1,0);
    } else if (nrows==3) {
        det = (*this)(0,0)*((*this)(1,1)*(*this)(2,2)-(*this)(2,1)*(*this)(1,2));
        det -= (*this)(1,0)*((*this)(0,1)*(*this)(2,2)-(*this)(2,1)*(*this)(0,2));
        det += (*this)(2,0)*((*this)(0,1)*(*this)(1,2)-(*this)(1,1)*(*this)(0,2));
    } else {
        //TODO not implemented yet
    }

    return det;
}

template<class T> void Matrix<T>::inverse(Matrix<T> *y) const
{
    assert(nrows==ncols);
    assert(nrows<4);

    if (nrows==1) {
        (*y)(0,0) = 1.0 / (*this)(0,0);
    } else {
        if (nrows==2) {
            (*y)(0,0) = (*this)(1,1);
            (*y)(0,1) = -(*this)(0,1);
            (*y)(1,0) = -(*this)(1,0);
            (*y)(1,1) = (*this)(0,0);
        } else if (nrows==3) {
            (*y)(0,0) = (*this)(1,1)*(*this)(2,2)-(*this)(2,1)*(*this)(1,2);
            (*y)(0,1) = (*this)(0,2)*(*this)(2,1)-(*this)(0,1)*(*this)(2,2);
            (*y)(0,2) =  (*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1);
            //
            (*y)(1,0) =  (*this)(1,2)*(*this)(2,0)-(*this)(2,2)*(*this)(1,0);
            (*y)(1,1) =  (*this)(0,0)*(*this)(2,2)-(*this)(2,0)*(*this)(0,2);
            (*y)(1,2) =  (*this)(0,2)*(*this)(1,0)-(*this)(1,2)*(*this)(0,0);
            //
            (*y)(2,0) =  (*this)(1,0)*(*this)(2,1)-(*this)(2,0)*(*this)(1,1);
            (*y)(2,1) =  (*this)(0,1)*(*this)(2,0)-(*this)(2,1)*(*this)(0,0);
            (*y)(2,2) =  (*this)(0,0)*(*this)(1,1)-(*this)(1,0)*(*this)(0,1);
        } else {
            //error
            std::cout << "***Error in Matrix::inverse(): nrows>3 is not supported." << std::endl;
        }
        double detA = determinant();
        (*y) /= detA;
    }
}

template<class T> Matrix<T>* Matrix<T>::inverse() const
{
    assert(nrows==ncols);
    assert(nrows<4);

    Matrix<T>* y(new Matrix<T> (nrows, ncols));
    inverse(y);

    return y;
}

template<class T> T& Matrix<T>::operator() (size_t row, size_t col)
    throw (std::range_error)
{
   if ( (row >= nrows) | ( col >= ncols) )
      throw std::range_error ("Matrix: op() const range error");
   return data [address(row,col)];
}


template<class T> T& Matrix<T>::operator() (size_t row, size_t col) const
    throw (std::range_error)
{
   if ( (row >= nrows) | ( col >= ncols) )
      throw std::range_error ("Matrix: op() const range error");
   return data [address(row,col)];
}

template <class T> void Matrix<T>::write (std::ostream &out) const
{
    for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
            out << data[address(i, j)] << "\t";
        }
        out << std::endl;
    }
}

template <class T> T sqrFrobNrm (const Matrix<T> &mat)
{
    T nrm ((T)(0));
    size_t i,j;
    for (j=0; j<mat.getNCols(); j++)
        for (i=0; i<mat.getNRows(); i++)
            nrm += mat(i,j) * mat(i,j);

    return nrm;
}

/** overload the output operator for class Matrix */
template <class T>
std::ostream& operator<< (std::ostream &os, const Matrix<T> &mat)
{
    mat.write (os);
    return os;
}

} // end namespace MathLib

#endif
