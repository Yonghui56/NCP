#ifndef EIGENINTERFACE_H
#define EIGENINTERFACE_H

#ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif

#include <iostream>
#include <fstream>
#include "LinAlg/Sparse/SparseTableCRS.h"
#include <Eigen>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Sparse>

namespace MathLib
{
class EigenTools
{
public:
static MathLib::CRSSigned* buildCRSMatrixFromEigenMatrix(Eigen::DynamicSparseMatrix<double, Eigen::RowMajor> &A)
{
    const size_t dimension = A.rows();

    size_t nonzero = 0;
    for (size_t i=0; i<A.outerSize(); ++i)
      for (Eigen::DynamicSparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A,i); it; ++it)
        nonzero++;
    //for (size_t i=0; i<A.rows(); i++) {
    //    for (size_t j=0; j<A.cols(); j++) {
    //        if (A.coeff(i,j)!=.0)
    //            nonzero++;
    //    }
    //}

    double *crs_data = new double [nonzero];
    size_t temp_cnt = 0;
    for (size_t i=0; i<A.outerSize(); ++i)
      for (Eigen::DynamicSparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A,i); it; ++it)
        crs_data[temp_cnt++] = it.value();
    //for (size_t i=0; i<A.rows(); i++) {
    //    for (size_t j=0; j<A.cols(); j++) {
    //        if (A.coeff(i,j)!=.0)
    //            crs_data[temp_cnt++] = A.coeff(i,j);
    //    }
    //}
    assert(temp_cnt==nonzero);
    int *ptr = new int[A.rows()+1];
    int *col_idx = new int[nonzero];
    long counter_ptr = 0, counter_col_idx = 0;
    long cnt_row = 0;

    for (size_t i=0; i<A.outerSize(); ++i)
    {
      ptr[cnt_row++] = counter_ptr;         // starting point of the row
      for (Eigen::DynamicSparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A,i); it; ++it)
      {
        col_idx[counter_col_idx] = it.col();
        ++counter_ptr;
        ++counter_col_idx;
      }
    }
    //for (size_t I=0; I<A.rows(); I++)
    //{
    //    ptr[cnt_row++] = counter_ptr;         // starting point of the row
    //    for (size_t J=0; J<A.cols(); J++)
    //    {
    //        if (A.coeff(I,J)==.0) continue;
    //        long K = counter_ptr;                // index in entry
    //        col_idx[counter_col_idx] = J;
    //        ++counter_ptr;
    //        ++counter_col_idx;
    //    }
    //}
    ptr[A.rows()] = counter_ptr;
#if 0
    //output CRS
    cout << "PTR:" << endl;
    for (size_t i=0; i<A.rows()+1; i++)
        cout << ptr[i] << ", ";
    cout << endl;
    cout << "ColID:" << endl;
    for (size_t i=0; i<nonzero; i++)
        cout << col_idx[i] << ", ";
    cout << endl;
    cout << "Data:" << endl;
    for (size_t i=0; i<nonzero; i++)
        cout << crs_data[i] << ", ";
    cout << endl;
#endif

    MathLib::CRSSigned *crs(new MathLib::CRSSigned);
    crs->dimension = dimension;
    crs->row_ptr = ptr;
    crs->col_idx = col_idx;
    crs->data = crs_data;

    return crs;
};


template<class EigenSpMatrix>
void outputEQS(const std::string &fileName, EigenSpMatrix &A, double *x, double *b)
{
    std::ofstream of(fileName.c_str());
    for (size_t i=0; i<A.rows(); i++) {
        for (size_t j=0; j<A.rows(); j++)
            of << A.coeff(i,j) << "\t";
        of << "\t" << x[i] << "\t" << b[i];
        of << std::endl;
    }

    of.close();
};

template<class EigenSpMatrix>
void setKnownXi(EigenSpMatrix &eqsA, double *eqsRHS, size_t id, double x)
{
    const size_t n_cols = eqsA.cols();
    //A(k, j) = 0.
    //for (size_t j=0; j<n_cols; j++)
    //    if (eqsA.coeff(id, j)!=.0)
    //        eqsA.coeffRef(id, j) = .0;
    for (typename EigenSpMatrix::InnerIterator it(eqsA,id); it; ++it)
    {
      eqsA.coeffRef(id, it.col()) = .0;
    }
    //A(k, k) = val,
    eqsA.coeffRef(id, id) = x;
    //b_i -= A(i,k)*val, i!=k
    //for (size_t j=0; j<eqsA.cols(); j++)
    //    eqsRHS[j] -= eqsA.coeff(j, id)*x;
    for (size_t i=0; i<eqsA.outerSize(); ++i)
    {
      if (i==id) continue;
      double v_i_k = .0;
      for (typename EigenSpMatrix::InnerIterator it(eqsA,i); it; ++it)
      {
        if (it.col()==id) {
          v_i_k = eqsA.coeff(i, id);
          eqsA.coeffRef(i, id) = .0;
          break;
        } else if (it.col()>id) {
          break;
        }
      }
      eqsRHS[i] -= v_i_k*x;
    }
    //b_k = A_kk*val
    eqsRHS[id] = eqsA.coeff(id, id)*x;
    ////A(i, k) = 0., i!=k
    //for (size_t j=0; j<eqsA.cols(); j++)
    //    if (eqsA.coeff(j, id)!=.0 && j!=id)
    //        eqsA.coeffRef(j, id) = .0;
}

} //namespace
} //namespace

#endif

