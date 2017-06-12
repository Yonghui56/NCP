/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparseTableCRS.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <cassert>

#include "Sparsity.h"

namespace MathLib 
{

/**
 * \brief Compressible row storage
 */
template <class T>
struct SparseTableCRS 
{
    size_t dimension; 
    size_t nonzero;
    T* row_ptr;
    T* col_idx;
    double* data;

    SparseTableCRS()
    : dimension(0), nonzero(0), row_ptr(NULL), col_idx(NULL), data(NULL)
    {};
};

typedef struct SparseTableCRS<unsigned> CRSUnsigned;
typedef struct SparseTableCRS<signed> CRSSigned;

/**
 * convert row-major sparsity to CRS data
 */
template<class INTTYPE>
void convertRowMajorSparsityToCRS(const RowMajorSparsity &row_major_entries, MathLib::SparseTableCRS<INTTYPE> &crs)
{
    const size_t n_rows = row_major_entries.size();
    assert(n_rows > 0);
    if (n_rows==0) return;
         

    //get number of nonzero 
    INTTYPE* ptr = new INTTYPE[n_rows+1];
    std::vector<INTTYPE> vec_col_idx;
    size_t counter_ptr = 0;
    size_t cnt_row = 0;

    for (size_t i=0; i<n_rows; i++) {
        ptr[cnt_row++] = counter_ptr;         // starting point of the row

        // entries at the i th row
        const std::set<size_t> &setConnection = row_major_entries[i];
        //
        for (std::set<size_t>::iterator it=setConnection.begin(); it!=setConnection.end(); it++) {
            vec_col_idx.push_back(static_cast<INTTYPE>(*it));
            ++counter_ptr;
        }
    }

    ptr[n_rows] = counter_ptr;

    crs.dimension = n_rows;
    crs.row_ptr = ptr;
    crs.col_idx = new INTTYPE[vec_col_idx.size()];
    for (size_t i=0; i<vec_col_idx.size(); i++)
        crs.col_idx[i] = vec_col_idx[i];
    crs.nonzero = vec_col_idx.size();
    crs.data = new double[vec_col_idx.size()];
    for (size_t i=0; i<vec_col_idx.size(); i++)
        crs.data[i] = .0;
}

template<class INTTYPE>
void outputSparseTableCRS(MathLib::SparseTableCRS<INTTYPE> *crs)
{
    //output CRS
    std::cout << "PTR:" << std::endl;
    for (size_t i=0; i<crs->dimension+1; i++)
        std::cout << crs->row_ptr[i] << ", "; 
    std::cout << std::endl;
    std::cout << "ColID:" << std::endl;
    for (size_t i=0; i<crs->nonzero; i++)
        std::cout << crs->col_idx[i] << ", "; 
    std::cout << std::endl;
    std::cout << "Data:" << std::endl;
    for (size_t i=0; i<crs->nonzero; i++)
        std::cout << crs->data[i] << ", "; 
    std::cout << std::endl;
}

} // end namespace MathLib

