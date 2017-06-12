/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DenseLinearEquations.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "AbstractDenseLinearEquation.h"

#include <algorithm>


namespace MathLib
{

void AbstractDenseLinearEquation::resize(size_t length)
{
    if (_A==0) {
        _A = new Matrix<double>(length, length);
    } else {
        _A->resize(length, length);
    }
    _b.resize(length);
    _x.resize(length);
    reset();
}

void AbstractDenseLinearEquation::applyKnownX()
{
    const size_t n_cols = _A->getNCols();
    for (std::size_t i=0; i<_vec_knownX_id.size(); i++) {
        std::size_t row_id = _vec_knownX_id[i];
        double x = _vec_knownX_x[i];
        //A(k, j) = 0.
        for (size_t j=0; j<n_cols; j++)
            (*_A)(row_id, j) = .0;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<n_cols; j++)
            _b[j] -= (*_A)(j, row_id)*x;
        //b_k = val
        _b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<n_cols; j++)
            (*_A)(j, row_id) = .0;
        //A(k, k) = 1.0
        (*_A)(row_id, row_id) = 1.0; //=x
    }
}

} //end namespace
