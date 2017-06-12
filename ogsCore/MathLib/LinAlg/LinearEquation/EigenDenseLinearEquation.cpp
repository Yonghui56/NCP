/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EigenDenseLinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "EigenDenseLinearEquation.h"


namespace MathLib
{

void EigenDenseLinearEquation::create(size_t length, RowMajorSparsity* /*sp*/)
{
    resize(length);
}

void EigenDenseLinearEquation::resize(size_t length)
{
    _A = MatrixType::Zero(length, length);
    _b = VectorType::Zero(length);
    _x = VectorType::Zero(length);
    //reset();
}

void EigenDenseLinearEquation::reset()
{
    _A *= .0;
    _b *= .0;
    _x *= .0;
}

void EigenDenseLinearEquation::applyKnownX()
{
    const size_t n_cols = _A.cols();
    for (std::size_t i=0; i<_vec_knownX_id.size(); i++) {
        std::size_t row_id = _vec_knownX_id[i];
        double x = _vec_knownX_x[i];
        //A(k, j) = 0.
        for (size_t j=0; j<n_cols; j++)
            _A(row_id, j) = .0;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<n_cols; j++)
            _b[j] -= _A(j, row_id)*x;
        //b_k = val
        _b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<n_cols; j++)
            _A(j, row_id) = .0;
        //A(k, k) = 1.0
        _A(row_id, row_id) = 1.0; //=x
    }
}

void EigenDenseLinearEquation::setKnownX(size_t row_id, double x)
{
    _vec_knownX_id.push_back(row_id);
    _vec_knownX_x.push_back(x);
}

void EigenDenseLinearEquation::setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
{
    _vec_knownX_id.assign(vec_id.begin(), vec_id.end());
    _vec_knownX_x.assign(vec_x.begin(), vec_x.end());
}

void EigenDenseLinearEquation::printout(std::ostream &os) const
{
    os << "A=" << _A << std::endl;
    os << "x=" << _x << std::endl;
    os << "b=" << _b << std::endl;
}

void setOption(const BaseLib::Options &/*option*/) {};

void EigenDenseLinearEquation::solve()
{
    this->applyKnownX();
    _x = _A.colPivHouseholderQr().solve(_b);
}

}
