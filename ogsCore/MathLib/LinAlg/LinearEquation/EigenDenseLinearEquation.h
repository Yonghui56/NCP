/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EigenDenseLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <valarray>
#include <vector>
#include <Eigen>

#include "ILinearEquation.h"


namespace MathLib
{

/**
 * \brief Dense linear equation using Eigen
 */
class EigenDenseLinearEquation : public ILinearEquation
{
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
    typedef Eigen::VectorXd VectorType;

    EigenDenseLinearEquation() {};

    virtual ~EigenDenseLinearEquation() {};

    void create(size_t length, RowMajorSparsity* /*sp*/=0);

    bool isCreated() const { return _A.rows()>0; };

    void resize(size_t length);

    void reset();

    size_t getDimension() const { return _A.rows(); }

    MatrixType* getA() { return &_A; }

    double getA(size_t rowId, size_t colId) const
    {
        return _A(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        _A(rowId, colId) = v;
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        _A(rowId, colId) += v;
    }

    double getRHS(size_t rowId) const
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    VectorType* getRHSAsVec() {return &_b;};

    void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    void addRHS(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    double* getX()
    {
        return &_x[0];
    }

    VectorType* getXAsVec() {return &_x;};


    void setKnownX(size_t row_id, double x);

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x);

    void printout(std::ostream &os=std::cout) const;

    void setOption(const BaseLib::Options &/*option*/) {};

    void solve();

private:
    void applyKnownX();

private:
    MatrixType _A;
    VectorType _b;
    VectorType _x;
    std::vector<size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;

};

}
