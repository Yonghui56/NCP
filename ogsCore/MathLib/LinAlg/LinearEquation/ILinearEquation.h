/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ILinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <iostream>
#include <vector>
#include "Eigen"
#include "BaseLib/Options.h"
#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Sparse/Sparsity.h"

namespace MathLib
{

/**
 * \brief Interface to all algebraic linear equation systems Ax=b
 *
 */
class ILinearEquation
{
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
    typedef Eigen::VectorXd LocalVector;

    /// 
    virtual ~ILinearEquation() {};

    /**
     * create a linear equation 
     *
     * \param dimension    dimension of the equation
     * \param sparsity     sparcity pattern
     */
    virtual void create(size_t dimension, RowMajorSparsity *sparsity=0) = 0;

    /**
     * return if this equation is already created
     */
    virtual bool isCreated() const = 0;

    /**
     * set properties
     */
    virtual void setOption(const BaseLib::Options &option) = 0;

    /**
     * reset the equation
     */
    virtual void reset() = 0;

    /// return dimension of this equation
    virtual size_t getDimension() const = 0;

    /// get an entriy in a matrix
    virtual double getA(size_t rowId, size_t colId) const = 0;

    /// set an entriy in a matrix
    virtual void setA(size_t rowId, size_t colId, double v) = 0;

    /// add a value into a matrix
    virtual void addA(size_t rowId, size_t colId, double v) = 0;

    /// add submatrix into a matrix
    virtual void addAsub(const std::vector<size_t> &vec_row_pos, const std::vector<size_t> &vec_col_pos, const LocalMatrix &sub_matrix, double fkt=1.0);

    /// add submatrix into a matrix
    virtual void addAsub(const std::vector<size_t> &vec_row_col_pos, const LocalMatrix &sub_matrix, double fkt=1.0);

    /// get RHS entriy
    virtual double getRHS(size_t rowId) const = 0;

    /// get RHS vector
    virtual double* getRHS() = 0;

    /// set RHS entry
    virtual void setRHS(size_t rowId, double v) = 0;

    /// add RHS entry
    virtual void addRHS(size_t rowId, double v) = 0;

    /// add a sub vector to RHS
    virtual void addRHSsub(const std::vector<size_t> &vec_row_pos, const double* sub_vector, double fkt=1.0);

    /// add a sub vector to RHS
    virtual void addRHSsub(const std::vector<size_t> &vec_row_pos, const LocalVector &sub_vector, double fkt=1.0);

    /// get a solution vector
    virtual double* getX() = 0;

    /// set prescrived unknown
    virtual void setKnownX(size_t row_id, double x) = 0;

    /// set prescrived unknown
    virtual void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x) = 0;

    /// solve
    virtual void solve() = 0;

    /// print out the equation for debugging
    virtual void printout(std::ostream &os=std::cout) const = 0;
};

}
