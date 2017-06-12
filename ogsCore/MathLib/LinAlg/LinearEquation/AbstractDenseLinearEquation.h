/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractDenseLinearEquation
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <valarray>
#include <vector>

#include "ILinearEquation.h"
#include "MathLib/LinAlg/Dense/Matrix.h"


namespace MathLib
{

/**
 * \brief Abstract class to any linear equation based on a dense matrix
 */
class AbstractDenseLinearEquation : public ILinearEquation
{
public:
    typedef Matrix<double> MatrixType;
    typedef std::valarray<double> VectorType;

    ///
    AbstractDenseLinearEquation() : _A(NULL) {};

    ///
    virtual ~AbstractDenseLinearEquation()
    {
        if (_A) delete _A;
    }

    ///
    virtual void create(size_t length, RowMajorSparsity* /*sp*/=0)
    {
        resize(length);
    }

    ///
    virtual bool isCreated() const { return _A!=0; };

    ///
    void resize(size_t length);

    ///
    virtual void reset()
    {
        (*_A) = .0;
        _b.resize(_b.size(), .0);
        _x.resize(_x.size(), .0);
    }

    virtual size_t getDimension() const
    {
        return _A->getNRows();
    }

    MatrixType* getA() 
    {
        return _A;
    }

    virtual double getA(size_t rowId, size_t colId) const
    {
        return (*_A)(rowId, colId);
    }

    virtual void setA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) = v;
    }

    virtual void addA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) += v;
    }

    virtual double getRHS(size_t rowId) const
    {
        return _b[rowId];
    }

    virtual double* getRHS()
    {
        return &_b[0];
    }

    VectorType* getRHSAsStdVec() {return &_b;};

    virtual void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    virtual void addRHS(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    virtual double* getX()
    {
        return &_x[0];
    }

    VectorType* getXAsStdVec() {return &_x;};


    virtual void setKnownX(size_t row_id, double x)
    {
        _vec_knownX_id.push_back(row_id);
        _vec_knownX_x.push_back(x);
    }

    virtual void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        _vec_knownX_id.assign(vec_id.begin(), vec_id.end());
        _vec_knownX_x.assign(vec_x.begin(), vec_x.end());
    }

    virtual void printout(std::ostream &os=std::cout) const
    {
        //os << "not implemented yet." << std::endl;
    	os << *_A << std::endl;
    	//return * _A;
    }

protected:
    void applyKnownX();

private:
    MatrixType *_A;
    VectorType _b;
    VectorType _x;
    std::vector<size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;

};

}
