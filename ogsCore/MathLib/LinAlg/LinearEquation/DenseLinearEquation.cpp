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

#include "DenseLinearEquation.h"

#include <algorithm>

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"


namespace MathLib
{

void DenseLinearEquation::setOption(const BaseLib::Options&)
{

}

void DenseLinearEquation::solve()
{
    this->applyKnownX();

    MathLib::GaussAlgorithm solver(*this->getA());
    AbstractDenseLinearEquation::VectorType *b = this->getRHSAsStdVec();
    AbstractDenseLinearEquation::VectorType *x = this->getXAsStdVec();
    //std::copy(b->begin(), b->end(), x->begin());
    *x = *b;
    solver.execute(&(*x)[0]);
}


} //end namespace
