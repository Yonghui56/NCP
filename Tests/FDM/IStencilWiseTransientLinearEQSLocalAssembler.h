/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IStencilWiseTransientLinearEQSLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"

namespace NumLib
{

class TimeStep;
}

namespace FdmLib
{
class IStencil;

class IStencilWiseTransientLinearEQSLocalAssembler
{
public:
    typedef MathLib::LocalEquation LocalEquationType;
    typedef LocalEquationType::MatrixType LocalMatrixType;
    typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~IStencilWiseTransientLinearEQSLocalAssembler() {};

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(const NumLib::TimeStep &time,  IStencil &nod, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs) = 0;
};

} //end
