/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DummyElementWiseTransientJacobianLocalAssembler.h
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#pragma once

#include "IElementWiseTransientJacobianLocalAssembler.h"

namespace NumLib
{

class DummyElementWiseTransientJacobianLocalAssembler : public IElementWiseTransientJacobianLocalAssembler
{
public:
    virtual ~DummyElementWiseTransientJacobianLocalAssembler() {};

    /// assemble a local Jacobian matrix for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_J        local Jacobian
    virtual void assembly(const TimeStep &,  const MeshLib::IElement &, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector &, const MathLib::LocalVector &, MathLib::LocalMatrix &) {};
};


} //end
