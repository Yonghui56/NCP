/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DummyElementWiseTransientResidualLocalAssembler.h
 *
 * Created on 2012-10-24 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"
#include "IElementWiseTransientResidualLocalAssembler.h"

namespace NumLib
{

/**
 * \brief Interface class for element local assembler for transient case
 */
class DummyElementWiseTransientResidualLocalAssembler : public IElementWiseTransientResidualLocalAssembler
{
public:

    virtual ~DummyElementWiseTransientResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_r        local residual
    virtual void assembly(
            const TimeStep &/*time*/,
            const MeshLib::IElement &/*e*/,
            const DiscreteLib::DofEquationIdTable &/*localDofManager*/,
            const MathLib::LocalVector &/*local_u_n1*/,
            const MathLib::LocalVector &/*local_u_n*/,
            MathLib::LocalVector &/*local_r*/) {};
};

} //end
