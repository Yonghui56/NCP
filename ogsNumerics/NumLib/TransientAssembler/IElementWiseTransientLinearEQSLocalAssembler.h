/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IElementWiseTransientLinearEQSLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"


namespace MeshLib
{
class IElement;
};

namespace DiscreteLib
{
class DofEquationIdTable;
}

namespace NumLib
{

class TimeStep;

/**
 * \brief Interface class for element local assembler for transient case
 */
class IElementWiseTransientLinearEQSLocalAssembler
{
public:
    virtual ~IElementWiseTransientLinearEQSLocalAssembler() {};

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(  const TimeStep &time,  const MeshLib::IElement &e, 
                            const DiscreteLib::DofEquationIdTable &localDofManager, 
                            const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, 
                            MathLib::LocalEquation &eqs) = 0;

};

} //end
