/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IElementWiseTimeODELocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"


namespace MeshLib
{
class IElement;
};

namespace NumLib
{

class TimeStep;

/**
 * \brief Interface class of element assembler for time ODE formulations
 */
class IElementWiseTimeODELocalAssembler
{
public:
    virtual ~IElementWiseTimeODELocalAssembler() {};

    /// assemble components of a local linear equation for the given element
    /// @param time        time step
    /// @param e        element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param M        mass matrix
    /// @param K        laplace matrix
    /// @param F        source/sink terms
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalMatrix &M, MathLib::LocalMatrix &K, MathLib::LocalVector &F)  = 0;
};

} //end
