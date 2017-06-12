/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticResidualLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientResidualLocalAssembler.h"

class FemLinearElasticResidualLocalAssembler : public NumLib::IElementWiseTransientResidualLocalAssembler
{
public:
    explicit FemLinearElasticResidualLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects)
    : _feObjects(*feObjects)
    {
    };

    virtual ~FemLinearElasticResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_r        local residual
    virtual void assembly(const NumLib::TimeStep &time,  const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalVector &local_r);

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};
