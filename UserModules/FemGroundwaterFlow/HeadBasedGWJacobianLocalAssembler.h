/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GWJacobianLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"


class HeadBasedGWJacobianLocalAssembler
: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    explicit HeadBasedGWJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects)
    : _feObjects(feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,  MathLib::LocalMatrix &localJ);

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};


