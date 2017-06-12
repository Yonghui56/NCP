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

#include "IElementWiseTransientLinearEQSLocalAssembler.h"

namespace NumLib
{

class DummyElementWiseTransientLinearEQSLocalAssembler : public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    virtual ~DummyElementWiseTransientLinearEQSLocalAssembler() {};

    virtual void assembly(const TimeStep &,  const MeshLib::IElement &, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector &, const MathLib::LocalVector &, MathLib::LocalEquation &) {};
};


} //end
