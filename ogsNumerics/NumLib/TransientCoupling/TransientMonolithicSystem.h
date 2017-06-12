/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientMonolithicSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cassert>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

class AbstractTransientMonolithicSystem : public NumLib::AbstractMonolithicSystem<ITransientCoupledSystem>
{
public:
    virtual ~AbstractTransientMonolithicSystem() {};
};


}
