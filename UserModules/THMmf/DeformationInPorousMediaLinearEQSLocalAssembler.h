/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DeformationInPorousMediaLinearEQSLocalAssembler.h
 *
 * Created on 2012-11-15 by Norihiro Watanabe
 */


#pragma once

#include "MeshLib/Core/CoordinateSystem.h"
#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"
#include "FemLib/Tools/IFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTransientCoupledLinearEQSLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"

namespace THMmf
{

class DeformationInPorousMediaLinearEQSLocalAssembler
: public NumLib::ElementWiseTransientCoupledLinearEQSLocalAssembler
{
public:

    DeformationInPorousMediaLinearEQSLocalAssembler(
                const FemLib::IFeObjectContainer* feObjects,
                const size_t n_var,
                const std::vector<size_t> &vec_order,
                const MeshLib::CoordinateSystem &problem_coordinates
                )
    : NumLib::ElementWiseTransientCoupledLinearEQSLocalAssembler(n_var, vec_order),
      _feObjects(feObjects->clone()), _problem_coordinates(problem_coordinates)
    {
    };

    virtual ~DeformationInPorousMediaLinearEQSLocalAssembler()
    {
        BaseLib::releaseObject(_feObjects);
    };

protected:
    virtual void assembleComponents(
                const NumLib::TimeStep &/*time*/,
                const MeshLib::IElement &e,
                const std::vector<size_t> &vec_order,
                const std::vector<MathLib::LocalVector> &vec_u0,
                const std::vector<MathLib::LocalVector> &vec_u1,
                std::vector<std::vector<MathLib::LocalMatrix> > &vec_K,
                std::vector<MathLib::LocalVector> &vec_F
                );

private:
    FemLib::IFeObjectContainer* _feObjects;
    const MeshLib::CoordinateSystem _problem_coordinates;
};

}
