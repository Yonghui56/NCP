/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTransientCoupledResidualLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "MathLib/DataType.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientResidualLocalAssembler.h"

namespace NumLib
{

class ElementWiseTransientCoupledResidualLocalAssembler
: public IElementWiseTransientResidualLocalAssembler
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    ElementWiseTransientCoupledResidualLocalAssembler(size_t n_var, const std::vector<size_t> &vec_order)
    : _n_var(n_var), _vec_order(vec_order)
    {
    };

    virtual ~ElementWiseTransientCoupledResidualLocalAssembler() {};

    virtual void assembly(  const NumLib::TimeStep &timestep,
                            const MeshLib::IElement &e,
                            const DiscreteLib::DofEquationIdTable &localDofManager,
                            const LocalVectorType &local_u_n1,
                            const LocalVectorType &local_u_n,
                            LocalVectorType &local_r
                            );


protected:
    virtual void assembleComponents(  const NumLib::TimeStep &/*time*/,
                            const MeshLib::IElement &e,
                            const std::vector<size_t> &vec_order,
                            const std::vector<LocalVectorType> &vec_u0,
                            const std::vector<LocalVectorType> &vec_u1,
                            std::vector<LocalVectorType> &vec_F
                            ) = 0;

private:
    size_t _n_var;
    std::vector<size_t> _vec_order;
};

} //end
