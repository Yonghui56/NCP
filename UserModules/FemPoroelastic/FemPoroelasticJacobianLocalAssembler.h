/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemPoroelasticJacobianLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/ElementWiseTransientCoupledJacobianLocalAssembler.h"

class FemPoroelasticJacobianLocalAssembler
: public NumLib::ElementWiseTransientCoupledJacobianLocalAssembler
{
public:
    FemPoroelasticJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, size_t n_var, const std::vector<size_t> &vec_order)
    : NumLib::ElementWiseTransientCoupledJacobianLocalAssembler(n_var, vec_order), _feObjects(*feObjects)
    {
    };

    virtual ~FemPoroelasticJacobianLocalAssembler() {};

    virtual void assembleComponents(  const NumLib::TimeStep &/*time*/,
                            const MeshLib::IElement &e,
                            const std::vector<size_t> &vec_order,
                            const std::vector<LocalVectorType> &vec_u0,
                            const std::vector<LocalVectorType> &vec_u1,
                            std::vector<std::vector<LocalMatrixType> > &vec_K
                            );

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};
