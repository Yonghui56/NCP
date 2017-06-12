/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemPoroelasticResidualLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/ElementWiseTransientCoupledResidualLocalAssembler.h"

class FemPoroelasticResidualLocalAssembler : public NumLib::ElementWiseTransientCoupledResidualLocalAssembler
{
public:
    FemPoroelasticResidualLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, size_t n_var, const std::vector<size_t> &vec_order)
    : NumLib::ElementWiseTransientCoupledResidualLocalAssembler(n_var, vec_order), _feObjects(*feObjects)
    {
    };

    virtual ~FemPoroelasticResidualLocalAssembler() {};

protected:
    virtual void assembleComponents(  const NumLib::TimeStep &/*time*/,  
                            const MeshLib::IElement &e, 
                            const std::vector<size_t> &vec_order, 
                            const std::vector<LocalVectorType> &vec_u0, 
                            const std::vector<LocalVectorType> &vec_u1, 
                            std::vector<LocalVectorType> &vec_F
                            );

private:
    FemLib::LagrangeFeObjectContainer _feObjects;
};
