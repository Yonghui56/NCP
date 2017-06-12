/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIC.cpp
 *
 * Created on 2012-09-22 by Norihiro Watanabe
 */

#include "FemIC.h"

#include "logog.hpp"
#include "FemLib/BC/IC2FEM.h"

namespace SolutionLib
{

///
void FemIC::addDistribution(const GeoLib::GeoObject* geo, const NumLib::ITXFunction* ic_func)
{
    _vec_geo.push_back(geo);
    _vec_func.push_back(ic_func);
}

/// setup 
void FemIC::setup(NumLib::ITXDiscreteFunction<double> &u0) const
{
    DiscreteLib::IDiscreteVector<double> *u0_array = u0.getDiscreteData();

    if (_vec_geo.size()==0)
        WARN("***WARN: IC not found.");

    for (size_t i=0; i<_vec_geo.size(); i++) {
        std::vector<size_t> vec_node_id;
        std::vector<double> vec_node_value;
        FemLib::IC2FEM ic2fem(*_msh, *_vec_geo[i], *_vec_func[i], vec_node_id, vec_node_value);


        for (size_t j=0; j<vec_node_id.size(); j++) {
            (*u0_array)[vec_node_id[j]] = vec_node_value[j]; 
        }
    }
    
}

}
