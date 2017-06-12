/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DecomposedDomain.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "DecomposedDomain.h"

#include "MeshLib/Core/IMesh.h"
#include "SubDomain.h"

namespace DiscreteLib
{

size_t DecomposedDomain::addSubDomain(SubDomain* sub)
{
    _list_dom.push_back(sub);
    sub->setDomainID(_list_dom.size()-1);
    if (_ddc_type==DecompositionType::Node) {
        _n_discrete_pt += sub->getLoalMesh()->getNumberOfNodes();
    } else {
        _n_discrete_pt += sub->getLoalMesh()->getNumberOfElements();
    }
    _n_discrete_pt -= sub->getNumberOfGhosts();

    _n_nodes += sub->getLoalMesh()->getNumberOfNodes();
    _n_eles += sub->getLoalMesh()->getNumberOfElements();
    if (_ddc_type==DecompositionType::Node) {
        _n_nodes -= sub->getNumberOfGhosts();
    } else {
        _n_eles -= sub->getNumberOfGhosts();
    }
    
    _id = sub->getLoalMesh()->getID();

    return sub->getDomainID();
}

size_t DecomposedDomain::findSubDomainID(size_t global_obj_id) const
{
    for (size_t i=0; i<_list_dom_start_obj.size(); ++i) {
        if (_list_dom_start_obj[i] <= global_obj_id)
            return i;
    }
    return 0; //error
}

} //end

