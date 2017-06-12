/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SubDomain.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "SubDomain.h"

namespace DiscreteLib
{

SubDomain::SubDomain(MeshLib::IMesh* msh, IGlobaLocalMappingTable* mapping, std::set<size_t>* list_ghosts)
: _local_msh(msh), _map_global_local(mapping)
{
    _dom_id = 0;
    if (list_ghosts) {
        _list_ghosts.insert(list_ghosts->begin(), list_ghosts->end());
    }
}

} //end

