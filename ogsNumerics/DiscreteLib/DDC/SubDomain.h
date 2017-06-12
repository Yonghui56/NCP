/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SubDomain.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <set>

#include "IGlobaLocalMappingTable.h"

namespace MeshLib
{
class IMesh;
}

namespace DiscreteLib
{

/**
 * \brief Sub-domain class
 *
 * This class owns
 * - sub-domain id
 * - local mesh
 * - mapping table for global and local id
 * - list of ghost objects
 */
class SubDomain
{
public:
    /**
     * 
     * @param msh
     * @param mapping
     * @param list_ghosts
     */
    SubDomain(MeshLib::IMesh* msh, IGlobaLocalMappingTable* mapping, std::set<size_t>* list_ghosts=0);

    /**
     * 
     */
    virtual ~SubDomain() {};

    /**
     * 
     * @param i
     */
    void setDomainID(size_t i) {_dom_id = i;};
    
    /**
     * 
     * @return
     */
    size_t getDomainID() const {return _dom_id;};

    /**
     * 
     * @return
     */
    MeshLib::IMesh* getLoalMesh() {return _local_msh;};

    /**
     * 
     * @return
     */
    IGlobaLocalMappingTable* getGlobalLocalIdMap() {return _map_global_local;};

    /**
     * 
     * @return
     */
    size_t getNumberOfGhosts() const {return _list_ghosts.size();};
    
    /**
     * 
     * @return
     */
    std::set<size_t>* getGhostList() {return &_list_ghosts;};

private:
    size_t _dom_id;
    MeshLib::IMesh* _local_msh;
    IGlobaLocalMappingTable* _map_global_local;
    std::set<size_t> _list_ghosts;
};

} //end

