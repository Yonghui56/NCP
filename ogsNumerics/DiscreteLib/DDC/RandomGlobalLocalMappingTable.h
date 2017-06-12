/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file RandomGlobalLocalMappingTable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/BidirectionalMap.h"
#include "IGlobaLocalMappingTable.h"

namespace DiscreteLib
{

/**
 * \brief Mapping table for random rule
 */
class RandomGlobalLocalMappingTable : public IGlobaLocalMappingTable
{
public:
    /**
     * 
     * @param map
     */
    explicit RandomGlobalLocalMappingTable(BaseLib::BidirectionalMap<size_t, size_t>* map)
    : _map_global_local(map)
    {
    }
    
    /**
     * 
     */
    virtual ~RandomGlobalLocalMappingTable()
    {
        BaseLib::releaseObject(_map_global_local);
    }
    
    /**
     * 
     * @param global_id
     * @return
     */
    virtual bool hasGlobal(size_t global_id) 
    {
        return _map_global_local->countInA(global_id)>0;
    }
    
    /**
     * 
     * @param global
     * @return
     */
    virtual size_t global2local(size_t global)
    {
        return _map_global_local->mapAtoB(global);
    }
    
    /**
     * 
     * @param local
     * @return
     */
    virtual size_t local2global(size_t local)
    {
        return _map_global_local->mapBtoA(local);
    }
private:
    BaseLib::BidirectionalMap<size_t, size_t>* _map_global_local;
};

} //end

