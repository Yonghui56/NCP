/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SequentialGlobaLocalMappingTable.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IGlobaLocalMappingTable.h"

namespace DiscreteLib
{

/**
 * \brief Mapping table for sequential rule
 * 
 * global id = local id + offset
 */
class SequentialGlobaLocalMappingTable : public IGlobaLocalMappingTable
{
public:
    /**
     * 
     * @param i_start   Valid global id range begin 
     * @param i_end     Valid global id range end
     * @param offset    offset
     */
    SequentialGlobaLocalMappingTable(size_t i_start, size_t i_end, size_t offset)
    : _i_start(i_start), _i_end(i_end), _offset(offset)
    { };

    ///
    virtual ~SequentialGlobaLocalMappingTable() {};

    /**
     * 
     * @param global_id
     * @return
     */
    virtual bool hasGlobal(size_t global_id) 
    {
        return _i_start<=global_id && global_id<_i_end;
    }

    /**
     * 
     * @param global
     * @return
     */
    virtual size_t global2local(size_t global)
    {
        return global - _offset;
    }

    /**
     * 
     * @param local
     * @return
     */
    virtual size_t local2global(size_t local)
    {
        return local + _offset;
    }
    
private:
    size_t _i_start;
    size_t _i_end;
    size_t _offset;
};

} //end

