/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SubVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "DiscreteLib/Core/IDiscreteVector.h"

namespace DiscreteLib
{
/**
 * \brief Local discrete vector for OpenMP
 */
template<typename T>
class SubVector : public DiscreteVector<T>
{     
 public:
     /// Constructor without ghost elements
     /// @param local_id    id of this local vector
     /// @param n_global    size of the global vector 
     /// @param i_begin     global index starting this local vector
     /// @param i_end       global index ending this local vector + 1 
    SubVector(IDiscreteSystem *sys, size_t local_id, size_t n_global, size_t i_begin, size_t i_end)
    : DiscreteVector<T>(i_end - i_begin, sys)
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    /// Constructor with ghost elements
    /// @param local_id    id of this local vector
    /// @param n_global    size of the global vector 
    /// @param i_begin     global index starting this local vector
    /// @param i_end       global index ending this local vector + 1 
    /// @param ghost_id    a list of ghost id
    SubVector(size_t local_id, size_t n_global, size_t i_begin, size_t i_end, std::vector<size_t> &ghost_id) : DiscreteVector<T>(i_end - i_begin + ghost_id.size()), _ghost_id(ghost_id)
    {
        _local_id = local_id;
        _global_n = n_global;
        _i_start = i_begin;
        _i_end = i_end;
    };

    /// destructor
    virtual ~SubVector() {};

    /// get the global vector size
    size_t getGlobalSize() const {return _global_n;};

    /// get the range begin
    size_t getRangeBegin() const {return _i_start;};
    /// get the range end
    size_t getRangeEnd() const {return _i_end;};

    /// get the number of ghost elements in this vector
    size_t getNumberOfGhostElements() const {return _ghost_id.size();};
    /// get the ghost element id (global) in this vector
    size_t getGhostElementId(size_t i) const {return _ghost_id[i];};
    /// check if the give id can be a ghost element in this vector. This function doesn't check if actually this vector contain the element.
    bool isGhost(size_t global_idx) const {return !(_i_start<=global_idx && global_idx<_i_end);};

    /// check if this vector contain an element with the given global id
    bool has(size_t global_idx) const
    {
        if (!isGhost(global_idx)) {
            return (_i_start<=global_idx && global_idx<_i_end);
        } else {
            return (_ghost_id.end()!=std::find(_ghost_id.begin(), _ghost_id.end(), global_idx));
        }
    }

    /// access with local index
    T& operator[] (size_t i_local) { return DiscreteVector<T>::_data[i_local]; }
    /// access with local index
    const T& operator[] (size_t i_local) const { return DiscreteVector<T>::_data[i_local]; }

    /// access with global index
    T& global (size_t i) { return DiscreteVector<T>::_data[access(i)]; }
    /// access with global index
   const T& global (size_t i) const { return DiscreteVector<T>::_data[access(i)]; }


private:
    size_t _global_n;
    size_t _i_start;
    size_t _i_end;
    size_t _local_id;
    /// a list of global id of ghost elements
    std::vector<size_t> _ghost_id;

    /// get local id
    inline size_t access(size_t global_idx) const 
    {
        assert (has(global_idx));
        if (!isGhost(global_idx)) {
            return global_idx - _i_start;
        } else {
            return _i_end + (std::find(_ghost_id.begin(), _ghost_id.end(), global_idx) - _ghost_id.begin());
        }
    };

};

} //end

