/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DecomposedVector.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "SubVector.h"

namespace DiscreteLib
{

/**
 * \brief Discrete vector container for OpenMP
 *
 * This vector container utilizes shared memory systems and makes local vectors keep actual memory
 */
template<typename T>
class DecomposedVector : public IDiscreteVector<T>
{
private:
    /// @param n_global the size of vector
    /// @param n_divide the number of local vectors
    DecomposedVector(IDiscreteSystem* sys, MeshLib::IMesh* msh, size_t n_global, std::vector<size_t> &list_range_start)
    {
        _sys = sys;
        _msh = msh;
        _global_n = n_global;
        _local_v.resize(list_range_start.size());
        for (size_t i=0; i<list_range_start.size(); i++) {
            size_t i_begin = list_range_start[i];
            size_t i_end = i+1<list_range_start.size() ? list_range_start[i+1] : n_global;
            _local_v[i] = new SubVector<T>(sys, i, _global_n, i_begin, i_end);
        }
    };

    DecomposedVector(IDiscreteSystem* sys, MeshLib::IMesh* msh, size_t n_global, const std::vector<SubVector<T>*> &local_v)
    {
        _sys = sys;
        _msh = msh;
        _global_n = n_global;
        _local_v.resize(local_v.size());
        for (size_t i=0; i<local_v.size(); i++) {
            _local_v[i] = local_v[i];
        }
    };


public:
    static DecomposedVector<T>* createInstance(IDiscreteSystem &sys, size_t n_global, std::vector<size_t> &list_range_start)
    {
        DecomposedVector<T>* v = new DecomposedVector<T>(&sys, sys.getMesh(), n_global, list_range_start);
        sys.addVector(v);
        return v;
    }

    virtual DecomposedVector<T>* clone() const
    {
        DecomposedVector<T>* v = new DecomposedVector<T>((IDiscreteSystem*)_sys, _sys->getMesh(), _global_n, _local_v);
        _sys->addVector(v);
        return v;
    }

    /// destructor
    virtual ~DecomposedVector()
    {
        BaseLib::releaseObjectsInStdVector(_local_v);
    };

    /// get the size of this vector
    size_t size() const { return _global_n; }

    /// access data
    T& operator[] (size_t idx) 
    {
        return local(idx).global(idx);
    }

    /// access data
    const T& operator[] (size_t idx) const
    {
        return local(idx).global(idx);
    }

    /// get the range begin
    size_t getRangeBegin() const {return 0;};
    /// get the range end
    size_t getRangeEnd() const {return _global_n;};

    /// This function should be called whenever local vectors are updated. It update ghost element from other local vectors.
    void finishUpdate()
    {
        for (size_t i=0; i<_local_v.size(); i++) {
            SubVector<T> *v = _local_v[i];
            for (size_t j=0; j<v->getNumberOfGhostElements(); i++) {
                size_t id = v->getGhostElementId(j);
                v->global(id) = (*this)[id];
            }
        }
    }
    
    /// construct
    virtual void construct(const DiscreteLib::DofEquationIdTable &dofManager, IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, dofManager, *this);
    }


private:
    /// return a local vector which contain the given element
    SubVector<T>& local(size_t idx)
    {
        return *_local_v[localId(idx)];
    };

    /// return a local vector which contain the given element
    const SubVector<T>& local(size_t idx) const
    {
        return *_local_v[localId(idx)];
    };

    /// return an index of a local vector which contain the given element
    size_t localId(size_t idx) const
    {
        size_t pos = 0;
        for (size_t i=0; i<_local_v.size(); ++i) {
            pos = i;
            if (idx < _local_v[i]->getRangeBegin()) {
                pos = i - 1;
                break;
            }
        }
        return pos;
    }

private:
    IDiscreteSystem* _sys;
    MeshLib::IMesh* _msh;
    /// vector length
    size_t _global_n;
    /// a list of local vectors
    std::vector<SubVector<T>*> _local_v;
};


} //end

