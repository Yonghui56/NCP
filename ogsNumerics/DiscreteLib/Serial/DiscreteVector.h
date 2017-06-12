/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteVector.h
 *
 * Created on 2012-08-20 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <cstddef>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 * \brief Vector container for single memory
 */
template<typename T>
class DiscreteVector : public IDiscreteVector<T>
{
protected:
    explicit DiscreteVector(MeshLib::IMesh* msh) : _data(0), _msh(msh), _sys(0) {};
    explicit DiscreteVector(size_t n, IDiscreteSystem* sys) : _data(n), _msh(sys->getMesh()), _sys(sys) {};
    DiscreteVector(MeshLib::IMesh* msh, size_t n) : _data(n), _msh(msh), _sys(0) {};

public:
    static DiscreteVector<T>* createInstance(IDiscreteSystem &sys, size_t n)
    {
        DiscreteVector<T>* vec = new DiscreteVector<T>(n, &sys);
        sys.addVector(vec);
        return vec;
    }
    virtual ~DiscreteVector() {};

    virtual DiscreteVector<T>* clone() const
    {
        DiscreteVector<T>* vec = DiscreteVector<T>::createInstance(*_sys, _data.size());
        *vec = (*this);
        return vec;
    }

    virtual void resize(size_t n) {_data.resize(n);};
    virtual size_t size() const {return _data.size();};

    virtual DiscreteVector<T>& operator= (T v)
    {
        for (size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }

    virtual T& operator[] (size_t idx) 
    {
        return _data[idx];
    }
    virtual const T& operator[] (size_t idx) const
    {
        return _data[idx];
    }

    typename std::vector<T>::iterator begin() 
    {
        return _data.begin();
    }
    typename std::vector<T>::iterator end() 
    {
        return _data.end();
    }
    virtual size_t getRangeBegin() const
    {
        return 0;
    }
    virtual size_t getRangeEnd() const
    {
        return _data.size();
    }

    /// construct
    virtual void construct(const DofEquationIdTable &dofEquationIdTable, IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, dofEquationIdTable, *this);
    }

protected:
    std::vector<T> _data;
    MeshLib::IMesh* _msh;
    IDiscreteSystem* _sys;
};

} // end
