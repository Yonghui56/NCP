/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteDataContainer.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"

namespace DiscreteLib
{

/**
 * \brief Data container for any discrete objects
 */
class DiscreteDataContainer
{
public:
    DiscreteDataContainer() {};
    virtual ~DiscreteDataContainer()
    {
        BaseLib::releaseObjectsInStdVector(_vec_vectors);
        BaseLib::releaseObjectsInStdVector(_vec_linear_sys);
    }

    size_t addVector(IDiscreteVectorBase* v)
    {
        _vec_vectors.push_back(v);
        v->setObjectID(_vec_vectors.size()-1);
        return _vec_vectors.size()-1;
    }

    void eraseVector(IDiscreteVectorBase* v)
    {
        const size_t i = v->getObjectID();
        if (_vec_vectors.size() > i) {
            _vec_vectors[i] = 0;
        }
    }

    size_t getNumberOfVectors() const {return _vec_vectors.size();};

    IDiscreteVectorBase* getVector(size_t i)
    {
        return _vec_vectors[i];
    }

    size_t addLinearEquation(IDiscreteLinearEquation* eq)
    {
        _vec_linear_sys.push_back(eq);
        eq->setObjectID(_vec_linear_sys.size()-1);
        return _vec_linear_sys.size()-1;
    }

    void eraseLinearEquation(IDiscreteLinearEquation* eq)
    {
        const size_t i = eq->getObjectID();
        if (_vec_linear_sys.size() > i) {
            _vec_linear_sys[i] = 0;
        }
    }

    size_t getNumberOfLinearEquations() const {return _vec_linear_sys.size();};

    IDiscreteLinearEquation* getLinearEquation(size_t i)
    {
        return _vec_linear_sys[i];
    }

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteDataContainer);

    std::vector<IDiscreteLinearEquation*> _vec_linear_sys;
    std::vector<IDiscreteVectorBase*> _vec_vectors;
};

} //end
