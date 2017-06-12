/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Utils/DiscreteDataContainer.h"

namespace DiscreteLib
{
class IDiscreteLinearEquation;
class IDiscreteVectorBase;

/**
 * \brief Interface for all kinds of discrete systems
 *
 *  Discrete system class contains the followings
 *  - discrete space (i.e. mesh)
 *  - discrete data (i.e. vector)
 *  - linear equations which are used to calculate discrete data
 */
class IDiscreteSystem
{
public:
    virtual ~IDiscreteSystem() {};
    
    /// return this mesh
    virtual MeshLib::IMesh* getMesh() const = 0;

    /// add equation object
    void addLinearEquation(IDiscreteLinearEquation *eqs)
    {
        _data.addLinearEquation(eqs);
    }

    /// delete equation object
    void deleteLinearEquation(IDiscreteLinearEquation* eqs)
    {
        if (eqs!=0) {
            _data.eraseLinearEquation(eqs);
            delete eqs;
        }
    }

    /// add vector object
    void addVector(IDiscreteVectorBase *vec)
    {
        _data.addVector(vec);
    }

    /// delete this vector object
    /// \param vector object
    void deleteVector(IDiscreteVectorBase* v)
    {
        if (v!=0) {
            _data.eraseVector(v);
            delete v;
        }
    };

private:
    DiscreteDataContainer _data;
};

} //end
