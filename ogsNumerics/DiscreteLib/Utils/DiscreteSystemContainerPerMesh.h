/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteSystemContainerPerMesh.h
 *
 * Created on 2012-08-17 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/DiscreteEnums.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"

namespace DiscreteLib
{

/**
 * \brief
 */
class DiscreteSystemContainerPerMesh
{
private:
    DiscreteSystemContainerPerMesh()
    {
    }

    static DiscreteSystemContainerPerMesh* _obj;

    typedef std::pair<size_t, DiscreteSystemType::type> MyKey;
public:
    static DiscreteSystemContainerPerMesh* getInstance();

    ~DiscreteSystemContainerPerMesh()
    {
        BaseLib::releaseObjectsInStdMap(_map_sys);
    }

    template <class T>
    T* createObject(const MeshLib::IMesh* msh)
    {
        MyKey key(msh->getID(), T::getSystemType());
        T* sys = 0;
        std::map<MyKey, IDiscreteSystem*>::iterator itr = _map_sys.find(key);
        if (itr == _map_sys.end()) {
            sys = new T(msh);
            _map_sys.insert(std::pair<MyKey, IDiscreteSystem*>(key, sys));
        } else {
            sys = static_cast<T*>(itr->second);
        }
        return sys;
    }

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystemContainerPerMesh);

private:
    std::map<MyKey, IDiscreteSystem*> _map_sys;
};

}
