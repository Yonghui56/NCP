/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectCachePerFeType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemElementCatalog.h"

namespace FemLib
{

/**
 * \brief Cache system for finite element objects 
 */
class FeObjectCachePerFeType
{
public:
    /**
     *
     * @param fe_catalog
     * @param msh
     */
    explicit FeObjectCachePerFeType(const FemElementCatalog* fe_catalog, MeshLib::IMesh* msh)
    : _fe_catalog(fe_catalog), _msh(msh) {};

    /**
     *
     * @param src
     */
    FeObjectCachePerFeType(const FeObjectCachePerFeType& src);

    /**
     *
     */
    virtual ~FeObjectCachePerFeType()
    {
        BaseLib::releaseObjectsInStdMap(_mapFeObj);
    }

    /**
     * get a finite element object
     * @param fe_type
     * @return
     */
    IFiniteElement* getFeObject(const int fe_type);

    /**
     * get a pointer to this mesh
     * @return
     */
    MeshLib::IMesh* getMesh() const {return _msh;};

    /**
     * get a pointer to FE catalog
     * @return
     */
    const FemElementCatalog* getCatalog() const {return _fe_catalog;};

private:
    const FemElementCatalog* _fe_catalog;
    MeshLib::IMesh* _msh;
    std::map<int, IFiniteElement*> _mapFeObj;
};

}
