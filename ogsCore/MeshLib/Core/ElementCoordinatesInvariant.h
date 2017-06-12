/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesInvariant
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElementCoordinatesMapping.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

/**
 * \brief This class keep original coordinates of element nodes
 */
class ElementCoordinatesInvariant
: public IElementCoordinatesMapping
{
public:
    /**
     *
     * @param msh
     * @param e
     */
    ElementCoordinatesInvariant(const IMesh* msh, const IElement* e) : _msh(msh), _e(e)
    {
    };
    
    ///
    virtual ~ElementCoordinatesInvariant() {};

    ///
    GeoLib::Point* getNodePoint(size_t local_id)
    {
        return (GeoLib::Point*)_msh->getNodeCoordinatesRef(_e->getNodeID(local_id));
    }

private:
    const IMesh* _msh;
    const IElement* _e;
};

}
