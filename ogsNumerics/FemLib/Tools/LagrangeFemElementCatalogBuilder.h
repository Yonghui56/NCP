/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFemElementCatalogBuilder.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include <map>
#include <cassert>

#include "FemLib/Core/Element/FiniteElementType.h"
#include "FemLib/Core/Element/C0IsoparametricElements.h"
#include "FemElementCatalog.h"
#include "MeshElementShapeToFemElementType.h"

namespace FemLib
{

/**
 * setup Lagrange FE type catalog
 */
class LagrangeFemElementCatalogBuilder
{
public:
    static void construct(FemElementCatalog &feCatalog, MeshElementShapeToFemElementType &shape2fetype)
    {
        feCatalog.registerFeType<LINE2>(FiniteElementType::LINE2);
        feCatalog.registerFeType<LINE3>(FiniteElementType::LINE3);
        feCatalog.registerFeType<QUAD4>(FiniteElementType::QUAD4);
        feCatalog.registerFeType<QUAD9>(FiniteElementType::QUAD9);
        feCatalog.registerFeType<TRI3>(FiniteElementType::TRI3);
        feCatalog.registerFeType<TRI6>(FiniteElementType::TRI6);
        feCatalog.registerFeType<TET4>(FiniteElementType::TET4);
        feCatalog.registerFeType<TET10>(FiniteElementType::TET10);

        shape2fetype.addFeType(MeshLib::ElementShape::LINE, 1, FiniteElementType::LINE2);
        shape2fetype.addFeType(MeshLib::ElementShape::LINE, 2, FiniteElementType::LINE3);
        shape2fetype.addFeType(MeshLib::ElementShape::QUAD, 1, FiniteElementType::QUAD4);
        shape2fetype.addFeType(MeshLib::ElementShape::QUAD, 2, FiniteElementType::QUAD9);
        shape2fetype.addFeType(MeshLib::ElementShape::TRIANGLE, 1, FiniteElementType::TRI3);
        shape2fetype.addFeType(MeshLib::ElementShape::TRIANGLE, 2, FiniteElementType::TRI6);
        shape2fetype.addFeType(MeshLib::ElementShape::TETRAHEDRON, 1, FiniteElementType::TET4);
        shape2fetype.addFeType(MeshLib::ElementShape::TETRAHEDRON, 2, FiniteElementType::TET10);
    }
};

}
