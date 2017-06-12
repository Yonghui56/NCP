/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DefaultFemElementCatalogBuilder.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include <map>
#include <cassert>

#include "FemLib/Core/Element/FiniteElementType.h"
#include "FemLib/Core/Element/C0IsoparametricElements.h"
#include "FemLib/Core/Element/TRI3CONST.h"
#include "FemElementCatalog.h"

namespace FemLib
{

class DefaultFemElementCatalogBuilder
{
public:
    static void construct(FemElementCatalog &feCatalog)
    {
        feCatalog.registerFeType<LINE2>(FiniteElementType::LINE2);
        feCatalog.registerFeType<LINE3>(FiniteElementType::LINE3);
        feCatalog.registerFeType<QUAD4>(FiniteElementType::QUAD4);
        feCatalog.registerFeType<QUAD8>(FiniteElementType::QUAD8);
        feCatalog.registerFeType<QUAD9>(FiniteElementType::QUAD9);
        feCatalog.registerFeType<TRI3CONST>(FiniteElementType::TRI3CONST);
        feCatalog.registerFeType<TRI3>(FiniteElementType::TRI3);
        feCatalog.registerFeType<TRI6>(FiniteElementType::TRI6);
        feCatalog.registerFeType<TET4>(FiniteElementType::TET4);
        feCatalog.registerFeType<TET10>(FiniteElementType::TET10);
    }
};

}
