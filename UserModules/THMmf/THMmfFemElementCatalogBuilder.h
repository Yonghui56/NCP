/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfFemElementCatalogBuilder.h
 *
 * Created on 2012-11-22 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Tools/FemElementCatalog.h"
#include "FemLib/Core/Element/C0IsoparametricElements.h"
#include "C0InterfaceElements.h"
#include "THMmfFiniteElementType.h"

namespace THMmf
{

class THMmfFemElementCatalogBuilder
{
public:
    static FemLib::FemElementCatalog* construct()
    {
        FemLib::FemElementCatalog* feCatalog = new FemLib::FemElementCatalog();
        feCatalog->registerFeType<FemLib::LINE2>(THMmfFiniteElementType::LINE2);
        feCatalog->registerFeType<FemLib::LINE3>(THMmfFiniteElementType::LINE3);
        feCatalog->registerFeType<FemLib::QUAD4>(THMmfFiniteElementType::QUAD4);
        feCatalog->registerFeType<FemLib::QUAD8>(THMmfFiniteElementType::QUAD8);
        feCatalog->registerFeType<FemLib::QUAD9>(THMmfFiniteElementType::QUAD9);
        feCatalog->registerFeType<FemLib::TRI3>(THMmfFiniteElementType::TRI3);
        feCatalog->registerFeType<FemLib::TRI6>(THMmfFiniteElementType::TRI6);
        feCatalog->registerFeType<FemLib::TET4>(THMmfFiniteElementType::TET4);
        feCatalog->registerFeType<FemLib::TET10>(THMmfFiniteElementType::TET10);
        feCatalog->registerFeType<IE_TRI3>(THMmfFiniteElementType::IE_TRI3);
        feCatalog->registerFeType<IE_QUAD4>(THMmfFiniteElementType::IE_QUAD4);

        return feCatalog;
    }
};

}
