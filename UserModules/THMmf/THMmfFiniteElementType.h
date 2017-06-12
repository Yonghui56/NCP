/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfFiniteElementType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace THMmf
{
/// Finite element type
struct THMmfFiniteElementType
{
    enum type {
        LINE2,
        LINE3,
        TRI3,
        TRI3CONST,
        TRI6,
        QUAD4,
        QUAD8,
        QUAD9,
        TET4,
        TET10,
        IE_QUAD4,
        IE_TRI3,
        INVALID
    };
};

}
