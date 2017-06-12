/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "Element.h"

namespace MeshLib
{

/**
 * \brief Mesh element factory
 */
class ElemenetFactory
{
public:
    /// create a new element object
    static IElement* createNewElement(const ElementShape::type t)
    {
        if (t == ElementShape::LINE)    return new Line();
        else if (t == ElementShape::QUAD) return new Quadrirateral();
        else if (t == ElementShape::HEXAHEDRON) return new Hexahedron();
        else if (t == ElementShape::TRIANGLE) return new Triangle();
        else if (t == ElementShape::TETRAHEDRON) return new Tetrahedron();
        else if (t == ElementShape::PRISM) return new Prism();
        else if (t == ElementShape::PYRAMID) return new Pyramid();
        return NULL;
    };
};

} // end namespace

