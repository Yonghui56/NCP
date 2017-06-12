/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemGaussIntegrationFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"
#include "FemIntegrationGaussLine.h"
#include "FemIntegrationGaussQuad.h"
#include "FemIntegrationGaussTriangle.h"
#include "FemIntegrationGaussTetra.h"

namespace FemLib
{

class FemGaussIntegrationFactory
{
public:
    static AbstractFemIntegrationGaussBase* create(MeshLib::IElement *e) 
    {
        switch (e->getShapeType()) {
            case MeshLib::ElementShape::QUAD:
                return new FemIntegrationGaussQuad();
            case MeshLib::ElementShape::TRIANGLE:
                return new FemIntegrationGaussQuad();
            case MeshLib::ElementShape::TETRAHEDRON:
                return new FemIntegrationGaussTetra();
            default:
                return 0;
        }
    }
};

}
