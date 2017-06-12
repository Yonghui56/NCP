/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file C0IsoparametricElements.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Extrapolation/FemExtrapolation.h"
#include "TemplateIsoparametric.h"

namespace FemLib
{
             
typedef TemplateIsoparametric<1, 2, 1, FemShapeLine2, FemIntegrationGaussLine, FeExtrapolationGaussLinear> LINE2;
typedef TemplateIsoparametric<1, 3, 2, FemShapeLine3, FemIntegrationGaussLine, FeExtrapolationGaussLinear> LINE3;
typedef TemplateIsoparametric<2, 4, 1, FemShapeQuad4, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD4;
typedef TemplateIsoparametric<2, 8, 2, FemShapeQuad8, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD8;
typedef TemplateIsoparametric<2, 9, 2, FemShapeQuad9, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD9;
typedef TemplateIsoparametric<2, 3, 1, FemShapeTriangle3, FemIntegrationGaussTriangle, FeExtrapolationGaussLinear> TRI3;
typedef TemplateIsoparametric<2, 6, 2, FemShapeTriangle6, FemIntegrationGaussTriangle, FeExtrapolationGaussLinear> TRI6;
typedef TemplateIsoparametric<3, 4, 1, FemShapeTetra4, FemIntegrationGaussTetra, FeExtrapolationGaussLinear> TET4;
typedef TemplateIsoparametric<3, 10, 2, FemShapeTetra10, FemIntegrationGaussTetra, FeExtrapolationGaussLinear> TET10;

}
