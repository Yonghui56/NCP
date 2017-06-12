/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file C0InterfaceElements.h
 *
 * Created on 2012-11-22 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Extrapolation/FemExtrapolation.h"
#include "THMmfFiniteElementType.h"
#include "TemplateInterfaceFe.h"
#include "FemShapeInterfaceGoodmanQuad4.h"
#include "FemShapeInterfaceGoodmanTri3.h"

namespace THMmf
{

typedef TemplateInterfaceFe<2, 4, 1, FemShapeInterfaceGoodmanQuad4, FemLib::FemIntegrationGaussQuad, FemLib::FeExtrapolationGaussLinear> IE_QUAD4;
typedef TemplateInterfaceFe<2, 3, 1, FemShapeInterfaceGoodmanTri3, FemLib::FemIntegrationGaussTriangle, FemLib::FeExtrapolationGaussLinear> IE_TRI3;

}
