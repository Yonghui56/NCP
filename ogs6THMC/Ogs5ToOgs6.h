/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5ToOgs6.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemIO/ogs5/Ogs5FemIO.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/PorousMedia.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "GeoLib/GEOObjects.h"
#include "Ogs6FemData.h"

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertSolidProperty(const ogs5::CSolidProperties &msp, MaterialLib::Solid &solid);

void convertPorousMediumProperty(const ogs5::Ogs5FemData &ogs5fem, const ogs5::CMediumProperties &mmp, MaterialLib::PorousMedia &pm);

bool convert(const ogs5::Ogs5FemData &ogs5fem, Ogs6FemData &ogs6fem, BaseLib::Options &option);
};
} //end
