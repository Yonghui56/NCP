/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Fluid.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"

namespace MaterialLib
{

struct Fluid
{
    NumLib::ITXFunction* density;
	NumLib::ITXFunction* molar_density;
    NumLib::ITXFunction* dynamic_viscosity;
	NumLib::ITXFunction* molar_mass;
    NumLib::ITXFunction* specific_heat;
    NumLib::ITXFunction* thermal_conductivity;
	NumLib::ITXFunction* drho_dp;

    Fluid()
    {
		BaseLib::zeroObject(
				density,
				molar_density,     //molar density of gas phase     //molar density of liquid phase
                dynamic_viscosity,
				molar_mass,
                specific_heat,
                thermal_conductivity,
				drho_dp
                );
    }
    ~Fluid()
    {
        BaseLib::releaseObject(
                density,
				molar_density,     //molar density of gas phase
                dynamic_viscosity,
				molar_mass,
                specific_heat,
                thermal_conductivity,
				drho_dp
                );
    }
};

} //end
