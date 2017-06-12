/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemcomp.cpp
 *
 * Created on 2012-08-23 by Haibing Shao
 */

#include "chemcomp.h"

namespace ogsChem
{

ChemComp::ChemComp(void):
    _name("NaN"), _composition("NaN"), _mole_mass(0.0),
    _charge(0.0), _G0(0.0)
{
    // set pointers to zero
    _mPhase    = 0;
    _mReaction = 0;
    _index = 0;
    _mole_volume = .0;
    neg_gamma_1 = .0;
    neg_gamma_2 = .0;
}

ChemComp::ChemComp(MaterialLib::Compound* mCompund):
    _name(mCompund->name), _composition("NaN"), _mole_mass(0.0),
    _charge(0.0), _G0(0.0)
{
    // set pointers to zero
    _mPhase    = 0;
    _mReaction = 0;
    _index = 0;
    _mole_volume = .0;
    neg_gamma_1 = .0;
    neg_gamma_2 = .0;
}

ChemComp::~ChemComp(void)
{
}

double ChemComp::get_mole_mass(void) const
{
return _mole_mass;
}

double ChemComp::get_mole_volume(void) const
{
return _mole_volume;
}

double ChemComp::get_charge(void) const
{
return _charge;
}

} // end of namespace
