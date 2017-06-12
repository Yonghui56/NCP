/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReaction.cpp
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#include "chemReaction.h"

namespace ogsChem 
{

chemReaction::chemReaction(void)
{
}

chemReaction::~chemReaction(void)
{
}


void chemReaction::addComp( ChemComp* pComp, double stoi )
{
    _vecComponents.push_back(pComp); 
    _vecStoi.push_back(stoi); 
}

} // end of namespace
