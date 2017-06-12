/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFeObjectContainer.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "LagrangeFeObjectContainer.h"


namespace FemLib
{

IFiniteElement* LagrangeFeObjectContainer::getFeObject(const MeshLib::IElement &e)
{
    e.setCurrentOrder(_order);
    IFiniteElement* fe = FeObjectContainerPerFeType::getFeObject(e);
    fe->configure(*const_cast<MeshLib::IElement*>(&e));
    return fe;
}

IFiniteElement* LagrangeFeObjectContainer::getFeObject(const MeshLib::IElement &e, size_t order)
{
    setPolynomialOrder(order);
    return getFeObject(e);
}


} //end
