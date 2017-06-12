/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFeObjectContainer.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FeObjectContainerPerFeType.h"
#include "FemElementCatalog.h"
#include "LagrangeFemElementCatalogBuilder.h"

namespace FemLib
{

/**
 * \brief Lagrangian finite element object containers
 */
class LagrangeFeObjectContainer: public FeObjectContainerPerFeType
{
public:
    /**
     *
     * @param msh
     */
    explicit LagrangeFeObjectContainer(MeshLib::IMesh* msh)
    : FeObjectContainerPerFeType(&_fe_catalog, &_shape2fetype, msh), _order(1)
    {
        LagrangeFemElementCatalogBuilder::construct(_fe_catalog, _shape2fetype);
    };

    /**
     * Copy constructor
     * @param src
     */
    LagrangeFeObjectContainer(const LagrangeFeObjectContainer &src)
    : FeObjectContainerPerFeType(src), _shape2fetype(src._shape2fetype), _order(src._order)
    {
    }

    /**
     *
     */
    virtual ~LagrangeFeObjectContainer() {};

    virtual LagrangeFeObjectContainer* clone() const
    {
        return new LagrangeFeObjectContainer(*this);
    }

    /**
     * set polynomial order
     * @param order
     */
    void setPolynomialOrder(size_t order) { _order = order; };

    /**
     * get FE object
     * @param e
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e);

    /**
     * get FE object
     * @param e
     * @param order
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, size_t order);

private:
    FemElementCatalog _fe_catalog;
    MeshElementShapeToFemElementType _shape2fetype;
    size_t _order;
};



}
