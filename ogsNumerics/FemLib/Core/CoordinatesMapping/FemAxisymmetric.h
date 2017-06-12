/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemAxisymmetric.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <iostream>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Mapping.h"
#include "FemLib/Core/ShapeFunction.h"

namespace FemLib
{

/**
 * Axisymmetric class
 * - coordinates (r, theta, z)
 * - dxdydz = 2pi*r*drdz
 */
class FemAxisymmetric
{
public:
    FemAxisymmetric() : _r(.0) {};
    virtual ~FemAxisymmetric() {};

    virtual void computeMappingFunctions(double* natural_pt)
    {
        //FemIsoparametricMapping::computeMappingFunctions(natural_pt, compType);
        _r = this->computeRadius();
    };


    /// return radius
    double getRadius() const {
        return _r;
    }

private:
    double _r;

    double computeRadius() {
        double r = .0;
        //MeshLib::IElement *e = getElement();
        //double *shapefct = getShapeFunction();
        //for (size_t i=0; i<e->getNumberOfNodes(); i++) {
        //    const GeoLib::Point *pt = e->getMappedGeometry()->getNodePoint(i);
        //    r += shapefct[i] * (*pt)[0];
        //}
        return r;
    }
};

}
