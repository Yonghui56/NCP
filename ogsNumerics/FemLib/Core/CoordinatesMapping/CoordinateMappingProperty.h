/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file CoordinateMappingProperty.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"

namespace FemLib
{

/**
 * \brief Properties of coordinate mapping at particular locations 
 */
struct CoordinateMappingProperty
{
    ///// mapping locations in computed coordinates
    //double r[3];
    /// shape function N(r)
    MathLib::LocalMatrix *shape_r;
    /// gradient of shape functions, dN(r)/dr
    MathLib::LocalMatrix *dshape_dr;
    /// gradient of shape functions, dN(r)/dx
    MathLib::LocalMatrix *dshape_dx;
    /// Jacobian matrix, J=dx/dr
    MathLib::LocalMatrix *jacobian_dxdr;
    /// determinant of the Jacobian
    double det_jacobian;
    /// inverse of the Jacobian
    MathLib::LocalMatrix *inv_jacobian_drdx;

    CoordinateMappingProperty()
    {
        shape_r = new MathLib::LocalMatrix();
        dshape_dr = new MathLib::LocalMatrix();
        dshape_dx = new MathLib::LocalMatrix();
        jacobian_dxdr = new MathLib::LocalMatrix();
        inv_jacobian_drdx = new MathLib::LocalMatrix();
        det_jacobian = .0;
    }

    ~CoordinateMappingProperty()
    {
        delete shape_r;
        delete dshape_dr;
        delete dshape_dx;
        delete jacobian_dxdr;
        delete inv_jacobian_drdx;
    }
};

}
