/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeTriangle3.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeTriangle3 : public TemplateShapeFunction<2, 3> 
{
    void computeShapeFunction(const double* u, double *N)
    {
        N[0] = 1. - u[0] - u[1];
        N[1] = u[0];
        N[2] = u[1];
    };
    void computeGradShapeFunction(const double*, double *dN3)
    {
        //   d()/dL_1
        dN3[0] = -1.0;
        dN3[1] =  1.0;
        dN3[2] =  0.0;
        //   d()/dL_2
        dN3[3] = -1.0;
        dN3[4] = 0.0;
        dN3[5] = 1.0;
    }
};

}
