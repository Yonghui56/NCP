/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeQuad4.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeQuad4 : public TemplateShapeFunction<2, 4> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = (1.0 + pt[0]) * (1.0 + pt[1]);
        N[1] = (1.0 - pt[0]) * (1.0 + pt[1]);
        N[2] = (1.0 - pt[0]) * (1.0 - pt[1]);
        N[3] = (1.0 + pt[0]) * (1.0 - pt[1]);
        for (int i = 0; i < 4; i++)
            N[i] *= 0.25;
    };

    void computeGradShapeFunction(const double* u, double *dN4)
    {
        dN4[0] = +(1.0 + u[1]);
        dN4[1] = -(1.0 + u[1]);
        dN4[2] = -(1.0 - u[1]);
        dN4[3] = +(1.0 - u[1]);
        dN4[4] = +(1.0 + u[0]);
        dN4[5] = +(1.0 - u[0]);
        dN4[6] = -(1.0 - u[0]);
        dN4[7] = -(1.0 + u[0]);
        for (int i = 0; i < 8; i++)
            dN4[i] *= 0.25;
    };
};

}
