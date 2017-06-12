/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeQuad9.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeQuad9 : public TemplateShapeFunction<2, 9> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[8] = (1.0 - pt[0] * pt[0]) * ( 1.0 - pt[1] * pt[1]);
        N[7] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 + pt[0]) - 0.5 * N[8];
        N[6] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 - pt[1]) - 0.5 * N[8];
        N[5] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 - pt[0]) - 0.5 * N[8];
        N[4] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 + pt[1]) - 0.5 * N[8];
        N[3] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]) - 0.5 * N[6] - 0.5 * N[7] - 0.25 * N[8];
        N[2] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]) - 0.5 * N[5] - 0.5 * N[6] - 0.25 * N[8];
        N[1] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]) - 0.5 * N[4] - 0.5 * N[5] - 0.25 * N[8];
        N[0] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]) - 0.5 * N[4] - 0.5 * N[7] - 0.25 * N[8];
    };

    void computeGradShapeFunction(const double* u, double *dN9)
    {
        dN9[8] = -2.0 * u[0] * (1.0 - u[1] * u[1]);
        dN9[7] = +0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
        dN9[6] = -1.0 * u[0] * (1.0 - u[1]) - 0.5 * dN9[8];
        dN9[5] = -0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
        dN9[4] = -1.0 * u[0] * (1.0 + u[1]) - 0.5 * dN9[8];
        dN9[3] = +0.25 * (1 - u[1]) - 0.5 * dN9[6] - 0.5 * dN9[7] - 0.25 * dN9[8];
        dN9[2] = -0.25 * (1 - u[1]) - 0.5 * dN9[5] - 0.5 * dN9[6] - 0.25 * dN9[8];
        dN9[1] = -0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[5] - 0.25 * dN9[8];
        dN9[0] = +0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[7] - 0.25 * dN9[8];

        dN9[17] = -2.0 * u[1] * (1.0 - u[0] * u[0]);
        dN9[16] = -1.0 * u[1] * (1.0 + u[0]) - 0.5 * dN9[17];
        dN9[15] = -0.5 * (1.0 - u[0] * u[0]) - 0.5 * dN9[17];
        dN9[14] = -1.0 * u[1] * (1.0 - u[0]) - 0.5 * dN9[17];
        dN9[13] = +0.5 * (1 - u[0] * u[0]) - 0.5 * dN9[17];
        dN9[12] = -0.25 * (1 + u[0]) - 0.5 * dN9[15] - 0.5 * dN9[16] - 0.25 * dN9[17];
        dN9[11] = -0.25 * (1 - u[0]) - 0.5 * dN9[14] - 0.5 * dN9[15] - 0.25 * dN9[17];
        dN9[10] = +0.25 * (1 - u[0]) - 0.5 * dN9[13] - 0.5 * dN9[14] - 0.25 * dN9[17];
        dN9[9] = +0.25 * (1 + u[0]) - 0.5 * dN9[13] - 0.5 * dN9[16] - 0.25 * dN9[17];
    };
};

}
