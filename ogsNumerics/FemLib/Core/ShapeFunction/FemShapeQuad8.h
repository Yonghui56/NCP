/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeQuad8.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeQuad8 : public TemplateShapeFunction<2, 8> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = -0.25 * (1.0 - pt[0]) * (1.0 - pt[1]) * (( 1.0 + pt[0] + pt[1]));
        N[1] =  0.25 * (1.0 + pt[0]) * (1.0 - pt[1]) * ((-1.0 + pt[0] - pt[1]));
        N[2] =  0.25 * (1.0 + pt[0]) * (1.0 + pt[1]) * ((-1.0 + pt[0] + pt[1]));
        N[3] = -0.25 * (1.0 - pt[0]) * (1.0 + pt[1]) * (( 1.0 + pt[0] - pt[1]));
        //
        N[4] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 - pt[1]);
        N[5] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 + pt[0]);
        N[6] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 + pt[1]);
        N[7] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 - pt[0]);
    };

    void computeGradShapeFunction(const double* pt, double *dN8)
    {
        double r = pt[0];
        double s = pt[1];

        //dN/dr
        dN8[0] = (1 - s) * (2 * r + s) * 0.25;
        dN8[1] = (1 - s) * (2 * r - s) * 0.25;
        dN8[2] = (1 + s) * (2 * r + s) * 0.25;
        dN8[3] = (1 + s) * (2 * r - s) * 0.25;
        dN8[4] = -r * (1 - s);
        dN8[5] = (1 - s * s) * 0.5;
        dN8[6] = -r * (1 + s);
        dN8[7] = -(1 - s * s) * 0.5;

        //dN/ds
        dN8[8] = (1 - r) * (r + 2 * s) * 0.25;
        dN8[9] = -(1 + r) * (r - 2 * s) * 0.25;
        dN8[10] = (1 + r) * (r + 2 * s) * 0.25;
        dN8[11] = -(1 - r) * (r - 2 * s) * 0.25;
        dN8[12] = -(1 - r * r) * 0.5;
        dN8[13] = -(1 + r) * s;
        dN8[14] = (1 - r * r) * 0.5;
        dN8[15] = -(1 - r) * s;
    };
};

}
