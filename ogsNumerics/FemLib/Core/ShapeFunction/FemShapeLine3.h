/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeLine3.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeLine3 : public TemplateShapeFunction<1, 3> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = 0.5 * pt[0] * (pt[0] - 1.0);
        N[1] = 0.5 * pt[0] * (pt[0] + 1.0);
        N[2] = 1.0 - pt[0] * pt[0];
    };

    void computeGradShapeFunction(const double* pt, double *dN)
    {
        dN[0] = pt[0] - 0.5;
        dN[1] = pt[0] + 0.5;
        dN[2] = -2.0 * pt[0];
    };
};

}
