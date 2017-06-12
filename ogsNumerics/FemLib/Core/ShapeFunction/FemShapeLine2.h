/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeLine2.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{


class FemShapeLine2 : public TemplateShapeFunction<1, 2> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = 1.0 - pt[0];
        N[1] = 1.0 + pt[0];
        for (int i = 0; i < 2; i++)
            N[i] *= 0.5;
    };

    void computeGradShapeFunction(const double*, double *dN)
    {
        dN[0] = -0.5;
        dN[1] = 0.5;
    };
};

}
