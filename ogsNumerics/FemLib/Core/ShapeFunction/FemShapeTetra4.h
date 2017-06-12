/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeTetra4.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeTetra4 : public TemplateShapeFunction<3, 4> 
{
    virtual void computeShapeFunction(const double* x, double* N)
    {
        N[0] = 1. - x[0] - x[1] - x[2];
        N[1] = x[0];
        N[2] = x[1];
        N[3] = x[2];
    }

    virtual void computeGradShapeFunction(const double* /*x*/, double* dNt4)
    {
        //dr
        dNt4[0] = -1.0;
        dNt4[1] = 1.0;
        dNt4[2] = 0.0;
        dNt4[3] = 0.0;

        //ds
        dNt4[4] = -1.0;
        dNt4[5] = 0.0;
        dNt4[6] = 1.0;
        dNt4[7] = 0.0;

        //dt
        dNt4[8] = -1.0;
        dNt4[9] = 0.0;
        dNt4[10] = 0.0;
        dNt4[11] = 1.0;
    }
};

}
