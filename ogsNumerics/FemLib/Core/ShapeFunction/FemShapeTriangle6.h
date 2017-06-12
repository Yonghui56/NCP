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

class FemShapeTriangle6 : public TemplateShapeFunction<2, 6> 
{
    void computeShapeFunction(const double* u, double *N6)
    {
        N6[0] = 2. * (1. - u[0] - u[1]) * (0.5 - u[0] - u[1]);
        N6[1] = u[0] * (2. * u[0] - 1.);
        N6[2] = u[1] * (2. * u[1] - 1.);
        N6[3] = 4. * u[0] * (1. - u[0] - u[1]);
        N6[4] = 4. * u[0] * u[1];
        N6[5] = 4. * u[1] * (1. - u[0] - u[1]);
    };
    void computeGradShapeFunction(const double* u, double *dN6)
    {
        dN6[0] = 4. * (u[0] + u[1]) - 3.;     // dN1/dL1
        dN6[6] = 4. * (u[0] + u[1]) - 3.;     // dN1/dL2

        dN6[1] = 4. * u[0] - 1.;              // dN2/dL1
        dN6[7] = 0.;                          // dN2/dL2

        dN6[2] = 0.;                          // dN3/dL1
        dN6[8] = 4. * u[1] - 1.;              // dN3/dL2

        dN6[3] =  4. * (1 - 2. * u[0] - u[1]); // dN4/dL1
        dN6[9] = -4. * u[0];                  // dN4/dL2

        dN6[4] = 4. * u[1];                   // dN5/dL1
        dN6[10] = 4. * u[0];                  // dN5/dL2

        dN6[5] = -4. * u[1];                  // dN6/dL1
        dN6[11] = 4. * (1 - u[0] - 2. * u[1]); // dN6/dL2
    }
};

}
