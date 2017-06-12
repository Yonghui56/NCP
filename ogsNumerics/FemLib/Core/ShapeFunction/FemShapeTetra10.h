/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeTetra10.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeTetra10 : public TemplateShapeFunction<3, 10>
{
    virtual void computeShapeFunction(const double* x, double* N10)
    {
        N10[0] = 2. * (1 - x[0] - x[1] - x[2]) * (0.5 - x[0] - x[1] - x[2]);
        N10[1] = x[0] * (2. * x[0] - 1);
        N10[2] = x[1] * (2. * x[1] - 1);
        N10[3] = x[2] * (2. * x[2] - 1);
        N10[4] = 4.0 * x[0] * (1.0 - x[0] - x[1] - x[2]);
        N10[5] = 4.0 * x[0] * x[1];
        N10[6] = 4.0 * x[1] * (1.0 - x[0] - x[1] - x[2]);
        N10[7] = 4.0 * x[0] * x[2];
        N10[8] = 4.0 * x[1] * x[2];
        N10[9] = 4.0 * x[2] * (1.0 - x[0] - x[1] - x[2]);
    }

    virtual void computeGradShapeFunction(const double* x, double* dN10)
    {
        dN10[0] = 4.0 * (x[0] + x[1] + x[2]) - 3.0;
        dN10[1] = 4. * x[0] - 1.;
        dN10[2] = 0.0;
        dN10[3] = 0.0;
        dN10[4] = 4.0 * (1.0 - 2.0 * x[0] - x[1] - x[2]);
        dN10[5] = 4.0 * x[1];
        dN10[6] = -4.0 * x[1];
        dN10[7] = 4.0 * x[2];
        dN10[8] = 0.0;
        dN10[9] = -4.0 * x[2];

        dN10[10] =  4. * (x[0] + x[1] + x[2]) - 3.;
        dN10[11] = 0.0;
        dN10[12] = 4. * x[1] - 1.;
        dN10[13] = 0.;
        dN10[14] = -4.0 * x[0];
        dN10[15] = 4.0 * x[0];
        dN10[16] = 4.0 * (1.0 - x[0] - 2.0 * x[1] - x[2]);
        dN10[17] = 0.0;
        dN10[18] = 4.0 * x[2];
        dN10[19] = -4.0 * x[2];

        dN10[20] = 4. * (x[0] + x[1] + x[2]) - 3.;
        dN10[21] = 0.;
        dN10[22] = 0.;
        dN10[23] = 4. * x[2] - 1.;
        dN10[24] = -4.0 * x[0];
        dN10[25] = 0.0;
        dN10[26] = -4.0 * x[1];
        dN10[27] = 4.0 * x[0];
        dN10[28] = 4.0 * x[1];
        dN10[29] = 4.0 * (1.0 - x[0] - x[1] - 2.0 * x[2]);
    }
};

}
