/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeHex8.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeHex8 : public TemplateShapeFunction<3, 8> 
{
    virtual void computeShapeFunction(const double* x, double* N8)
    {
        N8[0] = (1.0 + x[0]) * (1.0 + x[1]) * (1.0 + x[2]);
        N8[1] = (1.0 - x[0]) * (1.0 + x[1]) * (1.0 + x[2]);
        N8[2] = (1.0 - x[0]) * (1.0 - x[1]) * (1.0 + x[2]);
        N8[3] = (1.0 + x[0]) * (1.0 - x[1]) * (1.0 + x[2]);
        N8[4] = (1.0 + x[0]) * (1.0 + x[1]) * (1.0 - x[2]);
        N8[5] = (1.0 - x[0]) * (1.0 + x[1]) * (1.0 - x[2]);
        N8[6] = (1.0 - x[0]) * (1.0 - x[1]) * (1.0 - x[2]);
        N8[7] = (1.0 + x[0]) * (1.0 - x[1]) * (1.0 - x[2]);
        for (size_t i = 0; i < 8; i++)
            N8[i] *= 0.125;
    }

    virtual void computeGradShapeFunction(const double* x, double* dN8)
    {
        double r =  x[0];
        double s =  x[1];
        double t =  x[2];
        dN8[0] = +(1.0 + s) * (1.0 + t);
        dN8[1] = -(1.0 + s) * (1.0 + t);
        dN8[2] = -(1.0 - s) * (1.0 + t);
        dN8[3] = +(1.0 - s) * (1.0 + t);

        dN8[4] = +(1.0 + s) * (1.0 - t);
        dN8[5] = -(1.0 + s) * (1.0 - t);
        dN8[6] = -(1.0 - s) * (1.0 - t);
        dN8[7] = +(1.0 - s) * (1.0 - t);

        dN8[8] = +(1.0 + r) * (1.0 + t);
        dN8[9] = +(1.0 - r) * (1.0 + t);
        dN8[10] = -(1.0 - r) * (1.0 + t);
        dN8[11] = -(1.0 + r) * (1.0 + t);

        dN8[12] = +(1.0 + r) * (1.0 - t);
        dN8[13] = +(1.0 - r) * (1.0 - t);
        dN8[14] = -(1.0 - r) * (1.0 - t);
        dN8[15] = -(1.0 + r) * (1.0 - t);

        dN8[16] = +(1.0 + r) * (1.0 + s);
        dN8[17] = +(1.0 - r) * (1.0 + s);
        dN8[18] = +(1.0 - r) * (1.0 - s);
        dN8[19] = +(1.0 + r) * (1.0 - s);

        dN8[20] = -(1.0 + r) * (1.0 + s);
        dN8[21] = -(1.0 - r) * (1.0 + s);
        dN8[22] = -(1.0 - r) * (1.0 - s);
        dN8[23] = -(1.0 + r) * (1.0 - s);

        for (size_t i = 0; i < 24; i++)
            dN8[i] *= 0.125;

    }
};

}
