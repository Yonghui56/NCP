/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeHex20.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{

class FemShapeHex20 : public TemplateShapeFunction<3, 20> 
{
    virtual void computeShapeFunction(const double* x, double* N20)
    {
        double r = x[0];
        double s = x[1];
        double t = x[2];

        N20[0] = ShapeFunctionHexHQ_Corner(r,s,t);
        N20[1] = ShapeFunctionHexHQ_Corner(-r,s,t);
        N20[2] = ShapeFunctionHexHQ_Corner(-r,-s,t);
        N20[3] = ShapeFunctionHexHQ_Corner(r,-s,t);
        N20[4] = ShapeFunctionHexHQ_Corner(r,s,-t);
        N20[5] = ShapeFunctionHexHQ_Corner(-r,s,-t);
        N20[6] = ShapeFunctionHexHQ_Corner(-r,-s,-t);
        N20[7] = ShapeFunctionHexHQ_Corner(r,-s,-t);

        N20[8] = ShapeFunctionHexHQ_Middle(r,s,t);
        N20[10] = ShapeFunctionHexHQ_Middle(r,-s,t);
        N20[14] = ShapeFunctionHexHQ_Middle(r,-s,-t);
        N20[12] = ShapeFunctionHexHQ_Middle(r,s,-t);

        N20[11] = ShapeFunctionHexHQ_Middle(s,t,r);
        N20[15] = ShapeFunctionHexHQ_Middle(s,-t,r);
        N20[13] = ShapeFunctionHexHQ_Middle(s,-t,-r);
        N20[9]  =  ShapeFunctionHexHQ_Middle(s,t,-r);

        N20[16] = ShapeFunctionHexHQ_Middle(t,r,s);
        N20[17] = ShapeFunctionHexHQ_Middle(t,-r,s);
        N20[18] = ShapeFunctionHexHQ_Middle(t,-r,-s);
        N20[19] = ShapeFunctionHexHQ_Middle(t,r,-s);

    }

    virtual void computeGradShapeFunction(const double* x, double* dN20)
    {
        int co;
        double r = x[0];
        double s = x[1];
        double t = x[2];
        static double sign1[] = {-1.0, 1.0,1.0};
        static double sign2[] = { 1.0,-1.0,1.0};
        static double sign3[] = { 1.0, 1.0,-1.0};
        for(int i = 0; i < 3; i++)
        {
            dN20[20 * i + 0] = dShapeFunctionHexHQ_Corner(r,s,t,i);
            dN20[20 * i + 1] = sign1[i] * dShapeFunctionHexHQ_Corner(-r,s,t,i);
            dN20[20 * i + 2] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Corner(-r,-s,t,i);
            dN20[20 * i + 3] = sign2[i] * dShapeFunctionHexHQ_Corner(r,-s,t,i);
            dN20[20 * i + 4] = sign3[i] * dShapeFunctionHexHQ_Corner(r,s,-t,i);
            dN20[20 * i + 5] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r,s,-t,i);
            dN20[20 * i + 6] = sign1[i] * sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r,
                                                                                           -s,
                                                                                           -t,
                                                                                           i);
            dN20[20 * i + 7] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(r,-s,-t,i);

            dN20[20 * i + 8] =  dShapeFunctionHexHQ_Middle(r,s,t,i);
            dN20[20 * i + 10] = sign2[i] * dShapeFunctionHexHQ_Middle(r,-s,t,i);
            dN20[20 * i + 14] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Middle(r,-s,-t,i);
            dN20[20 * i + 12] = sign3[i] * dShapeFunctionHexHQ_Middle(r,s,-t,i);

            co = (i + 2) % 3;
            dN20[20 * i + 11] = dShapeFunctionHexHQ_Middle(s,t,r,co);
            dN20[20 * i + 15] = sign3[i] * dShapeFunctionHexHQ_Middle(s,-t,r,co);
            dN20[20 * i + 13] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Middle(s,-t,-r,co);
            dN20[20 * i + 9] =  sign1[i] * dShapeFunctionHexHQ_Middle(s,t,-r,co);

            co = (i + 1) % 3;
            dN20[20 * i + 16] = dShapeFunctionHexHQ_Middle(t,r,s,co);
            dN20[20 * i + 17] = sign1[i] * dShapeFunctionHexHQ_Middle(t,-r,s,co);
            dN20[20 * i + 18] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Middle(t,-r,-s,co);
            dN20[20 * i + 19] = sign2[i] * dShapeFunctionHexHQ_Middle(t,r,-s,co);
        }
    }

private:
    double ShapeFunctionHexHQ_Corner(const double r, const double s, const double t)
    {
        return 0.125 * (1 + r) * (1 + s) * (1 + t) * ( r + s + t - 2.0);
    }

    double ShapeFunctionHexHQ_Middle(const double r, const double s, const double t)
    {
        return 0.25 * (1 - r * r) * (1 + s) * (1 + t);
    }

    double dShapeFunctionHexHQ_Corner(const double r, const double s, const double t, const int ty)
    {
        switch(ty)
        {
        case 0:
            return 0.125 * (1 + s) * (1 + t) * (2.0 * r + s + t - 1.0);
            break;
        case 1:
            return 0.125 * (1 + t) * (1 + r) * (2.0 * s + r + t - 1.0);
            break;
        case 2:
            return 0.125 * (1 + r) * (1 + s) * (2.0 * t + s + r - 1.0);
            break;
        }
        return 0.0;
    }

    double dShapeFunctionHexHQ_Middle(const double r, const double s, const double t, const int ty)
    {
        switch(ty)
        {
        case 0:
            return -0.5 * r * (1 + s) * (1 + t);
            break;
        case 1:
            return 0.25 * (1 - r * r) * (1 + t);
            break;
        case 2:
            return 0.25 * (1 - r * r) * (1 + s);
            break;
        }
        return 0.0;
    }
};

}
