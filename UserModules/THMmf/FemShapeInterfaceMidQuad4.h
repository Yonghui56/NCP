/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceMidQuad4.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/ShapeFunction/TemplateShapeFunction.h"
#include "FemLib/Core/ShapeFunction/FemShapeLine2.h"

namespace THMmf
{

/**
 * \brief Shape function for mid-plane interface elements with Quad4
 *
 * mean(u) = 1/2* (u_top + u_bottom)
 *
 *  3                2
 *  *----------------*
 *  |                |
 *  *----------------*
 *  0                1
 */
class FemShapeInterfaceMidQuad4 : public FemLib::TemplateShapeFunction<2, 4>
{
public:
    virtual void computeShapeFunction(const double* pt, double *N4)
    {
        _fe_line2.computeShapeFunction(pt, N4);

        N4[2] = N4[1];
        N4[3] = N4[0];
        for (size_t i = 0; i < 4; i++)
            N4[i] *= 0.5;
    };

    virtual void computeGradShapeFunction(const double* pt, double *dN8)
    {
        //dN/dr
        _fe_line2.computeGradShapeFunction(pt, dN8);
        dN8[2] = dN8[1];
        dN8[3] = dN8[0];
        for (int i = 0; i < 4; i++)
            dN8[i] *= 0.5;
        //dN/ds
        for (int i = 4; i < 8; i++)
            dN8[i] = 0.0;
    };

private:
    FemLib::FemShapeLine2 _fe_line2;
};

}

