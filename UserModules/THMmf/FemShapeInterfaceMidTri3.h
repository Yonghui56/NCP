/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceMidTri3.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/ShapeFunction/TemplateShapeFunction.h"
#include "FemLib/Core/ShapeFunction/FemShapeLine2.h"

namespace THMmf
{

/**
 * \brief Shape function for mid-plane interface elements
 */
class FemShapeInterfaceMidTri3 : public FemLib::TemplateShapeFunction<2, 3>
{
public:
    virtual void computeShapeFunction(const double* pt, double *N3)
    {
        _fe_line2.computeShapeFunction(pt, N3);
        //*Notice: node order of an element
        //          2
        //       _ -*
        //   _ -    |
        //  *-------*
        //  0       1

        N3[2] = N3[1];
        for (size_t i = 1; i < 3; i++)
            N3[i] *= 0.5;
    };

    virtual void computeGradShapeFunction(const double* pt, double *dN6)
    {
        // dN/dr
        _fe_line2.computeGradShapeFunction(pt, dN6);
        dN6[2] = dN6[1];
        for (size_t i = 1; i < 3; i++)
            dN6[i] *= 0.5;
        // dN/ds
        for (size_t i = 3; i < 6; i++)
            dN6[i] = 0.0;
    };

    void GradShapeFunctionInterfaceMidTriEnd(double* dN6, const double* u)
    {
        //*Notice: node order of an element
        //  2
        //  *- _
        //  |    - _
        //  *-------*
        //  0       1
        // dN/dr
        double dN2[2] = {};
        GradShapeFunctionLine(dN2, u);
        dN6[0] = 0.5 * dN2[0];
        dN6[1] = dN2[1];
        dN6[2] = dN6[0];
        // dN/ds
        for (int i = 3; i < 6; i++)
            dN6[i] = 0.0;
    }

private:
    FemLib::FemShapeLine2 _fe_line2;
};

}

