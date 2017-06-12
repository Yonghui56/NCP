/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceGoodmanTri3.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/ShapeFunction/TemplateShapeFunction.h"
#include "FemLib/Core/ShapeFunction/FemShapeLine3.h"

namespace THMmf
{


class FemShapeInterfaceGoodmanTri5 : public FemLib::TemplateShapeFunction<2, 5>
{
public:
    virtual void computeShapeFunction(const double* pt, double* N5)
    {
        double N3[3] = {};
        _fe_line3.computeShapeFunction(pt, N3);

        //*Notice: node order of an element
        //          2
        //     4 _ -*
        //   _ -*   |
        //  *---*---*
        //  0   3    1

        N5[0] = 0.0;
        N5[1] = -N3[1];
        N5[3] = -N3[2];
        N5[2] = N3[1];
        N5[4] = N3[2];
    };

    void ShapeFunctionGoodmanTriEnd(double* N3, const double* u)
    {
        double N2[2] = {};
        ShapeFunctionLine(N2, u);

        //*Notice: node order of an element
        //  2
        //  *- _
        //  |    - _
        //  *-------*
        //  0       1

        N3[0] = -N2[1];
        N3[1] = 0.0;
        N3[2] = N2[1];
    }


    virtual void computeGradShapeFunction(const double*, double *)
    {
    };

private:
    FemLib::FemShapeLine3 _fe_line3;
};

}
