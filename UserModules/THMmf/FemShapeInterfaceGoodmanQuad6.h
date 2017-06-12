/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceGoodmanQuad6.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/ShapeFunction/TemplateShapeFunction.h"
#include "FemLib/Core/ShapeFunction/FemShapeLine3.h"

namespace THMmf
{

/**
 *
 */
class FemShapeInterfaceGoodmanQuad6 : public FemLib::TemplateShapeFunction<2, 6>
{
public:
    void computeShapeFunction(const double* pt, double *N6)
    {
        double N3[3] = {};
        _fe_line3.computeShapeFunction(pt, N3);

        //*Notice: node order of an element
        //  3        5        2
        //  *--------*--------*
        //  |                 |
        //  *--------*--------*
        //  0        4        1
        N6[0] = -N3[0];
        N6[1] = -N3[1];
        N6[4] = -N3[2];
        N6[3] = N3[0];
        N6[2] = N3[1];
        N6[5] = N3[2];
    };

    void computeGradShapeFunction(const double*, double *)
    {
        //TODO
    };

private:
    FemLib::FemShapeLine3 _fe_line3;
};

}

