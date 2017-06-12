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
#include "FemLib/Core/ShapeFunction/FemShapeLine2.h"

namespace THMmf
{

/**
 * \brief Shape function for Goodman joint elements
 *
 * [[u]] = u_top - u_bottom
 *
 *          2
 *       _ -*
 *   _ -    |
 *  *-------*
 *  0       1
 */
class FemShapeInterfaceGoodmanTri3 : public FemLib::TemplateShapeFunction<2, 3>
{
public:
    virtual void computeShapeFunction(const double* pt, double* N3)
    {
        double N2[2] = {};
        _fe_line2.computeShapeFunction(pt, N2);

        //*Notice: node order of an element
        //          2
        //       _ -*
        //   _ -    |
        //  *-------*
        //  0       1

        N3[0] = 0.0;
        N3[1] = -N2[1];
        N3[2] = N2[1];
    };

    void ShapeFunctionGoodmanTriEnd(double* N3, const double* pt)
    {
        double N2[2] = {};
        _fe_line2.computeShapeFunction(pt, N2);

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
    FemLib::FemShapeLine3 _fe_line2;
};

}
