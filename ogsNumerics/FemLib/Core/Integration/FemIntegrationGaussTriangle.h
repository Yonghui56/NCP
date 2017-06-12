/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationGaussTriangle.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "logog.hpp"

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

/**
 * \brief Gauss quadrature rule for triangles
 *
 * Gauss quadrature rule for triangles is originally given as
 * \f[
 *    \int F(x,y) dx dy = \int F(x(r, s), y(r, s)) j(r,s) dr ds \approx \frac{1}{2} \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w_i )
 * \f]
 *
 * To make it consistent with other elements, we rewrite the above formula as
 * \f[
 *    \int F(x,y) dx dy \approx \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w'_i )
 * \f]
 * by defining the new weight \f$ w'=\frac{1}{2} w \f$.
 */
class FemIntegrationGaussTriangle : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        const size_t n = getNumberOfSamplingPoints();
        switch (n)
        {
        case 1:
            getSamplePointTri1(igp, x);
             break;
        case 3:
            getSamplePointTri3(igp, x);
            break;
        default:
            ERR("***Error in FemIntegrationGaussTriangle::getSamplingPoint() - Given sampling points %d are not supported not supported.", n);
            break;
        }
    }

    double getWeight(size_t) const
    {
        double w = .0;
        switch (getNumberOfSamplingPoints())
        {
        case 1:
            w= 1.0;
            break;
        case 3:
            w= 0.333333333333333; // = 1/3
            break;
        }

        return w*0.5;
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        if (n_sampl_level==1) return 1;
        else if (n_sampl_level==2) return 3;
        else return 6;
    }

    void getSamplePointTri1(size_t, double *pt) const
    {
        pt[0] = 0.333333333333333 ;
        pt[1] = 0.333333333333333 ;
    }

    void getSamplePointTri3(size_t igp, double *pt) const
    {
        switch(igp)
        {
        case 0:
            pt[0] = 0.166666666666667 ;
            pt[1] = 0.166666666666667 ;
            break;
        case 1:
            pt[0] = 0.666666666666667 ;
            pt[1] = 0.166666666666667 ;
            break;
        case 2:
            pt[0] = 0.166666666666667 ;
            pt[1] = 0.666666666666667 ;
            break;
        }
    }
};

}
