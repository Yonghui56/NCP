/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIntegrationGaussTetra.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "logog.hpp"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

/**
 * \brief Gauss quadrature rule for tetrahedrals
 *
 */
class FemIntegrationGaussTetra : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        const size_t n = getNumberOfSamplingPoints();
        switch (n)
        {
        case 1:
            getSamplePointTet1(igp, x);
             break;
        case 5:
            getSamplePointTet5(igp, x);
            break;
        case 15:
            getSamplePointTet15(igp, x);
            break;
        default:
            ERR("***Error in FemIntegrationGaussTetra::getSamplingPoint() - Given sampling points %d are not supported not supported.", n);
            break;
        }
    }

    double getWeight(size_t igp) const
    {
        double w = .0;
        switch (getNumberOfSamplingPoints())
        {
        case 1:
            w= 1.0;
            break;
        case 5:
            switch(igp)
            {
            case 0:
                w = -0.133333333333333;
                break;
            case 1:
            case 2:
            case 3:
            case 4:
                w = 0.07500000000000;
                break;
            default:
                break;
            }

            break;
        case 15:
            switch(igp)
            {
            case 0:
                w = 0.019753086419753086;
                break;
            case 1:
            case 2:
            case 3:
            case 4:
                w = 0.011989513963169772;
                break;
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
                w = 0.008818342151675485;
                break;
            case 11:
            case 12:
            case 13:
            case 14:
                w = 0.011511367871045397;
                break;
            default:
                break;
            }
            break;
        }

        return w;
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        if (n_sampl_level==1) return 1;
        else if (n_sampl_level==2) return 5;
        else if (n_sampl_level==3) return 15;
        else return 0;
    }

    void getSamplePointTet1(size_t, double *pt) const
    {
        pt[0] = 0.333333333333333 ;
        pt[1] = 0.333333333333333 ;
        pt[2] = 0.333333333333333 ;
    }

    void getSamplePointTet5(size_t igp, double *x) const
    {
        switch(igp)
        {
        case 0:
            x[0] = 0.25;
            x[1] = 0.25;
            x[2] = 0.25;
            break;
        case 1:
            x[0] = 0.166666666666667;
            x[1] = 0.166666666666667;
            x[2] = 0.166666666666667;
            break;
        case 2:
            x[0] = 0.5;
            x[1] = 0.166666666666667;
            x[2] = 0.166666666666667;
            break;
        case 3:
            x[0] = 0.166666666666667;
            x[1] = 0.5;
            x[2] = 0.166666666666667;
            break;
        case 4:
            x[0] = 0.166666666666667;
            x[1] = 0.166666666666667;
            x[2] = 0.5;
            break;
        default: break;
        }
    }

    void getSamplePointTet15(size_t igp, double* x) const
    {
        switch(igp)
        {
        case 0:
            x[0] = 0.25;
            x[1] = 0.25;
            x[2] = 0.25;
            break;
        case 1:
            x[0] = 0.09197107805272303;
            x[1] = 0.09197107805272303;
            x[2] = 0.09197107805272303;
            break;
        case 2:
            x[0] = 0.72408676584183096;
            x[1] = 0.09197107805272303;
            x[2] = 0.09197107805272303;
            break;
        case 3:
            x[0] = 0.09197107805272303;
            x[1] = 0.72408676584183096;
            x[2] = 0.09197107805272303;
            break;
        case 4:
            x[0] = 0.09197107805272303;
            x[1] = 0.09197107805272303;
            x[2] = 0.72408676584183096;
            break;
        case 5:
            x[0] = 0.44364916731037080;
            x[1] = 0.05635083268962915;
            x[2] = 0.05635083268962915;
            break;
        case 6:
            x[0] = 0.05635083268962915;
            x[1] = 0.44364916731037080;
            x[2] = 0.05635083268962915;
            break;
        case 7:
            x[0] = 0.05635083268962915;
            x[1] = 0.05635083268962915;
            x[2] = 0.44364916731037080;
            break;
        case 8:
            x[0] = 0.05635083268962915;
            x[1] = 0.44364916731037080;
            x[2] = 0.44364916731037080;
            break;
        case 9:
            x[0] = 0.44364916731037080;
            x[1] = 0.05635083268962915;
            x[2] = 0.44364916731037080;
            break;
        case 10:
            x[0] = 0.44364916731037080;
            x[1] = 0.44364916731037080;
            x[2] = 0.05635083268962915;
            break;
        case 11:
            x[0] = 0.31979362782962989;
            x[1] = 0.31979362782962989;
            x[2] = 0.31979362782962989;
            break;
        case 12:
            x[0] = 0.04061911651111023;
            x[1] = 0.31979362782962989;
            x[2] = 0.31979362782962989;
            break;
        case 13:
            x[0] = 0.31979362782962989;
            x[1] = 0.04061911651111023;
            x[2] = 0.31979362782962989;
            break;
        case 14:
            x[0] = 0.31979362782962989;
            x[1] = 0.31979362782962989;
            x[2] = 0.04061911651111023;
            break;
        default: break;
        }
    }
};

}
