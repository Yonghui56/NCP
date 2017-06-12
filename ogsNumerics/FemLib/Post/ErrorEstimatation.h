/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ErrorEstimatation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace FemLib
{
/**
 * \brief Error estimator
 */    
class IErrorEstimator {};

/**
 * \brief Zienkiewicz-Zhu error estimator 
 */
class FemErrorZienkiewiczZhu : IErrorEstimator
{
public:
    void estimate(double* duh, double* dwh)
    {

    }
};

}
