/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFemShapeFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IFemShapeFunction.h"

namespace FemLib
{

template <size_t N_DIM, size_t N_NODES>
class TemplateShapeFunction : public IFemShapeFunction
{
public:
    virtual void computeShapeFunction(const double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(const double* pt, double* dN) = 0;
    virtual ~TemplateShapeFunction() {};
protected:
};
  
}
