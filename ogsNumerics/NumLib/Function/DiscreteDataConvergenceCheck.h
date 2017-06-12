/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DiscreteDataConvergenceCheck.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "MathLib/Vector.h"
#include "NumLib/Function/NormOfDiscreteDataFunction.h"
#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

namespace NumLib
{

/**
 * \brief Convergence check for iterative calculation of discrete data
 *
 */
class DiscreteDataConvergenceCheck : public NumLib::IConvergenceCheck
{
public:
    ///
    DiscreteDataConvergenceCheck() { };

    ///
    virtual ~DiscreteDataConvergenceCheck() {};

    ///
    virtual bool isConverged(NumLib::UnnamedParameterSet& vars_prev, NumLib::UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {

        for (size_t i=0; i<vars_prev.size(); i++) {
            v_diff = .0;
            bool isScalar = (vars_prev.getName(i).compare("h")==0 || vars_prev.getName(i).compare("c")==0);
            bool isIntegrationPointVector = (vars_prev.getName(i).compare("v")==0);

            if (isScalar) {
                v_diff = calculateDifference<double>(vars_prev, vars_current, i);
            } else if (isIntegrationPointVector) {
                v_diff = calculateDifference<MathLib::TemplateVectorX<MathLib::LocalVector> >(vars_prev, vars_current, i);
            }

            if (v_diff>eps) {
                return false;
            }
        }
        return true;
    }

    ///
    template <class T>
    double calculateDifference(NumLib::UnnamedParameterSet &vars_prev, NumLib::UnnamedParameterSet &vars_current, size_t i)
    {
        const NumLib::ITXDiscreteFunction<T>* f_fem_prev = vars_prev.get<NumLib::ITXDiscreteFunction<T> >(i);
        const NumLib::ITXDiscreteFunction<T>* f_fem_cur = vars_current.get<NumLib::ITXDiscreteFunction<T> >(i);
        NumLib::NormOfDiscreteDataFunction<T> _norm;
        return _norm(*f_fem_prev, *f_fem_cur);
    }
};

}
