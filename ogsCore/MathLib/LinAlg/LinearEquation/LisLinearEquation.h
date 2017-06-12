/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LisLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "lis.h"

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"
#include "AbstractCRSLinearEquation.h"
#include "LIS_option.h"

namespace MathLib
{

/**
 * \brief Linear equation using Lis solver based on CRS matrix
 */
class LisLinearEquation : public AbstractCRSLinearEquation<signed>
{
public:
    void initialize();
    void finalize();

    virtual ~LisLinearEquation();

    void setOption(const BaseLib::Options &option);

    void setOption(const LIS_option &option)
    {
        _option = option;
    }

    LIS_option &getOption() 
    {
        return _option;
    }

    void gatherX(std::vector<double> &x);

protected:
    void solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x);

private:
    LIS_option _option;
    LIS_VECTOR bb,xx;
};


}
