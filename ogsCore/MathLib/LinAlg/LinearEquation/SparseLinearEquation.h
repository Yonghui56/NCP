/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SparseLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */


#pragma once

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "AbstractCRSLinearEquation.h"

namespace MathLib
{

/**
 * \brief 
 */
class SparseLinearEquation : public AbstractCRSLinearEquation<unsigned>
{
public:
    enum SolverType
    {
        INVALID,
        SolverCG,
        SolverBiCGStab,
        SolverGMRes
    };

    enum PreconditionerType
    {
        NONE,
        DIAG,
        DIAGSCALE
    };

    struct SpLinearOptions
    {
        SolverType solver_type;
        PreconditionerType precon_type;
        double error_tolerance;
        size_t max_iteration_step;

        SpLinearOptions()
        {
            solver_type = SolverCG;
            precon_type = NONE;
            error_tolerance = 1.e-10;
            max_iteration_step = 500;
        };


        SolverType getSolverType(const std::string &str)
        {
            if (str.compare("CG")==0)
                return SolverCG;
            if (str.compare("BICGSTAB")==0)
                return SolverBiCGStab;
            if (str.compare("GMRES")==0)
                return SolverGMRes;

            return INVALID;
        }
        PreconditionerType getPreconType(const std::string &str)
        {
            RETURN_ENUM_IF_SAME_STRING(NONE, str);
            RETURN_ENUM_IF_SAME_STRING(DIAG, str);
            RETURN_ENUM_IF_SAME_STRING(DIAGSCALE, str);

            return NONE;
        }
    };

    SparseLinearEquation() {};

    virtual ~SparseLinearEquation()
    {
    }

    void initialize() {};
    void finalize() {};

    void setOption(const BaseLib::Options &option);

    void setOption(const SpLinearOptions &option)
    {
        _option = option;
    }

    SpLinearOptions &getOption()
    {
        return _option;
    }


protected:
    void solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x);

private:
    SpLinearOptions _option;

    DISALLOW_COPY_AND_ASSIGN(SparseLinearEquation);
};

}
