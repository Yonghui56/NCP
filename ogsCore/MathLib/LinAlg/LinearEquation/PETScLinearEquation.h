/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PETScLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#ifdef USE_PETSC

#include <string>

#include "petscmat.h"
#include "petscksp.h"

#include "BaseLib/Options.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"
#include "AbstractCRSLinearEquation.h"


namespace MathLib
{

struct PETSc_option
{
    typedef const char* SolverType;
    typedef const char* PreconType;
    SolverType solver_type;
    PreconType precon_type;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;

    PETSc_option()
    {
        solver_type = KSPBCGS;
        precon_type = PCNONE;
        ls_max_iterations = 500;
        ls_error_tolerance = 1.e-6;
    }

    SolverType getSolverType(const std::string &str)
    {
        if (str.compare("CG")==0)
            return KSPCG;
        if (str.compare("BICG")==0)
            return KSPBICG;
        if (str.compare("BICGSTAB")==0)
            return KSPBCGS;

        return "";
    }
    PreconType getPreconType(const std::string &str)
    {
        if (str.compare("NONE")==0)
            return PCNONE;
        if (str.compare("JACOBI")==0)
            return PCJACOBI;
        if (str.compare("ILU")==0)
            return PCILU;

        return PCNONE;
    }
};

class PETScLinearEquation : public AbstractCRSLinearEquation<signed>
{
public:
    virtual ~PETScLinearEquation();

    void initialize();
    void finalize();

    void setOption(const BaseLib::Options &option);
    void setOption(const PETSc_option &option)
    {
        _option = option;
    }
    PETSc_option &getOption() 
    {
        return _option;
    }


protected:
    void solveEqs(CRSMatrix<double, signed> *A, double *b, double *x);

private:
    PETSc_option _option;
    Vec bb;
    Vec xx;
    Mat  AA;
    KSP petsc_solver;
    PC petsc_prec;
    PetscInt i_start;
    PetscInt i_end;
};

} //end

#endif
