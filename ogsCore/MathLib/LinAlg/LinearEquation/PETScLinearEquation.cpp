/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PETScLinearEquation.cpp
 *
 * Created on 2012-08-24 by Norihiro Watanabe
 */

#include "PETScLinearEquation.h"

#include <iostream>

namespace MathLib
{

void PETScLinearEquation::initialize()
{
    char help[] = "Using PETSc package\n";
    int argc = 1;
    char *argv[] = {"test\n"};
    char **tmp = argv;
    PetscInitialize(&argc, &tmp,(char *)0,help);
}

void PETScLinearEquation::finalize()
{
    PetscFinalize();
}

PETScLinearEquation::~PETScLinearEquation()
{
}

void PETScLinearEquation::setOption(const BaseLib::Options &option)
{
    const BaseLib::Options *op = option.getSubGroup("LinearSolver");
    if (op==0) {
        return;
    }

    if (op->hasOption("solver_type"))
        _option.solver_type = _option.getSolverType(op->getOption("solver_type"));
    if (op->hasOption("precon_type"))
        _option.precon_type = _option.getPreconType(op->getOption("precon_type"));
    if (op->hasOption("error_tolerance"))
        _option.ls_error_tolerance = op->getOptionAsNum<double>("error_tolerance");
    if (op->hasOption("max_iteration_step"))
        _option.ls_max_iterations = op->getOption<int>("max_iteration_step");
}

void PETScLinearEquation::solveEqs(CRSMatrix<double, signed> *A, double *b, double *x)
{
    long dimension = static_cast<long>(A->getNRows());

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "*** PETSc solver computation" << std::endl;

//    std::cout << "A=" << std::endl;
//    A->printMat();
//    std::cout << "b=" << std::endl;
//    for (size_t i=0; i<A->getNRows(); i++)
//        std::cout << b[i] << " ";
//    std::cout << std::endl;

    const int dim = A->getNRows();

    //-------------------------------------------------------------------------
    // create vector and matrix
    VecCreate(PETSC_COMM_WORLD, &bb);
    VecSetSizes(bb, PETSC_DECIDE, dim);
    VecSetFromOptions(bb);
    VecDuplicate(bb, &xx);

    MatCreate(PETSC_COMM_WORLD, &AA);
    MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,dim,dim);
    MatSetFromOptions(AA);
  #ifdef USE_MPI
    MatSetType(AA,MATMPIAIJ);
    MatMPIAIJSetPreallocation(AA,d_nz,PETSC_NULL, o_nz,PETSC_NULL);
  #else
    MatSetType(AA,MATSEQAIJ);
  #endif
    MatSetUp(AA);
    MatGetOwnershipRange(AA,&i_start,&i_end);

    //-------------------------------------------------------------------------
    // set vector and matrix
    for (int i=0; i < dimension; ++i)
    {
        VecSetValues(bb,1,&i,&b[i], INSERT_VALUES);
        VecSetValues(xx,1,&i,&x[i], INSERT_VALUES);
    }
    VecAssemblyBegin(bb);
    VecAssemblyEnd(bb);
    VecAssemblyBegin(xx);
    VecAssemblyEnd(xx);

    //TODO can we directory pass pointer instead of copy all data?
    const int* crs_row_ptr = A->getRowPtrArray();
    const int* crs_col_id = A->getColIdxArray();
    const double* crs_data = A->getEntryArray();
    int k = 0;
    for (int i=0; i<dim; i++) {
        const int row_end = crs_row_ptr[i+1];
        for (int j=crs_row_ptr[i]; j<row_end; j++) {
            double val = crs_data[k++];
            const int col = crs_col_id[j];
            MatSetValue(AA, i, col , val, ADD_VALUES);
        }
    }
    MatAssemblyBegin(AA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(AA, MAT_FINAL_ASSEMBLY);

    //-------------------------------------------------------------------------
    // create solver
    KSPCreate(PETSC_COMM_WORLD, &petsc_solver);
    KSPSetOperators(petsc_solver, AA, AA, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(petsc_solver, _option.solver_type);
    KSPGetPC(petsc_solver, &petsc_prec);
    PCSetType(petsc_prec, _option.precon_type);
    KSPSetTolerances(petsc_solver, _option.ls_error_tolerance, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(petsc_solver);

    //-------------------------------------------------------------------------
    // solve
    PetscLogDouble v1,v2;

    PetscGetTime(&v1);

    KSPSolve(petsc_solver, bb, xx);

    PetscGetTime(&v2);

    //EQSV_Viewer("petsc.out.log");

    //-------------------------------------------------------------------------
    // check result
    PetscPrintf(PETSC_COMM_WORLD, "solver    : %s\n", _option.solver_type);
    PetscPrintf(PETSC_COMM_WORLD, "precon    : %s\n", _option.precon_type );
    KSPConvergedReason reason;
    KSPGetConvergedReason(petsc_solver,&reason); //CHKERRQ(ierr);
    if (reason==KSP_DIVERGED_INDEFINITE_PC)
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : Divergence because of indefinite preconditioner. Run the executable again but with -pc_factor_shift_positive_definite option.");
    }
    else if (reason<0)
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : Other kind of divergence: this should not happen.\n");
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : converged\n");
        int its;
        double res;
        KSPGetIterationNumber(petsc_solver,&its); //CHKERRQ(ierr);
        KSPGetResidualNorm(petsc_solver, &res);
        PetscPrintf(PETSC_COMM_WORLD, "iteration : %d/%d\n", (int)its, _option.ls_max_iterations);
        PetscPrintf(PETSC_COMM_WORLD, "residual  : %e\n", res);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");

//    VecAssemblyBegin(xx);
//    VecAssemblyEnd(xx);

    //-------------------------------------------------------------------------
    // Update the solution (answer) into the x vector
     PetscScalar *xp;
//     PetscInt low,high,otherlow;
//     VecGetOwnershipRange(xx, &low, &high);

     VecGetArray(xx, &xp);
     PetscInt count;
     VecGetLocalSize(xx, &count);
     for(int i=0; i<count; i++) //TODO this shouldn't work with DDC
       x[i] = xp[i];


    // Clear memory
     VecDestroy(&bb);
     VecDestroy(&xx);
     MatDestroy(&AA);
     KSPDestroy(&petsc_solver);
//     PCDestroy(&petsc_prec);

     std::cout << "------------------------------------------------------------------" << std::endl;
}


}
