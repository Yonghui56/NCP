/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestSolvers.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/LinearEquation/DenseLinearEquation.h"
#include "MathLib/LinAlg/LinearEquation/SparseLinearEquation.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquation/LisLinearEquation.h"
#endif
#ifdef USE_PETSC
#include "MathLib/LinAlg/LinearEquation/PETScMPILinearEquation.h"
#endif

#include "TestUtil.h"

using namespace MathLib;


struct Example1
{
    std::vector<double> mat;
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;
    RowMajorSparsity sparse;

    Example1()
    {
        double d_mat[] = {
            6.66667e-012, -1.66667e-012, 0, -1.66667e-012, -3.33333e-012, 0, 0, 0, 0, 
            -1.66667e-012, 1.33333e-011, -1.66667e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 0, 0, 0, 
            0, -1.66667e-012, 6.66667e-012, 0, -3.33333e-012, -1.66667e-012, 0, 0, 0, 
            -1.66667e-012, -3.33333e-012, 0, 1.33333e-011, -3.33333e-012, 0, -1.66667e-012, -3.33333e-012, 0, 
            -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 2.66667e-011, -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 
            0, -3.33333e-012, -1.66667e-012, 0, -3.33333e-012, 1.33333e-011, 0, -3.33333e-012, -1.66667e-012, 
            0, 0, 0, -1.66667e-012, -3.33333e-012, 0, 6.66667e-012, -1.66667e-012, 0, 
            0, 0, 0, -3.33333e-012, -3.33333e-012, -3.33333e-012, -1.66667e-012, 1.33333e-011, -1.66667e-012, 
            0, 0, 0, 0, -3.33333e-012, -1.66667e-012, 0, -1.66667e-012, 6.66667e-012
        };
        mat.assign(d_mat, d_mat+dim_eqs*dim_eqs);
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
        sparse.resize(dim_eqs);
        for (size_t i=0; i<dim_eqs; i++) {
            for (size_t j=0; j<dim_eqs; j++) {
                if (mat[i*dim_eqs+j]!=.0)
                    sparse[i].insert(j);
            }
        }

    }
};

TEST(Math, LinearSolverDirect)
{
    // set problem
    Example1 ex1;

    // construct discrete eqs
    DenseLinearEquation eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            eqs.addA(i, j, v);
        }
    }
    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}

#ifdef USE_BLAS_LAPACK
TEST(Math, LinearSolverSparse)
{
    // set problem
    Example1 ex1;

    // construct discrete eqs
    SparseLinearEquation eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);
    eqs.getOption().solver_type = SparseLinearEquation::SolverCG;
    eqs.getOption().precon_type = SparseLinearEquation::NONE;

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            if (v!=.0)
                eqs.addA(i, j, v);
        }
    }
    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}
#endif

#ifdef USE_LIS
TEST(Math, LinearSolverLis1)
{
#if 0
    {
        LIS_MATRIX        A;
        LIS_VECTOR        b,x,u;
        LIS_SOLVER        solver;
        int              nprocs,my_rank;
        int               i,n,gn,is,ie,iter;
        n  = 12;
        int argc=0;
        char **argv;
        lis_initialize(&argc, &argv);

    #ifdef USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    #else
        nprocs  = 1;
        my_rank = 0;
    #endif

        lis_matrix_create(LIS_COMM_WORLD,&A);
        lis_matrix_set_size(A,0,n);
        lis_matrix_get_size(A,&n,&gn);
        lis_matrix_get_range(A,&is,&ie);
        for(i=is;i<ie;i++)
        {
            if( i>0   )  lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if( i<gn-1 ) lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0,A);
        }
        lis_matrix_set_type(A,LIS_MATRIX_CRS);
        lis_matrix_assemble(A);

        lis_vector_duplicate(A,&u);
        lis_vector_duplicate(A,&b);
        lis_vector_duplicate(A,&x);
        lis_vector_set_all(1.0,u);
        lis_matvec(A,u,b);
        lis_solver_create(&solver);
        lis_solver_set_option("-print mem",solver);
        lis_solver_set_optionC(solver);
        lis_solve(A,b,x,solver);
        lis_solver_get_iters(solver,&iter);
        if (my_rank==0)
          {
        printf("iter = %d\n",iter);
        printf("\n");
          }
        lis_vector_print(x);

        lis_matrix_destroy(A);
        lis_vector_destroy(b);
        lis_vector_destroy(x);
        lis_vector_destroy(u);
        lis_solver_destroy(solver);
        lis_finalize();

    }
#endif
    // set problem
    Example1 ex1;

    // construct discrete eqs
    LisLinearEquation eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);
    eqs.getOption().ls_method = LIS_option::CG;
    eqs.getOption().ls_precond = LIS_option::NONE;

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            if (v!=.0)
                eqs.addA(i, j, v);
        }
    }

    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}

#endif

#if defined(USE_PETSC) && defined(USE_MPI)
TEST(Math, LinearSolverPETsc)
{
    // set problem
    Example1 ex1;
    std::vector<double> eqsX(ex1.dim_eqs, .0);
    {

        //
        char help[] = "Using PETSc package\n";
        int argc = 1;
        char exename[] = "test\n";
        char *argv[] = {"test\n"};
        char **tmp = argv;
        PetscInitialize(&argc, &tmp,(char *)0,help);

        int rank_p, size_p;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank_p);
        MPI_Comm_size(PETSC_COMM_WORLD, &size_p);


        PETScLinearSolver eqs(ex1.dim_eqs);
        eqs.Init();
        eqs.set_rank_size(rank_p, size_p);

        eqs.Config(1e-10, KSPCG,PCILU);
        //eqs.Config(1e-10, KSPCG,PCILU);
        eqs.Initialize(); 

        for (int i=0; i<ex1.dim_eqs; i++)
            for (int j=0; j<ex1.dim_eqs; j++)
                eqs.addMatrixEntry(i, j, ex1.mat[i*ex1.dim_eqs+j]);

        eqs.AssembleMatrixPETSc();
        eqs.AssembleRHS_PETSc();
        eqs.AssembleUnkowns_PETSc();

        size_t nrows = ex1.list_dirichlet_bc_value.size();   
        PetscInt *rows_toberemoved = new PetscInt[nrows];

        for (size_t i=0; i<nrows; i++) 
        {
            int id = ex1.list_dirichlet_bc_id[i];
            double val = ex1.list_dirichlet_bc_value[i];
            rows_toberemoved[i] = id;
            eqs.set_bVectorEntry(id, val);         
            eqs.set_xVectorEntry(id, val);   
        }

        eqs.zeroRows_in_Matrix(nrows, rows_toberemoved);

        eqs.AssembleRHS_PETSc();
        eqs.AssembleUnkowns_PETSc();
        eqs.AssembleMatrixPETSc();

        eqs.Solver();

        eqs.UpdateSolutions(&eqsX[0], &eqsX[0]);

        delete [] rows_toberemoved;
        rows_toberemoved = NULL;  
    }

    PetscFinalize();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], &eqsX[0], ex1.dim_eqs, 1.e-5);

}
#endif
