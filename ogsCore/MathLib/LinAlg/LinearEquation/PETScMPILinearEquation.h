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


namespace MathLib
{

typedef struct {
    int ls_method;
    int ls_precond;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;
} PETSc_option;


typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;

class PETScLinearSolver
{
  public:    
    PETScLinearSolver (const int size)
      :lsolver(NULL),prec(NULL) 
    {
       ltolerance = 1.e-10;
       m_size = size;
       time_elapsed = 0.0;           
    }
    ~PETScLinearSolver();

    void Config(const float tol, const KSPType lsol, const PCType prec_type);

    void Init()
    {

      VectorCreate(m_size);   
      MatrixCreate(m_size, m_size);
    }

    void Solver();
    void AssembleRHS_PETSc();
    void AssembleUnkowns_PETSc( );
    void AssembleMatrixPETSc();

    void UpdateSolutions(PetscScalar *u0, PetscScalar *u1);


    int Size() const {return m_size;}

    void set_xVectorEntry(const int i, const double value);
    void set_bVectorEntry(const int i, const double value);

    void add_xVectorEntry(const int i, const double value, InsertMode mode);
    void add_bVectorEntry(const int i, const double value, InsertMode mode);
    void addMatrixEntry(const int i, const int j, const double value);
 
    void Initialize();

    void zeroRows_in_Matrix(const int nrow, const  PetscInt *rows);
    void zeroMatrix()
    {
       MatZeroEntries(A);
    }


    void set_rank_size(const int m_rank, const int size)
    {
      mpi_size = size;
      rank = m_rank;  
    } 


    PetscInt getStartRow() const {return i_start;} 
    PetscInt getEndRow() const {return i_end;} 

    void EQSV_Viewer(std::string file_name, PetscViewer viewer);
    void EQSV_Viewer(std::string file_name);
   
private:
    PETSc_Mat  A;
    PETSc_Vec b;
    PETSc_Vec x;
    KSP lsolver;
    PC prec; 
    PetscInt i_start;
    PetscInt i_end;

    // Slover and preconditioner names, only for log
    std::string sol_type;
    std::string pc_type;

    PetscLogDouble time_elapsed;

    long m_size;
    float ltolerance;  

    int mpi_size;
    int rank;

    void VectorCreate(PetscInt m);
    void MatrixCreate(PetscInt m, PetscInt n);

};

} //end

#endif
