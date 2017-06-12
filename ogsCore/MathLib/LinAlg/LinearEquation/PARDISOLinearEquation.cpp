/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PARDISOLinearEquation.cpp
 *
 * Original work is by Chan-Hee Park
 * Moved the work here by Norihiro Watanabe on 2012-08-23
 */

#ifdef USE_PARDISO

#include "PARDISOLinearEquation.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "logog.hpp"

//#if USE_ORGINAL_PARDISO
//#define pardiso_ PARDISO
//#else
//#define PARDISO pardiso_
//#endif

#if defined(USE_ORGINAL_PARDISO)
/* PARDISO prototype. */
extern "C" int PARDISOINIT (void*, int*, int*, int*, double*, int*);
extern "C" int PARDISO (void*, int*, int*, int*, int*, int*,
                        double*, int*, int*, int*, int*, int*,
                        int*, double*, double*, int*, double*);

#elif defined(USE_MKL_PARDISO)
#include "mkl.h"
//#include "mkl_pardiso.h"
//#include "mkl_types.h"
///* PARDISO prototype. */
extern int PARDISO
        (int*, int*, int*, int*, int*, int*,
        double*, int*, int*, int*, int*, int*,
        int*, double*, double*, int*);
#endif

namespace MathLib
{

void PARDISOLinearEquation::solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x)
{

    //omp_set_num_threads (1);
    INFO("------------------------------------------------------------------");
    INFO("*** MKL PARDISO solver computation");
    // Assembling the matrix
    // Establishing CRS type matrix from GeoSys Matrix data storage type
    int nonzero = A->getNNZ();
    int dim = A->getNRows();
    double* value = (double*)A->getEntryArray();
    int* ptr = (int*)malloc((dim + 1) * sizeof( int));
    int* index = (int*)malloc((nonzero) * sizeof( int));

    // Reindexing ptr according to Fortran-based PARDISO
    for(int i = 0; i < dim; ++i)
        ptr[i] = A->getRowPtrArray()[i] + 1;
    //ptr needs one more storage
    ptr[dim] = A->getRowPtrArray()[dim] + 1;
    // Reindexing index according to Fortran-based PARDISO
    // and zonzero of Matrix A
    for(int i = 0; i < nonzero; ++i)
        index[i] = A->getColIdxArray()[i] + 1;

    int mtype = 11;           /* Real unsymmetric matrix */
    int nrhs = 1;             /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void* pt[64];
    /* Pardiso control parameters.*/
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;

    /* Auxiliary variables.*/
    double ddum;              /* Double dummy */
    int idum;                 /* Integer dummy. */

#ifdef USE_ORGINAL_PARDISO
    double dparm[64];
    int solver;
    // Check the license and initialize the solver
    {
        //static bool done = false;
        //if (!done) {
        PARDISOINIT (pt,  &mtype, &solver, iparm, dparm, &error);
        if (error != 0)
        {
            if (error == -10 )
                printf("->No license file found \n");
            if (error == -11 )
                printf("->License is expired \n");
            if (error == -12 )
                printf("->Wrong username or hostname \n");
            exit(1);
        }
        else
            printf("->PARDISO license check was successful ... \n");

        //  done = true;
        //}
    }
#endif

    /* --------------------------------------------------------------------*/
    /* .. Setup Pardiso control parameters.*/
    /* --------------------------------------------------------------------*/
    for (int i = 0; i < 64; i++)
        iparm[i] = 0;
    iparm[0] = 1;             /* No solver default */
    iparm[1] = 2;             /* Fill-in reordering from METIS */
    /* Numbers of processors, value of MKL_NUM_THREADS */
#ifdef USE_ORGINAL_PARDISO
    iparm[2] = omp_get_max_threads();
#else
    iparm[2] = mkl_get_max_threads();
    INFO("-> PARDISO runs with OpenMP threads: %d", iparm[2]);
#endif
    iparm[3] = 0;             /* No iterative-direct algorithm */
    iparm[4] = 0;             /* No user fill-in reducing permutation */
    iparm[5] = 0;             /* Write solution into x */
    iparm[6] = 0;             /* Not in use */
    iparm[7] = 2;             /* Max numbers of iterative refinement steps */
    iparm[8] = 0;             /* Not in use */
    iparm[9] = 13;            /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;            /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;            /* Not in use */
    iparm[12] = 0;            /* Not in use */
    iparm[13] = 0;            /* Output: Number of perturbed pivots */
    iparm[14] = 0;            /* Not in use */
    iparm[15] = 0;            /* Not in use */
    iparm[16] = 0;            /* Not in use */
    iparm[17] = -1;           /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;           /* Output: Mflops for LU factorization */
    iparm[19] = 0;            /* Output: Numbers of CG Iterations */
    maxfct = 1;               /* Maximum number of numerical factorizations. */
    mnum = 1;                 /* Which factorization to use. */
    msglvl = 0;               /* Print statistical information in file */
    error = 0;                /* Initialize error flag */

    /* --------------------------------------------------------------------*/
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* --------------------------------------------------------------------*/
    for (int i = 0; i < 64; i++)
        pt[i] = 0;

    /* --------------------------------------------------------------------*/
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* --------------------------------------------------------------------*/
    phase = 11;
#ifdef USE_ORGINAL_PARDISO
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
#endif

    if (error != 0)
    {
        ERR("***Error: PARDISO returned error code %d during symbolic factorization", error);
        return;
    }

    /* --------------------------------------------------------------------*/
    /* .. Numerical factorization.*/
    /* --------------------------------------------------------------------*/
    phase = 22;
#ifdef USE_ORGINAL_PARDISO
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
#endif
    if (error != 0)
    {
        ERR("***Error: PARDISO returned error code %d during numerical factorization", error);
        return;
    }

    /* --------------------------------------------------------------------*/
    /* .. Back substitution and iterative refinement. */
    /* --------------------------------------------------------------------*/
    phase = 33;
    iparm[7] = 2;             /* Max numbers of iterative refinement steps. */

    /* Set right hand side to one. */

#ifdef USE_ORGINAL_PARDISO
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, b, x, &error, dparm);
#else
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, rhs, x, &error);
#endif
    if (error != 0)
    {
        ERR("***Error: PARDISO returned error code %d during final solution", error);
        return;
    }

    phase = -1;               /* Release internal memory. */
#ifdef USE_ORGINAL_PARDISO
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &dim, value, ptr, index, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
#endif

    // Releasing the local memory
    //delete [] value;
    free(ptr);
    free(index);
    //      MKL_FreeBuffers();
    INFO("------------------------------------------------------------------");
}

}

#endif

