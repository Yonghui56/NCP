
#include "LisMPILinearEquation.h"

namespace MathLib
{

LisMPILinearEquation::~LisMPILinearEquation()
{
    lis_matrix_destroy(_A);
    lis_vector_destroy(_b);
    lis_vector_destroy(_x);
}

void LisMPILinearEquation::createDynamic(size_t local_n, size_t global_n)
{
    _dynamic = true;
    int err = 0;
    lis_matrix_create(LIS_COMM_WORLD, &_A);
    lis_matrix_set_size(_A,local_n,global_n);
    lis_matrix_get_size(_A, &_local_dim, &_global_dim);
    err = lis_matrix_get_range(_A,&_is,&_ie); CHKERR(err);
    err = lis_vector_duplicate(_A,&_b); CHKERR(err);
    err = lis_vector_duplicate(_b,&_x); CHKERR(err);

    reset();
}

void LisMPILinearEquation::assembleMatrix()
{
    int err;
    if (_dynamic) {
        err = lis_matrix_set_type(_A, LIS_MATRIX_CRS);
    } else {
        err = lis_matrix_set_crs(_crs.row_ptr[_ie-_is],_crs.row_ptr,_crs.col_idx,_crs.data,_A); CHKERR(err);
    }
    err = lis_matrix_assemble(_A); CHKERR(err);
}

void LisMPILinearEquation::setOption(const Base::Options &option)
{
    throw "LisMPILinearEquation::setOption() is not implemented.";
}

void LisMPILinearEquation::reset()
{
    lis_vector_set_all(0., _b);
    lis_vector_set_all(0., _x);
}

double LisMPILinearEquation::getA(size_t rowId, size_t colId)
{
    throw "not implemented.";
}

void LisMPILinearEquation::setA(size_t rowId, size_t colId, double v)
{
    lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _A);
}

void LisMPILinearEquation::addA(size_t rowId, size_t colId, double v)
{
    lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _A);
}

void LisMPILinearEquation::addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt)
{
    for (size_t i=0; i<vec_row_pos.size(); i++) {
        const size_t rowId = vec_row_pos[i];
        for (size_t j=0; j<vec_col_pos.size(); j++) {
            const size_t colId = vec_col_pos[j];
            addA(rowId, colId, fkt*sub_matrix(i,j));
        }
    }
}

void LisMPILinearEquation::addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt)
{
    addA(vec_pos, vec_pos, sub_matrix, fkt);
}

double LisMPILinearEquation::getRHS(size_t rowId)
{
    double v;
    lis_vector_get_value(_b, rowId, &v);
    return v;
}

double* LisMPILinearEquation::getRHS()
{
    return &_tmp_b[0];
}

void LisMPILinearEquation::setRHS(size_t rowId, double v)
{
    lis_vector_set_value(LIS_INS_VALUE, rowId, v, _b);
}

void LisMPILinearEquation::addRHS(size_t rowId, double v)
{
    lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _b);
}

void LisMPILinearEquation::addRHS(std::vector<size_t> &vec_pos, double *sub_vector, double fkt)
{
    for (size_t i=0; i<vec_pos.size(); i++) {
        const size_t rowId = vec_pos[i];
        addRHS(rowId, sub_vector[i]*fkt);
    }
}

double* LisMPILinearEquation::getX()
{
    return &_tmp_x[0];
}

void LisMPILinearEquation::setKnownX(size_t row_id, double x)
{

}

void LisMPILinearEquation::setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
{

}

void LisMPILinearEquation::solve()
{
    int err;

    // solve
    LIS_SOLVER solver;
    err = lis_solver_create(&solver); CHKERR(err);
    lis_solver_set_option("-print mem",solver);
    lis_solver_set_optionC(solver);

    err = lis_solve(_A,_b,_x,solver); CHKERR(err);

    int    iter,iter_double,iter_quad;
    double            times,itimes,ptimes,p_c_times,p_i_times;
    LIS_REAL        resid;
    int                nsol;
    char            solvername[128];
    lis_solver_get_itersex(solver,&iter,&iter_double,&iter_quad);
    lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
    lis_solver_get_residualnorm(solver,&resid);
    lis_solver_get_solver(solver,&nsol);
    lis_get_solvername(nsol,solvername);
    
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if( myrank==0 )
    {
        printf("%s: iter     = %d iter_double = %d iter_quad = %d\n",solvername,iter, iter_double, iter_quad);
        printf("%s: times    = %e\n",solvername,times);
        printf("%s: p_times  = %e (p_c = %e p_i = %e )\n",solvername, ptimes, p_c_times,p_i_times);
        printf("%s: i_times  = %e\n",solvername, itimes);
        printf("%s: Residual = %e\n\n",solvername,resid);
    }

    lis_solver_destroy(solver);

#if 0
    int iter = 0;
    double resid = 0.0;

    //assemble matrix
    lis_matrix_set_type(_A, LIS_MATRIX_CRS);
    lis_matrix_assemble(_A);




    // Matrix solver and Precondition can be handled better way.
    const size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

    int nthreads = omp_get_num_threads();
    //omp_set_num_threads (nthreads);

    sprintf(solver_options, "-i %d -p %d %s", _option.ls_method, _option.ls_precond, _option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d -omp_num_threads %d -initx_zeros 0", _option.ls_error_tolerance, _option.ls_max_iterations, nthreads);



    // Create solver
    LIS_SOLVER solver;
    lis_solver_create(&solver);
    lis_solver_set_option(solver_options, solver);
    lis_solver_set_option(tol_option, solver);
    lis_solver_set_option("-print mem", solver);

    // solve
    lis_solve(_A, _b, _x, solver);
    _tmp_x.resize(_global_dim);
    lis_vector_gather(_x, &_tmp_x[0]);

    //
    lis_solver_get_iters(solver, &iter);
    lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%d\n", iter, _option.ls_max_iterations);
    printf("\t residuals: %e\n", resid);

    lis_solver_destroy(solver);
#endif
}

} //end
