      implicit none
      
#include "lisf.h"

      LIS_INTEGER      my_rank,nprocs
      LIS_INTEGER      matrix_type,comm_world
      LIS_INTEGER      omp_get_num_procs,omp_get_max_threads
      LIS_INTEGER      i,n,gn,ln,is,ie,iter,ierr
      LIS_MATRIX       A
      LIS_VECTOR       b,x,u
      LIS_SOLVER       solver
      
      call lis_initialize(ierr)

#ifdef USE_MPI
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      n = 12
      ln = 0
      matrix_type = LIS_MATRIX_CRS
      comm_world = LIS_COMM_WORLD

      call lis_matrix_create(comm_world,A,ierr)
      call lis_matrix_set_size(A,ln,n,ierr)
      call lis_matrix_get_size(A,n,gn,ierr)
      call lis_matrix_get_range(A,is,ie,ierr)
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0d0,
     .                                        A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0d0,
     .                                        A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0d0,A,ierr)
      enddo
      call lis_matrix_set_type(A,matrix_type,ierr)
      call lis_matrix_assemble(A,ierr)
      call lis_vector_duplicate(A,u,ierr)
      call lis_vector_duplicate(A,b,ierr)
      call lis_vector_duplicate(A,x,ierr)
      call lis_vector_set_all(1.0d0,u,ierr)
      call lis_matvec(A,u,b,ierr)
      
      call lis_solver_create(solver,ierr)
      call lis_solver_set_option("-print mem",solver,ierr);
      call lis_solver_set_optionC(solver,ierr)
      call lis_solve(A,b,x,solver,ierr)
      call lis_solver_get_iters(solver,iter,ierr)
      call lis_vector_print(x,ierr)
      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(b,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_vector_destroy(u,ierr)
      call lis_solver_destroy(solver,ierr)
      call lis_finalize(ierr)

      stop
      end
      
