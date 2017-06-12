      implicit none
      
#include "lisf.h"

      LIS_MATRIX    A
      LIS_VECTOR    x
      LIS_REAL      evalue0
      LIS_ESOLVER   esolver
      LIS_INTEGER   ierr
      LIS_INTEGER   my_rank,nprocs
      LIS_INTEGER   matrix_type,comm_world 
      LIS_INTEGER   omp_get_num_procs,omp_get_max_threads
      LIS_INTEGER   n,gn,i
      LIS_INTEGER   nsol,iter,iter_double,iter_quad
      real*8        times,itimes,ptimes,p_c_times,p_i_times
      LIS_REAL      resid
      character*256 fname,solname,resname
      character*20  esolvername
      integer*4     iargc
      
      call lis_initialize(ierr) 

#ifdef USE_MPI
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      matrix_type = LIS_MATRIX_CRS
      comm_world = LIS_COMM_WORLD

      i = iargc()

      if( i.lt.3 ) then
        if( my_rank.eq.0 ) then
          write(*,*) 'etest1f matrix_filename ',
     .              'solution_filename residual_filename [options]'
          call lis_finalize(ierr)
        endif
        stop
      endif

      if (my_rank .eq. 0) then
         write(*,*) ''
         write(*,*) 'number of processes = ',nprocs
      endif


#ifdef _OPENMP
      write(*,*) 'max number of threads = ',omp_get_num_procs()
      write(*,*) 'number of threads = ', omp_get_max_threads()
#endif

! read matrix from file 
      call getarg(1,fname)
      call lis_matrix_create(comm_world,A,ierr)
      call CHKERR(ierr)
      call lis_matrix_set_type(A,matrix_type,ierr)
      call lis_input_matrix(A,fname,ierr)
      call CHKERR(ierr);
      call lis_matrix_get_size(A,n,gn,ierr)

      call lis_vector_duplicate(A,x,ierr)
      call lis_vector_set_all(1.0d0,x,ierr)

      call lis_esolver_create(esolver,ierr)
      call CHKERR(ierr)
      call lis_esolver_set_option('-eprint mem',esolver,ierr)
      call lis_esolver_set_optionC(esolver,ierr)
      
      call lis_esolve(A,x,evalue0,esolver,ierr)
      call CHKERR(ierr)

      call lis_esolver_get_itersex(esolver,iter,iter_double,iter_quad,
     .     ierr)
      call lis_esolver_get_timeex(esolver,times,itimes,ptimes,
     .     p_c_times,p_i_times,ierr)
      call lis_esolver_get_residualnorm(esolver,resid,ierr)
      call lis_esolver_get_esolver(esolver,nsol,ierr)
      call lis_get_esolvername(nsol,esolvername,ierr)
	
      If( my_rank.eq.0 ) then
        write(*,*) esolvername,': mode number              = ',0
        write(*,*) esolvername,': eigenvalue               = ',evalue0
        write(*,*) esolvername,': number of iterations     = ',iter
        write(*,*) esolvername,': elapsed time             = ',times
        write(*,*) esolvername,':   preconditioner         = ',ptimes
        write(*,*) esolvername,':     matrix creation      = ',p_c_times
        write(*,*) esolvername,': linear solver            = ',itimes
        write(*,*) esolvername,': relative residual 2-norm = ',resid
        write(*,*) 
      endif

! write solution 
      call getarg(2,solname)
      call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

! write residual 
      call getarg(3,resname)
      call lis_esolver_output_rhistory(esolver, resname,ierr);

      call lis_esolver_destroy(esolver,ierr)
      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(x,ierr)

      call lis_finalize(ierr)

      stop
      end
      
