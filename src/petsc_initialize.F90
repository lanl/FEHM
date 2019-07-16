!******************************************************************************
!              Initialize PETSc Solver, mainly for memory
!******************************************************************************
!
module petsc_initialize_package
    
#include <petsc/finclude/petscksp.h>            

    use petscksp
    use comco2

    public :: petsc_initialize
    
    Vec    x,b, nnz_seq
    Mat    A_matrix
    KSP    ksp
    PC     pc

    PetscInt,dimension(:),allocatable :: bp_index,x_index        ! index for right hand vector
    PetscReal,dimension(:),allocatable :: solution,sol_gather    ! index for right hand vector
    PetscReal,dimension(:),allocatable :: value                  ! value for Assemble A matrix
    PetscInt,dimension(:),allocatable :: recvcounts,displs       ! index for MPI_Gatherv

    PetscBool        flg
    PetscErrorCode   ierr_3

    PetscInt         A_size,i,j, nnz_start, nnz_end, nnz_size
    PetscInt,dimension(:),allocatable :: nnz,nnz_above            ! nnz: number of nonzero on each row
    PetscInt,dimension(:),allocatable :: nnz_global               ! nnz: number of nonzero on each row

    integer :: nnz_total,nnz_max, nnz_sum 
    integer :: Istart,Iend,rank_nonzero_num
    integer :: d_nz,o_nz, idf, neq_inc, idf_cnt, neq_adj, idf_adj   

    integer,dimension(:),allocatable :: col_id, test_id                    ! for assemble A_matrix
    integer,dimension(:),allocatable :: row_index, row_nelm       ! for assemble A_matrix
    integer,dimension(:),allocatable :: value_index               ! the index for value vector to assemble
    integer,dimension(:),allocatable :: rank_nonzero_array, global_nonzero        ! just for the order: from 1 to ...

contains

    subroutine petsc_initialize

!#include <petsc/finclude/petscksp.h>

       use petscksp

       use comai                                                   ! for neq, number of freedom per node
       use combi   

       implicit none


       call PetscInitialize(PETSC_NULL_CHARACTER,ierr_3)
       
       ! set the dof present in he a matrix
       if(idof_co2 .ne. 0)then
           idf = idof_co2
       else
           idf = idof
       endif
       
       if (ierr_3 .ne. 0) then
           print*,'Unable to initialize PETSc'
           stop
       endif
       
       nnz_total = 0
       A_size = neq * idf


       call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-A_size',A_size,flg,ierr_3)
       call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-neq',neq,flg,ierr_3)

       call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr_3)
       call MPI_Comm_size(PETSC_COMM_WORLD,rank_size,ierr_3)


       allocate(nnz(neq))
       allocate(nnz_global(neq*idf))
       allocate(nnz_above(neq*idf))

       allocate(bp_index(A_size))
       allocate(x_index(A_size))
       allocate(recvcounts(A_size))
       allocate(displs(A_size))
       allocate(solution(A_size))
       allocate(sol_gather(A_size))

       nnz          = 0
       nnz_global   = 0
       nnz_above    = 0
       bp_index     = 0
       x_index      = 0
       recvcounts   = 0
       displs       = 0
       solution     = 0.0
       sol_gather   = 0.0
       nnz_max      = 0
       Istart       = 0
       Iend         = 0

       d_nz         = 7*idf
       o_nz         = 7*idf

       do j = 1,A_size
           bp_index(j) = j-1      ! Assign index values from 0 to n-1
       end do
       
       
       


       ! obtain coefficient for each block, not the entire A matrix
      ! do j = 0, (idf-1)
		   do i = 1,neq
			   nnz(i) = nelm(i+1)-nelm(i)          ! decide how many non-zero values on each row in each block
			   nnz_total = nnz_total + nnz(i)         ! decide how many non-zero number in each block
					   
			   if (i == 1) then
				   nnz_above(i) = 0                   ! deccide how many non-zero number above that row
			   else
				   nnz_above(i) = nnz_above(i-1) + nnz(i-1)
			   endif

			 
			   !nnz_global(i+neq*j) = nnz(i)
			end do
		!end do
!       do i = 1,neq
!           nnz(i)    = nelm(i+1)-nelm(i)          ! decide how many non-zero values on each row in each block
!           nnz_total = nnz_total + nnz(i)         ! decide how many non-zero number in each block
!
!           if (i == 1) then
!               nnz_above(i) = 0                   ! deccide how many non-zero number above that row
!           else
!               nnz_above(i) = nnz_above(i-1) + nnz(i-1)
!           end if
!
!       end do
       nnz_global(1:neq) = nnz(1:neq)
       do i = 1, neq
       		do j = 0, (idf-1)
       			nnz_global(i+neq*j) = nnz(i)
       		end do
       end do
       do i = neq+1, neq*idf
       		nnz_above(i) = nnz_above(i-1) + nnz_global(i-1)
       end do


       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !         Compute the matrix A for the linear system, Ax = b
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !
       !  Create matrix.  When using MatCreate(), the matrix format can
       !  be specified at runtime.
       	   
       ! More robust algorithm for calculation of nz for unstructured grid !
       if (rank_size == 1) then

             d_nz = 7 * idf          ! 3 dimension * number of variables
             o_nz = 0                 ! o_nz: off diagonal    d_nz: diagonal non-zero

       else if (rank_size == 2) then

             d_nz = 7 * ceiling(idf/1.00/rank_size)
             o_nz = 7 * ceiling(idf/1.00/rank_size)

       else if (rank_size == 3) then

             d_nz = 7 * ceiling(idf/1.00/rank_size)
             o_nz = 7*idf - floor(idf/1.00/rank_size)

       else if (rank_size == 4) then

             d_nz = ceiling(7 * idf/1.00/rank_size)
             o_nz = 7*idf - floor(7 * idf/1.00/rank_size)

       else

             d_nz = ceiling(7 * idf/1.00/rank_size)
             o_nz = 7*idf - floor(7 * idf/1.00/rank_size)

       end if

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! - - - -
       !         Compute the matrix A for the linear system, Ax = b
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !
       !  Create matrix.  When using MatCreate(), the matrix format can
       !  be specified at runtime.
       !

       call MatCreate(PETSC_COMM_WORLD,A_matrix,ierr_3)        ! for entire large matrix
       call MatSetSizes(A_matrix,PETSC_DECIDE,PETSC_DECIDE,A_size,A_size,ierr_3)
       call MatSetFromOptions(A_matrix,ierr_3)
       call MatSetType(A_matrix,MATMPIAIJ,ierr_3)
       call MatSetUp(A_matrix,ierr_3)
       call MatMPIAIJSetPreallocation(A_matrix,d_nz,PETSC_NULL_INTEGER,o_nz,PETSC_NULL_INTEGER,ierr_3)


       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !         Compute the matrix A for the linear system, Ax = b
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !
       !  Create matrix.  When using MatCreate(), the matrix format can
       !  be specified at runtime.
       !

!       call MatCreate(PETSC_COMM_WORLD,A_matrix,ierr_3)        ! for entire large matrix
!       call MatSetSizes(A_matrix,PETSC_DECIDE,PETSC_DECIDE,A_size,A_size,ierr_3)
!       call MatSetFromOptions(A_matrix,ierr_3)
!       call MatSetUp(A_matrix,ierr_3)

!       call MatCreateSeqAIJ(PETSC_COMM_SELF,A_size,A_size,0,nnz_global,A_matrix,ierr_3)
!       call MatCreateAIJ(PETSC_COMM_WORLD,idf,idf,A_size,A_size,0,nnz,0,nnz,A_matrix,ierr_3)
!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global)

!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global,ierr_3)

!       call MatSetFromOptions(A_matrix,ierr_3)
!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global,ierr_3)

! -----------------------------------------------------------------------------
!      Decide from which row to which row for each processor
      
       call MatGetOwnershipRange(A_matrix,Istart,Iend,ierr_3)
        
       nnz_size = nnz_total*idf**2

       call VecCreate(PETSC_COMM_WORLD,nnz_seq,ierr_3)
       call VecSetSizes(nnz_seq,PETSC_DECIDE,nnz_size,ierr_3)
       call VecSetFromOptions(nnz_seq,ierr_3)
       call VecGetOwnershipRange(nnz_seq, nnz_start, nnz_end, ierr_3)

       rank_nonzero_num = nnz_end - nnz_start + 0.3*(nnz_end - nnz_start)
	   allocate(col_id(rank_nonzero_num))
	   allocate(value(rank_nonzero_num))
	   allocate(value_index(rank_nonzero_num))
	   allocate(rank_nonzero_array(rank_nonzero_num))

      ! do i =1, neq
      ! 		row_nelm = nelm(neq+i)
      ! end do

       allocate(row_index(Iend - Istart+1))

       col_id      = 0
       value       = 0.0
       row_index   = 0
       value_index = 0
       nnz_sum = 0
       neq_inc = neq
       neq_adj = 0
       idf_cnt = 0

       rank_nonzero_array = (/ (I, I = 1,rank_nonzero_num) /)


       call VecCreate(PETSC_COMM_WORLD,x,ierr_3)
       call VecSetSizes(x,PETSC_DECIDE,A_size,ierr_3)
       call VecSetFromOptions(x,ierr_3)
       call VecDuplicate(x,b,ierr_3)
             

    end subroutine petsc_initialize


end module petsc_initialize_package


