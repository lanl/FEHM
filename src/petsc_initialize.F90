!******************************************************************************
!              Initialize PETSc Solver, mainly for memory
!******************************************************************************
!
module petsc_initialize_package

#include <petsc/finclude/petscksp.h>

    use petscksp
    use comco2

    public :: petsc_initialize

    Vec    x,b, nnz_seq, global_asize
    Mat    A_matrix
    KSP    ksp
    PC     pc

    PetscInt,dimension(:),allocatable :: bp_index,x_index        ! index for right hand vector
    PetscReal,dimension(:),allocatable :: solution,sol_gather    ! index for right hand vector
    PetscReal,dimension(:),allocatable :: value                  ! value for Assemble A matrix
    PetscInt,dimension(:),allocatable :: recvcounts,displs       ! index for MPI_Gatherv

    PetscBool        flg
    integer(kind=4)   ierr_3

    PetscInt         A_size, nnz_start, nnz_end, nnz_size, local_asize, localstart, localend
    PetscInt,dimension(:),allocatable :: nnz,nnz_above            ! nnz: number of nonzero on each row
    PetscInt,dimension(:),allocatable :: nnz_global               ! nnz: number of nonzero on each row

    integer :: nnz_total,nnz_max, nnz_sum, row_nnz, i, j
    integer :: rank_nonzero_num, rank_rows, nnz_block_above
    integer :: idf, neq_inc, idf_cnt, neq_adj, idf_adj, neq_swc, nnz_cnt
    integer(kind=4) :: d_nz, o_nz, Istart, Iend
    integer(kind=4),dimension(:),allocatable :: col_id, d_nnz, o_nnz, row_tmp                    ! for assemble A_matrix
    REAL(8),dimension(:,:),allocatable :: test_id

    integer(kind=4),dimension(:),allocatable :: row_index, row_nelm       ! for assemble A_matrix
    integer,dimension(:),allocatable :: value_index               ! the index for value vector to assemble
    integer,dimension(:),allocatable :: rank_nonzero_array, global_nonzero        ! just for the order: from 1 to ...


    double precision info(MAT_INFO_SIZE)
    double precision mal, nz_a

contains

    subroutine petsc_initialize

!#include <petsc/finclude/petscksp.h>

       use petscksp

       use comai                                                   ! for neq, number of freedom per node
       use combi
       use comdti, only : nr
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
!       if (rank_size == 1) then
!
!             !d_nz = 9 * idf         ! 3 dimension * number of variables
!             !o_nz = 0                 ! o_nz: off diagonal    d_nz: diagonal non-zero
!             d_nz = nr * idf
!             o_nz = 0
!
!!       else if (rank_size == 2) then
!!
!!             d_nz = 7 * ceiling(idf/1.00/rank_size)
!!             o_nz = 7 * ceiling(idf/1.00/rank_size)
!!
!!       else if (rank_size == 3) then
!!
!!             d_nz = 7 * ceiling(idf/1.00/rank_size)
!!             o_nz = 7*idf - floor(idf/1.00/rank_size)
!!
!!       else if (rank_size == 4) then
!!
!!             d_nz = ceiling(7 * idf/1.00/rank_size)
!!             o_nz = 7*idf - floor(7 * idf/1.00/rank_size)
!!
!       else
!
!             !d_nz = ceiling(9 * idf/1.00/rank_size)
!             !o_nz = 9*idf - floor(9 * idf/1.00/rank_size)
!             d_nz = nr * idf
!             o_nz = nr * idf
!
!       end if
       d_nz = 0
       o_nz = 0


       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! - - - -
       !         Compute the matrix A for the linear system, Ax = b
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !
       !  Create matrix.  When using MatCreate(), the matrix format can
       !  be specified at runtime.
       !
       !call VecCreate(PETSC_COMM_WORLD,global_asize,ierr_3)
       !call VecSetSizes(global_asize,PETSC_DECIDE,A_size,ierr_3)
       !call VecSetFromOptions(global_asize,ierr_3)
       !call VecGetOwnershipRange(global_asize, localstart, localend, ierr_3)
       !local_asize = localend - localstart


       call MatCreate(PETSC_COMM_WORLD,A_matrix,ierr_3);CHKERRA(ierr_3)        ! for entire large matrix
       call MatSetSizes(A_matrix,PETSC_DECIDE,PETSC_DECIDE,A_size,A_size,ierr_3);CHKERRA(ierr_3)
       call MatSetFromOptions(A_matrix,ierr_3);CHKERRA(ierr_3)

       call MatSetType(A_matrix,MATMPIAIJ,ierr_3);CHKERRA(ierr_3)

       call MatSetUp(A_matrix,ierr_3);CHKERRA(ierr_3)
       !call MatSeqAIJSetPreallocation(A_matrix, d_nz, PETSC_NULL_INTEGER, ierr_3);CHKERRA(ierr_3)



       call MatGetInfo(A_matrix, MAT_LOCAL, info, ierr_3);CHKERRA(ierr_3)
       mal = info(MAT_INFO_MALLOCS)
       nz_a = info(MAT_INFO_NZ_ALLOCATED)

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

       ! size of diagonal nnz array for preallocation
       nnz_size = Iend - Istart

       call VecCreate(PETSC_COMM_WORLD,nnz_seq,ierr_3)
       call VecSetSizes(nnz_seq,PETSC_DECIDE,nnz_total*idf,ierr_3)
       call VecSetFromOptions(nnz_seq,ierr_3)
       call VecGetOwnershipRange(nnz_seq, nnz_start, nnz_end, ierr_3)

       rank_nonzero_num = nnz_total*idf**2
	   allocate(col_id(rank_nonzero_num))
	   allocate(value(rank_nonzero_num))
	   allocate(value_index(rank_nonzero_num))
	   allocate(rank_nonzero_array(rank_nonzero_num))

       !allocate(row_index(Iend - Istart+1))

       ! alloc for using MatSetValues
       allocate(row_index(rank_nonzero_num))
       allocate(row_tmp(rank_nonzero_num))
       allocate(d_nnz(Iend-Istart))
       allocate(o_nnz(Iend-Istart))

       row_index   = 0
       row_nnz = 1
       idf_cnt = 0
       neq_inc = 0
       neq_inc = neq
       neq_adj = 0
       idf_cnt = 0
       idf_adj = 0
       rank_nonzero_array = (/ (I, I = 1,rank_nonzero_num) /)
       value_index = 0
       nnz_block_above = 0
       d_nnz = 0
       o_nnz = 0
       nnz_cnt = 0
       row_tmp = 0

       do i  = 1, neq*idf
            neq_swc = 0

            row_index(row_nnz:(nnz_global(i)*idf+row_nnz-1)) = i - 1


            if(i > neq_inc)then
                    idf_adj = idf_adj + (idf-1)
                    idf_cnt = idf_cnt + 1
                    neq_inc = neq_inc + neq
                    neq_adj = neq_adj + neq
            end if

            row_tmp(row_nnz:(nnz_global(i)*idf+row_nnz-1)) = i - neq_adj

            do j = 0, idf - 1
                            col_id(row_nnz+nnz_global(i)*j:(nnz_global(i)*(j+1)+row_nnz-1)) = nelm(nelm(i-neq_adj)+1:nelm(i-neq_adj+1)) - 1	+ neq_swc
                            value_index(row_nnz+nnz_global(i)*j:(nnz_global(i)*(j+1)+row_nnz-1)) = rank_nonzero_array(1:nnz_global(i)) + nnz_above(i) + nnz_total*j + nnz_total*idf_adj
            neq_swc =  neq_swc + neq
            end do
            row_nnz = row_nnz + nnz_global(i) * idf
            nnz_block_above = nnz_block_above + nnz_global(i)
       end do


       do i = 1, rank_nonzero_num

       		if (row_index(i) >= Istart .AND. row_index(i) <= Iend - 1)then

       			if (col_id(i) >= Istart .AND. col_id(i) <= Iend - 1) then
					d_nnz(row_index(i) - Istart + 1) = d_nnz(row_index(i) - Istart + 1) + 1
				else
					o_nnz(row_index(i) - Istart + 1) = o_nnz(row_index(i) - Istart + 1) + 1
				end if

       		end if

       end do

 !      do while (i <= nnz_size - 1)
 !      		if (col_id(nnz_cnt) >= Istart .AND. col_id(nnz_cnt) <= Iend - 1) then
 !      			d_nnz(i+1) = d_nnz(i+1) + 1
 !      		else
 !      			o_nnz(i+1) = o_nnz(i+1) + 1
 !      		end if
 !      		nnz_cnt = nnz_cnt + 1
 !      		i = row_index(nnz_cnt)
 !      end do
!       do i = 1, Iend-Istart
!       		row_index(i) = Istart + i
!       end do

       !col_id      = 0
       value       = 0.0

       nnz_sum = 0
       neq_inc = neq
       neq_adj = 0
       idf_cnt = 0

       call MatMPIAIJSetPreallocation(A_matrix, d_nz, d_nnz, o_nz, o_nnz, ierr_3);CHKERRA(ierr_3)

       call VecCreate(PETSC_COMM_WORLD,x,ierr_3)
       call VecSetSizes(x,PETSC_DECIDE,A_size,ierr_3)
       call VecSetFromOptions(x,ierr_3)
       call VecDuplicate(x,b,ierr_3)

       call KSPCreate(PETSC_COMM_WORLD,ksp,ierr_3)

    end subroutine petsc_initialize


end module petsc_initialize_package