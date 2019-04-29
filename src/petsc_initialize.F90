!******************************************************************************
!              Initialize PETSc Solver, mainly for memory
!******************************************************************************
!
module petsc_initialize_package
    
#include <petsc/finclude/petscksp.h>            

    use petscksp

    public :: petsc_initialize
    
    Vec    x,b
    Mat    A_matrix
    KSP    ksp
    PC     pc

    PetscInt,dimension(:),allocatable :: bp_index,x_index        ! index for right hand vector
    PetscReal,dimension(:),allocatable :: solution,sol_gather    ! index for right hand vector
    PetscReal,dimension(:),allocatable :: value                  ! value for Assemble A matrix
    PetscInt,dimension(:),allocatable :: recvcounts,displs       ! index for MPI_Gatherv

    PetscBool        flg
    PetscErrorCode   ierr_3

    PetscInt         A_size,i,j
    PetscInt,dimension(:),allocatable :: nnz,nnz_above            ! nnz: number of nonzero on each row
    PetscInt,dimension(:),allocatable :: nnz_global               ! nnz: number of nonzero on each row

    integer :: nnz_total,nnz_max 
    integer :: Istart,Iend,rank_nonzero_num
    integer :: d_nz,o_nz

    integer,dimension(:),allocatable :: col_id                    ! for assemble A_matrix
    integer,dimension(:),allocatable :: row_index                 ! for assemble A_matrix
    integer,dimension(:),allocatable :: value_index               ! the index for value vector to assemble
    integer,dimension(:),allocatable :: rank_nonzero_array        ! just for the order: from 1 to ...

contains

    subroutine petsc_initialize

! #include <petsc/finclude/petscksp.h>
       use petscksp

       use comai                                                   ! for neq, number of freedom per node
       use combi   

       implicit none


       call PetscInitialize(PETSC_NULL_CHARACTER,ierr_3)

       if (ierr_3 .ne. 0) then
           print*,'Unable to initialize PETSc'
           stop
       endif

       nnz_total = 0
       A_size = neq * idof


       call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-A_size',A_size,flg,ierr_3)
       call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-neq',neq,flg,ierr_3)

       call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr_3)
       call MPI_Comm_size(PETSC_COMM_WORLD,rank_size,ierr_3)


       allocate(nnz(neq))
       allocate(nnz_global(idof*neq))
       allocate(nnz_above(neq))

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

       d_nz         = 7*idof
       o_nz         = 7*idof

       do j = 1,A_size
           bp_index(j) = j-1      ! Assign index values from 0 to n-1
       end do


       ! obtain coefficient for each block, not the entire A matrix
       do i = 1,neq
           nnz(i)    = nelm(i+1)-nelm(i)          ! decide how many non-zero values on each row in each block
           nnz_total = nnz_total + nnz(i)         ! decide how many non-zero number in each block

           if (i == 1) then
               nnz_above(i) = 0                   ! deccide how many non-zero number above that row
           else
               nnz_above(i) = nnz_above(i-1) + nnz(i-1)
           end if

        end do

        nnz_global(1:neq) = nnz(1:neq)


       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !         Compute the matrix A for the linear system, Ax = b
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !
       !  Create matrix.  When using MatCreate(), the matrix format can
       !  be specified at runtime.

       if (rank_size == 1) then

             d_nz = 7 * idof          ! 3 dimension * number of variables
             o_nz = 0                 ! o_nz: off diagonal    d_nz: diagonal non-zero

       else if (rank_size == 2) then

             d_nz = 7 * ceiling(idof/1.00/rank_size)
             o_nz = 7 * ceiling(idof/1.00/rank_size)

       else if (rank_size == 3) then

             d_nz = 7 * ceiling(idof/1.00/rank_size)
             o_nz = 7*idof - floor(idof/1.00/rank_size)

       else if (rank_size == 4) then

             d_nz = ceiling(7 * idof/1.00/rank_size)
             o_nz = 7*idof - floor(7 * idof/1.00/rank_size)

       else

             d_nz = ceiling(7 * idof/1.00/rank_size)
             o_nz = 7*idof - floor(7 * idof/1.00/rank_size)

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
!       call MatCreateAIJ(PETSC_COMM_WORLD,idof,idof,A_size,A_size,0,nnz,0,nnz,A_matrix,ierr_3)
!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global)

!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global,ierr_3)

!       call MatSetFromOptions(A_matrix,ierr_3)
!       call MatSeqAIJSetPreallocation(A_matrix,0,nnz_global,ierr_3)

! -----------------------------------------------------------------------------
!      Decide from which row to which row for each processor
      
       call MatGetOwnershipRange(A_matrix,Istart,Iend,ierr_3)

       if(Iend == neq) then

           rank_nonzero_num = (nnz_total-nnz_above(Istart+1))*2
           allocate(col_id(rank_nonzero_num))
           allocate(value(rank_nonzero_num))
           allocate(value_index(rank_nonzero_num))
           allocate(rank_nonzero_array(rank_nonzero_num))

       else if (Iend < neq) then    
    
           rank_nonzero_num = (nnz_above(Iend+1)-nnz_above(Istart+1))*2   
           allocate(col_id(rank_nonzero_num))
           allocate(value(rank_nonzero_num))
           allocate(value_index(rank_nonzero_num))
           allocate(rank_nonzero_array(rank_nonzero_num))

       else if (Iend == neq*idof) then

           if (Istart < neq) then
           
               rank_nonzero_num = (nnz_total-nnz_above(Istart+1)+nnz_total)*2
               allocate(col_id(rank_nonzero_num))
               allocate(value(rank_nonzero_num))
               allocate(value_index(rank_nonzero_num))
               allocate(rank_nonzero_array(rank_nonzero_num)) 

           else if (Istart >= neq) then

               rank_nonzero_num = (nnz_total-nnz_above(Istart+1-neq))*2
               allocate(col_id(rank_nonzero_num))
               allocate(value(rank_nonzero_num))
               allocate(value_index(rank_nonzero_num))
               allocate(rank_nonzero_array(rank_nonzero_num))

           end if 

       else if ((Iend < neq*idof) .AND. (Iend > neq)) then

           if (Istart < neq) then

               rank_nonzero_num = (nnz_total-nnz_above(Istart)+nnz_above(Iend+1-neq))*2 
               allocate(col_id(rank_nonzero_num))
               allocate(value(rank_nonzero_num))
               allocate(value_index(rank_nonzero_num))
               allocate(rank_nonzero_array(rank_nonzero_num))

           else if (Istart >= neq) then

               rank_nonzero_num = (nnz_above(Iend+1-neq)-nnz_above(Istart+1-neq))*2
               allocate(col_id(rank_nonzero_num))
               allocate(value(rank_nonzero_num))
               allocate(value_index(rank_nonzero_num))
               allocate(rank_nonzero_array(rank_nonzero_num))

           end if 

       end if 


       allocate(row_index(Iend - Istart+1))

       col_id      = 0
       value       = 0.0
       row_index   = 0
       value_index = 0

       rank_nonzero_array = (/ (I, I = 1,rank_nonzero_num) /)


       call VecCreate(PETSC_COMM_WORLD,x,ierr_3)
       call VecSetSizes(x,PETSC_DECIDE,A_size,ierr_3)
       call VecSetFromOptions(x,ierr_3)
       call VecDuplicate(x,b,ierr_3)


    end subroutine petsc_initialize


end module petsc_initialize_package


