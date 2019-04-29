!******************************************************************************
!                        PETSc solver for parallel computing
!                               Linear solver with KSP
!                                                 - by Jack Cheng An 06/27/2018
!******************************************************************************
!
module petsc_package

  public :: petsc_solver
contains


subroutine petsc_solver(a,bp,nmat,nrhs,nelm,tollr)
!
!,a_num,bp_num,nmat_num,nrhs_num,nelm_num)
!  neq: number of node
!   a : value of Jacobian vector
!  bp : residual vector
! nmat: number for each block in Jacobian, such as dm/dp, dm/dt, de/dp, de/dt
! nrhs: number for each block in residual
! nelm: pointer and index of each non-zero value in Jacobian
! tollr: defult 1E-5
! a_num,bp_num,nmat_num,nrhs_num,nelm_num: size for these input vector
!
#include <petsc/finclude/petscksp.h>

      use petscksp

      use petsc_initialize_package
      use comai

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
!        This problem is for two variables: P + T (so far)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     ksp     - linear solver context
!     ksp      - Krylov subspace method context
!     pc       - preconditioner context
!     x, b, u  - approx solution, right-hand-side, exact solution vectors
!     A_matrix - matrix that defines linear system
!     its      - iterations for convergence
!     norm     - norm of error in solution

      real(kind=8),intent(in) :: tollr
      real(kind=8),dimension(:),intent(in) :: a
      real(kind=8),dimension(:),intent(inout) :: bp
      integer,dimension(:),intent(in) :: nmat
      integer,dimension(:),intent(in) :: nrhs
      integer,dimension(:),intent(in) :: nelm

      integer :: b_row,b_col,k
      integer :: rank_rows,rank_cols
      integer :: nnz_block_above,value_block_above
      integer :: a_start_num,a_end_num
      
      PetscReal        norm,tol
      PetscInt         its,nnz_count          ! n: dimension of A matrix
      PetscInt         Istart_v,Iend_v        ! start and end row number for each rank
      PetscInt         bcore_size             ! vector size for each core/processor

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr_3)
      call MPI_Comm_size(MPI_COMM_WORLD,rank_size,ierr_3)

      col               = 0
      rank_rows         = 0
      rank_cols         = 0
      value_block_above = 0

! Create sequential AIJ (CSR) Sparse Matrices
! call MatCreateSeqAIJ(PETSC_COMM_SELF,neq,neq,nz,nnz,A_sub,ierr)
!
!                                                         dm_dp     dm_dT
! - - dm_dp:    insert value by looping over nodes  ----                   -------
!                                                         de_dp     de_dT
!
! -----------------------------------------------------------------------------

      rank_rows       = Iend - Istart
      row_index(1)    = 0
      col_id          = 0
      nnz_block_above = 0

      if(Iend .LE. neq) then 

            do i = Istart+1,Iend                  ! Matrix start [0][0], row[0] to row[n-1]; (i=1) is row[0]

                  row_index(i-Istart+1) = row_index(i-Istart) + nnz_global(i)*idof

                  do k = 1,idof

                      col_id(nnz_block_above*2+1+nnz_global(i)*(k-1):nnz_block_above*2+nnz_global(i)*k) =   &
                      &      nelm(nelm(i)+1:nelm(i+1)) -1 + (k-1)*neq

                      value_index(nnz_block_above*2+(k-1)*nnz_global(i)+1:nnz_block_above*2+(k-1)*nnz_global(i)+nnz_global(i)) &
                           &  =  rank_nonzero_array(1:nnz_global(i)) + nnz_above(Istart+1) + nnz_block_above +(k-1)*nnz_total

                  end do

                  nnz_block_above = nnz_block_above + nnz_global(i)

            end do

      else if ( (Iend > neq) .AND. (Istart .GE. neq) ) then

            do i = Istart+1,Iend                  ! Matrix start [0][0], row[0] to row[n-1]; (i=1) is row[0]

                  nnz_global(i) = nnz(i-neq)
                 
                  row_index(i-Istart+1) = row_index(i-Istart) + nnz_global(i)*idof

                  do k = 1,idof

                     col_id(nnz_block_above*2+1+nnz_global(i)*(k-1):nnz_block_above*2+nnz_global(i)*k) =   &
                     &      nelm(nelm(i-neq)+1:nelm(i+1-neq)) -1 + (k-1)*neq

                      value_index(nnz_block_above*2+(k-1)*nnz_global(i)+1:nnz_block_above*2+k*nnz_global(i))     &
                      &  = rank_nonzero_array(1:nnz_global(i)) +nnz_total*2 +(k-1)*nnz_total + nnz_above(i-neq)

                  end do

                  nnz_block_above = nnz_block_above + nnz_global(i)

            end do

      else if ( (Iend > neq) .AND. (Istart<neq) ) then

            do i = Istart+1,Iend                  ! Matrix start [0][0], row[0] to row[n-1]; (i=1) is row[0]

                  if (i > neq) then
                     nnz_global(i) = nnz(i-neq)
                  end if 

                  row_index(i-Istart+1) = row_index(i-Istart) + nnz_global(i)*idof

                  do k = 1,idof

                      if (i >  neq) then
                           col_id(nnz_block_above*2+1+nnz_global(i)*(k-1):nnz_block_above*2+nnz_global(i)*k)  =   &
                           &      nelm(nelm(i-neq)+1:nelm(i+1-neq)) -1 +(k-1)*neq

                      value_index(nnz_block_above*2+(k-1)*nnz_global(i)+1:nnz_block_above*2+k*nnz_global(i)) &
                      &  = rank_nonzero_array(1:nnz_global(i)) +nnz_total*2 +(k-1)*nnz_total + nnz_above(i-neq)

                      else 
                           col_id(nnz_block_above*2+1+nnz_global(i)*(k-1):nnz_block_above*2+nnz_global(i)*k) =   &
                           &      nelm(nelm(i)+1:nelm(i+1)) -1 + (k-1)*neq

                      value_index(nnz_block_above*2+(k-1)*nnz_global(i)+1:nnz_block_above*2+(k-1)*nnz_global(i)+nnz_global(i))&
                           &  =  rank_nonzero_array(1:nnz_global(i))+ nnz_above(Istart+1)+ nnz_block_above +(k-1)*nnz_total

                      end if 

                  end do

                  nnz_block_above = nnz_block_above + nnz_global(i)

            end do

       end if 

       value(1:nnz_block_above*idof) = a(value_index(1:nnz_block_above*idof))


       call MatCreateMPIAIJWithArrays(MPI_COMM_WORLD,rank_rows,PETSC_DECIDE,A_size,A_size, &
            &    row_index,col_id,value,A_matrix,ierr_3)

       call MatAssemblyBegin(A_matrix,MAT_FINAL_ASSEMBLY,ierr_3)
       call MatAssemblyEnd(A_matrix,MAT_FINAL_ASSEMBLY,ierr_3)

!       call MatView(A_matrix,PETSC_VIEWER_STDOUT_WORLD,ierr_3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Create the vector b for the linear system, Ax = b
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Set exact solution then assemble right-hand-side vector from bp and nrhs


      call VecGetOwnershipRange(b,Istart_v,Iend_v,ierr_3)

      bcore_size = Iend_v - Istart_v
         
      call VecSetValues(b,bcore_size,bp_index(Istart_v+1:Iend_v),bp(Istart_v+1:Iend_v),INSERT_VALUES,ierr_3) 
      ! insert values from bp to b, each core for one part


      call VecAssemblyBegin(b,ierr_3)

      call VecAssemblyEnd(b,ierr_3)

!      call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr_3)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Create linear solver context

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr_3)

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.

      call KSPSetOperators(ksp,A_matrix,A_matrix,ierr_3)

!  Set linear solver defaults for this problem (optional).
!   - By extracting the KSP and PC contexts from the KSP context,
!     we can then directly call any KSP and PC routines to set various options.
!   - The following four statements are optional; all of these
!     parameters could alternatively be specified at runtime via
!     KSPSetFromOptions();

      call KSPGetPC(ksp,pc,ierr_3)
      call KSPSetType(ksp,KSPGMRES,ierr_3)
!      call KSPGMRESSetRestart(ksp,20,ierr)

      call PCSetType(pc,PCASM,ierr_3)      ! set up preconditioners

      tol = tollr                           ! convergence tollrance

      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,       &
      &    PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr_3)

!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.

!      call KSPSetFromOptions(ksp,ierr)
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call KSPSolve(ksp,b,x,ierr_3)
!      call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr_3)

      x_index(1:bcore_size) = bp_index(Istart_v+1:Iend_v)     ! decide index of x for each processor

      call VecGetValues(x,bcore_size,x_index,solution,ierr_3)     ! Pass the values of x to solution
      
      ! gather how many unknowns for each processor / rank
      call MPI_GATHER(bcore_size,1,MPI_INT,recvcounts,1,MPI_INT,0,MPI_COMM_WORLD,ierr_3)

      ! displs: the index of start point to restore solution
      do i = 1,rank_size
          if (i == 1) then
              displs(i) = 0
          end if
          
          if (i > 1) then        
              displs(i) = displs(i-1) + recvcounts(i-1)
          end if   

      end do

      ! Gather values from each processor to rank 0
      call MPI_GATHERV(solution,bcore_size,MPI_REAL8,solution,recvcounts,displs,MPI_REAL8,0,MPI_COMM_WORLD,ierr_3)

      if (rank == 0 ) then
        bp(1:A_size) = 0.0
        bp(1:A_size) = bp(1:A_size) + solution(1:A_size)
      end if

      ! send the solver solution bp to all the processors
      call MPI_Bcast(bp,A_size,MPI_REAL8,0,MPI_COMM_WORLD,ierr_3)


end subroutine petsc_solver ! End of SUBROUTINE PETSc_solver
!
!
!
end module petsc_package
!
!*************************************************************************

