
!***********************************************************************
!BOP
! !IROUTINE: broadcast_array_dbl_1d
! !INTERFACE:

subroutine broadcast_array_dbl_1d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a vector dbl variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   use MPI


! !INPUT PARAMETERS:

   integer, intent(in) :: &
     root_pe           ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real *8, dimension(1000), intent(inout) :: &
     array             ! array to be broadcast

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   

!-----------------------------------------------------------------------
!EOC

end subroutine broadcast_array_dbl_1d

!***********************************************************************






! contains mpi vars
module rank_stor
	implicit none
	integer :: ierr_1
	real *8 petsca(1000)
	real *8 petscb(1000)	
end module rank_stor


program main	
	use MPI
!	use petsc_initialize_package
            
	use petscksp
	use rank_stor
	use comai
	use petsc_initialize_package
	implicit none
	! ----------------------------------------- for PETSC_solver 

	call MPI_INIT(ierr_1)
	
	call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr_1)
	call MPI_Comm_size(MPI_COMM_WORLD,rank_size,ierr_1)

	!call petsc_initialize
    
	print *, "calling para begin with rank ", rank, rank_size
    call parallel_begin(rank, .false.)
	
! -----------------------------------------------------------

      !call petsc_finalize
      call MPI_Finalize(ierr_1)          ! for MPI running
! -----------------------------------------------------------	
end program main

!This function is where all nonzero ranks will wait for rank 0 to call with solve = true 
subroutine parallel_begin(rank_ck, solve_ck)
    use MPI
    use rank_stor
    use comai
    use petsc_initialize_package
	implicit none
	integer rank_ck
	logical solve_ck
	
	if((rank_ck .eq. 0) .and. (solve_ck .eqv. .false.))then		
		call fehm_begin()
	end if
	if( (solve_ck .eqv. .true.) .or. (rank_ck .gt. 0) ) then
		
		print *, "waiting for other ranks: ", rank_ck 
		call MPI_Barrier(MPI_COMM_WORLD, ierr_1)
		call broadcast_array_dbl_1d(petsca, rank_ck)
		print *, "solving at rank: ", rank_ck, "with a, b: ", petsca(1), petscb(1)
	end if
	
end subroutine parallel_begin
