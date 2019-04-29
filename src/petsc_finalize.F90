!******************************************************************************
!          Finalize PETSc Solver, mainly for memory
!******************************************************************************
!
module petsc_finalize_package

#include <petsc/finclude/petscksp.h>
      public :: petsc_finalize
    contains

    subroutine petsc_finalize

         use petscksp
         use petsc_initialize_package
         use petsc_package

         implicit none 

         deallocate(nnz)
         deallocate(nnz_global)
         deallocate(nnz_above)
         deallocate(bp_index)
         deallocate(x_index)
         deallocate(solution)
         deallocate(sol_gather)
         deallocate(recvcounts)
         deallocate(displs)
         deallocate(col_id)
         deallocate(value)
         deallocate(value_index)
         deallocate(row_index)
         deallocate(rank_nonzero_array)

         call VecDestroy(x,ierr_3)
         call VecDestroy(b,ierr_3)
         call MatDestroy(A_matrix,ierr_3)
         call KSPDestroy(ksp,ierr_3)

         call PetscFinalize(ierr_3) 

    end subroutine petsc_finalize


end module petsc_finalize_package
