#ifdef DDFLAG
module decomp_unstruct

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private 

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: UGridRead, &
            UGridReadHDF5, &
            UGridReadSurfGrid, &
            UGridReadHDF5SurfGrid, &
            UGridDecompose, &
            UGridComputeInternConnect, &
            UGridPopulateConnection, &
            UGridComputeCoord, &
            UGridComputeVolumes, &
            UGridComputeAreas, &
            UGridComputeQuality, &
            UGridGetCellFromPoint, &
            UGridGetCellsInRectangle, &
            UGridEnsureRightHandRule, &
            UGridMapSideSet, &
            UGridMapSideSet2, &
            UGridGrowStencilSupport, &
            UGridMapBoundFacesInPolVol

contains

! ************************************************************************** !

function UGridCreate()
  ! 
  ! Creates an unstructured grid object
  ! 
  ! 

  implicit none
  
  type(grid_unstructured_type), pointer :: UGridCreate

  type(grid_unstructured_type), pointer :: unstructured_grid

  allocate(unstructured_grid)

  ! variables for all unstructured grids
  unstructured_grid%num_ghost_cells = 0
  unstructured_grid%global_offset = 0
  unstructured_grid%nmax = 0
  unstructured_grid%nlmax = 0
  unstructured_grid%ngmax = 0
  nullify(unstructured_grid%hash)
  unstructured_grid%num_hash = 100
  nullify(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  unstructured_grid%ao_natural_to_petsc = 0
  nullify(unstructured_grid%explicit_grid)
  nullify(unstructured_grid%polyhedra_grid)

  ! variables for implicit unstructured grids
  unstructured_grid%grid_type = THREE_DIM_GRID
  unstructured_grid%num_vertices_global = 0
  unstructured_grid%num_vertices_local = 0
  unstructured_grid%max_ndual_per_cell = 0
  unstructured_grid%max_nvert_per_cell = 0
  unstructured_grid%max_cells_sharing_a_vertex = 24
  nullify(unstructured_grid%cell_type)
  nullify(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%cell_to_face_ghosted)
  nullify(unstructured_grid%vertex_ids_natural)
  nullify(unstructured_grid%vertices)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
  nullify(unstructured_grid%connection_to_face)
  nullify(unstructured_grid%face_centroid)
  nullify(unstructured_grid%face_area)
  nullify(unstructured_grid%nat_ids_of_other_grid)

  UGridCreate => unstructured_grid
  
end function UGridCreate

! ************************************************************************** !
! ************************************************************************** !

subroutine UGridCreateUGDM(unstructured_grid,ugdm,ndof,option)
  ! 
  ! Constructs mappings / scatter contexts for PETSc DM
  ! object

  
#include "petsc/finclude/petscdm.h"
  use petscdm
  use Option_module
  use Utility_module, only: ReallocateArray
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(ugdm_type), pointer :: ugdm
  PetscInt :: ndof
  type(option_type) :: option
  
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  IS :: is_tmp
  Vec :: vec_tmp
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: ndof_word
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscViewer :: viewer

  PetscInt, allocatable :: int_array(:)
  
  ugdm => UGDMCreate()
  ugdm%ndof = ndof

#if UGRID_DEBUG
  write(ndof_word,*) ndof
  ndof_word = adjustl(ndof_word)
  ndof_word = '_' // trim(ndof_word)
  string = 'Vectors' // ndof_word
  call PrintMsg(option,string)
#endif

  ! create global vec
  !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
  !                  PETSC_DETERMINE,ugdm%global_vec,ierr)
  call VecCreate(option%mycomm,ugdm%global_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(ugdm%global_vec,unstructured_grid%nlmax*ndof, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(ugdm%global_vec,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(ugdm%global_vec,ierr);CHKERRQ(ierr)

  ! create local vec
  !call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax*ndof, &
  !                  ugdm%local_vec,ierr)
  call VecCreate(PETSC_COMM_SELF,ugdm%local_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(ugdm%local_vec,unstructured_grid%ngmax*ndof,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(ugdm%local_vec,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(ugdm%local_vec,ierr);CHKERRQ(ierr)
  
  ! IS for global numbering of local, non-ghosted cells
!geh  call VecGetOwnershipRange(ugdm%global_vec,istart,iend,ierr)
  ! ISCreateBlock requires block ids, not indices.  Therefore, istart should be
  ! the offset of the block from the beginning of the vector.
!geh  istart = istart / ndof
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax
 !geh   int_array(local_id) = (local_id-1)+istart
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo

  ! arguments for ISCreateBlock():
  ! option%mycomm  - the MPI communicator
  ! ndof  - number of elements in each block
  ! unstructured_grid%nlmax  - the length of the index set
  !                                      (the number of blocks
  ! int_array  - the list of integers, one for each block and count
  !              of block not indices
  ! PETSC_COPY_VALUES  - see PetscCopyMode, only PETSC_COPY_VALUES and
  !                      PETSC_OWN_POINTER are supported in this routine
  ! ugdm%is_local_petsc - the new index set
  ! ierr - PETScErrorCode
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_local_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_local_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do ghosted_id = 1, unstructured_grid%num_ghost_cells
    int_array(ghosted_id) = (ghosted_id+unstructured_grid%nlmax-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_local, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosts_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_ghosts_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
#if UGRID_DEBUG
  string = 'Index Sets' // ndof_word
  call PrintMsg(option,'Index Sets')
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do ghosted_id = 1, unstructured_grid%num_ghost_cells
    int_array(ghosted_id) = &
      (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosts_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_ghosts_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! IS for local numbering of local, non-ghosted cells
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax
    int_array(local_id) = (local_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_local, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_local_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_local_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  ! IS for ghosted numbering of local ghosted cells
  allocate(int_array(unstructured_grid%ngmax))
  do ghosted_id = 1, unstructured_grid%ngmax
    int_array(ghosted_id) = (ghosted_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_local, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosted_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_ghosted_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
             
  ! IS for petsc numbering of local ghosted cells
  allocate(int_array(unstructured_grid%ngmax))
  do local_id = 1, unstructured_grid%nlmax
!geh    int_array(local_id) = istart+(local_id-1)
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo
  do ghosted_id = 1,unstructured_grid%num_ghost_cells
    int_array(unstructured_grid%nlmax+ghosted_id) = &
      (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosted_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISView(ugdm%is_ghosted_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif    
                 
  ! create a local to global mapping
#if UGRID_DEBUG
  string = 'ISLocalToGlobalMapping' // ndof_word
  call PrintMsg(option,string)
#endif

  call ISLocalToGlobalMappingCreateIS(ugdm%is_ghosted_petsc, &
                                      ugdm%mapping_ltog,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  string = 'mapping_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call ISLocalToGlobalMappingView(ugdm%mapping_ltog,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
               
#if UGRID_DEBUG
  string = 'local to global' // ndof_word
  call PrintMsg(option,string)
#endif

  ! Create local to global scatter
  call VecScatterCreate(ugdm%local_vec,ugdm%is_local_local,ugdm%global_vec, &
                        ugdm%is_local_petsc,ugdm%scatter_ltog, &
                        ierr);CHKERRQ(ierr)
                        
#if UGRID_DEBUG
  string = 'scatter_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(ugdm%scatter_ltog,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

#if UGRID_DEBUG
  string = 'global to local' // ndof_word
  call PrintMsg(option,string)
#endif



  ! Set up global to natural scatter
  ! Create index set of local non-ghosted Petsc ordering
  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax, &
                    PETSC_DETERMINE,vec_tmp,ierr);CHKERRQ(ierr)
!geh  call VecGetOwnershipRange(vec_tmp,istart,iend,ierr)
  call VecDestroy(vec_tmp,ierr);CHKERRQ(ierr)
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax 
!geh    int_array(local_id) = (local_id-1)+istart
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo
  call ISCreateGeneral(option%mycomm,unstructured_grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp,ierr);CHKERRQ(ierr)
  deallocate(int_array)
  call AOPetscToApplicationIS(unstructured_grid%ao_natural_to_petsc, &
                              is_tmp,ierr);CHKERRQ(ierr)
  ! remap for ndof > 1  !geh: no longer need to accommodate ndof > 1, but leave
  ! alone for now.
  allocate(int_array(unstructured_grid%nlmax))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr);CHKERRQ(ierr)
  do local_id = 1, unstructured_grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
 VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
  !                  PETSC_DETERMINE,vec_tmp,ierr)
  call VecCreate(option%mycomm,vec_tmp,ierr);CHKERRQ(ierr)
  call VecSetSizes(vec_tmp,unstructured_grid%nlmax*ndof,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vec_tmp,ierr);CHKERRQ(ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_petsc,vec_tmp, &
                        ugdm%is_local_natural,ugdm%scatter_gton, &
                        ierr);CHKERRQ(ierr)
  call VecDestroy(vec_tmp,ierr);CHKERRQ(ierr)


  ! set the ao_natural_to_petsc pointer
  ugdm%ao_natural_to_petsc = unstructured_grid%ao_natural_to_petsc

end subroutine UGridCreateUGDM

! ************************************************************************** !

subroutine UGridDMCreateVector(unstructured_grid,ugdm,vec,vec_type,option)
  ! 
  ! Creates a global vector with PETSc ordering
  ! 

  use Option_module

  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  Vec :: vec
  PetscInt :: vec_type
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  select case(vec_type)
    case(GLOBAL)
      !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
      !                  ugdm%ndof, &
      !                  PETSC_DETERMINE,vec,ierr)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,unstructured_grid%nlmax*ugdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetLocalToGlobalMapping(vec,ugdm%mapping_ltog, &
                                      ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(LOCAL)
      !call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax* &
      !                  ugdm%ndof, &
      !                  vec,ierr)
      call VecCreate(PETSC_COMM_SELF,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,unstructured_grid%ngmax*ugdm%ndof, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
    case(NATURAL)
      !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
      !                  ugdm%ndof, &
      !                  PETSC_DETERMINE,vec,ierr)
      call VecCreate(option%mycomm,vec,ierr);CHKERRQ(ierr)
      call VecSetSizes(vec,unstructured_grid%nlmax*ugdm%ndof, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr);CHKERRQ(ierr)
      call VecSetFromOptions(vec,ierr);CHKERRQ(ierr)
  end select
    
end subroutine UGridDMCreateVector

! ************************************************************************** !

subroutine UGridMapIndices(unstructured_grid,ugdm,nG2L,nL2G,nG2A,option)
  ! 
  ! maps global, local and natural indices of cells to each other
  ! 


  use Option_module

  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nG2A(:)
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id

  allocate(nG2L(unstructured_grid%ngmax))
  allocate(nL2G(unstructured_grid%nlmax))
  allocate(nG2A(unstructured_grid%ngmax))
  
  ! initialize ghosted to 0
  !geh: any index beyond %nlmax will be 0 indicating that there is no local
  !     counterpart (i.e., it is a ghost cell)
  nG2L = 0

  !geh: Yes, it seems redundant that that we are setting both nL2G and nG2L to 
  !     the same index, but keep in mind that nG2L extends beyond %nlmax and
  !     we need these arrays to provide seemless integration for structured and
  !     unstructured
  do local_id = 1, unstructured_grid%nlmax
    nL2G(local_id) = local_id
    nG2L(local_id) = local_id
  enddo

  call ISGetIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    nG2A(ghosted_id) = int_ptr(ghosted_id)+1
  enddo
  call ISRestoreIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr);CHKERRQ(ierr)
  nG2A = nG2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%ngmax, &
                            nG2A,ierr);CHKERRQ(ierr)
  nG2A = nG2A + 1 ! 1-based


end subroutine UGridMapIndices

! ************************************************************************** !

subroutine UGridPartition(ugrid,option,Dual_mat,is_new, &
                          num_cells_local_new)
  ! 
  ! UGridGet_Dual_Part_IS: Given an adjacency matrix, calculates the dual
  ! partitions, and provides a new IS with the ids
  ! of the local cells on the processor
  ! 
  !

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  
  implicit none
  
  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  Mat :: Dual_mat
  IS :: is_new
  PetscInt :: num_cells_local_new

  MatPartitioning :: Part
  PetscInt, allocatable :: cell_counts(:)
  PetscInt :: iflag
  PetscInt :: tempint
  PetscViewer :: viewer
  PetscInt :: local_vertex_offset
  PetscErrorCode :: ierr

#if UGRID_DEBUG
  call PrintMsg(option,'Partitioning')
#endif

  ! create the partitioning
  call MatPartitioningCreate(option%mycomm,Part,ierr);CHKERRQ(ierr)
  ! MatPartitioningSetAdjacency sets the adjacency graph (matrix) of the 
  ! thing to be partitioned.  - petsc
  call MatPartitioningSetAdjacency(Part,Dual_mat,ierr);CHKERRQ(ierr)
  call MatPartitioningSetFromOptions(Part,ierr);CHKERRQ(ierr)
  ! MatPartitioningApply gets a partitioning for a matrix. For each local cell 
  ! this tells the processor number that that cell is assigned to. - petsc
  ! is_new holds this information
  call MatPartitioningApply(Part,is_new,ierr);CHKERRQ(ierr)
  call MatPartitioningDestroy(Part,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_new,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! calculate the number of local grid cells on each processor
  allocate(cell_counts(option%mycommsize))
  ! ISPartitioningCount takes a ISPartitioning and determines the number of  
  ! resulting elements on each (partition) process - petsc
  tempint = option%mycommsize
  call ISPartitioningCount(is_new,tempint,cell_counts,ierr);CHKERRQ(ierr)
  num_cells_local_new = cell_counts(option%myrank+1) 
  call MPI_Allreduce(num_cells_local_new,iflag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MIN,option%mycomm,ierr)
  deallocate(cell_counts)
  if (iflag < 1) then
    option%io_buffer = 'A processor core has been assigned zero cells.'
    call PrintErrMsg(option)
  endif
  
end subroutine UGridPartition  

! ************************************************************************** !

subroutine UGridDecompose(unstructured_grid,option)
  ! 
  ! Decomposes an unstructured grid across ranks
  ! 
  ! 
  
#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  use Utility_module, only: ReallocateArray, SearchOrderedArray
  
  implicit none
  
  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  
  PetscInt :: local_id, local_id2
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  PetscInt :: count, vertex_count
  PetscInt :: vertex_offset, global_vertex_offset
  PetscInt :: stride
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscInt :: num_rows, num_cols, istart, iend, icol
  PetscBool :: success
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  PetscViewer :: viewer
  Mat :: Adj_mat
  Mat :: Dual_mat
  MatPartitioning :: Part
  Vec :: elements_natural
  Vec :: elements_local
  Vec :: elements_old
  Vec :: vertices_old
  Vec :: vertices_new
  IS :: is_new
  IS :: is_scatter
  IS :: is_gather

  VecScatter :: vec_scatter
  
  PetscInt :: vertex_ids_offset
  PetscInt :: dual_offset
  PetscInt :: natural_id_offset

  PetscInt :: max_int_count
  PetscInt :: temp_int
  PetscInt :: min_value
  PetscInt :: num_cells_local_new
  PetscInt :: num_cells_local_old  
  PetscInt :: global_offset_old
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, pointer :: int_array_pointer(:)
  
  PetscInt :: idual, dual_id
  PetscInt :: iflag
  PetscBool :: found

!  cell distribution across processors (size = num_cores + 1)
!  core i owns cells cell_distribution(i):cell_distribution(i+1), note
!  the zero-based indexing
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%nlmax,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)

  num_cells_local_old = unstructured_grid%nlmax

  ! recalculate maximum number of vertices for any given cell
  temp_int = 0
  min_value = 2 ! min value should be either 0 or 1 after global reduction
  do local_id = 1, num_cells_local_old
    vertex_count = 0
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point, cell vertex can be 0
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) exit
      if (unstructured_grid%cell_vertices(ivertex,local_id) < min_value) then
        min_value = unstructured_grid%cell_vertices(ivertex,local_id)
      endif
      vertex_count = vertex_count+1
    enddo
    if (vertex_count > temp_int) temp_int = vertex_count
  enddo
  call MPI_Allreduce(temp_int,unstructured_grid%max_nvert_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(min_value,index_format_flag, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN,option%mycomm,ierr)

  ! let's make it Fortran indexing
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point we may be zero-based
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) then
        ! change no_value (UNINITIALIZED_INTEGER) to '0'
        unstructured_grid%cell_vertices(ivertex,local_id) = 0
      else
        if (index_format_flag == 0) then
          ! let's make it Fortran indexing
          unstructured_grid%cell_vertices(ivertex,local_id) = &
            unstructured_grid%cell_vertices(ivertex,local_id) + 1
        endif
      endif
    enddo
  enddo

#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_nvert_per_cell
  option%io_buffer = 'Maximum number of vertices per cell: ' // adjustl(string)
  call PrintMsg(option)
  write(string,*) index_format_flag
  option%io_buffer = 'Vertex indexing starts at: ' // adjustl(string)
  call PrintMsg(option)
  if (index_format_flag == 0) then
    option%io_buffer = 'Changing vertex indexing to 1-based.'
    call PrintMsg(option)
  endif
#endif

  num_cells_local_old = unstructured_grid%nlmax 
  allocate(local_vertices(unstructured_grid%max_nvert_per_cell* &
                          num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  local_vertices = 0
  local_vertex_offset = 0
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      if (unstructured_grid%cell_vertices(ivertex,local_id) == 0) exit
      count = count + 1
      ! local vertices must be zero-based for MatCreateMPIAdj; thus subtract 1
      local_vertices(count) = &
        unstructured_grid%cell_vertices(ivertex,local_id) - 1
    enddo
    local_vertex_offset(local_id+1) = count 
  enddo
    
  select case (unstructured_grid%grid_type)
    case(TWO_DIM_GRID)
      num_common_vertices = 2 ! cells must share at least this number of vertices
    case(THREE_DIM_GRID)
      num_common_vertices = 3 ! cells must share at least this number of vertices
    case default
        option%io_buffer = 'Grid type not recognized '
        call PrintErrMsg(option)
    end select

  ! determine the global offset from 0 for cells on this rank
  global_offset_old = 0
  call MPI_Exscan(num_cells_local_old,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! create an adjacency matrix for calculating the duals (connnections)
#if UGRID_DEBUG
  call PrintMsg(option,'Adjacency matrix')
#endif

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat, &
                       ierr);CHKERRQ(ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'Adj_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'Adj_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call MatView(Adj_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

#if UGRID_DEBUG
  call PrintMsg(option,'Dual matrix')
#endif

#if defined(PETSC_HAVE_PARMETIS)
  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat, &
                          ierr);CHKERRQ(ierr)
#endif
  call MatDestroy(Adj_mat,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'Dual_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'Dual_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call MatView(Dual_mat,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call UGridPartition(unstructured_grid,option,Dual_mat,is_new, &
                      num_cells_local_new)
  
  if (allocated(local_vertices)) deallocate(local_vertices)
  if (allocated(local_vertex_offset)) deallocate(local_vertex_offset)
  
  ! second argument of ZERO_INTEGER means to use 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  if (.not.success .or. num_rows /= num_cells_local_old) then
    print *, option%myrank, num_rows, success, num_cells_local_old
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call PrintErrMsg(option)
  endif

  ! calculate maximum number of connections for any given cell
  unstructured_grid%max_ndual_per_cell = 0
  do local_id = 1, num_cells_local_old
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) &
      unstructured_grid%max_ndual_per_cell = num_cols
  enddo
  temp_int = unstructured_grid%max_ndual_per_cell
  call MPI_Allreduce(temp_int,unstructured_grid%max_ndual_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  
#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_ndual_per_cell
  option%io_buffer = 'Maximum number of duals per cell: ' // adjustl(string)
  call PrintMsg(option)
#endif
  
  if (unstructured_grid%max_ndual_per_cell > 0) then
    call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                            num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  
  ! in order to redistributed vertex/cell data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  vertex_ids_offset = 1 + 1 ! +1 for -777
  dual_offset = vertex_ids_offset + unstructured_grid%max_nvert_per_cell + 1 ! +1 for -888
  stride = dual_offset+ unstructured_grid%max_ndual_per_cell + 1 ! +1 for -999999
  natural_id_offset = 1

  ! Information for each cell is packed in a strided petsc vec
  ! The information is ordered within each stride as follows:
  ! -cell_N   ! global cell id (negative indicates 1-based)
  ! -777      ! separator between cell id and vertex ids for cell_N
  ! vertex1   ! in cell_N
  ! vertex2
  ! ...
  ! vertexN   
  ! -888      ! separator between vertex and dual ids
  ! dual1     ! dual ids between cell_N and others
  ! dual2
  ! ...
  ! dualN     
  ! -999999   ! separator indicating end of information for cell_N
  
  ! the purpose of -777, -888, and -999999 is to allow one to use cells of 
  ! various geometry.  Currently, the max # vertices = 8 and max # duals = 6.
  ! But this will be generalized in the future.
  
  call UGridCreateOldVec(unstructured_grid,option,elements_old, &
                                num_cells_local_old, &
                                is_new,is_scatter,stride)

  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE,num_rows, &
                      ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(elements_old,vec_ptr,ierr);CHKERRQ(ierr)
  count = 0
  vertex_count = 0
  do local_id = 1, num_cells_local_old
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(global_offset_old+local_id)
    count = count + 1
    ! add the separator
    vec_ptr(count) = -777  ! help differentiate
    ! add the vertex ids
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      count = count + 1
      vertex_count = vertex_count + 1
      ! increment for 1-based ordering
      vec_ptr(count) = unstructured_grid%cell_vertices(ivertex,local_id)
    enddo


    count = count + 1 
    ! another vertex/dual separator
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Dual matrix is larger then max_ndual_per_cell.'
      call PrintErrMsgByRank(option)
    endif
    do icol = 1, unstructured_grid%max_ndual_per_cell
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo
    count = count + 1 
    ! final separator
    vec_ptr(count) = -999999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr);CHKERRQ(ierr)
  
  if (unstructured_grid%max_ndual_per_cell > 0) then
    call MatRestoreRowIJF90(Dual_mat,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                            num_rows,ia_ptr,ja_ptr,success,ierr);CHKERRQ(ierr)
  endif
  call MatDestroy(Dual_mat,ierr);CHKERRQ(ierr)
 
  
  call UGridNaturalToPetsc(unstructured_grid,option, &
                           elements_old,elements_local, &
                           num_cells_local_new,stride,dual_offset, &
                           natural_id_offset,is_scatter)
  
  ! make a list of local vertices
  max_int_count = 2*unstructured_grid%ngmax
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
  ! note that the vertices are still in natural numbering
  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id=1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride))
      if (vertex_id < 1) exit
      vertex_count = vertex_count + 1
      if (vertex_count > max_int_count) then
        call ReallocateArray(int_array_pointer,max_int_count)
      endif
      vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride) = vertex_count
      int_array_pointer(vertex_count) = vertex_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)

  ! sort the vertex ids
  allocate(int_array(vertex_count))
  int_array(1:vertex_count) = int_array_pointer(1:vertex_count)
  allocate(int_array2(vertex_count))
  do ivertex = 1, vertex_count
    int_array2(ivertex) = ivertex 
  enddo
  deallocate(int_array_pointer)
  nullify(int_array_pointer)
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(vertex_count,int_array,int_array2, &
                                   ierr);CHKERRQ(ierr)
  int_array2 = int_array2+1

  ! remove duplicates
  allocate(int_array3(vertex_count))
  allocate(int_array4(vertex_count))
  int_array3 = 0
  int_array4 = 0
  int_array3(1) = int_array(int_array2(1))
  count = 1
  int_array4(int_array2(1)) = count
  do ivertex = 2, vertex_count
    vertex_id = int_array(int_array2(ivertex))
    if (vertex_id > int_array3(count)) then
      count = count + 1
      int_array3(count) = vertex_id
    endif
    int_array4(int_array2(ivertex)) = count
  enddo
  vertex_count = count
  deallocate(int_array)

  allocate(unstructured_grid%vertex_ids_natural(vertex_count))
  unstructured_grid%vertex_ids_natural = int_array3(1:vertex_count)

  ! now load all the vertices needed to define all the local cells
  ! on the processor
  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  ! allocate the array that will store the vertex ids for each cell.
  ! remember that max_nvert_per_cell is the max # of vertices in a cell
  ! currently hardwired to 8.
  deallocate(unstructured_grid%cell_vertices)
  allocate(unstructured_grid%cell_vertices( &
             0:unstructured_grid%max_nvert_per_cell,unstructured_grid%ngmax))
  unstructured_grid%cell_vertices = 0
  
  ! permute the local ids calculated earlier in the int_array4
  call VecGetArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! extract the original vertex id
      vertex_id = int(vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride))
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices(0,ghosted_id)+1
      unstructured_grid%cell_vertices(count,ghosted_id) = &
        int_array4(vertex_id)
      unstructured_grid%cell_vertices(0,ghosted_id) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride) = &
        int_array4(vertex_id)
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(int_array2)
  deallocate(int_array3)
  deallocate(int_array4)

#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'elements_vert_local' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'elements_vert_local' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(elements_local,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  call VecDestroy(elements_local,ierr);CHKERRQ(ierr)

  ! now we need to work on aligning the original vertex coordinates with 
  ! the current ordering or permuted/rearranged ordering.

  ! IS for gather operation - need local numbering
  allocate(int_array(vertex_count))
  ! vertex_count = # of local vertices (I believe ghosted+non-ghosted)
  do ivertex = 1, vertex_count
    int_array(ivertex) = ivertex-1
  enddo

  ! include cell ids (use block ids, not indices)
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  ! create a parallel petsc vector with a stride of 3.
  !call VecCreateMPI(option%mycomm,unstructured_grid%num_vertices_local*3, &
  !                  PETSC_DETERMINE,vertices_old,ierr)
  call VecCreate(option%mycomm,vertices_old,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_old,unstructured_grid%num_vertices_local*3, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_old,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_old,ierr);CHKERRQ(ierr)

  ! create serial petsc vector with a stride of 3
  !call VecCreateSeq(PETSC_COMM_SELF,vertex_count*3,vertices_new,ierr)
  call VecCreate(PETSC_COMM_SELF,vertices_new,ierr);CHKERRQ(ierr)
  call VecSetSizes(vertices_new,vertex_count*3,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetBlockSize(vertices_new,3,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(vertices_new,ierr);CHKERRQ(ierr)

!  call VecCreate(option%mycomm,vertices_new,ierr)
!  call VecSetSizes(vertices_new, &
!                   vertex_count*3,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_new,ierr)
  
!  call VecCreate(option%mycomm,vertices_old,ierr)
!  call VecSetSizes(vertices_old, &
!                   3*unstructured_grid%num_vertices_local,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_old,ierr)
! load up the coordinates
  call VecGetArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    vec_ptr((ivertex-1)*3+1) = unstructured_grid%vertices(ivertex)%x
    vec_ptr((ivertex-1)*3+2) = unstructured_grid%vertices(ivertex)%y
    vec_ptr((ivertex-1)*3+3) = unstructured_grid%vertices(ivertex)%z
  enddo
  call VecRestoreArrayF90(vertices_old,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)

  ! IS for scatter - provide petsc global numbering
  allocate(int_array(vertex_count))
  do ivertex = 1, vertex_count
    int_array(ivertex) = (needed_vertices_petsc(ivertex)-1)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_scatter, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

  ! resize vertex array to new size
  unstructured_grid%num_vertices_natural = unstructured_grid%num_vertices_local
  unstructured_grid%num_vertices_local = vertex_count
  allocate(unstructured_grid%vertices(vertex_count))
  do ivertex = 1, vertex_count
    unstructured_grid%vertices(ivertex)%x = 0.d0
    unstructured_grid%vertices(ivertex)%y = 0.d0
    unstructured_grid%vertices(ivertex)%z = 0.d0
  enddo

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_scatter,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_vert_old_to_new_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_vert_old_to_new_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call ISView(is_gather,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecScatterCreate(vertices_old,is_scatter,vertices_new,is_gather, &
                        vec_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_gather,ierr);CHKERRQ(ierr)
  call VecScatterBegin(vec_scatter,vertices_old,vertices_new, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,vertices_old,vertices_new, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old_subsurf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old_surf.out',viewer, &
                              ierr);CHKERRQ(ierr)
  endif
  call VecView(vertices_old,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecDestroy(vertices_old,ierr);CHKERRQ(ierr)


  call VecGetArrayF90(vertices_new,vec_ptr,ierr);CHKERRQ(ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ivertex)%id = needed_vertices_petsc(ivertex)
    unstructured_grid%vertices(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    unstructured_grid%vertices(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    unstructured_grid%vertices(ivertex)%z = vec_ptr((ivertex-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_new,vec_ptr,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'vertex_coord_new' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'vertex_coord_new' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(vertices_new,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecDestroy(vertices_new,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  call PrintMsg(option,'Setting cell types')
#endif

  allocate(unstructured_grid%cell_type(unstructured_grid%ngmax))

  select case(unstructured_grid%grid_type)
    case(THREE_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        ! Determine number of faces and cell-type of the current cell
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(8)
            unstructured_grid%cell_type(ghosted_id) = HEX_TYPE
          case(6)
            unstructured_grid%cell_type(ghosted_id) = WEDGE_TYPE
          case(5)
            unstructured_grid%cell_type(ghosted_id) = PYR_TYPE
          case(4)
            unstructured_grid%cell_type(ghosted_id) = TET_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call PrintErrMsg(option)
        end select      
      enddo
    case(TWO_DIM_GRID)
      do ghosted_id = 1, unstructured_grid%ngmax
        select case(unstructured_grid%cell_vertices(0,ghosted_id))
          case(4)
            unstructured_grid%cell_type = QUAD_TYPE
          case(3)
            unstructured_grid%cell_type = TRI_TYPE
          case default
            option%io_buffer = 'Cell type not recognized: '
            call PrintErrMsg(option)
        end select
      end do
    case default
      option%io_buffer = 'Grid type not recognized: '
      call PrintErrMsg(option)
  end select
  
end subroutine UGridDecompose

! ************************************************************************** !

function UGridInternConnect(unstructured_grid,grid_x,grid_y,grid_z, &
                                   option)
  ! 
  ! divides connectivity of an
  ! unstructured grid
  ! 


  use Connection_module
  use Option_module
  use Utility_module, only : DotProduct, CrossProduct
  use Geometry_module  

  implicit none

  type(connection_set_type), pointer :: UGridComputeInternConnect
  type(option_type) :: option
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(grid_unstructured_type) :: unstructured_grid

  type(connection_set_type), pointer :: connections
  PetscInt :: nconn, iconn
  PetscInt :: idual, dual_id

  PetscInt, allocatable :: face_to_vertex(:,:)
  PetscInt, allocatable :: cell_to_face(:,:)
  PetscInt, allocatable :: face_to_cell(:,:)
  PetscInt, allocatable :: vertex_to_cell(:,:)
  PetscInt, allocatable :: temp_int(:)
  PetscInt, allocatable :: temp_int_2d(:,:)
  PetscBool, allocatable :: local_boundary_face(:)
  PetscInt :: num_match
  PetscInt :: found_count
  PetscBool :: found
  PetscBool :: match_found
  PetscInt :: face_count
  PetscInt :: count, i
  PetscInt :: iface, iface2, iside
  PetscInt :: face_id, face_id2
  PetscInt :: ghosted_id, ghosted_id2
  PetscInt :: local_id, local_id2
  PetscInt :: cell_id, cell_id2
  PetscInt :: dual_local_id
  PetscInt :: ivertex, ivertex2
  PetscInt :: vertex_id, vertex_id2
  PetscInt :: vertex_ids4(4)
  PetscInt :: nfaces, nfaces2, nvertices, nvertices2, cell_type, cell_type2
  PetscInt :: face_type, face_type2
  PetscBool :: face_found, vertex_found
  
  PetscReal :: v1(3), v2(3), v3(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: vcross(3), magnitude
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  PetscInt :: ivert
  
  type(plane_type) :: plane1, plane2
  type(point3d_type) :: point1, point2, point3, point4
  type(point3d_type) :: point_up, point_dn
  type(point3d_type) :: intercept1, intercept2, intercept

  character(len=MAXSTRINGLENGTH) :: string  
  
  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL* &
           unstructured_grid%ngmax))
  face_to_vertex = 0
  allocate(cell_to_face(MAX_FACE_PER_CELL, &
                        unstructured_grid%ngmax))
  cell_to_face = 0
  allocate(face_to_cell(2,MAX_FACE_PER_CELL* &
                        unstructured_grid%ngmax))
  face_to_cell = 0
  allocate(vertex_to_cell(0:unstructured_grid%max_cells_sharing_a_vertex, &
                          unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  allocate(unstructured_grid%face_to_vertex_natural(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL*unstructured_grid%ngmax))
  unstructured_grid%face_to_vertex_natural = 0

  face_count = 0
  do ghosted_id = 1, unstructured_grid%ngmax
    cell_type = unstructured_grid%cell_type(ghosted_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do iface = 1, nfaces
      face_count = face_count + 1
      cell_to_face(iface,ghosted_id) = face_count
      face_to_cell(1,face_count) = ghosted_id
      call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                      vertex_ids4)
      do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices(vertex_ids4(ivertex),ghosted_id)
          if (face_to_vertex(ivertex,face_count) > 0) then
            unstructured_grid%face_to_vertex_natural(ivertex,face_count) = &
              unstructured_grid%vertex_ids_natural(face_to_vertex(ivertex,face_count))
          endif
      enddo
    enddo
  enddo

  !
  ! Remove duplicate faces:
  !
  ! A cell (cell_id) and Neighboring-Cell (cell_id2) will only share ONE face.
  ! Find the face that cell_id ane cell_id2 share and remove it.
  !
  ! Method:
  !        - Pick i-th face (iface) of cell_id and check if ALL the vertices of
  !          the iface are present in cell_id2. If all the vertices of iface are
  !          not present in cell_id2, move to the next face.
  !        - After finding the iface, now find iface2 in cell_id2 that
  !          corresponds to iface.
  !        - Check to ensure that atleast on face of cell_id is shared
  !          with cell_id2.
  !
  !
  !
  ! NOTE: For a cell_type = WEDGE_TYPE, faces 1-3 have 4 vertices; while
  !       faces 4-5 have 3 vertices
  !
  do local_id = 1, unstructured_grid%nlmax
    ! Selet a cell and find number of vertices
    cell_id = local_id
    ! cell_type is ghosted, but local cells are in the first nlmax entries
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      ! Select a neighboring cell
      ! ghosted neighbors have a negative id
      cell_id2 = &
        abs(unstructured_grid%cell_neighbors_local_ghosted(idual,local_id))
      cell_type2 = unstructured_grid%cell_type(cell_id2)
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      ! Find the number of vertices for neighboring cell
      nfaces2 = UCellGetNFaces(cell_type2,option)
      ! Initialize
      face_found = PETSC_FALSE
      do iface = 1, nfaces
        ! Select a face and find number of vertices forming the face
        face_id = cell_to_face(iface,cell_id)
        nvertices = UCellGetNFaceVertices(cell_type,iface,option)
        do ivertex = 1, nvertices
          ! Select a vertex and initialize vertex_found
          vertex_id = face_to_vertex(ivertex,face_id) ! face_to_vertex is 1-based indexing
          vertex_found = PETSC_FALSE
          do ivertex2 = 1, unstructured_grid%cell_vertices(0,cell_id2)
            vertex_id2 = unstructured_grid%cell_vertices(ivertex2,cell_id2)
            if (vertex_id == vertex_id2) then
              vertex_found = PETSC_TRUE
              exit
            endif
  enddo

  UGridComputeInternConnect => connections

end function UGridComputeInternConnect

! ************************************************************************** !

subroutine UGridDivideVoronoi(unstructured_grid, connection, iface_cell, &
                                   iconn, ghosted_id, option)
  ! 
  ! Computes details of connection (area, dist, etc)
  ! 
  !

 
  
  implicit none
  
  
  
end subroutine UGridPopulateConnection

! ************************************************************************** !

subroutine UGridDivideCoord(unstructured_grid,option, &
                             grid_x,grid_y,grid_z, &
                             x_min,x_max,y_min,y_max,z_min,z_max)
  ! 
  ! Computes coordinates in x,y,z of unstructured grid cells
  !

  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_8(8)
  PetscReal :: centroid(3)
  PetscErrorCode :: ierr 

  do ghosted_id = 1, unstructured_grid%ngmax 
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    centroid = UCellComputeCentroid(unstructured_grid%cell_type(ghosted_id), &
                                    vertex_8,option)
    grid_x(ghosted_id) = centroid(1)
    grid_y(ghosted_id) = centroid(2)
    grid_z(ghosted_id) = centroid(3)
  enddo

  do ivertex = 1, unstructured_grid%num_vertices_local
    if (x_max < unstructured_grid%vertices(ivertex)%x) &
      x_max = unstructured_grid%vertices(ivertex)%x
    if (x_min > unstructured_grid%vertices(ivertex)%x) &
      x_min = unstructured_grid%vertices(ivertex)%x
    if (y_max < unstructured_grid%vertices(ivertex)%y) &
      y_max = unstructured_grid%vertices(ivertex)%y
    if (y_min > unstructured_grid%vertices(ivertex)%y) &
      y_min = unstructured_grid%vertices(ivertex)%y
    if (z_max < unstructured_grid%vertices(ivertex)%z) &
      z_max = unstructured_grid%vertices(ivertex)%z
    if (z_min > unstructured_grid%vertices(ivertex)%z) &
      z_min = unstructured_grid%vertices(ivertex)%z
  enddo
      
end subroutine UGridComputeCoord

! ************************************************************************** !


! ************************************************************************** !

subroutine UGridComputeQuality(unstructured_grid,option)
  ! 
  ! Computes quality of unstructured grid cells
  ! geh: Yes, this is very primitive as mesh quality can be based on any
  ! number of metrics (e.g., see http://cubit.sandia.gov/help-version8/
  ! Chapter_5/Mesh_Quality_Assessment.html).  However, the current edge
  ! length-based formula gives a ballpark estimate.
  ! 


  use Option_module
  use Geometry_module  
  
  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(option_type) :: option

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point3d_type) :: vertex_8(8)
  PetscReal :: quality, mean_quality, max_quality, min_quality
  PetscErrorCode :: ierr

  mean_quality = 0.d0
  max_quality = -1.d20
  min_quality = 1.d20
  
  do local_id = 1, unstructured_grid%nlmax
    ! ghosted_id = local_id on unstructured grids
    ghosted_id = local_id
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    quality = UCellQuality(unstructured_grid%cell_type( &
                           ghosted_id),vertex_8,option)
    if (quality < min_quality) min_quality = quality
    if (quality > max_quality) max_quality = quality
    mean_quality = mean_quality + quality
  enddo

  call MPI_Allreduce(MPI_IN_PLACE,mean_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
  mean_quality = mean_quality / unstructured_grid%nmax

  call MPI_Allreduce(MPI_IN_PLACE,max_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  call MPI_Allreduce(MPI_IN_PLACE,min_quality,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)

  if (OptionPrintToScreen(option)) then
    write(*,'(/," ---------- Mesh Quality ----------", &
            & /,"   Mean Quality: ",es10.2, &
            & /,"   Max Quality : ",es10.2, &
            & /,"   Min Quality : ",es10.2, &
            & /," ----------------------------------",/)') &
              mean_quality, max_quality, min_quality
  endif

end subroutine UGridComputeQuality

! ************************************************************************** !

subroutine UGridMapSideSet(unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)
  ! 
  ! A redirection function for boundary conditions
  ! ghosted cells
  ! 


#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  PetscInt :: face_vertices(:,:)
  PetscInt :: n_ss_faces
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)
  
  Mat :: Mat_vert_to_face
  Vec :: Vertex_vec, Face_vec
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscReal :: real_array4(4)
  PetscInt, allocatable :: boundary_faces(:)
  PetscInt, allocatable :: temp_int(:,:)
  PetscInt :: boundary_face_count
  PetscInt :: mapped_face_count
  PetscInt :: nfaces, nvertices
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: local_id
  PetscInt :: cell_type
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ivertex, cell_id, vertex_id_local
  PetscErrorCode :: ierr
  PetscReal :: min_verts_req
  PetscInt :: largest_vert_id, v_id_n
  Vec :: sideset_vert_vec
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  IS :: is_tmp1, is_tmp2
  VecScatter :: scatter_gton


  ! fill matrix with boundary faces of local cells
  ! count up the number of boundary faces
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    nfaces = UCellGetNFaces(unstructured_grid%cell_type(local_id),option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
      endif
    enddo
  enddo

  call MatCreateAIJ(option%mycomm, &
                       boundary_face_count, &
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       unstructured_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face, &
                       ierr);CHKERRQ(ierr)
  call MatZeroEntries(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  real_array4 = 1.d0

  offset=0
  call MPI_Exscan(boundary_face_count,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(boundary_faces(boundary_face_count))
  boundary_faces = 0
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type,option)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
        boundary_faces(boundary_face_count) = face_id
        call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                        int_array4)

        ! For this matrix:
        !   irow = local face id
        !   icol = natural (global) vertex id
        do ivertex = 1, nvertices
          vertex_id_local = &
            unstructured_grid%cell_vertices(int_array4(ivertex),local_id)
          int_array4_0(ivertex) = &
            unstructured_grid%vertex_ids_natural(vertex_id_local)-1
        enddo
        call MatSetValues(Mat_vert_to_face,1,boundary_face_count-1+offset, &
                          nvertices,int_array4_0,real_array4, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
  enddo

  call MatAssemblyBegin(Mat_vert_to_face,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_subsurf.out'
  else
    string = 'Mat_vert_to_face_' // trim(region_name) // '_global' // &
            '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecCreateMPI(option%mycomm,PETSC_DETERMINE, &
                    unstructured_grid%num_vertices_global, &
                    Vertex_vec,ierr);CHKERRQ(ierr)
  call VecZeroEntries(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyBegin(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(Vertex_vec,ierr);CHKERRQ(ierr)

  ! For this vector:
  !   irow = natural (global) vertex id
  nvertices = 0
  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        nvertices = nvertices + 1
      endif
    enddo
  enddo

  offset=0
  call MPI_Exscan(nvertices,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array(nvertices))
  do local_id = 1, nvertices 
    int_array(local_id) = (local_id-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nvertices, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'is_tmp1_' // trim(region_name) // '_subsurf.out'
  else
    string = 'is_tmp1_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call ISView(is_tmp1,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  nvertices = 0
  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        nvertices = nvertices + 1
        int_array(nvertices) = face_vertices(ivertex,iface)-1
      endif
    enddo
  enddo

  call ISCreateGeneral(option%mycomm,nvertices, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr);CHKERRQ(ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'is_tmp2_' // trim(region_name) // '_subsurf.out'
  else
    string = 'is_tmp2_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call ISView(is_tmp2,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call VecCreateMPI(option%mycomm,nvertices, PETSC_DETERMINE, &
                    sideset_vert_vec,ierr);CHKERRQ(ierr)
  call VecSet(sideset_vert_vec,1.d0,ierr);CHKERRQ(ierr)

  call VecScatterCreate(sideset_vert_vec,is_tmp1, &
                        Vertex_vec,is_tmp2,scatter_gton,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp1,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp2,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'scatter_gton_' // trim(region_name) // '_subsurf.out'
  else
    string = 'scatter_gton_' // trim(region_name) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(scatter_gton,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
  call VecScatterBegin(scatter_gton,sideset_vert_vec,Vertex_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter_gton,sideset_vert_vec,Vertex_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(scatter_gton,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Vertex_vec_' // trim(region_name) // '_global' // &
              '_subsurf.out'
  else
    string = 'Vertex_vec_' // trim(region_name) // '_global' // &
              '_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(Vertex_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  call VecCreateMPI(option%mycomm,boundary_face_count,PETSC_DETERMINE,Face_vec, &
                    ierr);CHKERRQ(ierr)
  call MatMult(Mat_vert_to_face,Vertex_vec,Face_vec,ierr);CHKERRQ(ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (unstructured_grid%grid_type == THREE_DIM_GRID) then
    string = 'Face_vec_' // trim(region_name) // '_global_subsurf.out'
  else
    string = 'Face_vec_' // trim(region_name) // '_global_surf.out'
  endif
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(Face_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  allocate(temp_int(MAX_FACE_PER_CELL,boundary_face_count))
  temp_int = 0
  
  mapped_face_count = 0
  if ( unstructured_grid%grid_type == THREE_DIM_GRID) then
    min_verts_req = 3.d0
  else
    min_verts_req = 2.d0
  endif
  
  call VecGetArrayF90(Face_vec,vec_ptr,ierr);CHKERRQ(ierr)
  ! resulting vec contains the number of natural vertices in the sideset that
  ! intersect a local face
  do iface = 1, boundary_face_count
    face_id = boundary_faces(iface)
    if (vec_ptr(iface) >= min_verts_req) then ! 3 or more vertices in sideset
      ! need to ensure that the right number of vertices are included
      cell_id = unstructured_grid%face_to_cell_ghosted(1,face_id)
      cell_type = unstructured_grid%cell_type(cell_id)
      nfaces = UCellGetNFaces(cell_type,option)
      nvertices = 0
      do iface2 = 1, nfaces
        face_id2 = unstructured_grid%cell_to_face_ghosted(iface2,cell_id)
        if (face_id == face_id2) then
          nvertices = UCellGetNFaceVertices(cell_type,iface2,option)
          exit
        endif
      enddo
      if (nvertices == 0) then ! the case if not found 
        option%io_buffer = 'Face not found in UGridMapSideSet'
        call PrintErrMsgByRank(option)
      endif
      if (abs(nvertices - vec_ptr(iface)) < 0.5d0) then
        mapped_face_count = mapped_face_count + 1
        temp_int(1,mapped_face_count) = cell_id
        temp_int(2,mapped_face_count) = iface2
      endif
    endif
  enddo
  call VecRestoreArrayF90(Face_vec,vec_ptr,ierr);CHKERRQ(ierr)
  deallocate(boundary_faces)
  
  allocate(cell_ids(mapped_face_count))
  allocate(face_ids(mapped_face_count))
  
  cell_ids(:) = temp_int(1,1:mapped_face_count)
  face_ids(:) = temp_int(2,1:mapped_face_count)

  call MatDestroy(Mat_vert_to_face,ierr);CHKERRQ(ierr)
  call VecDestroy(Face_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(Vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(sideset_vert_vec,ierr);CHKERRQ(ierr)
  
end subroutine UGridMapSideSet

! ************************************************************************** !

subroutine DivideVariableArray(data_array, unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)
  
  ! 1.1) Create an array containing cell-ids of 'exisiting' ghost cells + 
  !      ghost cells from 'cids_new'
  !
  count = ngmax_new-unstructured_grid%nlmax + &
          unstructured_grid%ngmax - unstructured_grid%nlmax
  allocate(int_array1(count))
  allocate(int_array2(count))

  count=0
  do ii=1,unstructured_grid%ngmax-unstructured_grid%nlmax
    count=count+1
    ghosted_id=ii+unstructured_grid%nlmax
    int_array1(count)=unstructured_grid%cell_ids_natural(ghosted_id)
    int_array2(count)=count
  enddo
  
  call VecGetArrayF90(cids_on_proc,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,ngmax_new
    if (vec_ptr(ii)/=option%myrank) then
      count=count+1
      int_array1(count)=cids_new(ii)
      int_array2(count)=count
    endif
  enddo
  call VecRestoreArrayF90(cids_on_proc,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(cids_on_proc,ierr);CHKERRQ(ierr)

  ! 1.2) Sort the array
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(count,int_array1, &
                                   int_array2,ierr);CHKERRQ(ierr)
  int_array2 = int_array2+1

  ! 1.3) Count the entries in the sorted array which appear only once.
  nghost_new=0
  ii=1
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) nghost_new=nghost_new+1
  
  do ii=2,count-1
    if ((int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))).and. &
       (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) ) nghost_new=nghost_new+1
  enddo

  ii=count
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))) nghost_new=nghost_new+1
  
  ! 1.4) Save the entries in the sorted array which appear only once.
  allocate(ghost_cids_new(nghost_new))
  nghost_new=0
  ii=1
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) then
    nghost_new=nghost_new+1
    ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
  endif
  
  do ii=2,count-1
    if ((int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))).and. &
       (int_array1(int_array2(ii)) /= int_array1(int_array2(ii+1))) ) then
      nghost_new=nghost_new+1
      ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
    endif
  enddo

  ii=count
  if (int_array1(int_array2(ii)) /= int_array1(int_array2(ii-1))) then
    nghost_new=nghost_new+1
    ghost_cids_new(nghost_new) = int_array1(int_array2(ii))
  endif
  
  deallocate(int_array1)
  deallocate(int_array2)

  ! Step-2: Find PETSc index of additional ghost cells
  call VecCreateMPI(option%mycomm, &
                    unstructured_grid%nlmax, &
                    PETSC_DETERMINE, &
                    cids_petsc,ierr);CHKERRQ(ierr)
  
  offset=0
  call MPI_Exscan(unstructured_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(int_array1(unstructured_grid%nlmax))
  allocate(tmp_scl_array(unstructured_grid%nlmax))
  
  do local_id=1,unstructured_grid%nlmax
    nat_id=unstructured_grid%cell_ids_natural(local_id)
    int_array1(local_id)=nat_id - 1
    tmp_scl_array(local_id)=local_id+offset+0.d0
  enddo
  
  call VecSetValues(cids_petsc,unstructured_grid%nlmax,int_array1,tmp_scl_array,INSERT_VALUES, &
                    ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  deallocate(tmp_scl_array)
  
  call VecAssemblyBegin(cids_petsc,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(cids_petsc,ierr);CHKERRQ(ierr)

  call VecCreateMPI(option%mycomm,nghost_new,PETSC_DETERMINE,ghosts_petsc, &
                    ierr);CHKERRQ(ierr)
  allocate(int_array1(nghost_new))

  int_array1=ghost_cids_new-1
  call ISCreateGeneral(option%mycomm,nghost_new, &
                       int_array1,PETSC_COPY_VALUES,is_from, &
                       ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(nghost_new,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii=1,nghost_new
    int_array1(ii)=ii-1+offset
  enddo
  call ISCreateGeneral(option%mycomm,nghost_new, &
                       int_array1,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
  deallocate(int_array1)
  
  call VecScatterCreate(cids_petsc,is_from,ghosts_petsc,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,cids_petsc,ghosts_petsc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,cids_petsc,ghosts_petsc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)
  
  allocate(ghost_cids_new_petsc(nghost_new))
  call VecGetArrayF90(ghosts_petsc,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,nghost_new
    ghost_cids_new_petsc(ii)=INT(vec_ptr(ii))
  enddo
  call VecRestoreArrayF90(ghosts_petsc,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(cids_petsc,ierr);CHKERRQ(ierr)
  call VecDestroy(ghosts_petsc,ierr);CHKERRQ(ierr)

end subroutine DivideVariableArray

! ************************************************************************** !

end module decomp_unstruct_module



! ***********
! TO CALL FROM FEHM_INITPETSC
!  cell distribution across processors (size = num_cores + 1)
!  core i owns cells cell_distribution(i):cell_distribution(i+1), note
!  the zero-based indexing
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%nlmax,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)

  num_cells_local_old = unstructured_grid%nlmax

  ! recalculate maximum number of vertices for any given cell
  temp_int = 0
  min_value = 2 ! min value should be either 0 or 1 after global reduction
  do local_id = 1, num_cells_local_old
    vertex_count = 0
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point, cell vertex can be 0
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) exit
      if (unstructured_grid%cell_vertices(ivertex,local_id) < min_value) then
        min_value = unstructured_grid%cell_vertices(ivertex,local_id)
      endif
      vertex_count = vertex_count+1
    enddo
    if (vertex_count > temp_int) temp_int = vertex_count
  enddo
  call MPI_Allreduce(temp_int,unstructured_grid%max_nvert_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(min_value,index_format_flag, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN,option%mycomm,ierr)

  ! let's make it Fortran indexing
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point we may be zero-based
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) then
        ! change no_value (UNINITIALIZED_INTEGER) to '0'
        unstructured_grid%cell_vertices(ivertex,local_id) = 0
      else
        if (index_format_flag == 0) then
          ! let's make it Fortran indexing
          unstructured_grid%cell_vertices(ivertex,local_id) = &
            unstructured_grid%cell_vertices(ivertex,local_id) + 1
        endif
      endif
    enddo
  enddo

#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_nvert_per_cell
  option%io_buffer = 'Maximum number of vertices per cell: ' // adjustl(string)
  call PrintMsg(option)
  write(string,*) index_format_flag
  option%io_buffer = 'Vertex indexing starts at: ' // adjustl(string)
  call PrintMsg(option)
  if (index_format_flag == 0) then
    option%io_buffer = 'Changing vertex indexing to 1-based.'
    call PrintMsg(option)
  endif
#endif

  num_cells_local_old = unstructured_grid%nlmax 
  allocate(local_vertices(unstructured_grid%max_nvert_per_cell* &
                          num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  local_vertices = 0
  local_vertex_offset = 0
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      if (unstructured_grid%cell_vertices(ivertex,local_id) == 0) exit
      count = count + 1
      ! local vertices must be zero-based for MatCreateMPIAdj; thus subtract 1
      local_vertices(count) = &
        unstructured_grid%cell_vertices(ivertex,local_id) - 1
    enddo
    local_vertex_offset(local_id+1) = count 
  enddo
    
  select case (unstructured_grid%grid_type)
    case(TWO_DIM_GRID)
      num_common_vertices = 2 ! cells must share at least this number of vertices
    case(THREE_DIM_GRID)
      num_common_vertices = 3 ! cells must share at least this number of vertices
    case default
        option%io_buffer = 'Grid type not recognized '
        call PrintErrMsg(option)
    end select

  ! determine the global offset from 0 for cells on this rank
  global_offset_old = 0
  call MPI_Exscan(num_cells_local_old,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! create an adjacency matrix for calculating the duals (connnections)
#if UGRID_DEBUG
  call PrintMsg(option,'Adjacency matrix')
#endif

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat, &
                       ierr);CHKERRQ(ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

  endif
  
  
!!!!!!!!!!!!!!!!!!!!!!  
! CALL UGridDecompose()
	!!!!!!!!!!!!!!!!!!
#endif