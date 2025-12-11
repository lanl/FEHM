      module comrlp
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global file for relative permeability and capillary pressure.
!D1
!***********************************************************************

      integer :: maxrp = 0, maxcp = 0, ntable =0, ntblines=0
      integer :: max_rp, max_rpf, max_cp
      integer :: ishisrlp, num_sat,phase_comp(30,3),nphases(10)
      integer, allocatable :: cap_coupling(:,:)
      integer, allocatable :: cap_pos(:,:)
      integer, allocatable :: cap_type(:), cap_type2(:,:)
      integer, allocatable :: rlp_group(:)
      integer, allocatable :: rlp_phase(:,:)
      integer, allocatable :: rlp_pos(:,:)
      integer, allocatable :: rlp_type(:), rlp_type2(:,:)
      integer, allocatable :: tblindx(:,:)

      real(8) :: delta_sat
      real(8), allocatable :: cap_param(:,:)      
      real(8), allocatable :: cap_fparam(:,:)      
      real(8), allocatable :: rlp_param(:,:)
      real(8), allocatable :: rlp_fparam(:,:)
      real(8), allocatable :: vg1(:,:), vg2(:,:), vg3(:,:), vg4(:,:)
      real(8), allocatable :: vg1f(:,:), vg2f(:,:), vg3f(:,:), vg4f(:,:)
      real(8), allocatable :: cp1(:,:), cp2(:,:), cp1f(:,:), cp2f(:,:)
      real(8), allocatable :: vg5(:,:), vg6(:,:), vg7(:,:), vg8(:,:)
      real(8), allocatable :: cp1s(:,:), cp2s(:,:), cp3s(:,:), cp4s(:,:)
      real(8), allocatable :: rlp_table(:,:)
      real(8), allocatable :: sat_out(:)
c gaz 070423
      character*30, allocatable :: table_file_names(:)
      logical :: rlpnew = .false.

      parameter(max_rp = 4, max_rpf = 8, max_cp = 6)


      end module comrlp
