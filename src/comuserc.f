      module comuserc
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
!D1 Global file for transport user subroutine
!D1
!***********************************************************************

      integer iuchk, iuserc, nu_zones, nu_points
      integer :: usroption = 0
      integer, allocatable :: prodzones(:)
      integer, allocatable :: inflzones(:)
      integer, allocatable :: nodeindex(:)
      integer :: iaunit = 0

      real(8) :: in1save
      real(8) :: erosion_factor = 0.d0
      real(8) :: timeCl36start
      real(8), allocatable :: userconc(:,:)
      real(8), allocatable :: userconc3(:,:,:)
      real(8), allocatable :: usertime(:)
      real(8), allocatable :: moles_recycle(:,:)
      real(8), allocatable :: recycle_factor(:)

      real(8), allocatable :: in(:)

      end module comuserc
