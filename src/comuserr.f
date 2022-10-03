      module comuserr
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
!D1 Global file for reaction user subroutine (userr)
!D1
!***********************************************************************

c      integer iuchk, iuserc, nu_zones, nu_points
c      integer :: usroption = 0
c      integer, allocatable :: prodzones(:)
c      integer, allocatable :: inflzones(:)
c      integer, allocatable :: nodeindex(:)
c      integer :: iaunit = 0

c      real(8) :: in1save
c      real(8) :: erosion_factor = 0.d0
c      real(8) :: timeCl36start
c      real(8), allocatable :: userconc(:,:)
c      real(8), allocatable :: userconc3(:,:,:)
c      real(8), allocatable :: moles_recycle(:,:)
c      real(8), allocatable :: recycle_factor(:)

c      real(8), allocatable :: in(:)

c      integer num_models, num_times
      integer, allocatable :: distcoeff_zones(:,:)
      integer, allocatable :: rxnon_u(:,:)
      real(8), allocatable :: distcoeffs(:,:)
      real(8), allocatable :: userrtime(:,:)
      real(8), allocatable :: distcoeffi(:)
      


      end module comuserr

