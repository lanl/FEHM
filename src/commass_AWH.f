      module commass_AWH
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

c 
c gaz 022521
c global variables for conservation phase change
c
      
      integer  imass_phase,neq_phase_chk,n_phase_chk, ivar_mass
c      real*8 r1min,r2min,r3min,r123min,r123tol
c      integer, allocatable :: nphase_chk(:)
      real*8, allocatable  :: wmass_awh(:)
      real*8, allocatable  :: energy_awh(:)
      real*8, allocatable  :: amass_awh(:)
      real*8, allocatable  :: r1_mass(:)
      real*8, allocatable  :: r2_energy(:)
      real*8, allocatable  :: r3_amass(:)
      
      end module commass_AWH
