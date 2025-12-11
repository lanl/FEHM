      module com_exphase
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
c gaz 122720
c global variables for explicit phase change
c
      integer i_ex_update,ieq_ex,iter_expa
      integer  neq_phase_chk,iphase_chk
      real*8 r1min,r2min,r3min,r123min,r123tol
      integer, allocatable :: nphase_chk(:)
      real*8, allocatable  :: r1(:)
      real*8, allocatable  :: r2(:)
      real*8, allocatable  :: r3(:)
      real*8, allocatable  :: r1_wo_acc(:)
      real*8, allocatable  :: r2_wo_acc(:)
      real*8, allocatable  :: r3_wo_acc(:)
      real*8, allocatable  :: r_ex(:)
      real*8, allocatable  :: a_ex(:,:)
      
      end module com_exphase
