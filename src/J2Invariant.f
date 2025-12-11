      subroutine J2Invariant(J2, dev_stress)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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
!
! Returns the J2-Invariant for a six-component deviatoric stress tensor
!
! Author : Sai Rapaka
!

      implicit none
      real*8,  dimension(6)        :: dev_stress
      real*8 J2
      real*8 pressure

      pressure = (dev_stress(1) + dev_stress(2) + dev_stress(3))/3.0d0
      dev_stress(1) = dev_stress(1) - pressure
      dev_stress(2) = dev_stress(2) - pressure
      dev_stress(3) = dev_stress(3) - pressure

      J2 = 0.0d0
      J2 = 0.5d0*sum(dev_stress(1:3)**2) + sum(dev_stress(4:6)**2)

      end subroutine J2Invariant
