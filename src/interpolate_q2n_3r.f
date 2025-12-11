      subroutine interpolate_q2n_3r(val, f, k)
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
! Extrapolates a variable defined at the gausspoints to the nodes
! 
! Author : Sai Rapaka
!

      use comfem

      implicit none
      real*8,  dimension(8)        :: f
      real*8                       :: val
      integer                      :: k
      
      real*8                       :: sq3, zeta, eta, mu
      real*8,  dimension(8)        :: N
      integer                      :: i

      sq3 = sqrt(3.0d0)
      
      zeta = gpcord(k,1)*sq3*sq3
      eta  = gpcord(k,2)*sq3*sq3
      mu   = gpcord(k,3)*sq3*sq3

      N(1) = 0.125d0*(1 - zeta)*(1 - eta)*(1 - mu)
      N(2) = 0.125d0*(1 + zeta)*(1 - eta)*(1 - mu)
      N(3) = 0.125d0*(1 + zeta)*(1 + eta)*(1 - mu)
      N(4) = 0.125d0*(1 - zeta)*(1 + eta)*(1 - mu)
      N(5) = 0.125d0*(1 - zeta)*(1 - eta)*(1 + mu)
      N(6) = 0.125d0*(1 + zeta)*(1 - eta)*(1 + mu)
      N(7) = 0.125d0*(1 + zeta)*(1 + eta)*(1 + mu)
      N(8) = 0.125d0*(1 - zeta)*(1 + eta)*(1 + mu)
      
      val = 0.0d0
      do i=1,8
        val = val + N(i)*f(i)
      enddo

      end subroutine interpolate_q2n_3r
 
