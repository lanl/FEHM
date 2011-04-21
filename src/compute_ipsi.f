      subroutine compute_ipsi()
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
! Computes the inverse of the shapefunction matrix for hax elements
!
! Author : Sai Rapaka
!
      use comfem, only: iPsi

      implicit none

      real*8   :: xc(8), yc(8), zc(8)
      real*8   :: zeta, eta, mu
      real*8   :: xmult, solved
      integer i,j,k
      real*8   :: Nmat(8, 8), rhs(8, 8)
      real*8   :: sol(8)

      xc = 1.0d0
      yc = 1.0d0
      zc = 1.0d0
      xc=(/-1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)
      yc=(/-1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0/)
      zc=(/-1.0d0, -1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
      do i=1,8
        do j=1,8
          zeta = xc(i)/dsqrt(3.0d0)
          eta = yc(i)/dsqrt(3.0d0)
          mu = zc(i)/dsqrt(3.0d0)
          Nmat(i, j) = (1.0d0/8.0d0)*(1 + xc(j)*zeta)*
     &                 (1 + yc(j)*eta)*(1 + zc(j)*mu)
        enddo
      enddo

      rhs = 0.d00
      do i=1,8
        rhs(i,i) = 1.0d0
      enddo

      do k=1,7
        do i=k+1,8
          xmult = Nmat(i,k)/Nmat(k,k)
          do j=1,8
            Nmat(i,j) = Nmat(i,j) - xmult*Nmat(k,j)
            rhs(i,j) = rhs(i,j) - xmult*rhs(k,j)
          enddo
        enddo
      enddo

      do k=1,8
        ! Solve for the k-th column of rhs
        sol(8) = rhs(8, k)/Nmat(8,8)
        do j=7,1,-1
          solved = 0.0d0
          do i=j+1,8
            solved = solved + sol(i)*Nmat(j,i)
          enddo
          sol(j) = (rhs(j, k) - solved)/Nmat(j, j)
        enddo
        do j=1,8
          iPsi(j, k) = sol(j)
        enddo
      enddo

      return
      end subroutine compute_ipsi

