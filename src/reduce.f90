subroutine reduce_to_tridiagonal(A)
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
!********************************************************************
! SUBROUTINE reduce_to_tridiagonal                                  *
!                                                                   *
! Accepts a 3x3 matrix A and uses a single Householder matrix to    *
! reduce A to tridiagonal form. The matrix A is expected to be a    *
! real matrix with symmetric components.                            *
!                                                                   *
! Author : Sai Rapaka                                               *
!                                                                   *
!********************************************************************


implicit none
double precision, dimension(3,3) :: A, P, B
double precision, dimension(2)   :: x
double precision                 :: norm_x, H
double precision, parameter      :: tol = 1.0d-12
integer i,j

x(1) = A(2,1)
x(2) = A(3,1)

norm_x = sqrt(x(1)*x(1) + x(2)*x(2))
x(1) = x(1) - norm_x

H = 0.5d0*(x(1)*x(1) + x(2)*x(2))

if(H.gt.tol) then
  P = 0.0d0
  P(1,1) = 1.0d0
  do i=2,3
    do j=2,3
      if(i.eq.j) then
        P(i,j) = 1.0d0
      endif
      P(i,j) = P(i,j) - x(i-1)*x(j-1)/H
    enddo
  enddo

  B = MATMUL(P,A)
  A = MATMUL(B,P)

endif

end subroutine reduce_to_tridiagonal

