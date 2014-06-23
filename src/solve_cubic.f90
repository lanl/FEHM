subroutine solve_cubic(a,b,c,lambda)
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
! SUBROUTINE solve_cubic(a,b,c,lambda)                              *
!                                                                   *
! Accepts the coefficients a,b and c of the monic cubic polynomial  *
! x**3 + a*x**2 + b*x + c = 0 and returns its three roots in the    *
! vector lambda. It is assumed that the roots are real. The         *
! intended application is the computation of eigenvalues of real,   *
! symmetric tensors - which justifies the assumption.               *
!                                                                   *
! Author : Sai Rapaka                                               *
!                                                                   *
!********************************************************************


implicit none
integer ii,j
double precision, dimension(3)           :: lambda
double precision                         :: a,b,c, tmp
double precision                         :: aa, bb
double precision                         :: P, Q, D
complex                                  :: i, zeta, t1, t2

aa = a/3.0d0
bb = b/3.0d0

P = bb - aa*aa
Q = c/2.0d0 + aa**3 - 1.5d0*aa*bb
D = Q**2 + P**3

i = cmplx(0.0d0, 1.0d0)
zeta = -0.5d0 + i*sqrt(3.0d0)/2.0d0
t1 = (Q - i*sqrt(abs(D)))**(1.0d0/3.0d0)
t2 = (Q + i*sqrt(abs(D)))**(1.0d0/3.0d0)

lambda(1) = -real(aa + t1 + t2)
lambda(2) = -real(aa + zeta*t1 + zeta*zeta*t2)
lambda(3) = -real(aa + zeta*zeta*t1 + zeta*t2)

! First sort the eigenvectors into increasing order
do ii=1,3
  do j=ii+1,3
    if(lambda(ii).gt.lambda(j)) then
      tmp = lambda(ii)
      lambda(ii) = lambda(j)
      lambda(j) = tmp
    endif
  enddo
enddo

end subroutine solve_cubic
