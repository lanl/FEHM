!********************************************************************
! SUBROUTINE eigenvectors                                           *
!                                                                   *
! Accepts as input a symmetric 3x3 tridiagonal matrix A and a       *
! vector of eigenvalues lambda and computes the eigenvectors        *
! corresponding to each eigenvalue. The eigenvalues are sorted      *
! during the computation in increasing order.                       *
!                                                                   *
!********************************************************************
subroutine eigenvectors(A, lambda, V)

implicit none
double precision, dimension(3,3)       :: A, V, E, B
double precision, dimension(3)         :: lambda
double precision                       :: tmp, theta, phi

double precision, parameter            :: cut_off = 1.0d-4
double precision, parameter            :: zero = 0.0d0
integer i,j
logical e1used, e3used

! E is the identity matrix
E = reshape((/1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/), SHAPE=(/3,3/))

V = 0.0d0
e1used = .false.
e3used = .false.

! First deal with the smallest eigenvalue
B = A - lambda(1)*E
! Make sure that numerically small values in B are truncated to 0.0d0
! Needed to ensure robustness of calls to atan2
do i=1,3
  do j=1,3
    if(abs(B(i,j)).lt.cut_off) then
      B(i,j) = zero
    endif
  enddo
enddo

if((B(2,3).eq.zero).and.(B(3,3).eq.zero)) then
  V(1,1) = 0.0d0
  V(2,1) = 0.0d0
  V(3,1) = 1.0d0
  e3used = .true.
elseif((B(1,1).eq.zero).and.(B(1,2).eq.zero)) then
  V(1,1) = 1.0d0
  V(2,1) = 0.0d0
  V(3,1) = 0.0d0
  e1used = .true.
else
  phi = datan2(-B(1,1), B(1,2))
  theta = datan2(-B(3,2)*dsin(phi), B(3,3))
  V(1,1) = dcos(theta)*dcos(phi)
  V(2,1) = dcos(theta)*dsin(phi)
  V(3,1) = dsin(theta)
endif

! Now deal with the largest eigenvalue
B = A - lambda(3)*E
! Make sure that numerically small values in B are truncated to 0.0d0
! Needed to ensure robustness of calls to atan2
do i=1,3
  do j=1,3
    if(abs(B(i,j)).lt.cut_off) then
      B(i,j) = zero
    endif
  enddo
enddo
if((B(2,3).eq.zero).and.(B(3,3).eq.zero)) then
  if(.not.e3used) then
    V(1,3) = 0.0d0
    V(2,3) = 0.0d0
    V(3,3) = 1.0d0
  else
    theta = datan2(-B(1,1), B(1,2))
    V(1,3) = dcos(theta)
    V(2,3) = dsin(theta)
    V(3,3) = 0.0d0
  endif
elseif((B(1,1).eq.zero).and.(B(1,2).eq.zero)) then
  if(.not.e1used) then
    V(1,3) = 1.0d0
    V(2,3) = 0.0d0
    V(3,3) = 0.0d0
  else
    theta = datan2(-B(3,2), B(3,3))
    V(1,3) = 0.0d0
    V(2,3) = dcos(theta)
    V(3,3) = dsin(theta)
  endif
else
  phi = datan2(-B(1,1), B(1,2))
  theta = datan2(-B(3,2)*dsin(phi), B(3,3))
  V(1,3) = dcos(theta)*dcos(phi)
  V(2,3) = dcos(theta)*dsin(phi)
  V(3,3) = dsin(theta)
endif

! Obtain the eigenvector for the intermediate eigenvalue using
! cross product of the other two eigenvectors
V(1,2) = V(2,1)*V(3,3) - V(3,1)*V(2,3)
V(2,2) = V(3,1)*V(1,3) - V(1,1)*V(3,3)
V(3,2) = V(1,1)*V(2,3) - V(2,1)*V(1,3)

end subroutine eigenvectors

