subroutine GetSubsWaveLength(SubsWaveLength)
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

use comdi, only: izone_top,z_plot
use comai, only: neq
use combi, only: cord
implicit none
real*8 :: SubsWaveLength,x_min,x_max,x_shift,x_scale,Integral, &
   a_max,DeltaChi,AvgSubs,MaxSubs,MinSubs
real*8, allocatable :: Chi(:),Subsid(:),a(:)
integer :: I,NumNodes,Counter,n,MaxCoeffNum
integer, parameter :: NumCoeffs = 100
real*8, parameter :: pi = 3.14159

! top boundary is indicated by the neq entries of izone_top:
! izone_top(I) = 3 indicates that node I is on the top boundary;
! izone_top(I) = 0 indicates that node I is not on the top 
! boundary

! cord(I,1) = x coordinates, cord(I,2) = y coordinatess, 
! cord(I,3) = z coordinates

! final heights are in the array z_plot, which also has neq 
! entries

! shift and rescale the x coords so that they lie between
! 0 and 1; convert resulting wavelength back to meters
! afterward

! find x_min, x_max, and # of nodes above the subsidence zone
x_min =  1.e30
x_max = -1.e30
NumNodes = 0
do I=1,neq
   if (izone_top(I)==3) then
      if (cord(I,1)<x_min) then
         x_min = cord(I,1)
      end if
      if (cord(I,1)>x_max) then
         x_max = cord(I,1)
      end if
      NumNodes = NumNodes + 1
   end if
end do
x_shift = -x_min
x_scale = (x_max - x_min)

allocate(Chi(NumNodes))
allocate(Subsid(NumNodes))
Counter = 1
AvgSubs = 0.0
MaxSubs = -1.e30
MinSubs =  1.e30
do I=1,neq
   if (izone_top(I)==3) then
      Chi(Counter) = (cord(I,1) + x_shift)/x_scale
      Subsid(Counter) = z_plot(I)
      if (Subsid(Counter)>MaxSubs) then
         MaxSubs = Subsid(Counter)
      end if
      if (Subsid(Counter)<MinSubs) then
         MinSubs = Subsid(Counter)
      end if
      !write(*,*) Counter,'XCoord: ',cord(I,1),'Chi: ',Chi(Counter)
      !write(*,*) Counter,'subsid: ',Subsid(Counter)
      Counter = Counter + 1
   end if
end do
AvgSubs = (MaxSubs + MinSubs)/2.0
! vertically shift the subsidence to center on AvgSubs
do I=1,Counter-1
   Subsid(I) = Subsid(I) - AvgSubs
end do

! calculate the fourier series coefficients a_n, find max | a_n |
a_max = 0.0
allocate(a(NumCoeffs))
do n=1,NumCoeffs
   Integral = 0.0
   do I=1,Counter-2
      DeltaChi = Chi(I+1)-Chi(I)
      Integral = Integral + sin(n*pi*Chi(I))*Subsid(I)*DeltaChi
   end do
   a(n) = 2.0*Integral
   !write(*,*) n,'| a |: ',dabs(a(n))
   if (dabs(a(n))>a_max) then
      a_max = dabs(a(n))
      MaxCoeffNum = n
   end if
end do
!write(*,*) 'max coeff number: ',MaxCoeffNum
!write(*,*) 'max coeff: ',a(MaxCoeffNum)
SubsWaveLength = 2.0*x_scale/MaxCoeffNum
end
