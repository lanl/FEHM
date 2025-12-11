function PorosCoeff(NodeNum,Temp,Pres,ToDrain)
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
!function PorosCoeff(NodeNum,Temp,Pres,ToDrain) result(Answer)
! Input:node number, temperature (C), pressure(MPa), increase in 
!       pressure over last timestep, initial and current porosity, 
!       and whether to allow for collapse after fluid drainage (1 
!       for yes and 0 for no)
! Output: factor by which to increase the porosity due to the 
!         melting of ice and drainage (result will be < 1)

use comdi, only: MaxDrainPart,InitPoros,InitPres
implicit none
real*8 Answer,VolumeRatio,TLower,TUpper,RelRampPos,Temp,PhasePart,   &
   DrainFormula,Pres,DrainCoeff,PorosCoeff
! RefPres is the pressure at which the porosity reaches a minimum
! when drainage is turned on; MinPoros is the minimum allowed
! porosity
real*8, parameter :: RefPres = 0.1D0, MinPoros = 0.3
integer :: ToDrain,NodeNum,DrainFlag

! liquid H2O to ice volume ratio
VolumeRatio = 0.9

! lower and upper bounds of the phase change ramp
Tlower = -0.25
Tupper =  0.25

! fraction of distance up ramp relative to Tlower
RelRampPos = (Temp-Tlower)/(Tupper-Tlower)
if (RelRampPos > 1.0) then
   RelRampPos = 1.0
elseif (RelRampPos < 0.0) then
   RelRampPos = 0.0
end if

if ((InitPres(NodeNum)-RefPres).ne.0.0) then
   DrainFormula = MinPoros + &
      (VolumeRatio*InitPoros(NodeNum) - MinPoros)* &
      (Pres-RefPres)/(InitPres(NodeNum)-RefPres)
else
   DrainFormula = VolumeRatio*InitPoros(NodeNum)
end if

PorosCoeff = (1.0-RelRampPos) + &
         RelRampPos*DrainFormula/InitPoros(NodeNum)
! for analytic solution only!
!Answer = (MinPoros + 3.0*(Pres-RefPres))/0.6
end
