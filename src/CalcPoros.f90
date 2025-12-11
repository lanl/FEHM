subroutine CalcPoros(NodeNum,Temp,Pres,Poros)
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
! Input:node number, temperature (C), pressure(MPa)
! Output: porosity, Poros

use comdi, only: InitPoros,phini,LiqEndTemp,IceEndTemp
implicit none
real*8 TLower,TUpper,Pres,Alpha,Poros,Temp,DeltaPres,DeltaPoros, &
   RelRampPos,MinPoros,MinPres,AquiferCompress
integer :: NodeNum

MinPoros = 0.3D0
MinPres = 0.1D0
AquiferCompress = 1.D-1

DeltaPres = phini(NodeNum) - MinPres
if (DeltaPres.ne.0.0) then
   Alpha = (InitPoros(NodeNum)-MinPoros)/DeltaPres
end if

RelRampPos = (Temp - IceEndTemp)/(LiqEndTemp - IceEndTemp)
RelRampPos = max(0.0d0,RelRampPos)
RelRampPos = min(1.0d0,RelRampPos)

if (InitPoros(NodeNum) > MinPoros) then
   Poros = (1.0 - RelRampPos)*(InitPoros(NodeNum) + &
      AquiferCompress*(Pres - phini(NodeNum))) +    &
      RelRampPos*(MinPoros + Alpha*(Pres - MinPres))
else
   Poros = InitPoros(NodeNum) + &
      AquiferCompress*(Pres - phini(NodeNum))
end if

Poros = max(MinPoros,Poros)
Poros = min(1.0d0,Poros)
end
