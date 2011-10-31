function DeltaHeight(NodeNum,InitHeight,Poros) result(Answer)
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
! Input: node number, initial control volume height (m),
!        and current porosity
! Output: increase in the control volume height (m) due to the 
!         melting of ice (result will be negative)

use comdi, only: InitPoros
implicit none
real*8 :: Answer,InitHeight,Poros
integer :: NodeNum

Answer = InitHeight*((1.0-InitPoros(NodeNum))/(1.0-Poros)-1.0)
end
