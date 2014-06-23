function KozenyCoeff(NodeNum,Poros) result(Answer)
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
! Input: node number and current porosity
! Output: multiplier M such that new_perm = M*old_perm

use comdi, only: InitPoros,phi,InitPerm
implicit none
real*8 :: Answer,Poros,Const,Beta,RefPerm,RefPres,Pres
integer :: NodeNum

Const = (1.0-InitPoros(NodeNum))**2.0/(InitPoros(NodeNum)**3.0)
Answer = Const*Poros**3.0/(1.0-Poros)**2.0

!for simplified solution only!
Answer = 1.0
end
