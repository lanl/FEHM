      subroutine stone_rlp(rlc, drlc, rw, rg, rl, drls, rv, drvs, r, dr)
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
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
!D1
!D1 PURPOSE
!D1
!D1 Calculate relative permeabilities and capillary pressures
!D1 and derivatives.
!D1
!***********************************************************************

      use comrlp
      implicit none

      real*8 rlc, drlc, rw, rg, rl, drls, rv, drvs, r, dr
      
      r = rlc * ((rl/rlc + rw)*(rv/rlc + rg) - rw -rg)
      

      end subroutine stone_rlp
