      subroutine rlp_cap_table(sw, itbl, iparam, cp, dpcp)
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

      integer i, itbl, i1, i2, iparam
      real*8 sw, cp, dpcp

      i1 = tblindx(itbl , 1)
      i2 = tblindx(itbl , 2)

      do i = i1, i2 - 1
         if (sw .le. rlp_table(i, 1) .and. i .eq. i1) then
            cp = rlp_table(i1, iparam)
            dpcp = 0.
            return
         else if (sw .ge. rlp_table(i + 1, 1) .and. i + 1 .eq. i2) then
            cp = rlp_table(i2, iparam)
            dpcp = 0.
            return
         else if (sw .ge. rlp_table(i, 1) .and. 
     &           sw .lt. rlp_table(i + 1, 1)) then
            dpcp = (rlp_table(i + 1, iparam) - rlp_table(i, iparam)) /
     &           (rlp_table(i + 1, 1) - rlp_table(i, 1))
            cp = rlp_table(i, iparam) + dpcp * (sw - rlp_table(i, 1))
            return
         end if
      end do

      end subroutine rlp_cap_table
