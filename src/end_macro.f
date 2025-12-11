      logical function end_macro(macro, last_macro, line)
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
!D1
!D1 PURPOSE
!D1
!D1 Check macro end designator
!D1
!***********************************************************************

      use comai, only : ierr, iptty
      implicit none

      integer l1
      character*4 macro, last_macro
      character*8 e1
      character*80 line

c     Check end macro

      if (len_trim(macro) .ne. 3) then
         e1 = trim(macro) // last_macro
      else
         e1 = macro // last_macro
      end if
      l1 = len_trim(e1)
      if (line(1:l1) .eq. e1(1:l1)) then
c     macro terminator is OK
         end_macro = .true.
      else
         end_macro = .false.
      end if

      end function end_macro
