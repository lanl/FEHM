      subroutine skip_macro(macro, inunit, found_end)
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
!D1 Don't use macro for current simulation
!D1
!***********************************************************************

      use comai, only : ierr, iout, iptty
      implicit none
      
      integer inunit, nwds, imsg(20), msg(20)

      real*8 xmsg(20)
      character*4 macro, macroend
      character*32 cmsg(20)
      character*80 line
      logical :: found_end, end_macro

      do
c     Find the end of input for this macro

         read (inunit, '(a80)', end = 100) line
         call parse_string2(line,imsg,msg,xmsg,cmsg,nwds)

         if (msg(1) .eq. 3) then
            macroend = cmsg(1)
            if (macroend(1:3) .eq. 'end' .or. 
     &           macroend(1:3) .eq.  'END') then
               found_end = end_macro (macroend, macro, line)
               if (found_end) return
            end if
         end if

      end do

 100  write (ierr, 200) trim(macro)
      if (iout .ne. 0) write (iout, 200) trim(macro)      
      if (iptty .ne. 0) write (iptty, 200) trim(macro)

      stop

 200  format ('**** Missing "end ', a, '" statement. STOP ****')

      end subroutine skip_macro
