      subroutine connections_list
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
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
!D1 Print out number of connections to a node after simplification 
!D1 of connectivity based on porosity and terminate execution after
!D1 list is output.
!D1
!***********************************************************************

      use comai, only : neq, ierr, iout, iptty
      use combi, only : nelm

      implicit none

      integer i, i1, i2, num, list_unit, open_file

      list_unit = open_file('connections_list.txt', 'unknown')

      do i = 1, neq
         i1 = nelm(i) + 1
         i2 = nelm(i + 1)
         num = i2 - i1
         write (list_unit, '(i8.8, i5)') i, num
      end do

      write (ierr, 10)
      if (iout .ne. 0) write (iout, 10)
      if (iptty .ne. 0) write (iptty, 10)

      stop

 10    format ('Number of nodal connections written to file: ',
     &     'connections_list.txt', /, ' **** STOPPING ****')

      end subroutine connections_list
