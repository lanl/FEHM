      program main
!***********************************************************************
!  Copyright, 1994, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Finite Element Heat and Mass Transfer in porous media.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 23-AUG-94    L. Trease      22      Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/mainrip.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:18   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   Tue Jan 30 15:27:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.0   08/23/94 15:25:30   llt
CD2 Added program main above fehmn, so that all arrays would be allocated
CD2 before being used - required for the ibm.
CD2
CD2 1980         G. Zyvoloski   22        Initial implementation.
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This routine mainly calls the driver routine of the code.
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.8 Provide Multiple Realization Option
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN fehmn
CPS 
CPS   call fehmn
CPS
CPS END fehmn
CPS
C***********************************************************************

c	PC Version
C       use ifport

C     not used V3.6
C     use compart, only : ripfehm
C     use comai, only : irun
      implicit none

      interface
         subroutine fehmn(method, state,in, out)
C!DEC$ ATTRIBUTES c :: fehmn
C!DEC$ ATTRIBUTES value :: method
C!DEC$ ATTRIBUTES reference :: state
C!DEC$ ATTRIBUTES reference :: in
C!DEC$ ATTRIBUTES reference :: out
         integer method, state
         real(8) in(4), out(3)
         end subroutine fehmn
      end interface
      integer return_flag
      integer method, state
      real(8) in(4), out(3)
      integer nsim, irun
      character(100) pre_string, post_string
      logical file_exists

c     temporary variables for testing
      real*8 oldtime
      integer its, nts, i, tindex

c zvd 09-Sep-2011 change size of in array to be consistent with iofile
c modification for GoldSim 
      in = 0.
      out = 0.
c      ripfehm = 0
c     Determine how many runs to perform
      call inmsim

c     Perform all runs
      do irun = 1, nsim

         if(file_exists) then
c	UNIX version
            write(pre_string,1000) 'sh fehmn.pre ',irun, nsim
            call system(pre_string)
c	PC version
C	     write(pre_string,1000) 'fehmn.pre ',irun, nsim
C            return_flag = system(pre_string)
         end if
c     First call for initialization, second call for calculation
         method = 0
         in(1) = irun
         call fehmn(method, state, in, out)
         in(1) = 0
         in(3) = 0
         method = 1
         call fehmn(method, state, in, out)
         if(file_exists) then
c	PC Version
C            write(post_string,1001) 'fehmn.post ',irun, nsim
C            return_flag = system(post_string)
c	UNIX Version
            write(post_string,1001) 'sh fehmn.post ',irun, nsim
            call system(post_string)
         end if
      end do


c	UNIX Version
 1000 format(a13, 1x, i10, 1x, i10)
 1001 format(a14, 1x, i10, 1x, i10)
c	PC Version
C 1000 format(a10, 1x, i10, 1x, i10)
C 1001 format(a11, 1x, i10, 1x, i10)

      stop

      contains

c     Routine called from main

      subroutine inmsim
      implicit none
      logical more
      character(80) single_line
      integer print_flag

c     Determine if multiple simulations are being performed

      inquire(file = 'fehmn.msim', exist = file_exists)

c     If there are multiple simulations, determine how many

      if(file_exists) then
         open(1,file='fehmn.msim')
         read(1,*) nsim
c       PC Version
C         open(2,file='fehmn.pre.bat')
C         open(3,file='fehmn.post.bat')
c       UNIX version
         open(2,file='fehmn.pre')
         open(3,file='fehmn.post')
         more = .true.
         print_flag = 2
         do while(more)
            read(1,'(a80)',end=1000) single_line
            if(single_line(1:4).ne.'post') then
               write(print_flag,'(a80)') single_line
            else
               print_flag = 3
            end if
         end do
 1000    more = .false.
         close(1)
         close(2)
         close(3)
      else
         nsim = 1
      end if

      return
      end subroutine inmsim

      end program main
