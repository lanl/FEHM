      subroutine infiles(simnum)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
CD1 Control reading of input data files.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 04-JAN-94    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/infiles.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Wed Jun 12 15:07:58 1996   zvd
CD2 Modified order of write to iout to compensate for writes from 
CD2 optional input file routine
CD2 
CD2    Rev 1.8   Tue Jan 30 11:41:28 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.7   08/16/95 11:42:46   zvd
CD2 Updated pseudocode to reflect change made by SEH.
CD2 
CD2    Rev 1.6   08/16/95 10:56:44   zvd
CD2  Added format for writing title to output files.
CD2 
CD2    Rev 1.5   08/08/95 08:56:38   awolf
CD2 SEH addition - rarng now called here to properly start xz problems
CD2 and to allow using zones for xz problems
CD2 
CD2    Rev 1.4   05/02/95 10:46:46   llt
CD2 title was getting written to wrong contour output file
CD2 
CD2    Rev 1.3   03/18/94 16:03:26   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.2   02/28/94 11:55:58   zvd
CD2 Corrected problem of writes to coefficient storage file when it should be
CD2 read only file.
CD2 
CD2    Rev 1.1   02/14/94 11:57:48   zvd
CD2 Added write for title in dual/dpdp contour data file.
CD2 
CD2    Rev 1.0   01/20/94 10:24:56   pvcs
CD2 original version in process of being certified
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
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
CD4   incoor          INT      faai   Unit number for coordinate input file
CD4   inpt            INT      faai   Unit number for input file
CD4   inzone          INT      faai   Unit number for zone input file
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   ischk           INT      faai   Unit number for input data check file
CD4   iscon           INT      faai   Unit number for contour data file
CD4   iscon1          INT      faai   Unit number for dual porosity or dpdp 
CD4                                     contour data file
CD4   ishis           INT      faai   Unit number for history data file
CD4   istrc           INT      faai   Unit number for tracer history data file
CD4   wdd             CHAR     faac   Character input string
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   zone                     Set zone information
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
CD5   cnum            INT      Number of times zone has been called
CD5   macro           CHAR     Control statement keyword
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
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN infiles
CPS 
CPS   rewind input file and read input title 
CPS   write title to all output files being used and tty if being used
CPS   call writeio to write files being used to general output file
CPS   
CPS   rewind the coordinate input file and call incoord to read . . .
CPS   . . . coordinate and element data
CPS   
CPS   IF this is a 2D problem
CPS      call rarng to rearrange coefficients if necessary
CPS   END IF
CPS
CPS   IF a zone input file exists
CPS      rewind the zone input file
CPS      REPEAT
CPS        read a line from the input file
CPS
CPS        IF the zone macro is read
CPS           write macro identifier and input file unit to output file . . .
CPS           . . . and tty if enabled
CPS           call zone to set zone information
CPS        ELSE IF stop macro is read
CPS           write macro identifier and input file unit to output file and . . .
CPS        . . . tty if enabled
CPS           EXIT repeat loop
CPS        END IF
CPS 
CPS      UNTIL stop is found
CPS   END IF
CPS      
CPS   rewind the input file   
CPS   reread the input title
CPS   call input to read data from the input file
CPS
CPS END infiles
CPS
C***********************************************************************

      use comai
      implicit none

      real*8 simnum
      integer cnum
c gaz 062823
      integer nzone_saved, icall_sv
      character*4 macro

      cnum = 0

      rewind inpt
      read (inpt, '(a80)') wdd

      if (iout .ne. 0) write(iout, 5005) wdd
      if (iptty  .gt. 0) write(iptty, 5005) wdd

      if (ishis  .gt. 0) write(ishis, '(a80)') wdd
      if (iscon  .gt. 0) write(iscon, '(a80)') wdd
      if (iscon1 .gt. 0) write(iscon1, '(a80)') wdd
      if (istrc  .gt. 0) write(istrc, '(a80)') wdd
      if (ischk  .gt. 0) write(ischk, 5005) wdd
 5005 format (/, 1x, a80, //)

      if (iout .ne. 0) call writeio (iout)

c****  coordinate input section ('incoor'-file) ****

      rewind incoor
      if(irun.eq.1) call incoord

c:: added by seh 8/95
c**** complete element information ****
      if(irun.eq.1) call rarng(0)

c****  zone input section ('inzone'-file) ****

      if (inzone .ne. inpt) then

         rewind inzone

 100     continue
         
         read (inzone, '(a4)') macro

c**** set zone information ****

         if (macro .eq. 'zone') then

            if (iout .ne. 0) write(iout, 6010) macro,  'inzone', inzone
            if (iptty .gt. 0) write(iptty, 6010) macro,  'inzone', 
     *           inzone
 6010       format(1x, '**** input title : ', a4, ' **** ', a6, ' = ',
     *           i3, ' ****')
c gaz 062823 check if zone file contains save zones
            nzone_saved = 0
            icall_sv = 0
            read(inzone,'(a30)') wdd(1:30)
            rewind inzone
c gaz 110123 added    'infi'
            call zone_saved(1,'infi',nzone_saved,icall_sv,inzone)
            rewind inzone
            cnum = cnum + 1
            call zone(cnum, inzone)
            goto 105

         else if (macro .eq.'stop'.or.macro .eq.'   ') then

            if (iout .ne. 0) write(iout, 6010) macro, 'inzone', inzone
            if (iptty .gt. 0) write(iptty, 6010) macro,  'inzone', 
     *           inzone

            go  to  110

         end if

 105     continue
         go  to  100
         
 110     continue

      end if

c**** input section ('inpt  '-file) ****

      rewind inpt
      call input(cnum, simnum)

      end
