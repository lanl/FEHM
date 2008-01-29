      subroutine close_files
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
CD1 Close all open files
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 16-NOV-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/close_files.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   03/18/94 15:53:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:38   pvcs
CD2 original version in process of being certified
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
CD3   See unit descriptors under Global Variables
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
CD4   iread           INT      faai   Unit number for restart file (to read)
CD4   isave           INT      faai   Unit number for restart file (to write)
CD4   ischk           INT      faai   Unit number for input data check file
CD4   iscon           INT      faai   Unit number for contour data file
CD4   iscon1          INT      faai   Unit number for dual porosity or dpdp
CD4                                     contour data file
CD4   ishis           INT      faai   Unit number for history data file
CD4   isstor          INT      faai   Unit number for element coefficient file
CD4   istrc           INT      faai   Unit number for tracer history data file
CD4
CD4 Global Subprograms
CD4
CD4   None
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 None
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
CD8 None
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 N/A
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
CPS BEGIN close_files
CPS 
CPS   [for each file (inpt, incoor, inzone, iout, iread, isave, ishis, . . .
CPS   . . . istrc, iscon, iscon1, isstor, ischk)] 
CPS   IF file unit number is not 0, 5, or 6 [undefined or terminal]
CPS      close the file
CPS   ENDIF
CPS   
CPS END close_files
CPS
C***********************************************************************
      use comai
      implicit none

      if (inpt .ne. 0 .and. inpt .ne. 5 .and. inpt .ne. 6)  
     *     close (inpt  , status = 'keep')
      if (incoor .ne. 0 .and. incoor .ne. 5 .and. incoor .ne. inpt)
     *     close (incoor, status = 'keep')
      if (inzone .ne. 0 .and. inzone .ne. 5 .and. inzone .ne. inpt)
     *     close (inzone, status = 'keep')
      if (iout .ne. 0 .and. iout .ne. 5 .and. iout .ne. 6)  
     *     close (iout  , status = 'keep')
      if (iread .ne. 0 .and. iread .ne. 5 .and. iread .ne. 6)
     *     close (iread , status = 'keep')
      if (isave .ne. 0 .and. isave .ne. 5 .and. isave .ne. 6)
     *     close (isave , status = 'keep')
      if (ishis .ne. 0 .and. ishis .ne. 5 .and. ishis .ne. 6) 
     *     close (ishis , status = 'keep')
      if (istrc .ne. 0 .and. istrc .ne. 5 .and. istrc .ne. 6)
     *     close (istrc , status = 'keep')
      if (iscon .ne. 0 .and. iscon .ne. 5 .and. iscon .ne. 6)
     *     close (iscon , status = 'keep')
      if (iscon1 .ne. 0 .and. iscon1 .ne. 5 .and. iscon1 .ne. 6)
     *     close (iscon1 , status = 'keep')
      if (isstor .ne. 0 .and. isstor .ne. 5 .and. isstor .ne. 6) 
     *     close (isstor, status = 'keep')
      if (ischk .ne. 0 .and. ischk .ne. 5 .and. ischk .ne. 6) 
     *     close (ischk , status = 'keep')

        end
