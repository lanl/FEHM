      subroutine setunits(numfil)
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
CD1 Set file unit numbers.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 25-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/setunits.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Thu Feb 01 16:08:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:55:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:32   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 None
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
CD4   isave           INT      faai   Unit number for restart file (to write)
CD4   ischk           INT      faai   Unit number for input data check file
CD4   iscon           INT      faai   Unit number for contour data file
CD4   iscon1          INT      faai   Unit number for dual porosity or dpdp
CD4                                     contour data file
CD4   ishis           INT      faai   Unit number for history data file
CD4   isstor          INT      faai   Unit number for element coefficient file
CD4   istrc           INT      faai   Unit number for tracer history data file
CD4   iowork          CHAR     faay   File use/type designator
CD4   isave           INT      faai   Unit number for restart file (to write)
CD4
CD4 Global Subprograms
CD4
CD4   None
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
CD5   numfil          INT ARY  I/O file unit numbers
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
CPS BEGIN setunits
CPS 
CPS   Assign file unit number for each I/O file [iocntl, inpt, incoor, . . .
CPS   . . . inzone, iout, iread, isave, ishis, istrc, iscon, iscon1, . . .
CPS   . . . isstor, ischk]
CPS   
CPS END setunits
CPS
C***********************************************************************

      use comai
      implicit none

      integer numfil(*)

C **** Set file unit numbers ****
      iocntl = numfil(1)
      inpt = numfil(2)
      incoor = numfil(3)
      inzone = numfil(4)
      iout = numfil(5)
      iread = numfil(6)
      isave = numfil(7)
      ishis = numfil(8)
      istrc = numfil(9)
      iscon = numfil(10)
      iscon1 = numfil(11)
      isstor = numfil(12)
      ischk = numfil(13)
      
      end
