      subroutine writeio(unit)
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
CD1 Write assigned file names, unit numbers, and file purpose to 
CD1 specified output unit.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 15-Nov-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/writeio.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Fri Feb 02 14:30:38 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:56:04   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:32   pvcs
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
CD3   nmfilc          CHAR     I    Terminal usage on/off
CD3   nmfild          CHAR     I    Terminal usage on/off
CD3   unit            INT      I    Output unit number
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
CD4   cuser           CHAR     faax   String written when not using a
CD4                                     particular input/output file
CD4   iatty           INT      faai   Unit number for all tty output
CD4   incoor          INT      faai   Unit number for coordinate input file
CD4   inpt            INT      faai   Unit number for input file
CD4   inzone          INT      faai   Unit number for zone input file
CD4   iocntl          INT      faai   Unit number for control file
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   iread           INT      faai   Unit number for restart file (to read)
CD4   isave           INT      faai   Unit number for restart file (to write)
CD4   ischk           INT      faai   Unit number for input data check file
CD4   iscon           INT      faai   Unit number for contour data file
CD4   ishis           INT      faai   Unit number for history data file
CD4   isstor          INT      faai   Unit number for element coefficient file
CD4   istrc           INT      faai   Unit number for tracer history data file
CD4   nmfil1          CHAR     faax   Main input file name
CD4   nmfil2          CHAR     faax   Coordinate input file name
CD4   nmfil3          CHAR     faax   Zone input file name
CD4   nmfil4          CHAR     faax   Main output file name
CD4   nmfil5          CHAR     faax   Restart input file name
CD4   nmfil6          CHAR     faax   Restart output file name
CD4   nmfil7          CHAR     faax   Simulation history output file name
CD4   nmfil8          CHAR     faax   Contour plot output file name
CD4   nmfil9          CHAR     faax   Solute output file name
CD4   nmfila          CHAR     faax   Coefficient storage output file name
CD4   nmfilb          CHAR     faax   Input check output file name
CD4   nmfils          CHAR     faax   Control file name
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
CPS BEGIN writeio
CPS 
CPS   write assigned files, unit numbers, and purpose to output unit
CPS   
CPS END writeio
CPS
C***********************************************************************

      use comai
      use comxi
      use comco2, only : icarb
      implicit none

      integer unit
c gaz 111621 added error file (unit 23 referenced in countour output)
      write(unit, 6100)  iocntl, nmfil(1),
     *     inpt  , nmfil(2),
     *     incoor, nmfil(3),
     *     inzone, nmfil(4),
     *     iout  , nmfil(5),
     *     iread , nmfil(6),
     *     isave , nmfil(7),
     *     ishis , nmfil(8),
     *     istrc , nmfil(9),
     *     iscon , nmfil(10),
     *     iscon1, nmfil(11),
     *     isstor, nmfil(12),
     *     ischk , nmfil(13),
     *     ierr  , nmfil(14),      
     *     cuser 

 6100 format(/,' File purpose - Variable - Unit number - File name',
     *     //,4x,'control       - iocntl -', i3, ' - ', a100,
     *     /,4x, 'input         - inpt   -', i3, ' - ', a100, 
     *     /,4x, 'geometry      - incoor -', i3, ' - ', a100, 
     *     /,4x, 'zone          - inzone -', i3, ' - ', a100, 
     *     /,4x, 'output        - iout   -', i3, ' - ', a100, 
     *     /,4x, 'initial state - iread  -', i3, ' - ', a100, 
     *     /,4x, 'final state   - isave  -', i3, ' - ', a100, 
     *     /,4x, 'time history  - ishis  -', i3, ' - ', a100, 
     *     /,4x, 'time his.(tr) - istrc  -', i3, ' - ', a100, 
     *     /,4x, 'contour plot  - iscon  -', i3, ' - ', a100, 
     *     /,4x, 'con plot (dp) - iscon1 -', i3, ' - ', a100, 
     *     /,4x, 'fe coef stor  - isstor -', i3, ' - ', a100, 
     *     /,4x, 'input check   - ischk  -', i3, ' - ', a100, 
     *     /,4x, 'error file    - ierr   -', i3, ' - ', a100,      
     *     /,1x, 'Value provided to subroutine user: ', a9, /)
  
      end