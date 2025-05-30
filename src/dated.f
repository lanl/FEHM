      subroutine  dated   ( idatex,jtimex )
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To determine the current date and time for the sun.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                     However, previous non-YMP
CD2                                     versions of FEHM exist, and
CD2                                     the current version differs
CD2                                     from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dated.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:06   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Tue May 07 15:59:40 1996   llt
CD2 added p for prototype to version label
CD2 
CD2    Rev 1.10   Tue May 07 12:58:44 1996   llt
CD2 new version number
CD2 
CD2    Rev 1.9   Mon Jan 29 14:25:50 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   09/13/95 09:33:34   llt
CD2 put "p" in version number - FEHMN 95-09-09p - for prototype
CD2 
CD2    Rev 1.7   09/11/95 09:11:14   llt
CD2 updated version number
CD2 
CD2    Rev 1.6   05/01/95 09:55:16   llt
CD2 protype version
CD2 
CD2    Rev 1.5   05/01/95 09:21:28   llt
CD2 new version
CD2 
CD2    Rev 1.4   04/25/95 09:26:38   llt
CD2 corrected log history information
CD2 
CD2    Rev 1.3   02/01/95 19:02:20   llt
CD2 added prototype to version number
CD2 
CD2    Rev 1.2   02/01/95 14:57:56   llt
CD2 new verno
CD2 
CD2    Rev 1.1   01/28/95 15:22:08   llt
CD2 updated to new version
CD2 
CD2    Rev 1.0   01/27/95 15:40:26   llt
CD2 remade archive to correct pvcs problem
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier     Type    Use      Description
CD3 
CD3 idatex         char     O       String containing the date
CD3                                    information
CD3 jtimex         char     O       String containing the time
CD3                                    information
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 None
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 verno    FEHMN Version number
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
CD4 idate     N/A        SUN system routine for determining the date
CD4 itime     N/A        SUN system routine for determining the time
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 jjdate       int         Returned date information in integer form
CD5 jjtime       int         Returned time information in integer form
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This routine is used for timing information only.
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 Not Applicable.  See special Comments. 
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN dated
CPS 
CPS idate - determine current date
CPS Perform internal write of date information into string
CPS itime - determine current time
CPS Perform internal write of time information into string
CPS set version number
CPS 
CPS END dated
CPS 
C**********************************************************************
      use comai
      implicit none

      character*8 current_date
      character*10 current_time
      character*11  idatex
      character*8   jtimex

c      integer*4  jjdate (3), jjtime (3)
c      idatex='        '
c      jtimex='        '

c      jjdate(1)=0
c      jjdate(2)=0
c      jjdate(3)=0
     
c      jjtime(1)=0
c      jjtime(2)=0
c      jjtime(3)=0

c      call  idate( jjdate )
c      jjdate(3) =  mod( jjdate(3),100 )
c      write(idatex,'(i2.2,1h/,i2.2,1h/,i2.2)')  jjdate(3) , jjdate(2)
c     *                                                   , jjdate(1)
c      call  itime( jjtime )
c      write(jtimex,'(i2.2,1h:,i2.2,1h:,i2.2)')  jjtime
      call date_and_time(current_date,current_time)
      jtimex(1:2) = current_time(1:2)
      jtimex(3:3) = ":"
      jtimex(4:5) = current_time(3:4)
      jtimex(6:6) = ":"
      jtimex(7:8) = current_time(5:6)
      idatex(1:2) = current_date(5:6)
      idatex(3:3) = "/"
      idatex(4:5) = current_date(7:8)
      idatex(6:6) = "/"
      idatex(7:8) = current_date(1:2)
      idatex(9:11) = current_date(3:4)
! Version number passed to GoldSim
!     vernum = 3.4 before gaz update
      vernum = 3.6
! Code version identifier
!     verno = "FEHM V3.4.2lbUbuntu16 25-02-20 QA:NA"
!     verno = "FEHM V3.6.2   DATE QA:REL"
      verno = "FEHM V3.6.3   DATE QA:DEV"

      end
