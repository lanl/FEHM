      function tyming( dummy )
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To calculate the elapsed cpu time for the ibm.
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
CD2 $Log:   /pvcs.config/fehm90/src/tyming.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:22   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Fri Feb 02 13:53:24 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   04/25/95 10:17:46   llt
CD2 corrected log history information
CD2 
CD2    Rev 1.0   01/27/95 15:40:56   llt
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
CD3 dummy          real*8   I       Dummy argument used for compatability 
CD3                                   with other versions
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
CD4 None
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
CD4 etime_    N/A        IBM system routine for computing cpu time
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
CD5 tb_type  Structure containing usrtime and systime
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 usrtime      real        Returned values from call to cpu time
CD5                             routine etime_
CD5 systime      real        Returned values from call to cpu time
CD5                             routine etime_
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
CD8 This is a general purpose routine used for timing the code 
CD8 execution which is useful information but does not aid any
CD8 high level code requirement.
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 Not Applicable.  See Special Comments.
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
CPS BEGIN tyming
CPS 
CPS etime_ - determine elapsed cpu time
CPS 
CPS END tyming
CPS 
C**********************************************************************

c     Note: for the Digital FORTRAN version, uncomment the use portlib 
c     line.

c	use portlib
      implicit none

c     Note: for the Digital FORTRAN version, comment the declaration
c     line for etime (already declared in portlib).

      real(4) etime
      real(8) tyming
      real(4) dummy(2)
c      tyming = etime(dummy)
      call cpu_time(etime)
      tyming = etime
      end

