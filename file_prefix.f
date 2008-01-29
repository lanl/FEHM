      subroutine file_prefix(namfil, length)
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
CD1 Determine prefix of input file name (everything before last dot or
CD1 first blank if no dot).
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/file_prefix.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 16:00:44 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:53:28   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:50   pvcs
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
CD3   namfil          CHAR     I    Filename for determining prefix
CD3   length          INT      O    Length of filename prefix
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 None
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
CD5   i               INT      Loop index
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
CPS BEGIN file_prefix
CPS   
CPS   FOR each character in the name (starting from the end)
CPS       search for the last . in the name
CPS   EXIT IF "." found
CPS   END FOR
CPS   
CPS   [if a '.' wasn't found]
CPS   FOR each character in the name (starting from the end)
CPS       search for the last character in the name
CPS   EXIT IF a character is found
CPS   END FOR
CPS   
CPS   set the prefix length to the postion of the character before "." 
CPS    or the last character
CPS   
CPS END file_prefix
CPS 
C***********************************************************************

      implicit none

      integer length

      character*100 namfil

c search for last "." in name
      do length = 100, 1, -1
         if (namfil(length:length) .eq. '.') go to 10
      end do
      
c if "." not found, search back for last character
      
      do length = 101, 2, -1
         if (namfil(length-1:length-1) .ne. ' ') go to 10
      end do
      
 10   continue

      length = length - 1

      end

