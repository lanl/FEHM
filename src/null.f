      logical function null1(wd)
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1 
CD1 Check for null lines or 0's in lines.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 04-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/null.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:34   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   03/18/94 15:55:40   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:02:08   zvd
CD2  
CD2 
CD2    Rev 1.0   01/20/94 10:26:04   pvcs
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
CD3   wd              CHAR     I   String whose content is to be checked
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
CD5   length          INT      Length of input string
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 Check whether a string has characters or is null or zero's. If the
CD6 string has characters return 'FALSE', otherwise if it is null or zero's
CD6 return 'TRUE'.
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
CPS BEGIN null1
CPS 
CPS   determine the length of the input string
CPS 
CPS   FOR each character in the string
CPS   
CPS       IF the character is not a ".", " ", ",", "0", or tab
CPS          set the value to .FALSE. and return
CPS       ENDIF   
CPS       
CPS   ENDFOR
CPS   
CPS   set the value to .TRUE. and return
CPS       
CPS END null1
CPS
C***********************************************************************

      implicit none

      integer i,length
      character*(*) wd

      length = len(wd)
      do i = 1, length
         if (wd(i:i) .ne. '.' .and. 
     *       wd(i:i) .ne. ' ' .and. 
     *       wd(i:i) .ne. ',' .and. 
     *       wd(i:i) .ne. '0' .and. 
     *       wd(i:i) .ne. char(9)) then
            null1 = .FALSE.
            return
         endif
      end do

      null1 = .TRUE.

      return
      end
