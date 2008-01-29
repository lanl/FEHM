      subroutine determ(det,dm,iorder)
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
CD1  PURPOSE
CD1
CD1  Evaluate determinates.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/determ.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Jan 18 09:08:32 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.3   Wed Jan 10 10:59:40 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.2   Tue Jan 09 15:42:26 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:57:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:58   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4  
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5
CD5   det             REAL     O    Value of determinate
CD5   dm              REAL     I    Input array node coordinates
CD5   iorder          INT      I    Equation order
CD5
CD5 Interface Tables
CD5
CD5   None 
CD5
CD5 Files
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6   None
CD6
CD6 Global Subprograms
CD6
CD6   None
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   det1            REAL     Intermediate subdeterminate value
CD7   det2            REAL     Intermediate subdeterminate value
CD7   det3            REAL     Intermediate subdeterminate value
CD7   det4            REAL     Intermediate subdeterminate value
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN determ
CPS 
CPS   IF the order is four
CPS      evaluate subdeterminates
CPS      evaluate determinate
CPS   END IF
CPS
CPS END determ
CPS
C***********************************************************************

      implicit none
      
      integer iorder
      real*8 det1,det2,det3,det4,det,dm(iorder,iorder)
      
      if(iorder.eq.4) then
c     evaluate subdeterminates
         det1=dm(2,2)*dm(3,3)*dm(4,4)
     &        +dm(3,2)*dm(4,3)*dm(2,4)
     &        +dm(4,2)*dm(2,3)*dm(3,4)
     &        -dm(2,4)*dm(3,3)*dm(4,2)
     &        -dm(2,3)*dm(3,2)*dm(4,4)
     &        -dm(2,2)*dm(3,4)*dm(4,3)
         det2=dm(1,2)*dm(3,3)*dm(4,4)
     &        +dm(3,2)*dm(4,3)*dm(1,4)
     &        +dm(4,2)*dm(1,3)*dm(3,4)
     &        -dm(1,4)*dm(3,3)*dm(4,2)
     &        -dm(1,3)*dm(3,2)*dm(4,4)
     &        -dm(1,2)*dm(3,4)*dm(4,3)
         det3=dm(1,2)*dm(2,3)*dm(4,4)
     &        +dm(2,2)*dm(4,3)*dm(1,4)
     &        +dm(4,2)*dm(1,3)*dm(2,4)
     &        -dm(1,4)*dm(2,3)*dm(4,2)
     &        -dm(1,3)*dm(2,2)*dm(4,4)
     &        -dm(1,2)*dm(2,4)*dm(4,3)
         det4=dm(1,2)*dm(2,3)*dm(3,4)
     &        +dm(2,2)*dm(3,3)*dm(1,4)
     &        +dm(3,2)*dm(1,3)*dm(2,4)
     &        -dm(1,4)*dm(2,3)*dm(3,2)
     &        -dm(1,3)*dm(2,2)*dm(3,4)
     &        -dm(1,2)*dm(2,4)*dm(3,3)
         det=dm(1,1)*det1-dm(2,1)*det2+dm(3,1)*det3-dm(4,1)*det4
      endif
      
      return
      end
      
