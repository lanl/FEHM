       subroutine pebi (x,y,cpb,vol)
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
CD1  To calculate internodal areas using perpendicular bisectors.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/pebi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:36   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Jan 18 10:06:00 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Wed Jan 10 15:20:50 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:57:50   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:10   pvcs
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
CD5   cpb             REAL*8   O    Perpindicular bisector
CD5   vol             REAL*8   O    Internodal area
CD5   x               REAL*8   I    X nodal coordinates
CD5   y               REAL*8   I    Y nodal coordinates
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
CD7   arx             REAL     Bisector calculation temporary variable
CD7   ary             REAL     Bisector calculation temporary variable
CD7   dsq             REAL     Bisector calculation temporary variable
CD7   i               INT      Loop index
CD7   i1              INT      Variable index
CD7   i2              INT      Variable index
CD7   i3              INT      Variable index
CD7   pneg            REAL     Area calculation temporary variable 
CD7   pos             REAL     Area calculation temporary variable
CD7   xc              REAL     Bisector X coordinate
CD7   yc              REAL     Bisector Y coordinate     
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN pebi
CPS 
CPS   FOR each coordinate position
CPS   
CPS       set indices
CPS       
CPS       IF first and second ycoordinate are equal
CPS          compute bisector xcoordinate using first and second y
CPS       ELSE IF first and third ycoordinate are equal
CPS          compute bisector xcoordinate using first and third y
CPS       ELSE
CPS          compute bisector xcoordinate using all y's
CPS       END IF
CPS       
CPS       IF first and second xcoordinate are equal
CPS          compute bisector ycoordinate using first and second x
CPS       ELSE IF first and third xcoordinate are equal
CPS          compute bisector ycoordinate using first and third x
CPS       ELSE
CPS          compute bisector ycoordinate using all x's
CPS       END IF
CPS       
CPS       calculate internodal area using bisectors
CPS       
CPS   END FOR
CPS   
CPS END pebi
CPS
C***********************************************************************

      implicit none
      
      integer i, i1, i2, i3
      real*8 cpb(*), vol(*), x(*), y(*) 
      real*8 arx, ary, dsq, pneg, pos, xc, yc
      
      do i = 0, 2
         i1 = i + 1
         i2 = mod(i+1,3) + 1
         i3 = mod(i+2,3) + 1
         if (y(i1).eq.y(i2)) then
	    xc = 0.5 * (x(i1)+x(i2))
         else if (y(i1).eq.y(i3)) then
	    xc = 0.5 * (x(i1)+x(i3))
         else
	    xc = 0.5 * ( (x(i3)*x(i3)-x(i1)*x(i1))/(y(i3)-y(i1)) -
     &           (x(i2)*x(i2)-x(i1)*x(i1))/(y(i2)-y(i1))+y(i3)-y(i2))/
     &           ( (x(i3)-x(i1))/(y(i3)-y(i1)) -
     &           (x(i2)-x(i1))/(y(i2)-y(i1)) )
         end if
         if (x(i1).eq.x(i2)) then
	    yc = 0.5 * (y(i1)+y(i2))
         else if (x(i1).eq.x(i3)) then
	    yc = 0.5 * (y(i1)+y(i3))
         else
	    yc = 0.5 * ( (y(i3)*y(i3)-y(i1)*y(i1))/(x(i3)-x(i1)) -
     &           (y(i2)*y(i2)-y(i1)*y(i1))/(x(i2)-x(i1))+x(i3)-x(i2))/
     &           ( (y(i3)-y(i1))/(x(i3)-x(i1)) -
     &           (y(i2)-y(i1))/(x(i2)-x(i1)) )
         end if
         
         arx = yc - 0.5*(y(i1)+y(i2))
         ary = xc - 0.5*(x(i1)+x(i2))
         dsq = (x(i2)-x(i1))**2 + (y(i2)-y(i1))**2
         cpb(i1) = sqrt( (arx*arx+ary*ary)/dsq )
         
         pos = x(i1)*0.5*(y(i1)+y(i2)) + 0.5*(x(i1)+x(i2))*yc +
     &        xc*0.5*(y(i1)+y(i3)) + 0.5*(x(i1)+x(i3))*y(i1)
         pneg = y(i1)*0.5*(x(i1)+x(i2)) + 0.5*(y(i1)+y(i2))*xc +
     &        yc*0.5*(x(i1)+x(i3)) + 0.5*(y(i1)+y(i3))*x(i1)
         vol(i1) = 0.5 * ( pos - pneg )
      enddo
      
      return
      end
