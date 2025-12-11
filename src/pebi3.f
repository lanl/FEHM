      subroutine pebi3 (px,py,pz, cpb,vol)
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
CD1  To calculate internodal volumes using perpendicular bisectors.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/pebi3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:36   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:50   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Jan 18 10:06:16 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Wed Jan 10 15:26:58 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:57:52   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:12   pvcs
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
CD5   px              REAL*8   I    X nodal coordinates
CD5   py              REAL*8   I    Y nodal coordinates
CD5   pz              REAL*8   I    Z nodal coordinates
CD5   vol             REAL*8   O    Internodal volume
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
CD6   Identifier      Type     Description
CD6   
CD6   lubksub0                 Solves the set of n linear equations A*x = b.
CD6   ludcmp0                  Given an n by n matrix A, with physical 
CD6                              dimension np, this routine replaces it
CD6                              by the LU decompostion of a row wise 
CD6                              permutation of itself.
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
CD7   amat            REAL*8   Array
CD7   b               REAL*8   Vector, right-hand side
CD7   cx              REAL*8
CD7   cy              REAL*8
CD7   cz              REAL*8
CD7   d               REAL*8   Flag denoting even or odd number of
CD7                              row interchanges
CD7   i               INT      Loop index
CD7   i1              INT      Array index
CD7   i2              INT      Array index
CD7   i3              INT      Array index
CD7   i4              INT      Array index
CD7   iar             INT
CD7   indx            INT      Permutation vector
CD7   iperm           INT
CD7   j               INT      Loop index
CD7   k               INT      Array index
CD7   l               INT      Array index
CD7   n3              INT      Array dimension
CD7   qx              REAL*8
CD7   qy              REAL*8
CD7   qz              REAL*8
CD7   vol1            REAL*8
CD7   vol2            REAL*8
CD7   vol3            REAL*8
CD7   vol4            REAL*8
CD7   vol5            REAL*8
CD7   x1              REAL*8
CD7   x2              REAL*8
CD7   x3              REAL*8
CD7   x4              REAL*8
CD7   xm              REAL*8
CD7   y1              REAL*8
CD7   y2              REAL*8
CD7   y3              REAL*8
CD7   y4              REAL*8
CD7   ym              REAL*8
CD7   z1              REAL*8
CD7   z2              REAL*8
CD7   z3              REAL*8
CD7   z4              REAL*8
CD7   zm              REAL*8
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7   
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN pebi3
CPS
CPS   FOR each coordinate position
CPS       FOR each incremental coordinate position
CPS           calculate coordinate midpoints
CPS       END FOR
CPS   END FOR
CPS
CPS   FOR each coordinate position
CPS       set indices and calculate internodal lengths, etc.
CPS       call ludcmp0 and lubksb0 to solve for bisector coordinates
CPS   END FOR
CPS   
CPS   calculate internodal lengths, etc. from first coordinate to each . . .
CPS   . . . other coordinate
CPS   call ludcmp0 and lubksb0 to solve for coordinates
CPS   
CPS   FOR each coordinate position
CPS       FOR each incremental coordinate position
CPS           IF x internode length is greater than y and z internode length
CPS              use y and z bisector coordinates to calculate . . .
CPS              . . . perpindicular bisector     
CPS           ELSE IF y internode length is greater than z internode length
CPS              use x and z bisector coordinates to calculate . . .
CPS              . . . perpindicular bisector     
CPS           ELSE
CPS              use x and y bisector coordinates to calculate . . .
CPS              . . . perpindicular bisector     
CPS           END IF
CPS           set perpindicular bisector
CPS       END FOR
CPS   END FOR
CPS
CPS   FOR each coordinate position
CPS       set indices
CPS       calculate internodal volumes
CPS   END FOR
CPS   
CPS END pebi3
CPS
C***********************************************************************

      implicit none

      integer n3, i, i1, i2, i3, i4, j, k, l
      parameter (n3 = 3)
      integer indx(n3), iperm(4, 4), iar(4, 4)
      real*8 cpb(4, 4), vol(4), px(4), py(4), pz(4) 
      real*8 cx(4), cy(4), cz(4), xm(4, 4), ym(4, 4), zm(4, 4)
      real*8 amat(n3, n3), b(n3), d
      real*8 arx, ary, arz, qx, qy, qz
      real*8 vol1, vol2, vol3, vol4, vol5
      real*8 x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
      
      data iperm/1,2,3,4, 2,1,4,3, 3,1,2,4, 4,1,3,2/
      data iar /0,4,2,3, 3,0,4,1, 4,1,0,2, 2,3,1,0/
      
      do i = 1, 3
         do j = i+1, 4
	    xm(i,j) = 0.5 * (px(i)+px(j))
	    ym(i,j) = 0.5 * (py(i)+py(j))
	    zm(i,j) = 0.5 * (pz(i)+pz(j))
	    xm(j,i) = xm(i,j)
	    ym(j,i) = ym(i,j)
	    zm(j,i) = zm(i,j)
         enddo
      enddo
      
      do i = 1, 4
         i1 = iperm(1,i)
         i2 = iperm(2,i)
         i3 = iperm(3,i)
         i4 = iperm(4,i)
         
         amat(1,1) = px(i2)-px(i3)
         amat(1,2) = py(i2)-py(i3)
         amat(1,3) = pz(i2)-pz(i3)
         amat(2,1) = px(i2)-px(i4)
         amat(2,2) = py(i2)-py(i4)
         amat(2,3) = pz(i2)-pz(i4)
         amat(3,1) = (py(i2)-py(i3))*(pz(i2)-pz(i4)) -
     &        (pz(i2)-pz(i3))*(py(i2)-py(i4))
         amat(3,2) = (pz(i2)-pz(i3))*(px(i2)-px(i4)) -
     &        (px(i2)-px(i3))*(pz(i2)-pz(i4))
         amat(3,3) = (px(i2)-px(i3))*(py(i2)-py(i4)) -
     &        (py(i2)-py(i3))*(px(i2)-px(i4))

         b(1) = xm(i2,i3)*amat(1,1)+ym(i2,i3)*amat(1,2)+zm(i2,i3)*
     &        amat(1,3)
         b(2) = xm(i2,i4)*amat(2,1)+ym(i2,i4)*amat(2,2)+zm(i2,i4)*
     &        amat(2,3)
         b(3) = px(i2)*amat(3,1)+py(i2)*amat(3,2)+pz(i2)*amat(3,3)

         call ludcmp0 (amat,n3,n3,indx,d)
         call lubksb0 (amat,n3,n3,indx,b)
         cx(i1) = b(1)
         cy(i1) = b(2)
         cz(i1) = b(3)
      enddo

      amat(1,1) = px(1)-px(2)
      amat(1,2) = py(1)-py(2)
      amat(1,3) = pz(1)-pz(2)
      amat(2,1) = px(1)-px(3)
      amat(2,2) = py(1)-py(3)
      amat(2,3) = pz(1)-pz(3)
      amat(3,1) = px(1)-px(4)
      amat(3,2) = py(1)-py(4)
      amat(3,3) = pz(1)-pz(4)

      b(1) = xm(1,2)*amat(1,1) + ym(1,2)*amat(1,2) + zm(1,2)*amat(1,3)
      b(2) = xm(1,3)*amat(2,1) + ym(1,3)*amat(2,2) + zm(1,3)*amat(2,3)
      b(3) = xm(1,4)*amat(3,1) + ym(1,4)*amat(3,2) + zm(1,4)*amat(3,3)

      call ludcmp0 (amat,n3,n3,indx,d)
      call lubksb0 (amat,n3,n3,indx,b)
      qx = b(1)
      qy = b(2)
      qz = b(3)

      do i = 1, 3
         do j = i+1, 4
	    k = iar(i,j)
	    l = iar(j,i)
	    if (abs(px(i)-px(j)) .ge. abs(py(i)-py(j)) .and.
     &           abs(px(i)-px(j)) .ge. abs(pz(i)-pz(j))) then
               arx = (qy-ym(i,j))*(cz(k)-cz(l)) -
     &              (qz-zm(i,j))*(cy(k)-cy(l))
               cpb(i,j) = 0.5 * arx / (px(i)-px(j))
	    else if (abs(py(i)-py(j)) .ge. abs(pz(i)-pz(j))) then
               ary = (qz-zm(i,j))*(cx(k)-cx(l)) -
     &              (qx-xm(i,j))*(cz(k)-cz(l))
               cpb(i,j) = 0.5 * ary / (py(i)-py(j))
	    else
               arz = (qx-xm(i,j))*(cy(k)-cy(l)) -
     &              (qy-ym(i,j))*(cx(k)-cx(l))
               cpb(i,j) = 0.5 * arz / (pz(i)-pz(j))
	    end if
	    cpb(j,i) = cpb(i,j)
         enddo
      enddo

      do i = 1, 4
         i1 = iperm(1,i)
         i2 = iperm(2,i)
         i3 = iperm(3,i)
         i4 = iperm(4,i)
         x1 = qx
         y1 = qy
         z1 = qz
         x2 = cx(i4)
         y2 = cy(i4)
         z2 = cz(i4)
         x3 = cx(i3)
         y3 = cy(i3)
         z3 = cz(i3)
         x4 = cx(i2)
         y4 = cy(i2)
         z4 = cz(i2)
         vol1 = (y1*z2-y2*z1)*(x3-x4) + (y1*z3-y3*z1)*(x4-x2) +
     &        (y1*z4-y4*z1)*(x2-x3) + (y2*z3-y3*z2)*(x1-x4) +
     &        (y2*z4-y4*z2)*(x3-x1) + (y3*z4-y4*z3)*(x1-x2)

         x1 = px(i1)
         y1 = py(i1)
         z1 = pz(i1)
         x3 = cx(i2)
         y3 = cy(i2)
         z3 = cz(i2)
         x4 = cx(i3)
         y4 = cy(i3)
         z4 = cz(i3)
         vol2 = (y1*z2-y2*z1)*(x3-x4) + (y1*z3-y3*z1)*(x4-x2) +
     &        (y1*z4-y4*z1)*(x2-x3) + (y2*z3-y3*z2)*(x1-x4) +
     &        (y2*z4-y4*z2)*(x3-x1) + (y3*z4-y4*z3)*(x1-x2)

         x2 = xm(i1,i3)
         y2 = ym(i1,i3)
         z2 = zm(i1,i3)
         x3 = cx(i2)
         y3 = cy(i2)
         z3 = cz(i2)
         x4 = cx(i4)
         y4 = cy(i4)
         z4 = cz(i4)
         vol3 = (y1*z2-y2*z1)*(x3-x4) + (y1*z3-y3*z1)*(x4-x2) +
     &        (y1*z4-y4*z1)*(x2-x3) + (y2*z3-y3*z2)*(x1-x4) +
     &        (y2*z4-y4*z2)*(x3-x1) + (y3*z4-y4*z3)*(x1-x2)

         x2 = xm(i1,i4)
         y2 = ym(i1,i4)
         z2 = zm(i1,i4)
         x3 = cx(i3)
         y3 = cy(i3)
         z3 = cz(i3)
         x4 = cx(i2)
         y4 = cy(i2)
         z4 = cz(i2)
         vol4 = (y1*z2-y2*z1)*(x3-x4) + (y1*z3-y3*z1)*(x4-x2) +
     &        (y1*z4-y4*z1)*(x2-x3) + (y2*z3-y3*z2)*(x1-x4) +
     &        (y2*z4-y4*z2)*(x3-x1) + (y3*z4-y4*z3)*(x1-x2)

         x2 = xm(i1,i2)
         y2 = ym(i1,i2)
         z2 = zm(i1,i2)
         x3 = cx(i4)
         y3 = cy(i4)
         z3 = cz(i4)
         x4 = cx(i3)
         y4 = cy(i3)
         z4 = cz(i3)
         vol5 = (y1*z2-y2*z1)*(x3-x4) + (y1*z3-y3*z1)*(x4-x2) +
     &        (y1*z4-y4*z1)*(x2-x3) + (y2*z3-y3*z2)*(x1-x4) +
     &        (y2*z4-y4*z2)*(x3-x1) + (y3*z4-y4*z3)*(x1-x2)

         vol(i1) = (vol1+vol2+vol3+vol4+vol5) / 6.
      enddo

      return
      end
