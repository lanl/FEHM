      subroutine shap3r(nga)
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
CD1  Evaluate 3-D finite element shape functions at quadrature points.
CD1  Calculates 8-node quad.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/shap3r.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Jan 18 10:50:16 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Thu Jan 11 11:13:34 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:58:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:46   pvcs
CD2 original version in process of being certified
CD2 
************************************************************************
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
CD5   nga             INT      I    Integration point
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
CD6                            COMMON
CD6   Identifier      Type     Block  Description
CD6
CD6   dr              REAL*8   fbs    Contains weights for integration
CD6                                     points (bricks, rectangles)
CD6   eta             REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   exci            REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   intg            INT      faai   Indicates integration type used
CD6   si              REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   xd              REAL*8   fbs    Global coordinates of the nodes in a
CD6                                     finite element
CD6   yd              REAL*8   fbs    Global coordinates of the nodes in a
CD6                                     finite element
CD6   wr              REAL*8   fbs    Finite element shape functions
CD6                                     (rectangles)
CD6   wxr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to x (rectangles)
CD6   wyr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y (rectangles)
CD6   wzr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to z (rectangles)
CD6   zd              REAL*8   fbs    Global coordinates of the nodes in a
CD6                                     finite element
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
CD7   i               INT      Loop index
CD7   nsl             INT      Number of nodes in a quad element
CD7   tintg           REAL*8
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN shap3r
CPS
CPS   set nodes in element to 8
CPS
CPS   IF Gauss quadrature is used
CPS      set integration constant to 1/(3^.5)
CPS   ELSE
CPS      set integration constant to 1
CPS   END IF
CPS
CPS   set local X, Y, Z coordinates
CPS   calculate local coordinates for integration point
CPS   FOR each node of element
CPS       calculate shape function and derivatives for integration point
CPS   END FOR
CPS
CPS END shap3r
CPS
C***********************************************************************


      use combi
      use comdti
      use comai
      implicit none


      integer i,nga,nsl
      real*8 tintg


      nsl=8
c define local coordinates,integration types
      if ( intg.eq.1 ) then
         tintg=1.0/sqrt(3.d0)
      else
         tintg=1.0
      endif
c local coordinates
      xd(1)=-1.0
      xd(2)=1.0
      xd(3)=1.0
      xd(4)=-1.0
      xd(5)=-1.0
      xd(6)=1.0
      xd(7)=1.0
      xd(8)=-1.0
      yd(1)=-1.0
      yd(2)=-1.0
      yd(3)=1.0
      yd(4)=1.0
      yd(5)=-1.0
      yd(6)=-1.0
      yd(7)=1.0
      yd(8)=1.0
      zd(1)=-1.0
      zd(2)=-1.0
      zd(3)=-1.0
      zd(4)=-1.0
      zd(5)= 1.0
      zd(6)= 1.0
      zd(7)= 1.0
      zd(8)= 1.0
c define integration points
      si(nga)=xd(nga)*tintg
      eta(nga)=yd(nga)*tintg
      exci(nga)=zd(nga)*tintg
      do 20 i=1,nsl
      wr(nga,i)=0.125*(1.+xd(i)*si(nga))*(1.+yd(i)*eta(nga))
     **(1.+zd(i)*exci(nga))
      wxr(nga,i)=0.125*xd(i)*(1.+yd(i)*eta(nga))*(1.+zd(i)*exci(nga))
      wyr(nga,i)=0.125*yd(i)*(1.+xd(i)*si(nga))*(1.+zd(i)*exci(nga))
      wzr(nga,i)=0.125*zd(i)*(1.+xd(i)*si(nga))*(1.+yd(i)*eta(nga))
      dr(nga)=1.
   20 continue


      return 
      end