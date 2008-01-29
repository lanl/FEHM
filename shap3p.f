      subroutine shap3p(nga)
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
CD1  Evaluate 3-D prism element shape functions at quadrature points.
CD1  Calculates 6-node prism.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/shap3p.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Jan 18 10:50:24 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Thu Jan 11 11:09:44 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:58:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:44   pvcs
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
CD6   dp              REAL*8   fbs    Contains weights for integration points
CD6                                     (prisms, triangles)
CD6   eta             REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   exci            REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   si              REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   wp              REAL*8   fbs    Finite element shape functions
CD6                                     (prisms)
CD6   wxp             REAL*8   fbs    Derivative of shape functions with
CD6                                      respect to x (prisms)
CD6   wyp             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y (prisms)
CD6   wzp             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to z (prisms)
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
CD7   None
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN shap3p
CPS 
CPS   set local coordinates for integration point
CPS   calculate shape function and derivatives for integration point of
CPS    each node
CPS   
CPS END shap3p
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      integer nga

      si(1)=0.5
      si(2)=0.5
      si(3)=0.0
      si(4)=0.5
      si(5)=0.5
      si(6)=0.0
      eta(1)=0.0
      eta(2)=0.5
      eta(3)=0.5
      eta(4)=0.0
      eta(5)=0.5
      eta(6)=0.5
      exci(1)=1.0
      exci(2)=1.0
      exci(3)=1.0
      exci(4)=-1.0
      exci(5)=-1.0
      exci(6)=-1.0
      wp(nga,1)=(1.0-si(nga)-eta(nga))*(1.0+exci(nga))/2.0
      wxp(nga,1)=-(1.0+exci(nga))/2.0
      wyp(nga,1)=-(1.0+exci(nga))/2.0
      wzp(nga,1)=(1.0-si(nga)-eta(nga))/2.0
      wp(nga,2)=si(nga)*(1.0+exci(nga))/2.0
      wxp(nga,2)=(1.0+exci(nga))/2.0
      wyp(nga,2)=0.0
      wzp(nga,2)=si(nga)/2.0
      wp(nga,3)=eta(nga)*(1.0+exci(nga))/2.0
      wxp(nga,3)=0.0
      wyp(nga,3)=(1.0+exci(nga))/2.0
      wzp(nga,3)=eta(nga)/2.0
      wp(nga,4)=(1.0-si(nga)-eta(nga))*(1.0-exci(nga))/2.0
      wxp(nga,4)=-(1.0-exci(nga))/2.0
      wyp(nga,4)=-(1.0-exci(nga))/2.0
      wzp(nga,4)=-(1.0-si(nga)-eta(nga))/2.0
      wp(nga,5)=si(nga)*(1.0-exci(nga))/2.0
      wxp(nga,5)=(1.0-exci(nga))/2.0
      wyp(nga,5)=0.0
      wzp(nga,5)=-si(nga)/2.0
      wp(nga,6)=eta(nga)*(1.0-exci(nga))/2.0
      wxp(nga,6)=0.0
      wyp(nga,6)=(1.0-exci(nga))/2.0
      wzp(nga,6)=-eta(nga)/2.0
      dp(nga)=1./6.

      return
      end
