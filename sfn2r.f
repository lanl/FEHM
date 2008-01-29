      subroutine sfn2r (sid, etad)
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
CD1 Evaluate shape functions for 2D calculations.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 06-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sfn2r.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Fri Feb 02 10:37:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:55:46   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:40   pvcs
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
CD3   sid             REAL*8   I    Local coordinates in a finite element of
CD3                                   the numerical integration points
CD3   etad            REAL*8   I    Local coordinates in a finite element of
CD3                                   the numerical integration points
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
CD4   w               REAL*8   fbs    Finite element shape functions
CD4   wx              REAL*8   fbs    Derivative of shape functions with
CD4                                     respect to x
CD4   wy              REAL*8   fbs    Derivative of shape functions with
CD4                                     respect to y
CD4   xd              REAL*8   fbs    Global coordinates of the nodes in a
CD4                                     finite element
CD4   yd              REAL*8   fbs    Global coordinates of the nodes in a
CD4                                     finite element
CD4
CD4 Global Subprograms
CD4
CD4   None
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
CD5   nsl             INT      Number of coordinate positions for 2D
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
CD9 2.2 Finite-Element Coefficient Generation
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
CPS BEGIN sfn2r
CPS 
CPS   set number of coordinate postions to 4
CPS   initialize node global coordinates for 2D
CPS   
CPS   FOR each coordinate position
CPS       calculate the shape function and its derivatives 
CPS   ENDFOR
CPS   
CPS END sfn2r
CPS
C***********************************************************************

      use combi
      use comdti
      implicit none

      integer i, nsl
      real*8 etad, sid

      nsl    =  4
c**** define local coordinates, integration types ****
c**** local coordinates ****
      xd(1)  = -1.0
      xd(2)  =  1.0
      xd(3)  =  1.0
      xd(4)  = -1.0
      yd(1)  = -1.0
      yd(2)  = -1.0
      yd(3)  =  1.0
      yd(4)  =  1.0

      do i = 1, nsl
         w (i,1) = 0.25 * (1.0 + xd(i) * sid) * (1.0 + yd(i) * etad)
         wx(i,1) = 0.25 * xd(i) * (1.0 + yd(i) * etad)
         wy(i,1) = 0.25 * yd(i) * (1.0 + xd(i) * sid)
      end do

      end
