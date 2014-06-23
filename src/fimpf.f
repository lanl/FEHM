      subroutine fimpf
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
CD1 Calculate fraction of variables over a given tolerance.
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
CD2 $Log:   /pvcs.config/fehm90/src/fimpf.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   03/18/94 15:53:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:54   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 None
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
CD4   delat           REAL*8   faar   Given tolerance for gas pressure
CD4   delpt           REAL*8   faar   Given tolerance for pressure
CD4   delst           REAL*8   faar   Given tolerance for saturation
CD4   deltt           REAL*8   faar   Given tolerance for temperature
CD4   fimp            REAL*8   faar   Fraction of variables over a given
CD4                                     tolerance
CD4   ico2            INT      faai   Indicates if noncondensible gas solution
CD4                                     is enabled
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   pci             REAL*8          Gas pressure at each node
CD4   pcio            REAL*8          Last time step gas pressure at each node
CD4   phi             REAL*8          Pressure at each node
CD4   pho             REAL*8          Last time step pressure at each node
CD4   s               REAL*8          Liquid saturation at each node
CD4   so              REAL*8          Last time step saturation at each node
CD4   t               REAL*8          Temperature at each node
CD4   to              REAL*8          Last time step temperature at each node
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
CD5   ampt            REAL*8   Counter for number of variables over given
CD5                              tolerance
CD5   dela            REAL*8   Absolute value of the difference between the
CD5                              current gas pressure at a node and its value
CD5                              at the last time step
CD5   delp            REAL*8   Absolute value of the difference between the
CD5                              current pressure at a node and its value at
CD5                              the last time step 
CD5   dels            REAL*8   Absolute value of the difference between the
CD5                              current saturation at a node and its value at
CD5                              the last time step  
CD5   delt            REAL*8   Absolute value of the difference between the
CD5                              current temperature at a node and its value
CD5                              at the last time step  
CD5   i               INT      Loop variable
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
CPS BEGIN fimpf
CPS 
CPS   initialize count for variables out of tolerance to zero
CPS   
CPS   FOR each node
CPS       compute the difference between the current value of pressure, 
CPS        temperature, and saturation at the node and its value at 
CPS        the last time step
CPS       
CPS       IF noncondensible gas solution is enabled
CPS          compute the difference between the current value of gas 
CPS           pressure at the node and its value at the last time step
CPS       ELSE
CPS          set the difference value within tolerance
CPS       ENDIF
CPS       
CPS       IF any difference exceeds the given tolerance
CPS          increment count for variables out of tolerance by 1
CPS       ENDIF
CPS       
CPS   ENDFOR
CPS   
CPS   compute the fraction of variables out of tolerance
CPS       
CPS END fimpf
CPS
C***********************************************************************

      use comci
      use combi
      use comdi
      use comfi
      use comdti
      use comai
      use davidi
      implicit none

      integer i
      real*8 xdela, xdelp, xdels, xdelt 
      real*8 xdelam, xdelpm, xdelsm, xdeltm 
      real*8 xdelmax, tolv
      parameter(tolv=1.d-6)

      if(impf.eq.0) return

      xdelpm = 0.0d00
      xdeltm = 0.0d00
      xdelsm = 0.0d00
      xdelam = 0.0d00
      do i = 1, n
         xdelp = abs(phi(i)-pho(i))
         xdelt = abs(t(i)-to(i))
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            xdels = abs(s(i)-so(i))
         else
            xdels = 0.d0
         end if

         if (ico2 .gt. 0) then
            xdela = abs(pci(i)-pcio(i))
         else
            xdela = 0.99d00*delat
         endif

         xdelpm = max(xdelp/delpt,xdelpm)
         xdeltm = max(xdelt/deltt,xdeltm)
         xdelsm = max(xdels/delst,xdelsm)
         xdelam = max(xdela/delat,xdelam)
         xdelmax = max(xdelpm,xdeltm,xdelsm,xdelam)

      end do

      if(abs(xdelpm-xdelmax).le.tolv) then
         impf = 1
         fimp = xdelpm
      else if(abs(xdeltm-xdelmax).le.tolv) then
         impf = 2
         fimp = xdeltm
      else if(abs(xdelsm-xdelmax).le.tolv) then
         impf = 3
         fimp = xdelsm
      else if(abs(xdelam-xdelmax).le.tolv) then
         impf = 4
         fimp = xdelam
      endif
      
      end
