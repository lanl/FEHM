      subroutine  dualex
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
CD1 To backsubstitute to get heat and mass transfer
CD1 solution for dual porosity nodes.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-4-93      G. Zyvoloski   N/A     Initial implementation, but
CD2                                        previous non-YMP versions
CD2                                        of FEHM exist, and the
CD2                                        current version may differ
CD2                                        from these
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dualex.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 15:30:40 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:50:04   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:22   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 None
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
CD4 neq, bp, rb2mf, a21mpf, a21mef, rb2ef, a21epf, a21eef, wb11, wb12,
CD4 wb21, wb22, r3mf, a32mpf, a32mef, r3ef, a32epf, a32eef, tb21, tb22,
CD4 ico2, ieosc, ieos, phi, t, s, 
CD4 
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 None
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
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 a1           real*8      Parameter used in function definition
CD5 a2           real*8      Parameter used in function definition
CD5 a3           real*8      Parameter used in function definition
CD5 a4           real*8      Parameter used in function definition
CD5 b1           real*8      Parameter used in function definition
CD5 b2           real*8      Parameter used in function definition
CD5 id           int         Do loop parameter over all nodes
CD5 x1m          real*8      Corrections for first degree of freedom
CD5                              variable at each primary node
CD5 x1e          real*8      Corrections for second degree of freedom
CD5                              variable at each primary node
CD5 idp          int         Index of first matrix node
CD5 idpp         int         Index of second matrix node
CD5 rtm          real*8      Parameter used in calculating corrections
CD5 rte          real*8      Parameter used in calculating corrections
CD5 x2m          real*8      Parameter used in pressure correction
CD5                              calculation at first matrix node
CD5 x2e          real*8      Parameter used in temperature or
CD5                               saturation correction calculation at
CD5                               first matrix node
CD5 x3m          real*8      Parameter used in pressure correction
CD5                              calculation at second matrix node
CD5 x3e          real*8      Parameter used in temperature or
CD5                               saturation correction calculation at
CD5                               second matrix node
CD5 ieosd        int         Equation of state flag at current node
CD5 
CD5 
CD5 
CD5 Local Subprograms
CD5
CD5 Identifier     Type        Description
CD5 
CD5 alm            real*8      Function used in calculating corrections
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
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
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.4.7 Dual-porosity formulation
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN dualex
CPS 
CPS Perform preliminary calculations
CPS 
CPS FOR each node
CPS   Calculate indexes for matrix nodes
CPS   Calculate corrections for first set of matrix nodes
CPS   Calculate corrections for second set of matrix nodes
CPS   IF this is a heat and mass transfer solution
CPS     Set equation of state parameter
CPS     IF it is a single phase fluid
CPS       Compute pressure and temperature in first set of matrix nodes
CPS     ENDIF
CPS     IF it is a two-phase fluid
CPS       Compute pressure and saturation at first set of matric nodes
CPS     ENDIF
CPS     IF it is a superheated vapor
CPS       Compute pressure and temperature at first set of matrix nodes
CPS     ENDIF
CPS     IF it is a heat transfer only solution
CPS       Compute temperature at first set of matrix nodes
CPS     ENDIF
CPS     Set equation of state parameter at second set of matrix nodes
CPS     IF it is a single phase fluid
CPS       Compute pressure and temperature in second set of matrix nodes
CPS     ENDIF
CPS     IF it is a two-phase fluid
CPS       Compute pressure and saturation at second set of matric nodes
CPS     ENDIF
CPS     IF it is a superheated vapor
CPS       Compute pressure and temperature at second set of matrix nodes
CPS     ENDIF
CPS     IF it is a heat transfer only solution
CPS       Compute temperature at second set of matrix nodes
CPS     ENDIF
CPS   ELSEIF this is a air water nonisothermal problem
CPS     Compute pressures and saturations at matrix nodes
CPS   ENDIF
CPS   
CPS ENDFOR each node
CPS 
CPS END dualex
CPS 
C**********************************************************************

      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use combi
      use comdti
      use comai
      implicit none

      real*8 a1,a2,b1,b2,x1m,x1e,rtm,rte,x2m,x2e,x3m,x3e,alm
      integer id,idp,idpp,ieosd

c***     linear algebra     ***
      alm(a1,a2,b1,b2) =  a1*b1+a2*b2
      
      do id=1,neq
         idp =  id+neq
         idpp =  id+neq+neq
c     
c     calculate node idp solution
c     x2=ab22i*(rb2-a21*x1)
c     find rt=rb2-a21*x1
c     
         x1m =  bp(id)
         x1e =  bp(idp)
         rtm =  rb2mf(id)-alm(a21mpf(id),a21mef(id),x1m,x1e)
         rte =  rb2ef(id)-alm(a21epf(id),a21eef(id),x1m,x1e)
         x2m =  alm(wb11(id),wb12(id),rtm,rte)
         x2e =  alm(wb21(id),wb22(id),rtm,rte)
c     
c     calculate node idpp solution
c     x3=a33i*(r3-a32*x2)
c     find rt=r3-a32*x2
c     
         rtm =  r3mf(id)-alm(a32mpf(id),a32mef(id),x2m,x2e)
         rte =  r3ef(id)-alm(a32epf(id),a32eef(id),x2m,x2e)
         x3m =  alm(tb11(id),tb12(id),rtm,rte)
         x3e =  alm(tb21(id),tb22(id),rtm,rte)
c     
c     update appropriate variable set
c     n-r corrections
c     
         if(ico2.eq.0) then
c     
c     corrections for heat / mass solution
c     
            strd =  1.0
            if ( ieosc .ne. 0) strd =  0.5
            ieosd =  ieos(idp)
            if ( ieosd .eq. 1 ) then
               phi(idp) =  phi(idp)-x2m*strd
               t(idp) =  t(idp)-x2e*strd
            endif
            if ( ieosd .eq. 2 ) then
               phi(idp)=phi(idp)-x2m*strd
               s(idp)=s(idp)-x2e*strd
            endif
            if ( ieosd .eq. 3 )  then
               phi(idp)=phi(idp)-x2m*strd
               t(idp)=t(idp)-x2e*strd
            endif
            if ( ieosd .eq. 4 )  then
               t(idp)=t(idp)-x2e*strd
            endif
            strd=1.0
            if ( ieosc .ne. 0) strd=0.5
            ieosd=ieos(idpp)
            if ( ieosd .eq. 1 )  then
               phi(idpp)=phi(idpp)-x3m*strd
               t(idpp)=t(idpp)-x3e*strd
            endif
            if ( ieosd .eq. 2 )  then
               phi(idpp)=phi(idpp)-x3m*strd
               s(idpp)=s(idpp)-x3e*strd
            endif
            if ( ieosd .eq. 3 )  then
               phi(idpp)=phi(idpp)-x3m*strd
               t(idpp)=t(idpp)-x3e*strd
            endif
            if ( ieosd .eq. 4 )  then
               t(idpp)=t(idpp)-x3e*strd
            endif
         else if(ico2.lt.0) then
c     
c     corrections for air/water isthermal solution
c     
            phi(idp)=phi(idp)-x2m
            s(idp)=s(idp)-x2e
            phi(idpp)=phi(idpp)-x3m
            s(idpp)=s(idpp)-x3e
         endif
      enddo
      
      return
      end
