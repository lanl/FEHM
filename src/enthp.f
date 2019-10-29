          real*8 function enthp(mi,td)
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
CD1 To calculate enthalpy at a node as a function of temperature and
CD1 pressure.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/enthp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 15:57:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:11:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:44   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 mi                      int       I   Node number
CD3 td                      int       I   Temperature
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
CD4 ps, phi, cpr, iieos, cel, idof
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
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 tl           real*8      Temperature
CD5 psd          real*8      Porosity
CD5 pl           real*8      Pressure
CD5 cprd         real*8      Rock heat capacity
CD5 iieosd       int         Equation of state parameter
CD5 ela0         real*8      Coefficient in water enthalpy correlation
CD5 elpa1        real*8      Coefficient in water enthalpy correlation
CD5 elpa2        real*8      Coefficient in water enthalpy correlation
CD5 elpa3        real*8      Coefficient in water enthalpy correlation
CD5 elta1        real*8      Coefficient in water enthalpy correlation
CD5 elta2        real*8      Coefficient in water enthalpy correlation
CD5 elta3        real*8      Coefficient in water enthalpy correlation
CD5 elpta        real*8      Coefficient in water enthalpy correlation
CD5 elp2ta       real*8      Coefficient in water enthalpy correlation
CD5 elpt2a       real*8      Coefficient in water enthalpy correlation
CD5 elb0         real*8      Coefficient in water enthalpy correlation
CD5 elpb1        real*8      Coefficient in water enthalpy correlation
CD5 elpb2        real*8      Coefficient in water enthalpy correlation
CD5 elpb3        real*8      Coefficient in water enthalpy correlation
CD5 eltb1        real*8      Coefficient in water enthalpy correlation
CD5 eltb2        real*8      Coefficient in water enthalpy correlation
CD5 eltb3        real*8      Coefficient in water enthalpy correlation
CD5 elptb        real*8      Coefficient in water enthalpy correlation
CD5 elp2tb       real*8      Coefficient in water enthalpy correlation
CD5 elpt2b       real*8      Coefficient in water enthalpy correlation
CD5 x            real*8      Pressure
CD5 x2           real*8      x squared
CD5 x3           real*8      x cubed
CD5 tl2          real*8      tl squared
CD5 tl3          real*8      tl cubed
CD5 enwn1        real*8      Intermediate term in enthalpy calculation
CD5 enwn2        real*8      Intermediate term in enthalpy calculation
CD5 enwn3        real*8      Intermediate term in enthalpy calculation
CD5 enwn         real*8      Intermediate term in enthalpy calculation
CD5 enwd1        real*8      Intermediate term in enthalpy calculation
CD5 enwd2        real*8      Intermediate term in enthalpy calculation
CD5 enwd3        real*8      Intermediate term in enthalpy calculation
CD5 enwd         real*8      Intermediate term in enthalpy calculation
CD5 enw          real*8      Water enthalpy
CD5 
CD5 Local Subprograms
CD5 
CD5 None
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
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN enthp
CPS 
CPS Set parameter values
CPS 
CPS IF the porosity is non-zero and this is a heat and mass solution
CPS   Set coefficients for water enthalpy
CPS   Compute intermediate parameters
CPS   Compute water enthalpy
CPS ENDIF
CPS 
CPS IF the porosity is zero or this is a heat transfer only solution
CPS   Compute rock enthalpy
CPS ENDIF
CPS 
CPS END enthp
CPS 
C**********************************************************************

      use comdi
      use comfi
      use comii
      use comdti
      use comai
      use combi
      implicit none

      integer mi,iieosd
      real*8 td,tl,psd,pl,cprd,ela0,elpa1,elpa2,elpa3,elta1,elta2
      real*8 elta3,elpta,elp2ta,elpt2a,elb0,elpb1,elpb2,elpb3,eltb1
      real*8 eltb2,eltb3,elptb,elp2tb,elpt2b,x,x2,x3,tl2,tl3,enwn1,enwn2
      real*8 enwn3,enwn,enwd1,enwd2,enwd3,enwd,enw
      real*8 p_energy
c gaz 110715
      real*8 dum1,dumb,dumc,value(9)
      integer istate 
c
c calculates enthalpy as a function of t and p
c
      if(igrav.ne.0) then
       p_energy = -grav*cord(mi,igrav)
      else
       p_energy = 0.0d0
      endif
      psd=ps(mi)
      cprd=cpr(mi)
c gaz 081317      
      cprd=1.0d0
      tl=td
      if(psd.ne.0.0.and.idof.ge.2) then
         pl=phi(mi) - phi_inc
        if(iwater_table.ne.1) then
         iieosd=iieos(mi)
c liquid enthalpy
c numerator coefficients
         ela0=cel(1,iieosd)
         elpa1=cel(2,iieosd)
         elpa2=cel(3,iieosd)
         elpa3=cel(4,iieosd)
         elta1=cel(5,iieosd)
         elta2=cel(6,iieosd)
         elta3=cel(7,iieosd)
         elpta=cel(8,iieosd)
         elp2ta=cel(9,iieosd)
         elpt2a=cel(10,iieosd)
c denomenator coefficients
         elb0=cel(11,iieosd)
         elpb1=cel(12,iieosd)
         elpb2=cel(13,iieosd)
         elpb3=cel(14,iieosd)
         eltb1=cel(15,iieosd)
         eltb2=cel(16,iieosd)
         eltb3=cel(17,iieosd)
         elptb=cel(18,iieosd)
         elp2tb=cel(19,iieosd)
         elpt2b=cel(20,iieosd)
         x=pl
         x2=x*x
         x3=x2*x
         tl2=tl*tl
         tl3=tl2*tl
         enwn1=ela0+elpa1*x+elpa2*x2+elpa3*x3
         enwn2=elta1*tl+elta2*tl2+elta3*tl3
         enwn3=elpta*tl*x+elpt2a*tl2*x+elp2ta*tl*x2
         enwn=enwn1+enwn2+enwn3
         enwd1=elb0+elpb1*x+elpb2*x2+elpb3*x3
         enwd2=eltb1*tl+eltb2*tl2+eltb3*tl3
         enwd3=elptb*tl*x+elpt2b*tl2*x+elp2tb*tl*x2
         enwd=enwd1+enwd2+enwd3
         enw=enwn/enwd
         enthp=enw + p_energy
       else
c gaz 110915 (calculates too many variables)
          call h2o_properties_new(4,1,pl,tl,dum1,1,
     &                 dumb,value,dumc)
          enthp = value(4) + p_energy
       endif
      endif
      if(psd.eq.0.0.or.idof.le.1) then
         enthp=cprd*tl
      endif

      return
      end
