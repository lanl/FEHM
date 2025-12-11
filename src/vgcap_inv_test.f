      subroutine vgcap_inv( iflg,sl, slr, smr, alpha, beta, 
     &       slcut, sucut, hlcut, hucut, ac3, ac4, bc3, bc4, hp, dslh)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract with the U.S. Department of Energy (DOE). 
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
CD1 To compute the saturation from capillary pressure and derivatives for the 
CD1 van Genuchten model. The inverse of subroutine vg_cap.f .
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2 4-04-07      george Zyvoloski initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vgcap_inv.f_a  $                                     
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 sl           real*8   I/O    Liquid saturation on input,
CD3                                  normalized liquid saturation on
CD3                                  output
CD3 slr          real*8   I      Residual liquid saturation
CD3 smr          real*8   I      Maximum liquid saturation
CD3 alpha        real*8   I      Van Ganuchten parameter
CD3 beta         real*8   I      Van Ganuchten parameter
CD3 ac1          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac2          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 smcut        real*8   I      Lower cutoff saturation below which a
CD3                                  spline fit is used for capillary
CD3                                  pressure versus saturation
CD3 sucut        real*8   I      Upper cutoff saturation above which a
CD3                                  linear fit is used for capillary
CD3                                  pressure versus saturation
CD3 bc3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 bc4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 hp           real*8   O      Capillary pressure
CD3 dhp          real*8   O      Derivative of capillary pressure with
CD3                                  respect to saturation
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 NONE
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4 
CD4 NONE
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4 
CD4 NONE
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 NONE
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 alamda       real*8      Exponent in correlation
CD5 ds           real*8      Reciprocal of demoninator in normalized
CD5                              saturation expression
CD5 denom        real*8      Denominator in normalized saturation
CD5                              expression
CD5 star         real*8      Normalized saturation
CD5 alpi         real*8      Exponent in correlation
CD5 hp           real*8      Capillary pressure before unit conversion
CD5 dhp          real*8      Derivative of hp with respect to saturation
CD5 termstar1    real*8      Term used in capillary pressure calculation
CD5 termstar2    real*8      Term used in capillary pressure calculation
CD5 termb1       real*8      Term used in capillary pressure calculation
CD5 termb2       real*8      Term used in capillary pressure calculation
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7 
CD7  Requirements from SDN: 10086-RD-2.20-00
CD7    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD7    FEHM Application Version 2.20
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8 
CD8 2.4.4 Relative-permeability and capillary-pressure functions
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN vgcap_inv
CPS 
CPS 
CPS END vgcap_inv
CPS
C**********************************************************************

      implicit none

      integer iflg
      real*8 sl
      real*8 slr
      real*8 smr
      real*8 alpha
      real*8 beta
      
      real*8 ac3
      real*8 ac4
      real*8 bc3
      real*8 bc4
      real*8 hlcut
      real*8 hucut
	real*8 slcut
	real*8 sucut
  
      real*8 hp
      real*8 dhp
      real*8 star
  
      real*8 dslh      
      real*8 alamda
      real*8 denom
      real*8 alpi
      real*8 termh1
      real*8 termh2
      real*8 dtermh2h
      real*8 dstarh
      real*8 hpa

      hpa = abs(hp)
      if(iflg.eq.1) then
       if(hpa.gt.hlcut.and.hpa.lt.hucut) then
c      check if within spline cutoffs    
         alamda = 1.0-1.0/beta
         denom = smr-slr     
         termh1 = (alpha*hpa)**beta
         termh2 = 1.d0/(1.d0+termh1)
         star = (termh2)**alamda
         sl =  star*denom +slr
c  assume beta (vg=n) is > 1         
         dtermh2h = -termh2/(1.d0+termh1)*beta*(alpha*hpa)**(beta-1.d0)
     &   *alpha    
         dstarh = alamda*(1.d0/termh2)**(1.d0-alamda)*dtermh2h        
         dslh = dstarh*denom
	   
       else  if(hpa.le.hlcut) then
c      near sl = 1.0  
        sl = (hpa-bc4)/bc3
        dslh = 1.d0/bc3
       else  if(hpa.ge.hucut) then
c      near sl = 0.0   
        termh1 = slcut/(hucut-ac4) 
	  sl = termh1*(hpa-ac4) 
        dslh = termh1
       endif
      else if(iflg.eq.2) then      
      endif     
      return
      end
