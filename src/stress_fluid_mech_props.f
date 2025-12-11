      subroutine stress_fluid_mech_props(iflg,ndummy)
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
CD1 To calculate coefficients and derivatives for fluid mech interaction
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 01-20-2007     G. Zyvoloski   00022   Initial implementation.
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
CD9 2.3.3 Noncondensible gas flow equations
CD9 2.3.7 Sources and sinks
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9       Stress and displacement calculations
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
CPS 
C**********************************************************************
c

c
      use davidi
	use comai
      use comii
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comwt
      use comsi
      implicit none

      integer ndummy,mid,mi,ieosd,kq,iflg
      real*8 pl,tl,dtin
c     
c     
c     dependent variables vap p and sl
c     
c     misc. constants
 
c     
      dtin=1.0/dtot
c     
c     generate coef and derivatives
c     
      do mid=1,neq
         mi=mid+ndummy 
         pl = phi(mi) 
         tl = t(mi) 
c     

      enddo
      
      
      return
      end
      subroutine ps_disp(iflg,ndummy)
c
c  displacement-porosity relationship
c
	use comai
	implicit none
	integer ndummy,mid,mi,ieosd,kq,iflg
      do mid = 1,neq
       
      enddo
      end
