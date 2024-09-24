      subroutine  air_cp(temp_celsius, heatcap, dheatcapt)
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
CD1 To compute the heat capacity of air (J/kg-K) and the
CD1 derivative with respect to temperature.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 8-14-94      B. Robinson    22      Initial Implementation
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/air_cp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 13:10:42 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.1   Mon Jan 29 09:52:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.0   08/22/94 11:21:04   llt
CD2 Original version
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name          Type        Description
CD3 
CD3 temp_celsius  real*8      Temperature in deg. C
CD3 heatcap       real*8      Heat capacity of air
CD3 dheatcapt     real*8      Derivative of heat capacity
CD3                              with temperature
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
CD4 None
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4 None
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
CD5 a_air        real*8      Constants in heat capacity
CD5                             calculation
CD5 a_aird       real*8      Constants in heat capacity
CD5                             derivative calculation
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
CD5 temp_celsius real*8      Temperature (deg. C)
CD5 tplus        real*8      Running product in heat capacity
CD5                             calculation
CD5 tplusd       real*8      Running product in heat capacity
CD5                             derivative calculation
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
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
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA Sychev, V. V. et al., Thermodynamic Properties of Air, Hemisphere
CDA Publishing Corp., 1988.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN air_cp
CPS 
CPS   Initialize heat capacity and derivative to 0
CPS   Initialize running products
CPS
CPS   FOR each term in series
CPS     Compute running product terms
CPS     Compute next term in series
CPS   ENDFOR
CPS
CPS END air_cp
CPS 
C**********************************************************************
      implicit none
      integer i
      real*8  a_air(4),a_aird(4),temp_celsius,heatcap,dheatcapt,
     &     tplus, tplusd
      data a_air / 1003.7, .025656, .00045457, -2.7107e-7 /
      data a_aird / 0., .025656, .00090914, -8.1321e-7 /
      
      heatcap = a_air(1)
      dheatcapt = 0.0
      tplus = 1.
      
      do i = 2,4
         tplusd = tplus
         tplus = tplus*temp_celsius
         heatcap = heatcap + a_air(i)*tplus
         dheatcapt = dheatcapt + a_aird(i)*tplusd
      enddo

      return
      end
c gaz 121418
c need to add solubility as function of temperature      
      subroutine air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
c calculate air solubility in water 
       use com_prop_data, only :   xnl_max  
       implicit none
       real*8 tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt, alpha0
       real*8 xtol,alpha, dalpca,dalphat,tsolmax, alpha_tol       
       integer imod_sol
        parameter (imod_sol = 1)
        parameter (alpha0 = 1.6111d-04)     
        parameter(xtol=1.d-16, tsolmax = 300., alpha_tol = 1.d-9)  
        xnl_max=0.1      
c gaz 121618 disable temperature dependant Henry's law   
        if(imod_sol.ne.0) then           
         if(tl.le.tsolmax) then
          alpha=alpha0*(1.-tl/tsolmax)**2 + alpha_tol
          dalpca=0.0
          dalphat = -alpha0*2.*(1.-tl/tsolmax)/tsolmax
          xnl=alpha*pcl + xtol
          dxnlp=0.0
          dxnlpc=alpha
          dxnlt=dalphat*pcl
         else
          xnl=alpha_tol*pcl
          dxnlp=0.0
          dxnlpc=alpha_tol
          dxnlt=0.0     
         endif
        else
          xnl=alpha0*pcl + alpha_tol
          if(xnl.lt.xnl_max) then
           dxnlp=0.0   
           dxnlpc=alpha0
           dxnlt=0.0
          else
           xnl=xnl_max
           dxnlp=0.0
           dxnlpc=alpha0
           dxnlt=0.0           
          endif
        endif
        return
c gaz 112018 account for pcl < 0    (this produces the old model)    
          if(pcl.gt.-100.) then
           xnl=alpha0*pcl
           dxnlp=0.0
           dxnlpc=alpha0
           dxnlt=0.0
          else
           xnl=0.0
           dxnlp=0.0
           dxnlpc=alpha0
           dxnlt=0.0              
          endif
c
c       derivatives of liquid mass fraction
c
         
      return
      end
      subroutine gas_frac_2phase_iso(iflg)
c
c gaz 040224 
c convert total mass fraction to gas and dissolved gas      
c
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comxi
      use comwt
      use davidi
      use comflow 
      use comsplitts 
      use com_prop_data, only : xnl_ngas, ihenryiso, xnl_max, xnl_ini,
     &   xnl_chng_low,xnl_chng_high
      implicit none

      integer iflg,ndummy,ico2d,ndum,nndum,ieoss,iieoss,iwelbs,i,mid
c gaz 121821 added n0dum for iso water properties using fluid_control_prop.f
      integer n0dum
      integer mi,i1,i2,ilev,mlev,il,md,irlpsv,irlptsv,icapsv,nr1
      integer nr2,irdofsv 
      integer  ii,ij,im,inode,iwm,j,ja,k,kb
c gaz 110819 removed tref, pref (now global)       
      real*8 pssv,ssv,phisv,dmpfd,dmefd,dqd,rqd,qcd
      real*8 inflow_thstime,inen_thstime,denht,deneht
      real*8 dels,delp,dfdum11,dfdum12,dfdum21,dfdum22
      real*8 dfdum11i,dfdum12i,dfdum21i,dfdum22i,detdf
      real*8 fdum01,fdum02,sx1d,phidum,phi_dif,phi_1,phi_2
      real*8 hmax, hmin, hmid
      real*8 cden_correction
      character*80 form_string
      character*80 dum_air
      real*8 pref_1_2,pref_2_1,s_1_2
      real*8 t_low
c gaz 103019 added strd_satneg for under relaxation when neg saturations 
c gaz 081623   
      real*8 strd_satneg, strd_old, phi_unsat_to_sat, p_uzmin
      real*8 strd1, strd2, phi_old
c gaz 040424
      real*8 mass_tot, gas_tot
c gaz 081823 ieosd
      integer ieosd,nr_test
      integer i_t_bad 
c gaz debug 112119      
      integer ieos_c
c gaz 071223
      integer ico2_sv
c gaz 092723
      integer isotherup  
      if(iflg.eq.1) then
c gaz 040424 for now, just isothermal with dissolved gas
c gaz 0402024 added  cnlof(i)
       do i = 1, n0 
        if(ihenryiso.ne.0.and.s(i).gt.0.0d0) then
         cnlf(i) = frac_gas_iso(i)
         cnlof(i) =  cnlf(i)
        else 
         cnlf(i) =  0.0d0
         cnlof(i) =  cnlf(i)
        endif     
       enddo
25    continue
      else
      endif
          
      end


