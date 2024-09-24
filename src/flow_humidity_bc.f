          subroutine flow_humidity_bc(iflg,tl,pl,h,qin_ng,enth_avg)    
!***********************************************************************
! Copyright 2016. Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los
! Alamos National Laboratory (LANL), which is operated by Los Alamos
! National Security, LLC for the U.S. Department of Energy. The U.S.
! Government has rights to use, reproduce, and distribute this software.
! Neither the U.S. Government nor Los Alamos National Security, LLC or
! persons acting on their behalf, make any warranty, express or implied,
! or assumes any liability for the accuracy, completeness, or usefulness
! of the software, any information pertaining to the software, or
! represents that its use would not infringe privately owned rights.
!
! The software being licensed may be Export Controlled.   It may not be
! distributed or used by individuals or entities prohibited from having
! access to the software package, pursuant to United States export
! control laws and regulations.  An export control review and
! determination must be completed before LANS will provide access to the
! identified Software.
!***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To calculate the inlet flowrates for dry air and water vapor given
CD1 total flowrate of air(or ngas), T, P, and humidity. Will also output the energy
CD1 inflow.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 12/232015    G. Zyvoloski   N/A     Initial implementation
CD2
C**********************************************************************

      use comii
      use comdti
      use comai
      use combi
      implicit none

c tl-input temperature (input variable)
c pl-input pressure (input variable)
c h-input relative humidity (input variable)
c qin-total input flowrate (input variable)
c qin_ng- dry air inflow rate (output variable)
c qin_h2o- h2o inflow rate (output variable)
c qin_enth- enthalpy inflow rate (output variable)
c enth_avg- average enthapy inflow (output variable)
c calculated enthalpy does not include grav term

      integer iflg
      real*8 tl,pl,pl_in,h,qin,qin_ng,qin_h2o,qin_enth,enth_avg
      real*8 w_vap_pres,h_vap_pres,w_vap_den,a_vap_density
      real*8 xnv_a,xnv_h2o,dcpat,w_vap_enth,xtol
      real*8 a_vap_enth,cpa,drocpc,pcl0,v_den,v_enth
      real*8 dumv
      real*8 psatl
      real*8 dum1,dumb,dumc,value(9),value_a(9)
      integer istate, ifail  
      parameter(xtol=1.d-16)      
c      
       pcl0=0.101325
       roc0=1.292864
           
       if(iflg.eq.0) then
c
c initialization
c
       else if(iflg.eq.1) then
c calculate input water vapor density and enthalpy
         w_vap_pres = psatl(tl,0.d0,0.d0,dumv,dumv,0,0.d0)
         h_vap_pres = h*w_vap_pres
c gaz 060820 use table when available  
         if(iwater_table.eq.0) then  
          call h2o_v_den(1,h_vap_pres,tl,w_vap_den)
         else
          call h2o_properties_new(4,3,h_vap_pres,tl,dum1,istate,
     &                 dumb,value,dumc) 
c both w_vap_den  and  ew_vap_enth available from one call to  h2o_properties_new    
          w_vap_den = value(1)
          w_vap_enth = value(4)
         endif
c         
c calculate input air fraction
c
c gaz 070820 added air table EOS
      if(iair_table.eq.1) then 
       if(pl.le.pmin_air_tabl(1)) then
        pl_in = pmin_air_tabl(1)
       else
        pl_in = pl
       endif    
       call air_properties_new(4,3,pl_in,tl,dum1,istate,
     &                 dumb,value_a,dumc)            
       a_vap_density = value_a(1) 
       a_vap_enth = value_a(4)
      else             
**** density of air eqn (31)
       pcl0=0.101325
       roc0=1.292864    
       drocpc=roc0*(273./(tl+273.))/pcl0
       a_vap_density=drocpc*pl
      endif
c      
c air mass fraction 
c
      v_den = a_vap_density + w_vap_den
      xnv_a = a_vap_density / v_den + xtol
      xnv_h2o = 1.0-xnv_a
c
c
c       enthalpy of air and mixture
c
c
c air heat capacity and derivative, (BAR - 8-15-94)
c
        if(iair_table.eq.0) then
         call air_cp(tl, cpa, dcpat)
         a_vap_enth=1.e-6*cpa*tl
        endif
c gaz 060820 use table (enthapy available from previous h2o_properties_new call 
        if(iwater_table.eq.0) 
     &  call h2o_v_enth(1,h_vap_pres,tl,w_vap_enth)
        v_enth=xnv_a*a_vap_enth+xnv_h2o*w_vap_enth

c       
c
c calculate flowrates
c
        qin = 1.0
        qin_ng= xnv_a*qin
        qin_h2o= xnv_h2o*qin
        qin_enth= a_vap_enth*qin_ng+w_vap_enth*qin_h2o
        enth_avg = qin_enth/qin
       endif 
      return
      end
      
      subroutine h2o_v_den(iflg,pl,tl,denv)
c
c calculates water vapor density as a function of t and p
c vapor density
c
      use comai 
      use combi
      use comii
      implicit none

      integer iieosd,iflg
      real*8 dva0,pl,tl,denv
      real*8 dvpa1,dvpa2,dvpa3,dvta1,dvta2,dvta3,dvpta,dvp2ta,dvpt2a
      real*8 dvpt2b
      real*8 dvb0,dvpb1,dvpb2,dvpb3,dvtb1,dvtb2,dvtb3,dvptb,dvp2tb
      real*8 rnsd3,rnsn,rnsd,rns,rov,rnsn1,rnsn2,rnsn3,rnsd1,rnsd2
      real*8 tlx,tl2x,tlx2,tl2,tl3,tlxv,tlxv2,tl2xv
      real*8 xv,xv2,xv3 
      parameter (iieosd = 1)
     
c numerator coefficients
      dva0=crv(1,iieosd)
      dvpa1=crv(2,iieosd)
      dvpa2=crv(3,iieosd)
      dvpa3=crv(4,iieosd)
      dvta1=crv(5,iieosd)
      dvta2=crv(6,iieosd)
      dvta3=crv(7,iieosd)
      dvpta=crv(8,iieosd)
      dvp2ta=crv(9,iieosd)
      dvpt2a=crv(10,iieosd)
c denomenator coefficients
      dvb0=crv(11,iieosd)
      dvpb1=crv(12,iieosd)
      dvpb2=crv(13,iieosd)
      dvpb3=crv(14,iieosd)
      dvtb1=crv(15,iieosd)
      dvtb2=crv(16,iieosd)
      dvtb3=crv(17,iieosd)
      dvptb=crv(18,iieosd)
      dvp2tb=crv(19,iieosd)
      dvpt2b=crv(20,iieosd)
c
      xv = pl
      xv2 = xv*xv
      xv3 = xv2*xv
      tl2=tl*tl
      tl3=tl2*tl
      tlxv=xv*tl
      tl2xv=tl2*xv
      tlxv2=tl*xv2
      rnsn1=dva0+dvpa1*xv+dvpa2*xv2+dvpa3*xv3
      rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
      rnsn3=dvpta*tlxv+dvpt2a*tl2xv+dvp2ta*tlxv2
      rnsn=rnsn1+rnsn2+rnsn3
      rnsd1=dvb0+dvpb1*xv+dvpb2*xv2+dvpb3*xv3
      rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
      rnsd3=dvptb*tlxv+dvpt2b*tl2xv+dvp2tb*tlxv2
      rnsd=rnsd1+rnsd2+rnsd3
      rns=rnsn/rnsd
      denv=rns                 
      return
      end
      
      subroutine h2o_v_enth(iflg,pl,tl,env)
      use comai, only : phi_inc
      use combi
      use comii
      implicit none

      integer iieosd,iflg
      real*8 tl,pl,cprd,eva0,evpa1,evpa2,evpa3,evta1,evta2
      real*8 evta3,evpta,evp2ta,evpt2a,evb0,evpb1,evpb2,evpb3,evtb1
      real*8 evtb2,evtb3,evptb,evp2tb,evpt2b,x,x2,x3,tl2,tl3,ensn1,ensn2
      real*8 ensn3,ensn,ensd1,ensd2,ensd3,ensd,ens,env
      real*8 xv,xv2,xv3,tl2xv,tlxv2,tlxv
      parameter (iieosd = 1)
c
c calculates water vapor enthalpy as a function of t and p
c
      if(iflg.eq.1)then
         pl=pl - phi_inc
c vapor enthalpy
c numerator coefficients

      eva0=cev(1,iieosd)
      evpa1=cev(2,iieosd)
      evpa2=cev(3,iieosd)
      evpa3=cev(4,iieosd)
      evta1=cev(5,iieosd)
      evta2=cev(6,iieosd)
      evta3=cev(7,iieosd)
      evpta=cev(8,iieosd)
      evp2ta=cev(9,iieosd)
      evpt2a=cev(10,iieosd)
c denomenator coefficients
      evb0=cev(11,iieosd)
      evpb1=cev(12,iieosd)
      evpb2=cev(13,iieosd)
      evpb3=cev(14,iieosd)
      evtb1=cev(15,iieosd)
      evtb2=cev(16,iieosd)
      evtb3=cev(17,iieosd)
      evptb=cev(18,iieosd)
      evp2tb=cev(19,iieosd)
      evpt2b=cev(20,iieosd)
         xv=pl
         xv2=xv*xv
         xv3=xv2*xv
         tl2=tl*tl
         tl3=tl2*tl
         tlxv = tl*xv
         tl2xv =tl2*xv
         tlxv2 =tl*xv2
      ensn1=eva0+evpa1*xv+evpa2*xv2+evpa3*xv3
      ensn2=evta1*tl+evta2*tl2+evta3*tl3
      ensn3=evpta*tlxv+evpt2a*tl2xv+evp2ta*tlxv2
      ensn=ensn1+ensn2+ensn3
      ensd1=evb0+evpb1*xv+evpb2*xv2+evpb3*xv3
      ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
      ensd3=evptb*tlxv+evpt2b*tl2xv+evp2tb*tlxv2
      ensd=ensd1+ensd2+ensd3
      ens=ensn/ensd
      env=ens 
      endif
      
      return
      end      
