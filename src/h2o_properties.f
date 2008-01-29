      subroutine h2o_properties(iflg,iphase,var1,var2,var3,
     &                              var4,prop,der1,der2,der3,der4)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1 
!D1 To calculate water properties and derivatives (new subroutine form).
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Date ?, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/h2o_properties.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 ?
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      implicit none
      integer iflg,iphase
      real*8 var1,var2,var3,var4
      real*8 prop,der1,der2,der3,der4
      real*8 tm0,pm0
      real*8 cpms,cpml,cpml1,cpmv
      real*8 lhmsl,lhmlv,enthl0
      real*8 denms,ddenmsp,ddenmst
      real*8 denml,ddenmlp,ddenmlt
      real*8 denmv,prefmv,pterm,tterm   
c sold phase is immobile (set visms very large)
      real*8 visms  
      real*8 visml,pvisml1,pvisml2,tvisml
      real*8 vismv,dvismvp,dvismvt    
      real*8 tphasesl,tphaselv        
      real*8 dtphaseslp,dtphaselvp        
      real*8 pphasesl,pphaselv        
      real*8 cap_max                  
      real*8 fracw_min,fracw_max            
c     remove following declarations later
      real*8 rolref, comw, rcomd,pref
      real*8 vla0,vlpa1,vlpa2,vlpa3,vlta1,vlta2,vlta3
      real*8 vlpta,vlp2ta,vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3
      real*8 vlptb,vlp2tb,vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tref,tlx,tl2x,tlx2
      real*8 viln1,viln2,viln3,viln,vil
      real*8 vild1,vild2,vild3,vild
c cpms - heat capacity for water solid
c cpml - heat capacity for water liquid
c cpmv - heat capacity for water vapor 
c tm0 - reference temperature for water
c pm0 - reference pressure for water
c lhmsl - reference value for liquid enthalpy
c lhmlv - latent heat for water liquid to vapor phase change
      parameter(cpms=4.0d-3,cpml=4.178d-3,cpml1=9.897e-4,cpmv=0.002)
      parameter (pm0=7.0,tm0=5.0)
      parameter(lhmsl=0.1,lhmlv=2.0,enthl0=2.796e-3)                      
c denms - reference density for water solid
c ddenmsp - derivative of solid density wrt pressure 
c ddenmst - derivative of solid density wrt temperature
      parameter(denms=1006.6,ddenmsp=0.3558,ddenmst=-5.371e-2)
c denml - reference density for water liquid
c ddenmlp - derivative of liquid density wrt pressure 
c ddenmlt - derivative of liquid density wrt temperature
      parameter(denml=1006.6,ddenmlp=0.3558,ddenmlt=-5.371e-2)
c denmv - reference density for water gas    
c prefmv - reference pressure for perfect gas law   
c pterm,tterm - intermediate calculations
      parameter(denmv=1.0,prefmv=0.1)       
c visms  reference viscosity for water solid  
      parameter(visms=1.d10)                                     
c visml - reference viscosity for water liquid
c pvisml1 - parameter used with liquid viscosity 
c pvisml2 - parameter used with liquid viscosity 
c tvisml - cutoff temperature for constant liquid viscosity 
c      parameter(visml=0.005,pvisml1=5.0)
c      parameter(pvisml2=0.1,tvisml=0.0001)
      parameter(visml=0.0011,pvisml1=42.368)
      parameter(pvisml2=27.368,tvisml=-10000.)
c vismv- reference  viscosity for water vapor  
c dvismvp - derivative of vapor viscosity wrt pressure 
c dvismvt - derivative of vapor viscosity wrt temperature
      parameter(vismv=1.d-5,dvismvp=-8.e-7,dvismvt=4.d-8)
c tphasesl - phase change temperature for solid to liquid 
c tphaselv - phase change temperature for liquid to vapor 
c dtphaseslp - derivative phase change temperature (s-l) wrt pressure
c dtphaselvp - derivative phase change temperature (l-v) wrt pressure 
c     parameter(tphasesl=0.0,tphaselv=100.0)                   
      parameter(tphasesl=-100.0,tphaselv=800.0)                   
c     parameter(dtphaseslp=1.0,dtphaselvp=250.0)                   
      parameter(dtphaseslp=1.0e-5,dtphaselvp=250.0e-5)                   
c pphasesl - phase change pressure for solid to liquid 
c pphaselv - phase change pressure for liquid to vapor 
      parameter(pphasesl=0.1,pphaselv=0.1)                   
c cap_max - maximum capillary pressure wrt water_frac 
c      parameter(cap_max=1e-5)                  
c      parameter(cap_max=2.0)                  
      parameter(cap_max=0.1)                  
c fracw_min - minimum water fraction for methane production
      parameter(fracw_min=0.01)               
c fracw_max - maximum water fraction for methane production
      parameter(fracw_max=1.01)               

c input ----------
c iflag - designator for type of calculation
c iphase - phase state                                   
c var1 - variable 1(pressure)
c var2 - variable 2(temperature)
c var3 - variable 3
c var4 - variable 4
c output ----------
c prop - property of material (methane)
c der1 - derivative of property wrt variable 1
c der2 - derivative of property wrt variable 2
c der3 - derivative of property wrt variable 3
c der4 - derivative of property wrt variable 4
c -----------------

c iflg=1, enthalpy and derivative wrt pressure temperature
c iflg=2, density and derivative wrt pressure temperature
c iflg=3, viscosity and derivative wrt pressure temperature
c iflg=4, saturation (or melting) temperature and derivative wrt pressure
c iflg=5, saturation (or melting) pressure and derivative wrt temperature
c iflg=10, capillary pres (water frac) and derivative wrt water frac

      if(iflg.eq.1.and.iphase.eq.1) then
c  solid enthalpy and derivative wrt pressure and temperature
        prop=cpms*(var2-tm0)
        der1 = 0.0
        der2 = cpms
      else if(iflg.eq.1.and.iphase.eq.2) then
c  liquid enthalpy and derivative wrt pressure and temperature
        prop = cpml1*(var1-pm0) + cpml*(var2-tm0) + enthl0
        der1 = cpml1
        der2 = cpml 
      else if(iflg.eq.1.and.iphase.eq.3) then
c  vapor enthalpy and derivative wrt pressure and temperature
        prop=cpmv*(var2-tm0) + lhmlv
        der1 = 0.0
        der2 = cpmv
      else if(iflg.eq.2.and.iphase.eq.1) then
c  solid density and derivative wrt pressure and temperature
        prop=denms + ddenmsp*(var1-pm0) + ddenmst*(var2-tm0)
        der1 = ddenmsp
        der2 = ddenmst
      else if(iflg.eq.2.and.iphase.eq.2) then
c  liquid density and derivative wrt pressure and temperature
        prop=denml + ddenmlp*(var1-pm0) + ddenmlt*(var2-tm0)
        der1 = ddenmlp
        der2 = ddenmlt
      else if(iflg.eq.2.and.iphase.eq.3) then
c  vapor density and derivative wrt pressure and temperature
c  assume perfect gas law
        tterm=273.0/(var2+273.0)
        pterm=denmv*tterm/prefmv
        prop=pterm*var1
        der1 = pterm
        der2 = (denmv*var1/prefmv)*(-273.0/(var2+273)**2)
      else if(iflg.eq.3.and.iphase.eq.1) then
c  solid viscosity and derivative wrt pressure and temperature
c  solid is immobile 
        prop=visms
        der1 = 0.0
        der2 = 0.0
      else if(iflg.eq.3.and.iphase.eq.2) then
c  liquid viscosity and derivative wrt pressure and temperature
        if(var2.ge.tvisml) then
          prop= pvisml1*visml/(var2+pvisml2)
          der2 = -pvisml1*visml/(var2+pvisml2)**2
          der1 = 0.0
        else
          prop=pvisml1*visml/(tvisml+pvisml2)
          der1 = 0.0
          der2 = 0.0
       endif
      else if(iflg.eq.3.and.iphase.eq.3) then
c  vapor viscosity and derivative wrt pressure and temperature
        prop=vismv + dvismvp*(var1-pm0) + dvismvt*(var1-tm0)
        der1 = dvismvp
        der2 = dvismvt
      else if(iflg.eq.4.and.iphase.eq.1) then
c  phase-change(s-l) temperature and derivative wrt pressure
        prop = tphasesl + dtphaseslp*(var1-pphasesl)
        der1 = dtphaseslp
        der2 = 0.0
      else if(iflg.eq.4.and.iphase.eq.2) then
c  phase-change(l-v) temperature and derivative wrt pressure
        prop = tphaselv + dtphaselvp*(var1-pphaselv)
        der1 = dtphaselvp
        der2 = 0.0
      else if(iflg.eq.5.and.iphase.eq.1) then
c  phase-change(s-l) pressure and derivative wrt temperature
        prop = pphasesl + 1.0/dtphaseslp*(var2-tphasesl)
        der1 = 0.0
        der2 = 1.0/dtphaseslp
      else if(iflg.eq.5.and.iphase.eq.2) then
c  phase-change(l-v) pressure and derivative wrt temperature
        prop = pphaselv + 1.0/dtphaselvp*(var2-tphaselv)
        der1 = 0.0
        der2 = 1.0/dtphaselvp
      else if(iflg.eq.6.and.iphase.eq.1) then
c  liquid-vapor capillary pressure and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if(iflg.eq.6.and.iphase.eq.2) then
c  solid-liquid capillary pressure and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if(iflg.eq.7.and.iphase.eq.1) then
c  vapor relative permeability and derivative wrt liq. saturation
        prop = 1.0-var3      
        der1 = 0.0
        der2 = -1.0
      else if(iflg.eq.7.and.iphase.eq.2) then
c  liquid relative permeability and derivative wrt liq. saturation
        prop = var3      
        der1 = 0.0
        der2 = 1.0
      else if(iflg.eq.8.and.iphase.eq.1) then
c  solid relative permeability and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if(iflg.eq.8.and.iphase.eq.2) then
c  liquid relative permeability and derivative wrt liq. saturation
        prop = var3      
        der1 = 0.0
        der2 = 1.0
      else if(iflg.eq.9) then
c  source enthalpy returned when a temperature is specified            
        if(var2.le.tphasesl) then
          prop=cpms*(var2-tm0)
        else if(var2.ge.tphaselv) then
          prop=cpmv*(var2-tm0) + lhmlv 
        else
          prop = cpml1*(var1-pm0) + cpml*(var2-tm0) + enthl0
        endif
      else if(iflg.eq.10.and.iphase.eq.1) then
c  liquid-vapor capillary pressure and derivative wrt water fraction
        prop = 0.0                   
        der1 = 0.0       
        der2 = 0.0
c     updated by Rajesh to take care of negative saturations 11/22/2003
c     prop = cap_max*(1.0-var1)
c     der1 = -cap_max
c     der2 = 0.0
        prop = 1.0-var1
        if (prop.le.0.0) then
           prop = 0.0
           der1 = 0.0
           der2 = 0.0
        else
           prop = cap_max*prop
           der1 = -cap_max
           der2 = 0.0
        endif
      else if(iflg.eq.11) then
c  water-methane relative permeability based on fractions 
c  for liquid water
c  indepedent of fluid state 
         if(var3.lt.fracw_max) then
            prop = var3                    
            der1 = 0.0        
            der2 = 0.0
            der3 = 1.0
         else
            prop = fracw_max               
            der1 = 0.0        
            der2 = 0.0
            der3 = 0.0
         endif
c     Modified by Rajesh 11/25/03 
c        if(var3.lt.fracw_max) then
c         prop = var3**4                    
c         der1 = 0.0        
c         der2 = 0.0
c         der3 = 4*var3**3
c        else
c         prop = fracw_max               
c         der1 = 0.0        
c         der2 = 0.0
c         der3 = 0.0
c	endif

      else if(iflg.eq.12) then
c  water-methane production factor 
c  for liquid water
c  indepedent of fluid state 
       if(var3.gt.fracw_min) then
        prop = (1.0-var3-var4)/var3             
        der1 = 0.0        
        der3 = (-1.0)/var3-(1.0-var3-var4)/(var3*var3)
        der4 = -1.0/var3
       else
        prop = 0.0                              
        der1 = 0.0        
        der3 = 0.0  
        der4 = 0.0 
       endif
      endif
      return 
      end
