      subroutine methane_properties(iflg,iphase,var1,var2,var3,
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
!D1 To calculate methane properties and derivatives.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/methane_properties.f_a  $
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
      real*8 cpms,cpml,cpmv
      real*8 lhmsl,lhmlv
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
c     remove following declarations after code varification RJP 11/25/03
      real*8 roc0,tref,tempc,pcl0,drocp0
c     remove following declarations later
      real*8 rolref, comw, rcomd,pref
      real*8 vla0,vlpa1,vlpa2,vlpa3,vlta1,vlta2,vlta3
      real*8 vlpta,vlp2ta,vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3
      real*8 vlptb,vlp2tb,vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tlx,tl2x,tlx2
      real*8 viln1,viln2,viln3,viln,vil
      real*8 vild1,vild2,vild3,vild
c cpms - heat capacity for methane solid
c cpml - heat capacity for methane liquid
c cpmv - heat capacity for methane vapor 
c tm0 - reference temperature for methane
c lhmsl - latent heat for methane solid to liquid phase change
c lhmlv - latent heat for methane liquid to vapor phase change
      parameter(cpms=4.0d-3,cpml=4.0d-3,cpmv=0.002,pm0=0.0,tm0=0.0)
      parameter(lhmsl=0.1,lhmlv=0.84)                      
c denms - reference density for methane solid
c ddenmsp - derivative of solid density wrt pressure 
c ddenmst - derivative of solid density wrt temperature
      parameter(denms=1000.,ddenmsp=0.5,ddenmst=-0.5)
c denml - reference density for methane liquid
c ddenmlp - derivative of liquid density wrt pressure 
c ddenmlt - derivative of liquid density wrt temperature
      parameter(denml=1000.,ddenmlp=0.5,ddenmlt=-0.5)
c denmv - reference density for methane gas    
c prefmv - reference pressure for perfect gas law   
c pterm,tterm - intermediate calculations
      parameter(denmv=1.0,prefmv=0.1)       
c visms  reference viscosity for methane solid  
      parameter(visms=1.d10)                                     
c visml - reference viscosity for methane liquid
c pvisml1 - parameter used with liquid viscosity 
c pvisml2 - parameter used with liquid viscosity 
c tvisml - cutoff temperature for constant liquid viscosity 
      parameter(visml=0.005,pvisml1=5.0)
      parameter(pvisml2=0.1,tvisml=0.0001)
c vismv- reference  viscosity for methane vapor  
c dvismvp - derivative of vapor viscosity wrt pressure 
c dvismvt - derivative of vapor viscosity wrt temperature
      parameter(vismv=1.d-5,dvismvp=-8.e-7,dvismvt=4.d-8)
c tphasesl - phase change temperature for solid to liquid 
c tphaselv - phase change temperature for liquid to vapor 
c dtphaseslp - derivative phase change temperature (s-l) wrt pressure
c dtphaselvp - derivative phase change temperature (l-v) wrt pressure 
c   adjust phase temepratures so methane always stays as a gas
c     parameter(tphasesl=0.0,tphaselv=-100.0)                   
      parameter(tphasesl=-100.0,tphaselv=0.0)                   
      parameter(dtphaseslp=1.0d-5,dtphaselvp=1.0d-5)                   
c pphasesl - phase change pressure for solid to liquid 
c pphaselv - phase change pressure for liquid to vapor 
      parameter(pphasesl=0.1,pphaselv=0.1)                   

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

      if((iflg.eq.1).and.(iphase.eq.1)) then
c  solid enthalpy and derivative wrt pressure and temperature
        prop=cpms*(var2-tm0)
        der1 = 0.0
        der2 = cpms
      else if((iflg.eq.1).and.(iphase.eq.2)) then
c  liquid enthalpy and derivative wrt pressure and temperature
        prop=cpml*(var2-tm0) + lhmsl
        der1 = 0.0
        der2 = cpml
      else if((iflg.eq.1).and.(iphase.eq.3)) then
c  vapor enthalpy and derivative wrt pressure and temperature
         prop = 852.18-(13.163*var1)+(2.65667*(var2+273.15d0-270.d0))
         prop = prop*1d-3
         der1 = -13.163d0
         der1 = der1*1d-3
         der2 = 2.265667d0
         der2 = der2*1d-3
      else if((iflg.eq.2).and.(iphase.eq.1)) then
c  solid density and derivative wrt pressure and temperature
        prop=denms + ddenmsp*(var1-pm0) + ddenmst*(var2-tm0)
        der1 = ddenmsp
        der2 = ddenmst
      else if((iflg.eq.2).and.(iphase.eq.2)) then
c  liquid density and derivative wrt pressure and temperature
        prop=denml + ddenmlp*(var1-pm0) + ddenmlt*(var2-tm0)
        der1 = ddenmlp
        der2 = ddenmlt
      else if((iflg.eq.2).and.(iphase.eq.3)) then
c    vapor density and derivative wrt pressure and temperature
c         prop=-0.1548d0*var2-0.729d0+(7.8143d0*var1)
c         der1=7.8143d0
c         der2=-0.1548d0
         prop=-0.6534d0*var2+89.95d0+(10.321d0*(var1-10.d0))
         der1=10.321d0
         der2=-0.6534d0
c     Modified RJP 11/25/2003 to air for code varification
c     remove and change to top one later
c         roc0 = 1.292864
c         tref= 4.
c         tempc=273.0/(tref+273.0)
c         pcl0 = 0.101325
c         drocp0=roc0*tempc/pcl0
c         prop=drocp0*var1
c         der1=drocp0
c         der2=0.
      else if((iflg.eq.3).and.(iphase.eq.1)) then
c  solid viscosity and derivative wrt pressure and temperature
c  solid is immobile
        prop=visms
        der1 = 0.0
        der2 = 0.0
      else if((iflg.eq.3).and.(iphase.eq.2)) then
c  liquid viscosity and derivative wrt pressure and temperature
        if(var2.ge.tvisml) then
         prop= pvisml1*visml/(var2+pvisml2)
         der1 = 0.0
         der2 = -pvisml1*visml/(var2+pvisml2)**2
        else
         prop=pvisml1*visml/(tvisml+pvisml2)
         der1 = 0.0
         der2 = 0.0
        endif
      else if((iflg.eq.3).and.(iphase.eq.3)) then
c  vapor viscosity and derivative wrt pressure and temperature
c        prop=vismv + dvismvp*(var1-pm0) + dvismvt*(var1-tm0)
c        der1 = dvismvp
c        der2 = dvismvt
c        prop = 1.0613d-5
c        der1 = 1.8236d-7
c        der2 = 3.286d-8
         prop = 9.6965d0+(0.3751d0*var1)+
     &        (0.032233d0*(var2+273.15d0-270.d0))
         prop = prop*1d-6
         der1 = 0.3751d0
         der1 = der1*1d-6
         der2 = 0.032233d0
         der2 = der2*1d-6
c     Modified RJP 11/25/2003 to air for code varification
c     remove and change to top one later
c         prop=  182.e-7
c         der1 = 0.d0
c         der2 = 0.d0
      else if((iflg.eq.4).and.(iphase.eq.1)) then
c  phase-change(s-l) temperature and derivative wrt pressure
        prop = tphasesl + dtphaseslp*(var1-pphasesl)
        der1 = dtphaseslp
        der2 = 0.0
      else if((iflg.eq.4).and.(iphase.eq.2)) then
c  phase-change(l-v) temperature and derivative wrt pressure
c        prop = tphaselv + dtphaselvp*(var1-pphaselv)
c        der1 = dtphaselvp
c        der2 = 0.0
c gaz 073004 this insures gas phase for methane
         prop = -100000.
         der1 = 0.
         der2 = 0.
      else if((iflg.eq.5).and.(iphase.eq.1)) then
c  phase-change(s-l) pressure and derivative wrt temperature
        prop = pphasesl + 1.0/dtphaseslp*(var2-tphasesl)
        der1 = 0.0
        der2 = 1.0/dtphaseslp
      else if((iflg.eq.5).and.(iphase.eq.2)) then
c  phase-change(l-v) pressure and derivative wrt temperature
        prop = pphaselv + 1.0/dtphaselvp*(var2-tphaselv)
        der1 = 0.0
        der2 = 1.0/dtphaselvp
      else if((iflg.eq.6).and.(iphase.eq.1)) then
c  liquid-vapor capillary pressure and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if((iflg.eq.6).and.(iphase.eq.2)) then
c  solid-liquid capillary pressure and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if((iflg.eq.7).and.(iphase.eq.1)) then
c  vapor relative permeability and derivative wrt liq. saturation
        prop = 1.0-var3      
        der1 = 0.0
        der2 = -1.0
      else if((iflg.eq.7).and.(iphase.eq.2)) then
c  liquid relative permeability and derivative wrt liq. saturation
        prop = var3      
        der1 = 0.0
        der2 = 1.0
      else if((iflg.eq.8).and.(iphase.eq.1)) then
c  solid relative permeability and derivative wrt liq. saturation
        prop = 0.0       
        der1 = 0.0
        der2 = 0.0
      else if((iflg.eq.8).and.(iphase.eq.2)) then
c  liquid relative permeability and derivative wrt liq. saturation
        prop = var3      
        der1 = 0.0
        der2 = 1.0
      else if(iflg.eq.9) then
c  source enthalpy returned when a temperature is specified            
         prop = 852.18-(13.163*var1)+(2.65667*(var2+273.15d0-270.d0))
         prop = prop*1d-3
c        if(var2.le.tphasesl) then
c          prop=cpms*(var2-tm0)
c        else if(var2.ge.tphaselv) then
c          prop=cpmv*(var2-tm0) + lhmlv
c        else
c          prop=cpml*(var2-tm0) + lhmsl
c        endif
      else if(iflg.eq.11) then
c  water-methane relative permeability based on fractions    	    
c  for methane gas
c  indepedent of fluid state
        prop = (1.0 - var3 - var4)
	der1 = 0.0
	der2 = 0.0
	der3 = -1.0
	der4 = -1.0        
        if (prop.le.0.0) then
           prop = 0.0
           der1 = 0.0
           der2 = 0.0
           der3 = 0.0
           der4 = 0.0
        endif
c     Modified by Rajesh 11/25/03
c        prop = (1.0 - var3 - var4)**4
c	der1 = 0.0
c	der2 = 0.0
c	der3 = -4*(1.0-var3-var4)**3
c	der4 = -4*(1.0-var3-var4)**3        
c        if (prop.le.0.0) then
c           prop = 0.0
c           der1 = 0.0
c           der2 = 0.0
c           der3 = 0.0
c           der4 = 0.0
c        endif
      endif
      return 
      end
