      subroutine hydrate_properties(iflg,iphase,var1,var2,var3,
     &                              var4,prop,der1,der2,der3,der4)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************

!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To calculate methane hydrate properties and derivatives.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-feb-03, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/hydrate_properties.f_a  $
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

      use comai
      use commeth
      implicit none

      integer iflg,iphase
      real*8 var1,var2,var3,var4
      real*8 prop,der1,der2,der3,der4
      real*8 tm0,pm0,aip_dis
      real*8 cpms,cpml,cpml1,cpmv
      real*8 lhmsl,lhmlv,enthl0
      real*8 denms,ddenmsp,ddenmst
      real*8 denml,ddenmlp,ddenmlt
      real*8 denmv,prefmv,pterm,tterm   
c solid phase is immobile (set visms very large)
      real*8 visms  
      real*8 visml,pvisml1,pvisml2,tvisml
      real*8 vismv,dvismvp,dvismvt    
      real*8 tphasesl,tphaselv        
      real*8 dtphaseslp,dtphaselvp        
      real*8 pphasesl,pphaselv        
      real*8 pt,diss1,pt0,pt1,dptvar2,tp,dtp                 
      real*8 hdiss,hdiss0,hdiss1                 
      real*8 prop1,ddiss1dvar2,diss
      real*8 propv,propl
      real*8 derv1,derv2,derl1,derl2
      real*8 der11,der21,der22

      real*8 gas_const,temp_conv,v_tol,vard
      parameter (gas_const = 8.314, temp_conv = 273.15, v_tol = 1.e-6)

      real*8 eact2,aterm,hterm,dhterm4,wterm,dwterm3,dpt2
      real*8 gvar,gterm,dgterm3,dgterm4,pdif,dpterm1,dpterm2
      real*8 actterm,dactterm2,porterm,mol_hyd,den_hyd,dhdiss2
      real*8 den_hyd0,mol_gas 
c mol_hyd is kg/mol, den_hyd is mol/m**3
      parameter (mol_gas = 0.016d0, den_hyd0 = 7866.)
      parameter (mol_hyd = 0.1195d0)

 
c cpms - heat capacity for methane solid
c cpml - heat capacity for methane liquid
c cpmv - heat capacity for methane vapor 
c tm0 - reference temperature for methane
c lhmsl - latent heat for methane solid to liquid phase change
c lhmlv - latent heat for methane liquid to vapor phase change
      parameter(cpms=4.0d-3,cpml=4.178d-3,cpml1=9.897e-4,cpmv=0.002)
      parameter (pm0=7.0,tm0=5.0)
      parameter(lhmsl=0.1,lhmlv=2.0,enthl0=2.796e-3)                     
c denms - reference density for methane solid
c ddenmsp - derivative of solid density wrt pressure 
c ddenmst - derivative of solid density wrt temperature
      parameter(denms=910.,ddenmsp=0.,ddenmst=0.)
c denml - reference density for methane liquid
c ddenmlp - derivative of liquid density wrt pressure 
c ddenmlt - derivative of liquid density wrt temperature
      parameter(denml=1006.6,ddenmlp=0.3558,ddenmlt=-5.371e-2)
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
      parameter(tphasesl=0.0,tphaselv=100.0)                   
      parameter(dtphaseslp=1.0,dtphaselvp=250.0)                   
c pphasesl - phase change pressure for solid to liquid 
c pphaselv - phase change pressure for liquid to vapor 
      parameter(pphasesl=0.1,pphaselv=0.1)                   
c 
c pt - phase change pressure (solid to gas)     

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
c iflg=3, dissociation rate of methane hydrate solid                    

      if(iflg.eq.1.and.iphase.eq.1) then
c  solid enthalpy and derivative wrt pressure and temperature
c       prop=cpms*(var2-tm0)
c       der1 = 0.0
c       der2 = cpms
c  vapor enthalpy and derivative wrt pressure and temperature
         propv = 852.18-(13.163*var1)+(2.65667*(var2+temp_conv-270.d0))
         propv = propv*1d-3
         derv1 = -13.163d0
         derv1 = derv1*1d-3
         derv2 = 2.265667d0
         derv2 = derv2*1d-3
c  liquid enthalpy and derivative wrt pressure and temperature
c  note this is liquid water (for use in weighted average below)
         propl = cpml1*(var1-pm0) + cpml*(var2-tm0) + enthl0
         derl1 = cpml1
         derl2 = cpml 
c solid should contain weighted values of individual phases

c subtract dissociation entahlpy 
         if (var2 .lt. 0.d0) then
c    low temperature (var2 is temperature)
            hdiss = hdiss01 + (hdiss11*(var2 + temp_conv))
            dhdiss2 = hdiss11
         else
c    higher temperature (var2 is temperature)
            hdiss = hdiss02 + (hdiss12*(var2 + temp_conv))
            dhdiss2 = hdiss12
         endif
         prop = fracv*propv + fracl*propl 
     &        - hdiss
         der1 = fracv*derv1 + fracl*derl1 
         der2 = fracv*derv2 + fracl*derl2 
     &        - dhdiss2


      else if(iflg.eq.1.and.iphase.eq.2) then
c  liquid enthalpy and derivative wrt pressure and temperature
         prop=cpml*(var2-tm0) + lhmsl
         der1 = 0.0
         der2 = cpml
      else if(iflg.eq.1.and.iphase.eq.3) then
c  vapor enthalpy and derivative wrt pressure and temperature
         prop=cpmv*(var2-tm0) + lhmlv
         der1 = 0.0
         der2 = cpmv
      else if(iflg.eq.2.and.iphase.eq.1) then
c solid density and derivative wrt pressure and temperature
c liquid properties for solid hydrate (ie water)
         propl=denml + ddenmlp*(var1-pm0) + ddenmlt*(var2-tm0)
         derl1 = ddenmlp
         derl2 = ddenmlt 
c gas properties for solid hydrate
         propv=-0.6534d0*var2+89.95d0+(10.321d0*(var1-10.d0))
         derv1=10.321d0
         derv2=-0.6534d0
c solid should contain weighted values of individual phases
         prop = fracv*propv + fracl*propl
         der1 = fracv*derv1 + fracl*derl1
         der2 = fracv*derv2 + fracl*derl2
c gaz 080704  fixed density of hydrate so it has no compressibility
c        weighted values are 60 and 10006 
         prop = 900.
         der1 = 0.
         der2 = 0.
      else if(iflg.eq.2.and.iphase.eq.2) then
c  liquid density and derivative wrt pressure and temperature
         prop=denml + ddenmlp*(var1-pm0) + ddenmlt*(var2-tm0)
         der1 = ddenmlp
         der2 = ddenmlt
      else if(iflg.eq.2.and.iphase.eq.3) then
c  vapor density and derivative wrt pressure and temperature
c  assume perfect gas law
         tterm=temp_conv/(var2+temp_conv)
         pterm=denmv*tterm/prefmv
         prop=pterm*var1
         der1 = pterm
         der2 = (denmv*var1/prefmv)*(-temp_conv/(var2+temp_conv)**2)

      else if(iflg.eq.3.and.idof_meth.ne.7) then
c
c      dissociation rate 
c      var1 is phimeth=phi
c      var2 is t       
c      var3 is water fraction
c      var4 is hydrate fraction
c      pt is the phase transition line
c      da - grain area term
c      oe - exponent of grain area term
c      pe - exponent of fugacity difference
c      le - exponent of hydrate saturation term
c      me - exponent of water saturation term
c      ne - exponent of gas saturation term
c 
c      afhyd   - term for activation energy(input) 
c      porh    - porosity (assume = constant but spatially variable)
c      porh is in module commeth(set to ps(i) in calling routine)
c      volh    - cell volume (assume = constant but spatially variable)
c      volh is in module commeth(set to sx1(i) in calling routine)
c      den_hyd - hydrate density (assume = constant) 
c      e_act   - activation energy
c      gas_const - universal gas constant
c      mol_hyd - molecular weight of hydrate

c  
c  hydrate transition line
c   
c  ptc1, ptc2 are defined in declarations and parameter statements
c
c
c Modified 012005 by GAZ to simplify rate term
c and make function of temperature differance
c
         if(ihyd_diss_type.eq.7) then
c Equilibrium model
            tp = (log(var1)-ptc2)/ptc1 - temp_conv
            dtp = 1.0/(ptc1*var1)
            me =(var2-tp)
            if(me.lt.0) then
c growing (need gas and water)
c    gas fraction + water term
            else
	         
c dissociating (need hydrate)

            endif

         else if(ihyd_diss_type.eq.4) then 
c
c    calculate pressure dependent dissociation temperature
c 
c dar is input (fraction of max available mass)
c da is input (multiplier)
c mer is input (fraction of maximum rate,used in otehr iterations)
            tp = (log(var1)-ptc2)/ptc1 - temp_conv
            dtp = 1.0/(ptc1*var1)
            me =(var2-tp)
            if(me.lt.0) then
c growing (need gas and water)
c    gas fraction + water term
               gvar = (1.0-var3-var4)*var3
               dgterm3 =  -1.0*var3 + (1.0-var3-var4)
               dgterm4 =  -1.0*var3
               oe = volh*da
               prop = oe*gvar*me
               der1 = -oe*gvar*dtp
               der2 = oe*gvar
               der3 = oe*me*dgterm3
               der4 = oe*me*dgterm4
            else
c dissociation (need hydrate)
               gvar = var4
               dgterm3 =  0.0
               dgterm4 =  1.0
               oe = volh*da
               prop = oe*gvar*me
               der1 = -oe*gvar*dtp
               der2 = oe*gvar
               der3 = oe*me*dgterm3
               der4 = oe*me*dgterm4
              
            end if
             
            return
         endif
c
c Modified 072804 by GAZ to generalize rate term
c
         pt = exp((ptc1*(var2+temp_conv))+ptc2)
         dpt2 = pt*ptc1

         if (pt.le.var1) then
c
c    hydrate grows (or reforms)
c
c    activation energy should be divided by gas constant
c
            if(ihyd_grow_type.eq.1) then 
c Kim-Bishnoi type
               eact2 = e_act1r
               den_hyd = 1.0
c    activity term
               actterm   =  exp(-eact2/(var2+temp_conv))
               dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
            else  if(ihyd_grow_type.eq.2) then 
c Sakamoto type
               den_hyd = den_hyd0
               afhydr = 1.0
c    activity term
               actterm   =  exp(e_act1r/(var2+temp_conv)+e_act2r)
               dactterm2 =  -actterm*(e_act1r/(var2+temp_conv)**2)
            else
c return    
               prop = 0.0
               der1 = 0.0
               der2 = 0.0
               der3 = 0.0
               der4 = 0.0
               return
            endif
c      den_hyd = 0.0
c    area term (no variable dependence)
            aterm  =  dar**oer
c    hydrate fraction term (included hydrate frac from start of eq.)
            if(ler+1.0.gt.0.0) then
               hterm   =  var4**(ler+1.0)
               if(ler+1.0.eq.1) then
                  dhterm4 = 1.d0
               else
                  dhterm4 =  (ler+1.0)*var4**(ler)
               endif
            else
               hterm = 1.0
               dhterm4 = 0.0
            endif
c    water fraction term
            if(mer.gt.0.0) then
               if(var3.lt.0) then
                  vard = 0.0
               else
                  vard = var3
               endif
               wterm   =  vard**mer 
               if(mer.ne.1) then
                  dwterm3 =  mer*vard**(mer-1.0)
               else
                  dwterm3 = 1.0d0
               endif
            else
               wterm = 1.0
               dwterm3 = 0.0
            endif
c    gas fraction term 
            gvar = 1.0-var3-var4
            if(ner.gt.0.0) then
               gterm   =  gvar**ner
               if(ner.ne.1) then
                  dgterm3 =  -ner*gvar**(ner-1.0)
                  dgterm4 =  dgterm3
               else
                  dgterm3 = 1.0d0
                  dgterm4 = 1.0d0
               endif
            else
               gterm = 1.0
               dgterm3 = 0.0
               dgterm4 = 0.0
            endif
c    pressure (fugacity) term
c    pdif need to be able to change sign
c    might restrict pe to greater than or equal to 1
            pdif   =  pt-var1
            pterm   =  pdif**per
            dpterm1 =  -per*pdif**(per-1.0)
            dpterm2 =  -dpterm1*dpt2

c    porosity, perm, density term
            porterm = mol_gas*porh*den_hyd*afhydr*volh*aterm

c    rate term is made up of individual terms
            prop = porterm*actterm*hterm*wterm*gterm*pterm
c    pressure derivative  (der1)
            der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c    temperature derivative  (der2)
            der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &           + dactterm2*pterm)
c    water fraction derivative  (der3)
            der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &           + wterm*dgterm3)
c    hydrate fraction derivative  (der4)   
            der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &           + hterm*dgterm4)   

         else

c
c    activation energy should be divided by gas constant
c
            if(ihyd_diss_type.eq.1) then 
c Kim-Bishnoi type
               eact2 = e_act1
               den_hyd = 1.0
c    activity term
               actterm   =  exp(-eact2/(var2+temp_conv))
               dactterm2 =  actterm*(eact2/(var2+temp_conv)**2)
            else
c Sakamoto type
               den_hyd = den_hyd0
               afhyd = 1.0
c    activity term
               actterm   =  exp(e_act1/(var2+temp_conv)+e_act2)
               dactterm2 =  -actterm*(e_act1/(var2+temp_conv)**2)
            endif

c    area term (no variable dependence)
            aterm  =  da**oe
c    hydrate fraction term (included hydrate frac from start of eq.)
            if(le+1.0.gt.0.0) then
               hterm   =  var4**(le+1.0)
               dhterm4 =  (le+1.0)*var4**(le)
            else
               hterm = 1.0
               dhterm4 = 0.0
            endif
c    water fraction term
            if(me.gt.0.0) then
               wterm   =  var3**me
               dwterm3 =  me*var3**(me-1.0)
            else
               wterm = 1.0
               dwterm3 = 0.0
            endif
c    gas fraction term
            gvar = 1.0-var3-var4
            if(ne.gt.0.0) then
               gterm   =  gvar**ne
               dgterm3 =  -ne*gvar**(ne-1.0)
               dgterm4 =  dgterm3
            else
               gterm = 1.0
               dgterm3 = 0.0
               dgterm4 = 0.0
            endif
c    pressure (fugacity) term
c    pdif need to be able to change sign
c    might restrict pe to greater than or equal to 1
            pdif   =  pt-var1
            pterm   =  pdif**pe
            dpterm1 =  -pe*pdif**(pe-1.0)
            dpterm2 =  -dpterm1*dpt2
c    activity term

c    porosity, perm, density term
            porterm = mol_gas*porh*den_hyd*afhyd*volh*aterm

c    rate term is made up of individual terms
            prop = porterm*actterm*hterm*wterm*gterm*pterm
c    pressure derivative  (der1)
            der1 = porterm*actterm*hterm*wterm*gterm*dpterm1
c    temperature derivative  (der2)
            der2 = porterm*hterm*wterm*gterm*(actterm*dpterm2
     &           + dactterm2*pterm)
c    water fraction derivative  (der3)
            der3 = porterm*actterm*hterm*pterm*(dwterm3*gterm 
     &           + wterm*dgterm3)
c    hydrate fraction derivative  (der4)   
            der4 = porterm*actterm*wterm*pterm*(dhterm4*gterm 
     &           + hterm*dgterm4)   

         endif
     
      else if(iflg.eq.4.and.idof_meth.ne.7) then

c    heat dissociation rate
c    the following variable are passed through module commeth
c    from last call to hydrate properties (iflg eq 3)
c      source_hyd_temp = skhyd(mi)
c      ds1_tmp = dskhydp
c      ds2_tmp = dskhydt
c      ds3_tmp = dskhydw
c      ds4_tmp = dskhydm
         if (var2 .lt. 0.d0) then
c    low temperature (var2 is temperature)
            hdiss = hdiss01 + (hdiss11*(var2 + temp_conv))
            dhdiss2 = hdiss11
         else
c    higher temperature (var2 is temperature)
            hdiss = hdiss02 + (hdiss12*(var2 + temp_conv))
            dhdiss2 = hdiss12
         endif
         prop = hdiss*source_hyd_temp
         der1 = hdiss*ds1_tmp
         der2 = hdiss*ds2_tmp + dhdiss2*source_hyd_temp
         der3 = hdiss*ds3_tmp
         der4 = hdiss*ds4_tmp
c zero out enthalpy here, put in enthalpy of hydrate
         prop = 0.d0
         der1 = 0.d0 
         der2 = 0.d0
         der3 = 0.d0
         der4 = 0.d0
      elseif (iflg.eq.5) then 
c    calculate temp. dependent dissociation pressure
         prop = exp((ptc1*(var2+temp_conv))+ptc2)
      elseif (iflg.eq.6) then
c    calculate pressure dependent dissociation temperature
         if(var1.ge.0.1) then
            prop = (log(var1)-ptc2)/ptc1 - temp_conv
            der1 = 1.0/(ptc1*var1)
         else
            prop = 0.0
            der1 = 0.0
         endif
      elseif (iflg.eq.8) then
         if(iphase.eq.-1) then
c    Calculate hydrate fraction
c var3 is water frac 
c var4 is gas frac
            prop = 1.0d0/fracv*(1.d0-var3)
	    der3 = -1.0d0/fracv
         else if(iphase.eq.-2) then
            prop = 1.0d0/fracl*(1.d0-var4)
	    der4 = -1.0d0/fracl
         endif
      elseif (iflg.eq.9) then
c    Calculate hydrate fraction
c var3 is water frac 
c var4 is gas frac
         derl1 = 1.0d0/fracv*(1.d0-var3)
         derl2 = 1.0d0/fracl*(1.d0-var4)
	  
         if(derl1.le.derl2) then
c using all the water
            prop =  max(derl1,0.d00)
            der3 = -1.0d0/fracv
            iphase = -1
         else
c using all the gas
            prop =  max(derl2,0.d00)
            der4 = -1.0d0/fracl
            iphase = -2
         endif

      endif
      return 
      end
