      subroutine rlperm_co2(ndummy,iz,i,prop1,
     &           dprop11,dprop12,prop2,dprop21,dprop22,
     &           prop3,dprop31,dprop32)
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

C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To compute the relative permeabilities for CO2 problem
CD1 Two phases present if single phase CO2 
CD1 Three phases present if two phase CO2 
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 04-04-2007   R. Pawar   N/A     Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/rlperm_hyd.f_a  $
CD2 
C**********************************************************************
c
c calculates relative permiabilities and derivatives
c
      use comhi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comki
      use comki
      use comco2
      implicit none

      integer iz,ndummy,i,irlpd,mi,it,ir,j,num_models,ireg,icesd,iflg
      real*8 prop1,dprop11,dprop12,prop2,dprop21,dprop22,alpha,beta
      real*8 prop3,dprop31,dprop32,rp1,rp2,rp3,rp4,rp5,rp6
      real*8 rp7,rp8,rp9,rp10,rp11,rp111,rp12,rp13,rp14,krow,krog
      real*8 dkrow1,dkrow2,dkrog1,dkrog2,sw,sl,sg,denom
      real*8 rp31,rp32,sw_star,dsw_stardw,sc_star,sl_star
      real*8 dsl_stardw,dsl_stardg,dsc_stardw,dsc_stardg,dsw_stardg
      real*8 slope,dhp,ds,hcut,hmin,hmax,slcut,smcut,fac,fac_use
      real*8 fac_min,sucut,smcutm,rp15,rp16
      real*8 s_hat,s_star,b1,b2,cr,d,ddds,dcds,sl1,sl2,tol_l,tol_u,rl
      real*8 drls,rv,drvs,hp,smr,slr
      parameter(hmin=0.d0)
      parameter(fac_min=2.d0)
      
      if(idof_co2.le.2) then
         if(iz.eq.0) then
            if(idof_co2.eq.1) then
               prop1 = 1
               dprop11 = 0
               dprop12 = 0
            else
               if(ices(i).ne.2) then
                  prop2 = 1
                  dprop21 = 0
                  dprop22 = 0
               else
                  prop2 = fl(i)
                  dprop21 = -1.d0
                  dprop22 = -1.d0
                  prop3 = fg(i)
                  dprop31 = 0.d0
                  dprop32 = 1.d0
               endif
            endif
         else
            prop1 = 0
            prop2 = 0
            dprop11 = 0
            dprop12 = 0
            dprop21 = 0
            dprop22 = 0
         endif
      else
         mi = i+ndummy
         icesd = ices(mi)       
         it = irlp(mi)
         if(it.eq.0) then
            irpd=0
         else
            irpd = irlpt(it)
         endif
c     rp1 minimum water saturation (connate water)
c     rp2 max water saturation
c     rp3 water rel perm curve exponent (1 equals linear)
c     rp4 minimun CO2 (l/s single phase) saturation, it is assumed to be l in 3 phase calculation
c     rp5 max CO2 (l/s single phase) saturation
c     rp6 CO2 (l/s single phase)-water rel perm curve exponent
c     rp7 CO2-l CO2-g rel perm curve exponent
c     rp8 minimun CO2 (gas phase) saturation
c     rp9 max CO2 (gas phase) saturation
c     rp10 CO2 (gas phase) rel perm curve exponent
         rp1 = rp1f(it)
         rp2 = rp2f(it)
         rp3 = rp3f(it)
         rp4 = rp4f(it)
         rp5 = rp5f(it)
         rp6 = rp6f(it)
         rp7 = rp7f(it)
         rp8 = rp8f(it)
         rp9 = rp9f(it)
         rp10 = rp10f(it)
c     rp11 exponent for co2(l)-water capillary pressures (if rp11 is -ve
c     Brooks-Corey capillary pressure relationship is used) then rp11
c     is exponent lambda.
c     rp12 max co2(l)-water capillary pressures (Entry pressure for 
c     Brooks-Corey)
c     rp13 exponent for co2(l)-co2(g) capillary pressures
c     rp14 max co2(l)-co2(g) capillary pressures
c     rp15 is the low saturation fitting parameter for Brooks Corey
c     rp16 is the cutoff saturation used in fits described for rp15 must be 
c     greater than rp1
         rp11 = rp11f(it)
         rp12 = rp12f(it)
         rp13 = rp13f(it)
         rp14 = rp14f(it)
         rp15 = rp15f(it)
         rp16 = rp16f(it)
c     Perform a check on connate water saturation before proceeding
c     Max liquid and gas saturation should be 1-rp1, if not overwrite specified values
         if(rp5.ne.(1-rp1)) rp5=1-rp1
         if(rp9.ne.(1-rp1)) rp9=1-rp1		
         sw = fw(mi)
         sl = fl(mi)
         sg = fg(mi)
         if(icesd.eq.3) sl = fg(mi)
         prop1 = 0.0
         prop2 = 0.0
         prop3 = 0.d0
         dprop11 = 0.d0
         dprop12 = 0.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
         dprop31 = 0.d0
         dprop32 = 0.d0
         

c     
c     calculate relative perm and derivatives
c     note: must be called for each node i
c     
c     
         if(irpd.eq.19) then
c test for co2 phase
c this formulation is not general for 3phase (water/co2-gas/co2-liquid)
	 if(sg>sl) sl=sg
c     assumes liquid water and single phase-CO2
c     van Genuchten for water and Corey for CO2 (regardless of phase)
c     equivalent to Pruess 2002; requires 6 parameters 
c     rp1 = irreducible water saturation (for rel perm)
c     rp2 = irreducable CO2 saturation 
c     rp3 = exponent (lambda), equivalent to Van Genuchten (1980) m
c     rp4 = irreducable water saturation (for cap pressure))
c     rp5 = lower cuttoff for linear interpolation  (for rel perm)
c     rp6 = upper cuttoff for  linear interpolation  (for rel perm)
c     rp7 = Strength coefficient (Mpa)
c note:  rp5 and rp6 are applied to s-star for rel perm interpolation
c and to saturation for capillary pressure interpolation

            rp1=rp1f(it)
            rp2=rp2f(it)
            rp3=rp3f(it)
            rp4=rp4f(it)
            rp5=rp5f(it)
            rp6=rp6f(it)
            rp7=rp7f(it)
            b1=1.D0/(1.d0-rp1)
            s_star=(sw-rp1)*b1
            prop3=0.
            dprop31=0.
            dprop32=0.
            alpha=1./rp7
c alpha is 1/Po
            beta=1.d0/(1.d0-rp3)
c beta is Van Genuchten's n; rp3 is Van Genuchten's m = Preuss's lambda
            if(iz.eq.0) then
c we're calculating rel perms

c first we calculate the water part
               iflg=1
c	star=s_star
               sl1=rp1
               sl2=1.
c sl2 is maximum liquid saturation
               tol_l=0.
               tol_u=1.
c rel perm calculations for water and single-phase CO2
               call vgrlps(iflg, sw, sl1,sl2, rp3, rp5,
     2              rp6, rl,drls, rv, drvs )
               prop1=rl
               dprop11=drls
               dprop12=0.

c		cr=1-s_star**(1.d0/rp3)
c		d=1-cr**rp3
c		prop1=sqrt(s_star)*(1.D0-d)**2.D0
c		dcds=(-1.d0/rp3)*s_star**(1.d0/rp3-1.d0)
c		ddds=-rp3*cr**(rp3-1.d0)*dcds
c		dprop11=b1*(0.5d0*s_star**-0.5d0*(1.d0-d)**2+
c   +	sqrt(s_star)*2.d0*(1.d0-d)*ddds)
c	dprop12=0.
c
c  now we're doing the co2 part
               if(sl.ge.rp2) then
c co2 properties, if we're above min co2 sat
                  b2=1.D0/(1.d0-rp1-rp2)
                  s_hat=min(1.d0,(sw-rp1)*b2)
                  prop2=(1.d0-s_hat)**2.d0*(1.d0-s_hat**2.d0)
                  if (prop2 .gt. 1.0d0) then
                     prop2 = 1.d0
                     dprop21 = 0.d0
                  else if (prop2 .lt. 0.0d0) then
                     prop2 = 0.d0
                     dprop21 = 0.d0
                  else
                     dprop21=(2.0*(1.d0-s_hat)*(1.d0-s_hat**2.d0) +
     &                    (1.d0-s_hat)**2.d0*2.0*(-s_hat))*(-b2)
c                  dprop21=-b2
                  end if
                  dprop22=0.
               else
c co2 properties, if we're at or below min co2 sat
                  prop2=0.
                  dprop21=0.
                  dprop22=0.
               endif
            else
c     cap pressure
c these 2 are set according to Preuss(2004)
               slr=0.;smr=1.
c these are flexible
               smcut=rp5;sucut=1.0-rp6
               call vgcap_ek	( sw, slr, smr,alpha, beta,smcut, 
     2              sucut,hp, dhp)
      
     
     
               prop1=hp
               dprop11=dhp
               dprop12=0.d0
c end of iz =1 or iz = 0
            endif
            return
c done wtih model 19
         endif

         if(iz.eq.0) then
c     CO2-water 2-phase
c     This is always between sw and sl
c     two-phase mixture
            if(irpd.eq.18) then
c     Brooks-Corey relationship
c     rp3 is Brooks-Corey exponent
c     rp1 is residual water saturation
c     rp4 is residual CO2 liquid/gas saturation for 2-phase, residual gas saturation for 3-phase
c     rp7 is residual CO2 liquid saturation for 3-phase
               rp31=((2.d0+3.d0*rp3)/rp3)
               rp32=(2.d0+rp3)/rp3
               sw_star = (sw-rp1)/(1.d0-rp1)
               dsw_stardw = 1.d0/(1.d0-rp1)
               if(sw_star.le.0.d0) then
                  prop1 = 0.d0
                  dprop11 = 0.d0
                  dprop12 = 0.d0
               elseif(sw_star.ge.1.d0) then
                  prop1 = 1.d0
                  dprop11 = 0.d0
                  dprop12 = 0.d0
               else
                  prop1 = sw_star**rp31
                  dprop11=rp31*((sw-rp1)**(rp31-1))/((1-rp1)**rp31)
               endif

               if(icesd.eq.2) then
                  sl_star = (sw+sl-rp1)/(1.d0-rp1)
c     dsl_stardw = 1.d0/(1.d0-rp1)
                  dsl_stardw = 0.d0
                  dsl_stardg = -1.d0/(1.d0-rp1)
                  sc_star = (sl-rp7)/(1-rp1-rp7)
                  dsc_stardw = -1.d0/(1-rp1-rp7)
                  dsc_stardg = -1.d0/(1-rp1-rp7)
                  dsw_stardg = 0.d0
                  if(sl_star.le.0.d0) then
                     prop3 = 1.d0
                     dprop31 = 0.d0
                     dprop32 = 0.d0
                  elseif(sl_star.ge.1.d0) then
                     prop3 = 0.d0
                     dprop31 = 0.d0
                     dprop32 = 0.d0
                  else
                     prop3=((1.d0-sl_star)**2.d0)*(1.d0-(sl_star**rp32))
                     dprop31 = (-2.d0*(1.d0-sl_star)*
     &                    (1.d0-(sl_star**rp32))+((1.d0-sl_star)**2.d0)*
     &                    (-rp32*sl_star**(rp32-1.d0)))*(dsl_stardw)
                     dprop32 = (-2.d0*(1.d0-sl_star)*
     &                    (1.d0-(sl_star**rp32))+((1.d0-sl_star)**2.d0)*
     &                    (-rp32*sl_star**(rp32-1.d0)))*(dsl_stardg)
                  endif
                  if(sc_star.le.0.d0) then
                     prop2 = 1.d0
                     dprop21 = 0.d0
                     dprop22 = 0.d0
                  elseif(sc_star.ge.1.d0) then
                     prop2 = 0.d0
                     dprop21 = 0.d0
                     dprop22 = 0.d0
                  else
                     prop2=(sc_star**2.d0)*((sl_star**rp32)-
     &                    (sw_star**rp32))
                     dprop21=((sl_star**rp32)-(sw_star**rp32))*
     &                    2.d0*sc_star*dsc_stardw+(sc_star**2.d0)*rp32*
     &                    (sl_star**(rp32-1.d0)*dsl_stardw-
     &                    sc_star**(rp32-1.d0)*dsc_stardw)
                     dprop22=((sl_star**rp32)-(sw_star**rp32))*
     &                    2.d0*sc_star*dsc_stardg+(sc_star**2.d0)*rp32*
     &                    (sl_star**(rp32-1.d0)*dsl_stardg-
     &                    sc_star**(rp32-1.d0)*dsc_stardg)
                  endif
               else
                  if(sw_star.le.0.d0) then
                     prop2 = 1.d0
                     dprop21 = 0.d0
                     dprop22 = 0.d0
                  elseif(sw_star.ge.1.d0) then
                     prop2 = 0.d0
                     dprop21 = 0.d0
                     dprop22 = 0.d0
                  else
                     prop2=((1.d0-sw_star)**2.d0)*
     &                    (1.d0-(sw_star**rp32))
                     dprop21 = (-2.d0*(1.d0-sw_star)*
     &                    (1.d0-(sw_star**rp32))+((1.d0-sw_star)**2.d0)*
     &                    (-rp32*sw_star**(rp32-1.d0)))*(dsw_stardw)
                  endif
               endif
            else
               if(rp3.eq.1) then
c     linear function of saturations
                  if(sw.ge.rp1) then
                     prop1 = 1.d0
                     dprop11 = 0.d0
                     if(sw.le.rp2) then
                        denom = rp2-rp1
                        prop1=(sw-rp1)/denom
                        dprop11 = 1.d0/denom
                     endif
                  else
                     prop1 = 0.d0
                     dprop11 = 0.d0
                  endif
               else
c     Corey relationship
                  prop1 = ((sw-rp1)/(1.d0-rp1))**rp3
                  dprop11 = rp3*((sw-rp1)**(rp3-1))/((1-rp1)**rp3)
c     if(icesd.eq.2) dprop12=-rp3*((sw-rp1)**(rp3-1))/
c     &		((1-rp1)**rp3)
                  if ((sw-rp1).eq.0.d0) then
                     dprop11 = 0.d0
                     dprop12 = 0.d0
                  endif
               endif
c     
               if(rp6.eq.1) then
                  if(sl.ge.rp4) then
                     krow = 1.d0
                     prop2 = krow
                     dprop21 = 0.d0
                     if(sl.le.rp5) then
                        denom=rp5-rp4
                        krow=(sl-rp4)/denom
                        dkrow1= -1.d0/denom
                        if(icesd.eq.2) then
                           dkrow2=-1.d0/denom
                        else
                           prop2 = krow
                           dprop21 = dkrow1
                        endif
                     endif
                  else
                     prop2 = 0.d0
                     dprop21 = 0.d0
                  endif
               else
                  krow= ((sl-rp4)/(1.d0-rp4))**rp6
                  dkrow1 = rp6*((sl-rp4)**(rp6-1))/((1-rp4)**rp6)
                  if((sl-rp4).eq.0.d0) then
                     dkrow1 = 0.d0
                     dkrow2 = 0.d0
                  endif
                  if(icesd.eq.2) then
                     dkrow2 =-rp6*((sl-rp4)**(rp6-1))/
     &                    ((1-rp4)**rp6)
                  else
                     prop2 = krow
                     dprop21 = dkrow1
                  endif
               endif
c     
               if(icesd.eq.2) then
                  if(rp7.eq.1) then
                     krog = krow
                     dkrog1 = dkrow1
                     dkrog2 = dkrow2
                  else
                     krog= ((sl-rp4)/(1.d0-rp4))**rp7
                     dkrog1 = rp7*((sl-rp4)**(rp7-1))/((1-rp4)**rp7)
                     if(icesd.eq.2) dkrog2 =-rp7*((sl-rp4)**(rp7-1))/
     &                    ((1-rp4)**rp7)
                     if((sl-rp4).eq.0.d0) then
                        dkrog1 = 0.d0
                        dkrog2 = 0.d0
                     endif
                  endif
c     
                  if(rp10.eq.1) then
                     if(sg.gt.rp8) then
                        prop3 = 1.d0
                        if(sg.le.rp9) then
                           denom = rp9-rp8
                           prop3=(sg-rp8)/denom
c     dprop31=-1.d0/denom
                           dprop32 = 1.d0/denom
                        endif
                     endif
                  else
                     prop3 = ((sg-rp8)/(1.d0-rp8))**rp10
c     dprop31 = -rp10*((sg-rp8)**(rp10-1))/((1-rp8)**rp10)
                     dprop32 = rp10*((sg-rp8)**(rp10-1))/((1-rp8)**rp10)
                     if ((sg-rp8).eq.0.d0) then
                        dprop31 = 0.d0
                        dprop32 = 0.d0
                     endif
                  endif
c     modify prop2 using the default 3-phase formulation in Eclipse
                 if(irpd.eq.17) then
                     prop2 = sl
                     dprop21 = -1.d0
                     dprop22 = -1.d0
                  else
                     prop2 = (sg*krog+(sw-rp1)*krow)/(sg+sw-rp1)
                     dprop21 = sg*krow/(sg+sw-rp1)**2.d0
                     dprop22 = (sw-rp1)*krog/(sg+sw-rp1)**2.d0
                  endif
               endif
            endif
*** end of model 18 or not ****************


         else
c     capillary pressures
c     rp11 exponent for co2(g/l/sc)-water capillary pressures
c     rp12 max co2(l)-water capillary pressures
c     rp13 exponent for co2(l)-co2(g) capillary pressures
c     rp14 max co2(l)-co2(g) capillary pressures
            prop1 = 0.d0
            dprop11 = 0.d0
            dprop12 = 0.d0
            prop2 = 0.d0
            dprop21 = 0.d0
            dprop22 = 0.d0
            prop3 = 0.d0
            dprop31 = 0.d0
            dprop32 = 0.d0
            fac=rp15
            if(irpd.eq.18) then
c     Use Brooks-Corey capillary pressure relationship
c     rp11 is exponent, rp12 is entry pressure
c     sucut = 1.d0
               smcut=(rp16-rp1)/(rp2-rp1)
               smcutm=max(smcut,1.d-3)
               slcut=smcutm*(rp2-rp1)+rp1
               if(sw.gt.0.d0) then
                  if(sw.le.rp16) then
                     hcut=rp12*(smcutm**(-1.d0/rp11))
                     if(rp15.gt.0.0) then
c     a multiple of value of
                        fac_use=max(rp15,fac_min)
                        hmax=hcut*fac
                        slope=-(hmax-hcut)/slcut
                     else
c     linear fit from saturation curve
                        ds=1.d0/(rp2-rp1)
                        dhp=rp12*(-1.d0/rp11)*
     &                       smcutm**((-rp11-1.d0)/rp11)*ds
                        hmax=hcut-dhp*smcutm
                        slope=-(hmax-hcut)/smcutm
                     endif
c     prop1=slope*((sw-rp1)/(rp2-rp1))+hmax
                     prop1=slope*sw+hmax
                     dprop11=slope
c     elseif(sw.ge.sucut) then
c     slcut=(sucut-rp1)/(rp2-rp1)
c     hcut=rp12*(slcut**(-1.d0/rp11))
c     slope=(hcut-hmin)/(sucut-1.d0)
c     prop1=slope*sw+hmin-slope
c     dprop11=slope
                  else
                     prop1 = rp12*((sw-rp1)/(rp2-rp1))**(-1.d0/rp11)
                     dprop11 = rp12*(-1.d0/rp11)*
     &                    ((sw-rp1)**((-rp11-1.d0)/rp11))/
     &                    ((rp2-rp1)**(-1.d0/rp11))
                  endif
               elseif(sw.le.0.d0) then
                  hcut=rp12*(smcutm**(-1.d0/rp11))
                  if(rp15.gt.0.0) then
c     a multiple of value of
                     fac_use=max(rp15,fac_min)
                     hmax=hcut*fac
                  else
c     linear fit from saturation curve
                     ds=1.d0/(rp2-rp1)
                     dhp=rp12*(-1.d0/rp11)*
     &                    smcutm**((-rp11-1.d0)/rp11)*ds
                     hmax=hcut-dhp*smcutm
                  endif
                  prop1=hmax
                  dprop11=0.d0
c     else
c     slope=(hcut-hmin)/(1.d0-slcut)
c     prop1=-slope*((sw-rp1)/(rp2-rp1))+slope-hmin
c     dprop11=-slope
c     slcut=(sucut-rp1)/(rp2-rp1)
c     hcut=rp12*(slcut**(-1.d0/rp11))
c     slope=(hcut-hmin)/(sucut-1.d0)
c     prop1=slope*sw+hmin-slope
c     dprop11=slope

c     prop1=0.d0
c     dprop11=0.d0
               endif
            else					
               if(sw.le.rp1) then
                  prop1 = rp12
                  dprop11 = 0.d0
                  dprop12 = 0.d0
               elseif(sw.le.rp2) then
                  if(rp11.eq.1) then
                     prop1 = rp12*(rp2-sw)/(rp2-rp1)
                     dprop11 = -rp12/(rp2-rp1)
                     dprop12 = 0.d0
                  else
                     prop1 = rp12*(rp2-sw)**rp11/(rp2-rp1)
                     dprop11 = -rp12*rp11*((rp2-sw)**(rp11-1))/(rp2-rp1)
                     dprop12 = 0.d0
                  endif
               else
                  prop1 = 0.d0
                  dprop11 = 0.d0
                  dprop12 = 0.d0
               endif
            endif
            if(icesd.eq.2) then
               if(sg.le.rp8) then
                  prop2 = 0.d0
                  dprop21 = 0.d0
                  dprop22 = 0.d0
               elseif(sg.le.rp9) then
                  if(rp13.eq.1) then
                     prop2 = rp14*(sg-rp8)/(rp9-rp8)
                     dprop21 = 0.d0
                     dprop22 = rp14/(rp9-rp8)
                  else
                     prop2 = rp14*(sg-rp8)**rp13/(rp9-rp8)
                     dprop21 = 0.0
                     dprop22 = rp14*rp13*((sg-rp8)**(rp13-1))/(rp9-rp8)
                  endif
               else
                  prop2 = rp14
                  dprop21 = 0.d0
                  dprop22 = 0.d0
               endif
            endif
         endif
      endif
      return
      end
