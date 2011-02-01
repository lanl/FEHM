      subroutine stressperm_6(i)
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

c     same as model 5 but added pore pressure term 
c     
c     perm model 2 - volume strains
c     lagged permeability only (stress-based) (only after each time step)
c     only good for 2D or 3D models
c     
c     need to subtract biot and thermal stresses to get 
c     lithostatic component
      use comai
      use combi
      use comdi
      use comsi
      
      implicit none
      integer i, iispmd,ipchk

      real*8 biot,ctherm,e1i,e2i,e3i,erat,efac,epi,eti,dpd,dt
      real*8 shti,shpi,lith_damage_fac,str_eff,pary1,pary2
      real*8 pary3,parsy3,parx1,parx2,parx3,parsx3,stress_min_prin
      real*8 lith_min,pterm,parxy3,parxz3,paryz3,parz1,parz2
      real*8 parsz3,parz3

      parameter (lith_min = 1.0)      
c     
c     
      
      biot=bulk(i)
      ctherm=alp(i)
      e1i = e1(i)
      e2i = e2(i)
      e3i = e3(i)
      erat=e2i/e1i
      efac=3.d0*e2i+2.d0*e3i
c     stress due to temp and pore pressure changes
      epi=efac*biot
      eti=efac*ctherm
      dpd=phi(i)-phini(i)
      dt=t(i)-tini(i)
      shti=(eti*dt)
      shpi=(epi*dpd)
      
c     spm1f is strx_min, min x-tensile stressfor damage to occur
c     spm2f is stry_min, min y-tensile stressfor damage to occur
c     spm3f is stry_min, min z-tensile stressfor damage to occur
c     spm4f is e10_facx, damage factor (maximum x) for elastic modulus
c     spm5f is e10_facy, damage factor (maximum y) for elastic modulus
c     spm6f is e10_facz, damage factor (maximum z) for elastic modulus
c     spm7f is str_multx, maximum change in x-permeability allowed 
c     spm8f is str_multy, maximum change in y-permeability  allowed 
c     spm9f is str_multz, maximum change in z-permeability  allowed 
c     model 3 and model 5 are fully coupled
c     model 6 is simple directional plasticity
c     
      if(icnl.ne.0) then   
c     2D x-y version  (y can be vertical) 
         strx_min = spm1f(iispmd)
         stry_min = spm2f(iispmd)   
         e10_facx = spm4f(iispmd) 
         e10_facy = spm5f(iispmd) 
         perx_m = spm7f(iispmd)
         pery_m = spm8f(iispmd)
         lith_damage_fac = spm10f(iispmd)
         str_eff = spm11f(iispmd)
c     
         ipchk = 0
c     changes occur 	
c     x direction changes 
         if(str_y(i).lt.-stry_min) then
c     pary1 is stress diffrence in tension	 
            pary1 = abs((str_y(i))+stry_min) 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            pary2 = stry_min*(tensile_str_fac-1.)
            pary3 = (perx_m-1.)/pary2      
            pnx(i) = pnx0(i)*(pary3*pary1+1.0)
c     rock strength never increases	 
            parsy3 =  (e10_facy-1.)/pary2
            e1(i) = min(e1(i),e10(i)*(parsy3*pary1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsy3*pary1+1.0))
            e3(i) = min(e3(i),e10(i)*(parsy3*pary1+1.0))
            e1(i) = max(e1(i), e10_facy*e10(i))
            e2(i) = max(e2(i), e10_facy*e20(i))
            e3(i) = max(e3(i), e10_facy*e30(i))
         endif
         if(str_x(i).lt.-strx_min) then
c     parx1 is stress diffrence in tension	 
            parx1 = abs((str_x(i))+strx_min) 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            parx2 = stry_min*(tensile_str_fac-1.)
            parx3 = (pery_m-1.)/parx2
            pny(i) = pny0(i)*(parx3*parx1+1.0)
c     rock strength never increases	 
            parsx3 =  (e10_facx-1.)/parx2
            e1(i) = min(e1(i),e10(i)*(parsx3*parx1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsx3*parx1+1.0))
            e3(i) = min(e3(i),e30(i)*(parsx3*parx1+1.0))
            e1(i) = max(e1(i), e10_facx*e10(i))
            e2(i) = max(e2(i), e10_facx*e20(i))
            e3(i) = max(e3(i), e10_facx*e30(i))
         endif
c     
c     pore pressure term (greater than the minimum earth stress)
c     only effects y direction (assumed the maximum stress direction)
c     
         if(lith_damage_fac.gt.0.0) then
            stress_min_prin = 
     &           max((str_y(i))*lith_damage_fac,lith_min)
            pterm = phi(i)-pres_elev
            pary1 = (pterm-stress_min_prin)/stress_min_prin
            if(pary1.gt.0.0) then
               pary3 = (pery_m-1.)/str_eff 
               pny(i) = pny0(i)*(pary3*pary1+1.0)
            endif 
         endif
         
c     
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
c     
      else if(icnl.eq.0) then   
c     
c     3D  version  (z is always the vertical direction) 
c     
         strx_min = spm1f(iispmd)
         stry_min = spm2f(iispmd) 
         strz_min = spm3f(iispmd) 
         e10_facx = spm4f(iispmd) 
         e10_facy = spm5f(iispmd) 
         e10_facz = spm6f(iispmd)
         perx_m = spm7f(iispmd)
         pery_m = spm8f(iispmd) 
         perz_m = spm9f(iispmd) 
         lith_damage_fac = spm10f(iispmd)
         str_eff = spm11f(iispmd)
c     
         ipchk = 0
c     
c     x direction tensile stress affects y and z perms 
c     
c     determine maximum changes affecting x direction
         parx1 = 0.
         if(str_x(i).lt.-strx_min) then
c     parx1 is stress difference in tension	- x direction 
            parx1 = abs((str_x(i))+strx_min) 
         endif     
         if(parx1.gt.0.0) then	
c     parx1 is stress diffrence in tension calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            parx2 = strx_min*(tensile_str_fac-1.)
            parxy3 = (pery_m-1.)/parx2
            pny(i) = pny0(i)*(parxy3*parx1+1.0)
            parxz3 = (perz_m-1.)/parx2
            pnz(i) = pnz0(i)*(parxz3*parx1+1.0)	  
c     rock strength never increases	 
            parsx3 =  (e10_facx-1.)/parx2
            e1(i) = min(e1(i),e10(i)*(parsx3*parx1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsx3*parx1+1.0))
            e3(i) = min(e3(i),e30(i)*(parsx3*parx1+1.0))
            e1(i) = max(e1(i), e10_facx*e10(i))
            e2(i) = max(e2(i), e10_facx*e20(i))
            e3(i) = max(e3(i), e10_facx*e30(i))
         endif	 
c     
c     y direction tensile stress affects x and z perms 
c     
c     determine maximum changes affecting y direction
         pary1 = 0.
         if(str_y(i).lt.-stry_min) then
c     pary1 is stress difference in tension	- y direction 
            pary1 = abs((str_y(i))+stry_min) 
         endif     
         if(pary1.gt.0.0) then	
c     pary1 is stress diffrence in tension	calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            pary2 = stry_min*(tensile_str_fac-1.)
            parxy3 = (perx_m-1.)/pary2
            pnx(i) = pnx0(i)*(parxy3*pary1+1.0)
            paryz3 = (perz_m-1.)/pary2
            pnz(i) = pnz0(i)*(paryz3*pary1+1.0)	  
c     rock strength never increases	 
            parsy3 =  (e10_facy-1.)/pary2
            e1(i) = min(e1(i),e10(i)*(parsy3*pary1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsy3*pary1+1.0))
            e3(i) = min(e3(i),e30(i)*(parsy3*pary1+1.0))
            e1(i) = max(e1(i), e10_facy*e10(i))
            e2(i) = max(e2(i), e10_facy*e20(i))
            e3(i) = max(e3(i), e10_facy*e30(i))
         endif
c     
c     z direction tensile stress affects x and y perms 
c     
c     determine maximum changes affecting z direction
         parz1 = 0.0
         if(str_z(i).lt.-strz_min) then
c     parz1 is stress difference in tension - z direction 
            parz1 = abs((str_z(i))+strz_min) 
         endif     
         if(parz1.gt.0.0) then	
c     parz1 is stress diffrence in tension	calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and y dirextions
            parz2 = strz_min*(tensile_str_fac-1.)
            parxz3 = (perx_m-1.)/parz2
            pnx(i) = pnx0(i)*(parxz3*parz1+1.0)	
            paryz3 = (pery_m-1.)/parz2
            pny(i) = pny0(i)*(paryz3*parz1+1.0)	    
c     rock strength never increases	 
            parsz3 =  (e10_facz-1.)/parz2
            e1(i) = min(e1(i),e10(i)*(parsz3*parz1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsz3*parz1+1.0))
            e3(i) = min(e3(i),e30(i)*(parsz3*parz1+1.0))
            e1(i) = max(e1(i), e10_facz*e10(i))
            e2(i) = max(e2(i), e10_facz*e20(i))
            e3(i) = max(e3(i), e10_facz*e30(i))
         endif
c     
c     effective stress permeability enhancement
c     assume z direction is the maximum lithostatic stress direction
c     
         if(lith_damage_fac.gt.0.0) then
            stress_min_prin = 
     &           max((str_z(i))*lith_damage_fac,lith_min)
            pterm = phi(i)-pres_elev
            parz1 = (pterm-stress_min_prin)/stress_min_prin
            if(parz1.gt.0.0) then
               parz3 = (perz_m-1.)/str_eff 
               pnz(i) = pnz0(i)*(parz3*parz1+1.0)
            endif 
         endif 
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
         pnz(i) = min(perx_m*pnz0(i),pnz(i))	 
         
      endif   	
      
      return
      
      end
c....................................................................
