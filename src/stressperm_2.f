      subroutine stressperm_2(i)
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

c     perm model 2 - volume strains
c     lagged permeability only(stress-based,only after each time step)
c     only good for 2D or 3D models
      use comai
      use combi
      use comdi
      use comsi
      implicit none

      integer i,iispmd,ipchk

      real*8 dt,alpv,pary1,pary2,pary3,parsy3,parx1,parx2,parx3,parsx3
      real*8 parxy3,parxz3,paryz3,parz1,parz2,parsz3
      real*8 p_eff, perm_fac, str_max, fac
      
c     
c     calculate components of volume strain
c     
      str_max = 1.
      iispmd = ispm(i)       
      dt = t(i) - tini(i)
      p_eff = phi(i)-phini(i)
      alpv = alp(i)
      vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
c     spm1f is strx_min, min x-tensile stress for damage to occur
c     spm2f is stry_min, min y-tensile stress for damage to occur
c     spm3f is stry_min, min z-tensile stress for damage to occur
c     spm4f is e10_facx, damage factor (maximum x) for elastic modulus
c     spm5f is e10_facy, damage factor (maximum y) for elastic modulus
c     spm6f is e10_facz, damage factor (maximum z) for elastic modulus
c     spm7f is str_multx, maximum change in x-permeability allowed 
c     spm8f is str_multy, maximum change in y-permeability  allowed 
c     spm9f is str_multz, maximum change in z-permeability  allowed 
c     
      if(icnl.ne.0) then   
c     2D x-y version  (y can be vertical) 
         strx_min = spm1f(iispmd)
         stry_min = spm2f(iispmd)   
         e10_facx = spm4f(iispmd) 
         e10_facy = spm5f(iispmd) 
         perx_m = spm7f(iispmd)
         pery_m = spm8f(iispmd)
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
c     
         ipchk = 0
c s kelkar 5 june 2012. If strx_min <or= 0 model is interpreted to be that
c for local failure driven by pore pressure only. 
c  If strx_min>0 , it is taken to be the usual model 2
         if(strx_min.gt.0.) then     
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
         else
            fac = p_eff + strx_min
            if(fac.gt.0.0) then
               fac=fac/str_max
               perm_fac = fac*(perx_m - 1.0d0) + 1.0d0
               pnx(i) = pnx0(i)*perm_fac
               pny(i) = pny0(i)*perm_fac
               pnz(i) = pnz0(i)*perm_fac
            endif
            
         endif
         
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
         pnz(i) = min(perx_m*pnz0(i),pnz(i))	 
         
      endif
      
      return

      end
c.....................................................................
