      subroutine stressperm_222(jpt)

c s kelkar 4/29/2010 simillar to ispmd=2 except
c using effective stress instead of the total stress
c
c     only good for 2D or 3D models
c     
      use comai
      use combi
      use comdi
      use comsi

      implicit none
      integer jpt,ipchk,iispmd
      real*8 pary1,pary2,pary3,parsy3,parx1,parx2,parx3
      real*8 parsx3,parxy3,parxz3,paryz3,parz1,parz2,parz3,parsz3
      real*8 dt,alpv,e1i,e2i,e3i,efac,epi,eti,dpd,biot,ctherm,shti
      real*8 shpi

      iispmd = ispm(jpt)    
               
      biot=bulk(jpt)
      ctherm=alp(jpt)
      e1i = e1(jpt)
      e2i = e2(jpt)
      e3i = e3(jpt)
      efac=3.d0*e2i+2.d0*e3i
c     stress due to temp and pore pressure changes
      epi=efac*biot
      eti=efac*ctherm
      dpd=phi(jpt)-phini(jpt)
      dt=t(jpt)-tini(jpt)
      shti=(eti*dt)
      shpi=(epi*dpd)
      
      vol_temp(jpt) = alpv*dt*sx1(jpt) - vol_strain(jpt)*sx1(jpt) 
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
         pary1 = (str_y(jpt)-shpi-shti)+stry_min  
         if(pary1.lt.0.) then
c     pary1 is stress diffrence in tension	 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            pary1 = abs(pary1) 
            pary2 = stry_min*(tensile_str_fac-1.)
            pary3 = (perx_m-1.)/pary2      
            pnx(jpt) = pnx0(jpt)*(pary3*pary1+1.0)
c     rock strength never increases	 
            parsy3 =  (e10_facy-1.)/pary2
            e1(jpt) = min(e1(jpt),e10(jpt)*(parsy3*pary1+1.0))
            e2(jpt) = min(e2(jpt),e20(jpt)*(parsy3*pary1+1.0))
            e3(jpt) = min(e3(jpt),e10(jpt)*(parsy3*pary1+1.0))
            e1(jpt) = max(e1(jpt), e10_facy*e10(jpt))
            e2(jpt) = max(e2(jpt), e10_facy*e20(jpt))
            e3(jpt) = max(e3(jpt), e10_facy*e30(jpt))
         endif
         parx1 = (str_x(jpt)-shpi-shti)+strx_min
         if(parx1.lt.0.) then
c     parx1 is stress diffrence in tension	 
            parx1 = abs(parx1) 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            parx2 = stry_min*(tensile_str_fac-1.)
            parx3 = (pery_m-1.)/parx2
            pny(jpt) = pny0(jpt)*(parx3*parx1+1.0)
c     rock strength never increases	 
            parsx3 =  (e10_facx-1.)/parx2
            e1(jpt) = min(e1(jpt),e10(jpt)*(parsx3*parx1+1.0))
            e2(jpt) = min(e2(jpt),e20(jpt)*(parsx3*parx1+1.0))
            e3(jpt) = min(e3(jpt),e30(jpt)*(parsx3*parx1+1.0))
            e1(jpt) = max(e1(jpt), e10_facx*e10(jpt))
            e2(jpt) = max(e2(jpt), e10_facx*e20(jpt))
            e3(jpt) = max(e3(jpt), e10_facx*e30(jpt))
         endif
c     
         pnx(jpt) = min(perx_m*pnx0(jpt),pnx(jpt))
         pny(jpt) = min(pery_m*pny0(jpt),pny(jpt))
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
c     
c     x direction tensile stress affects y and z perms    
c     determine maximum changes affecting x direction
         parx1 = (str_x(jpt)-shpi-shti)+strx_min
c     parx1 is effective stress difference in tension: x direction 
         if(parx1.lt.0.0) then	
c     parx1 is stress diffrence in tension calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            parx1 = abs(parx1) 
            parx2 = strx_min*(tensile_str_fac-1.)
            parxy3 = (pery_m-1.)/parx2
            pny(jpt) = pny0(jpt)*(parxy3*parx1+1.0)
            parxz3 = (perz_m-1.)/parx2
            pnz(jpt) = pnz0(jpt)*(parxz3*parx1+1.0)	  
c     rock strength never increases	 
            parsx3 =  (e10_facx-1.)/parx2
            e1(jpt) = min(e1(jpt),e10(jpt)*(parsx3*parx1+1.0))
            e2(jpt) = min(e2(jpt),e20(jpt)*(parsx3*parx1+1.0))
            e3(jpt) = min(e3(jpt),e30(jpt)*(parsx3*parx1+1.0))
            e1(jpt) = max(e1(jpt), e10_facx*e10(jpt))
            e2(jpt) = max(e2(jpt), e10_facx*e20(jpt))
            e3(jpt) = max(e3(jpt), e10_facx*e30(jpt))
         endif	 
c     
c     y direction tensile stress affects x and z perms 
c     
c     determine maximum changes affecting y direction
         pary1 = (str_y(jpt)-shpi-shti)+stry_min
         if(pary1.lt.0.0) then
c     pary1 is stress difference in tension: y direction 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            pary1 = abs(pary1) 
            pary2 = stry_min*(tensile_str_fac-1.)
            parxy3 = (perx_m-1.)/pary2
            pnx(jpt) = pnx0(jpt)*(parxy3*pary1+1.0)
            paryz3 = (perz_m-1.)/pary2
            pnz(jpt) = pnz0(jpt)*(paryz3*pary1+1.0)	  
c     rock strength never increases	 
            parsy3 =  (e10_facy-1.)/pary2
            e1(jpt) = min(e1(jpt),e10(jpt)*(parsy3*pary1+1.0))
            e2(jpt) = min(e2(jpt),e20(jpt)*(parsy3*pary1+1.0))
            e3(jpt) = min(e3(jpt),e30(jpt)*(parsy3*pary1+1.0))
            e1(jpt) = max(e1(jpt), e10_facy*e10(jpt))
            e2(jpt) = max(e2(jpt), e10_facy*e20(jpt))
            e3(jpt) = max(e3(jpt), e10_facy*e30(jpt))
         endif
c     
c     z direction tensile stress affects x and y perms 
c     
c     determine maximum changes affecting z direction
         parz1 = (str_z(jpt)-shti-shpi)+strz_min
         if(parz1.lt.0.0) then
c     parz1 is stress difference in tension : z direction 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and y dirextions
            parz1 = abs(parz1) 
            parz2 = strz_min*(tensile_str_fac-1.)
            parxz3 = (perx_m-1.)/parz2
            pnx(jpt) = pnx0(jpt)*(parxz3*parz1+1.0)	
            paryz3 = (pery_m-1.)/parz2
            pny(jpt) = pny0(jpt)*(paryz3*parz1+1.0)	    
c     rock strength never increases	 
            parsz3 =  (e10_facz-1.)/parz2
            e1(jpt) = min(e1(jpt),e10(jpt)*(parsz3*parz1+1.0))
            e2(jpt) = min(e2(jpt),e20(jpt)*(parsz3*parz1+1.0))
            e3(jpt) = min(e3(jpt),e30(jpt)*(parsz3*parz1+1.0))
            e1(jpt) = max(e1(jpt), e10_facz*e10(jpt))
            e2(jpt) = max(e2(jpt), e10_facz*e20(jpt))
            e3(jpt) = max(e3(jpt), e10_facz*e30(jpt))
         endif
         
         pnx(jpt) = min(perx_m*pnx0(jpt),pnx(jpt))
         pny(jpt) = min(pery_m*pny0(jpt),pny(jpt))
         pnz(jpt) = min(perx_m*pnz0(jpt),pnz(jpt))	 
         
      endif
      
      return
      end
      
c.....................................................................
