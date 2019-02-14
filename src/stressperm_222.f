      subroutine stressperm_222(i)
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
c s kelkar Nov 2011
c combined tensil and shear stress(2 and 22) model
      use comai
      use combi
      use comdi
      use comsi
      implicit none

      integer i,iispmd,ipchk
      integer fail_flag

      real*8 dt,alpv,pary1,pary2,pary3,parsy3,parx1,parx2,parx3,parsx3
      real*8 parxy3,parxz3,paryz3,parz1,parz2,parsz3
      real*8 pp_fac, xpermfac,ypermfac,zpermfac
      real*8 str_x_eff, str_y_eff, str_z_eff
      real*8 xperm_stry,xperm_strz, yperm_strx,yperm_strz
      real*8 zperm_strx,zperm_stry
      real*8 fac(3,3),fac_E(3)
c      real*8 eigenvec(3,3),alambda(3),rm(3,3)
      real*8 rm(3,3)
c     
c     calculate components of volume strain
c     
      fac = 1.
      fac_E=1.

      iispmd = ispm(i)  
      pp_fac=spm3f(iispmd)

      dt = t(i) - tini(i)
      alpv = alp(i)
      vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
c     spm1f222 is strx_min, min x-tensile stress for damage to occur
c     spm2f222 is stry_min, min y-tensile stress for damage to occur
c     spm3f222 is stry_min, min z-tensile stress for damage to occur
c     spm4f222 is e10_facx, damage factor (maximum x) for elastic modulus
c     spm5f222 is e10_facy, damage factor (maximum y) for elastic modulus
c     spm6f222 is e10_facz, damage factor (maximum z) for elastic modulus
c     spm7222f is str_multx, maximum change in x-permeability allowed 
c     spm8f222 is str_multy, maximum change in y-permeability  allowed 
c     spm9f222 is str_multz, maximum change in z-permeability  allowed 
c     
      if(icnl.ne.0) then   
c     2D x-y version  (y can be vertical) 
         str_x_eff=str_x(i)-pp_fac*phi(i)
         str_y_eff=str_y(i)-pp_fac*phi(i)
         strx_min = spm1f222(iispmd)
         stry_min = spm2f222(iispmd)   
         e10_facx = spm4f222(iispmd) 
         e10_facy = spm5f222(iispmd) 
         perx_m = spm7f222(iispmd)
         pery_m = spm8f222(iispmd)
c     
         ipchk = 0
c     changes occur 	
c     x direction changes 
         if(str_y_eff.lt.-stry_min) then
c     pary1 is stress diffrence in tension	 
            pary1 = abs(str_y_eff+stry_min) 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            pary2 = stry_min*(tensile_str_fac-1.)
            pary3 = (perx_m-1.)/pary2      
            xperm_stry = (pary3*pary1+1.0)
c     rock strength never increases	 
            parsy3 =  (e10_facy-1.)/pary2
            e1(i) = min(e1(i),e10(i)*(parsy3*pary1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsy3*pary1+1.0))
            e3(i) = min(e3(i),e10(i)*(parsy3*pary1+1.0))
            e1(i) = max(e1(i), e10_facy*e10(i))
            e2(i) = max(e2(i), e10_facy*e20(i))
            e3(i) = max(e3(i), e10_facy*e30(i))
         endif
         if(str_x_eff.lt.-strx_min) then
c     parx1 is stress diffrence in tension	 
            parx1 = abs(str_x_eff+strx_min) 
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac (=10)
            parx2 = stry_min*(tensile_str_fac-1.)
            parx3 = (pery_m-1.)/parx2
            yperm_strx = (parx3*parx1+1.0)
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
         pnx(i)=pnx0(i)*xperm_stry
         pny(i)=pny0(i)*yperm_strx
     
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
c     
      else if(icnl.eq.0) then   
c     
c     3D  version  (z is always the vertical direction) 
c     
         str_x_eff=str_x(i)-pp_fac*phi(i)
         str_y_eff=str_y(i)-pp_fac*phi(i)
         str_z_eff=str_z(i)-pp_fac*phi(i)
         strx_min = spm1f222(iispmd)
         stry_min = spm2f222(iispmd) 
         strz_min = spm3f222(iispmd) 
         e10_facx = spm4f222(iispmd) 
         e10_facy = spm5f222(iispmd) 
         e10_facz = spm6f222(iispmd)
         perx_m = spm7f222(iispmd)
         pery_m = spm8f222(iispmd) 
         perz_m = spm9f222(iispmd) 
c     
         ipchk = 0
c     
c     x direction tensile stress affects y and z perms 
c     
c     determine maximum changes affecting x direction
         parx1 = 0.
         if(str_x_eff.lt.-strx_min) then
c     parx1 is stress difference in tension	- x direction 
            parx1 = abs((str_x_eff)+strx_min) 
         endif     
         if(parx1.gt.0.0) then	
c     parx1 is stress diffrence in tension calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            parx2 = strx_min*(tensile_str_fac-1.)
            parxy3 = (pery_m-1.)/parx2
            yperm_strx = (parxy3*parx1+1.0)
            parxz3 = (perz_m-1.)/parx2
            zperm_strx = (parxz3*parx1+1.0)	  
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
         if(str_y_eff.lt.-stry_min) then
c     pary1 is stress difference in tension	- y direction 
            pary1 = abs((str_y_eff)+stry_min) 
         endif     
         if(pary1.gt.0.0) then	
c     pary1 is stress diffrence in tension	calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and z dirextions
            pary2 = stry_min*(tensile_str_fac-1.)
            parxy3 = (perx_m-1.)/pary2
            xperm_stry = (parxy3*pary1+1.0)
            paryz3 = (perz_m-1.)/pary2
            zperm_stry = (paryz3*pary1+1.0)	  
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
         if(str_z_eff.lt.-strz_min) then
c     parz1 is stress difference in tension - z direction 
            parz1 = abs((str_z_eff)+strz_min) 
         endif     
         if(parz1.gt.0.0) then	
c     parz1 is stress diffrence in tension	calculated above
c     assume max perm changes occur at tensile stress of 
c     tesile_stess_fac(=10)
c     change x and y dirextions
            parz2 = strz_min*(tensile_str_fac-1.)
            parxz3 = (perx_m-1.)/parz2
            xperm_strz = (parxz3*parz1+1.0)	
            paryz3 = (pery_m-1.)/parz2
            yperm_strz = (paryz3*parz1+1.0)	    
c     rock strength never increases	 
            parsz3 =  (e10_facz-1.)/parz2
            e1(i) = min(e1(i),e10(i)*(parsz3*parz1+1.0))
            e2(i) = min(e2(i),e20(i)*(parsz3*parz1+1.0))
            e3(i) = min(e3(i),e30(i)*(parsz3*parz1+1.0))
            e1(i) = max(e1(i), e10_facz*e10(i))
            e2(i) = max(e2(i), e10_facz*e20(i))
            e3(i) = max(e3(i), e10_facz*e30(i))
         endif

c now do permeability enhanced by shear
         fail_flag = 0
         call stressperm_22_failure(i,fail_flag,rm)
         if(fail_flag.eq.1) then
            if(itemp_perm22(i).eq.0) then
               write(91,*)i
               itemp_perm22(i) = 1
            endif
            call stressperm_22_perm(i,rm, fac)
            call stressperm_22_emod(i,rm,fac_E)
         endif

         xpermfac = max(xperm_stry,xperm_strz, fac(1,1))
         ypermfac = max(yperm_strx,yperm_strz, fac(2,2))
         zpermfac = max(zperm_stry,zperm_strx, fac(3,3))

         pnx(i) = pnx0(i)*xpermfac
         pny(i) = pny0(i)*ypermfac
         pnz(i) = pnz0(i)*zpermfac
         
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
         pnz(i) = min(perx_m*pnz0(i),pnz(i))	 
         
         e1(i)=e10(i)*fac_E(1)
         e2(i)=e20(i)*fac_E(1)
         e3(i)=e30(i)*fac_E(1)
      
      endif
      
      return

      end
c.....................................................................
