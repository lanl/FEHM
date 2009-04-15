      subroutine stress_perm(iflg,ndummy)         
c
c directional_perm_information
c max and min directional quantities
c    
c      use comai
c      use combi
c      use comdti
c	 use comji
c      use comsi 
	
	use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comji
      use comki
      use comxi
      use davidi
      use comsi    
      implicit none 
	real*8 xdmin,xdmax,ydmin,ydmax,zdmin,zdmax,xi,yi,zi
	real*8 xkb,ykb,zkb,xd1,yd1,zd1,dis,disx,disy,disz
	real*8 dt,alpv,kx_fac,ky_fac 
 	real*8 drlsu1,drlsu2,drlsv1,drlsv2,drlsw1,drlsw2,cperm  
      real*8 pary1,pary2,pary3,parsy3
	real*8 parx1,parx2,parx3,parsx3
	real*8 parz1,parz2,parz3,parsz3
	real*8  paryz1, paryz2, paryz3, parsyz3, stryz_min
	real*8  parxz1, parxz2, parxz3, parsxz3, strxz_min
	real*8  parxy1, parxy2, parxy3, parsxy3, strxy_min
	real*8 permsum,permsum0,permchng,permchng_tol
	real*8  permx_max, permy_max, permz_max
      real*8  sxx_min, syy_min, szz_min
      real*8  coorxz_max,cooryz_max,coorzz_max
 
	integer ipermx_max, ipermy_max, ipermz_max
	integer iflg,i,kb,kc,i1,i2,jj,jj1,ndummy,ipchk,ipr_str
	integer ikbxmin,ikbxmax,ikbymin,ikbymax,ikbzmin,ikbzmax
	integer ikcxmin,ikcxmax,ikcymin,ikcymax,ikczmin,ikczmax
	integer kbxmin,kbxmax,kbymin,kbymax,kbzmin,kbzmax
	integer kbx1,kbx2,kby1,kby2,kbz1,kbz2,idir,iispmd
	parameter (cperm=1.0)
	parameter(permchng_tol = 0.01)
c 
c
      if(istrs.eq.0) return
      if(ipermstr.eq.0) return
c
      if(iflg.eq.0) then
c      
c model 2 and model 4 require an initial setup      
c
c
      if(ipermstr2.ne.0) then
c allocate space for parameters      
       if(.not.allocated(pnx0)) then
	  allocate(pnx0(neq))
	  allocate(pny0(neq))
	  allocate(pnz0(neq))
	  allocate(e10(neq))
	  allocate(e20(neq))
	  allocate(e30(neq))	 
	  pnx0 = pnx
	  pny0 = pny
	  pnz0 = pnz
	  e10 = e1
	  e20 = e2
	  e30 = e3
	 endif

	endif
c      
c model 3 and model 5 require an initial setup      
c
c only allocate if there is a model 3 and model 5
c	
      if(ipermstr3.ne.0) then
       if(.not.allocated(ipermx)) then
	  allocate(ipermx(n0,2))
 	  allocate(ipermy(n0,2))
	  allocate(ipermz(n0,2))
	   ipermx = 0
	   ipermy = 0
	   ipermz = 0
	 endif  
c
c only calculate for model 3 and model 5
c initial setup calcs node neighbor information
c
	   do i = 1,n0
	    iispmd = ispm(i)
	    ispmd = ispmt(iispmd)
	    if(ispmd.eq.3.or.ispmd.eq.5) then
	     xi = cord(i,1)
	     yi = cord(i,2)
	     zi = cord(i,3)
	     xdmax = 1.d30
	     xdmin = 1.d30
	     ydmax = 1.d30
	     ydmin = 1.d30
	     zdmax = 1.d30
	     zdmin = 1.d30
           kbxmax = i
           kbxmin = i
           kbymax = i
           kbymin = i
           kbzmax = i
           kbzmin = i
	     i1 = nelm(i)+1
	     i2 = nelm(i+1)
	     do jj = i1,i2
	      kb = nelm(jj)
  	       xkb = cord(kb,1)
  	       ykb = cord(kb,2)
	       zkb = cord(kb,3)
               if(xkb-xi.gt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.xdmax) then
                   kbxmax = kb
                   xdmax = dis
	           endif
               else if(xkb-xi.lt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.xdmin) then
                   kbxmin = kb
                   xdmin = dis
	           endif
	         endif
               if(ykb-yi.gt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.ydmax) then
                   kbymax = kb
                   ydmax = dis
	           endif
               else if(ykb-yi.lt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.ydmin) then
                   kbymin = kb
                   ydmin = dis
	           endif
	         endif 
               if(zkb-zi.gt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.zdmax) then
                   kbzmax = kb
                   zdmax = dis
	           endif
               else if(zkb-zi.lt.0) then
	          dis = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
	           if(dis.lt.zdmin) then
                   kbzmin = kb
                   zdmin = dis
	           endif
	          endif
		  enddo                              
	    
c
c  store max and min of each direction (may create larger stencil)
c	     
	     ipermx(i,1) = kbxmin
	     ipermx(i,2) = kbxmax
	     ipermy(i,1) = kbymin
	     ipermy(i,2) = kbymax	     
	     ipermz(i,1) = kbzmin
	     ipermz(i,2) = kbzmax	  
	    endif 	  
	   enddo
       	endif
      endif
      if (iflg.eq.-1) then
c
c  allocate memory for stress derivatives for fully coupled solution
c  just before call to generate equations
c  
c
       if(ihms.gt.0)then  
        allocate(rlxs(n0))
 	  allocate(rlys(n0))
	  allocate(rlzs(n0))
	  allocate(drlxs(n0,4))
	  allocate(drlys(n0,4))
	  allocate(drlzs(n0,4))	  
	  allocate (idum_str1(n0))
	  idum_str1 = 0 
	 endif
      else if (iflg.eq.-2) then
c
c  deallocate memory for stress derivatives for fully coupled solution
c  
c
       if(ihms.gt.0)then       
        deallocate(rlxs,rlys,rlzs)
        deallocate(drlxs,drlys,drlzs)
	  deallocate(idum_str1)
	 endif
c	 
	else if (iflg.eq.1) then
c 
c  This is the general loop to calculate permeability models and derivatives
c	
c skip loop if permmacro in read
      if(.not.allocated(ispm)) return
c  
      do i = 1,n0
c
c identify model associated with node i
c		    
      iispmd = ispm(i)    
	ispmd = ispmt(iispmd) 
c
	if(ispmd.eq.1.and.ihms.eq.1) then
c perm model 1 - not implemented yet (default)
c not needed unless fully coupled
       if(allocated(rlxs)) then
        rlxs(i) = 1.0
	  rlys(i) = 1.0 
	  rlzs(i) = 1.0
        drlxs(i,1) = 0.0
	  drlys(i,1) = 0.0
	  drlzs(i,1) = 0.0
	  drlxs(i,2) = 0.0
	  drlys(i,2) = 0.0
	  drlzs(i,2) = 0.0
        drlxs(i,3) = 0.0
	  drlys(i,3) = 0.0
	  drlzs(i,3) = 0.0
	  drlxs(i,4) = 0.0
	  drlys(i,4) = 0.0
	  drlzs(i,4) = 0.0	  
	 endif	 
	else if(ispmd.eq.3) then
c	
c perm model 3 fully coupled model
c
c  perm changes for individual node
c  assumes call has been made to allocate memory
c      
c  x direction (orthogonal pieces)
c
      kbx1 = ipermx(i,1)
      kbx2 = ipermx(i,2)
      kby1 = ipermy(i,1)
      kby2 = ipermy(i,2)
	kbz1 = ipermz(i,1)
      kbz2 = ipermz(i,2)
c
      disx = cord(kbx2,1)-cord(kbx1,1)
      disy = cord(kby2,2)-cord(kby1,2)
	disz = cord(kbz2,3)-cord(kbz1,3)
c
      rlxs(i) = cperm*(1. + (dv(kby2)-dv(kby1))/disy)**2*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2
      drlxs(i,2) = 2*cperm*(1. + (dv(kby2)-dv(kby1))/disy)*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2/disy
      drlxs(i,1) = -2*cperm*(1. + (dv(kby2)-dv(kby1))/disy)*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2/disy
      drlxs(i,4) = 2*cperm*(1. + (dw(kbz2)-dw(kbz1))/disz)*
     &          (1. + (dv(kby2)-dv(kby1))/disy)**2/disz
      drlxs(i,3) = -2*cperm*(1. + (dw(kbz2)-dw(kbz1))/disz)*
     &          (1. + (dv(kby2)-dv(kby1))/disy)**2/disz

      rlys(i) = cperm*(1. + (du(kbx2)-du(kbx1))/disx)**2*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2
      drlys(i,2) = 2*cperm*(1. + (du(kbx2)-du(kbx1))/disx)*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2/disx
      drlys(i,1) = -2*cperm*(1. + (du(kbx2)-du(kbx1))/disx)*
     &          (1. + (dw(kbz2)-dw(kbz1))/disz)**2/disx
      drlys(i,4) = 2*cperm*(1. + (dw(kbz2)-dw(kbz1))/disz)*
     &          (1. + (du(kbx2)-du(kbx1))/disx)**2/disz
      drlys(i,3) = -2*cperm*(1. + (dw(kbz2)-dw(kbz1))/disz)*
     &          (1. + (du(kbx2)-du(kbx1))/disx)**2/disz

      rlzs(i) = cperm*(1. + (dv(kby2)-dv(kby1))/disy)**2*
     &          (1. + (du(kbx2)-du(kbx1))/disx)**2
      drlzs(i,2) = 2*cperm*(1. + (du(kbx2)-du(kbx1))/disx)*
     &          (1. + (dv(kby2)-dv(kby1))/disy)**2/disx
      drlzs(i,1) = -2*cperm*(1. + (du(kbx2)-du(kbx1))/disx)*
     &          (1. + (dv(kby2)-dv(kby1))/disy)**2/disx
      drlzs(i,4) = 2*cperm*(1. + (dv(kby2)-dv(kby1))/disy)*
     &          (1. + (du(kbx2)-du(kbx1))/disx)**2/disy
      drlzs(i,3) = -2*cperm*(1. + (dv(kby2)-dv(kby1))/disy)*
     &          (1. + (du(kbx2)-du(kbx1))/disx)**2/disy


	else if(ispmd.eq.2) then
c perm model 2 - volume strains
c lagged permeability only (stress-based) (only after each time step)
c only good for 2D or 3D models
c      
c calculate components of volume strain
c   
      
	 dt = t(i) - tini(i)
	 alpv = alp(i)
	 vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
c spm1f is strx_min, min tensile stress (x direction) for damage to occur
c spm2f is stry_min, min tensile stress (y direction) for damage to occur
c spm3f is stry_min, min tensile stress (z direction) for damage to occur
c spm4f is e10_facx, damage factor (maximum x) for elastic modulus
c spm5f is e10_facy, damage factor (maximum y) for elastic modulus
c spm6f is e10_facz, damage factor (maximum z) for elastic modulus
c spm7f is str_multx, maximum change in permeability (x direction) allowed 
c spm8f is str_multy, maximum change in permeability (y direction) allowed 
c spm9f is str_multz, maximum change in permeability (z direction) allowed 
c model 3 and model 5 are fully coupled
c model 6 is simple directional plasticity
c
       if(icnl.ne.0) then   
c 2D x-y version  (y can be vertical) 
       strx_min = spm1f(iispmd)
       stry_min = spm2f(iispmd)   
       e10_facx = spm4f(iispmd) 
       e10_facy = spm5f(iispmd) 
       perx_m = spm7f(iispmd)
       pery_m = spm8f(iispmd)
c              
	 ipchk = 0
c changes occur 	
c x direction changes 
	 if(str_y(i).lt.-stry_min) then
c  pary1 is stress diffrence in tension	 
        pary1 = abs((str_y(i))+stry_min) 
c assume max perm changes occur at tensile stress of tesile_stess_fac (=10)
        pary2 = stry_min*(tensile_str_fac-1.)
        pary3 = (perx_m-1.)/pary2      
	  pnx(i) = pnx0(i)*(pary3*pary1+1.0)
c rock strength never increases	 
        parsy3 =  (e10_facy-1.)/pary2
	  e1(i) = min(e1(i),e10(i)*(parsy3*pary1+1.0))
	  e2(i) = min(e2(i),e20(i)*(parsy3*pary1+1.0))
	  e3(i) = min(e3(i),e10(i)*(parsy3*pary1+1.0))
	  e1(i) = max(e1(i), e10_facy*e10(i))
	  e2(i) = max(e2(i), e10_facy*e20(i))
	  e3(i) = max(e3(i), e10_facy*e30(i))
	 endif
	 if(str_x(i).lt.-strx_min) then
c  parx1 is stress diffrence in tension	 
        parx1 = abs((str_x(i))+strx_min) 
c assume max perm changes occur at tensile stress of tesile_stess_fac (=10)
        parx2 = stry_min*(tensile_str_fac-1.)
        parx3 = (pery_m-1.)/parx2
	  pny(i) = pny0(i)*(parx3*parx1+1.0)
c rock strength never increases	 
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
c 3D  version  (z is always the vertical direction) 
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
c x direction tensile stress affects y and z perms 
c 
c determine maximum changes affecting x direction
       parx1 = 0.
	 if(str_x(i).lt.-strx_min) then
c  parx1 is stress difference in tension	- x direction 
        parx1 = abs((str_x(i))+strx_min) 
       endif     
       if(parx1.gt.0.0) then	
c parx1 is stress diffrence in tension calculated above
c assume max perm changes occur at tensile stress of tesile_stess_fac(=10)
c change x and z dirextions
        parx2 = strx_min*(tensile_str_fac-1.)
        parxy3 = (pery_m-1.)/parx2
	  pny(i) = pny0(i)*(parxy3*parx1+1.0)
        parxz3 = (perz_m-1.)/parx2
	  pnz(i) = pnz0(i)*(parxz3*parx1+1.0)	  
c rock strength never increases	 
        parsx3 =  (e10_facx-1.)/parx2
	  e1(i) = min(e1(i),e10(i)*(parsx3*parx1+1.0))
	  e2(i) = min(e2(i),e20(i)*(parsx3*parx1+1.0))
	  e3(i) = min(e3(i),e30(i)*(parsx3*parx1+1.0))
	  e1(i) = max(e1(i), e10_facx*e10(i))
	  e2(i) = max(e2(i), e10_facx*e20(i))
	  e3(i) = max(e3(i), e10_facx*e30(i))
	 endif	 
c 	
c y direction tensile stress affects x and z perms 
c 
c determine maximum changes affecting y direction
       pary1 = 0.
	 if(str_y(i).lt.-stry_min) then
c  pary1 is stress difference in tension	- y direction 
        pary1 = abs((str_y(i))+stry_min) 
       endif     
       if(pary1.gt.0.0) then	
c pary1 is stress diffrence in tension	calculated above
c assume max perm changes occur at tensile stress of tesile_stess_fac(=10)
c change x and z dirextions
        pary2 = stry_min*(tensile_str_fac-1.)
        parxy3 = (perx_m-1.)/pary2
	  pnx(i) = pnx0(i)*(parxy3*pary1+1.0)
        paryz3 = (perz_m-1.)/pary2
	  pnz(i) = pnz0(i)*(paryz3*pary1+1.0)	  
c rock strength never increases	 
        parsy3 =  (e10_facy-1.)/pary2
	  e1(i) = min(e1(i),e10(i)*(parsy3*pary1+1.0))
	  e2(i) = min(e2(i),e20(i)*(parsy3*pary1+1.0))
	  e3(i) = min(e3(i),e30(i)*(parsy3*pary1+1.0))
	  e1(i) = max(e1(i), e10_facy*e10(i))
	  e2(i) = max(e2(i), e10_facy*e20(i))
	  e3(i) = max(e3(i), e10_facy*e30(i))
	 endif
c 	
c z direction tensile stress affects x and y perms 
c 
c determine maximum changes affecting z direction
       parz1 = 0.0
	 if(str_z(i).lt.-strz_min) then
c  parz1 is stress difference in tension - z direction 
        parz1 = abs((str_z(i))+strz_min) 
       endif     
       if(parz1.gt.0.0) then	
c parz1 is stress diffrence in tension	calculated above
c assume max perm changes occur at tensile stress of tesile_stess_fac(=10)
c change x and y dirextions
        parz2 = strz_min*(tensile_str_fac-1.)
        parxz3 = (perx_m-1.)/parz2
	  pnx(i) = pnx0(i)*(parxz3*parz1+1.0)	
        paryz3 = (pery_m-1.)/parz2
	  pny(i) = pny0(i)*(paryz3*parz1+1.0)	    
c rock strength never increases	 
        parsz3 =  (e10_facz-1.)/parz2
	  e1(i) = min(e1(i),e10(i)*(parsz3*parz1+1.0))
	  e2(i) = min(e2(i),e20(i)*(parsz3*parz1+1.0))
	  e3(i) = min(e3(i),e30(i)*(parsz3*parz1+1.0))
	  e1(i) = max(e1(i), e10_facz*e10(i))
	  e2(i) = max(e2(i), e10_facz*e20(i))
	  e3(i) = max(e3(i), e10_facz*e30(i))
	 endif

	   pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))
         pnz(i) = min(perx_m*pnz0(i),pnz(i))	 

	 endif
	 
      	else if(ipermstr.eq.4) then
c perm model 4 - keita model
c differs from perm model 2 by allowing cubic variation 
c of perm with stress  
c lagged permeability (only after each time step)
c only good for 2-D models
c      
c calculate components of volume strain
c      
	 dt = t(i) - tini(i)
	 alpv = alp(i)
	 vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
	 ipchk = 0

c        kx_fac = 1. + stry_min*abs(str_y(i)**3)

        
	 endif
      enddo
c      
c check for damage zone (permeability changes)
c and maximum allowable changes
c  
       if(ipermstr2.ne.0) then
        permx_max = 0.0
        permy_max = 0.0
        permz_max = 0.0
        ipermx_max = 1
        ipermy_max = 1
        ipermz_max = 1
        sxx_min = 0.0
        syy_min = 0.0
        szz_min = 0.0
        ipr_str = 0 
        do i = 1,n0
        if(icnl.eq.0) then
         permsum = pnx(i)+pny(i)+pnz(i)
         permsum0 = pnx0(i)+pny0(i)+pnz0(i)
        else
         pnx(i) = min(perx_m*pnx0(i),pnx(i))
         pny(i) = min(pery_m*pny0(i),pny(i))        
         permsum = pnx(i)+pny(i)
         permsum0 = pnx0(i)+pny0(i)       
        endif
         permchng = (permsum - permsum0)/permsum0
         if(permchng.gt.permchng_tol) then
c count damaged zone nodes
           ipr_str = ipr_str +1
c find minimum stresses in the damaged zone           
           if(str_x(i).lt.sxx_min) then
             sxx_min  = str_x(i) 
             ipermx_max = i
           endif
           if(str_y(i).lt.syy_min) then
             syy_min  = str_y(i) 
             ipermy_max = i
           endif  
           if(icnl.eq.0.and.str_z(i).lt.szz_min) then
             szz_min  = str_z(i) 
             ipermz_max = i
           endif                 
         endif
        enddo
c        
c write out information for the damaged zones    
c
        if(icnl.eq.0) then
         coorxz_max = cord(ipermx_max,3)
         cooryz_max = cord(ipermy_max,3)
         coorzz_max = cord(ipermz_max,3)
         else
         coorxz_max = 0.0
         cooryz_max = 0.0
         coorzz_max = 0.0        
        endif
        if(iout.ne.0) write(iout,99) l,days 
	  if(iptty.ne.0) write(iptty,99) l,days  
        if(iout.ne.0) write(iout,100) ipr_str 
	  if(iptty.ne.0) write(iptty,100) ipr_str
        if(iout.ne.0) write(iout,101) 
	  if(iptty.ne.0) write(iptty,101)   
      if(iout.ne.0)
     &  write(iout,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &       cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max 
	if(iptty.ne.0)
     &  write(iptty,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &       cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max 
      if(iout.ne.0) 
     &  write(iout,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &       cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
	if(iptty.ne.0)
     & write(iptty,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &       cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
        if(icnl.eq.0) then
      if(iout.ne.0)
     & write(iout,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &       cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max  
	if(iptty.ne.0)
     & write(iptty,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &       cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max          
        endif     
       endif         
      endif
99    format(/,1x'Time step ',i6,' Days'1x,f9.2)   
100   format(1x,'Number of damaged gridbolcks', 1x,i6)
101   format(1x,'Largest tensile stresses in damaged zone')
102   format(1x,'Node ',i6,1x,'Sxx ',g12.4,1x,'Kxx ',
     &       g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
103   format(1x,'Node ',i6,1x,'Syy ',g12.4,1x,'Kyy ',
     &       g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
104   format(1x,'Node ',i6,1x,'Szz ',g12.4,1x,'Kzz ',
     &       g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)            
	return 
	end
	 
	 
	 
	 
	 