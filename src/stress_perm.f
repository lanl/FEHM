      subroutine stress_perm(iflg,ndummy)         
c
c directional_perm_information
c max and min directional quantities
c    
c      use comai
c      use combi
c      use comdti
c	use comji
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
	integer iflg,i,kb,kc,i1,i2,jj,jj1,ndummy,ipchk,ipr_str
	integer ikbxmin,ikbxmax,ikbymin,ikbymax,ikbzmin,ikbzmax
	integer ikcxmin,ikcxmax,ikcymin,ikcymax,ikczmin,ikczmax
	integer kbxmin,kbxmax,kbymin,kbymax,kbzmin,kbzmax
	integer kbx1,kbx2,kby1,kby2,kbz1,kbz2
	parameter (cperm=1.0)
c
c
      if(istrs.eq.0) return
      if(ipermstr.eq.0) return
c
      if(iflg.eq.0) then
c
      if(.not.allocated(ipermx)) then
	 allocate(ipermx(n0,2))
 	 allocate(ipermy(n0,2))
	 allocate(ipermz(n0,2))
	   ipermx = 0
	   ipermy = 0
	   ipermz = 0
	endif
c
	   do i = 1,neq
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
	   enddo

      else if (iflg.eq.-1) then
c
c  allocate memory for stress derivatives
c
  
       allocate(rlxs(n0))
 	 allocate(rlys(n0))
	 allocate(rlzs(n0))
	 allocate(drlxs(n0,4))
	 allocate(drlys(n0,4))
	 allocate(drlzs(n0,4))
	  
	 allocate (idum_str1(n0))
	 idum_str1 = 0 
      else if (iflg.eq.-2) then
       deallocate(rlxs,rlys,rlzs)
       deallocate(drlxs,drlys,drlzs)
	 deallocate(idum_str1)
	else if (iflg.eq.1) then
	if(ipermstr.eq.3) then
c perm model 1 
c
c  perm changes for individual node
c  assumes call has been made to allocate memory
c      
      do i = 1,neq
c  x direction (orthogonal pieces)
      kbx1 = ipermx(i,1)
      kbx2 = ipermx(i,2)
      kby1 = ipermy(i,1)
      kby2 = ipermy(i,2)
	kbz1 = ipermz(i,1)
      kbz2 = ipermz(i,2)

      disx = cord(kbx2,1)-cord(kbx1,1)
      disy = cord(kby2,2)-cord(kby1,2)
	disz = cord(kbz2,3)-cord(kbz1,3)

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
      enddo
	else if(ipermstr.eq.1) then
c perm model 1 - not implemented yet (default)
      if(allocated(rlxs)) then
       rlxs = 1.0
	 rlys = 1.0 
	 rlzs = 1.0
       drlxs = 0.0
	 drlys = 0.0
	 drlzs = 0.0
	endif

	else if(ipermstr.eq.2) then
c perm model 2 - volume strains
c lagged permeability (only after each time step)
c only good for 2-D models
      if(icnl.eq.0) then
        if(iout.ne.0) then
         write(iout,*) '>>> perm model 2 not valid for 3D(stopping)'
        endif
        if(iptty.ne.0) then
         write(iptty,*) '>>>perm model 2 not valid for 3D(stopping)'
        endif
        stop
      endif
      if(.not.allocated(pnx0)) then
	 allocate(pnx0(neq))
	 allocate(pny0(neq))
	 allocate(e10(neq))
	 allocate(e20(neq))
	 allocate(e30(neq))	 
	 pnx0 = pnx
	 pny0 = pny
c pnz is a counter here	 
	 pnz = 1.e20
	 e10 = e1
	 e20 = e2
	 e30 = e3
	endif
	ipr_str = 0
      do i = 1,neq
c      
c calculate components of volume strain
c      
	 dt = t(i) - tini(i)
	 alpv = alp(i)
c	 stry_min = 0.05
c	 strx_min = 0.05
c	 e10_fac = 0.1
	 vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
	 ipchk = 0
	 if(str_y(i).lt.-stry_min) then
        kx_fac = abs((str_y(i)/stry_min))
	  kx_fac = min(str_mult,kx_fac)
	  pnx(i) = max(pnx0(i)*kx_fac,pnx(i))
	  e1(i) = e10(i)*e10_fac
	  e2(i) = e20(i)*e10_fac
	  e3(i) = e30(i)*e10_fac	
	  if(pnz(i).gt.3.) pnz(i) = 1. 
	  if(pnz(i).eq.2) pnz(i) = 3. 	  
	  if(pnz(i).le.3.)ipchk = 1 
	 endif
	 if(str_x(i).lt.-strx_min) then
        ky_fac = abs((str_x(i)/strx_min))
	  ky_fac = min(str_mult,ky_fac)
	  pny(i) = max(pny0(i)*ky_fac,pny(i))
	  e1(i) = e10(i)*e10_fac
	  e2(i) = e20(i)*e10_fac
	  e3(i) = e30(i)*e10_fac	 
	  if(pnz(i).gt.3.) pnz(i) = 2. 	  
	  if(pnz(i).eq.1.) pnz(i)=3.
	  if(pnz(i).le.3.)ipchk = 1 	  
	 endif
        if(ipchk.ne.0) ipr_str = ipr_str +1
	enddo
	 write(*,100) l,days,ipr_str
	  write(iout,100) l,days,ipr_str
    
      	else if(ipermstr.eq.4) then
c perm model 3 - keita model
c lagged permeability (only after each time step)
c only good for 2-D models
      if(icnl.eq.0) then
        if(iout.ne.0) then
         write(iout,*) '>>> perm model 4 not valid for 3D(stopping)'
        endif
        if(iptty.ne.0) then
         write(iptty,*) '>>>perm model 4 not valid for 3D(stopping)'
        endif
        stop
      endif
      if(.not.allocated(pnx0)) then
	 allocate(pnx0(neq))
	 allocate(pny0(neq))
	 allocate(e10(neq))
	 allocate(e20(neq))
	 pnx0 = pnx
	 pny0 = pny
c pnz is a counter here	 
	 pnz = 1.e20
	 e10 = 1.0
	 e20 = 1.0
	endif
	ipr_str = 0
      do i = 1,neq
c      
c calculate components of volume strain
c      
	 dt = t(i) - tini(i)
	 alpv = alp(i)
	 vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
	 ipchk = 0
	 if(str_y(i).lt.0.0) then
        kx_fac = 1. + stry_min*abs(str_y(i)**3)
	  kx_fac = min(str_mult,kx_fac)
	  pnx(i) = max(pnx0(i)*kx_fac,e10(i)*pnx0(i))
	  if(kx_fac.gt.e10_fac) e10(i) = e10_fac
	  if(pnz(i).gt.3.) pnz(i) = 1. 
	  if(pnz(i).eq.2) pnz(i) = 3. 	  
	  if(pnz(i).le.3.)ipchk = 1 
	 endif
	 if(str_x(i).lt.0.0) then
        ky_fac = 1. + strx_min*abs(str_x(i)**3)
	  ky_fac = min(str_mult,ky_fac)
	  pny(i) = max(pny0(i)*ky_fac,e20(i)*pny0(i))
	  if(kx_fac.gt.e10_fac) e20(i) = e10_fac
	  if(pnz(i).gt.3.) pnz(i) = 2. 	  
	  if(pnz(i).eq.1.) pnz(i)=3.
	  if(pnz(i).le.3.)ipchk = 1 	  
	 endif
        if(ipchk.ne.0) ipr_str = ipr_str +1
	enddo
	 write(*,100) l,days,ipr_str
	  write(iout,100) l,days,ipr_str
      endif
	
	endif
100   format(/,'time step', 1x,i6,' days'1x,f9.2,' damaged nodes ',i8) 	
	return 
	end
	 
	 
	 
	 
	 