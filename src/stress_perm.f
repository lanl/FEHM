      subroutine stress_perm(iflg,ndummy)         
c     
c     directional_perm_information
c     max and min directional quantities
c     
      
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
      real*8 xdmint,xdmaxt,ydmint,ydmaxt,zdmint,zdmaxt
      real*8 xkb,ykb,zkb,xd1,yd1,zd1,dis,dist,disx,disy,disz
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
      real*8 frac_bx, frac_by, frac_bz
      real*8 amultx, amulty, amultz
      real*8 amultxy, amultxz, amultyz
      real*8 dispx12, dispy12, dispz12
      real*8 disptx12, dispty12, disptz12
      real*8 dispxy12, dispxz12, dispyx12, dispyz12, dispzx12, dispzy12
      real*8 alpkbx1, alpkbx2, alpkby1, alpkby2, alpkbz1, alpkbz2, alpi
      real*8 tx12, ty12, tz12
      real*8 biot,erat,efac,epi,dpd,shpi,stress_min_prin,lith_min
      real*8 ctherm,eti,shti,e1i,e2i,e3i,lith_damage_fac,str_eff,pterm
      
c****************************local parameters used in perm model 4      
      real*8  ipmd4_p1,ipmd4_p2,ipmd4_p3,ipmd4_p4
      real*8  ipmd4_k
      
c****************************local parameters used in perm model 7 Bai 
      real*8  ipmd7_dsx,ipmd7_dsy,ipmd7_dsz
      real*8  ipmd7_dsxy,ipmd7_dsxz,ipmd7_dsyz
      real*8  ipmd7_knx,ipmd7_kny,ipmd7_knz
      real*8  ipmd7_nmx,ipmd7_nmy,ipmd7_nmz
      real*8  ipmd7_shx,ipmd7_shy,ipmd7_shz
      real*8  pi,frac_tol
c****************************local parameters used in perm model 8  Bai 
      real*8  ipmd8_dsx,ipmd8_dsy,ipmd8_dsz
      real*8  ipmd8_dsxy
      real*8  ipmd8_knx,ipmd8_kny,ipmd8_knz
      real*8  ipmd8_nmx,ipmd8_nmy,ipmd8_nmz
      real*8  ipmd8_shx,ipmd8_shy,ipmd8_shz
      
      integer ipermx_max, ipermy_max, ipermz_max
      integer iflg,i,kb,kc,i1,i2,jj,jj1,ndummy,ipchk,ipr_str
      integer ikbxmin,ikbxmax,ikbymin,ikbymax,ikbzmin,ikbzmax
      integer ikcxmin,ikcxmax,ikcymin,ikcymax,ikczmin,ikczmax
      integer kbxmin,kbxmax,kbymin,kbymax,kbzmin,kbzmax
      integer kbx1,kbx2,kby1,kby2,kbz1,kbz2,idir,iispmd
c................
      real*8 gk0,gpmod,gmexp,gn,sigy_eff
      real*8 dis_tol
c.........................................
      parameter (cperm=1.0,lith_min = 1.0,pi = 3.14159)
      parameter (permchng_tol = 0.01,frac_tol=0.00001)
      parameter (dis_tol = 1.d-8)
c     
c     
      if(idoff .eq. -1) return      
      if(istrs.eq.0) return
      if(ipermstr.eq.0) return
c     
      if(iflg.eq.0) then
c     
c     model 2 and model 4 require an initial setup      
c     
c     
         if(ipermstr2.ne.0.or.ipermstr6.ne.0) then
c     allocate space for parameters      
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
            if(ipermstr6.ne.0) then
c     
c     effective stress damage model - need lowest pressure 
c     at highest elevation  
c     
               elev_high = -1.e30
               do i = 1,n0
                  if(cord(i,igrav).gt.elev_high) then
                     elev_high = cord(i,igrav)
                     pres_elev = pho(i)
                  endif
               enddo  
               pres_elev = pres_elev-0.1   
            endif
         endif	 
         
         if(ipermstr7.ne.0.or.ipermstr8.ne.0.or.ipermstr5.ne.0) then
            if(.not.allocated(pnx0)) then
               allocate(pnx0(neq))
               allocate(pny0(neq))
               allocate(pnz0(neq))
               pnx0 = pnx
               pny0 = pny
               pnz0 = pnz
            endif
            
            if(ipermstr8.ne.0) then
               allocate(frac_flg(neq))
               do i = 1,n0
                  frac_flg(i) = 0
               enddo
            endif
         endif
c     
         
c     
c     model 3 and model 5 require an initial setup      
c     
c     only allocate if there is a model 3 and model 5
c     
         if(ipermstr3.ne.0.or.ipermstr5.ne.0) then
            if(.not.allocated(ipermx)) then
               allocate(ipermx(n0,2))
               allocate(ipermy(n0,2))
               allocate(ipermz(n0,2))
               ipermx = 0
               ipermy = 0
               ipermz = 0
            endif  
c     
c     only calculate for model 3 and model 5
c     initial setup calcs node neighbor information
c     
            jj1 =0 
            do i = 1,n0
               iispmd = ispm(i)
               ispmd = ispmt(iispmd)
               if(ispmd.eq.3.or.ispmd.eq.5) then
                  xi = cord(i,1)
                  yi = cord(i,2)
                  zi = cord(i,3)
                  xdmax = -1.d30
                  xdmin = -1.d30
                  ydmax = -1.d30
                  ydmin = -1.d30
                  zdmax = -1.d30
                  zdmin = -1.d30
                  xdmaxt = -1.d30
                  xdmint = 1.d30
                  ydmaxt = -1.d30
                  ydmint = 1.d30
                  zdmaxt = -1.d30
                  zdmint = 1.d30                  
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
                        dis = abs(xkb-xi)
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        if(dis.gt.xdmax) then
                           kbxmax = kb
                           xdmax = dis
                           xdmaxt = dist
                        else if(abs(dis-xdmax).lt.dis_tol) then
                           if(dist.lt.xdmaxt) then
                              xdmaxt = dist
                              kbxmax = kb
                           endif
                        endif
                     else if(xkb-xi.lt.0) then
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        dis = abs(xkb-xi)
                        if(dis.gt.xdmin) then
                           kbxmin = kb
                           xdmin = dis
                           xdmint = dist
                        else if(abs(dis-xdmin).lt.dis_tol) then
                           if(dist.lt.xdmint) then
                              xdmint = dist
                              kbxmin = kb
                           endif
                        endif
                     endif
                     if(ykb-yi.gt.0) then
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        dis = abs(ykb-yi)
                        if(dis.gt.ydmax) then
                           kbymax = kb
                           ydmax = dis
                           ydmaxt = dist
                        else if(abs(dis-ydmax).lt.dis_tol) then
                           if(dist.lt.ydmaxt) then
                              ydmaxt = dist
                              kbymax = kb
                           endif
                        endif
                     else if(ykb-yi.lt.0) then
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        dis = abs(ykb-yi)
                        if(dis.gt.ydmin) then
                           kbymin = kb
                           ydmin = dis
                           ydmint = dist
                        else if(abs(dis-ydmin).lt.dis_tol) then
                           if(dist.lt.ydmint) then
                              ydmint = dist
                              kbymin = kb
                           endif
                        endif
                     endif 
                     if(zkb-zi.gt.0) then
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        dis = abs(zkb-zi)
                        if(dis.gt.zdmax) then
                           kbzmax = kb
                           zdmax = dis
                           zdmaxt = dist
                        else if(abs(dis-zdmax).lt.dis_tol) then
                           if(dist.lt.zdmaxt) then
                              zdmaxt = dist
                              kbzmax = kb
                           endif
                        endif
                     else if(zkb-zi.lt.0) then
                        dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
                        dis = abs(zkb-zi)
                        if(dis.gt.zdmin) then
                           kbzmin = kb
                           zdmin = dis
                           zdmint = dist
                        else if(abs(dis-zdmin).lt.dis_tol) then
                           if(dist.lt.zdmint) then
                              zdmint = dist
                              kbzmin = kb
                           endif
                        endif
                     endif
                  enddo                       
c     
c     store max and min of each direction (may create (slightly)larger stencil)
c     
                  ipermx(i,1) = kbxmin
                  ipermx(i,2) = kbxmax
                  ipermy(i,1) = kbymin
                  ipermy(i,2) = kbymax	     
                  ipermz(i,1) = kbzmin
                  ipermz(i,2) = kbzmax	  
                  if(cord(kbxmax,1)-cord(kbxmin,1).le.0.0) then
                     jj1 =1
                     write(ierr,*) 'dis(x) failed, node ',
     &                    i,' model 3 or 5 sub stres_perm'
                  endif 
                  if(cord(kbymax,2)-cord(kbymin,2).le.0.0) then
                     jj1 =1
                     write(ierr,*) 'dis(y) failed, node ',
     &                    i,' model 3 or 5 sub stres_perm' 
                  endif 
                  if(cord(kbzmax,3)-cord(kbzmin,3).le.0.0) then
                     jj1 =1
                     write(ierr,*) 'dis(y) failed, node ',
     &                    i,' model 3 or 5 sub stres_perm'
                  endif                             
c     stop for zero distances               
                  if(jj1.ne.0) stop     
               endif 	                 
            enddo
         endif
      endif
      go to 2001
c     
c     test code
c     
      open(99,file='stress_perm',status = 'unknown')
      do i = 1,n0
         iispmd = ispm(i)
         ispmd = ispmt(iispmd)
         if(ispmd.eq.3.or.ispmd.eq.5) then
            write(99,*) 'node = ',i
            write(99,*) 
     &           ' x node pair', ipermx(i,2), ipermx(i,1),
     &           ' dis ', cord(ipermx(i,2),1)-cord(ipermx(i,1),1)
            if(cord(ipermx(i,2),1)-cord(ipermx(i,1),1).le.0.001)
     &           write (99,*) '>>> zero distance x <<<<'
            write(99,*) 
     &           ' y node pair', ipermy(i,2), ipermy(i,1),
     &           ' dis ', cord(ipermy(i,2),2)-cord(ipermy(i,1),2)
            if(cord(ipermy(i,2),2)-cord(ipermy(i,1),2).le.0.001)
     &           write (99,*) '>>> zero distance y <<<<'
            write(99,*) 
     &           ' z node pair', ipermz(i,2), ipermz(i,1),
     &           ' dis ', cord(ipermz(i,2),3)-cord(ipermz(i,1),3)
            if(cord(ipermz(i,2),3)-cord(ipermz(i,1),3).le.0.001)
     &           write (99,*) '>>> zero distance z <<<<'    
         endif
      enddo
      close(99)
      
 2001 continue
      if (iflg.eq.-1) then
c     
c     allocate memory for stress derivatives for fully coupled solution
c     just before call to generate equations
c     
c     
         if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &        .not.allocated(rlxs))then  
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
c     deallocate memory for stress derivatives for fully coupled solution
c     
c     
         if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &        allocated(rlxs))then       
            deallocate(rlxs,rlys,rlzs)
            deallocate(drlxs,drlys,drlzs)
            deallocate(idum_str1)
         endif
c     
      else if (iflg.eq.1) then
c     
c     general loop to calculate permeability models and derivatives
c     
c     skip loop if perm macro is read
         if(.not.allocated(ispm)) return
c     
         do i = 1,n0
c     
c     identify model associated with node i
c     
            iispmd = ispm(i)    
            ispmd = ispmt(iispmd) 
c     
            if(ispmd.eq.1.and.ihms.eq.1) then
c     perm model 1 - not implemented yet (default)
c     not needed unless fully coupled
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
c     perm model 3 fully coupled model
c     
c     perm changes for individual node
c     assumes call has been made to allocate memory
c     
c     x direction (orthogonal pieces)
c     
c     identify parameters
c     
c     
               strx_min = spm1f(iispmd)
               stry_min = spm2f(iispmd) 
               strz_min = spm3f(iispmd) 
               frac_bx = spm4f(iispmd) 
               frac_by = spm5f(iispmd) 
               frac_bz = spm6f(iispmd)
               perx_m = spm7f(iispmd)
               pery_m = spm8f(iispmd) 
               perz_m = spm9f(iispmd) 
               
               kbx1 = ipermx(i,1)
               kbx2 = ipermx(i,2)
               kby1 = ipermy(i,1)
               kby2 = ipermy(i,2)
               kbz1 = ipermz(i,1)                       
               kbz2 = ipermz(i,2)
c     
c     identify displacements
c     
               disx = (cord(kbx2,1)-cord(kbx1,1))/2.
               disy = (cord(kby2,2)-cord(kby1,2))/2.
               disz = (cord(kbz2,3)-cord(kbz1,3))/2.
c     
               amultx = frac_bx**3
               amulty = frac_by**3
               amultz = frac_bz**3
               amultxy = 1./(amultx + amulty)
               amultxz = 1./(amultx + amultz)
               amultyz = 1./(amulty + amultz)
c     
c     displacement terms (orthogonal)
c     
               dispx12 = du(kbx2)-du(kbx1)
               dispy12 = dv(kby2)-dv(kby1)  
               dispz12 = dw(kbz2)-dw(kbz1)  
c     
c     displacement terms (thermal)
c     
               alpi = alp(i)     
               tx12 = (t(i)-tini(i))*alpi*disx/2. 
               ty12 = (t(i)-tini(i))*alpi*disy/2. 
               tz12 = (t(i)-tini(i))*alpi*disz/2. 
c     
c     determine the net contribution  
c     
               disptx12 = max(dispx12-tx12,0.0d0)
               dispty12 = max(dispy12-ty12,0.0d0) 
               disptz12 = max(dispz12-tz12,0.0d0)              
c     
c     displacement terms (shear)
c     
               dispxy12 = du(kby2)-du(kby1)
               dispxz12 = du(kbz2)-du(kbz1)  
               dispyx12 = dv(kbx2)-dv(kbx1)
               dispyz12 = dv(kbz2)-dv(kbz1)
               dispzx12 = dw(kbx2)-dw(kbx1)
               dispzy12 = dw(kby2)-dw(kby1)                         
               
               
               rlxs(i) = amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3
               drlxs(i,2) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))
     &              **2*(1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
               drlxs(i,1) =-3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))
     &              **2*(1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
               drlxs(i,4) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))
     &              **3*(1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
               drlxs(i,3) =-3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))
     &              **3*(1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
               
               rlys(i) = amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3
               drlys(i,2) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))
     &              **2*(1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
               drlys(i,1) =-3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))
     &              **2*(1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
               drlys(i,4) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))
     &              **3*(1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
               drlys(i,3) =-3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))
     &              **3*(1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz  
               
               rlzs(i) = amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**3
               drlzs(i,2) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))
     &              **2*(1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
               drlzs(i,1) =-3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))
     &              **2*(1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
               drlzs(i,4) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))
     &              **3*(1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty
               drlzs(i,3) =-3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))
     &              **3*(1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty     
               
c     
c     thermal contribution based formulation
c     
c     
c     
            else if(ispmd.eq.5) then
c     
c     perm model 5 displacement based model-explicitly coupled
c     
c     perm changes for individual node
c     assumes call has been made to allocate memory
c     
c     x direction (orthogonal pieces)
c     
c     identify parameters (model 5 uses only the rock strengths)
c     
c     
               strx_min = spm1f(iispmd)
               stry_min = spm2f(iispmd) 
               strz_min = spm3f(iispmd) 
c     
c     the parameters in the fully coupled model are replaced 
c     by initial permeabilities
c     
               
               frac_bx = max(spm4f(iispmd),frac_tol) 
               frac_by = max(spm4f(iispmd),frac_tol) 
               frac_bz = max(spm4f(iispmd),frac_tol) 
c     perx_m = spm7f(iispmd)
c     pery_m = spm8f(iispmd) 
c     perz_m = spm9f(iispmd) 
               
               kbx1 = ipermx(i,1)
               kbx2 = ipermx(i,2)
               kby1 = ipermy(i,1)
               kby2 = ipermy(i,2)
               kbz1 = ipermz(i,1)
               kbz2 = ipermz(i,2)
c     
c     identify displacements
c     
               disx = (cord(kbx2,1)-cord(kbx1,1))/2.
               disy = (cord(kby2,2)-cord(kby1,2))/2.
               disz = (cord(kbz2,3)-cord(kbz1,3))/2.
c     
c     
c     different definitions from model 3
c     
               amultx = 1./(disx*frac_bx)
               amulty = 1./(disy*frac_by)
               amultz = 1./(disz*frac_bz)
               amultxy = 1.
               amultxz = 1.
               amultyz = 1.              
c     
c     displacement terms (orthogonal)
c     
               dispx12 = du(kbx2)-du(kbx1)-(du_ini(kbx2)-du_ini(kbx1))
               dispy12 = dv(kby2)-dv(kby1)-(dv_ini(kby2)-dv_ini(kby1))  
               dispz12 = dw(kbz2)-dw(kbz1)-(dw_ini(kbz2)-dw_ini(kbz1)) 
c     
c     displacement terms (thermal)
c     
               alpi = alp(i)     
               tx12 = (t(i)-tini(i))*alpi*disx
               ty12 = (t(i)-tini(i))*alpi*disy
               tz12 = (t(i)-tini(i))*alpi*disz 
c     
c     determine the net contribution  
c     
               disptx12 = max(dispx12-tx12,0.0d0)
               dispty12 = max(dispy12-ty12,0.0d0) 
               disptz12 = max(dispz12-tz12,0.0d0)              
c     
c     displacement terms (shear)
c     
               dispxy12 = du(kby2)-du(kby1)
               dispxz12 = du(kbz2)-du(kbz1)  
               dispyx12 = dv(kbx2)-dv(kbx1)
               dispyz12 = dv(kbz2)-dv(kbz1)
               dispzx12 = dw(kbx2)-dw(kbx1)
               dispzy12 = dw(kby2)-dw(kby1)                         
               
               
               rlxs(i) = amultyz*(1. + amulty*dispty12)**3*
     &              (1. + amultz*disptz12)**3
               drlxs(i,2) = 3.*amultyz*(1. + amulty*dispty12)**2*
     &              (1. + amultz*disptz12)**3*amulty
               drlxs(i,1) = -3.*amultyz*(1. + amulty*dispty12)**2*
     &              (1. + amultz*disptz12)**3*amulty
               drlxs(i,4) = 3.*amultyz*(1. + amulty*dispty12)**3*
     &              (1. + amultz*disptz12)**2*amultz
               drlxs(i,3) = -3.*amultyz*(1. + amulty*dispty12)**3*
     &              (1. + amultz*disptz12)**2*amultz
               
               rlys(i) = amultxz*(1. + amultx*disptz12)**3*
     &              (1. + amultz*disptz12)**3
               drlys(i,2) = 3.*amultxz*(1. + amultx*disptz12)**2*
     &              (1. + amultz*disptz12)**3*amultx
               drlys(i,1) = -3.*amultxz*(1. + amultx*disptz12)**2*
     &              (1. + amultz*disptz12)**3*amultx
               drlys(i,4) = 3.*amultxz*(1. + amultx*disptz12)**3*
     &              (1. + amultz*disptz12)**2*amultz
               drlys(i,3) = -3.*amultxz*(1. + amultx*disptz12)**3*
     &              (1. + amultz*disptz12)**2*amultz  
               
               rlzs(i) = amultxy*(1. + amultx*disptz12)**3*
     &              (1. + amulty*dispty12)**3
               drlzs(i,2) = 3.*amultxy*(1. + amultx*disptz12)**2*
     &              (1. + amulty*dispty12)**3*amultx
               drlzs(i,1) = -3.*amultxy*(1. + amultx*disptz12)**2*
     &              (1. + amulty*dispty12)**3*amultx
               drlzs(i,4) = 3.*amultxy*(1. + amultx*disptz12)**3*
     &              (1. + amulty*dispty12)**2*amulty
               drlzs(i,3) = -3.*amultxy*(1. + amultx*disptz12)**3*
     &              (1. + amulty*dispty12)**2*amulty     
               
c     
c     thermal contribution based formulation
c     
c     
c     now change absolute permeabilities         

               pnx(i) = pnx0(i)*rlxs(i)
               pny(i) = pny0(i)*rlys(i)
               pnz(i) = pnz0(i)*rlzs(i)
               
c     
            else if(ispmd.eq.2) then
c     perm model 2 - volume strains
c     lagged permeability only (stress-based) (only after each time step)
c     only good for 2D or 3D models
c     
c     calculate components of volume strain
c     
               
               dt = t(i) - tini(i)
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
                  
                  pnx(i) = min(perx_m*pnx0(i),pnx(i))
                  pny(i) = min(pery_m*pny0(i),pny(i))
                  pnz(i) = min(perx_m*pnz0(i),pnz(i))	 
                  
               endif
            else if(ispmd.eq.6) then	 
c     
c     same as above but added pore pressure term 
c     
c     perm model 2 - volume strains
c     lagged permeability only (stress-based) (only after each time step)
c     only good for 2D or 3D models
c     
c     need to subtract biot and thermal stresses to get 
c     lithostatic component
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
     &                    max((str_y(i))*lith_damage_fac,lith_min)
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
     &                    max((str_z(i))*lith_damage_fac,lith_min)
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
               
c     Min model
            else if(ispmd.eq.4) then
               
c     perm model 4 - keita model
c     differs from perm model 2 by allowing cubic variation 
c     of perm with stress  
c     lagged permeability (only after each time step)
c     only good for 2-D models
c     
c     calculate components of volume strain
c     
c     dt = t(i) - tini(i)
c     alpv = alp(i)
c     vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
c     ipchk = 0
               
c     kx_fac = 1. + stry_min*abs(str_y(i)**3)
               ipmd4_fx = spm1f(iispmd)
               ipmd4_br = spm2f(iispmd) 
               ipmd4_bmx = spm3f(iispmd) 
               ipmd4_alx = spm4f(iispmd) 
               ipmd4_aly = spm5f(iispmd) 
               ipmd4_fdx = spm6f(iispmd)
               ipmd4_dmx = spm7f(iispmd)
               ipmd4_gmx = spm8f(iispmd) 
               ipmd4_kc = spm9f(iispmd)
               
               ipmd4_fy = spm10f(iispmd)
               ipmd4_btx = spm11f(iispmd)
               ipmd4_bty = spm12f(iispmd)
               ipmd4_fdy = spm13f(iispmd)
               ipmd4_gmy = spm14f(iispmd)
               
c********************have to decide how to do with perx_m etc. input?   
               perx_m = 100
               pery_m = 100     
               
c***********************************NORMAL
               ipmd4_p1 = ipmd4_alx*(str_x(i)-pho(i))+ipmd4_aly*
     &              (str_y(i)-pho(i))
c     maximum aperture dilation	   
               if(ipmd4_p1.lt.-5) ipmd4_p1 = -5
               
c     normal dilation component	   
               ipmd4_p2 = ipmd4_fx/12*(ipmd4_br+ipmd4_bmx*
     &              exp(-ipmd4_p1))**3
               
               ipmd4_p3 = ipmd4_btx*(str_x(i)-pho(i))+ipmd4_bty*
     &              (str_y(i)-pho(i))
               if(ipmd4_p3.lt.-5) ipmd4_p4 = -5
               
               ipmd4_p4 = ipmd4_fy/12*(ipmd4_br+ipmd4_bmx*
     &              exp(-ipmd4_p3))**3
               
c***********************************SHEAR
c     ipmd4_k is the ratio of min stress to max stress
c     if(str_x(i).gt.str_y(i)) then
c     ipmd4_k = str_x(i)/str_y(i)
c     else 
c     ipmd4_k = str_y(i)/str_x(i)   
c     endif
c     if(ipmd4_k.gt.0) then
c     ipmd4_p3 = ipmd4_dmx*(1-exp(-ipmd4_gmx*(ipmd4_k-ipmd4_kc)))
c     else
c     ipmd4_p3 = 0
c     endif
c***********************************SHEAR
               
               pnx(i) = ipmd4_p2*1e6
               pny(i) = ipmd4_p4*1e6
               
               pnx(i) = min(perx_m*pnx0(i),pnx(i))
               pny(i) = min(pery_m*pny0(i),pny(i))
               
c     
c     ******************************Bai et al model (stress based)
c     
            else if(ispmd.eq.7) then
               if(icnl.ne.0) then 
                  ipmd7_bx = spm1f(iispmd)
                  ipmd7_Jx = spm2f(iispmd)
                  ipmd7_sx = spm3f(iispmd)
                  ipmd7_by = spm4f(iispmd)
                  ipmd7_Jy = spm5f(iispmd)
                  ipmd7_sy = spm6f(iispmd)
                  ipmd7_phid = spm10f(iispmd)
                  ipmd7_ksh = spm11f(iispmd)
                  
                  ipmd7_nmx = 0
                  ipmd7_nmy = 0
                  ipmd7_shx = 0
                  ipmd7_shy = 0
                  
                  perx_m=100
                  pery_m=100
c     initial stress - current stress	   
                  ipmd7_dsx = estr_x0(i) - (str_x(i)-phi(i))
                  ipmd7_dsy = estr_y0(i) - (str_y(i)-phi(i))
                  ipmd7_dsxy = abs(str_xy0(i) - str_xy(i))
                  
c     fracture toughness can be a function of stress as well 
c     (not implemented yet)
                  ipmd7_knx = ipmd7_Jx
                  ipmd7_kny = ipmd7_Jy
c     normal dilation 
                  ipmd7_nmx = ipmd7_dsx/(ipmd7_knx*ipmd7_bx) + 
     &                 (ipmd7_sx-ipmd7_bx)/(elastic_mod(i)*ipmd7_bx)
     &                 *(ipmd7_dsx-poisson(i)*ipmd7_dsy) 
                  
                  ipmd7_nmy = ipmd7_dsy/(ipmd7_kny*ipmd7_by) + 
     &                 (ipmd7_sy-ipmd7_by)/(elastic_mod(i)*ipmd7_by)
     &                 *(ipmd7_dsy-poisson(i)*ipmd7_dsx)
                  
                  
                  
c     shear dilation
                  ipmd7_shx = (2*(1+poisson(i))*ipmd7_sx/elastic_mod(i)
     &                 +1/ipmd7_ksh)
     &                 *ipmd7_dsxy/ipmd7_bx*tan(ipmd7_phid/180*pi)
                  
                  ipmd7_shy = (2*(1+poisson(i))*ipmd7_sy/elastic_mod(i)
     &                 +1/ipmd7_ksh)
     &                 *ipmd7_dsxy/ipmd7_by*tan(ipmd7_phid/180*pi)
                  
                  
                  
                  pnx(i) = pnx0(i)*(1+ipmd7_nmy+ipmd7_shy)**3
                  pny(i) = pny0(i)*(1+ipmd7_nmx+ipmd7_shx)**3
                  
c     pnx(i) = pnx0(i)*(1+ipmd7_nmx)**3
c     pny(i) = pny0(i)*(1+ipmd7_nmy)**3
                  
c     pnx(i) = pnx0(i)*(1+ipmd7_shx)**3
c     pny(i) = pny0(i)*(1+ipmd7_shy)**3
                  
                  pnx(i) = min(perx_m*pnx0(i),pnx(i))
                  pny(i) = min(pery_m*pny0(i),pny(i))
                  pnx(i) = max(1.*pnx0(i),pnx(i))
                  pny(i) = max(1.*pny0(i),pny(i))  
               endif        
c     
c************************************3D
               if(icnl.eq.0) then 
                  
                  ipmd7_bx = spm1f(iispmd)
                  ipmd7_Jx = spm2f(iispmd)
                  ipmd7_sx = spm3f(iispmd)
                  ipmd7_by = spm4f(iispmd)
                  ipmd7_Jy = spm5f(iispmd)
                  ipmd7_sy = spm6f(iispmd)
                  ipmd7_bz = spm7f(iispmd)
                  ipmd7_Jz = spm8f(iispmd)
                  ipmd7_sz = spm9f(iispmd)
                  ipmd7_phid = spm10f(iispmd)
                  ipmd7_ksh = spm11f(iispmd)
                  
                  ipmd7_nmx = 0
                  ipmd7_nmy = 0
                  ipmd7_nmz = 0
                  ipmd7_shx = 0
                  ipmd7_shy = 0
                  ipmd7_shz = 0
                  
                  perx_m=100
                  pery_m=100
                  perz_m=100
c     initial stress - current stress	   
                  ipmd7_dsx = estr_x0(i) - (str_x(i)-phi(i))
                  ipmd7_dsy = estr_y0(i) - (str_y(i)-phi(i))
                  ipmd7_dsz = estr_z0(i) - (str_z(i)-phi(i))
                  ipmd7_dsxy = abs(str_xy0(i) - str_xy(i))
                  ipmd7_dsxz = abs(str_xz0(i) - str_xz(i))
                  ipmd7_dsyz = abs(str_yz0(i) - str_yz(i))
                  
c     fracture toughness can be a function of stress as well 
c     (not implemented yet)
                  ipmd7_knx = ipmd7_Jx
                  ipmd7_kny = ipmd7_Jy
                  ipmd7_knz = ipmd7_Jz
                  
c     normal dilation         	   
                  ipmd7_nmx = ipmd7_dsx/(ipmd7_knx*ipmd7_bx) + 
     &                 (ipmd7_sx-ipmd7_bx)/(elastic_mod(i)*ipmd7_bx)
     &                 *(ipmd7_dsx-poisson(i)*
     &                 (ipmd7_dsy+ipmd7_dsz)) 
                  
                  ipmd7_nmy = ipmd7_dsy/(ipmd7_kny*ipmd7_by) + 
     &                 (ipmd7_sy-ipmd7_by)/(elastic_mod(i)*ipmd7_by)
     &                 *(ipmd7_dsy-poisson(i)*(ipmd7_dsx+ipmd7_dsz))
                  
                  ipmd7_nmz = ipmd7_dsz/(ipmd7_knz*ipmd7_bz) + 
     &                 (ipmd7_sz-ipmd7_bz)/
     &                 (elastic_mod(i)*ipmd7_bz)
     &                 *(ipmd7_dsz-poisson(i)*
     &                 (ipmd7_dsx+ipmd7_dsy))
                  
c     shear dilation
                  ipmd7_shx = (2*(1+poisson(i))*ipmd7_sx/
     &                 elastic_mod(i)
     &                 +1/ipmd7_ksh)
     &                 *(ipmd7_dsxy+ipmd7_dsxz)/
     &                 ipmd7_bx*tan(ipmd7_phid/180*pi)
                  
                  ipmd7_shy = (2*(1+poisson(i))*ipmd7_sy/
     &                 elastic_mod(i)
     &                 +1/ipmd7_ksh)
     &                 *(ipmd7_dsxy+ipmd7_dsyz)/
     &                 ipmd7_by*tan(ipmd7_phid/180*pi)
                  
                  ipmd7_shz = (2*(1+poisson(i))*
     &                 ipmd7_sz/elastic_mod(i)
     &                 +1/ipmd7_ksh)
     &                 *(ipmd7_dsxz+ipmd7_dsyz)/
     &                 ipmd7_bz*tan(ipmd7_phid/180*pi)
                  
                  pnx(i) = pnx0(i)/(ipmd7_by**3+ipmd7_bz**3)*
     &                 (ipmd7_by**3*(1+ipmd7_nmy+ipmd7_shy)**3 
     &                 +ipmd7_bz**3*(1+ipmd7_nmz+ipmd7_shz)**3)
                  
                  pny(i) = pny0(i)/(ipmd7_bx**3+ipmd7_bz**3)*
     &                 (ipmd7_bx**3*(1+ipmd7_nmx+ipmd7_shx)**3
     &                 +ipmd7_bz**3*(1+ipmd7_nmz+ipmd7_shz)**3)
                  
                  pnz(i) = pnz0(i)/(ipmd7_bx**3+ipmd7_by**3)*
     &                 (ipmd7_bx**3*(1+ipmd7_nmx+ipmd7_shx)**3
     &                 +ipmd7_by**3*(1+ipmd7_nmy+ipmd7_shy)**3)
                  
                  
                  pnx(i) = min(perx_m*pnx0(i),pnx(i))
                  pny(i) = min(pery_m*pny0(i),pny(i))
                  pnz(i) = min(perz_m*pnz0(i),pnz(i))
                  pnx(i) = max(1.*pnx0(i),pnx(i))
                  pny(i) = max(1.*pny0(i),pny(i))  
                  pnz(i) = max(1.*pnz0(i),pnz(i))
               endif   
               
               
c     **************failure criteria + stress-perm (for intact rock)
            else if(ispmd.eq.8) then
               
               
               if(icnl.ne.0) then 
                  ipmd8_bx = spm1f(iispmd)
                  ipmd8_Jx = spm2f(iispmd)
                  ipmd8_sx = spm3f(iispmd)
                  ipmd8_by = spm4f(iispmd)
                  ipmd8_Jy = spm5f(iispmd)
                  ipmd8_sy = spm6f(iispmd)
                  ipmd8_phid = spm10f(iispmd)
                  ipmd8_ksh = spm11f(iispmd)
c     Tensile strength and Coulomb criteria
                  ipmd8_tns = spm12f(iispmd)
                  ipmd8_clb = spm13f(iispmd)	   
                  
                  perx_m=100
                  pery_m=100
c     2D	    
c     frac_flg = 0 (no fracture)
c     = 1 (frac - x direction)
c     = 2 (frac - y direction)
c     = 3 (frac - z direction)
c     = 4 (frac - x and y direction)
c     = 5 (frac - x and z direction)
c     = 6 (frac - y and z direction)
c     = 7 (frac - all direction)   
                  
c     Tensile criteria
                  if((str_y(i)-phi(i)).lt.ipmd8_tns) then
                     if(frac_flg(i).eq.0) then
                        frac_flg(i) = 1
                        pnz(i) = 1.e5
                        estr_x0(i) = str_x(i)-phi(i)
                        pnx0(i) = 1.2*pnx0(i)
                     endif
                     if(frac_flg(i).eq.2) then
                        frac_flg(i) = 4
                        pnz(i) = 1.e2
                        estr_x0(i) = str_x(i)-phi(i)
                        pnx0(i) = 1.2*pnx0(i)
                     endif
                  endif
                  if((str_x(i)-phi(i)).lt.ipmd8_tns) then
                     if(frac_flg(i).eq.0) then
                        frac_flg(i) = 2
                        pnz(i) = 1.e4
                        estr_y0(i) = str_y(i)-phi(i)
                        pny0(i) = 1.2*pny0(i)
                     endif
                     if(frac_flg(i).eq.1) then
                        frac_flg(i) = 4
                        pnz(i) = 1.e2
                        estr_y0(i) = str_y(i)-phi(i)
                        pny0(i) = 1.2*pny0(i)
                     endif
                  endif
                  
                  
c     Shear criteria
c     if( (str_x(i)-phi(i)).gt.0.and.(str_y(i)-phi(i)).gt.0 ) then
c     if((str_x(i)-phi(i))/(str_y(i)-phi(i)).gt.ipmd8_clb) then
c     if(frac_flg(i).eq.0) then
c     frac_flg(i) = 1
c     pnz(i) = 1.e5
c     estr_x0(i) = str_x(i)-phi(i)
c     pnx0(i) = 1.2*pnx0(i)
c     endif
c     if(frac_flg(i).eq.2) then
c     frac_flg(i) = 4
c     pnz(i) = 1.e2
c     estr_x0(i) = str_x(i)-phi(i)
c     pnx0(i) = 1.2*pnx0(i)
c     endif
c     endif
                  
c     if((str_y(i)-phi(i))/(str_x(i)-phi(i)).gt.ipmd8_clb) then
c     if(frac_flg(i).eq.0) then
c     frac_flg(i) = 2
c     pnz(i) = 1.e4
c     estr_y0(i) = str_y(i)-phi(i)
c     pny0(i) = 1.2*pny0(i)
c     endif
c     if(frac_flg(i).eq.1) then
c     frac_flg(i) = 4
c     pnz(i) = 1.e2
c     estr_y0(i) = str_y(i)-phi(i)
c     pny0(i) = 1.2*pny0(i)
c     endif
c     endif
c     endif       
                  
c     normal dilation
                  ipmd8_nmx = 0
                  ipmd8_nmy = 0     	   
                  
c     normail dilation in x-direction	   
                  if(frac_flg(i).eq.1.or.frac_flg(i).eq.4) then
                     ipmd8_dsx = ipmd8_tns - (str_x(i)-phi(i))
                     ipmd8_dsy = estr_y0(i) - (str_y(i)-phi(i))
                     ipmd8_nmx = ipmd8_dsy/(ipmd8_kny*ipmd8_by) + 
     &                    (ipmd8_sy-ipmd8_by)/(elastic_mod(i)*ipmd8_by)
     &                    *(ipmd8_dsy-poisson(i)*ipmd8_dsx)
                  endif
c     normail dilation in y-direction         
                  if(frac_flg(i).eq.2.or.frac_flg(i).eq.4) then
                     ipmd8_dsx = estr_x0(i) - (str_y(i)-phi(i))
                     ipmd8_dsy = ipmd8_tns - (str_y(i)-phi(i))
                     ipmd8_nmy = ipmd8_dsx/(ipmd8_knx*ipmd8_bx) +
     &                    (ipmd8_sx-ipmd8_bx)/(elastic_mod(i)*ipmd8_bx)
     &                    *(ipmd8_dsx-poisson(i)*ipmd8_dsy)
                  endif
                  pnx(i) = pnx0(i)*(1+ipmd7_nmx)**3
                  pny(i) = pny0(i)*(1+ipmd7_nmy)**3	
                  
                  
c     pnx(i) = pnx0(i)*10.0
                  pnx(i) = min(perx_m*pnx0(i),pnx(i))
                  pny(i) = min(pery_m*pny0(i),pny(i))
                  pnx(i) = max(1.*pnx0(i),pnx(i))
                  pny(i) = max(1.*pny0(i),pny(i)) 
               endif
c     
c     3D
c     
               if(icnl.eq.0) then 
                  ipmd8_bx = spm1f(iispmd)
                  ipmd8_Jx = spm2f(iispmd)
                  ipmd8_sx = spm3f(iispmd)
                  ipmd8_by = spm4f(iispmd)
                  ipmd8_Jy = spm5f(iispmd)
                  ipmd8_sy = spm6f(iispmd)
                  ipmd8_bz = spm7f(iispmd)
                  ipmd8_Jz = spm8f(iispmd)
                  ipmd8_sz = spm9f(iispmd)
                  ipmd8_phid = spm10f(iispmd)
                  ipmd8_ksh = spm11f(iispmd)
c     Tensile strength and Coulomb criteria
                  ipmd8_tns = spm12f(iispmd)
                  ipmd8_clb = spm13f(iispmd)	   
                  
                  perx_m=100
                  pery_m=100
                  
c     frac_flg = 0 (no fracture)
c     = 1 (frac - xz plane, vertical   -> kx)
c     = 2 (frac - yz plane, vertical   -> ky)
c     = 3 (frac - xy plane, horizontal -> kz)
c     = 4 (frac - xz and xy plane 1&3  -> kx)
c     = 5 (frac - yz and xy plane 2&3  -> ky)
c     = 6 (frac - xz and yz plane 1&2  -> kz)
c     = 7 (frac - all direction        -> just for record)	   
c     Tensile criteria
                  
                  if((str_y(i)-phi(i)).lt.ipmd8_tns) then
                     if(frac_flg(i).eq.0) then
                        frac_flg(i) = 1
                        pnx(i) = 10*pnx0(i)
                        pnz(i) = 10*pnz0(i)
                     elseif(frac_flg(i).eq.3) then
                        frac_flg(i) = 4
                        pnx(i) = 10*pnx(i)
                     elseif(frac_flg(i).eq.2) then
                        frac_flg(i) = 6
                        pnz(i) = 10*pnz(i)
                     endif
                  endif
                  if((str_x(i)-phi(i)).lt.ipmd8_tns) then
                     if(frac_flg(i).eq.0) then
                        frac_flg(i) = 2
                        pny(i) = 10*pny0(i)
                        pnz(i) = 10*pnz0(i)
                     elseif(frac_flg(i).eq.3) then
                        frac_flg(i) = 5
                        pny(i) = 10*pny(i)
                     elseif(frac_flg(i).eq.1) then
                        frac_flg(i) = 6
                        pnz(i) = 10*pnz(i)
                     endif
                  endif
                  if((str_z(i)-phi(i)).lt.ipmd8_tns) then
                     if(frac_flg(i).eq.0) then
                        frac_flg(i) = 3
                        pnx(i) = 10*pnx0(i)
                        pny(i) = 10*pny0(i)
                     elseif(frac_flg(i).eq.1) then
                        frac_flg(i) = 4
                        pnx(i) = 10*pnx(i)
                     elseif(frac_flg(i).eq.2) then
                        frac_flg(i) = 5
                        pny(i) = 10*pny(i)
                     endif         
                  endif 
                  
               endif
c...........................................................
            elseif(ispmd.eq.11) then
c     s kelkar June 9 2009, simple gangi (1978, Int.J.Rock Mech)
c     model in 2-D
c     assumed model in x-y plane and fracture in x-z plane, so only 
c     x-perm is modified. This formulation is for fracture faces
c     in contact, so str_y = effective stress in y-dir has to be
c     compressive
               if(icnl.ne.0) then
                  gk0  =spm1f(iispmd)
                  gpmod=spm2f(iispmd)
                  gmexp =spm3f(iispmd)
                  gn=1./gmexp
                  sigy_eff=str_y(i)
                  if(sigy_eff.gt.0.d0) then
                     pnx(i)=gk0*(1.-(sigy_eff/gpmod)**gmexp)**3.
                     e2(i) =  gn*sigy_eff*(gpmod/sigy_eff)**gn
                  endif
               endif
c...........................................................

            endif   
            
         enddo
         
         
c     
c     check for damage zone (permeability changes)
c     and maximum allowable changes
c     
         if(ipermstr2.ne.0.or.ipermstr5.ne.0.or.ipermstr7.ne.0
     &        .or.ipermstr8.ne.0) then
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
c     count damaged zone nodes
                  ipr_str = ipr_str +1
c     find minimum stresses in the damaged zone           
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
c     write out information for the damaged zones    
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
            ipermx_max = 308
            if(iout.ne.0) write(iout,102)
     &           ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max
*******************************
            
c     &  write(iout,201)ipermx_max,str_y(ipermx_max)-phi(ipermx_max),
c     &        frac_flg(ipermx_max)
c*******************************
            if(iptty.ne.0) write(iptty,102)
     &           ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max 
            if(iout.ne.0) write(iout,103)
     &           ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            
c     &  write(iout,106)ipermy_max,
c     &       estr_x0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       estr_y0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       pny(ipermx_max)*1.e-6,
c     &       cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            if(iptty.ne.0) write(iptty,103)
     &           ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
            if(icnl.eq.0) then
               if(iout.ne.0) write(iout,104)
     &              ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &              cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max  
               if(iptty.ne.0) write(iptty,104)
     &              ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &              cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max    
            endif     
         endif         
      endif
 99   format(/,1x'Time step ',i6,' Days'1x,f9.2)   
 100  format(1x,'Number of damaged gridbolcks (gt 0.01 k/k0 ) ', 1x,i6)
 101  format(1x,'Largest tensile stresses in damaged zone')
 102  format(1x,'Node ',i6,1x,'Sxx ',g12.4,1x,'Kxx ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 103  format(1x,'Node ',i6,1x,'Syy ',g12.4,1x,'Kyy ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 104  format(1x,'Node ',i6,1x,'Szz ',g12.4,1x,'Kzz ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)            
      
 105  format(1x,'Node ',i6,1x,'Sxx ',g12.4,1x,'Kxo',g12.4,1x
     &     'Kxx ',g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 106  format(1x,'Node ',i6,1x,'Syy ',g12.4,1x,'Kyo ',g12.4,1x
     &     'Kyy ',g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)	 
 201  format(1x,'Node ',i6,1x,'eSyy',g12.4,'frac_flg ',i6)   	 
      
      return 
      end	 
      
