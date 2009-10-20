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
      real*8 tx12, ty12, tz12, phid, phib
      real*8 biot,erat,efac,epi,dpd,shpi,stress_min_prin,lith_min
      real*8 ctherm,eti,shti,e1i,e2i,e3i,lith_damage_fac,str_eff,pterm
      
c****************************local parameters used in perm model 4      
      real*8  ipmd4_p1,ipmd4_p2,ipmd4_p3,ipmd4_p4
      real*8  ipmd4_k
c**************************** parameters added for perm model 5, may need to be merged with 7&8
      real*8  knx, kny, knz, ksx, ksy, ksz
      real*8  epsxx12, epsyy12, epszz12 
      real*8  epsxy12, epsxz12, epsyz12
      real*8  sisjx, sisjy, sisjz           
c****************************local parameters used in perm model 7 Bai   
      real*8  ipmd7_dsx,ipmd7_dsy,ipmd7_dsz
      real*8  ipmd7_dsxy,ipmd7_dsxz,ipmd7_dsyz
      real*8  ipmd7_knx,ipmd7_kny,ipmd7_knz
      real*8  ipmd7_nmx,ipmd7_nmy,ipmd7_nmz
      real*8  ipmd7_shx,ipmd7_shy,ipmd7_shz
      real*8  pi
      real*8  frac_tol
c****************************local parameters used in perm model 8  Bai 
      real*8  ipmd8_dsx,ipmd8_dsy,ipmd8_dsz
      real*8  ipmd8_dsxy,ipmd8_dsyz,ipmd8_dsxz
      real*8  ipmd8_tns,ipmd8_clb,cohes,sh_cohe
      real*8  ipmd8_nmx,ipmd8_nmy,ipmd8_nmz
      real*8  ipmd8_shx,ipmd8_shy,ipmd8_shz
      real*8  esxi, esyi, eszi
      real*8  dsx12,dsy12,dsz12,dsxy12,dsxz12,dsyz12
      real*8  L11,L12,L13,L21,L22,L23,L31,L32,L33
      
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
c     only allocate if there is a model 3, model 5 and model 8
c     
         if(ipermstr3.ne.0.or.ipermstr5.ne.0.or.ipermstr7.ne.0.
     &      .or.ipermstr8.ne.0) then
            if(.not.allocated(ipermx)) then
               allocate(ipermx(n0,2))
               allocate(ipermy(n0,2))
               allocate(ipermz(n0,2))
               ipermx = 0
               ipermy = 0
               ipermz = 0
            endif  
c     
c     only calculate for model 3, model 5, model 7 and model 8
c     initial setup calcs node neighbor information
c    
            jj1 =0 
            do i = 1,n0
               iispmd = ispm(i)
               ispmd = ispmt(iispmd)
               if(ispmd.eq.3.or.ispmd.eq.5.or.ispmd.eq.7
     &            .or.ispmd.eq.8) then
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
                write(ierr,*) 
     &          'dis(x) failed, node ',i,' model 3 or 5 sub stres_perm'
               endif 
               if(cord(kbymax,2)-cord(kbymin,2).le.0.0) then
                jj1 =1
                write(ierr,*) 
     &          'dis(y) failed, node ',i,' model 3 or 5 sub stres_perm' 
               endif 
               if(cord(kbzmax,3)-cord(kbzmin,3).le.0.0) then
                jj1 =1
                write(ierr,*) 
     &          'dis(z) failed, node ',i,' model 3 or 5 sub stres_perm'
               endif                             
c stop for zero distances               
              if(jj1.ne.0) stop     
             endif 	                 
            enddo
         endif
      endif
      go to 2001
c
c test code
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
     &             write (99,*) '>>> zero distance x <<<<'
                    write(99,*) 
     &           ' y node pair', ipermy(i,2), ipermy(i,1),
     &           ' dis ', cord(ipermy(i,2),2)-cord(ipermy(i,1),2)
                 if(cord(ipermy(i,2),2)-cord(ipermy(i,1),2).le.0.001)
     &             write (99,*) '>>> zero distance y <<<<'
                    write(99,*) 
     &           ' z node pair', ipermz(i,2), ipermz(i,1),
     &           ' dis ', cord(ipermz(i,2),3)-cord(ipermz(i,1),3)
                  if(cord(ipermz(i,2),3)-cord(ipermz(i,1),3).le.0.001)
     &             write (99,*) '>>> zero distance z <<<<'    
               endif
              enddo
              close(99)
    
2001  continue   
      if (iflg.eq.-1) then
c     
c     allocate memory for stress derivatives for fully coupled solution
c     just before call to generate equations
c     
c     
         if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &         .not.allocated(rlxs))then  
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
     &                    allocated(rlxs))then       
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
               disptx12 = max(dispx12-tx12,0.0)
               dispty12 = max(dispy12-ty12,0.0) 
               disptz12 = max(dispz12-tz12,0.0)              
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
          drlxs(i,2) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**2*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
          drlxs(i,1) = -3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**2*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
          drlxs(i,4) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
          drlxs(i,3) = -3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
               
               rlys(i) = amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3
          drlys(i,2) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
          drlys(i,1) = -3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
          drlys(i,4) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
          drlys(i,3) = -3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz  
               
               rlzs(i) = amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**3
          drlzs(i,2) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
          drlzs(i,1) = -3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
          drlzs(i,4) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty
          drlzs(i,3) = -3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &              (1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty     
               
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
c  the parameters in the fully coupled model are replaced 
c  by initial permeabilities
c               
               
               frac_bx = max(spm4f(iispmd),frac_tol) 
               frac_by = max(spm4f(iispmd),frac_tol) 
               frac_bz = max(spm4f(iispmd),frac_tol) 
c               perx_m = spm7f(iispmd)
c               pery_m = spm8f(iispmd) 
c               perz_m = spm9f(iispmd) 
               
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
c different definitions from model 3
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
               dispx12 = 
     &          (du(kbx2)-du(kbx1)-(du_ini(kbx2)-du_ini(kbx1)))/2.
               dispy12 = 
     &          (dv(kby2)-dv(kby1)-(dv_ini(kby2)-dv_ini(kby1)))/2.  
               dispz12 = 
     &          (dw(kbz2)-dw(kbz1)-(dw_ini(kbz2)-dw_ini(kbz1)))/2. 
               epsxx12 = dispx12/disx
               epsyy12 = dispy12/disy
               epszz12 = dispz12/disz 
c     
c     displacement terms (thermal)
c     
               alpi = alp(i)     
               dt = t(i)-tini(i)
               efac = 3.d0*e2(i)+2.d0*e3(i)        
c               tx12 = (t(i)-tini(i))*alpi*disx
c               ty12 = (t(i)-tini(i))*alpi*disy
c               tz12 = (t(i)-tini(i))*alpi*disz 

c     
c     determine the net contribution  
c     

               dispx12= dispx12
     &         +(e1(i)*epsxx12+e2(i)*(epsyy12+epszz12)
     &         -alpi*efac*dt)/knx     

               dispy12 = dispy12
     &         +(e1(i)*epsyy12+e2(i)*(epsxx12+epszz12)
     &         -alpi*efac*dt)/kny
               
               dispz12 = dispz12
     &         +(e1(i)*epszz12+e2(i)*(epsxx12+epsyy12)
     &         -alpi*efac*dt)/knz
               
               disptx12 = max(dispx12,0.0)
               dispty12 = max(dispy12,0.0) 
               disptz12 = max(dispz12,0.0)              
c     
c     displacement terms (shear)
c     
               
               epsxy12 = abs((du(kby2)-du(kby1))/2)/disy
     &                  +abs((dv(kbx2)-dv(kbx1))/2)/disx          
               epsxz12 = abs((du(kbz2)-du(kbz1))/2)/disz
     &                  +abs((dw(kbx2)-dw(kbx1))/2)/disx
               epsyz12 = abs((dv(kbz2)-dv(kbz1))/2)/disz
     &                  +abs((dw(kby2)-dw(kby1))/2)/disy
               
            if(str_x(i)-phi(i).gt.0) then
                dispyx12 = (disy+e3(i)/ksy)*epsxy12
                dispzx12 = (disz+e3(i)/ksz)*epsxz12
            endif
            
            if(str_y(i)-phi(i).gt.0) then
               dispxy12 = (disx+e3(i)/ksx)*epsxy12
               dispzy12 = (disz+e3(i)/ksz)*epsyz12  
            endif

            if(str_z(i)-phi(i).gt.0) then
                dispxz12 = (disx+e3(i)/ksx)*epsxz12
                dispyz12 = (disy+e3(i)/ksy)*epsyz12
            endif

c               dispxy12 = (disx+e3(i)/ksx)*epsxy12
c               dispxz12 = (disx+e3(i)/ksx)*epsxz12
c               dispyx12 = (disy+e3(i)/ksy)*epsxy12
c               dispyz12 = (disy+e3(i)/ksy)*epsyz12
c               dispzx12 = (disz+e3(i)/ksz)*epsxz12
c               dispzy12 = (disz+e3(i)/ksz)*epsyz12  

c
c     perm enhancement                               
c

          if(frac_by.eq.0.and.frac_bz.eq.0) then
             rlxs(i) = 1.
          else
             rlxs(i) = 1/(frac_by**3+frac_bz**3)*((frac_by+dispty12
     &        +(dispxy12+dispzy12)*tan(phid/180*pi))**3+
     &       (frac_bz+disptz12+(dispyz12+dispxz12)*tan(phid/180*pi))**3)
          endif
          
          if(frac_bx.eq.0.and.frac_bz.eq.0) then
             rlys(i) = 1.
          else
             rlys(i) = 1/(frac_bx**3+frac_bz**3)*((frac_bx+disptx12
     &       +(dispyx12+dispzx12)*tan(phid/180*pi))**3+
     &       (frac_bz+disptz12+(dispyz12+dispxz12)*tan(phid/180*pi))**3)
          endif
          
          if(frac_bx.eq.0.and.frac_by.eq.0) then
             rlzs(i) = 1.
          else
             rlzs(i) = 1/(frac_bx**3+frac_by**3)*((frac_bx+disptx12
     &       +(dispyx12+dispzx12)*tan(phid/180*pi))**3+
     &       (frac_by+dispty12+(dispxy12+dispzy12)*tan(phid/180*pi))**3)
          endif
          
ccccccccccccccccccccccccccccccccccccccccccccccccccccCHECK          
             check(i)= (dispxy12+dispzy12)*tan(phid/180*pi)
               
c               rlxs(i) = amultyz*(1. + amulty*dispty12)**3*
c     &              (1. + amultz*disptz12)**3
          drlxs(i,2) = 3.*amultyz*(1. + amulty*dispty12)**2*
     &              (1. + amultz*disptz12)**3*amulty
          drlxs(i,1) = -3.*amultyz*(1. + amulty*dispty12)**2*
     &              (1. + amultz*disptz12)**3*amulty
          drlxs(i,4) = 3.*amultyz*(1. + amulty*dispty12)**3*
     &              (1. + amultz*disptz12)**2*amultz
          drlxs(i,3) = -3.*amultyz*(1. + amulty*dispty12)**3*
     &              (1. + amultz*disptz12)**2*amultz
               
c               rlys(i) = amultxz*(1. + amultx*disptz12)**3*
c     &              (1. + amultz*disptz12)**3
          drlys(i,2) = 3.*amultxz*(1. + amultx*disptz12)**2*
     &              (1. + amultz*disptz12)**3*amultx
          drlys(i,1) = -3.*amultxz*(1. + amultx*disptz12)**2*
     &              (1. + amultz*disptz12)**3*amultx
          drlys(i,4) = 3.*amultxz*(1. + amultx*disptz12)**3*
     &              (1. + amultz*disptz12)**2*amultz
          drlys(i,3) = -3.*amultxz*(1. + amultx*disptz12)**3*
     &              (1. + amultz*disptz12)**2*amultz  
               
c               rlzs(i) = amultxy*(1. + amultx*disptz12)**3*
c     &              (1. + amulty*dispty12)**3
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
      ipmd4_p1 = ipmd4_alx*(str_x(i)-pho(i))+ipmd4_aly*(str_y(i)-pho(i))
c     maximum aperture dilation	   
               if(ipmd4_p1.lt.-5) ipmd4_p1 = -5
               
c     normal dilation component	   
      ipmd4_p2 = ipmd4_fx/12*(ipmd4_br+ipmd4_bmx*exp(-ipmd4_p1))**3
               
      ipmd4_p3 = ipmd4_btx*(str_x(i)-pho(i))+ipmd4_bty*(str_y(i)-pho(i))
               if(ipmd4_p3.lt.-5) ipmd4_p4 = -5
               
      ipmd4_p4 = ipmd4_fy/12*(ipmd4_br+ipmd4_bmx*exp(-ipmd4_p3))**3
               
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
c     ************************ 2D is not implemented
	  if(icnl.ne.0) then 
	     write(ierr,*) 'Bai model not implemented in 2D'
	     write(ierr,*) 'stopping in stress_perm.f'

        endif        
c
c **************************** 3D
c
	  if(icnl.eq.0) then 

	   frac_bx = spm1f(iispmd)
	   frac_by = spm2f(iispmd)
	   frac_bz = spm3f(iispmd)
	   knx = spm4f(iispmd)
	   ksx = spm5f(iispmd)
	   phid = spm6f(iispmd)

         kny = knx
         knz = knx
         ksy = ksx
         ksz = ksx	     
	   ipmd7_nmx = 0.
	   ipmd7_nmy = 0.
	   ipmd7_nmz = 0.
	   ipmd7_shx = 0.
	   ipmd7_shy = 0.
	   ipmd7_shz = 0.
	     
	   perx_m=1.e4
	   pery_m=1.e4
	   perz_m=1.e4
c initial stress - current stress	   
	   ipmd7_dsx = estr_x0(i) - (str_x(i)-phi(i))
	   ipmd7_dsy = estr_y0(i) - (str_y(i)-phi(i))
	   ipmd7_dsz = estr_z0(i) - (str_z(i)-phi(i))
	   ipmd7_dsxy = abs(str_xy0(i) - str_xy(i))
	   ipmd7_dsxz = abs(str_xz0(i) - str_xz(i))
	   ipmd7_dsyz = abs(str_yz0(i) - str_yz(i))
	   
	   if(ipmd7_dsx.lt.0) ipmd7_dsx = 0.
	   if(ipmd7_dsy.lt.0) ipmd7_dsy = 0.
	   if(ipmd7_dsz.lt.0) ipmd7_dsz = 0.

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
c       fracture dilation (normal and shear)        
c
         if(frac_bx.ne.0.) then
           ipmd7_nmx = ipmd7_dsx/(knx*frac_bx) + disx/frac_bx
     &     *((ipmd7_dsx-poisson(i)*(ipmd7_dsy+ipmd7_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(str_x(i)-phi(i).lt.0) then
             ipmd7_shx = 0.
           else     
             ipmd7_shx = ((disy/e3(i)+1/ksy)*ipmd7_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd7_dsxz)*tan(phid/180*pi)/frac_bx
           endif
           
         endif 
         
         if(frac_by.ne.0.) then
     	     ipmd7_nmy = ipmd7_dsy/(kny*frac_by) + disy/frac_by
     &     *((ipmd7_dsy-poisson(i)*(ipmd7_dsx+ipmd7_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(str_y(i)-phi(i).lt.0) then
             ipmd7_shy = 0.
           else
             ipmd7_shy = ((disx/e3(i)+1/ksx)*ipmd7_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd7_dsyz)*tan(phid/180*pi)/frac_by
           endif
           
         endif
     
         if(frac_bz.ne.0.) then
	     ipmd7_nmz = ipmd7_dsz/(knz*frac_bz) + disz/frac_bz
     &     *((ipmd7_dsz-poisson(i)*(ipmd7_dsx+ipmd7_dsy))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(str_z(i)-phi(i).lt.0) then
             ipmd7_shy = 0.
           else
             ipmd7_shz = ((disx/e3(i)+1/ksx)*ipmd7_dsxz
     &       +(disy/e3(i)+1/ksy)*ipmd7_dsyz)*tan(phid/180*pi)/frac_bz
           endif
           
         endif

cccccccccccccccccccccccccccccccccccCHECK            
          check(i) = ipmd7_shy*frac_by    
c
c         perm enhancement    
c

          if(frac_by.ne.0.or.frac_bz.ne.0) then
             pnx(i) = pnx0(i)/(frac_by**3+frac_bz**3)*
     &     ( frac_by**3*(1+ipmd7_nmy+ipmd7_shy)**3 
     &       + frac_bz**3*(1+ipmd7_nmz+ipmd7_shz)**3 ) 
          endif
          
          if(frac_bx.ne.0.or.frac_bz.ne.0) then
            pny(i) = pny0(i)/(frac_bx**3+frac_bz**3)*
     &     ( frac_bx**3*(1+ipmd7_nmx+ipmd7_shx)**3
     &       + frac_bz**3*(1+ipmd7_nmz+ipmd7_shz)**3 )                 
          endif
     
          if(frac_bx.ne.0.or.frac_by.ne.0) then
            pnz(i) = pnz0(i)/(frac_bx**3+frac_by**3)*
     &     ( frac_bx**3*(1+ipmd7_nmx+ipmd7_shx)**3
     &       + frac_by**3*(1+ipmd7_nmy+ipmd7_shy)**3 )
          endif     

          
c         pnx(i) = min(perx_m*pnx0(i),pnx(i))
c         pny(i) = min(pery_m*pny0(i),pny(i))
c         pnz(i) = min(perz_m*pnz0(i),pnz(i))
c         pnx(i) = max(0.01*pnx0(i),pnx(i))
c         pny(i) = max(0.01*pny0(i),pny(i))  
c         pnz(i) = max(0.01*pnz0(i),pnz(i))

        endif   
c     **************failure criteria + stress-perm (for intact rock)
            else if(ispmd.eq.8) then
               
               
               if(icnl.ne.0) then 
               
c   ******** 2D is not implemented yet

               endif

c   ******** 3D 
c     
               if(icnl.eq.0) then 
                  frac_bx = spm1f(iispmd)
                  frac_by = spm2f(iispmd)
                  frac_bz = spm3f(iispmd)
                  knx = spm4f(iispmd)
                  ksx = spm5f(iispmd)
                  phid = spm6f(iispmd)
                  ipmd8_tns = spm7f(iispmd)
                  phib = spm8f(iispmd)
                  cohes = spm9f(iispmd)
                  
                  ipmd8_clb = (1+sin(phib/180*pi))/(1-sin(phib/180*pi))
                  sh_cohe=2.*cohes*cos(phib/180*pi)/(1-sin(phib/180*pi))
                  

	    perx_m=1.e4
	    pery_m=1.e4
	    perz_m=1.e4
c Fracture toughness is input (can be a funciton of stress)	    
	    kny = knx
	    knz = knx
	    ksy = ksx
	    ksz = ksx
	    
	    ipmd8_nmx = 0.
	    ipmd8_nmy = 0.
	    ipmd8_nmz = 0.
	    ipmd8_shx = 0.
	    ipmd8_shy = 0.

c Distance	    
          kbx1 = ipermx(i,1)
          kbx2 = ipermx(i,2)
          kby1 = ipermy(i,1)
          kby2 = ipermy(i,2)
	    kbz1 = ipermz(i,1)
          kbz2 = ipermz(i,2)

          disx = (cord(kbx2,1)-cord(kbx1,1))/2
          disy = (cord(kby2,2)-cord(kby1,2))/2
          disz = (cord(kbz2,3)-cord(kbz1,3))/2
        
c frac_flg = 0 (no fracture)
c          = 1 (frac - xz plane, vertical   -> kx,kz)
c          = 2 (frac - yz plane, vertical   -> ky,kz)
c          = 3 (frac - xy plane, horizontal -> kx,ky)
c          = 4 (frac - xz and xy plane 1&3  -> kx*,ky,kz)
c          = 5 (frac - yz and xy plane 2&3  -> kx,ky*,kz)
c          = 6 (frac - xz and yz plane 1&2  -> kx,ky,kz*)
c          = 7 (frac - all planes           -> kx*,ky*,kz*)	   

c********************************** check rock failure
          esxi = str_x(i)-phi(i)
          esyi = str_y(i)-phi(i)
          eszi = str_z(i)-phi(i)

c********************************** Frac on xz-plane, frac_flg = 1
          if(esyi.lt.ipmd8_tns.or.
     &    ((esyi.gt.esxi).and.(esyi.gt.eszi).and.
     &     ((esyi.gt.sh_cohe+eszi*ipmd8_clb).or.
     &     (esyi.gt.sh_cohe+esxi*ipmd8_clb)))) then
c          if(esyi.lt.ipmd8_tns) then
c Save stresses when xz-plane frac, normal to y(=2), is initiated for futher dilation calc            
            if(frac_flg(i).eq.0.or.frac_flg(i).eq.3.or.frac_flg(i).eq.2
     &         .or.frac_flg(i).eq.5) then
              if(esyi.lt.ipmd8_tns) then
                es_f_y0(i,2) = ipmd8_tns
                frc_zen(i,2) = 0.
                frc_azm(i,2) = 0.
              elseif(esyi.gt.sh_cohe+eszi*ipmd8_clb) then
                es_f_y0(i,2) = esyi
                frc_zen(i,2) = cohes
                frc_azm(i,2) = 0.
              elseif(esyi.gt.sh_cohe+esxi*ipmd8_clb) then
                es_f_y0(i,2) = esyi
                frc_zen(i,2) = 0.
                frc_azm(i,2) = cohes
              endif               
              es_f_x0(i,2) = esxi
              es_f_z0(i,2) = eszi
              s_f_xy0(i,2) = str_xy(i)
              s_f_yz0(i,2) = str_yz(i)
              s_f_xz0(i,2) = str_xz(i)
c save frac directions in frac_flg           
              if(frac_flg(i).eq.0) then
                frac_flg(i) = 1
              elseif(frac_flg(i).eq.3) then
                frac_flg(i) = 4
              elseif(frac_flg(i).eq.2) then
                frac_flg(i) = 6
              elseif(frac_flg(i).eq.5) then
                frac_flg(i) = 7
              endif
            endif
          endif

c********************************** Frac on yz-plane, frac_flg = 2         
          if(esxi.lt.ipmd8_tns.or.
     &    ((esxi.gt.esyi).and.(esxi.gt.eszi).and.
     &     ((esxi.gt.sh_cohe+eszi*ipmd8_clb)
     &     .or.(esxi.gt.sh_cohe+esyi*ipmd8_clb)))) then
c          if(esxi.lt.ipmd8_tns) then
c Save stresses when yz-plane frac, normal to x(=1), is initiated for futher dilation calc            
            if(frac_flg(i).eq.0.or.frac_flg(i).eq.3.or.frac_flg(i).eq.1
     &        .or.frac_flg(i).eq.4) then
              if(esxi.lt.ipmd8_tns) then
                es_f_x0(i,1) = ipmd8_tns
                frc_zen(i,1) = 0.
                frc_azm(i,2) = 0.
              elseif(esxi.gt.sh_cohe+eszi*ipmd8_clb) then
                es_f_x0(i,1) = esxi
                frc_zen(i,1) = cohes
                frc_azm(i,1) = 0.
              elseif(esxi.gt.sh_cohe+esyi*ipmd8_clb) then
                es_f_x0(i,1) = esxi
                frc_zen(i,1) = 0.
                frc_azm(i,1) = cohes
              endif               
              es_f_y0(i,1) = esyi
              es_f_z0(i,1) = eszi
              s_f_xy0(i,1) = str_xy(i)
              s_f_yz0(i,1) = str_yz(i)
              s_f_xz0(i,1) = str_xz(i)
c save frac directions in frac_flg            
              if(frac_flg(i).eq.0) then
                frac_flg(i) = 2
              elseif(frac_flg(i).eq.3) then
                frac_flg(i) = 5
              elseif(frac_flg(i).eq.1) then
                frac_flg(i) = 6
              elseif(frac_flg(i).eq.4) then
                frac_flg(i) = 7
              endif
            endif
          endif
          
c********************************** Frac on xy-plane, frac_flg = 3        
          if(eszi.lt.ipmd8_tns.or.
     &    ((eszi.gt.esyi).and.(eszi.gt.esxi).and.
     &     ((eszi.gt.sh_cohe+esyi*ipmd8_clb)
     &     .or.(eszi.gt.sh_cohe+esxi*ipmd8_clb)))) then
c          if(eszi.lt.ipmd8_tns) then
c Save stresses when xy-plane frac, normal to z(=3), is initiated for futher dilation calc 
            if(frac_flg(i).eq.0.or.frac_flg(i).eq.1.or.frac_flg(i).eq.2
     &        .or.frac_flg(i).eq.6) then            
              if(eszi.lt.ipmd8_tns) then
                es_f_z0(i,3) = ipmd8_tns
                frc_zen(i,3) = 0.
                frc_azm(i,3) = 0.
              elseif(eszi.gt.sh_cohe+esyi*ipmd8_clb) then
                es_f_z0(i,3) = eszi
                frc_zen(i,3) = cohes
                frc_azm(i,3) = 90.
              elseif(eszi.gt.sh_cohe+esxi*ipmd8_clb) then
                es_f_z0(i,3) = eszi
                frc_zen(i,3) = cohes
                frc_azm(i,3) = 0.
              endif               
              es_f_x0(i,3) = esxi
              es_f_y0(i,3) = esyi
              s_f_xy0(i,3) = str_xy(i)
              s_f_yz0(i,3) = str_yz(i)
              s_f_xz0(i,3) = str_xz(i)
            
              if(frac_flg(i).eq.0) then
                frac_flg(i) = 3
              elseif(frac_flg(i).eq.1) then
                frac_flg(i) = 4
              elseif(frac_flg(i).eq.2) then
                frac_flg(i) = 5
              elseif(frac_flg(i).eq.6) then
                frac_flg(i) = 7
              endif
            endif         
          endif 

c frac_flg = 1, xz-plane -> kx, kz          
          if(frac_flg(i).eq.1) then
            dsx12  = es_f_x0(i,2) - esxi
            dsy12  = es_f_y0(i,2) - esyi
            dsz12  = es_f_z0(i,2) - eszi
            dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
            dsyz12 = abs(s_f_yz0(i,2)-str_yz(i))
            dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
            L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L13 = -sin(frc_zen(i,2)/180*pi)
            L21 = -sin(frc_azm(i,2)/180*pi)
            L22 = cos(frc_azm(i,2)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L33 = cos(frc_zen(i,2)/180*pi)
			            
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
            
     	     ipmd8_nmy = ipmd8_dsy/kny + disy
     &     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
     
           if(esyi.lt.0) then
             ipmd8_shy = 0.
           else
             ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
           endif
            
            ipmd8_nmy = max(0.,ipmd8_nmy)
            ipmd8_shy = max(0.,ipmd8_shy)
            
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)

c frac_flg = 2, yz-plane -> ky, kz            
          elseif(frac_flg(i).eq.2) then
            dsx12  = es_f_x0(i,1) - esxi
            dsy12  = es_f_y0(i,1) - esyi
            dsz12  = es_f_z0(i,1) - eszi
            dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
            dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))	
            L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L13 = -sin(frc_zen(i,1)/180*pi)
            L21 = -sin(frc_azm(i,1)/180*pi)
            L22 = cos(frc_azm(i,1)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L33 = cos(frc_zen(i,1)/180*pi)
            
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12     
                               
           ipmd8_nmx = ipmd8_dsx/knx + disx
     &     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
       
           if(esxi.lt.0) then
             ipmd8_shx = 0.
           else
             ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
           endif
           
            ipmd8_nmx = max(0.,ipmd8_nmx)
            ipmd8_shx = max(0.,ipmd8_shx)
            
            pny(i) = pny0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)

          
c frac_flg = 3, xy-plane -> kx, ky     
          elseif(frac_flg(i).eq.3) then
            dsx12  = es_f_x0(i,3) - esxi
            dsy12  = es_f_y0(i,3) - esyi
            dsz12  = es_f_z0(i,3) - eszi         
            dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
            dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
            
  			L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L13 = -sin(frc_zen(i,3)/180*pi)
            L21 = -sin(frc_azm(i,3)/180*pi)
            L22 = cos(frc_azm(i,3)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                       
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
            L33 = cos(frc_zen(i,3)/180*pi)          
	
	     ipmd8_nmz = ipmd8_dsz/knz + disz
     &     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
      
           if(eszi.lt.0) then
             ipmd8_shz = 0.
           else
             ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz
     &       +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
           endif

            ipmd8_nmz = max(0.,ipmd8_nmz)
            ipmd8_shz = max(0.,ipmd8_shz)
                       
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
            pny(i) = pnz0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)

c frac_flg = 4, xz and xy plane -> kx*, ky, kz   
          elseif(frac_flg(i).eq.4) then
c
c frac on xz
c         
            dsx12  = es_f_x0(i,2) - esxi
            dsy12  = es_f_y0(i,2) - esyi
            dsz12  = es_f_z0(i,2) - eszi
            dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
            dsyz12 = abs(s_f_yz0(i,2)-str_yz(i))     
            dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
            L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L13 = -sin(frc_zen(i,2)/180*pi)
            L21 = -sin(frc_azm(i,2)/180*pi)
            L22 = cos(frc_azm(i,2)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L33 = cos(frc_zen(i,2)/180*pi)
            	           
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
     
	     ipmd8_nmy = ipmd8_dsy/kny + disy
     &     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(esyi.lt.0) then
             ipmd8_shy = 0.
           else
             ipmd8_shy = ((disx/e3(i)+1/ksy)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
           endif
           
c
c frac on xy plane
c 
            dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
            dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
            dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
            dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
            dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
			L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L13 = -sin(frc_zen(i,3)/180*pi)
            L21 = -sin(frc_azm(i,3)/180*pi)
            L22 = cos(frc_azm(i,3)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L33 = cos(frc_zen(i,3)/180*pi)
            
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12   
               	
	     ipmd8_nmz = ipmd8_dsz/knz + disz
     &     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
      
           if(eszi.lt.0) then
             ipmd8_shz = 0.
           else
             ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz
     &       +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
           endif

            ipmd8_nmy = max(0.,ipmd8_nmy)
            ipmd8_shy = max(0.,ipmd8_shy)
            ipmd8_nmz = max(0.,ipmd8_nmz)
            ipmd8_shz = max(0.,ipmd8_shz)
                          
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
     &               +1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
            pny(i) = pny0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
  
c frac_flg = 5, yz and xy plane -> kx, ky*, kz   
          elseif(frac_flg(i).eq.5) then
c          
c frac on yz plane
c
            dsx12  = ipmd8_tns - esxi
            dsy12  = es_f_y0(i,1) - esyi
            dsz12  = es_f_z0(i,1) - eszi
            dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
            dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
            L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L13 = -sin(frc_zen(i,1)/180*pi)
            L21 = -sin(frc_azm(i,1)/180*pi)
            L22 = cos(frc_azm(i,1)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L33 = cos(frc_zen(i,1)/180*pi)
            
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12                         
	
           ipmd8_nmx = ipmd8_dsx/knx + disx
     &     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
       
           if(esxi.lt.0) then 
             ipmd8_shx = 0.
           else
             ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksy)*ipmd8_dsxz)*tan(phid/180*pi)
           endif

c
c frac on xy plane
c     
            dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
            dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
            dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
            dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
            dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
			L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L13 = -sin(frc_zen(i,3)/180*pi)
            L21 = -sin(frc_azm(i,3)/180*pi)
            L22 = cos(frc_azm(i,3)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L33 = cos(frc_zen(i,3)/180*pi)
            
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
            	
           ipmd8_nmz = ipmd8_dsz/knz + disz
     &     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
      
           if(eszi.lt.0) then
             ipmd8_shz = 0.
           else
             ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz
     &       +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)	
           endif
                 
            ipmd8_nmx = max(0.,ipmd8_nmx)
            ipmd8_shx = max(0.,ipmd8_shx)
            ipmd8_nmz = max(0.,ipmd8_nmz)
            ipmd8_shz = max(0.,ipmd8_shz)
      
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
            pny(i) = pny0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
     &               +1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)

c frac_flg = 6, xz and yz plane -> kx, ky, kz*   
          elseif(frac_flg(i).eq.6) then 
c
c frac on yz plane
c                     
            dsx12  = es_f_x0(i,1) - esxi
            dsy12  = es_f_y0(i,1) - esyi
            dsz12  = es_f_z0(i,1) - eszi
            dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
            dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))	
            dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
            L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L13 = -sin(frc_zen(i,1)/180*pi)
            L21 = -sin(frc_azm(i,1)/180*pi)
            L22 = cos(frc_azm(i,1)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L33 = cos(frc_zen(i,1)/180*pi)   
                     
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12    
                          
           ipmd8_nmx = ipmd8_dsx/knx + disx
     &     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
       
           if(esxi.lt.0) then
             ipmd8_shx = 0.
           else
             ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
           endif

c
c frac on xy plane
c            
            dsx12  = es_f_x0(i,2) - (str_x(i)-phi(i))
            dsy12  = es_f_y0(i,2) - (str_y(i)-phi(i))
            dsz12  = es_f_z0(i,2) - (str_z(i)-phi(i))
            dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
            dsyz12 = abs(s_f_yz0(i,2)-str_yz(i)) 
            dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
            L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L13 = -sin(frc_zen(i,2)/180*pi)
            L21 = -sin(frc_azm(i,2)/180*pi)
            L22 = cos(frc_azm(i,2)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L33 = cos(frc_zen(i,2)/180*pi)
                       
            ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12
                  
	     ipmd8_nmy = ipmd8_dsy/kny + disy
     &     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(esyi.lt.0) then
             ipmd8_shy = 0.
           else
             ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
           endif
            
            ipmd8_nmx = max(0.,ipmd8_nmx)
            ipmd8_shx = max(0.,ipmd8_shx)
            ipmd8_nmy = max(0.,ipmd8_nmy)
            ipmd8_shy = max(0.,ipmd8_shy)
    
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
            pny(i) = pny0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
     &               +1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                 
c frac_flg = 7, all plane -> kx*, ky*, kz*   
          elseif(frac_flg(i).eq.7) then  
c
c frac on yz plane
c          
            dsx12  = es_f_x0(i,1) - esxi
            dsy12  = es_f_y0(i,1) - esyi
            dsz12  = es_f_z0(i,1) - eszi
            dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
            dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))	
            dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
            L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L13 = -sin(frc_zen(i,1)/180*pi)
            L21 = -sin(frc_azm(i,1)/180*pi)
            L22 = cos(frc_azm(i,1)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
            L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
            L33 = cos(frc_zen(i,1)/180*pi)
            
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
            
           ipmd8_nmx = ipmd8_dsx/knx + disx
     &     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
       
           if(esxi.lt.0) then
             ipmd8_shx = 0.
           else
             ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy
     &       +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
           endif

c
c frac on xz plane
c     
            
            dsx12  = es_f_x0(i,2) - (str_x(i)-phi(i))
            dsy12  = es_f_y0(i,2) - (str_y(i)-phi(i))
            dsz12  = es_f_z0(i,2) - (str_z(i)-phi(i))
            dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
            dsyz12 = abs(s_f_yz0(i,2)-str_yz(i)) 
            dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
            L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L13 = -sin(frc_zen(i,2)/180*pi)
            L21 = -sin(frc_azm(i,2)/180*pi)
            L22 = cos(frc_azm(i,2)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
            L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
            L33 = cos(frc_zen(i,2)/180*pi)
            
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
                 	           
	     ipmd8_nmy = ipmd8_dsy/kny + disy
     &     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
           
           if(esyi.lt.0) then
             ipmd8_shy = 0.
           else
             ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy
     &      +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
           endif

c
c frac on xy plane
c
            dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
            dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
            dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
            dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
            dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
            dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
			L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L13 = -sin(frc_zen(i,3)/180*pi)
            L21 = -sin(frc_azm(i,3)/180*pi)
            L22 = cos(frc_azm(i,3)/180*pi)
            L23 = 0.
            L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
            L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
            L33 = cos(frc_zen(i,3)/180*pi)   
                    
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12
     &      + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
     
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12
     &      + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
     
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12
     &      + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
     
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12
     &      + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12
     &      + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12
     &      + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12
     &      + (L21*L33+L23*L31)*dsxz12
     
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12
     &      + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12
     &      + (L11*L33+L13*L31)*dsxz12  
     	
           ipmd8_nmz = ipmd8_dsz/knz + disz
     &     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i)
     &     +alp(i)*(t(i)-tini(i)))
      
           if(eszi.lt.0) then
             ipmd8_shz = 0.
           else
             ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz
     &       +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
           endif

            ipmd8_nmx = max(0.,ipmd8_nmx)
            ipmd8_shx = max(0.,ipmd8_shx)
            ipmd8_nmy = max(0.,ipmd8_nmy)
            ipmd8_shy = max(0.,ipmd8_shy)
            ipmd8_nmz = max(0.,ipmd8_nmz)
            ipmd8_shz = max(0.,ipmd8_shz)                
            
            pnx(i) = pnx0(i)+
     &               1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
     &               +1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
            pny(i) = pny0(i)+
     &               1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
     &               +1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
            pnz(i) = pnz0(i)+
     &               1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
     &               +1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy) 
     

          endif          
            pnx(i) = min(perx_m*pnx0(i),pnx(i))
            pny(i) = min(pery_m*pny0(i),pny(i))
            pnz(i) = min(perz_m*pnz0(i),pnz(i))

c            pnx(i) = max(pnx0(i),pnx(i))
c            pny(i) = max(pny0(i),pny(i))  
c            pnz(i) = max(pnz0(i),pnz(i))
 
                            
                  endif 
                  

c...........................................................
               elseif(ispmd.eq.11) then
c s kelkar June 9 2009, simple gangi (1978, Int.J.Rock Mech)
c  model in 2-D
c assumed model in x-y plane and fracture in x-z plane, so only 
c x-perm is modified. This formulation is for fracture faces
c in contact, so str_y = effective stress in y-dir has to be
c compressive
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
     &    .or.ipermstr8.ne.0) then
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
            if(iout.ne.0)
     &          write(iout,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max
*******************************
            
c     &  write(iout,201)ipermx_max,str_y(ipermx_max)-phi(ipermx_max),
c     &        frac_flg(ipermx_max)
c*******************************
            if(iptty.ne.0)
     &         write(iptty,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max 
            if(iout.ne.0)
     &         write(iout,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            
c     &  write(iout,106)ipermy_max,
c     &       estr_x0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       estr_y0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       pny(ipermx_max)*1.e-6,
c     &       cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            if(iptty.ne.0)
     &         write(iptty,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
            if(icnl.eq.0) then
               if(iout.ne.0)
     &          write(iout,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &              cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max  
               if(iptty.ne.0)
     &         write(iptty,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
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
      
