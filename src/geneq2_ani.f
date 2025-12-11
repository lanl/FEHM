      subroutine geneq2_ani(i)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To generate equations isothermal air-water solution at each node.
!D1 anisotropic permeability
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    E!D
!D2 Date         Programmer     Number  Comments
!D2
!D2 Initial implementation: 01-May-04, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/geneq2_ani.f_a  $
!D2
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
!D3
!**********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************
c
c generate equations for 3-d schemes'isothermal air-water ,
c full derivatives, anisotropic permeability'
c
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i
      integer icd
      integer ii1
      integer ii2
      integer idg 
      integer iq
      integer jmi
      integer jml
      integer jmia
      integer jm
      integer neqp1
      integer ij
      integer ij1
      integer ij2
      integer iz
      integer kb
      integer neighc
      integer iau
      integer ial
      integer kz
      integer nmatavw
      real*8 sx1d
      real*8 axi
      real*8 ayi
      real*8 azi
      real*8 alxi
      real*8 alyi
      real*8 alzi
      real*8 avxi
      real*8 avyi
      real*8 avzi
      real*8 pvii
      real*8 phii
      real*8 swi
      real*8 dlpi
      real*8 dlpkb
      real*8 dvpi
      real*8 dvpkb
      real*8 dlei
      real*8 dvei
      real*8 dlekb
      real*8 dvekb
      real*8 dlapi
      real*8 dvapi
      real*8 dlapkb
      real*8 dvapkb
      real*8 dili
      real*8 dilkb
      real*8 divi
      real*8 divkb
      real*8 axkb
      real*8 aykb
      real*8 azkb
      real*8 alxkb
      real*8 alykb
      real*8 alzkb
      real*8 sx2c
      real*8 sx4d
c     real*8 sxzc
      real*8 pvikb
      real*8 phikb
      real*8 radi
      real*8 radkb
      real*8 fid
      real*8 fid1
      real*8 axyd
      real*8 axy
      real*8 axyf
      real*8 heatc
      real*8 pvxy
      real*8 pxy
      real*8 pxyh
      real*8 pxyi
      real*8 sx4h
      real*8 vxyd
      real*8 vxy
      real*8 vxyf
      real*8 dpvti
      real*8 dilpi
      real*8 dilei
      real*8 divpi
      real*8 divei
      real*8 dilpkb
      real*8 divpkb
      real*8 divekb
      real*8 dilekb
      real*8 sx2t
      real*8 sx3t
      real*8 sxzt
      real*8 dlaei
      real*8 dlaekb
      real*8 dvaei
      real*8 dvaekb
      real*8 ti
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8 grav_air
      real*8 grav_term

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    
      integer i1,i2,j,kb1,kb2,kz1,kz2,iface
      integer kb_adv_x,kb_adv_y,kb_adv_z
      integer ij_x,ij_y,ij_z,ji_x,ji_y,ji_z 
      integer iqx,iqy,iqz,kc,kd,jj,i3,i4 
      integer iqxm,iqym,iqzm     
      integer iposx,iposy,iposz

      real*8 flux_x,flux_y,flux_z,flux_dir(3)                
      real*8 dfxpi,dfxpkb,dfxpkb1,dfxpkb2
      real*8 dfypi,dfypkb,dfypkb1,dfypkb2
      real*8 dfzpi,dfzpkb,dfzpkb1,dfzpkb2
      real*8 dfxpkb1a,dfxpkb2a,dfypkb1a,dfypkb2a
      real*8 dfzpkb1a,dfzpkb2a
      real*8 pvikb1,phikb1,pvikb2,phikb2 
      real*8 daxyfpi ,daxyfpkb 
      real*8 sumx,sumy,sumz,termx,termy,termz   
      real*8 presd1,presd2,cordd1,cordd2  
      real*8 sumxg,sumyg,sumzg,dsumxp,dsumyp,dsumzp 
      real*8 dend,ddendpi,ddendpkb  
c changed by avw -- entered here by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif
      if(icons.le.abs(maxit)) then
       grav_air=0.0
      else
       grav_air=grav
      endif

c
c zero out temporary storage
c
      if(.not.allocated(dum_ani)) then
       allocate(dum_ani(neq))
       allocate(dumx_ani(neq))
       allocate(dumy_ani(neq))
       allocate(dumz_ani(neq))
       allocate(axy_ani(neq))

       allocate(it4a(200))
       allocate(it5a(200))
       allocate(it6a(200))     
      endif
c
c storage for upwind
c
      sx1d=sx1(i)
      axi=pnx(i)
      ayi=pny(i)
      azi=pnz(i)
      alxi=axi
      avxi=axi
      alyi=ayi
      avyi=ayi
      alzi=azi
      avzi=azi
      pvii=phi(i)
      dili=dil(i)
      dilpi=dilp(i)
      if(irdof.ne.13) then
        phii=pvii-pcp(i)
        dpvti=dpcef(i)
        divi=div(i)
        divpi=divp(i)
        divei=dive(i)
        dilei=dile(i)
      else
        phii=pvii
      endif
      ti=t(i)
c      swi=s(i)
c
c form constants for i>neq
c
      if(i.gt.neq.and.idualp.eq.0) then
         icd=neq
      else
         icd=0
      endif
      iz=i-icd
c
      iz4m1 = 4*(iz-1)+1
c
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
c define diagonal term for connectivity
      jmi=nelmdg(i-icd)
c define diagonal term for jacobian matrix 
      jmia=jmi-neqp1

c Take care of source/sink term by putting it into empty slot (from 
c same node to same node). 
c If this is an isothermal air-water simulation
      a_axy(jmia+nmatavw)=sk(i)
c Take care of source/sink term 
c If this is an isothermal air-water simulation
c      a_vxy(jmia+nmatavw)=qh(i)

c
c   anisotropy coding
c
c
c create correspondence with array connectivity
c
         kb_adv_x = ncon_adv(iz,1)
         kb_adv_y = ncon_adv(iz,2)
         kb_adv_z = ncon_adv(iz,3)
c         
c   x face terms
      ii1=icxani(i-icd)+1
      ii2=icxani(i-icd+1)
      iqx=0
      do  jm=ii1,ii2   
       iqx=iqx+1
       kb1 = ncon_x1(jm)
       kb2 = ncon_x2(jm)
       termx = sx_x(jm)
       it4(iqx)= kb1
       it4a(iqx)= kb2
       t1(iqx)=termx
      enddo
c   y face terms
      ii1=icyani(i-icd)+1
      ii2=icyani(i-icd+1)
      iqy=0
      do  jm=ii1,ii2   
       iqy=iqy+1
       kb1 = ncon_y1(jm)
       kb2 = ncon_y2(jm)
       termy = sx_y(jm)
       it5(iqy)=kb1
       it5a(iqy)=kb2
       t2(iqy)=termy
      enddo
c   z face terms
c
c   note using -termz
c
      ii1=iczani(i-icd)+1
      ii2=iczani(i-icd+1)
      iqz=0
      do  jm=ii1,ii2   
       iqz=iqz+1
       kb1 = ncon_z1(jm)
       kb2 = ncon_z2(jm)
       termz = sx_z(jm)
       it6(iqz)=kb1
       it6a(iqz) = kb2
       t3(iqz)=termz
      enddo

c note must search below diagonal because + face node can be < iz 
c figure this out once above
           i1 = nelm(iz)+1
           i2 = nelm(iz+1)
           iq = 0
           do jj= i1,i2
            kb = nelm(jj)
            iq = iq +1
             it8(iq) = kb
             it9(iq) = jj - neqp1
           enddo
c  find -x positions
          iqxm = 0
          if(kb_adv_x.gt.0) then
           i3 = nelm(kb_adv_x)+1
           i4 = nelm(kb_adv_x+1) 
            do kd = i3,i4
             iqxm = iqxm +1
             kc = nelm(kd)
              it11(iqxm) = kc
              it11a(iqxm) = kd - neqp1
              if(kc.eq.iz) iposx = kd - neqp1
           enddo
          endif
c  find -y positions
          iqym = 0
          if(kb_adv_y.gt.0) then
           i3 = nelm(kb_adv_y)+1
           i4 = nelm(kb_adv_y+1) 
            do kd = i3,i4
             iqym = iqym +1
             kc = nelm(kd)
              it12(iqym) = kc
              it12a(iqym) = kd - neqp1
              if(kc.eq.iz) iposy = kd - neqp1
           enddo
          endif
c  find -z positions
          iqzm = 0
          if(kb_adv_z.gt.0) then
           i3 = nelm(kb_adv_z)+1
           i4 = nelm(kb_adv_z+1) 
            do kd = i3,i4
             iqzm = iqzm +1
             kc = nelm(kd)
              it13(iqzm) = kc
              it13a(iqzm) = kd - neqp1
              if(kc.eq.iz) iposz = kd - neqp1
           enddo
          endif

       do jm = 1,iq
        kb = it8(jm)
        dum_ani(kb)= 0.0
        axy_ani(kb)= 0.0
       enddo
       do jm = 1,iqxm
         kb = it11(jm)
         dumx_ani(kb) = 0.0
       enddo
       do jm = 1,iqym
        kb = it12(jm)     
        dumy_ani(kb) = 0.0
       enddo
       do jm = 1,iqzm
        kb = it13(jm)
        dumz_ani(kb) = 0.0
       enddo
c
c 3-d geometry
c
      if(icnl.eq.0) then
c
      if(irdof.ne.11) then
c     
c liquid phase calculations
c
         flux_dir = 0.0
c  check x direction upwinding
        if(kb_adv_x.gt.0) then
         sumx = 0.0
         sumxg = 0.0
         do jm=1,iqx
            kb1=it4(jm)
            kz1=kb1-icd
            kb2=it4a(jm)
            kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)            
            cordd1 = cord(kz1,igrav)
            cordd2 = cord(kz2,igrav)
            sumx = sumx + t1(jm)*(presd1-presd2) 
            sumxg = sumxg + t1(jm)*(cordd1-cordd2) 
         enddo
         flux_dir(1) = sumx - 
     &        0.5*grav*(rolf(iz)+rolf(kb_adv_x))*sumxg
        endif
c  check y direction upwinding
        if(kb_adv_y.gt.0) then
         sumy = 0.0
         sumyg = 0.0
         do jm=1,iqy
            kb1=it5(jm)
            kz1=kb1-icd
            kb2=it5a(jm)
            kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)  
            cordd1 = cord(kz1,igrav)
            cordd2 = cord(kz2,igrav)
            sumy = sumy + t2(jm)*(presd1-presd2) 
            sumyg = sumyg + t2(jm)*(cordd1-cordd2) 
         enddo
         flux_dir(2) = sumy - 
     &        0.5*grav*(rolf(iz)+rolf(kb_adv_y))*sumyg
        endif
c  check z direction upwinding
        if(kb_adv_z.gt.0) then
         sumz = 0.0
         sumzg = 0.0
         do jm=1,iqz
            kb1=it6(jm)
            kz1=kb1-icd
            kb2=it6a(jm)
            kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)  
            cordd1 = cord(kz1,igrav)
            cordd2 = cord(kz2,igrav)
            sumz = sumz + t3(jm)*(presd1-presd2) 
            sumzg = sumzg + t3(jm)*(cordd1-cordd2) 
         enddo
         flux_dir(3) = sumz -
     &        0.5*grav*(rolf(iz)+rolf(kb_adv_z))*sumzg
        endif
c
c determine upwind nodes and if liquid phase exists
c
         isl=1
         do iface=1,3
c
c coding to save upwind position
c 
          if(iad.le.iad_up) then
            fid=0.5
            axyd=flux_dir(iface)
            if(axyd.lt.0.0) fid=dnwgt
            if(axyd.gt.0.0) fid=upwgt
            t9(iface)=fid
c
            call setbit(nbits,iface,upwind_l(iz4m1),fid)
c
          else
           if(bit(nbits,iface,upwind_l(iz4m1))) then 
            t9(iface )=1.0

           else
            t9(iface )=0.0

           endif
          endif
         enddo
c
c form equations
c
         if(isl.ne.0) then
c     + x face
          if(kb_adv_x.ne.0) then
           dilkb=dil(kb_adv_x)
           dilpkb=dilp(kb_adv_x)
           fid=t9(1)
           fid1=1.0-fid
           axyf=(fid*dilkb+fid1*dili)
           dend=-(rolf(iz)+rolf(kb_adv_x))*0.5*grav
           ddendpi = -0.5*grav*dglp(iz)
           ddendpkb = -0.5*grav*dglp(kb_adv_x)
           daxyfpi = fid1*dilpi
           daxyfpkb = fid*dilpkb
            do  jm=1,iqx
               kb1=it4(jm)
               kz1=kb1-icd
               kb2=it4a(jm)
               kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)  
               dsumxp = t1(jm)
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumxp
               dumx_ani(kb1) = dumx_ani(kb1) - axyf*dsumxp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumxp
               dumx_ani(kb2) = dumx_ani(kb2) + axyf*dsumxp
            enddo
           flux_x = axyf*(sumx + dend*sumxg)
           dfxpi = daxyfpi*(sumx + dend*sumxg) + 
     &                  axyf*ddendpi*sumxg
           dfxpkb = daxyfpkb*(sumx + dend*sumxg) + 
     &                   axyf*ddendpkb*sumxg
           dum_ani(iz) = dum_ani(iz) + dfxpi
           dum_ani(kb_adv_x) = dum_ani(kb_adv_x) + dfxpkb 
           dumx_ani(iz) = dumx_ani(iz) - dfxpi
           dumx_ani(kb_adv_x) = dumx_ani(kb_adv_x) - dfxpkb 
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_x
           bp(kb_adv_x+nrhs(1))=bp(kb_adv_x+nrhs(1))-flux_x
           axy_ani(kb_adv_x) = axy_ani(kb_adv_x) + flux_x
           a_axy(iposx) = a_axy(iposx) - flux_x
          endif
c     + y face
          if(kb_adv_y.ne.0) then
           dilkb=dil(kb_adv_y)
           dilpkb=dilp(kb_adv_y)
           fid=t9(2)
           fid1=1.0-fid
           axyf=(fid*dilkb+fid1*dili)
           dend= -(rolf(iz)+rolf(kb_adv_y))*0.5*grav
           ddendpi = -0.5*grav*dglp(iz)
           ddendpkb = -0.5*grav*dglp(kb_adv_y)
           daxyfpi = fid1*dilpi
           daxyfpkb = fid*dilpkb
            do  jm=1,iqy
               kb1=it5(jm)
               kz1=kb1-icd
               kb2=it5a(jm)
               kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)  
               dsumyp = t2(jm)
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumyp
               dumy_ani(kb1) = dumy_ani(kb1) - axyf*dsumyp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumyp
               dumy_ani(kb2) = dumy_ani(kb2) + axyf*dsumyp
            enddo
           flux_y = axyf*(sumy + dend*sumyg)
           dfypi = daxyfpi*(sumy + dend*sumyg) + 
     &                  axyf*ddendpi*sumyg
           dfypkb = daxyfpkb*(sumy + dend*sumyg) + 
     &                   axyf*ddendpkb*sumyg
           dum_ani(iz) = dum_ani(iz) + dfypi
           dum_ani(kb_adv_y) = dum_ani(kb_adv_y) + dfypkb 
           dumy_ani(iz) = dumy_ani(iz) - dfypi
           dumy_ani(kb_adv_y) = dumy_ani(kb_adv_y) - dfypkb                 
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_y
           bp(kb_adv_y+nrhs(1))=bp(kb_adv_y+nrhs(1))-flux_y
           axy_ani(kb_adv_y) = axy_ani(kb_adv_y) + flux_y
           a_axy(iposy) = a_axy(iposy) - flux_y
          endif
c     + z face
          if(kb_adv_z.ne.0) then
           dilkb=dil(kb_adv_z)
           dilpkb=dilp(kb_adv_z)
           fid=t9(3)
           fid1=1.0-fid
           axyf=(fid*dilkb+fid1*dili)
           dend= -(rolf(iz)+rolf(kb_adv_z))*0.5*grav
           ddendpi = -0.5*grav*dglp(iz)
           ddendpkb = -0.5*grav*dglp(kb_adv_z)
           daxyfpi = fid1*dilpi
           daxyfpkb = fid*dilpkb
            do  jm=1,iqz
               kb1=it6(jm)
               kz1=kb1-icd
               kb2=it6a(jm)
               kz2=kb2-icd
c            presd1 = phi(kb1)-pcp(kb1)
c            presd2 = phi(kb2)-pcp(kb2)
            presd1 = phi(kb1)
            presd2 = phi(kb2)  
               dsumzp = t3(jm)
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumzp
               dumz_ani(kb1) = dumz_ani(kb1) - axyf*dsumzp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumzp
               dumz_ani(kb2) = dumz_ani(kb2) + axyf*dsumzp
            enddo
           flux_z = axyf*(sumz + dend*sumzg)
           dfzpi = daxyfpi*(sumz + dend*sumzg) + 
     &                  axyf*ddendpi*sumzg
           dfzpkb = daxyfpkb*(sumz + dend*sumzg) + 
     &                   axyf*ddendpkb*sumzg
           dum_ani(iz) = dum_ani(iz) + dfzpi
           dum_ani(kb_adv_z) = dum_ani(kb_adv_z) + dfzpkb 
           dumz_ani(iz) = dumz_ani(iz) - dfzpi
           dumz_ani(kb_adv_z) = dumz_ani(kb_adv_z) - dfzpkb   
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_z
           bp(kb_adv_z+nrhs(1))=bp(kb_adv_z+nrhs(1))-flux_z
           axy_ani(kb_adv_z) = axy_ani(kb_adv_z) + flux_z
           a_axy(iposz) = a_axy(iposz) - flux_z
          endif
         endif
c
c  now load jacobian matrix (node i)
c
       do jm = 1,iq
        kb = it8(jm)
        jj = it9(jm)
        a(jj) = a(jj) + dum_ani(kb)
        if(kb.ne.iz) a_axy(jj) = a_axy(jj) + axy_ani(kb)
       enddo
c
c  now load jacobian matrix (node kb_adv_x) (negative contribution)
c
       do jm = 1,iqxm
        kb = it11(jm)
        jj = it11a(jm)
         a(jj) = a(jj) + dumx_ani(kb)
       enddo
c
c  now load jacobian matrix (node kb_adv_y) (negative contribution)
c
       do jm = 1,iqym
        kb = it12(jm)
        jj = it12a(jm)
        a(jj) = a(jj) + dumy_ani(kb)
       enddo
c
c  now load jacobian matrix (node kb_adv_z) (negative contribution)
c
       do jm = 1,iqzm
        kb = it13(jm)
        jj = it13a(jm)
        a(jj) = a(jj) + dumz_ani(kb)
       enddo
c
      endif
      endif
      if(irdof.ne.13) then
c     
c vapour phase calculations
c

      endif
      return
      end


