      subroutine geneq1_ani(i)
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
!D1 To generate equations heat and mass transfer solution at each node.
!D1 anisotropic permeability
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    E!D
!D2 Date         Programmer     Number  Comments
!D2
!D2 Initial implementation: 01-april-2010, Programmer: G. Zyvoloski
!D2
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
!D4   FEHM Application Version 3.06
!D4
!**********************************************************************
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
      real*8 sx3c
      real*8 dlaei
      real*8 dlaekb
      real*8 dvaei
      real*8 dvaekb
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8 grav_air
      real*8 grav_term
      real*8  thxi, thxkb, thyi, thykb, thzi, thzkb, ti

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
      
      integer neq1,neq2,neq3,neq4

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
      
      real*8 dfxei,dfxekb,daxyfei,daxyfekb
      real*8 ddendei,ddendekb
      real*8 dfyei,dfyekb
      real*8 dfzei,dfzekb
      
      real*8  enli, enlkb, envi, envkb
      real*8  deli,delei,devi,devei,delkb,delekb,aexyf
      real*8 flux_xe, flux_ye, flux_ze
      real*8 daexyfpi,daexyfpkb,dfxepi,dfxepkb,daexyfei,daexyfekb
      real*8 dfxeei,dfxeekb

c changed by avw -- entered here by seh
      neqp1=neq+1
      neq1 = neq
      neq2 = neq1 + neq
      neq3 = neq2 + neq
      neq4 = neq3 + neq
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
       allocate(dum_ani(neq4))
       allocate(dumx_ani(neq4))
       allocate(dumy_ani(neq4))
       allocate(dumz_ani(neq4))      
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
        phii=pvii-pcp(i)
        dpvti=dpcef(i)
        divi=div(i)
        divpi=divp(i)
        divei=dive(i)
        dilei=dile(i)
      enli=enlf(i)
      deli=delf(i)
      delei=delef(i)
      envi=envf(i)
      devi=devf(i)
      devei=devef(i)
      ti=t(i)
      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      
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
c
c create correspondence with array connectivity
c
         kb_adv_x = ncon_adv(iz,1)
         kb_adv_y = ncon_adv(iz,2)
         kb_adv_z = ncon_adv(iz,3)

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
           enddo
          endif

       do jm = 1,iq
        kb = it8(jm)
        dum_ani(kb)= 0.0
        dum_ani(kb+neq1)= 0.0
        dum_ani(kb+neq2)= 0.0
        dum_ani(kb+neq3)= 0.0
       enddo
        do jm = 1,iqxm
         kb = it11(jm)
         dumx_ani(kb) = 0.0
         dumx_ani(kb+neq1) = 0.0
         dumx_ani(kb+neq2) = 0.0
         dumx_ani(kb+neq3) = 0.0
       enddo
       do jm = 1,iqym
        kb = it12(jm)     
         dumy_ani(kb) = 0.0
         dumy_ani(kb+neq1) = 0.0
         dumy_ani(kb+neq2) = 0.0
         dumy_ani(kb+neq3) = 0.0        
       enddo
       do jm = 1,iqzm
        kb = it13(jm)
         dumz_ani(kb) = 0.0
         dumz_ani(kb+neq1) = 0.0
         dumz_ani(kb+neq2) = 0.0
         dumz_ani(kb+neq3) = 0.0
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
            if (irdof .ne. 13) then
             presd1 = phi(kb1)-pcp(kb1)
             presd2 = phi(kb2)-pcp(kb2)
            else
             presd1 = phi(kb1)
             presd2 = phi(kb2) 
            endif
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
            if (irdof .ne. 13) then
             presd1 = phi(kb1)-pcp(kb1)
             presd2 = phi(kb2)-pcp(kb2)
            else
             presd1 = phi(kb1)
             presd2 = phi(kb2) 
            endif  
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
            if (irdof .ne. 13) then
             presd1 = phi(kb1)-pcp(kb1)
             presd2 = phi(kb2)-pcp(kb2)
            else
             presd1 = phi(kb1)
             presd2 = phi(kb2) 
            endif  
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
           enlkb=enlf(kb_adv_x)
           delkb=delf(kb_adv_x)
	     delekb=delef(kb_adv_x)
           dilekb=dile(kb_adv_x)
           dilekb=dile(kb_adv_x)
           ddendei = -0.5*grav*dgle(iz)
           ddendekb = -0.5*grav*dgle(kb_adv_x)
           aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            do  jm=1,iqx
               kb1=it4(jm)
               kz1=kb1-icd
               kb2=it4a(jm)
               kz2=kb2-icd
               dsumxp = t1(jm)
c mass eq (primary) pressure derivative               
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumxp
               dumx_ani(kb1) = dumx_ani(kb1) - axyf*dsumxp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumxp
               dumx_ani(kb2) = dumx_ani(kb2) + axyf*dsumxp
c energy eq (primary) pressure derivative                              
               dum_ani(kb1+neq2) = dum_ani(kb1+neq2) + aexyf*dsumxp
               dumx_ani(kb1+neq2) = dumx_ani(kb1+neq2) - aexyf*dsumxp
               dum_ani(kb2+neq2) = dum_ani(kb2+neq2) - aexyf*dsumxp
               dumx_ani(kb2+neq2) = dumx_ani(kb2+neq2) + aexyf*dsumxp
            enddo
c mass equation terms            
           flux_x = axyf*(sumx + dend*sumxg)
c flow terms for mass equation, derivatives wrt pressure           
           dfxpi = daxyfpi*(sumx + dend*sumxg) + 
     &                  axyf*ddendpi*sumxg
           dfxpkb = daxyfpkb*(sumx + dend*sumxg) + 
     &                   axyf*ddendpkb*sumxg
           dum_ani(iz) = dum_ani(iz) + dfxpi
           dum_ani(kb_adv_x) = dum_ani(kb_adv_x) + dfxpkb 
           dumx_ani(iz) = dumx_ani(iz) - dfxpi
           dumx_ani(kb_adv_x) = dumx_ani(kb_adv_x) - dfxpkb
c flow terms for mass equation, derivatives wrt energy variable  
            daxyfei = fid1*dilei
            daxyfekb = fid*dilekb            
            dfxei = daxyfei*(sumx + dend*sumxg) +  
     &                  axyf*ddendei*sumxg
            dfxekb = daxyfekb*(sumx + dend*sumxg) + 
     &                   axyf*ddendekb*sumxg  
            dum_ani(iz+neq1) = dum_ani(iz+neq1) + dfxei
            dum_ani(kb_adv_x+neq1) = dum_ani(kb_adv_x+neq1) + dfxekb
            dumx_ani(iz+neq1) = dumx_ani(iz+neq1) - dfxei
            dumx_ani(kb_adv_x+neq1) = dumx_ani(kb_adv_x+neq1) - dfxekb
c add to balance equation
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_x
           bp(kb_adv_x+nrhs(1))=bp(kb_adv_x+nrhs(1))-flux_x
c 
c energy equation terms  
c  
c           enlkb=enlf(kb_adv_x)
c           delkb=delf(kb_adv_x)
c	      delekb=delef(kb_adv_x)
c           dilekb=dile(kb_adv_x)
c           aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
           flux_xe = aexyf*(sumx + dend*sumxg)
           daexyfpi = fid1*(dilpi*enli+dili*deli)
           daexyfpkb = fid*(dilpkb*enlkb+dilkb*delkb)
c flow terms for energy equation, derivatives wrt pressure           
           dfxepi = daexyfpi*(sumx + dend*sumxg) + 
     &                  aexyf*ddendpi*sumxg
           dfxepkb = daexyfpkb*(sumx + dend*sumxg) + 
     &                   aexyf*ddendpkb*sumxg
           dum_ani(iz+neq2) = dum_ani(iz+neq2) + dfxepi
           dum_ani(kb_adv_x+neq2) = dum_ani(kb_adv_x+neq2) + dfxepkb 
           dumx_ani(iz+neq2) = dumx_ani(iz+neq2) - dfxepi
           dumx_ani(kb_adv_x+neq2) = dumx_ani(kb_adv_x+neq2) - dfxepkb
c flow terms for energy equation, derivatives wrt energy variable  
            dilekb=dile(kb_adv_x)
            ddendei = -0.5*grav*dgle(iz)
            ddendekb = -0.5*grav*dgle(kb_adv_x)
            daexyfei = fid1*(dilei*enli+dili*delei)
            daexyfekb = fid*(dilekb*enlkb+dilkb*delekb)       
            dfxeei = daexyfei*(sumx + dend*sumxg) +  
     &                  aexyf*ddendei*sumxg
            dfxeekb = daexyfekb*(sumx + dend*sumxg) + 
     &                   aexyf*ddendekb*sumxg  
            dum_ani(iz+neq3) = dum_ani(iz+neq3) + dfxeei
            dum_ani(kb_adv_x+neq3) = dum_ani(kb_adv_x+neq3) + dfxeekb
            dumx_ani(iz+neq3) = dumx_ani(iz+neq3) - dfxeei
            dumx_ani(kb_adv_x+neq3) = dumx_ani(kb_adv_x+neq3) - dfxeekb
c add to balance equation
           bp(iz+nrhs(2))=bp(iz+nrhs(2))+flux_xe
           bp(kb_adv_x+nrhs(2))=bp(kb_adv_x+nrhs(2))-flux_xe
          endif
c     + y face
          if(kb_adv_y.ne.0) then
           dilkb=dil(kb_adv_y)
           dilpkb=dilp(kb_adv_y)
           fid=t9(2)
           fid1=1.0-fid
           axyf=(fid*dilkb+fid1*dili)
           dend=-(rolf(iz)+rolf(kb_adv_y))*0.5*grav
           ddendpi = -0.5*grav*dglp(iz)
           ddendpkb = -0.5*grav*dglp(kb_adv_y)
           daxyfpi = fid1*dilpi
           daxyfpkb = fid*dilpkb
           enlkb=enlf(kb_adv_y)
           delkb=delf(kb_adv_y)
	     delekb=delef(kb_adv_y)
           dilekb=dile(kb_adv_y)
           dilekb=dile(kb_adv_y)
           ddendei = -0.5*grav*dgle(iz)
           ddendekb = -0.5*grav*dgle(kb_adv_y)
           aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            do  jm=1,iqy
               kb1=it5(jm)
               kz1=kb1-icd
               kb2=it5a(jm)
               kz2=kb2-icd
               dsumyp = t2(jm)
c mass eq (primary) pressure derivative               
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumyp
               dumy_ani(kb1) = dumy_ani(kb1) - axyf*dsumyp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumyp
               dumy_ani(kb2) = dumy_ani(kb2) + axyf*dsumyp
c energy eq (primary) pressure derivative                              
               dum_ani(kb1+neq2) = dum_ani(kb1+neq2) + aexyf*dsumyp
               dumy_ani(kb1+neq2) = dumy_ani(kb1+neq2) - aexyf*dsumyp
               dum_ani(kb2+neq2) = dum_ani(kb2+neq2) - aexyf*dsumyp
               dumy_ani(kb2+neq2) = dumy_ani(kb2+neq2) + aexyf*dsumyp  
            enddo
c mass equation terms            
           flux_y = axyf*(sumy + dend*sumyg)
c flow terms for mass equation, derivatives wrt pressure           
           dfxpi = daxyfpi*(sumy + dend*sumyg) + 
     &                  axyf*ddendpi*sumyg
           dfxpkb = daxyfpkb*(sumy + dend*sumyg) + 
     &                   axyf*ddendpkb*sumyg
           dum_ani(iz) = dum_ani(iz) + dfxpi
           dum_ani(kb_adv_y) = dum_ani(kb_adv_y) + dfxpkb 
           dumy_ani(iz) = dumy_ani(iz) - dfxpi
           dumy_ani(kb_adv_y) = dumy_ani(kb_adv_y) - dfxpkb
c flow terms for mass equation, derivatives wrt energy variable  
c            ddendei = -0.5*grav*dgle(iz)
c            ddendekb = -0.5*grav*dgle(kb_adv_y)   
            daxyfei = fid1*dilei 
            daxyfekb = fid*dilekb            
            dfxei = daxyfei*(sumy + dend*sumyg) +  
     &                  axyf*ddendei*sumyg
            dfxekb = daxyfekb*(sumy + dend*sumyg) + 
     &                   axyf*ddendekb*sumyg  
            dum_ani(iz+neq1) = dum_ani(iz+neq1) + dfxei
            dum_ani(kb_adv_y+neq1) = dum_ani(kb_adv_y+neq1) + dfxekb
            dumy_ani(iz+neq1) = dumy_ani(iz+neq1) - dfxei
            dumy_ani(kb_adv_y+neq1) = dumy_ani(kb_adv_y+neq1) - dfxekb
c add to balance equation
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_y
           bp(kb_adv_y+nrhs(1))=bp(kb_adv_y+nrhs(1))-flux_y
c 
c energy equation terms  
c  
           flux_ye = aexyf*(sumy + dend*sumyg)
           daexyfpi = fid1*(dilpi*enli+dili*deli)
           daexyfpkb = fid*(dilpkb*enlkb+dilkb*delkb)
c flow terms for energy equation, derivatives wrt pressure           
           dfxepi = daexyfpi*(sumy + dend*sumyg) + 
     &                  aexyf*ddendpi*sumyg
           dfxepkb = daexyfpkb*(sumy + dend*sumyg) + 
     &                   aexyf*ddendpkb*sumyg
           dum_ani(iz+neq2) = dum_ani(iz+neq2) + dfxepi
           dum_ani(kb_adv_y+neq2) = dum_ani(kb_adv_y+neq2) + dfxepkb 
           dumy_ani(iz+neq2) = dumy_ani(iz+neq2) - dfxepi
           dumy_ani(kb_adv_y+neq2) = dumy_ani(kb_adv_y+neq2) - dfxepkb
c flow terms for energy equation, derivatives wrt energy variable  
c            dilekb=dile(kb_adv_y)
c            ddendei = -0.5*grav*dgle(iz)
c            ddendekb = -0.5*grav*dgle(kb_adv_y)
            daexyfei = fid1*(dilei*enli+dili*delei)
            daexyfekb = fid*(dilekb*enlkb+dilkb*delekb)       
            dfxeei = daexyfei*(sumy + dend*sumyg) +  
     &                  aexyf*ddendei*sumyg
            dfxeekb = daexyfekb*(sumy + dend*sumyg) + 
     &                   aexyf*ddendekb*sumyg  
            dum_ani(iz+neq3) = dum_ani(iz+neq3) + dfxeei
            dum_ani(kb_adv_y+neq3) = dum_ani(kb_adv_y+neq3) + dfxeekb
            dumy_ani(iz+neq3) = dumy_ani(iz+neq3) - dfxeei
            dumy_ani(kb_adv_y+neq3) = dumy_ani(kb_adv_y+neq3) - dfxeekb
c add to balance equation
           bp(iz+nrhs(2))=bp(iz+nrhs(2))+flux_ye
           bp(kb_adv_y+nrhs(2))=bp(kb_adv_y+nrhs(2))-flux_ye
          endif
c     + z face
          if(kb_adv_z.ne.0) then
           dilkb=dil(kb_adv_z)
           dilpkb=dilp(kb_adv_z)
           fid=t9(3)
           fid1=1.0-fid
           axyf=(fid*dilkb+fid1*dili)
           dend=-(rolf(iz)+rolf(kb_adv_z))*0.5*grav
           ddendpi = -0.5*grav*dglp(iz)
           ddendpkb = -0.5*grav*dglp(kb_adv_z)
           daxyfpi = fid1*dilpi
           daxyfpkb = fid*dilpkb
           enlkb=enlf(kb_adv_z)
           delkb=delf(kb_adv_z)
	     delekb=delef(kb_adv_z)
           dilekb=dile(kb_adv_z)
           dilekb=dile(kb_adv_z)
           ddendei = -0.5*grav*dgle(iz)
           ddendekb = -0.5*grav*dgle(kb_adv_z)
           aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            do  jm=1,iqz
               kb1=it6(jm)
               kz1=kb1-icd
               kb2=it6a(jm)
               kz2=kb2-icd
               dsumzp = t3(jm)
c mass eq (primary) pressure derivative               
               dum_ani(kb1) = dum_ani(kb1) + axyf*dsumzp
               dumz_ani(kb1) = dumz_ani(kb1) - axyf*dsumzp
               dum_ani(kb2) = dum_ani(kb2) - axyf*dsumzp
               dumz_ani(kb2) = dumz_ani(kb2) + axyf*dsumzp
c energy eq (primary) pressure derivative                              
               dum_ani(kb1+neq2) = dum_ani(kb1+neq2) + aexyf*dsumzp
               dumz_ani(kb1+neq2) = dumz_ani(kb1+neq2) - aexyf*dsumzp
               dum_ani(kb2+neq2) = dum_ani(kb2+neq2) - aexyf*dsumzp
               dumz_ani(kb2+neq2) = dumz_ani(kb2+neq2) + aexyf*dsumzp  
            enddo
c mass equation terms            
           flux_z = axyf*(sumz + dend*sumzg)
c flow terms for mass equation, derivatives wrt pressure           
           dfxpi = daxyfpi*(sumz + dend*sumzg) + 
     &                  axyf*ddendpi*sumzg
           dfxpkb = daxyfpkb*(sumz + dend*sumzg) + 
     &                   axyf*ddendpkb*sumzg
           dum_ani(iz) = dum_ani(iz) + dfxpi
           dum_ani(kb_adv_z) = dum_ani(kb_adv_z) + dfxpkb 
           dumz_ani(iz) = dumz_ani(iz) - dfxpi
           dumz_ani(kb_adv_z) = dumz_ani(kb_adv_z) - dfxpkb
c flow terms for mass equation, derivatives wrt energy variable  
c            dilekb=dile(kb_adv_z)
c            ddendei = -0.5*grav*dgle(iz)
c            ddendekb = -0.5*grav*dgle(kb_adv_z)
            daxyfei = fid1*dilei
            daxyfekb = fid*dilekb            
            dfxei = daxyfei*(sumz + dend*sumzg) +  
     &                  axyf*ddendei*sumzg
            dfxekb = daxyfekb*(sumz + dend*sumzg) + 
     &                   axyf*ddendekb*sumzg  
            dum_ani(iz+neq1) = dum_ani(iz+neq1) + dfxei
            dum_ani(kb_adv_z+neq1) = dum_ani(kb_adv_z+neq1) + dfxekb
            dumz_ani(iz+neq1) = dumz_ani(iz+neq1) - dfxei
            dumz_ani(kb_adv_z+neq1) = dumz_ani(kb_adv_z+neq1) - dfxekb
c add to balance equation
           bp(iz+nrhs(1))=bp(iz+nrhs(1))+flux_z
           bp(kb_adv_z+nrhs(1))=bp(kb_adv_z+nrhs(1))-flux_z
c 
c energy equation terms  
c  
c           enlkb=enlf(kb_adv_z)
c           delkb=delf(kb_adv_z)
c	      delekb=delef(kb_adv_z)
c           dilekb=dile(kb_adv_z)
c           aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
           flux_ze = aexyf*(sumz + dend*sumzg)
           daexyfpi = fid1*(dilpi*enli+dili*deli)
           daexyfpkb = fid*(dilpkb*enlkb+dilkb*delkb)
c flow terms for energy equation, derivatives wrt pressure           
           dfxepi = daexyfpi*(sumz + dend*sumzg) + 
     &                  aexyf*ddendpi*sumzg
           dfxepkb = daexyfpkb*(sumz + dend*sumzg) + 
     &                   aexyf*ddendpkb*sumzg
           dum_ani(iz+neq2) = dum_ani(iz+neq2) + dfxepi
           dum_ani(kb_adv_z+neq2) = dum_ani(kb_adv_z+neq2) + dfxepkb 
           dumz_ani(iz+neq2) = dumz_ani(iz+neq2) - dfxepi
           dumz_ani(kb_adv_z+neq2) = dumz_ani(kb_adv_z+neq2) - dfxepkb
c flow terms for energy equation, derivatives wrt energy variable  
            dilekb=dile(kb_adv_z)
            ddendei = -0.5*grav*dgle(iz)
            ddendekb = -0.5*grav*dgle(kb_adv_z)
            daexyfei = fid1*(dilei*enli+dili*delei)
            daexyfekb = fid*(dilekb*enlkb+dilkb*delekb)       
            dfxeei = daexyfei*(sumz + dend*sumzg) +  
     &                  aexyf*ddendei*sumzg
            dfxeekb = daexyfekb*(sumz + dend*sumzg) + 
     &                   aexyf*ddendekb*sumzg  
            dum_ani(iz+neq3) = dum_ani(iz+neq3) + dfxeei
            dum_ani(kb_adv_z+neq3) = dum_ani(kb_adv_z+neq3) + dfxeekb
            dumz_ani(iz+neq3) = dumz_ani(iz+neq3) - dfxeei
            dumz_ani(kb_adv_z+neq3) = dumz_ani(kb_adv_z+neq3) - dfxeekb
c add to balance equation
           bp(iz+nrhs(2))=bp(iz+nrhs(2))+flux_ze
           bp(kb_adv_z+nrhs(2))=bp(kb_adv_z+nrhs(2))-flux_ze
          endif
         endif
c
c  now load jacobian matrix (node i)
c
       do jm = 1,iq
        kb = it8(jm)
        jj = it9(jm)
         a(jj+nmat(1)) = a(jj+nmat(1)) + dum_ani(kb)
         a(jj+nmat(2)) = a(jj+nmat(2)) + dum_ani(kb+neq1)
         a(jj+nmat(3)) = a(jj+nmat(3)) + dum_ani(kb+neq2)
         a(jj+nmat(4)) = a(jj+nmat(4)) + dum_ani(kb+neq3)
       enddo
c
c  now load jacobian matrix (node kb_adv_x) (negative contribution)
c
       do jm = 1,iqxm
        kb = it11(jm)
        jj = it11a(jm)
         a(jj+nmat(1)) = a(jj+nmat(1)) + dumx_ani(kb)
	   a(jj+nmat(2)) = a(jj+nmat(2)) + dumx_ani(kb+neq1)
	   a(jj+nmat(3)) = a(jj+nmat(3)) + dumx_ani(kb+neq2)
         a(jj+nmat(4)) = a(jj+nmat(4)) + dumx_ani(kb+neq3)
       enddo
c
c  now load jacobian matrix (node kb_adv_y) (negative contribution)
c
       do jm = 1,iqym
        kb = it12(jm)
        jj = it12a(jm)
         a(jj+nmat(1)) = a(jj+nmat(1)) + dumy_ani(kb)
	   a(jj+nmat(2)) = a(jj+nmat(2)) + dumy_ani(kb+neq1)
	   a(jj+nmat(3)) = a(jj+nmat(3)) + dumy_ani(kb+neq2)
         a(jj+nmat(4)) = a(jj+nmat(4)) + dumy_ani(kb+neq3)        
       enddo
c
c  now load jacobian matrix (node kb_adv_z) (negative contribution)
c
       do jm = 1,iqzm
        kb = it13(jm)
        jj = it13a(jm)
         a(jj+nmat(1)) = a(jj+nmat(1)) + dumz_ani(kb)
	   a(jj+nmat(2)) = a(jj+nmat(2)) + dumz_ani(kb+neq1)
	   a(jj+nmat(3)) = a(jj+nmat(3)) + dumz_ani(kb+neq2)
         a(jj+nmat(4)) = a(jj+nmat(4)) + dumz_ani(kb+neq3)        
       enddo
c       
c now load accumulation terms
c
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i)
      a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*deef(i)+deqh(i)

      endif
      endif
c
c set up traditional pointers
c
        ii2=nelm(i-icd+1)
        jmi=nelmdg(i-icd)
        jmia=jmi-neqp1
        iq = 0
        do j=jmi+1,ii2
         kb = nelm(j)
         if(kb.le.neq_primary) then
          iq = iq+1
          i3 = nelmdg(kb)+1
          i4 = nelm(kb+1)
          it8(iq) = kb
          it10(iq) = istrw(j-neqp1)
          it11(iq)=j-neqp1
          ij1=nelm(kb)+1
          ij2=nelmdg(kb)-1
          do ij=ij1,ij2
           if(nelm(ij).eq.iz) then
            it12(iq)=ij-neqp1
            go to 1500
           endif
          enddo
         endif
1500    continue 
        enddo
c
c add heat conduction
c
      do jm=1,iq
       iw = it10(jm)
       if(iw.gt.0)then
        kb=it8(jm)
        kz=kb-icd
        neighc=it9(jm)
        iau=it11(jm)
        ial=it12(jm)
        jml=nelmdg(kb)-neqp1
        thxkb=thx(kb)
        thykb=thy(kb)
        thzkb=thz(kb)
        sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
        sx2t=2.*thxi*thxkb/(thxi+thxkb)
        sx3t=2.*thyi*thykb/(thyi+thykb)
        sxzt=2.*thzi*thzkb/(thzi+thzkb)       
        delx2=(cord(kz,1)-cord(iz,1))**2
        dely2=(cord(kz,2)-cord(iz,2))**2
        delz2=(cord(kz,3)-cord(iz,3))**2 
        dis2=delx2+dely2+delz2
        if(dis2.gt.dis_tol.and.iw.gt.0) then
          sx3c=sx2c*dis2/
     &    (delx2/sx2t+dely2/sx3t+
     &    delz2/sxzt)
        else
          sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
        endif       
        heatc=sx3c
        bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatc*(t(kb)-ti)
        bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatc*(t(kb)-ti)
        a(jmia+nmat(3))=a(jmia+nmat(3))-heatc*dtpa(i)
        a(jmia+nmat(4))=a(jmia+nmat(4))-heatc*dtpae(i)
        a(ial+nmat(3))=a(ial+nmat(3))+heatc*dtpa(i)
        a(ial+nmat(4))=a(ial+nmat(4))+heatc*dtpae(i)
        a(iau+nmat(3))=a(iau+nmat(3))+heatc*dtpa(kb)
        a(iau+nmat(4))=a(iau+nmat(4))+heatc*dtpae(kb)
        a(jml+nmat(3))=a(jml+nmat(3))-heatc*dtpa(kb)
        a(jml+nmat(4))=a(jml+nmat(4))-heatc*dtpae(kb)
       endif
      enddo              
      if(irdof.ne.13) then
c     
c vapour phase calculations
c

      endif
      return
      end


