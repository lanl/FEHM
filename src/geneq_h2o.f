      subroutine  geneq_h2o  ( i )
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
CD1
CD1  PURPOSE
CD1
CD1  This subroutine generates the equations for 3-dimensional heat
CD1  and mass transfer for one componet of a mixture(water)
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/geneq_h2o.f_a  $
!D2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use comflow
      use davidi
      use comji
      use comhi
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use commeth
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
      integer imd,iwd
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8  aexy, aexyf, alxi, alxkb, alyi, alykb, alzi, alzkb
      real*8  avxi, avyi, avzi, axi, axkb, axy, axyd, axyf
      real*8  ayi, aykb, azi, azkb
      real*8  delei, delekb, deli, delkb, devei, devekb, devi, devkb 
      real*8  dilei, dilekb, dili, dilkb, dilpi, dilpkb 
      real*8  divei, divekb, divi, divkb, divpi, divpkb
      real*8  dlaei, dlaekb, dlapi, dlapkb, dleei, dleekb 
      real*8  dlei, dlekb, dlepi, dlepkb, dlpi, dlpkb, dpvti  
      real*8  dvaei, dvaekb, dvapi, dvapkb, dveei, dveekb
      real*8  dvei, dvekb, dvepi, dvepkb, dvpi, dvpkb
      real*8  fid, fid1, grav_air, heatc, enli, enlkb, envi, envkb
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyi, pxyh
      real*8  radi, radkb
      real*8  swi, sx1d, sx2c, sx2t, sx3c, sx3t, sx4d, sx4h, sxzt
      real*8  thxi, thxkb, thyi, thykb, thzi, thzkb, ti
      real*8  vexy, vexyf, vxy, vxyd, vxyf

      real*8  divmi, divmkb, dilmi, dilmkb                      
      real*8  dlami, dlamkb, dlemi, dlemkb                      
      real*8  dvami, dvamkb, dvemi, dvemkb                      
      real*8  dleim, dlekbm                              

      real*8  divwi, divwkb, dilwi, dilwkb                      
      real*8  dlawi, dlawkb, dlewi, dlewkb                      
      real*8  dvawi, dvawkb, dvewi, dvewkb, dleihyd                      
      real*8  dleiw, dlekbw, dpvtw, dpvthyd, dlekbhyd                             

      parameter(dis_tol=1.d-12)
c
c
! Guessed at value for grav_air
      grav_air = 0.
c changed by avw 4/95 -- entered into new version by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif

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
c gaz 1-13-04
      phii=pvii-pcp(i)
      dpvti=dpcef(i)
      dpvtw=dpcp3(i)
      dpvthyd=dpcp4(i) 
      enli=enlf(i)
      deli=delf(i)
      delei=delef(i)
      envi=envf(i)
      devi=devf(i)
      devei=devef(i)
      dili=dil(i)
      divi=div(i)
      dilpi=dilp(i)
      dilei=dile(i)
      dilwi=dilw(i)
      dilmi=dilm(i)
      divpi=divp(i)
      divei=dive(i)
      divwi=divw(i)
      divmi=divm(i)
      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      ti=t(i)
      swi=s(i)
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
      iz4m1 = 4*(iz-1) +1
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1

c Take care of source/sink term by putting it into empty slot (from
c same node to same node).
c If this is an isothermal air-water simulation
      a_axy(jmia+nmatavw)=sk(i)*s(i)
c Take care of souce/sink term
c If this is an isothermal air-water simulation
      a_vxy(jmia+nmatavw)=sk(i)*(1-s(i))

      do 58 jm=jmi+1,ii2
      iq=iq+1
      kb=nelm(jm)+icd
      it8(iq)=kb
      it9(iq)=jm-ii1+1
      it10(iq)=istrw(jm-neqp1)
      it11(iq)=jm-neqp1
      ij1=nelm(nelm(jm))+1
      ij2=nelmdg(nelm(jm))-1
      do 68 ij=ij1,ij2
      if(nelm(ij).eq.iz) then
      it12(iq)=ij-neqp1
      endif
68    continue
58    continue
      if(icnl.eq.0) then
c
c 3-d geometry
c
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw=it10(jm)
            iwd=it10(jm)
            iw =abs(iwd)
            axkb=pnx(kb)
            aykb=pny(kb)
            azkb=pnz(kb)
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
            reduction_factor = red_factor(istrw_itfc(it11(jm)))
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
            thxkb=thx(kb)
            thykb=thy(kb)
            thzkb=thz(kb)
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            sxzt=2.*thzi*thzkb/(thzi+thzkb)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)-pcpmeth(kb)
c           pxy=sx2c*perml(1)+sx3c*perml(2)+sxzc*perml(3)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
            endif
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t+
     &           delz2/sxzt)
            else
              sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c                    
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue
      else if ( icnl.ne.0 )  t h e n
c
c 2-d geometry
c
      radi=cord(iz,3)
         do 69 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iwd=it10(jm)
            iw = abs(iwd)
            neighc=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
            reduction_factor = red_factor(istrw_itfc(it11(jm)))
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
            thxkb=thx(kb)
            thykb=thy(kb)
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
c           pxy=sx2c*perml(1)+sx3c*perml(2)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
            else
              sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 69      continue
c
      endif
c
c liquid phase calculations
c
      do 60 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      pxyi=t1(neighc)
      sx4d=t6(neighc)
      axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     **(cord(kz,igrav)-cord(iz,igrav))
      t8(neighc)=axyd
60    continue
c
c determine upwind nodes and if liquid phase exists
c
      isl=1
      do 61 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      fid=.5
c add coding to save upwind position
         if(iad.le.iad_up) then
      axyd=t8(neighc)
      if(axyd.lt.0.0) fid=dnwgt
      if(axyd.gt.0.0) fid=upwgt
      if(t3(neighc).le.0.0) t9(neighc)=fid
      if(t3(neighc).gt.0.0) t9(neighc)=fid
c
            call setbit(nbits,neighc,upwind_l(iz4m1),fid)
c
         else
           if(bit(nbits,neighc,upwind_l(iz4m1))) then 
            t9(neighc)=1.0
           else
            t9(neighc)=0.0
           endif
         endif
61    continue
c
c form equations
c
      if(isl.ne.0) then
      do 62 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      iau=it11(jm)
      ial=it12(jm)
      jml=nelmdg(kz)-neqp1
      axyd=t8(neighc)
      fid=t9(neighc)
      fid1=1.0-fid
      pxyi=t1(neighc)
      pxy=t3(neighc)
      sx4d=t6(neighc)
      dilkb=dil(kb)
      dilpkb=dilp(kb)
      dilekb=dile(kb)
      dilwkb=dilw(kb)
      dilmkb=dilm(kb)
      enlkb=enlf(kb)
      delkb=delf(kb)
      delekb=delef(kb)
      dilekb=dile(kb)
      dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-cord(iz,igrav))
      dlpkb=pxy+0.5*sx4d*dglp(kb)*(cord(kz,igrav)-cord(iz,igrav))
      dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-cord(iz,igrav))
      dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     & *(cord(kz,igrav)-cord(iz,igrav))
      axyf=(fid*dilkb+fid1*dili)
      aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
      axy=axyd*axyf
      aexy=axyd*aexyf
      dlapi=dlpi*axyf+axyd*fid1*dilpi
      dlapkb=dlpkb*axyf+axyd*fid*dilpkb
      dlaei=dlei*axyf+axyd*fid1*dilei
      dlaekb=dlekb*axyf+axyd*fid*dilekb
      dlepi=dlpi*aexyf+axyd*fid1*(dilpi*enli+dili*deli)
      dlepkb=dlpkb*aexyf+axyd*fid*(dilpkb*enlkb+dilkb*delkb)
      dleei=dlei*aexyf+axyd*fid1*(dilei*enli+dili*delei)
      dleekb=dlekb*aexyf+axyd*fid*(dilekb*enlkb+dilkb*delekb)

c these just contain derivatives wrt water fraction (fracw)


c gaz 1-12-04
c     dlekbw = -pxy*dpcefmeth(kb)
      dleiw = pxy*dpvtw
      dlekbw = -pxy*dpcp3(kb)
      dlawi=dleiw*axyf+axyd*fid1*dilwi
      dlawkb=dlekbw*axyf+axyd*fid*dilwkb
      dlewi=dleiw*aexyf+axyd*fid1*(dilwi*enli)
      dlewkb=dlekbw*aexyf+axyd*fid*(dilwkb*enlkb)

c these just contain derivatives wrt methane hydrate fraction (fracm)
c gaz 1-12-04
      dleihyd = pxy*dpvthyd
      dlekbhyd = -pxy*dpcp4(kb)
      dlami=dleihyd*axyf+axyd*fid1*dilmi
      dlamkb=dlekbhyd*axyf+axyd*fid*dilmkb
      dlemi=dleihyd*aexyf+axyd*fid1*(dilmi*enli)
      dlemkb=dlekbhyd*aexyf+axyd*fid*(dilmkb*enlkb)
c
c changed by avw -- entered here by seh
      a_axy(iau+nmatavw)=axy
      a_axy(ial+nmatavw)=-axy

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
      bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
      a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
      a(jmia+nmat(2))=a(jmia+nmat(2))+dlaei
      a(ial+nmat(1))=a(ial+nmat(1))-dlapi
      a(ial+nmat(2))=a(ial+nmat(2))-dlaei
      a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
      a(iau+nmat(2))=a(iau+nmat(2))+dlaekb
      a(jml+nmat(1))=a(jml+nmat(1))-dlapkb
      a(jml+nmat(2))=a(jml+nmat(2))-dlaekb
c
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+aexy
      bp(kz+nrhs(2))=bp(kz+nrhs(2))-aexy
      a(jmia+nmat(7))=a(jmia+nmat(7))+dlepi
      a(jmia+nmat(8))=a(jmia+nmat(8))+dleei
      a(ial+nmat(7))=a(ial+nmat(7))-dlepi
      a(ial+nmat(8))=a(ial+nmat(8))-dleei
      a(iau+nmat(7))=a(iau+nmat(7))+dlepkb
      a(iau+nmat(8))=a(iau+nmat(8))+dleekb
      a(jml+nmat(7))=a(jml+nmat(7))-dlepkb
      a(jml+nmat(8))=a(jml+nmat(8))-dleekb
c
      a(jmia+nmat(5))=a(jmia+nmat(5))+dlawi
      a(ial+nmat(5))=a(ial+nmat(5))-dlawi
      a(iau+nmat(5))=a(iau+nmat(5))+dlawkb
      a(jml+nmat(5))=a(jml+nmat(5))-dlawkb
      a(jmia+nmat(11))=a(jmia+nmat(11))+dlewi
      a(ial+nmat(11))=a(ial+nmat(11))-dlewi
      a(iau+nmat(11))=a(iau+nmat(11))+dlewkb
      a(jml+nmat(11))=a(jml+nmat(11))-dlewkb
c
      a(jmia+nmat(6))=a(jmia+nmat(6))+dlami
      a(ial+nmat(6))=a(ial+nmat(6))-dlami
      a(iau+nmat(6))=a(iau+nmat(6))+dlamkb
      a(jml+nmat(6))=a(jml+nmat(6))-dlamkb
      a(jmia+nmat(12))=a(jmia+nmat(12))+dlemi
      a(ial+nmat(12))=a(ial+nmat(12))-dlemi
      a(iau+nmat(12))=a(iau+nmat(12))+dlemkb
      a(jml+nmat(12))=a(jml+nmat(12))-dlemkb
   62 continue
      e n d i f
c
c vapour phase calculations
c
      do 63 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      pxyh=t2(neighc)
      sx4h=t7(neighc)
      vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     **(cord(kz,igrav)-cord(iz,igrav))
      t8(neighc)=vxyd
63    continue
c
c determine upwind nodes and if vapour phase exists
c
      isl=1
      do 64 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      fid=.5
c add coding to save upwind position
         if(iad.le.iad_up) then
      vxyd=t8(neighc)
      if(vxyd.lt.0.0) fid=dnwgt
      if(vxyd.gt.0.0) fid=upwgt
      if(t4(neighc).le.0.0) t9(neighc)=fid
      if(t4(neighc).gt.0.0) t9(neighc)=fid
c
            call setbit(nbits,neighc,upwind_v(iz4m1),fid)
c
         else
           if(bit(nbits,neighc,upwind_v(iz4m1))) then
            t9(neighc)=1.0
           else
            t9(neighc)=0.0
           endif
         endif
64    continue
c
c form equations
c
      if(isl.ne.0) then
      do 65 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      iau=it11(jm)
      ial=it12(jm)
      jml=nelmdg(kz)-neqp1
      fid=t9(neighc)
      fid1=1.0-fid
      pxyh=t2(neighc)
      pvxy=t4(neighc)
      sx4h=t7(neighc)
      vxyd=t8(neighc)
      divkb=div(kb)
      divpkb=divp(kb)
      divekb=dive(kb)
      divwkb=divw(kb)
      divmkb=divm(kb)
      envkb=envf(kb)
      devkb=devf(kb)
      devekb=devef(kb)
      divekb=dive(kb)
      dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-cord(iz,igrav))
      dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-cord(iz,igrav))
      dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
      dvekb=0.5*sx4h*dgve(kb)
     & *(cord(kz,igrav)-cord(iz,igrav))
      vxyf=(fid*divkb+fid1*divi)
      vexyf=(fid*divkb*envkb+fid1*divi*envi)
      vxy=vxyd*vxyf
      vexy=vxyd*vexyf
      dvapi=dvpi*vxyf+vxyd*fid1*divpi
      dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
      dvaei=dvei*vxyf+vxyd*fid1*divei
      dvaekb=dvekb*vxyf+vxyd*fid*divekb
      dvepi=dvpi*vexyf+vxyd*fid1*(divpi*envi+divi*devi)
      dvepkb=dvpkb*vexyf+vxyd*fid*(divpkb*envkb+divkb*devkb)
      dveei=dvei*vexyf+vxyd*fid1*(divei*envi+divi*devei)
      dveekb=dvekb*vexyf+vxyd*fid*(divekb*envkb+divkb*devekb)
c changed by avw -- entered here by seh
      a_vxy(iau+nmatavw)=vxy
      a_vxy(ial+nmatavw)=-vxy

c these just contain derivatives wrt water fraction

      dvawi=vxyd*fid1*divwi
      dvawkb=vxyd*fid*divwkb
      dvewi=vxyd*fid1*(divwi*envi)
      dvewkb=vxyd*fid*(divwkb*envkb)

c these just contain derivatives wrt methane hydrate fraction

      dvami=vxyd*fid1*divmi
      dvamkb=vxyd*fid*divmkb
      dvemi=vxyd*fid1*(divmi*envi)
      dvemkb=vxyd*fid*(divmkb*envkb)

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+vxy
      bp(kz+nrhs(1))=bp(kz+nrhs(1))-vxy
      a(jmia+nmat(1))=a(jmia+nmat(1))+dvapi
      a(jmia+nmat(2))=a(jmia+nmat(2))+dvaei
      a(ial+nmat(1))=a(ial+nmat(1))-dvapi
      a(ial+nmat(2))=a(ial+nmat(2))-dvaei
      a(iau+nmat(1))=a(iau+nmat(1))+dvapkb
      a(iau+nmat(2))=a(iau+nmat(2))+dvaekb
      a(jml+nmat(1))=a(jml+nmat(1))-dvapkb
      a(jml+nmat(2))=a(jml+nmat(2))-dvaekb
c
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+vexy
      bp(kz+nrhs(2))=bp(kz+nrhs(2))-vexy
      a(jmia+nmat(7))=a(jmia+nmat(7))+dvepi
      a(jmia+nmat(8))=a(jmia+nmat(8))+dveei
      a(ial+nmat(7))=a(ial+nmat(7))-dvepi
      a(ial+nmat(8))=a(ial+nmat(8))-dveei
      a(iau+nmat(7))=a(iau+nmat(7))+dvepkb
      a(iau+nmat(8))=a(iau+nmat(8))+dveekb
      a(jml+nmat(7))=a(jml+nmat(7))-dvepkb
      a(jml+nmat(8))=a(jml+nmat(8))-dveekb
c
      a(jmia+nmat(5))=a(jmia+nmat(5))+dvawi
      a(ial+nmat(5))=a(ial+nmat(5))-dvawi
      a(iau+nmat(5))=a(iau+nmat(5))+dvawkb
      a(jml+nmat(5))=a(jml+nmat(5))-dvawkb
      a(jmia+nmat(11))=a(jmia+nmat(11))+dvewi
      a(ial+nmat(11))=a(ial+nmat(11))-dvewi
      a(iau+nmat(11))=a(iau+nmat(11))+dvewkb
      a(jml+nmat(11))=a(jml+nmat(11))-dvewkb
c
      a(jmia+nmat(6))=a(jmia+nmat(6))+dvami
      a(ial+nmat(6))=a(ial+nmat(6))-dvami
      a(iau+nmat(6))=a(iau+nmat(6))+dvamkb
      a(jml+nmat(6))=a(jml+nmat(6))-dvamkb
      a(jmia+nmat(12))=a(jmia+nmat(12))+dvemi
      a(ial+nmat(12))=a(ial+nmat(12))-dvemi
      a(iau+nmat(12))=a(iau+nmat(12))+dvemkb
      a(jml+nmat(12))=a(jml+nmat(12))-dvemkb
   65 continue
      e n d i f
c
c add heat conduction
c
      do 66 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      heatc=t5(neighc)
      iau=it11(jm)
      ial=it12(jm)
      jml=nelmdg(kz)-neqp1
      heatc=t5(neighc)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatc*(t(kb)-ti)
      bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatc*(t(kb)-ti)
      a(jmia+nmat(7))=a(jmia+nmat(7))-heatc*dtpa(i)
      a(jmia+nmat(8))=a(jmia+nmat(8))-heatc*dtpae(i)
      a(ial+nmat(7))=a(ial+nmat(7))+heatc*dtpa(i)
      a(ial+nmat(8))=a(ial+nmat(8))+heatc*dtpae(i)
      a(iau+nmat(7))=a(iau+nmat(7))+heatc*dtpa(kb)
      a(iau+nmat(8))=a(iau+nmat(8))+heatc*dtpae(kb)
      a(jml+nmat(7))=a(jml+nmat(7))-heatc*dtpa(kb)
      a(jml+nmat(8))=a(jml+nmat(8))-heatc*dtpae(kb)
66    continue

c add accumulation terms

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i) + skwhyd(i) 
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i) + qhwhyd(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i) + dskhyd1(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i) + dskhyd2(i)
c      a(jmia+nmat(7))=a(jmia+nmat(7))+sx1d*depf(i)+dqh(i) + dqhhyd1(i)
c      a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*deef(i)+deqh(i) + dqhhyd2(i)
      a(jmia+nmat(7))=a(jmia+nmat(7))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*deef(i)+deqh(i)

c     derivative wrt hydrate fraction

      a(jmia+nmat(6))=a(jmia+nmat(6))+sx1d*dmmf(i) + dskhyd4(i) + dq4(i)
c      a(jmia+nmat(12))=a(jmia+nmat(12))+sx1d*demf(i) + dqhhyd4(i)
      a(jmia+nmat(12))=a(jmia+nmat(12))+sx1d*demf(i)

c     derivative wrt water fraction

      a(jmia+nmat(5))=a(jmia+nmat(5))+sx1d*dmwf(i) + dq3(i) + dskhyd3(i)
      a(jmia+nmat(11))=a(jmia+nmat(11))+sx1d*dewf(i) + deqw(i) 

      r e t u r n
      e    n    d
