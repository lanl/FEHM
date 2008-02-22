      subroutine  geneq_meth ( i )
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
!D1  PURPOSE
!D1
!D1  This subroutine generates the equations for 3-dimensional heat
!D1  and mass transfer for methane componet of a mixture 
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 10-OCT-02    G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/geneq_meth.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.3.2 Heat- and mass-transfer equations
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

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

      real*8  divwi, divwkb, dilwi, dilwkb                      
      real*8  dlawi, dlawkb, dlewi, dlewkb                      
      real*8  dvawi, dvawkb, dvewi, dvewkb                      
      real*8  dleiw, dlekbw, dpvtw                              

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
      pvii=phimeth(i)
      phii=pvii-pcpmeth(i)
      dpvti=dpcpmeth2(i)
      enli=enlfmeth(i)
      deli=delf(i)
      delei=delef(i)
      envi=envfmeth(i)
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
      ti=tmeth(i)
      swi=smeth(i)
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
      a_axy(jmia+nmatavw)=skmeth(i)*smeth(i)
c Take care of souce/sink term
c If this is an isothermal air-water simulation
      a_vxy(jmia+nmatavw)=skmeth(i)*(1-smeth(i))

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
            pvikb=phimeth(kb)
            phikb=pvikb-pcpmeth(kb)
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
            pvikb=phimeth(kb)
            phikb=pvikb-pcpmeth(kb)
c           pxy=sx2c*perml(1)+sx3c*perml(2)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2))
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
      axyd=pxyi+0.5*sx4d*(rolfmeth(i)+rolfmeth(kb))
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
      enlkb=enlfmeth(kb)
      delkb=delf(kb)
      delekb=delef(kb)
      dilekb=dile(kb)
      dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-cord(iz,igrav))
      dlpkb=pxy+0.5*sx4d*dglp(kb)*(cord(kz,igrav)-cord(iz,igrav))
      dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-cord(iz,igrav))
      dlekb=-pxy*dpcpmeth2(kb)+0.5*sx4d*dgle(kb)
     **(cord(kz,igrav)-cord(iz,igrav))
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

c these just contain derivatives wrt methane hydrate fraction
c note: for methane equations water faction dependence is only
c in accumulation term

      dlami=axyd*fid1*dilmi
      dlamkb=axyd*fid*dilmkb
      dlemi=axyd*fid1*(dilmi*enli)
      dlemkb=axyd*fid*(dilmkb*enlkb)
c Modified by RJP 12/01/2003. Added dependence on water fraction
      dlawi=axyd*fid1*dilwi
      dlawkb=axyd*fid*dilwkb
      dlewi=axyd*fid1*(dilwi*enli)
      dlewkb=axyd*fid*(dilwkb*enlkb)
c
c changed by avw -- entered here by seh
      a_axy(iau+nmatavw)=axy
      a_axy(ial+nmatavw)=-axy

      bp(iz+nrhs(3))=bp(iz+nrhs(3))+axy
      bp(kz+nrhs(3))=bp(kz+nrhs(3))-axy
      a(jmia+nmat(15))=a(jmia+nmat(15))+dlapi
      a(jmia+nmat(16))=a(jmia+nmat(16))+dlaei
      a(ial+nmat(15))=a(ial+nmat(15))-dlapi
      a(ial+nmat(16))=a(ial+nmat(16))-dlaei
      a(iau+nmat(15))=a(iau+nmat(15))+dlapkb
      a(iau+nmat(16))=a(iau+nmat(16))+dlaekb
      a(jml+nmat(15))=a(jml+nmat(15))-dlapkb
      a(jml+nmat(16))=a(jml+nmat(16))-dlaekb
c
      bp(iz+nrhs(4))=bp(iz+nrhs(4))+aexy
      bp(kz+nrhs(4))=bp(kz+nrhs(4))-aexy
      a(jmia+nmat(21))=a(jmia+nmat(21))+dlepi
      a(jmia+nmat(22))=a(jmia+nmat(22))+dleei
      a(ial+nmat(21))=a(ial+nmat(21))-dlepi
      a(ial+nmat(22))=a(ial+nmat(22))-dleei
      a(iau+nmat(21))=a(iau+nmat(21))+dlepkb
      a(iau+nmat(22))=a(iau+nmat(22))+dleekb
      a(jml+nmat(21))=a(jml+nmat(21))-dlepkb
      a(jml+nmat(22))=a(jml+nmat(22))-dleekb
c
      a(jmia+nmat(18))=a(jmia+nmat(18))+dlami
      a(ial+nmat(18))=a(ial+nmat(18))-dlami
      a(iau+nmat(18))=a(iau+nmat(18))+dlamkb
      a(jml+nmat(18))=a(jml+nmat(18))-dlamkb
      a(jmia+nmat(24))=a(jmia+nmat(24))+dlemi
      a(ial+nmat(24))=a(ial+nmat(24))-dlemi
      a(iau+nmat(24))=a(iau+nmat(24))+dlemkb
      a(jml+nmat(24))=a(jml+nmat(24))-dlemkb
c Modified by RJP 12/01/2003. Added dependence on water fraction
c
      a(jmia+nmat(17))=a(jmia+nmat(17))+dlawi
      a(ial+nmat(17))=a(ial+nmat(17))-dlawi
      a(iau+nmat(17))=a(iau+nmat(17))+dlawkb
      a(jml+nmat(17))=a(jml+nmat(17))-dlawkb
      a(jmia+nmat(23))=a(jmia+nmat(23))+dlewi
      a(ial+nmat(23))=a(ial+nmat(23))-dlewi
      a(iau+nmat(23))=a(iau+nmat(23))+dlewkb
      a(jml+nmat(23))=a(jml+nmat(23))-dlewkb
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
      vxyd=pxyh+0.5*sx4h*(rovfmeth(i)+rovfmeth(kb))
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
      envkb=envfmeth(kb)
      devkb=devf(kb)
      devekb=devef(kb)
      divekb=dive(kb)
      dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-cord(iz,igrav))
      dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-cord(iz,igrav))
      dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
      dvekb=0.5*sx4h*dgve(kb)
     **(cord(kz,igrav)-cord(iz,igrav))
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

c cross derivatives for methane
c these just contain derivatives wrt methane fraction

      dvami=vxyd*fid1*divmi
      dvamkb=vxyd*fid*divmkb
      dvemi=vxyd*fid1*(divmi*envi)
      dvemkb=vxyd*fid*(divmkb*envkb)
c Modified by RJP 12/01/2003. Added dependence on water fraction
c     

      dvawi=vxyd*fid1*divwi
      dvawkb=vxyd*fid*divwkb
      dvewi=vxyd*fid1*(divwi*envi)
      dvewkb=vxyd*fid*(divwkb*envkb)
      bp(iz+nrhs(3))=bp(iz+nrhs(3))+vxy
      bp(kz+nrhs(3))=bp(kz+nrhs(3))-vxy
      a(jmia+nmat(15))=a(jmia+nmat(15))+dvapi
      a(jmia+nmat(16))=a(jmia+nmat(16))+dvaei
      a(ial+nmat(15))=a(ial+nmat(15))-dvapi
      a(ial+nmat(16))=a(ial+nmat(16))-dvaei
      a(iau+nmat(15))=a(iau+nmat(15))+dvapkb
      a(iau+nmat(16))=a(iau+nmat(16))+dvaekb
      a(jml+nmat(15))=a(jml+nmat(15))-dvapkb
      a(jml+nmat(16))=a(jml+nmat(16))-dvaekb
c
      bp(iz+nrhs(4))=bp(iz+nrhs(4))+vexy
      bp(kz+nrhs(4))=bp(kz+nrhs(4))-vexy
      a(jmia+nmat(21))=a(jmia+nmat(21))+dvepi
      a(jmia+nmat(22))=a(jmia+nmat(22))+dveei
      a(ial+nmat(21))=a(ial+nmat(21))-dvepi
      a(ial+nmat(22))=a(ial+nmat(22))-dveei
      a(iau+nmat(21))=a(iau+nmat(21))+dvepkb
      a(iau+nmat(22))=a(iau+nmat(22))+dveekb
      a(jml+nmat(21))=a(jml+nmat(21))-dvepkb
      a(jml+nmat(22))=a(jml+nmat(22))-dveekb
c
      a(jmia+nmat(18))=a(jmia+nmat(18))+dvami
      a(ial+nmat(18))=a(ial+nmat(18))-dvami
      a(iau+nmat(18))=a(iau+nmat(18))+dvamkb
      a(jml+nmat(18))=a(jml+nmat(18))-dvamkb
      a(jmia+nmat(24))=a(jmia+nmat(24))+dvemi
      a(ial+nmat(24))=a(ial+nmat(24))-dvemi
      a(iau+nmat(24))=a(iau+nmat(24))+dvemkb
      a(jml+nmat(24))=a(jml+nmat(24))-dvemkb
c Modified by RJP 12/01/2003. Added water fraction dependence
c
      a(jmia+nmat(17))=a(jmia+nmat(17))+dvawi
      a(ial+nmat(17))=a(ial+nmat(17))-dvawi
      a(iau+nmat(17))=a(iau+nmat(17))+dvawkb
      a(jml+nmat(17))=a(jml+nmat(17))-dvawkb
      a(jmia+nmat(23))=a(jmia+nmat(23))+dvewi
      a(ial+nmat(23))=a(ial+nmat(23))-dvewi
      a(iau+nmat(23))=a(iau+nmat(23))+dvewkb
      a(jml+nmat(23))=a(jml+nmat(23))-dvewkb

   65 continue
      e n d i f
c
c add heat conduction
c
c  heat conduction is only in water equation
c

c add accumulation terms

      bp(iz+nrhs(3))=bp(iz+nrhs(3))+sx1d*denmethi(i)+skmeth(i)
     &              + skmhyd(i)
      bp(iz+nrhs(4))=bp(iz+nrhs(4))+sx1d*denemethi(i)+qhmeth(i)
     &              + qhmhyd(i)
c deleted dq(i) on next line
      a(jmia+nmat(15))=a(jmia+nmat(15))+sx1d*dmpf(i)
     &              + dskmethw1(i) + dskhyd1(i)
      a(jmia+nmat(16))=a(jmia+nmat(16))+sx1d*dmef(i)+dqt(i)
     &              + dskhyd2(i)
      a(jmia+nmat(21))=a(jmia+nmat(21))+sx1d*depf(i)+dqmh(i)
      a(jmia+nmat(22))=a(jmia+nmat(22))+sx1d*deef(i)+deqmh(i)

      a(jmia+nmat(17))=a(jmia+nmat(17))+sx1d*dmwf(i) 
     &              + dskmethw3(i) + dskhyd3(i)
      a(jmia+nmat(18))=a(jmia+nmat(18))+sx1d*dmmf(i)
     &              + dskmethw4(i) + dskhyd4(i)
      a(jmia+nmat(23))=a(jmia+nmat(23))+sx1d*dewf(i)+deqw(i)
      a(jmia+nmat(24))=a(jmia+nmat(24))+sx1d*demf(i)+deqm(i)


      r e t u r n
      e    n    d
