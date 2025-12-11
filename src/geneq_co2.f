      subroutine  geneq_co2 ( i )
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

C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine generates the equations for 3-dimensional heat
CD1  and mass transfer for co2 componet of a mixture 
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-OCT-02    G. Zyvoloski           Initial implementation.
CD2
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
CD4 RJP 12/02/04 Major revision for CO2 only application. Hydrate
CD4 formation is ignored.
CD4 RJP 12/03/04 Reduced size of a, similar to geneq_h2o subroutine
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
      use comco2
      use comriv
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
      integer imd,iwd, icesd
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8  aexy, aexyf, alxi, alxkb, alyi, alykb, alzi, alzkb
      real*8  avxi, avyi, avzi, axi, axkb, axy, axyd, axyf
      real*8  ayi, aykb, azi, azkb
      real*8  delekb, delkb, devei, devekb, devi, devkb 
      real*8  dilei, dilekb, dili, dilkb, dilpi, dilpkb 
      real*8  divei, divekb, divi, divkb, divpi, divpkb
      real*8  dlaei, dlaekb, dlapi, dlapkb, dleei, dleekb 
      real*8  dlei, dlekb, dlepi, dlepkb, dlpi, dlpkb, dpcpwi  
      real*8  dvaei, dvaekb, dvapi, dvapkb, dveei, dveekb
      real*8  dvei, dvekb, dvepi, dvepkb, dvpi, dvpkb
      real*8  fid, fid1, grav_air, heatc, enlkb, envi, envkb
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyi, pxyh
      real*8  radi, radkb, plii
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
      
      real*8  dilyci,dilyai,divyci,divyai,dilh2omi   
      real*8  dilh2okb,dilh2opkb,dilh2oekb,dilh2owkb,dilh2omkb

      real*8 rowi, dgwpi, dgwei, dgwwi, dgwyci, dgwyai
      real*8 roci, dgcpi, dgcei, dgcwi, dgcyci, dgcyai
      real*8 roli, dglpi, dglei, dglwi, dglyci, dglyai
      real*8 enwi, dewi, dewei, dewwi, dewyci, dewyai
      real*8 enci, deci, decei, decwi, decyci, decyai
      real*8 enli, deli, delei, delwi, delyci, delyai

      real*8 diwi, diwpi, diwei, diwwi, diwyci, diwyai

      real*8 rowkb, dgwpkb, dgwekb, dgwwkb, dgwyckb, dgwyakb
      real*8 rockb, dgcpkb, dgcekb, dgcwkb, dgcyckb, dgcyakb
      real*8 rolkb, dglpkb, dglekb, dglwkb, dglyckb, dglyakb
      real*8 enwkb, dewkb, dewekb, dewwkb, dewyckb, dewyakb
      real*8 enckb, deckb, decekb, decwkb, decyckb, decyakb
      real*8 delyckb, delyakb

      real*8 diwkb, diwpkb, diwekb, diwwkb, diwyckb, diwyakb
      real*8 plikb, pxyl, dwpi, dwpkb, dwei, dwekb
      real*8 dwapi, dwapkb, dwaei, dwaekb, dwepi, dwepkb
      real*8 dweei, dweekb, dwwi, dwwkb, dwawi
      real*8 dwawkb, dwewkb, dwyci, dwyckb, dwayci,dwayckb
      real*8 dweyci, dweyckb, dwyai, dwyakb, dwayai, dwayakb
      real*8 dweyai, dweyakb, dwewi, divyckb, divyakb
      real*8 dvyci, dvyckb, dvayci, dvayckb, dveyci, dveyckb
      real*8 dvyai, dvyakb,dvayakb,dveyai,dveyakb
      real*8 dilyckb, dilyakb, dvayai,delwkb,lxyf,lxyd,lexyf,lxy,lexy
      
      real*8 dlyci, dlyckb, dlayci, dlayckb, dleyckb
      real*8 dlyai, dlyakb, dlayai, dlayakb, dleyakb
      real*8 dleyci, dleyai, yci, yckb, diffc
      real*8 dfxd, hfc, dhfcp, dhfct, dhfcw, dhfcyc, dhfcya

      real*8 diffi,diffkb,df2t,df3t,dfzt,df3c

      parameter(dis_tol=1.d-12)
c     
c     
                                ! Guessed at value for grav_air
      grav_air = 0.
      
      hfc = 0.0
      dhfcp = 0.0
      dhfct = 0.0
      
c     changed by avw 4/95 -- entered into new version by seh
      neqp1=neq+1
      ldna=nelm(neqp1)-neqp1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif
      
      icesd=ices(i)
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
c     phico2 is pressure in gas phase for 2-phase co2 or pressure in 
c     co2-rich phase for single phase co2
      pvii=phico2(i)
c     pcp is cap. pr. between CO2-liq/gas for single phase CO2 and water.
c     pcp is cap. pr. between CO2-gas phase and water for 2-phase CO2
c     pcpl is cap. pr. between CO2-liq phase and CO2-gas phase for 2-phase CO2
      phii=phi(i)
      plii=pvii-pcg(i)
      dpcpwi=dpcpw(i)

      rowi=wat_prop(i)
      dgwpi=wat_prop(neq+i)
      dgwei=wat_prop(2*neq+i)
      dgwwi=0.0
      dgwyci=wat_prop(3*neq+i)
      dgwyai=wat_prop(4*neq+i)	
      enwi=wat_prop(5*neq+i)
      dewi=wat_prop(6*neq+i)
      dewei=wat_prop(7*neq+i)
      dewwi=0.0
      dewyci=wat_prop(11*neq+i)
      dewyai=0.0


      roci=co2_prop(9*neq+i)
      dgcei=co2_prop(10*neq+i)
      dgcpi=co2_prop(11*neq+i)
      dgcwi=0.0
      dgcyci=0.0
      dgcyai=0.0
      enci=co2_prop(12*neq+i)
      deci=co2_prop(14*neq+i)
      decei=co2_prop(13*neq+i)
      decwi=0.0
      decyci=0.0
      decyai=0.0

      roli=co2_prop(i)
      dglei=co2_prop(neq+i)
      dglpi=co2_prop(2*neq+i)
      dglwi=0.0
      dglyci=0.0
      dglyai=0.0
      enli=co2_prop(3*neq+i)
      deli=co2_prop(5*neq+i)
      delei=co2_prop(4*neq+i)
      delwi=0.0
      delyci=0.0
      delyai=0.0

      diwi=diw(i)
      diwpi=diwp(i)
      diwei=diwe(i)
      diwwi=diww(i)
      diwyci=diwyc(i)
      diwyai=diwya(i)
      dili=dil(i)
      divi=div(i)
      dilpi=dilp(i)
      dilei=dile(i)
      dilwi=dilw(i)
      dilyci=dilyc(i)
      dilyai=dilya(i)
      divpi=divp(i)
      divei=dive(i)
      divwi=divw(i)
      divyci=divyc(i)
      divyai=divya(i)

      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      ti=tco2(i)

      yci = yc(i)
      diffi = diff(i)*ps(i)
c     
c     form constants for i>neq
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

c     Take care of source/sink term by putting it into empty slot (from
c     same node to same node).
c     If this is an isothermal air-water simulation
c     RJP 6/29/04 Took out a_axy and a_vxy assignment for now
      c_axy(jmia+nmatavw)=skco2(i)
c     Take care of souce/sink term
c     If this is an isothermal air-water simulation
c     a_vxy(jmia+nmatavw)=skco2(i)*(1-sco2(i))

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
 68      continue
 58   continue
      if(icnl.eq.0) then
c     
c     3-d geometry
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
c     RJP 08/02/07 modified below. Took out reduction_factor from earlier version
c     instead multiplied with it later
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
c     
c     RJP 05/18/08 diffusivity related terms
c     
            diffkb=diff(kb)*ps(kb)
            if (diffi .eq. 0. .or. diffkb .eq. 0) then
               df2t = 0.d0
               df3t = 0.d0
               dfzt = 0.d0
            else
               df2t=2.*diffi*diffkb/(diffi+diffkb)
               df3t=2.*diffi*diffkb/(diffi+diffkb)
               dfzt=2.*diffi*diffkb/(diffi+diffkb)
            end if
            pvikb=phico2(kb)
            phikb=phi(kb)
            plikb=pvikb-pcg(kb)
c     pxy=sx2c*perml(1)+sx3c*perml(2)+sxzc*perml(3)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            if(iriver.ne.0) then
               if(delz2.eq.0.d0.and.kb.lt.(neq_primary-npoint_riv)) then
                  if((mdnodes_riv(i).ne.0).or.
     &                 (mdnodes_riv(kb).ne.0)) then
                     delx2=mod_dis(iw-nic_old,1)**2
                     dely2=mod_dis(iw-nic_old,2)**2
                  endif
               endif
            endif
            dis2=delx2+dely2+delz2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
            endif
c     RJP 08/02/07 added below
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            pxyl=pxy*(plikb-plii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/(delx2/sx2t+dely2/sx3t+delz2/sxzt)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
c     
c     RJP 05/18/08 diffusivity related terms
c     
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               if(diffi.eq.0.d0.or.diffkb.eq.0) then
                  df3c = 0.d0
               else
                  df3c=sx2c*dis2/(delx2/df2t+dely2/df3t+delz2/dfzt)
               endif
            else
               df3c=sx2c*sx_mult*max(df2t,df3t,dfzt)
            endif

            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c                    
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
            t10(neighc)=pxyl
            t13(neighc)=-grav*t4(neighc)
c     
c     RJP 05/18/08 added new array for diffusivity related terms storage
c     
            t17(neighc)=df3c

 59      continue
      else if ( icnl.ne.0 )  t h e n
c     
c     2-d geometry
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
c     RJP 08/02/07 modified below. Took out reduction_factor from earlier version
c     instead multiplied with it later
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
            thxkb=thx(kb)
            thykb=thy(kb)
c     
c     RJP 05/18/08 added diffusivity terms
c     
            diffkb=diff(kb)*ps(kb)
            if (diffi .eq. 0. .or. diffkb .eq. 0) then
                df2t = 0.d0
                df3t = 0.d0
            else
                df2t=2.*diffi*diffkb/(diffi+diffkb)
                df3t=2.*diffi*diffkb/(diffi+diffkb)
            end if
c     
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            pvikb=phico2(kb)
            phikb=phi(kb)
            plikb=pvikb-pcg(kb)
c     pxy=sx2c*perml(1)+sx3c*perml(2)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+dely2/perml(2))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
c     RJP 08/02/07 added below
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            pxyl=pxy*(plikb-plii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/(delx2/sx2t+dely2/sx3t)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
c     
c     RJP 05/18/08 added diffusivity terms
c     
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               if(diffi.eq.0.d0.or.diffkb.eq.0) then
                  df3c = 0.d0
               else
                  df3c=sx2c*dis2/(delx2/df2t+dely2/df3t)
               endif
            else
               sx3c=sx2c*sx_mult*max(df2t,df3t)
            endif
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
            t10(neighc)=pxyl
            t13(neighc)=-grav*t4(neighc)
c     
c     RJP 05/18/08 added new array for diffusivity related terms storage
c     
            t17(neighc)=df3c

 69      continue
c     
      endif
c     
c     water-rich phase calculations
c     
      do 60 jm=1,iq
         kb=it8(jm)
         rowkb=wat_prop(kb)
         kz=kb-icd
         neighc=it9(jm)
         pxyi=t1(neighc)
         sx4d=t6(neighc)
         axyd=pxyi+0.5*sx4d*(rowi+rowkb)
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=axyd
 60   continue
c     
c     determine upwind nodes and if liquid phase exists
c     
      isl=1
      do 61 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         fid=.5
c     add coding to save upwind position
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
 61   continue
c     
c     form equations
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
            diwkb=diw(kb)
            diwpkb=diwp(kb)
            diwekb=diwe(kb)
            diwwkb=diww(kb)
            diwyckb=diwyc(kb)
            diwyakb=diwya(kb)
            dgwyckb=wat_prop(3*neq+kb)
            dgwyakb=wat_prop(4*neq+kb)
            enwkb=wat_prop(5*neq+kb)
            dewkb=wat_prop(6*neq+kb)
            dewekb=wat_prop(7*neq+kb)
c            enckb=co2_prop(12*neq+kb)
c            deckb=co2_prop(14*neq+kb)
c            decekb=co2_prop(13*neq+kb)
c            decwkb= 0.d0
            dewwkb=0.d0
            dewyckb=wat_prop(11*neq+kb)
            dewyakb=0.d0
c            decyckb = 0.0
            dwpi=-pxy+0.5*sx4d*wat_prop(neq+i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dwpkb=pxy+0.5*sx4d*wat_prop(neq+kb)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dwei=0.5*sx4d*wat_prop(2*neq+i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dwekb=0.5*sx4d*wat_prop(2*neq+kb)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            if(icesd.eq.2) dwei=0.d0
            if(icesd.eq.2) dwekb=0.d0
            axyf=(fid*diwkb+fid1*diwi)
            aexyf=(fid*diwkb*enwkb+fid1*diwi*enwi)
            axy=axyd*axyf
            aexy=axyd*aexyf
            dwapi=dwpi*axyf+axyd*fid1*diwpi
            dwapkb=dwpkb*axyf+axyd*fid*diwpkb
            dwaei=dwei*axyf+axyd*fid1*diwei
            dwaekb=dwekb*axyf+axyd*fid*diwekb
            dwepi=dwpi*aexyf+axyd*fid1*(diwpi*enwi+diwi*dewi)
            dwepkb=dwpkb*aexyf+axyd*fid*(diwpkb*enwkb+diwkb*dewkb)
            dweei=dwei*aexyf+axyd*fid1*(diwei*enwi+diwi*dewei)
            dweekb=dwekb*aexyf+axyd*fid*(diwekb*enwkb+diwkb*dewekb)

c     derivatives wrt water-rich phase fraction

            dwwi = pxy*dpcpwi
            dwwkb = -pxy*dpcpw(kb)
            dwawi=dwwi*axyf+axyd*fid1*diwwi
            dwawkb=dwwkb*axyf+axyd*fid*diwwkb
            dwewi=dwwi*aexyf+axyd*fid1*(diwwi*enwi+diwi*dewwi)
            dwewkb=dwwkb*aexyf+axyd*fid*(diwwkb*enwkb+diwkb*dewwkb)

c     derivatives wrt co2 mass fraction in water-rich phase

            dwyci = 0.5*sx4d*dgwyci*(cord(kz,igrav)-cord(iz,igrav))
            dwyckb = 0.5*sx4d*dgwyckb*(cord(kz,igrav)-cord(iz,igrav))
            dwayci = dwyci*axyf+axyd*fid1*diwyci
            dwayckb = dwyckb*axyf+axyd*fid*diwyckb
            dweyci = dwyci*aexyf+axyd*fid1*(diwyci*enwi+diwi*dewyci)
            dweyckb = dwyckb*aexyf+axyd*fid*(diwyckb*enwkb+diwkb*
     &           dewyckb)

            if(iprtype.eq.4) then
               if(ico2dis(i).eq.0) then
                  dwawi=dwayci
                  dwewi=dweyci
               else
                  if(ices(i).eq.2) then
                     dwapi=dwapi+dwayci*dmol(i)
     &                    +dwayci*dmol(i+neq)
                     dwepi=dwepi+dweyci*dmol(i)
     &                    +dweyci*dmol(i+neq)
                  else
                     dwapi=dwapi+dwayci*dmol(i)
                     dwaei=dwaei+dwayci*dmol(i+neq)
                     dwepi=dwepi+dweyci*dmol(i)
                     dweei=dweei+dweyci*dmol(i+neq)
                  endif
               endif
               if(ico2dis(kb).eq.0) then
                  dwawkb=dwayckb
                  dwewkb=dweyckb
               else
                  if(ices(kb).eq.2) then
                     dwapkb=dwapkb+dwayckb*dmol(kb)
     &                    +dwayckb*dmol(kb+neq)
                     dwepkb=dwepkb+dweyckb*dmol(kb)
     &                    +dweyckb*dmol(kb+neq)
                  else
                     dwapkb=dwapkb+dwayckb*dmol(kb)
                     dwaekb=dwaekb+dwayckb*dmol(kb+neq)
                     dwepkb=dwepkb+dweyckb*dmol(kb)
                     dweekb=dweekb+dweyckb*dmol(kb+neq)
                  endif
               endif
            endif
            
c     derivatives wrt air mass fraction in water-rich phase

            dwyai = 0.5*sx4d*dgwyai*(cord(kz,igrav)-cord(iz,igrav))
            dwyakb = 0.5*sx4d*dgwyakb*(cord(kz,igrav)-cord(iz,igrav))
            dwayai = dwyai*axyf+axyd*fid1*diwyai
            dwayakb = dwyakb*axyf+axyd*fid*diwyakb
            dweyai = dwyai*aexyf+axyd*fid1*(diwyai*enwi+diwi*dewyai)
            dweyakb = dwyakb*aexyf+axyd*fid*(diwyakb*enwkb+diwkb*
     &           dewyakb)

c     
c     RJP 07/05/07 added new co2 flux arrays c_axy and c_vxy
c     
            c_axy(iau+nmatavw)=axy
            c_axy(ial+nmatavw)=-axy

            bp(iz+nrhs(3))=bp(iz+nrhs(3))+axy
            bp(kz+nrhs(3))=bp(kz+nrhs(3))-axy
            a(jmia+nmat(11))=a(jmia+nmat(11))+dwapi
            a(jmia+nmat(12))=a(jmia+nmat(12))+dwaei
            a(ial+nmat(11))=a(ial+nmat(11))-dwapi
            a(ial+nmat(12))=a(ial+nmat(12))-dwaei
            a(iau+nmat(11))=a(iau+nmat(11))+dwapkb
            a(iau+nmat(12))=a(iau+nmat(12))+dwaekb
            a(jml+nmat(11))=a(jml+nmat(11))-dwapkb
            a(jml+nmat(12))=a(jml+nmat(12))-dwaekb

c     
            bp(iz+nrhs(4))=bp(iz+nrhs(4))+aexy
            bp(kz+nrhs(4))=bp(kz+nrhs(4))-aexy

            a(jmia+nmat(16))=a(jmia+nmat(16))+dwepi
            a(jmia+nmat(17))=a(jmia+nmat(17))+dweei
            a(ial+nmat(16))=a(ial+nmat(16))-dwepi
            a(ial+nmat(17))=a(ial+nmat(17))-dweei
            a(iau+nmat(16))=a(iau+nmat(16))+dwepkb
            a(iau+nmat(17))=a(iau+nmat(17))+dweekb
            a(jml+nmat(16))=a(jml+nmat(16))-dwepkb
            a(jml+nmat(17))=a(jml+nmat(17))-dweekb
c     
            a(jmia+nmat(13))=a(jmia+nmat(13))+dwawi
            a(ial+nmat(13))=a(ial+nmat(13))-dwawi
            a(iau+nmat(13))=a(iau+nmat(13))+dwawkb
            a(jml+nmat(13))=a(jml+nmat(13))-dwawkb
            a(jmia+nmat(18))=a(jmia+nmat(18))+dwewi
            a(ial+nmat(18))=a(ial+nmat(18))-dwewi
            a(iau+nmat(18))=a(iau+nmat(18))+dwewkb
            a(jml+nmat(18))=a(jml+nmat(18))-dwewkb
c     
            a(jmia+nmat(14))=a(jmia+nmat(14))+dwayci
            a(ial+nmat(14))=a(ial+nmat(14))-dwayci
            a(iau+nmat(14))=a(iau+nmat(14))+dwayckb
            a(jml+nmat(14))=a(jml+nmat(14))-dwayckb
            a(jmia+nmat(19))=a(jmia+nmat(19))+dweyci
            a(ial+nmat(19))=a(ial+nmat(19))-dweyci
            a(iau+nmat(19))=a(iau+nmat(19))+dweyckb
            a(jml+nmat(19))=a(jml+nmat(19))-dweyckb
c     
            a(jmia+nmat(15))=a(jmia+nmat(15))+dwayai
            a(ial+nmat(15))=a(ial+nmat(15))-dwayai
            a(iau+nmat(15))=a(iau+nmat(15))+dwayakb
            a(jml+nmat(15))=a(jml+nmat(15))-dwayakb
            a(jmia+nmat(20))=a(jmia+nmat(20))+dweyai
            a(ial+nmat(20))=a(ial+nmat(20))-dweyai
            a(iau+nmat(20))=a(iau+nmat(20))+dweyakb
            a(jml+nmat(20))=a(jml+nmat(20))-dweyakb

 62      continue
      endif
c     
c     co2-rich phase calculations
c     
      do 63 jm=1,iq
         kb=it8(jm)
         rockb = co2_prop(9*neq+kb)
         kz=kb-icd
         neighc=it9(jm)
         pxyh=t2(neighc)
         sx4h=t7(neighc)
         vxyd=pxyh+0.5*sx4h*(roci+rockb)
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=vxyd
 63   continue
c     
c     determine upwind nodes and if vapour phase exists
c     
      isl=1
      do 64 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         fid=.5
c     add coding to save upwind position
         if(iad.le.iad_up) then
            vxyd=t8(neighc)
            if(vxyd.lt.0.0) fid=dnwgt
            if(vxyd.gt.0.0) fid=upwgt
            if(t4(neighc).le.0.0) t9(neighc)=fid
            if(t4(neighc).gt.0.0) t9(neighc)=fid
            call setbit(nbits,neighc,upwind_v(iz4m1),fid)
         else
            if(bit(nbits,neighc,upwind_v(iz4m1))) then
               t9(neighc)=1.0
            else
               t9(neighc)=0.0
            endif
         endif
 64   continue
c     
c     form equations
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
            divyckb=divyc(kb)
            divyakb=divya(kb)
            dgcekb=co2_prop(10*neq+kb)
            dgcpkb=co2_prop(11*neq+kb)
            dgcwkb=0.0
            dgcyckb=0.0
            dgcyakb=0.0
            enckb=co2_prop(12*neq+kb)
            decekb=co2_prop(13*neq+kb)
            deckb=co2_prop(14*neq+kb)
            decwkb=0.0
            decyckb=0.0
            decyakb=0.0

            dvpi=-pvxy+0.5*sx4h*dgcpi*(cord(kz,igrav)-cord(iz,igrav))
            dvpkb=pvxy+0.5*sx4h*dgcpkb*(cord(kz,igrav)-cord(iz,igrav))
            dvei=0.5*sx4h*dgcei*(cord(kz,igrav)-cord(iz,igrav))
            dvekb=0.5*sx4h*dgcekb*(cord(kz,igrav)-cord(iz,igrav))
            if(icesd.eq.2) dvei=0.d0
            if(icesd.eq.2) dvekb=0.d0
            vxyf=(fid*divkb+fid1*divi)
            vexyf=(fid*divkb*enckb+fid1*divi*enci)
            vxy=vxyd*vxyf
            vexy=vxyd*vexyf
            dvapi=dvpi*vxyf+vxyd*fid1*divpi
            dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
            dvaei=dvei*vxyf+vxyd*fid1*divei
            dvaekb=dvekb*vxyf+vxyd*fid*divekb
            dvepi=dvpi*vexyf+vxyd*fid1*(divpi*enci+divi*deci)
            dvepkb=dvpkb*vexyf+vxyd*fid*(divpkb*enckb+divkb*deckb)
            dveei=dvei*vexyf+vxyd*fid1*(divei*enci+divi*decei)
            dveekb=dvekb*vexyf+vxyd*fid*(divekb*enckb+divkb*decekb)
c     
c     RJP 07/05/07 added new co2 flux arrays c_axy and c_vxy
c     
            c_vxy(iau+nmatavw)=vxy
            c_vxy(ial+nmatavw)=-vxy

c     derivatives wrt water-rich phase fraction (fw)

            dvawi=vxyd*fid1*divwi
            dvawkb=vxyd*fid*divwkb
            dvewi=vxyd*fid1*(divwi*enci)
            dvewkb=vxyd*fid*(divwkb*enckb)

c     derivatives wrt co2 mass fraction in water-rich phase (yc)
            
            dvyci = 0.5*sx4h*dgcyci*(cord(kz,igrav)-cord(iz,igrav))
            dvyckb = 0.5*sx4h*dgcyckb*(cord(kz,igrav)-cord(iz,igrav))
            dvayci=vxyd*fid1*divyci
            dvayckb=vxyd*fid*divyckb
            dveyci=vxyd*fid1*(divyci*enci)
            dveyckb=vxyd*fid*(divyckb*enckb)

            if(iprtype.eq.4) then
               if(ico2dis(i).eq.0) then
                  dvawi=dvayci
                  dvewi=dveyci
               else
                  if(ices(i).eq.2) then
                     dvapi=dvapi+dvayci*dmol(i)
     &                    +dvayci*dmol(i+neq)
                     dvepi=dvepi+dveyci*dmol(i)
     &                    +dveyci*dmol(i+neq)
                  else
                     dvapi=dvapi+dvayci*dmol(i)
                     dvaei=dvaei+dvayci*dmol(i+neq)
                     dvepi=dvepi+dveyci*dmol(i)
                     dveei=dveei+dveyci*dmol(i+neq)
                  endif
               endif
               if(ico2dis(kb).eq.0) then
                  dvawkb=dvayckb
                  dvewkb=dveyckb
               else
                  if(ices(kb).eq.2) then
                     dvapkb=dvapkb+dvayckb*dmol(kb)
     &                    +dvayckb*dmol(kb+neq)
                     dvepkb=dvepkb+dveyckb*dmol(kb)
     &                    +dveyckb*dmol(kb+neq)
                  else
                     dvapkb=dvapkb+dvayckb*dmol(kb)
                     dvaekb=dvaekb+dvayckb*dmol(kb+neq)
                     dvepkb=dvepkb+dveyckb*dmol(kb)
                     dveekb=dveekb+dveyckb*dmol(kb+neq)
                  endif
               endif
            endif

c     derivatives wrt air mass fraction in water-rich phase (ya)

            dvyai = 0.5*sx4h*dgcyai*(cord(kz,igrav)-cord(iz,igrav))
            dvyakb = 0.5*sx4h*dgcyakb*(cord(kz,igrav)-cord(iz,igrav))
            dvayai=vxyd*fid1*divyai
            dvayakb=vxyd*fid*divyakb
            dveyai=vxyd*fid1*(divyai*enci)
            dveyakb=vxyd*fid*(divyakb*enckb)

            bp(iz+nrhs(3))=bp(iz+nrhs(3))+vxy
            bp(kz+nrhs(3))=bp(kz+nrhs(3))-vxy

            a(jmia+nmat(11))=a(jmia+nmat(11))+dvapi
            a(jmia+nmat(12))=a(jmia+nmat(12))+dvaei
            a(ial+nmat(11))=a(ial+nmat(11))-dvapi
            a(ial+nmat(12))=a(ial+nmat(12))-dvaei
            a(iau+nmat(11))=a(iau+nmat(11))+dvapkb
            a(iau+nmat(12))=a(iau+nmat(12))+dvaekb
            a(jml+nmat(11))=a(jml+nmat(11))-dvapkb
            a(jml+nmat(12))=a(jml+nmat(12))-dvaekb

c     
            bp(iz+nrhs(4))=bp(iz+nrhs(4))+vexy
            bp(kz+nrhs(4))=bp(kz+nrhs(4))-vexy
c            if(iz.eq.1.and.kz.eq.2.and.iad.le.10) then
c             write(ierr,*) 'l ', l,'iad ',iad,'connection = ' , iz,kz
c             write(ierr,*)'vxy,vexy vexy/vxy 1-2 co2 rich-co2'
c             write(ierr,*) vxy, vexy,vexy/vxy
c             write(ierr,*) 'fid,fid1,divkb,divi,enckb,enci'
c             write(ierr,*) fid,fid1,divkb,divi,enckb,enci
c            endif 
            a(jmia+nmat(16))=a(jmia+nmat(16))+dvepi
            a(jmia+nmat(17))=a(jmia+nmat(17))+dveei
            a(ial+nmat(16))=a(ial+nmat(16))-dvepi
            a(ial+nmat(17))=a(ial+nmat(17))-dveei
            a(iau+nmat(16))=a(iau+nmat(16))+dvepkb
            a(iau+nmat(17))=a(iau+nmat(17))+dveekb
            a(jml+nmat(16))=a(jml+nmat(16))-dvepkb
            a(jml+nmat(17))=a(jml+nmat(17))-dveekb
c     
            a(jmia+nmat(13))=a(jmia+nmat(13))+dvawi
            a(ial+nmat(13))=a(ial+nmat(13))-dvawi
            a(iau+nmat(13))=a(iau+nmat(13))+dvawkb
            a(jml+nmat(13))=a(jml+nmat(13))-dvawkb
            a(jmia+nmat(18))=a(jmia+nmat(18))+dvewi
            a(ial+nmat(18))=a(ial+nmat(18))-dvewi
            a(iau+nmat(18))=a(iau+nmat(18))+dvewkb
            a(jml+nmat(18))=a(jml+nmat(18))-dvewkb
c     
            a(jmia+nmat(14))=a(jmia+nmat(14))+dvayci
            a(ial+nmat(14))=a(ial+nmat(14))-dvayci
            a(iau+nmat(14))=a(iau+nmat(14))+dvayckb
            a(jml+nmat(14))=a(jml+nmat(14))-dvayckb
            a(jmia+nmat(19))=a(jmia+nmat(19))+dveyci
            a(ial+nmat(19))=a(ial+nmat(19))-dveyci
            a(iau+nmat(19))=a(iau+nmat(19))+dveyckb
            a(jml+nmat(19))=a(jml+nmat(19))-dveyckb
c     
            a(jmia+nmat(15))=a(jmia+nmat(15))+dvayai
            a(ial+nmat(15))=a(ial+nmat(15))-dvayai
            a(iau+nmat(15))=a(iau+nmat(15))+dvayakb
            a(jml+nmat(15))=a(jml+nmat(15))-dvayakb
            a(jmia+nmat(20))=a(jmia+nmat(20))+dveyai
            a(ial+nmat(20))=a(ial+nmat(20))-dveyai
            a(iau+nmat(20))=a(iau+nmat(20))+dveyakb
            a(jml+nmat(20))=a(jml+nmat(20))-dveyakb

 65      continue
      endif
c     
c     
c     co2-rich liquid phase calculations
c     
      do 66 jm=1,iq
         kb=it8(jm)
         rolkb = co2_prop(kb)
         kz=kb-icd
         neighc=it9(jm)
         pxyh=t10(neighc)
         sx4h=t13(neighc)
         lxyd=pxyh+0.5*sx4h*(roli+rolkb)
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=lxyd
 66   continue
c     
c     determine upwind nodes and if vapour phase exists
c     
      isl=1
      do 67 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         fid=.5
c     add coding to save upwind position
         if(iad.le.iad_up) then
            lxyd=t8(neighc)
            if(lxyd.lt.0.0) fid=dnwgt
            if(lxyd.gt.0.0) fid=upwgt
            if(t4(neighc).le.0.0) t9(neighc)=fid
            if(t4(neighc).gt.0.0) t9(neighc)=fid
            call setbit(nbits,neighc,upwind_v(iz4m1),fid)
         else
            if(bit(nbits,neighc,upwind_v(iz4m1))) then
               t9(neighc)=1.0
            else
               t9(neighc)=0.0
            endif
         endif
 67   continue
c     
c     form equations
c     
      if(isl.ne.0) then
         do jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iau=it11(jm)
            ial=it12(jm)
            jml=nelmdg(kz)-neqp1
            fid=t9(neighc)
            fid1=1.0-fid
            pxyl=t10(neighc)
            pvxy=t4(neighc)
            sx4h=t13(neighc)
            lxyd=t8(neighc)

            dilkb=dil(kb)
            dilpkb=dilp(kb)
            dilekb=dile(kb)
            dilwkb=dilw(kb)
            dilyckb=dilyc(kb)
            dilyakb=dilya(kb)
            dglpkb=co2_prop(2*neq+kb)
            dglekb=co2_prop(neq+kb)
            dglwkb=0.0
            dglyckb=0.0
            dglyakb=0.0
            enlkb=co2_prop(3*neq+kb)
            delkb=co2_prop(5*neq+kb)
            delekb=co2_prop(4*neq+kb)
            delwkb=0.0
            delyckb=0.0
            delyakb=0.0

            dlpi=-pvxy+0.5*sx4h*dglpi*(cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pvxy+0.5*sx4h*dglpkb*(cord(kz,igrav)-cord(iz,igrav))
            dlei=0.5*sx4h*dglei*(cord(kz,igrav)-cord(iz,igrav))
            dlekb=0.5*sx4h*dglekb*(cord(kz,igrav)-cord(iz,igrav))
            if(icesd.eq.2) dlei=0.d0
            if(icesd.eq.2) dlekb=0.d0
            lxyf=(fid*dilkb+fid1*dili)
            lexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            lxy=lxyd*lxyf
            lexy=lxyd*lexyf
            dlapi=dlpi*lxyf+lxyd*fid1*dilpi
            dlapkb=dlpkb*lxyf+lxyd*fid*dilpkb
            dlaei=dlei*lxyf+lxyd*fid1*dilei
            dlaekb=dlekb*lxyf+lxyd*fid*dilekb
            dlepi=dlpi*lexyf+lxyd*fid1*(dilpi*enli+dili*deli)
            dlepkb=dlpkb*lexyf+lxyd*fid*(dilpkb*enlkb+dilkb*delkb)
            dleei=dlei*lexyf+lxyd*fid1*(dilei*enli+dili*delei)
            dleekb=dlekb*lexyf+lxyd*fid*(dilekb*enlkb+dilkb*delekb)

c     RJP 07/25/07 added flux arrays
            c_axy(iau+nmatavw)=lxy
            c_axy(ial+nmatavw)=-lxy

c     derivatives wrt water-rich phase fraction (fw)

            dlawi=lxyd*fid1*dilwi
            dlawkb=lxyd*fid*dilwkb
            dlewi=lxyd*fid1*(dilwi*enli)
            dlewkb=lxyd*fid*(dilwkb*enlkb)

c     derivatives wrt co2 mass fraction in water-rich phase (yc)
            
            dlyci = 0.5*sx4h*dglyci*(cord(kz,igrav)-cord(iz,igrav))
            dlyckb = 0.5*sx4h*dglyckb*(cord(kz,igrav)-cord(iz,igrav))
            dlayci=lxyd*fid1*dilyci
            dlayckb=lxyd*fid*dilyckb
            dleyci=lxyd*fid1*(dilyci*enli)
            dleyckb=lxyd*fid*(dilyckb*enlkb)

            if(iprtype.eq.4) then
               if(ico2dis(i).eq.0) then
                  dlawi=dlayci
                  dlewi=dleyci
               else
                  if(ices(i).eq.2) then
                     dlapi=dlapi+dlayci*dmol(i)
     &                    +dlayci*dmol(i+neq)
                     dlepi=dlepi+dleyci*dmol(i)
     &                    +dleyci*dmol(i+neq)
                  else
                     dlapi=dlapi+dlayci*dmol(i)
                     dlaei=dlaei+dlayci*dmol(i+neq)
                     dlepi=dlepi+dleyci*dmol(i)
                     dleei=dleei+dleyci*dmol(i+neq)
                  endif
               endif
               if(ico2dis(kb).eq.0) then
                  dlawkb=dlayckb
                  dlewkb=dleyckb
               else
                  if(ices(kb).eq.2) then
                     dlapkb=dlapkb+dlayckb*dmol(kb)
     &                    +dlayckb*dmol(kb+neq)
                     dlepkb=dlepkb+dleyckb*dmol(kb)
     &                    +dleyckb*dmol(kb+neq)
                  else
                     dlapkb=dlapkb+dlayckb*dmol(kb)
                     dlaekb=dlaekb+dlayckb*dmol(kb+neq)
                     dlepkb=dlepkb+dleyckb*dmol(kb)
                     dleekb=dleekb+dleyckb*dmol(kb+neq)
                  endif
               endif
            endif

c     derivatives wrt air mass fraction in water-rich phase (ya)

            dlyai = 0.5*sx4h*dglyai*(cord(kz,igrav)-cord(iz,igrav))
            dlyakb = 0.5*sx4h*dglyakb*(cord(kz,igrav)-cord(iz,igrav))
            dlayai=lxyd*fid1*dilyai
            dlayakb=lxyd*fid*dilyakb
            dleyai=lxyd*fid1*(dilyai*enli)
            dleyakb=lxyd*fid*(dilyakb*enlkb)

            bp(iz+nrhs(3))=bp(iz+nrhs(3))+lxy
            bp(kz+nrhs(3))=bp(kz+nrhs(3))-lxy
c            if(iz.eq.1.and.kz.eq.2.and.iad.le.10) then
c            write(ierr,*) 'l ', l,'iad ',iad,'connection = ' , iz,kz
c             write(ierr,*)'lxy,lexy lexy/lxy 1-2 co2 liquid-co2'
c             write(ierr,*) lxy, lexy,lexy/lxy
c             write(ierr,*) 'fid,fid1,dilkb,dili,enlkb,enli'
c             write(ierr,*) fid,fid1,dilkb,dili,enlkb,enli
c            endif 
            a(jmia+nmat(11))=a(jmia+nmat(11))+dlapi
            a(jmia+nmat(12))=a(jmia+nmat(12))+dlaei
            a(ial+nmat(11))=a(ial+nmat(11))-dlapi
            a(ial+nmat(12))=a(ial+nmat(12))-dlaei
            a(iau+nmat(11))=a(iau+nmat(11))+dlapkb
            a(iau+nmat(12))=a(iau+nmat(12))+dlaekb
            a(jml+nmat(11))=a(jml+nmat(11))-dlapkb
            a(jml+nmat(12))=a(jml+nmat(12))-dlaekb

c     
            bp(iz+nrhs(4))=bp(iz+nrhs(4))+lexy
            bp(kz+nrhs(4))=bp(kz+nrhs(4))-lexy
            
            a(jmia+nmat(16))=a(jmia+nmat(16))+dlepi
            a(jmia+nmat(17))=a(jmia+nmat(17))+dleei
            a(ial+nmat(16))=a(ial+nmat(16))-dlepi
            a(ial+nmat(17))=a(ial+nmat(17))-dleei
            a(iau+nmat(16))=a(iau+nmat(16))+dlepkb
            a(iau+nmat(17))=a(iau+nmat(17))+dleekb
            a(jml+nmat(16))=a(jml+nmat(16))-dlepkb
            a(jml+nmat(17))=a(jml+nmat(17))-dleekb
c     
            a(jmia+nmat(13))=a(jmia+nmat(13))+dlawi
            a(ial+nmat(13))=a(ial+nmat(13))-dlawi
            a(iau+nmat(13))=a(iau+nmat(13))+dlawkb
            a(jml+nmat(13))=a(jml+nmat(13))-dlawkb
            a(jmia+nmat(18))=a(jmia+nmat(18))+dlewi
            a(ial+nmat(18))=a(ial+nmat(18))-dlewi
            a(iau+nmat(18))=a(iau+nmat(18))+dlewkb
            a(jml+nmat(18))=a(jml+nmat(18))-dlewkb
c     
            a(jmia+nmat(14))=a(jmia+nmat(14))+dlayci
            a(ial+nmat(14))=a(ial+nmat(14))-dlayci
            a(iau+nmat(14))=a(iau+nmat(14))+dlayckb
            a(jml+nmat(14))=a(jml+nmat(14))-dlayckb
            a(jmia+nmat(19))=a(jmia+nmat(19))+dleyci
            a(ial+nmat(19))=a(ial+nmat(19))-dleyci
            a(iau+nmat(19))=a(iau+nmat(19))+dleyckb
            a(jml+nmat(19))=a(jml+nmat(19))-dleyckb
c     
            a(jmia+nmat(15))=a(jmia+nmat(15))+dlayai
            a(ial+nmat(15))=a(ial+nmat(15))-dlayai
            a(iau+nmat(15))=a(iau+nmat(15))+dlayakb
            a(jml+nmat(15))=a(jml+nmat(15))-dlayakb
            a(jmia+nmat(20))=a(jmia+nmat(20))+dleyai
            a(ial+nmat(20))=a(ial+nmat(20))-dleyai
            a(iau+nmat(20))=a(iau+nmat(20))+dleyakb
            a(jml+nmat(20))=a(jml+nmat(20))-dleyakb

         enddo
      endif
c     
c     add heat conduction
c     
      do jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         heatc=t5(neighc)
         iau=it11(jm)
         ial=it12(jm)
         jml=nelmdg(kz)-neqp1
         heatc=t5(neighc)
         bp(iz+nrhs(4))=bp(iz+nrhs(4))+heatc*(t(kb)-ti)
         bp(kz+nrhs(4))=bp(kz+nrhs(4))-heatc*(t(kb)-ti)

         a(jmia+nmat(16))=a(jmia+nmat(16))-heatc*dtpaco2(i)
         a(jmia+nmat(17))=a(jmia+nmat(17))-heatc*dtpaeco2(i)
         a(ial+nmat(16))=a(ial+nmat(16))+heatc*dtpaco2(i)
         a(ial+nmat(17))=a(ial+nmat(17))+heatc*dtpaeco2(i)
         a(iau+nmat(16))=a(iau+nmat(16))+heatc*dtpaco2(kb)
         a(iau+nmat(17))=a(iau+nmat(17))+heatc*dtpaeco2(kb)
         a(jml+nmat(16))=a(jml+nmat(16))-heatc*dtpaco2(kb)
         a(jml+nmat(17))=a(jml+nmat(17))-heatc*dtpaeco2(kb)
      enddo
c     
c     CO2 diffusion determine upwind nodes based on the dissolved CO2 mass fraction
c     
      do jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         diffc=t17(neighc)
         yckb=yc(kb)
         dfxd=diffc*(yckb-yci)
         t8(neighc)=dfxd
      enddo

      isl=1
      do jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         fid=0.5
         if (iad.le.iad_up) then
            dfxd=t8(neighc)
            if(dfxd.lt.0.d0) fid=dnwgt
            if(dfxd.gt.0.d0) fid=upwgt
            t9(neighc)=fid
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
      enddo

      do 70 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         diffc=t17(neighc)
         iau=it11(jm)
         ial=it12(jm)
         jml=nelmdg(kz)-neqp1
         fid=t9(neighc)
         fid1=1.0-fid
c     mass balance terms
         yckb=yc(kb)
         dgwpkb=wat_prop(neq+kb)
         dgwekb=wat_prop(2*neq+kb)
         dgwwkb=0.0
         dgwyckb=wat_prop(3*neq+kb)
         dgwyakb=wat_prop(4*neq+kb)

         bp(iz+nrhs(3))=bp(iz+nrhs(3))+diffc*(yckb*rowkb-yci*rowi)
         bp(kz+nrhs(3))=bp(kz+nrhs(3))-diffc*(yckb*rowkb-yci*rowi)
         if(iprtype.eq.4) then
            if(ico2dis(i).eq.0) then
               a(jmia+nmat(11))=a(jmia+nmat(11))-diffc*yci*dgwpi
               a(ial+nmat(11))=a(ial+nmat(11))+diffc*yci*dgwpi
               a(jmia+nmat(12))=a(jmia+nmat(12))-diffc*yci*dgwei
               a(ial+nmat(12))=a(ial+nmat(12))+diffc*yci*dgwei
               a(jmia+nmat(13))=a(jmia+nmat(13))-diffc*(yci*dgwyci+rowi)
               a(ial+nmat(13))=a(ial+nmat(13))+diffc*(yci*dgwyci+rowi)
               a(jmia+nmat(14))=a(jmia+nmat(14))-diffc*(yci*dgwyci+rowi)
               a(ial+nmat(14))=a(ial+nmat(14))+diffc*(yci*dgwyci+rowi)
               a(jmia+nmat(15))=a(jmia+nmat(15))-diffc*yci*dgwyai
               a(ial+nmat(15))=a(ial+nmat(15))+diffc*yci*dgwyai
            else
               if(ices(i).eq.2) then
                  a(jmia+nmat(11))=a(jmia+nmat(11))-diffc*yci*(dgwpi+
     &                 dgwyci*dmol(i)+dgwyci*dmol(i+neq))
                  a(ial+nmat(11))=a(ial+nmat(11))+diffc*yci*(dgwpi+
     &                 dgwyci*dmol(i)+dgwyci*dmol(i+neq))
                  a(jmia+nmat(13))=a(jmia+nmat(13))-diffc*yci*dgwwi
                  a(ial+nmat(13))=a(ial+nmat(13))+diffc*yci*dgwwi
               a(jmia+nmat(14))=a(jmia+nmat(14))-diffc*(yci*dgwyci+rowi)
                a(ial+nmat(14))=a(ial+nmat(14))+diffc*(yci*dgwyci+rowi)
                  a(jmia+nmat(15))=a(jmia+nmat(15))-diffc*yci*dgwyai
                  a(ial+nmat(15))=a(ial+nmat(15))+diffc*yci*dgwyai
               else
                  a(jmia+nmat(11))=a(jmia+nmat(11))-diffc*yci*(dgwpi+
     &                 dgwyci*dmol(i))
                  a(ial+nmat(11))=a(ial+nmat(11))+diffc*yci*(dgwpi+
     &                 dgwyci*dmol(i))
                  a(jmia+nmat(12))=a(jmia+nmat(12))-diffc*yci*(dgwei+
     &                 dgwyci*dmol(i+neq))
                  a(ial+nmat(12))=a(ial+nmat(12))+diffc*yci*(dgwei+
     &                 dgwyci*dmol(i+neq))
                  a(jmia+nmat(13))=a(jmia+nmat(13))-diffc*yci*dgwwi
                  a(ial+nmat(13))=a(ial+nmat(13))+diffc*yci*dgwwi
                  a(jmia+nmat(14))=a(jmia+nmat(14))-
     &                 diffc*(yci*dgwyci+rowi)
                  a(ial+nmat(14))=a(ial+nmat(14))+
     &                 diffc*(yci*dgwyci+rowi)
                  a(jmia+nmat(15))=a(jmia+nmat(15))-diffc*yci*dgwyai
                  a(ial+nmat(15))=a(ial+nmat(15))+diffc*yci*dgwyai
               endif
            endif
         else
            a(jmia+nmat(11))=a(jmia+nmat(11))-diffc*yci*dgwpi
            a(ial+nmat(11))=a(ial+nmat(11))+diffc*yci*dgwpi
            a(jmia+nmat(12))=a(jmia+nmat(12))-diffc*yci*dgwei
            a(ial+nmat(12))=a(ial+nmat(12))+diffc*yci*dgwei
            a(jmia+nmat(13))=a(jmia+nmat(13))-diffc*yci*dgwwi
            a(ial+nmat(13))=a(ial+nmat(13))+diffc*yci*dgwwi
            a(jmia+nmat(14))=a(jmia+nmat(14))-diffc*(yci*dgwyci+rowi)
            a(ial+nmat(14))=a(ial+nmat(14))+diffc*(yci*dgwyci+rowi)
            a(jmia+nmat(15))=a(jmia+nmat(15))-diffc*yci*dgwyai
            a(ial+nmat(15))=a(ial+nmat(15))+diffc*yci*dgwyai
         endif

         if(iprtype.eq.4) then
            if(ico2dis(kb).eq.0) then
               a(iau+nmat(11))=a(iau+nmat(11))+diffc*yckb*dgwpkb
               a(jml+nmat(11))=a(jml+nmat(11))-diffc*yckb*dgwpkb
               a(iau+nmat(12))=a(iau+nmat(12))+diffc*yckb*dgwekb
               a(jml+nmat(12))=a(jml+nmat(12))-diffc*yckb*dgwekb
               a(iau+nmat(13))=a(iau+nmat(13))+diffc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(13))=a(jml+nmat(13))-diffc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(14))=a(iau+nmat(14))+diffc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(14))=a(jml+nmat(14))-diffc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(15))=a(iau+nmat(15))+diffc*yckb*dgwyakb
               a(jml+nmat(15))=a(jml+nmat(15))-diffc*yckb*dgwyakb
            else
			if(ices(kb).eq.2) then
               a(iau+nmat(11))=a(iau+nmat(11))+diffc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb)+dgwyckb*dmol(kb+neq))
               a(jml+nmat(11))=a(jml+nmat(11))-diffc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb)+dgwyckb*dmol(kb+neq))
               a(iau+nmat(13))=a(iau+nmat(13))+diffc*yckb*dgwwkb
               a(jml+nmat(13))=a(jml+nmat(13))-diffc*yckb*dgwwkb
               a(iau+nmat(14))=a(iau+nmat(14))+diffc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(14))=a(jml+nmat(14))-diffc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(15))=a(iau+nmat(15))+diffc*yckb*dgwyakb
               a(jml+nmat(15))=a(jml+nmat(15))-diffc*yckb*dgwyakb
			else
	         a(iau+nmat(11))=a(iau+nmat(11))+diffc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))
               a(jml+nmat(11))=a(jml+nmat(11))-diffc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))
               a(iau+nmat(12))=a(iau+nmat(12))+diffc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))
               a(jml+nmat(12))=a(jml+nmat(12))-diffc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))
               a(iau+nmat(13))=a(iau+nmat(13))+diffc*yckb*dgwwkb
               a(jml+nmat(13))=a(jml+nmat(13))-diffc*yckb*dgwwkb
               a(iau+nmat(14))=a(iau+nmat(14))+diffc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(14))=a(jml+nmat(14))-diffc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(15))=a(iau+nmat(15))+diffc*yckb*dgwyakb
               a(jml+nmat(15))=a(jml+nmat(15))-diffc*yckb*dgwyakb
			endif
            endif
         else
            a(iau+nmat(11))=a(iau+nmat(11))+diffc*yckb*dgwpkb
            a(jml+nmat(11))=a(jml+nmat(11))-diffc*yckb*dgwpkb
            a(iau+nmat(12))=a(iau+nmat(12))+diffc*yckb*dgwekb
            a(jml+nmat(12))=a(jml+nmat(12))-diffc*yckb*dgwekb
            a(iau+nmat(13))=a(iau+nmat(13))+diffc*yckb*dgwwkb
            a(jml+nmat(13))=a(jml+nmat(13))-diffc*yckb*dgwwkb
            a(iau+nmat(14))=a(iau+nmat(14))+diffc*(yckb*dgwyckb+rowkb)
            a(jml+nmat(14))=a(jml+nmat(14))-diffc*(yckb*dgwyckb+rowkb)
            a(iau+nmat(15))=a(iau+nmat(15))+diffc*yckb*dgwyakb
            a(jml+nmat(15))=a(jml+nmat(15))-diffc*yckb*dgwyakb
            
         endif

c     energy balance terms
         if(ices(kb).eq.3) then
            enlkb=co2_prop(12*neq+kb)
            delkb=co2_prop(14*neq+kb)
            delekb=co2_prop(13*neq+kb)
         elseif(ices(kb).eq.1) then
            enlkb=co2_prop(3*neq+kb)
            delkb=co2_prop(5*neq+kb)
            delekb=co2_prop(4*neq+kb)
         elseif(ices(kb).eq.4) then
            enlkb=co2_prop(3*neq+kb)
            delkb=co2_prop(5*neq+kb)
            delekb=co2_prop(4*neq+kb)
         endif
         delwkb=0.d0
         delyckb=0.d0
         delyakb=0.d0
         if(ices(i).eq.3) then
            hfc=fid*enlkb+fid1*enci
            dhfcp=fid*delkb+fid1*deci
            dhfct=fid*delekb+fid1*decei
         elseif(ices(i).eq.1) then
            hfc=fid*enlkb+fid1*enli
            dhfcp=fid*delkb+fid1*deli
            dhfct=fid*delekb+fid1*delei
         elseif(ices(i).eq.4) then
            hfc=fid*enlkb+fid1*enli
            dhfcp=fid*delkb+fid1*deli
            dhfct=fid*delekb+fid1*delei
         endif
C     hfc=0.d0
C     dhfcp=0.d0
C     dhfct=0.d0	
         dhfcw=0.d0
         dhfcyc=0.d0
         dhfcya=0.d0
         bp(iz+nrhs(4))=bp(iz+nrhs(4))+diffc*hfc*(rowkb*yckb-yci*rowi)
         bp(kz+nrhs(4))=bp(kz+nrhs(4))-diffc*hfc*(rowkb*yckb-yci*rowi)
         if(iprtype.eq.4) then
            if(ico2dis(i).eq.0) then
               a(jmia+nmat(16))=a(jmia+nmat(16))-diffc*(hfc*yci*dgwpi
     &              +dhfcp*yci*rowi)
               a(ial+nmat(16))=a(ial+nmat(16))+diffc*(hfc*yci*dgwpi
     &              +dhfcp*yci*rowi)
               a(jmia+nmat(17))=a(jmia+nmat(17))-diffc*(hfc*yci*dgwei
     &              +dhfct*yci*rowi)
               a(ial+nmat(17))=a(ial+nmat(17))+diffc*(hfc*yci*dgwei
     &              +dhfct*yci*rowi)
               a(jmia+nmat(18))=a(jmia+nmat(18))-diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(ial+nmat(18))=a(ial+nmat(18))+diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(jmia+nmat(19))=a(jmia+nmat(19))-diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(ial+nmat(19))=a(ial+nmat(19))+diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(jmia+nmat(20))=a(jmia+nmat(20))-diffc*hfc*yci*dgwyai
               a(ial+nmat(20))=a(ial+nmat(20))+diffc*hfc*yci*dgwyai
            else
			if(ices(i).eq.2) then
               a(jmia+nmat(16))=a(jmia+nmat(16))-diffc*(hfc*yci*(dgwpi+
     &              dgwyci*dmol(i))+dhfcp*yci*rowi)
               a(ial+nmat(16))=a(ial+nmat(16))+diffc*(hfc*yci*(dgwpi+
     &              dgwyci*dmol(i))+dhfcp*yci*rowi)
               a(jmia+nmat(17))=a(jmia+nmat(17))-diffc*(hfc*yci*(dgwei+
     &              dgwyci*dmol(i+neq))+dhfct*yci*rowi)
               a(ial+nmat(17))=a(ial+nmat(17))+diffc*(hfc*yci*(dgwei+
     &              dgwyci*dmol(i+neq))+dhfct*yci*rowi)
               a(jmia+nmat(18))=a(jmia+nmat(18))-diffc*hfc*yci*dgwwi
               a(ial+nmat(18))=a(ial+nmat(18))+diffc*hfc*yci*dgwwi
               a(jmia+nmat(19))=a(jmia+nmat(19))-diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(ial+nmat(19))=a(ial+nmat(19))+diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(jmia+nmat(20))=a(jmia+nmat(20))-diffc*hfc*yci*dgwyai
               a(ial+nmat(20))=a(ial+nmat(20))+diffc*hfc*yci*dgwyai
			else
	         a(jmia+nmat(16))=a(jmia+nmat(16))-diffc*(hfc*yci*(dgwpi+
     &              dgwyci*dmol(i))+dhfcp*yci*rowi)
               a(ial+nmat(16))=a(ial+nmat(16))+diffc*(hfc*yci*(dgwpi+
     &              dgwyci*dmol(i))+dhfcp*yci*rowi)
               a(jmia+nmat(17))=a(jmia+nmat(17))-diffc*(hfc*yci*(dgwei+
     &              dgwyci*dmol(i+neq))+dhfct*yci*rowi)
               a(ial+nmat(17))=a(ial+nmat(17))+diffc*(hfc*yci*(dgwei+
     &              dgwyci*dmol(i+neq))+dhfct*yci*rowi)
               a(jmia+nmat(18))=a(jmia+nmat(18))-diffc*hfc*yci*dgwwi
               a(ial+nmat(18))=a(ial+nmat(18))+diffc*hfc*yci*dgwwi
               a(jmia+nmat(19))=a(jmia+nmat(19))-diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(ial+nmat(19))=a(ial+nmat(19))+diffc*hfc*(yci*dgwyci+
     &              rowi)
               a(jmia+nmat(20))=a(jmia+nmat(20))-diffc*hfc*yci*dgwyai
               a(ial+nmat(20))=a(ial+nmat(20))+diffc*hfc*yci*dgwyai
			endif
            endif
         else
            a(jmia+nmat(16))=a(jmia+nmat(16))-diffc*(hfc*yci*dgwpi
     &           +dhfcp*yci*rowi)
            a(ial+nmat(16))=a(ial+nmat(16))+diffc*(hfc*yci*dgwpi
     &           +dhfcp*yci*rowi)
            a(jmia+nmat(17))=a(jmia+nmat(17))-diffc*(hfc*yci*dgwei
     &           +dhfct*yci*rowi)
            a(ial+nmat(17))=a(ial+nmat(17))+diffc*(hfc*yci*dgwei
     &           +dhfct*yci*rowi)
            a(jmia+nmat(18))=a(jmia+nmat(18))-diffc*hfc*yci*dgwwi
            a(ial+nmat(18))=a(ial+nmat(18))+diffc*hfc*yci*dgwwi
            a(jmia+nmat(19))=a(jmia+nmat(19))-diffc*hfc*(yci*dgwyci+
     &              rowi)
            a(ial+nmat(19))=a(ial+nmat(19))+diffc*hfc*(yci*dgwyci+
     &              rowi)
            a(jmia+nmat(20))=a(jmia+nmat(20))-diffc*hfc*yci*dgwyai
            a(ial+nmat(20))=a(ial+nmat(20))+diffc*hfc*yci*dgwyai
         endif

         if(iprtype.eq.4) then
            if(ico2dis(kb).eq.0) then
               a(iau+nmat(16))=a(iau+nmat(16))+diffc*(hfc*yckb*dgwpkb
     &              +dhfcp*yckb*rowkb)
               a(iau+nmat(17))=a(iau+nmat(17))+diffc*(hfc*yckb*dgwekb
     &              +dhfct*yckb*rowkb)
               a(jml+nmat(16))=a(jml+nmat(16))-diffc*(hfc*yckb*dgwpkb
     &              +dhfcp*yckb*rowkb)
               a(jml+nmat(17))=a(jml+nmat(17))-diffc*(hfc*yckb*dgwekb
     &              +dhfct*yckb*rowkb)
               a(iau+nmat(18))=a(iau+nmat(18))+diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(18))=a(jml+nmat(18))-diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(19))=a(iau+nmat(19))+diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(19))=a(jml+nmat(19))-diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(20))=a(iau+nmat(20))+diffc*hfc*yckb*dgwyakb
               a(jml+nmat(20))=a(jml+nmat(20))-diffc*hfc*yckb*dgwyakb
            else
			if(ices(kb).eq.2) then
               a(iau+nmat(16))=a(iau+nmat(16))+diffc*(hfc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))+dhfcp*yckb*rowkb)
               a(iau+nmat(17))=a(iau+nmat(17))+diffc*(hfc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))+dhfct*yckb*rowkb)
               a(jml+nmat(16))=a(jml+nmat(16))-diffc*(hfc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))+dhfcp*yckb*rowkb)
               a(jml+nmat(17))=a(jml+nmat(17))-diffc*(hfc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))+dhfct*yckb*rowkb)
               a(iau+nmat(18))=a(iau+nmat(18))+diffc*hfc*yckb*dgwwkb
               a(jml+nmat(18))=a(jml+nmat(18))-diffc*hfc*yckb*dgwwkb
               a(iau+nmat(19))=a(iau+nmat(19))+diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(19))=a(jml+nmat(19))-diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(20))=a(iau+nmat(20))+diffc*hfc*yckb*dgwyakb
               a(jml+nmat(20))=a(jml+nmat(20))-diffc*hfc*yckb*dgwyakb
			else
               a(iau+nmat(16))=a(iau+nmat(16))+diffc*(hfc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))+dhfcp*yckb*rowkb)
               a(iau+nmat(17))=a(iau+nmat(17))+diffc*(hfc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))+dhfct*yckb*rowkb)
               a(jml+nmat(16))=a(jml+nmat(16))-diffc*(hfc*yckb*(dgwpkb+
     &              dgwyckb*dmol(kb))+dhfcp*yckb*rowkb)
               a(jml+nmat(17))=a(jml+nmat(17))-diffc*(hfc*yckb*(dgwekb+
     &              dgwyckb*dmol(kb+neq))+dhfct*yckb*rowkb)
               a(iau+nmat(18))=a(iau+nmat(18))+diffc*hfc*yckb*dgwwkb
               a(jml+nmat(18))=a(jml+nmat(18))-diffc*hfc*yckb*dgwwkb
               a(iau+nmat(19))=a(iau+nmat(19))+diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(jml+nmat(19))=a(jml+nmat(19))-diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
               a(iau+nmat(20))=a(iau+nmat(20))+diffc*hfc*yckb*dgwyakb
               a(jml+nmat(20))=a(jml+nmat(20))-diffc*hfc*yckb*dgwyakb
			endif
            endif
         else
            a(iau+nmat(16))=a(iau+nmat(16))+diffc*(hfc*yckb*dgwpkb
     &           +dhfcp*yckb*rowkb)
            a(jml+nmat(16))=a(jml+nmat(16))-diffc*(hfc*yckb*dgwpkb
     &           +dhfcp*yckb*rowkb)
            a(iau+nmat(17))=a(iau+nmat(17))+diffc*(hfc*yckb*dgwekb
     &           +dhfct*yckb*rowkb)
            a(jml+nmat(17))=a(jml+nmat(17))-diffc*(hfc*yckb*dgwekb
     &           +dhfct*yckb*rowkb)
            a(iau+nmat(18))=a(iau+nmat(18))+diffc*hfc*yckb*dgwwkb
            a(jml+nmat(18))=a(jml+nmat(18))-diffc*hfc*yckb*dgwwkb
            a(iau+nmat(19))=a(iau+nmat(19))+diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
            a(jml+nmat(19))=a(jml+nmat(19))-diffc*hfc*(yckb*dgwyckb+
     &              rowkb)
            a(iau+nmat(20))=a(iau+nmat(20))+diffc*hfc*yckb*dgwyakb
            a(jml+nmat(20))=a(jml+nmat(20))-diffc*hfc*yckb*dgwyakb
         endif

 70   continue

c     add accumulation terms
      bp(iz+nrhs(3))=bp(iz+nrhs(3))+sx1d*denco2i(i)+skco2(i)
      bp(iz+nrhs(4))=bp(iz+nrhs(4))+sx1d*deneco2i(i)+qhco2(i)
      a(jmia+nmat(11))=a(jmia+nmat(11))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(12))=a(jmia+nmat(12))+sx1d*dmef(i)+dqt(i)
      a(jmia+nmat(13))=a(jmia+nmat(13))+sx1d*dmwf(i)+dqw(i) 
      a(jmia+nmat(14))=a(jmia+nmat(14))+sx1d*dmycf(i)+dqyc(i)
      a(jmia+nmat(15))=a(jmia+nmat(15))+sx1d*dmyaf(i)+dqya(i)
      a(jmia+nmat(16))=a(jmia+nmat(16))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(17))=a(jmia+nmat(17))+sx1d*deef(i)+deqh(i)
      a(jmia+nmat(18))=a(jmia+nmat(18))+sx1d*dewf(i)+dqhw(i)
      a(jmia+nmat(19))=a(jmia+nmat(19))+sx1d*deycf(i)+dqhyc(i)
      a(jmia+nmat(20))=a(jmia+nmat(20))+sx1d*deyaf(i)+dqhya(i)
      
      
      return
      end
