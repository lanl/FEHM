      subroutine  geneq_h2o_co2  ( i )
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
CD1  and mass transfer for one componet of a mixture(water) for CO2
CD1  problem. This subroutine is modified version of geneq_h2o.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 21-FEB-07    R. Pawar               Initial implementation
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
CD4  
CD4  
CD4
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
      integer neighc, neqp1, nmatavw, icesd
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
      real*8  dpvtco2,dleico2,dlekbco2

      real*8 rowi, dgwpi, dgwei, dgwwi, dgwyci, dgwyai
      real*8 roci, dgcpi, dgcei, dgcwi, dgcyci, dgcyai
      real*8 enwi, dewei, dewwi, dewyci, dewyai
      real*8 enci, decei, decwi, decyci, decyai
      real*8 dilyci, dilyai, dewi, deci
      real*8 divyci, divyai
      real*8 rowkb, dilyckb, dilyakb, dgwpkb, dgwekb, dgwwkb
      real*8 dgwyckb, dgwyakb, enwkb, dewkb, dewekb, dewwkb
      real*8 dewyckb, dewyakb, dlwi, dlwkb, dlyci, dlyckb
      real*8 dlyakb, dleyci, dleyckb, dlyai, dlayai
      real*8 dlayakb, dleyai, delyakb, rockb, divyckb
      real*8 dlayci, dlayckb, dleyakb, divyakb, dgcpkb, dgcekb
      real*8 dgcwkb, dgcycka, dgcyakb, enckb, deckb, decekb
      real*8 decwkb, decyckb, decyakb, dvyci, dvyckb, dvayci
      real*8 dvayckb, dveyci, dveyckb, dvyai, dvyakb, dvayai
      real*8 dvayakb, dveyai, dveyakb, dgcyckb

      parameter(dis_tol=1.d-12)
c     
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
      pvii=phico2(i)
      phii=phi(i)
      dpvtw=dpcpw(i)
      dpvtco2=dpcgw(i) 

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

      if(icesd.ne.2) then
         roci=co2_prop(i)
         dgcei=co2_prop(neq+i)
         dgcpi=co2_prop(2*neq+i)
         enci=co2_prop(3*neq+i)
         deci=co2_prop(5*neq+i)
         decei=co2_prop(4*neq+i)
      else
         roci=co2_prop(9*neq+i)
         dgcei=co2_prop(10*neq+i)
         dgcpi=co2_prop(11*neq+i)
         enci=co2_prop(12*neq+i)
         deci=co2_prop(14*neq+i)
         decei=co2_prop(13*neq+i)
         dgwei=0.d0
         dewei=0.d0
      endif
      dgcwi=0.0
      dgcyci=0.0
      dgcyai=0.0
      decwi=0.0
      decyci=0.0
      decyai=0.0

      dili=diw(i)
      dilpi=diwp(i)
      dilei=diwe(i)
      dilwi=diww(i)
      dilyci=diwyc(i)
      dilyai=diwya(i)

      divi=div(i)
      divpi=divp(i)
      divei=dive(i)
      divwi=divw(i)
      divyci=divyc(i)
      divyai=divya(i)

      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      ti=t(i)
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
      if(ico2.lt.0.and.ice.ne.0) then
         a_axy(jmia+nmatavw)=sk(i)*fracw(i)
c     Take care of souce/sink term
c     If this is an isothermal air-water simulation
         a_vxy(jmia+nmatavw)=sk(i)*(1-fracw(i))
      else
         a_axy(jmia+nmatavw)=sk(i)
         a_vxy(jmia+nmatavw)=sk(i)
      endif

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
            pvikb=phico2(kb)
            phikb=phi(kb)
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
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/(delx2/sx2t+dely2/sx3t+delz2/sxzt)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c                    
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
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
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            pvikb=phico2(kb)
            phikb=phi(kb)
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
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/(delx2/sx2t+dely2/sx3t)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
 69      continue
c     
      endif
c     
c     liquid phase calculations
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
            dilkb=diw(kb)
            dilpkb=diwp(kb)
            dilekb=diwe(kb)
            dilwkb=diww(kb)
            dilyckb=diwyc(kb)
            dilyakb=diwya(kb)

            dgwpkb=wat_prop(neq+kb)
            dgwekb=wat_prop(2*neq+kb)
            dgwwkb=0.0
            dgwyckb=wat_prop(3*neq+kb)
            dgwyakb=wat_prop(4*neq+kb)
            enwkb=wat_prop(5*neq+kb)
            dewkb=wat_prop(6*neq+kb)
            dewekb=wat_prop(7*neq+kb)
            if(icesd.eq.2) dgwekb=0.d0
            if(icesd.eq.2) dewekb=0.d0
            dewwkb=0.0
            dewyckb=wat_prop(11*neq+kb)
            dewyakb=0.0

            dlpi=-pxy+0.5*sx4d*dgwpi*(cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pxy+0.5*sx4d*dgwpkb*(cord(kz,igrav)-cord(iz,igrav))
c     For single phase CO2 derivative of cap. pr. wrt. temp, dpvti, is zero
c     For double phase CO2 derivative of cap. pr. wrt. temp, dpvti, is derivative
c     of cap. pr. wrt. CO2 gas saturation.
            dlei=0.5*sx4d*dgwei*(cord(kz,igrav)-cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgwekb
     &           *(cord(kz,igrav)-cord(iz,igrav))
            if(icesd.eq.2) dlei=0.d0
            if(icesd.eq.2) dlekb=0.d0
            axyf=(fid*dilkb+fid1*dili)
            aexyf=(fid*dilkb*enwkb+fid1*dili*enwi)
            axy=axyd*axyf
            aexy=axyd*aexyf
            dlapi=dlpi*axyf+axyd*fid1*dilpi
            dlapkb=dlpkb*axyf+axyd*fid*dilpkb
            dlaei=dlei*axyf+axyd*fid1*dilei
            dlaekb=dlekb*axyf+axyd*fid*dilekb
            dlepi=dlpi*aexyf+axyd*fid1*(dilpi*enwi+dili*dewi)
            dlepkb=dlpkb*aexyf+axyd*fid*(dilpkb*enwkb+dilkb*dewkb)
            dleei=dlei*aexyf+axyd*fid1*(dilei*enwi+dili*dewei)
            dleekb=dlekb*aexyf+axyd*fid*(dilekb*enwkb+dilkb*dewekb)

c     derivatives wrt water-rich phase fraction (fracw)

            dlwi = pxy*dpvtw
            dlwkb = -pxy*dpcpw(kb) 
            dlawi=dlwi*axyf+axyd*fid1*dilwi
            dlawkb=dlwkb*axyf+axyd*fid*dilwkb
            dlewi=dlwi*aexyf+axyd*fid1*(dilwi*enwi+dili*dewi)
            dlewkb=dlwkb*aexyf+axyd*fid*(dilwkb*enwkb+dilkb*dewkb)

c     derivatives wrt co2 mass fraction in water-rich phase (yc)

            dlyci = 0.5*sx4d*dgwyci*(cord(kz,igrav)-cord(iz,igrav))
            dlyckb = 0.5*sx4d*dgwyckb*(cord(kz,igrav)-cord(iz,igrav))
            dlayci = dlyci*axyf+axyd*fid1*dilyci
            dlayckb = dlyckb*axyf+axyd*fid*dilyckb
            dleyci = dlyci*aexyf+axyd*fid1*(dilyci*enwi+dili*dewyci)
            dleyckb = dlyckb*aexyf+axyd*fid*(dilyckb*enwkb+
     &           dilkb*dewyckb)

            if(iprtype.eq.4) then
               if(ico2dis(i).eq.0) then
                  dlawi=dlayci
                  dlewi=dleyci
               else
				if(ices(i).eq.2) then
					dlapi=dlapi+dlayci*dmol(i)
     &				+dlayci*dmol(i+neq)
					dlepi=dlepi+dleyci*dmol(i)
     &				+dleyci*dmol(i+neq)
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
     &				+dlayckb*dmol(kb+neq)
					dlepkb=dlepkb+dleyckb*dmol(kb)
     &				+dleyckb*dmol(kb+neq)
				else
					dlapkb=dlapkb+dlayckb*dmol(kb)
					dlaekb=dlaekb+dlayckb*dmol(kb+neq)
					dlepkb=dlepkb+dleyckb*dmol(kb)
					dleekb=dleekb+dleyckb*dmol(kb+neq)
				endif
               endif
            endif
c     derivatives wrt air mass fraction in water-rich phase (ya)

            dlyai = 0.5*sx4d*dgwyai*(cord(kz,igrav)-cord(iz,igrav))
            dlyakb = 0.5*sx4d*dgwyakb*(cord(kz,igrav)-cord(iz,igrav))
            dlayai = dlyai*axyf+axyd*fid1*dilyai
            dlayakb = dlyakb*axyf+axyd*fid*dilyakb
            dleyai = dlyai*aexyf+axyd*fid1*(dilyai*enwi+dili*dewyai)
            dleyakb = dlyakb*aexyf+axyd*fid*(dilyakb*enwkb+dilkb*
     &           dewyakb)

c     
c     changed by avw -- entered here by seh
c     
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
c            if(iz.eq.1.and.kz.eq.2.and.iad.le.10) then
c             write(ierr,*) 'l ', l,'iad ',iad,'connection = ' , iz,kz
c             write(ierr,*)'axy,aexy aexy/axy 1-2 water rich-h2o'
c             write(ierr,*) axy, aexy,aexy/axy
c             write(ierr,*) 'fid,fid1,dilkb,dili,enwkb,enwi'
c             write(ierr,*) fid,fid1,dilkb,dili,enwkb,enwi
c            endif
            a(jmia+nmat(6))=a(jmia+nmat(6))+dlepi
            a(jmia+nmat(7))=a(jmia+nmat(7))+dleei
            a(ial+nmat(6))=a(ial+nmat(6))-dlepi
            a(ial+nmat(7))=a(ial+nmat(7))-dleei
            a(iau+nmat(6))=a(iau+nmat(6))+dlepkb
            a(iau+nmat(7))=a(iau+nmat(7))+dleekb
            a(jml+nmat(6))=a(jml+nmat(6))-dlepkb
            a(jml+nmat(7))=a(jml+nmat(7))-dleekb
c     
            a(jmia+nmat(3))=a(jmia+nmat(3))+dlawi
            a(ial+nmat(3))=a(ial+nmat(3))-dlawi
            a(iau+nmat(3))=a(iau+nmat(3))+dlawkb
            a(jml+nmat(3))=a(jml+nmat(3))-dlawkb
            a(jmia+nmat(8))=a(jmia+nmat(8))+dlewi
            a(ial+nmat(8))=a(ial+nmat(8))-dlewi
            a(iau+nmat(8))=a(iau+nmat(8))+dlewkb
            a(jml+nmat(8))=a(jml+nmat(8))-dlewkb
c     
            a(jmia+nmat(4))=a(jmia+nmat(4))+dlayci
            a(ial+nmat(4))=a(ial+nmat(4))-dlayci
            a(iau+nmat(4))=a(iau+nmat(4))+dlayckb
            a(jml+nmat(4))=a(jml+nmat(4))-dlayckb
            a(jmia+nmat(9))=a(jmia+nmat(9))+dleyci
            a(ial+nmat(9))=a(ial+nmat(9))-dleyci
            a(iau+nmat(9))=a(iau+nmat(9))+dleyckb
            a(jml+nmat(9))=a(jml+nmat(9))-dleyckb
c     
            a(jmia+nmat(5))=a(jmia+nmat(5))+dlayai
            a(ial+nmat(5))=a(ial+nmat(5))-dlayai
            a(iau+nmat(5))=a(iau+nmat(5))+dlayakb
            a(jml+nmat(5))=a(jml+nmat(5))-dlayakb
            a(jmia+nmat(10))=a(jmia+nmat(10))+dleyai
            a(ial+nmat(10))=a(ial+nmat(10))-dleyai
            a(iau+nmat(10))=a(iau+nmat(10))+dleyakb
            a(jml+nmat(10))=a(jml+nmat(10))-dleyakb
c     
 62      continue
      end if
c     
c     vapour phase calculations
c     
      do 63 jm=1,iq
         kb=it8(jm)
         if(ices(kb).ne.2) then
            rockb=co2_prop(kb)
         else
            rockb=co2_prop(9*neq+kb)
         endif
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

            if(ices(kb).ne.2) then
               dgcpkb=co2_prop(2*neq+kb)
               dgcekb=co2_prop(neq+kb)
               dgcwkb=0.0
               dgcyckb=0.0
               dgcyakb=0.0
               enckb=co2_prop(3*neq+kb)
               deckb=co2_prop(5*neq+kb)
               decekb=co2_prop(4*neq+kb)
c               enwkb=wat_prop(5*neq+kb)
c               dewkb=wat_prop(6*neq+kb)
c               dewekb=wat_prop(7*neq+kb)
c               dewwkb=0.0
c               dewyckb=0.0
c               dewyakb=0.0
               decwkb=0.0
               decyckb=0.0
               decyakb=0.0
            else
               dgcpkb=co2_prop(11*neq+kb)
               dgcekb=co2_prop(10*neq+kb)
               dgcwkb=0.0
               dgcyckb=0.0
               dgcyakb=0.0
               enckb=co2_prop(12*neq+kb)
               deckb=co2_prop(14*neq+kb)
               decekb=co2_prop(13*neq+kb)
               decwkb=0.0
               decyckb=0.0
               decyakb=0.0
               enwkb=wat_prop(5*neq+kb)
               dewkb=wat_prop(6*neq+kb)
               dewekb=wat_prop(7*neq+kb)
               dewwkb=0.0
               dewyckb=0.0
               dewyakb=0.0
            endif

            dvpi=-pvxy+0.5*sx4h*dgcpi*(cord(kz,igrav)-cord(iz,igrav))
            dvpkb=pvxy+0.5*sx4h*dgcpkb*(cord(kz,igrav)-cord(iz,igrav))
            dvei=0.5*sx4h*dgcei*(cord(kz,igrav)-cord(iz,igrav))
            dvekb=0.5*sx4h*dgcekb*(cord(kz,igrav)-cord(iz,igrav))
            if(icesd.eq.2) dvei = 0.d0
            if(icesd.eq.2) dvekb = 0.d0
            vxyf=(fid*divkb+fid1*divi)
            vexyf=(fid*divkb*enckb+fid1*divi*enci)
c Raj 040111       vexyf=(fid*divkb*enwkb+fid1*divi*enwi)
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
c     changed by avw -- entered here by seh
            a_vxy(iau+nmatavw)=vxy
            a_vxy(ial+nmatavw)=-vxy

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
     &				+dvayci*dmol(i+neq)
					dvepi=dvepi+dveyci*dmol(i)
     &				+dveyci*dmol(i+neq)
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
     &				+dvayckb*dmol(kb+neq)
					dvepkb=dvepkb+dveyckb*dmol(kb)
     &				+dveyckb*dmol(kb+neq)
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

            a(jmia+nmat(6))=a(jmia+nmat(6))+dvepi
            a(jmia+nmat(7))=a(jmia+nmat(7))+dveei
            a(ial+nmat(6))=a(ial+nmat(6))-dvepi
            a(ial+nmat(7))=a(ial+nmat(7))-dveei
            a(iau+nmat(6))=a(iau+nmat(6))+dvepkb
            a(iau+nmat(7))=a(iau+nmat(7))+dveekb
            a(jml+nmat(6))=a(jml+nmat(6))-dvepkb
            a(jml+nmat(7))=a(jml+nmat(7))-dveekb
c     
            a(jmia+nmat(3))=a(jmia+nmat(3))+dvawi
            a(ial+nmat(3))=a(ial+nmat(3))-dvawi
            a(iau+nmat(3))=a(iau+nmat(3))+dvawkb
            a(jml+nmat(3))=a(jml+nmat(3))-dvawkb
            a(jmia+nmat(8))=a(jmia+nmat(8))+dvewi
            a(ial+nmat(8))=a(ial+nmat(8))-dvewi
            a(iau+nmat(8))=a(iau+nmat(8))+dvewkb
            a(jml+nmat(8))=a(jml+nmat(8))-dvewkb
c     
            a(jmia+nmat(4))=a(jmia+nmat(4))+dvayci
            a(ial+nmat(4))=a(ial+nmat(4))-dvayci
            a(iau+nmat(4))=a(iau+nmat(4))+dvayckb
            a(jml+nmat(4))=a(jml+nmat(4))-dvayckb
            a(jmia+nmat(9))=a(jmia+nmat(9))+dveyci
            a(ial+nmat(9))=a(ial+nmat(9))-dveyci
            a(iau+nmat(9))=a(iau+nmat(9))+dveyckb
            a(jml+nmat(9))=a(jml+nmat(9))-dveyckb
c     
            a(jmia+nmat(5))=a(jmia+nmat(5))+dvayai
            a(ial+nmat(5))=a(ial+nmat(5))-dvayai
            a(iau+nmat(5))=a(iau+nmat(5))+dvayakb
            a(jml+nmat(5))=a(jml+nmat(5))-dvayakb
            a(jmia+nmat(10))=a(jmia+nmat(10))+dveyai
            a(ial+nmat(10))=a(ial+nmat(10))-dveyai
            a(iau+nmat(10))=a(iau+nmat(10))+dveyakb
            a(jml+nmat(10))=a(jml+nmat(10))-dveyakb
c     
 65      continue
      end if
c     add heat conduction for water-only problem
      if(idof_co2.eq.1) then
         do jm=1,iq
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

            a(jmia+nmat(6))=a(jmia+nmat(6))+0.d0
            a(jmia+nmat(7))=a(jmia+nmat(7))-heatc
            a(ial+nmat(6))=a(ial+nmat(6))+0.d0
            a(ial+nmat(7))=a(ial+nmat(7))+heatc
            a(iau+nmat(6))=a(iau+nmat(6))+0.d0
            a(iau+nmat(7))=a(iau+nmat(7))+heatc
            a(jml+nmat(6))=a(jml+nmat(6))+0.d0
            a(jml+nmat(7))=a(jml+nmat(7))-heatc
         enddo
      endif

c     add accumulation terms
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i)
      a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*dmwf(i)+dqw(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*dmycf(i)+dqyc(i)
      a(jmia+nmat(5))=a(jmia+nmat(5))+sx1d*dmyaf(i)+dqya(i)
      a(jmia+nmat(6))=a(jmia+nmat(6))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(7))=a(jmia+nmat(7))+sx1d*deef(i)+deqh(i)
      a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*dewf(i)+deqw(i)
      a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*deycf(i)+deqyc(i)
      a(jmia+nmat(10))=a(jmia+nmat(10))+sx1d*deyaf(i)+deqya(i)
      r e t u r n
      e    n    d
