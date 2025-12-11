      subroutine  dpdpfh
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  To load the double porosity/double permeability solution into the
CD1  solution matrix for a heat, water, water vapor simulation.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dpdpfh.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:22   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:50   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Feb 15 09:50:46 1996   zvd
CD2 Corrected purpose and requirements.
CD2 
CD2    Rev 1.4   Wed Jan 10 11:10:16 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.3   Wed Jan 10 09:34:44 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   03/23/94 14:41:04   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:12   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.3.2 Heat- and mass-transfer equations
CD3 2.4.9 Double-porosity/double-permeability formulation
CD3 2.5.2  Solve nonlinear equation set at each time step
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

      use davidi
      use comji
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i, ial, iau, id, idg, idl, idum, isl, iz, jm, jml
      integer kb, kz, neighc, neq2, neqp1, nmat2avw
      real*8  aexy, aexyf, alen, al0, al1, area, athkb
      real*8  axkb, axy, axyd, axyf
      real*8  coef1, coefmx 
      real*8  delci, delckb, delei, delekb, deli, delkb 
      real*8  devci, devckb, devei, devekb, devi, devkb 
      real*8  devmci, devmei, devmpi 
      real*8  dfmlis, dfmlkbs, dfmvis, dfmvkbs 
      real*8  dilci, dilckb, dili, dilei, dilekb, dilkb, dilpi, dilpkb
      real*8  divci, divckb, divei, divekb, divi, divkb, divpi, divpkb
      real*8  dist01
      real*8  dlaei, dlaekb, dlapi, dlapkb
      real*8  dleei, dleekb, dlei, dlekb, dlepi, dlepkb
      real*8  dlpi, dlpkb, dpvti 
      real*8  dvaei, dvaekb, dvapi, dvapkb
      real*8  dveei, dveekb, dvei, dvekb, dvepi, dvepkb, dvpi, dvpkb
      real*8  fid, fid1, fmkb, fmli, fmlkb, fmvi, fmvkb, frac0, frac1
      real*8  heatc, enli, enlkb, envi, envkb
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb 
      real*8  swi, sx2c, sx2t, sx4d, sx4h
      real*8  thxi, thxkb, thyi, thzi, ti, tot 
      real*8  vexy, vexyf, vxy, vxyd, vxyf 
      parameter (coefmx=1000000000.0)

      neq2   =  neq+neq
      neqp1  =  neq+1
c
c zero out transfer terms
c
      do i=1,nmat(2)
         a(i+nmat(3))=0.
         a(i+nmat(4))=0.
         a(i+nmat(7))=0.
         a(i+nmat(8))=0.
         a(i+nmat(9))=0.
         a(i+nmat(10))=0.
         a(i+nmat(13))=0.
         a(i+nmat(14))=0.
      enddo
c     
c     compute the fracture-matrix transfer terms
c     
c     loop on nodes
      do id=1,neq
c     identify diagonal member
         idg =  nelmdg(id)
         idum =  idg-neqp1
c     identify matrix node
         i=id
         idl=id+neq
         tot=sx1(id)+sx1(idl)
         frac0=sx1(id)/tot
         frac1=sx1(idl)/tot
         alen=apuv1(id)
         al0=frac0*alen
         al1=frac1*alen
c     
c     determine contribution to node id
c     
c     
         dist01=0.5*(al0+al1)
         if(dist01.le.1.e-15) dist01=1.e-15
         area=tot/alen
         pvii=phi(i)
         phii=pvii-pcp(i)
         dpvti=dpcef(i)
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
         divpi=divp(i)
         divei=dive(i)
         thxi=thx(i)
         thyi=thy(i)
         thzi=thz(i)
         ti=t(i)
         swi=s(i)
c     
c     form constants for i>neq
c     
         iz=i
         coef1=-area/dist01
         if(coef1.lt.-coefmx) coef1=coefmx
         if(iupk.ne.0) then
            axkb=max(pny(idl),zero_t)*pnx(idl)
         else
          if(idpdp.eq.1) then
            axkb =  max( pnx(idl),zero_t)       
          else if(idpdp.eq.2) then
            axkb =  max( pny(idl),zero_t)       
          else if(idpdp.eq.3) then
            axkb =  max( pnz(idl),zero_t)       
          else
            axkb =  max( pnx(idl  ),pny(idl  ),
     *           pnz(idl  ),zero_t )
          endif
         endif
          if(idpdp.eq.1) then
            athkb =  max( thx(idl),zero_t)     
          else if(idpdp.eq.2) then
            athkb =  max( thy(idl),zero_t)     
          else if(idpdp.eq.3) then
            athkb =  max( thz(idl),zero_t)     
          else
            athkb =  max( thx(idl  ),thy(idl  ),
     *        thz(idl  ),zero_t )
          endif
c     
c     2-d geometry
c     
         if (icnl.ne.0) then
            radi=cord(iz,3)
         elseif (icnl.eq.0) then
            radi=1.0
         endif
         jm=1
         kb=idl
         kz=kb
         neighc=1
         perml(1)=axkb
         permv(1)=perml(1)
         radkb=radi
         sx2c=radkb*coef1
c     sx4d=-radkb*perml(igrav)*sx(iw,igrav)*grav
c     sx4h=-radkb*permv(igrav)*sx(iw,igrav)*grav
         thxkb=athkb
         sx2t=sx2c*thxkb
         pvikb=phi(kb)
         phikb=pvikb-pcp(kb)
         pxy=sx2c*perml(1)
         pvxy=sx2c*permv(1)
         pxyi=pxy*(phikb-phii)
         pxyh=pvxy*(pvikb-pvii)
         t1(neighc)=pxyi
         t2(neighc)=pxyh
         t3(neighc)=pxy
         t4(neighc)=pvxy
         t5(neighc)=sx2t
         t6(neighc)=0.0
         t7(neighc)=0.0
c     
c     liquid phase calculations
c     
         jm=1
         kb=idl
         kz=kb
         neighc=1
         pxyi=t1(neighc)
         sx4d=t6(neighc)
         axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=axyd
c     
c     determine upwind nodes and if liquid phase exists
c     
         isl=0
         jm=1
         kb=idl
         kz=kb
         neighc=1
         fid=.5
         axyd=t8(neighc)
         if(axyd.lt.0.0) fid=dnwgt
         if(axyd.gt.0.0) fid=upwgt
         if(t3(neighc).lt.0.0) t9(neighc)=fid
         if(t3(neighc).gt.0.0) t9(neighc)=fid
         if(swi+s(kb).ne.0.0) isl=1
c     
c     form equations
c     
         if(isl.ne.0) then
            jm=1
            kb=idl
            kz=kb
            neighc=1
            axyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            pxyi=t1(neighc)
            pxy=t3(neighc)
            sx4d=t6(neighc)
            dilkb=dil(kb)
            dilpkb=dilp(kb)
            dilekb=dile(kb)
            enlkb=enlf(kb)
            delkb=delf(kb)
            delekb=delef(kb)
            dilekb=dile(kb)
            dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pxy+0.5*sx4d*dglp(kb)*(cord(kz,igrav)-cord(iz,igrav))
            dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
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

            bp(i+nrhs(1))=bp(i+nrhs(1))+axy
            bp(i+nrhs(3))=bp(i+nrhs(3))-axy
            a(idum+nmat(1))=a(idum+nmat(1))+dlapi
            a(idum+nmat(2))=a(idum+nmat(2))+dlaei
            a(idum+nmat(9))=a(idum+nmat(9))-dlapi
            a(idum+nmat(10))=a(idum+nmat(10))-dlaei
            a(idum+nmat(3))=a(idum+nmat(3))+dlapkb
            a(idum+nmat(4))=a(idum+nmat(4))+dlaekb
            a(idum+nmat(11))=a(idum+nmat(11))-dlapkb
            a(idum+nmat(12))=a(idum+nmat(12))-dlaekb

            bp(i+nrhs(2))=bp(i+nrhs(2))+aexy
            bp(i+nrhs(4))=bp(i+nrhs(4))-aexy
            a(idum+nmat(5))=a(idum+nmat(5))+dlepi
            a(idum+nmat(6))=a(idum+nmat(6))+dleei
            a(idum+nmat(13))=a(idum+nmat(13))-dlepi
            a(idum+nmat(14))=a(idum+nmat(14))-dleei
            a(idum+nmat(7))=a(idum+nmat(7))+dlepkb
            a(idum+nmat(8))=a(idum+nmat(8))+dleekb
            a(idum+nmat(15))=a(idum+nmat(15))-dlepkb
            a(idum+nmat(16))=a(idum+nmat(16))-dleekb
         endif
c     
c     vapour phase calculations
c     
         jm=1
         kb=idl
         kz=kb
         neighc=1
         pxyh=t2(neighc)
         sx4h=t7(neighc)
         vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=vxyd
c     
c     determine upwind nodes and if vapour phase exists
c     
c     isl=0
         jm=1
         kb=idl
         kz=kb
         neighc=1
         fid=.5
         vxyd=t8(neighc)
         if(vxyd.lt.0.0) fid=dnwgt
         if(vxyd.gt.0.0) fid=upwgt
         if(t4(neighc).lt.0.0) t9(neighc)=fid
         if(t4(neighc).gt.0.0) t9(neighc)=fid
c     if(swi+s(kb).ne.2.0) isl=1
c     
c     form equations
c     
         if(isl.ne.0) then
            jm=1
            kb=idl
            kz=kb
            neighc=1
            fid=t9(neighc)
            fid1=1.0-fid
            pxyh=t2(neighc)
            pvxy=t4(neighc)
            sx4h=t7(neighc)
            vxyd=t8(neighc)
            divkb=div(kb)
            divpkb=divp(kb)
            divekb=dive(kb)
            envkb=envf(kb)
            devkb=devf(kb)
            devekb=devef(kb)
            divekb=dive(kb)
            dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-cord(iz,igrav))
            dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-cord(iz,igrav))
            dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
            dvekb=0.5*sx4h*dgve(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
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
            bp(i+nrhs(1))=bp(i+nrhs(1))+vxy
            bp(i+nrhs(3))=bp(i+nrhs(3))-vxy
            a(idum+nmat(1))=a(idum+nmat(1))+dvapi
            a(idum+nmat(2))=a(idum+nmat(2))+dvaei
            a(idum+nmat(9))=a(idum+nmat(9))-dvapi
            a(idum+nmat(10))=a(idum+nmat(10))-dvaei
            a(idum+nmat(3))=a(idum+nmat(3))+dvapkb
            a(idum+nmat(4))=a(idum+nmat(4))+dvaekb
            a(idum+nmat(11))=a(idum+nmat(11))-dvapkb
            a(idum+nmat(12))=a(idum+nmat(12))-dvaekb

            bp(i+nrhs(2))=bp(i+nrhs(2))+vexy
            bp(i+nrhs(4))=bp(i+nrhs(4))-vexy
            a(idum+nmat(5))=a(idum+nmat(5))+dvepi
            a(idum+nmat(6))=a(idum+nmat(6))+dveei
            a(idum+nmat(13))=a(idum+nmat(13))-dvepi
            a(idum+nmat(14))=a(idum+nmat(14))-dveei
            a(idum+nmat(7))=a(idum+nmat(7))+dvepkb
            a(idum+nmat(8))=a(idum+nmat(8))+dveekb
            a(idum+nmat(15))=a(idum+nmat(15))-dvepkb
            a(idum+nmat(16))=a(idum+nmat(16))-dveekb
         endif
c     
c     add heat conduction
c     
         jm=1
         kb=idl
         kz=kb
         neighc=1
         heatc=t5(neighc)
         iau=it11(jm)
         ial=it12(jm)
         jml=nelmdg(kb)-neqp1
         heatc=t5(neighc)
         bp(i+nrhs(2))=bp(i+nrhs(2))+heatc*(t(kb)-ti)
         bp(i+nrhs(4))=bp(i+nrhs(4))-heatc*(t(kb)-ti)
         a(idum+nmat(5))=a(idum+nmat(5))-heatc*dtpa(i)
         a(idum+nmat(6))=a(idum+nmat(6))-heatc*dtpae(i)
         a(idum+nmat(13))=a(idum+nmat(13))+heatc*dtpa(i)
         a(idum+nmat(14))=a(idum+nmat(14))+heatc*dtpae(i)
         a(idum+nmat(7))=a(idum+nmat(7))+heatc*dtpa(kb)
         a(idum+nmat(8))=a(idum+nmat(8))+heatc*dtpae(kb)
         a(idum+nmat(15))=a(idum+nmat(15))-heatc*dtpa(kb)
         a(idum+nmat(16))=a(idum+nmat(16))-heatc*dtpae(kb)
      enddo
      
      return
      end
