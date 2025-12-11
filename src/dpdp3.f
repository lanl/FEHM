      subroutine  dpdp3
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
CD1  solution matrix for a three degree of freedom system (air,water,
CD1  and heat).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dpdp3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:20   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Feb 15 09:50:46 1996   zvd
CD2 Corrected purpose and requirements.
CD2 
CD2    Rev 1.4   Wed Jan 10 11:07:20 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.3   Wed Jan 10 08:38:30 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   03/23/94 14:41:00   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:08   pvcs
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
c 18-nov-92
c got rid of initialization of na and nb
c must be passed through
c 25-july-92
c started programming implementation of dpdp 3dof
c 1-13-92
c finished initial programming

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comhi
      use comji
      use davidi
      use comdti
      use comai
      use comflow
      implicit none

      integer i, ial, iau, icd, id, idg, idl, idum, isl, iz, jm, jml 
      integer kb, kz, neighc, neq2, neqp1, nmat2avw
      real*8  acxy, acxyf, aexy, aexyf, alen, al0, al1, area, athkb 
      real*8  axkb, axy, axyd, axyf
      real*8  cnli, cnlkb, cnvi, cnvkb, coef1, coefmx 
      real*8  dclci, dclckb, dclei, dclekb, dcli, dclkb 
      real*8  dcvci, dcvckb, dcvei, dcvekb, dcvi, dcvkb 
      real*8  delci, delckb, delei, delekb, deli, delkb 
      real*8  devci, devckb, devei, devekb, devi, devkb 
      real*8  devmci, devmei, devmpi 
      real*8  dfmlis, dfmlkbs, dfmvis, dfmvkbs 
      real*8  dilci, dilckb, dili, dilei, dilekb, dilkb, dilpi, dilpkb
      real*8  divci, divckb, divei, divekb, divi, divkb, divpi, divpkb
      real*8  dist01
      real*8  dlaci, dlackb, dlaei, dlaekb, dlapi, dlapkb
      real*8  dlcci, dlcckb, dlci, dlcei, dlcekb, dlckb, dlcpi, dlcpkb 
      real*8  dleci, dleckb, dleei, dleekb, dlei, dlekb, dlepi, dlepkb
      real*8  dlpi, dlpkb, dpvti 
      real*8  dvaci, dvackb, dvaei, dvaekb, dvai, dvakb, dvapi, dvapkb
      real*8  dvam, dvame
      real*8  dvcci, dvcckb, dvcei, dvcekb, dvci, dvckb, dvcpi, dvcpkb 
      real*8  dveci, dveckb, dveei, dveekb, dvei, dvekb, dvepi, dvepkb
      real*8  dvmckb, dvmekb, dvmpkb, dvpi, dvpkb 
      real*8  heatc, enli, enlkb, envi, envkb, envmi, envmkb
      real*8  fid, fid1, fmkb, fmli, fmlkb, fmvi, fmvkb, frac0, frac1 
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb 
      real*8  swi, sx2c, sx2t, sx4d, sx4h
      real*8  thxi, thxkb, thyi, thzi, ti, tot 
      real*8  vcxy, vcxyf, vexy, vexyf, vxy, vxyd, vxyf 
      parameter (coefmx=1000000000.0)
      
      neq2   =  neq+neq
      neqp1  =  neq+1
! Guessed at values nmat2avw, icd, fmkb
      nmat2avw = 0
      icd = 0
      fmkb = 1.0
      if(irdof.eq.0) then
         do  i=1,nmat(2)
            a(i+nmat(4))=0.
            a(i+nmat(5))=0.
            a(i+nmat(6))=0.
            a(i+nmat(10))=0.
            a(i+nmat(11))=0.
            a(i+nmat(12))=0.
            a(i+nmat(16))=0.
            a(i+nmat(17))=0.
            a(i+nmat(18))=0.
            a(i+nmat(19))=0.
            a(i+nmat(20))=0.
            a(i+nmat(21))=0.
            a(i+nmat(25))=0.
            a(i+nmat(26))=0.
            a(i+nmat(27))=0.
            a(i+nmat(31))=0.
            a(i+nmat(32))=0.
            a(i+nmat(33))=0.
         enddo
      else
c     
c     zero out transfer terms
c     
c     note that transfer terms are neq long ********************
c     this means that if we normalize with the 6dof system we will
c     miss some fill on the transfer terms
c     
         do i=1,neq
            a(i+nmat(4))=0.
            a(i+nmat(5))=0.
            a(i+nmat(6))=0.
            a(i+nmat(10))=0.
            a(i+nmat(11))=0.
            a(i+nmat(12))=0.
            a(i+nmat(16))=0.
            a(i+nmat(17))=0.
            a(i+nmat(18))=0.
            a(i+nmat(19))=0.
            a(i+nmat(20))=0.
            a(i+nmat(21))=0.
            a(i+nmat(25))=0.
            a(i+nmat(26))=0.
            a(i+nmat(27))=0.
            a(i+nmat(31))=0.
            a(i+nmat(32))=0.
            a(i+nmat(33))=0.
         enddo
      endif
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
         kb=idl
         tot=sx1(id)+sx1(idl)
         frac0=sx1(id)/tot
         frac1=sx1(idl)/tot
         alen=apuv1(id)
         al0=frac0*alen
         al1=frac1*alen
c     
c     determine contribution to node id
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
         cnli=cnlf(i)
         dcli=dclf(i)
         dclei=dclef(i)
         dclci=dclcf(i)
         cnvi=cnvf(i)
         dcvi=dcvf(i)
         dcvei=dcvef(i)
         dcvci=dcvcf(i)
         delci=delcf(i)
         devci=devcf(i)
         if(irlpt(irlp(i)).eq.7) then
c     add fracture-matrix factors(Sandia)
            fmvi=fmvf(i)
            dfmvis=dfmvf(i)
            fmvkb=fmvf(kb)
            dfmvkbs=dfmvf(kb)
            fmli=fmlf(i)
            dfmlis=dfmlf(i)
            fmlkb=fmlf(kb)
            dfmlkbs=dfmlf(kb)
         else
            fmvi=1.0
            dfmvis=0.0
            fmvkb=1.0
            dfmvkbs=0.0
            fmli=1.0
            dfmlis=0.0
            fmlkb=1.0
            dfmlkbs=0.0
         endif
         if(iupk.ne.0) then
c     
c     even with iupk set pny contains the correct permeability
c     this will allow for harmonic weighting of the transfer perms
c     
            dili=dil(i)/pny(i)*fmli
            divi=div(i)/pny(i)*fmvi
            dilpi=dilp(i)/pny(i)*fmli
            dilei=dile(i)/pny(i)*fmli + dil(i)/pny(i)*dfmlis
            divpi=divp(i)/pny(i)*fmvi
            divei=dive(i)/pny(i)*fmvi + div(i)/pny(i)*dfmvis
            dilkb=dil(kb)/pny(kb)*fmlkb
            divkb=div(kb)/pny(kb)*fmvkb
            dilpkb=dilp(kb)/pny(kb)*fmkb
            dilekb=dile(kb)/pny(kb)*fmlkb+ dil(kb)/pny(kb)*dfmlkbs
            divpkb=divp(kb)/pny(kb)*fmvkb
            divekb=dive(kb)/pny(kb)*fmvkb+ div(kb)/pny(kb)*dfmvkbs
            dilci=dilc(i)/pny(i)*fmli
            divci=divc(i)/pny(i)*fmvi
            dilckb=dilc(kb)/pny(kb)*fmlkb
            divckb=divc(kb)/pny(kb)*fmvkb
         else
            dili=dil(i)*fmli
            divi=div(i)*fmvi
            dilpi=dilp(i)*fmli
            dilei=dile(i)*fmli + dil(i)*dfmlis
            divpi=divp(i)*fmvi
            divei=dive(i)*fmvi + div(i)*dfmvis
            dilkb=dil(kb)*fmlkb
            divkb=div(kb)*fmvkb
            dilpkb=dilp(kb)*fmkb
            dilekb=dile(kb)*fmlkb+ dil(kb)*dfmlkbs
            divpkb=divp(kb)*fmvkb
            divekb=dive(kb)*fmvkb+ div(kb)*dfmvkbs
            dilci=dilc(i)*fmli
            divci=divc(i)*fmvi
            dilckb=dilc(kb)*fmlkb
            divckb=divc(kb)*fmvkb
         endif
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
         if(icnl.ne.0) then
            radi=cord(iz,3)
         else if(icnl.eq.0) then
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
         t10(neighc)=sx2c
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
            enlkb=enlf(kb)
            cnlkb=cnlf(kb)
            dilpkb=dilp(kb)
            dilekb=dile(kb)
            dilckb=dilc(kb)
            delkb=delf(kb)
            delekb=delef(kb)
            delckb=delcf(kb)
            dclkb=dclf(kb)
            dclekb=dclef(kb)
            dclckb=dclcf(kb)
            dilkb=dil(kb)
            dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pxy+0.5*sx4d*dglp(kb)*(cord(kz,igrav)-cord(iz,igrav))
            dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
            dlci=0.5*sx4d*dglc(i)*(cord(kz,igrav)-cord(iz,igrav))
            dlckb=0.5*sx4d*dglc(kb)*(cord(kz,igrav)-cord(iz,igrav))
            aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            axyf=(fid*dilkb+fid1*dili)
            acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
            aexy=axyd*aexyf
            axy=axyd*axyf
            acxy=axyd*acxyf
            dlapi=dlpi*axyf+axyd*fid1*dilpi
            dlapkb=dlpkb*axyf+axyd*fid*dilpkb
            dlaei=dlei*axyf+axyd*fid1*dilei
            dlaekb=dlekb*axyf+axyd*fid*dilekb
            dlaci=dlci*axyf+axyd*fid1*dilci
            dlackb=dlckb*axyf+axyd*fid*dilckb
            dlepi=dlpi*aexyf+axyd*fid1*(deli*dili+enli*dilpi)
            dlepkb=dlpkb*aexyf+axyd*fid*(delkb*dilkb+enlkb*dilpkb)
            dleei=dlei*aexyf+axyd*fid1*(delei*dili+enli*dilei)
            dleekb=dlekb*aexyf+axyd*fid*(delekb*dilkb+enlkb*dilekb)
            dleci=dlci*aexyf+axyd*fid1*(delci*dili+enli*dilci)
            dleckb=dlckb*aexyf+axyd*fid*(delckb*dilkb+enlkb*dilckb)
            dlcpi=dlpi*acxyf+axyd*fid1*(dcli*dili+cnli*dilpi)
            dlcpkb=dlpkb*acxyf+axyd*fid*(dclkb*dilkb+cnlkb*dilpkb)
            dlcei=dlei*acxyf+axyd*fid1*(dclei*dili+cnli*dilei)
            dlcekb=dlekb*acxyf+axyd*fid*(dclekb*dilkb+cnlkb*dilekb)
            dlcci=dlci*acxyf+axyd*fid1*(dclci*dili+cnli*dilci)
            dlcckb=dlckb*acxyf+axyd*fid*(dclckb*dilkb+cnlkb*dilckb)
            a_axy(nmat2avw+i) = axy
            
            bp(i+nrhs(1))=bp(i+nrhs(1))+axy
            bp(i+nrhs(4))=bp(i+nrhs(4))-axy
            a(idum+nmat(1))=a(idum+nmat(1))+dlapi
            a(idum+nmat(2))=a(idum+nmat(2))+dlaei
            a(idum+nmat(3))=a(idum+nmat(3))+dlaci
            a(idum+nmat(19))=a(idum+nmat(19))-dlapi
            a(idum+nmat(20))=a(idum+nmat(20))-dlaei
            a(idum+nmat(21))=a(idum+nmat(21))-dlaci
            a(idum+nmat(4))=a(idum+nmat(4))+dlapkb
            a(idum+nmat(5))=a(idum+nmat(5))+dlaekb
            a(idum+nmat(6))=a(idum+nmat(6))+dlackb
            a(idum+nmat(22))=a(idum+nmat(22))-dlapkb
            a(idum+nmat(23))=a(idum+nmat(23))-dlaekb
            a(idum+nmat(24))=a(idum+nmat(24))-dlackb
            
            bp(i+nrhs(2))=bp(i+nrhs(2))+aexy
            bp(i+nrhs(5))=bp(i+nrhs(5))-aexy
            a(idum+nmat(7))=a(idum+nmat(7))+dlepi
            a(idum+nmat(8))=a(idum+nmat(8))+dleei
            a(idum+nmat(9))=a(idum+nmat(9))+dleci
            a(idum+nmat(25))=a(idum+nmat(25))-dlepi
            a(idum+nmat(26))=a(idum+nmat(26))-dleei
            a(idum+nmat(27))=a(idum+nmat(27))-dleci
            a(idum+nmat(10))=a(idum+nmat(10))+dlepkb
            a(idum+nmat(11))=a(idum+nmat(11))+dleekb
            a(idum+nmat(12))=a(idum+nmat(12))+dleckb
            a(idum+nmat(28))=a(idum+nmat(28))-dlepkb
            a(idum+nmat(29))=a(idum+nmat(29))-dleekb
            a(idum+nmat(30))=a(idum+nmat(30))-dleckb
            
            bp(i+nrhs(3))=bp(i+nrhs(3))+acxy
            bp(i+nrhs(6))=bp(i+nrhs(6))-acxy
            a(idum+nmat(13))=a(idum+nmat(13))+dlcpi
            a(idum+nmat(14))=a(idum+nmat(14))+dlcei
            a(idum+nmat(15))=a(idum+nmat(15))+dlcci
            a(idum+nmat(31))=a(idum+nmat(31))-dlcpi
            a(idum+nmat(32))=a(idum+nmat(32))-dlcei
            a(idum+nmat(33))=a(idum+nmat(33))-dlcci
            a(idum+nmat(16))=a(idum+nmat(16))+dlcpkb
            a(idum+nmat(17))=a(idum+nmat(17))+dlcekb
            a(idum+nmat(18))=a(idum+nmat(18))+dlcckb
            a(idum+nmat(34))=a(idum+nmat(34))-dlcpkb
            a(idum+nmat(35))=a(idum+nmat(35))-dlcekb
            a(idum+nmat(36))=a(idum+nmat(36))-dlcckb
            
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
         isl=1
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
            
            envkb=envf(kb)
            cnvkb=cnvf(kb)
            divpkb=divp(kb)
            divckb=divc(kb)
            divekb=dive(kb)
            devkb=devf(kb)
            devekb=devef(kb)
            devckb=devcf(kb)
            dcvkb=dcvf(kb)
            dcvekb=dcvef(kb)
            dcvckb=dcvcf(kb)
            divkb=div(kb)
            
            dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
            dvekb=0.5*sx4h*dgve(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
            dvci=0.5*sx4h*dgvc(i)*(cord(kz,igrav)-cord(iz,igrav))
            dvckb=0.5*sx4h*dgvc(kb)*(cord(kz,igrav)-cord(iz,igrav))
            
            vexyf=(fid*divkb*envkb+fid1*divi*envi)
            vxyf=(fid*divkb+fid1*divi)
            vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
            vexy=vxyd*vexyf
            vxy=vxyd*vxyf
            vcxy=vxyd*vcxyf
            dvapi=dvpi*vxyf+vxyd*fid1*divpi
            dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
            dvaei=dvei*vxyf+vxyd*fid1*divei
            dvaekb=dvekb*vxyf+vxyd*fid*divekb
            dvaci=dvci*vxyf+vxyd*fid1*divci
            dvackb=dvckb*vxyf+vxyd*fid*divckb
            dvepi=dvpi*vexyf+vxyd*fid1*(devi*divi+envi*divpi)
            dvepkb=dvpkb*vexyf+vxyd*fid*(devkb*divkb+envkb*divpkb)
            dveei=dvei*vexyf+vxyd*fid1*(devei*divi+envi*divei)
            dveekb=dvekb*vexyf+vxyd*fid*(devekb*divkb+envkb*divekb)
            dveci=dvci*vexyf+vxyd*fid1*(devci*divi+envi*divci)
            dveckb=dvckb*vexyf+vxyd*fid*(devckb*divkb+envkb*divckb)
            dvcpi=dvpi*vcxyf+vxyd*fid1*(dcvi*divi+cnvi*divpi)
            dvcpkb=dvpkb*vcxyf+vxyd*fid*(dcvkb*divkb+cnvkb*divpkb)
            dvcei=dvei*vcxyf+vxyd*fid1*(dcvei*divi+cnvi*divei)
            dvcekb=dvekb*vcxyf+vxyd*fid*(dcvekb*divkb+cnvkb*divekb)
            dvcci=dvci*vcxyf+vxyd*fid1*(dcvci*divi+cnvi*divci)
            dvcckb=dvckb*vcxyf+vxyd*fid*(dcvckb*divkb+cnvkb*divckb)
            a_vxy(nmat2avw+i) = vxy
            
            bp(i+nrhs(1))=bp(i+nrhs(1))+vxy
            bp(i+nrhs(4))=bp(i+nrhs(4))-vxy
            a(idum+nmat(1))=a(idum+nmat(1))+dvapi
            a(idum+nmat(2))=a(idum+nmat(2))+dvaei
            a(idum+nmat(3))=a(idum+nmat(3))+dvaci
            a(idum+nmat(19))=a(idum+nmat(19))-dvapi
            a(idum+nmat(20))=a(idum+nmat(20))-dvaei
            a(idum+nmat(21))=a(idum+nmat(21))-dvaci
            a(idum+nmat(4))=a(idum+nmat(4))+dvapkb
            a(idum+nmat(5))=a(idum+nmat(5))+dvaekb
            a(idum+nmat(6))=a(idum+nmat(6))+dvackb
            a(idum+nmat(22))=a(idum+nmat(22))-dvapkb
            a(idum+nmat(23))=a(idum+nmat(23))-dvaekb
            a(idum+nmat(24))=a(idum+nmat(24))-dvackb
            
            bp(i+nrhs(2))=bp(i+nrhs(2))+vexy
            bp(i+nrhs(5))=bp(i+nrhs(5))-vexy
            a(idum+nmat(7))=a(idum+nmat(7))+dvepi
            a(idum+nmat(8))=a(idum+nmat(8))+dveei
            a(idum+nmat(9))=a(idum+nmat(9))+dveci
            a(idum+nmat(25))=a(idum+nmat(25))-dvepi
            a(idum+nmat(26))=a(idum+nmat(26))-dveei
            a(idum+nmat(27))=a(idum+nmat(27))-dveci
            a(idum+nmat(10))=a(idum+nmat(10))+dvepkb
            a(idum+nmat(11))=a(idum+nmat(11))+dveekb
            a(idum+nmat(12))=a(idum+nmat(12))+dveckb
            a(idum+nmat(28))=a(idum+nmat(28))-dvepkb
            a(idum+nmat(29))=a(idum+nmat(29))-dveekb
            a(idum+nmat(30))=a(idum+nmat(30))-dveckb
            
            bp(i+nrhs(3))=bp(i+nrhs(3))+vcxy
            bp(i+nrhs(6))=bp(i+nrhs(6))-vcxy
            a(idum+nmat(13))=a(idum+nmat(13))+dvcpi
            a(idum+nmat(14))=a(idum+nmat(14))+dvcei
            a(idum+nmat(15))=a(idum+nmat(15))+dvcci
            a(idum+nmat(31))=a(idum+nmat(31))-dvcpi
            a(idum+nmat(32))=a(idum+nmat(32))-dvcei
            a(idum+nmat(33))=a(idum+nmat(33))-dvcci
            a(idum+nmat(16))=a(idum+nmat(16))+dvcpkb
            a(idum+nmat(17))=a(idum+nmat(17))+dvcekb
            a(idum+nmat(18))=a(idum+nmat(18))+dvcckb
            a(idum+nmat(34))=a(idum+nmat(34))-dvcpkb
            a(idum+nmat(35))=a(idum+nmat(35))-dvcekb
            a(idum+nmat(36))=a(idum+nmat(36))-dvcckb
            if(iadif.ne.0) then
c     
c     add air-water vapor diffusion
c     
c     
c     determine upwind nodes and if vapour phase exists
c     
               jm=1
               kb=idl    
               neighc=1
               fid=.5
               vxyd=t10(neighc)*(cnvf(kb)-cnvi)
               if ( vxyd.lt.0.0) fid=dnwgt
               if ( vxyd.gt.0.0) fid=upwgt
c              t9(neighc)=0.5
c use matrix values
               t9(neighc)=1.0
               dvai=dva(i) 
               envmi=enva(i)
               devmpi=denvap(i)
               devmei=denvae(i)
               devmci=denvac(i)
               devmpi=devmpi*dvai+envmi*ddvap(i)
               devmei=devmei*dvai+envmi*ddvae(i)
               devmci=devmci*dvai+envmi*ddvac(i)
               jm=1
               kb=idl    
               kz=kb-icd
               neighc=1
               heatc=t10(neighc)
               fid=t9(neighc)
               fid1=1.0-fid
c     energy equation
               envmkb=enva(kb)
               dvmpkb=denvap(kb)
               dvmekb=denvae(kb)
               dvmckb=denvac(kb)
               dvakb=dva(kb)
               dvmckb=dvmckb*dvakb+envmkb*ddvac(kb)
               dvmekb=dvmekb*dvakb+envmkb*ddvae(kb)
               dvmpkb=dvmpkb*dvakb+envmkb*ddvap(kb)
               dvame=(fid*dvakb*envmkb+fid1*dvai*envmi)
               bp(i+nrhs(2))=bp(i+nrhs(2))+heatc*dvame*
     &              (cnvf(kb)-cnvi)
               bp(i+nrhs(5))=bp(i+nrhs(5))-heatc*dvame*
     &              (cnvf(kb)-cnvi)
               a(idum+nmat(7))=a(idum+nmat(7))+heatc*(devmpi*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dcvf(i))
               a(idum+nmat(8))=a(idum+nmat(8))+heatc*(devmei*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dcvef(i))
               a(idum+nmat(9))=a(idum+nmat(9))+heatc*(devmci*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dcvcf(i))
               a(idum+nmat(25))=a(idum+nmat(25))-heatc*(devmpi*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dcvf(i))
               a(idum+nmat(26))=a(idum+nmat(26))-heatc*(devmei*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dcvef(i))
               a(idum+nmat(27))=a(idum+nmat(27))-heatc*(devmci*fid1
     &              *(cnvf(kb)-cnvi)-dvame*dvam*dcvef(i))
               a(idum+nmat(10))=a(idum+nmat(10))+heatc*(dvmpkb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvf(kb))
               a(idum+nmat(11))=a(idum+nmat(11))+heatc*(dvmekb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvef(kb))
               a(idum+nmat(12))=a(idum+nmat(12))+heatc*(dvmckb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvcf(kb))
               a(idum+nmat(28))=a(idum+nmat(28))-heatc*(dvmpkb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvf(kb))
               a(idum+nmat(29))=a(idum+nmat(29))-heatc*(dvmekb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvef(kb))
               a(idum+nmat(30))=a(idum+nmat(30))-heatc*(dvmckb*fid
     &              *(cnvf(kb)-cnvi)+dvame*dcvcf(kb))
c     air equation
               dvakb=dva(kb)
               dvam=(fid*dvakb+fid1*dvai)
               bp(i+nrhs(3))=bp(i+nrhs(3))+heatc*dvam*(cnvf(kb)-cnvi)
               bp(i+nrhs(6))=bp(i+nrhs(6))-heatc*dvam*(cnvf(kb)-cnvi)
               a(idum+nmat(13))=a(idum+nmat(13))+heatc*(ddvap(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvf(i))
               a(idum+nmat(14))=a(idum+nmat(14))+heatc*(ddvae(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat(15))=a(idum+nmat(15))+heatc*(ddvac(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvcf(i))
               a(idum+nmat(31))=a(idum+nmat(31))-heatc*(ddvap(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvf(i))
               a(idum+nmat(32))=a(idum+nmat(32))-heatc*(ddvae(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat(33))=a(idum+nmat(33))-heatc*(ddvac(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat(16))=a(idum+nmat(16))+heatc*(ddvap(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvf(kb))
               a(idum+nmat(17))=a(idum+nmat(17))+heatc*(ddvae(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvef(kb))
               a(idum+nmat(18))=a(idum+nmat(18))+heatc*(ddvac(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvcf(kb))
               a(idum+nmat(34))=a(idum+nmat(34))-heatc*(ddvap(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvf(kb))
               a(idum+nmat(35))=a(idum+nmat(35))-heatc*(ddvae(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvef(kb))
               a(idum+nmat(36))=a(idum+nmat(36))-heatc*(ddvac(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvcf(kb))
c    water equation
               dvakb=dva(kb)
               dvam=(fid*dvakb+fid1*dvai)
               bp(i+nrhs(1))=bp(i+nrhs(1))-heatc*dvam*(cnvf(kb)-cnvi)
               bp(i+nrhs(4))=bp(i+nrhs(4))+heatc*dvam*(cnvf(kb)-cnvi)
               a(idum+nmat( 1))=a(idum+nmat( 1))-heatc*(ddvap(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvf(i))
               a(idum+nmat( 2))=a(idum+nmat( 2))-heatc*(ddvae(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat( 3))=a(idum+nmat( 3))-heatc*(ddvac(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvcf(i))
               a(idum+nmat(19))=a(idum+nmat(19))+heatc*(ddvap(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvf(i))
               a(idum+nmat(20))=a(idum+nmat(20))+heatc*(ddvae(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat(21))=a(idum+nmat(21))+heatc*(ddvac(i)*fid1
     &              *(cnvf(kb)-cnvi)-dvam*dcvef(i))
               a(idum+nmat( 4))=a(idum+nmat( 4))-heatc*(ddvap(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvf(kb))
               a(idum+nmat( 5))=a(idum+nmat( 5))-heatc*(ddvae(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvef(kb))
               a(idum+nmat( 6))=a(idum+nmat( 6))-heatc*(ddvac(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvcf(kb))
               a(idum+nmat(22))=a(idum+nmat(22))+heatc*(ddvap(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvf(kb))
               a(idum+nmat(23))=a(idum+nmat(23))+heatc*(ddvae(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvef(kb))
               a(idum+nmat(24))=a(idum+nmat(24))+heatc*(ddvac(kb)*fid
     &              *(cnvf(kb)-cnvi)+dvam*dcvcf(kb))
            endif
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
         bp(i+nrhs(5))=bp(i+nrhs(5))-heatc*(t(kb)-ti)
         a(idum+nmat(7))=a(idum+nmat(7))-heatc*dtpa(i)
         a(idum+nmat(8))=a(idum+nmat(8))-heatc*dtpae(i)
         a(idum+nmat(9))=a(idum+nmat(9))-heatc*dtpac(i)
         a(idum+nmat(25))=a(idum+nmat(25))+heatc*dtpa(i)
         a(idum+nmat(26))=a(idum+nmat(26))+heatc*dtpae(i)
         a(idum+nmat(27))=a(idum+nmat(27))+heatc*dtpac(i)
         a(idum+nmat(10))=a(idum+nmat(10))+heatc*dtpa(kb)
         a(idum+nmat(11))=a(idum+nmat(11))+heatc*dtpae(kb)
         a(idum+nmat(12))=a(idum+nmat(12))+heatc*dtpac(kb)
         a(idum+nmat(28))=a(idum+nmat(28))-heatc*dtpa(kb)
         a(idum+nmat(29))=a(idum+nmat(29))-heatc*dtpae(kb)
         a(idum+nmat(30))=a(idum+nmat(30))-heatc*dtpac(kb)
      enddo
      
      return
      end
