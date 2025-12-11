      subroutine  dpdpfa
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
CD1  solution matrix for an isothermal air-water simulation.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dpdpfa.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:32   pvcs
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
CD2    Rev 1.8   Thu Feb 15 09:50:46 1996   zvd
CD2 Corrected purpose and requirements.
CD2 
CD2    Rev 1.7   Wed Jan 10 11:08:50 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.6   Wed Jan 10 09:33:16 1996   hend
CD2 Revised Prolog
CD2 
CD2    Rev 1.5   Wed Jan 10 08:44:48 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   11/15/95 16:19:42   gaz
CD2 fixes so upwind perms except fracture matrix interaction
CD2 
CD2    Rev 1.3   08/07/95 11:09:04   awolf
CD2 Loads last part of a_axy with connection fluxes between
CD2 matrix and fracture for dpdp transport use later.
CD2 
CD2    Rev 1.2   03/23/94 14:41:02   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:22   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:10   pvcs
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

      use comflow
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

      integer i, id, idg, idl, idum, isl, iz, jm
      integer kb, kz, neighc, neq2, neqp1, nmat2avw
      real*8  alen, al0, al1, area
      real*8  axkb, axy, axyd, axyf
      real*8  coef1, coefmx 
      real*8  dfmlis, dfmlkbs, dfmvis, dfmvkbs 
      real*8  dilci, dilckb, dili, dilei, dilekb, dilkb, dilpi, dilpkb
      real*8  divci, divckb, divei, divekb, divi, divkb, divpi, divpkb
      real*8  dist01
      real*8  dlaei, dlaekb, dlapi, dlapkb
      real*8  dlei, dlekb, dlpi, dlpkb, dpvti 
      real*8  dvaei, dvaekb, dvapi, dvapkb
      real*8  dvei, dvekb, dvpi, dvpkb
      real*8  fid, fid1, fmkb, fmli, fmlkb, fmvi, fmvkb, frac0, frac1 
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb 
      real*8  swi, sx2c, sx2t, sx4d, sx4h
      real*8  ti, tot 
      real*8  vxy, vxyd, vxyf 
      parameter (coefmx=1000000000.0)

      neq2   =  neq+neq
      neqp1  =  neq+1
      nmat2avw = 2*(nelm(neqp1)-neqp1)
                                ! Guessed at value for fmkb
      fmkb = 1.0
c     
c     zero out transfer terms
c     
      do i=1,nmat(2)
         if (jswitch .ne. 0) then
            a(i+nmat(3))=0.
            a(i+nmat(4))=0.
            a(i+nmat(5))=0.
            a(i+nmat(6))=0.            
         else
            a(i+nmat(3))=0.
            a(i+nmat(4))=0.
            a(i+nmat(7))=0.
            a(i+nmat(8))=0.
            a(i+nmat(9))=0.
            a(i+nmat(10))=0.
            a(i+nmat(13))=0.
            a(i+nmat(14))=0.
         end if
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
         dist01=0.5*(al0+al1)
         if(dist01.le.1.e-15) dist01=1.e-15
         area=tot/alen
         pvii=phi(i)
         phii=pvii-pcp(i)
         dpvti=dpcef(i)
         kb=idl
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
         endif
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
     *              pnz(idl  ),zero_t )
            endif
         endif
c     
c     2-d geometry
c     
         if (icnl.ne.0)  then
            radi=cord(iz,3)
         elseif (icnl.eq.0)  then
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
         t5(neighc)=0.0
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
            dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pxy+0.5*sx4d*dglp(kb)*(cord(kz,igrav)-cord(iz,igrav))
            dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
            axyf=(fid*dilkb+fid1*dili)
            axy=axyd*axyf
            dlapi=dlpi*axyf+axyd*fid1*dilpi
            dlapkb=dlpkb*axyf+axyd*fid*dilpkb
            dlaei=dlei*axyf+axyd*fid1*dilei
            dlaekb=dlekb*axyf+axyd*fid*dilekb
            a_axy(nmat2avw+i) = axy
            
            bp(i+nrhs(1))=bp(i+nrhs(1))+axy
            bp(i+nrhs(3))=bp(i+nrhs(3))-axy
            a(idum+nmat(1))=a(idum+nmat(1))+dlapi
            a(idum+nmat(2))=a(idum+nmat(2))+dlaei
            a(idum+nmat(9-joff))=a(idum+nmat(9-joff))-dlapi
            a(idum+nmat(10-joff))=a(idum+nmat(10-joff))-dlaei
            a(idum+nmat(3))=a(idum+nmat(3))+dlapkb
            a(idum+nmat(4))=a(idum+nmat(4))+dlaekb
            a(idum+nmat(11-joff))=a(idum+nmat(11-joff))-dlapkb
            a(idum+nmat(12-joff))=a(idum+nmat(12-joff))-dlaekb
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
         
c     form equations
c     
         if(isl.ne.0 .and. jswitch .eq. 0) then
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
            dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
            dvekb=0.5*sx4h*dgve(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
            vxyf=(fid*divkb+fid1*divi)
            vxy=vxyd*vxyf
            dvapi=dvpi*vxyf+vxyd*fid1*divpi
            dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
            dvaei=dvei*vxyf+vxyd*fid1*divei
            dvaekb=dvekb*vxyf+vxyd*fid*divekb
            a_vxy(nmat2avw+i) = vxy
            
            bp(i+nrhs(2))=bp(i+nrhs(2))+vxy
            bp(i+nrhs(4))=bp(i+nrhs(4))-vxy
            a(idum+nmat(5))=a(idum+nmat(5))+dvapi
            a(idum+nmat(6))=a(idum+nmat(6))+dvaei
            a(idum+nmat(13))=a(idum+nmat(13))-dvapi
            a(idum+nmat(14))=a(idum+nmat(14))-dvaei
            a(idum+nmat(7))=a(idum+nmat(7))+dvapkb
            a(idum+nmat(8))=a(idum+nmat(8))+dvaekb
            a(idum+nmat(15))=a(idum+nmat(15))-dvapkb
            a(idum+nmat(16))=a(idum+nmat(16))-dvaekb
         endif
      enddo
      
      return
      end
