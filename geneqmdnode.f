      subroutine geneqmdnode
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To generate equations fo multiply defined nodes                  
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 04-04-96     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/geneqmdnode.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Jun 06 15:10:48 1996   robinson
CD2 Fixed implementation for dpdp
CD2 
CD2    Rev 1.2   Wed May 22 08:55:58 1996   gaz
CD2 corrections and additions
CD2 
C**********************************************************************
CD3
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 2.3.2 Heat- and mass-transfer equations
CD3 2.3.3 Noncondensible gas flow equations
CD3
C**********************************************************************
CD4
CD4 SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

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

      integer i,i1,i2,icd,ii,jj,ii1,ii2,idg,iq,jmi,jml,jmia,jm,neqp1
      integer ij,ij1,ij2,iz,kb,neighc,iau,ial,kz,nmatavw,iz4m1
      real*8 sx1d,pvii,phii,dlpi,dlpkb,dvpi,dvpkb,dlei,dvei,dlekb
      real*8 dvekb,dlapi,dvapi,dlapkb,dvapkb,pvikb,phikb,axyd,axy,axyf
      real*8 pvxy,pxy,pxyh,pxyi,vxyd,vxy,vxyf,sx4d,sx4h,dpvti
      real*8 dlaei,dlaekb,dvaei,dvaekb,sx_md,fid,fid1,dili,dilpi,dilei
      real*8 dilkb,dilpkb,dilekb,divi,divpi,divei,divkb,divpkb,divekb
      real*8 heatc,ti,cnvi,dcvi,dcvei,dcvci,cnvkb,dcvkb,dcvekb,dcvckb
      real*8 dlci,dlckb,dlaci,dlackb,dvci,dvckb,dvaci,dvackb 
      logical bit

      neqp1=neq+1
      
c     loop on nodes
c     
      do i=1,n
c     
c     form constants for i>neq
c     
         if(i.gt.neq.and.idualp.eq.0) then
            icd=neq
         else
            icd=0
         endif
c     
c     loop on mdnodes
c     
         do jj=1,abs(mdnodes(i-icd)) 
            ii=mdnode(i-icd,jj)
            if(ii.gt.i-icd) then
c     
c     storage locations for fluxes
c     
               if(i.gt.neq) then
                  nmatavw=ldna
               else
                  nmatavw=0
               endif
c     
c     node i value of parameters
c     
               sx1d=sx1(i)
               if(i.gt.neq.and.idpdp.ne.0) then
c     means double permeability
                  if(ico2.le.0) then
                     nmat(1)=ldna*(11-1)
                     nmat(2)=ldna*(12-1)
                     nmat(3)=ldna*(15-1)
                     nmat(4)=ldna*(16-1)
                     nrhs(1)=(3-1)*neq
                     nrhs(2)=(4-1)*neq
                  else if(ico2.gt.0) then 
                     nmat(1)=ldna*(22-1)
                     nmat(2)=ldna*(23-1)
                     nmat(3)=ldna*(24-1)
                     nmat(4)=ldna*(28-1)
                     nmat(5)=ldna*(29-1)
                     nmat(6)=ldna*(30-1)
                     nmat(7)=ldna*(34-1)
                     nmat(8)=ldna*(35-1)
                     nmat(9)=ldna*(36-1)
                     nrhs(1)=(4-1)*neq
                     nrhs(2)=(5-1)*neq
                     nrhs(3)=(6-1)*neq
                  endif
               else if(idpdp.ne.0) then
                  if(ico2.le.0) then
                     nmat(1)=0
                     nmat(2)=ldna
                     nmat(3)=ldna*(5-1)
                     nmat(4)=ldna*(6-1)
                     nrhs(1)=0    
                     nrhs(2)=neq
                  else if(ico2.gt.0) then 
                     nmat(1)=0
                     nmat(2)=ldna
                     nmat(3)=ldna*2
                     nmat(4)=ldna*(7-1)
                     nmat(5)=ldna*(8-1)
                     nmat(6)=ldna*(9-1)
                     nmat(7)=ldna*(13-1)
                     nmat(8)=ldna*(14-1)
                     nmat(9)=ldna*(15-1)
                     nrhs(1)=0
                     nrhs(2)=neq
                     nrhs(3)=neq+neq
                  endif
               endif
               iz=i-icd
c     
c     calculate flow coefficient
c     
               sx_md=0.0
               i1=nelm(iz)+1
               i2=nelm(iz+1)
               do jm=i1,i2
                  iw=istrw(jm-neqp1)
                  sx_md=sx_md+(sx(iw,1)+sx(iw,2)+sx(iw,3))
               enddo
               i1=nelm(ii)+1
               i2=nelm(ii+1)
               do jm=i1,i2
                  iw=istrw(jm-neqp1)
                  sx_md=sx_md+(sx(iw,1)+sx(iw,2)+sx(iw,3))
               enddo
               sx_md=sx_md*sx_mult
               iz4m1 = 4*(iz-1)+1
               ii1=nelm(iz)+1
               ii2=nelm(iz+1)
               idg=nelmdg(iz)-ii1+1
               neqp1=neq+1
               jmi=nelmdg(iz)
               jmia=jmi-neqp1
               
c     Take care of source/sink term by putting it into empty slot (from 
c     same node to same node). 
c     If this is an isothermal air-water simulation
               a_axy(jmia+nmatavw)=sk(i)
c     Take care of source/sink term 
c     If this is an isothermal air-water simulation
               a_vxy(jmia+nmatavw)=qh(i)
               
               iq=0
               do  jm=jmi+1,ii2
                  if(nelm(jm).eq.ii) then
                     iq=iq+1
                     it8(iq)=nelm(jm)+icd
                     it9(iq)=jm-ii1+1
                     it10(iq)=istrw(jm-neqp1)
                     it11(iq)=jm-neqp1
                     ij1=nelm(nelm(jm))+1
                     ij2=nelmdg(nelm(jm))-1
                     do  ij=ij1,ij2
                        if(nelm(ij).eq.iz) then
                           it12(iq)=ij-neqp1
                        endif
                     enddo
                  endif
               enddo
c     
c     if block on degree of freedom
c     
               if(idof.gt.1) then
c     non heat conduction only
c     
c     ifblock on equation type
c     
                  pvii=phi(i)
                  phii=pvii-pcp(i)
                  dpvti=dpcef(i)
                  if(ico2.lt.0) then

c     3-d and 2-d geometry
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     pvikb=phi(kb)
                     phikb=pvikb-pcp(kb)
                     pxy=sx_md*
     &               max(pnx(i),pny(i),pnz(i),
     &               pnx(kb),pny(kb),pnz(kb))
                     pvxy=pxy    
                     pxyi=pxy*(phikb-phii)
                     pxyh=pvxy*(pvikb-pvii)
                     t1(neighc)=pxyi
                     t2(neighc)=pxyh
                     t3(neighc)=pxy
                     t4(neighc)=pvxy
                     t5(neighc)=sx_md*max(thx(i),thy(i),thz(i),
     &               thx(kb),thy(kb),thz(kb))
                     t6(neighc)=-grav*t3(neighc)
                     t7(neighc)=-grav*t4(neighc)
c     
c     liquid phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyi=t1(neighc)
                     sx4d=t6(neighc)
                     axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=axyd
c     
c     determine upwind nodes and if liquid phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        axyd=t8(neighc)
                        if(axyd.lt.0.0) fid=dnwgt
                        if(axyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_l(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_l(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     dili=s(i)
                     dilpi=0.0
                     dilei=1.0
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyi=t1(neighc)
                     pxy=t3(neighc)
                     sx4d=t6(neighc)
                     axyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     dilkb=s(kb)
                     dilpkb=0.0
                     dilekb=1.0
                     dlpi=-pxy+0.5*sx4d*dglp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlpkb=pxy+0.5*sx4d*dglp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlei=pxy*dpvti+0.5*sx4d*dgle(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))

                     axyf=(fid*dilkb+fid1*dili)
                     axy=axyd*axyf
                     dlapi=dlpi*axyf+axyd*fid1*dilpi
                     dlapkb=dlpkb*axyf+axyd*fid*dilpkb
                     dlaei=dlei*axyf+axyd*fid1*dilei
                     dlaekb=dlekb*axyf+axyd*fid*dilekb

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
c     vapour phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyh=t2(neighc)
                     sx4h=t7(neighc)
                     vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=vxyd
c     
c     determine upwind nodes and if vapour phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        vxyd=t8(neighc)
                        if(vxyd.lt.0.0) fid=dnwgt
                        if(vxyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_v(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_v(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     divi=1.0-s(i)
                     divpi=0.0
                     divei=-1.0
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyh=t2(neighc)
                     pvxy=t4(neighc)
                     sx4h=t7(neighc)
                     vxyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     divkb=1.0-s(kb)
                     divpkb=0.0
                     divekb=-1.0
                     dvpi=-pvxy+0.5*sx4h*dgvp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvpkb=pvxy+0.5*sx4h*dgvp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-
     &                    cord(iz,igrav))
                     dvekb=0.5*sx4h*dgve(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     vxyf=(fid*divkb+fid1*divi)
                     vxy=vxyd*vxyf
                     dvapi=dvpi*vxyf+vxyd*fid1*divpi
                     dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
                     dvaei=dvei*vxyf+vxyd*fid1*divei
                     dvaekb=dvekb*vxyf+vxyd*fid*divekb

                     a_vxy(iau+nmatavw)=vxy
                     a_vxy(ial+nmatavw)=-vxy
                     
                     bp(iz+nrhs(2))=bp(iz+nrhs(2))+vxy
                     bp(kz+nrhs(2))=bp(kz+nrhs(2))-vxy
                     a(jmia+nmat(3))=a(jmia+nmat(3))+dvapi
                     a(jmia+nmat(4))=a(jmia+nmat(4))+dvaei
                     a(ial+nmat(3))=a(ial+nmat(3))-dvapi
                     a(ial+nmat(4))=a(ial+nmat(4))-dvaei
                     a(iau+nmat(3))=a(iau+nmat(3))+dvapkb
                     a(iau+nmat(4))=a(iau+nmat(4))+dvaekb
                     a(jml+nmat(3))=a(jml+nmat(3))-dvapkb
                     a(jml+nmat(4))=a(jml+nmat(4))-dvaekb
                  else if(ico2.eq.0.and.idof.ne.1) then
c     
c     ifblock for heat/water (no air) equation type
c     
c     3-d and 2-d geometry
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     pvikb=phi(kb)
                     phikb=pvikb-pcp(kb)
                     pxy=sx_md*
     &               max(pnx(i),pny(i),pnz(i),
     &               pnx(kb),pny(kb),pnz(kb))
                     pvxy=pxy  
                     pxyi=pxy*(phikb-phii)
                     pxyh=pvxy*(pvikb-pvii)
                     t1(neighc)=pxyi
                     t2(neighc)=pxyh
                     t3(neighc)=pxy
                     t4(neighc)=pvxy
                     t5(neighc)=sx_md*max(thx(i),thy(i),thz(i),
     &               thx(kb),thy(kb),thz(kb))
                     t6(neighc)=-grav*t3(neighc)
                     t7(neighc)=-grav*t4(neighc)
c     
c     liquid phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyi=t1(neighc)
                     sx4d=t6(neighc)
                     axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=axyd
c     
c     determine upwind nodes and if liquid phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        axyd=t8(neighc)
                        if(axyd.lt.0.0) fid=dnwgt
                        if(axyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_l(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_l(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     dili=s(i)
                     dilpi=0.0
                     if(ieos(i).ne.2) then
                        dilei=0.0
                     else
                        dilei=1.0
                     endif
c     
c     enli=0.01*t(i)
c     deli=0.01*dtpa(i)
c     delei=0.01*dtpae(i)
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyi=t1(neighc)
                     pxy=t3(neighc)
                     sx4d=t6(neighc)
                     axyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     dilkb=s(kb)
                     dilpkb=0.0
                     if(ieos(kb).ne.2) then
                        dilekb=0.0
                     else
                        dilekb=1.0
                     endif
                     dlpi=-pxy+0.5*sx4d*dglp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlpkb=pxy+0.5*sx4d*dglp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlei=pxy*dpvti+0.5*sx4d*dgle(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     axyf=(fid*dilkb+fid1*dili)
                     axy=axyd*axyf
                     dlapi=dlpi*axyf+axyd*fid1*dilpi
                     dlapkb=dlpkb*axyf+axyd*fid*dilpkb
                     dlaei=dlei*axyf+axyd*fid1*dilei
                     dlaekb=dlekb*axyf+axyd*fid*dilekb

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
c     
c     vapour phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyh=t2(neighc)
                     sx4h=t7(neighc)
                     vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=vxyd
c     
c     determine upwind nodes and if vapour phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        vxyd=t8(neighc)
                        if(vxyd.lt.0.0) fid=dnwgt
                        if(vxyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_v(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_v(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     divi=1.0-s(i)
                     divpi=0.0
                     if(ieos(i).ne.2) then
                        divei=0.0
                     else
                        divei=-1.0
                     endif
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyh=t2(neighc)
                     pvxy=t4(neighc)
                     sx4h=t7(neighc)
                     vxyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     divkb=1.0-s(kb)
                     divpkb=0.0
                     if(ieos(kb).ne.2) then
                        divekb=0.0
                     else
                        divekb=-1.0
                     endif
                     dvpi=-pvxy+0.5*sx4h*dgvp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvpkb=pvxy+0.5*sx4h*dgvp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-
     &                    cord(iz,igrav))
                     dvekb=0.5*sx4h*dgve(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     vxyf=(fid*divkb+fid1*divi)
                     vxy=vxyd*vxyf
                     dvapi=dvpi*vxyf+vxyd*fid1*divpi
                     dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
                     dvaei=dvei*vxyf+vxyd*fid1*divei
                     dvaekb=dvekb*vxyf+vxyd*fid*divekb

                     a_vxy(iau+nmatavw)=vxy
                     a_vxy(ial+nmatavw)=-vxy

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
c     add heat conduction
c     
                     ti=t(i)
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     heatc=t5(neighc)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
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
c     end ifblock for heat/water (no air) equation type
                  else if(ico2.gt.0) then
c     
c     ifblock for heat/water/air equation type
c     
c     3-d and 2-d geometry
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iw=it10(jm)
                     pvikb=phi(kb)
                     phikb=pvikb-pcp(kb)
                     pxy=sx_md*
     &               max(pnx(i),pny(i),pnz(i),
     &               pnx(kb),pny(kb),pnz(kb))
                     pvxy=pxy   
                     pxyi=pxy*(phikb-phii)
                     pxyh=pvxy*(pvikb-pvii)
                     t1(neighc)=pxyi
                     t2(neighc)=pxyh
                     t3(neighc)=pxy
                     t4(neighc)=pvxy
                     t5(neighc)=sx_md*max(thx(i),thy(i),thz(i),
     &               thx(kb),thy(kb),thz(kb))
                     t6(neighc)=-grav*t3(neighc)
                     t7(neighc)=-grav*t4(neighc)
c     
c     liquid phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyi=t1(neighc)
                     sx4d=t6(neighc)
                     axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=axyd
c     
c     determine upwind nodes and if liquid phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        axyd=t8(neighc)
                        if(axyd.lt.0.0) fid=dnwgt
                        if(axyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_l(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_l(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     dili=s(i)
                     dilpi=0.0
                     if(ieos(i).ne.2) then
                        dilei=0.0
                     else
                        dilei=1.0
                     endif
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyi=t1(neighc)
                     pxy=t3(neighc)
                     sx4d=t6(neighc)
                     axyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     dilkb=s(kb)
                     dilpkb=0.0
                     if(ieos(kb).ne.2) then
                        dilekb=0.0
                     else
                        dilekb=1.0
                     endif
                     dlpi=-pxy+0.5*sx4d*dglp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlpkb=pxy+0.5*sx4d*dglp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlei=pxy*dpvti+0.5*sx4d*dgle(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlci=0.5*sx4d*dglc(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dlckb=0.5*sx4d*dglc(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     axyf=(fid*dilkb+fid1*dili)
                     axy=axyd*axyf
                     dlapi=dlpi*axyf+axyd*fid1*dilpi
                     dlapkb=dlpkb*axyf+axyd*fid*dilpkb
                     dlaei=dlei*axyf+axyd*fid1*dilei
                     dlaekb=dlekb*axyf+axyd*fid*dilekb
                     dlaci=dlci*axyf
                     dlackb=dlckb*axyf

                     a_axy(iau+nmatavw)=axy
                     a_axy(ial+nmatavw)=-axy

                     bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
                     bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
                     a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
                     a(jmia+nmat(2))=a(jmia+nmat(2))+dlaei
                     a(jmia+nmat(3))=a(jmia+nmat(3))+dlaci
                     a(ial+nmat(1))=a(ial+nmat(1))-dlapi
                     a(ial+nmat(2))=a(ial+nmat(2))-dlaei
                     a(ial+nmat(3))=a(ial+nmat(3))-dlaci
                     a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
                     a(iau+nmat(2))=a(iau+nmat(2))+dlaekb
                     a(iau+nmat(3))=a(iau+nmat(3))+dlackb
                     a(jml+nmat(1))=a(jml+nmat(1))-dlapkb
                     a(jml+nmat(2))=a(jml+nmat(2))-dlaekb
                     a(jml+nmat(3))=a(jml+nmat(3))-dlackb
c     
c     
c     vapour phase calculations
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     pxyh=t2(neighc)
                     sx4h=t7(neighc)
                     vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *                    *(cord(kz,igrav)-cord(iz,igrav))
                     t8(neighc)=vxyd
c     
c     determine upwind nodes and if vapour phase exists
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
c     add coding to save upwind position
                     if(iad.le.iad_up) then
                        fid=0.5
                        vxyd=t8(neighc)
                        if(vxyd.lt.0.0) fid=dnwgt
                        if(vxyd.gt.0.0) fid=upwgt
                        t9(neighc)=fid
                        call setbit(nbits,neighc,upwind_v(iz4m1),fid)
                     else
                        if(bit(nbits,neighc,upwind_v(iz4m1))) then
                           t9(neighc)=1.0
                        else
                           t9(neighc)=0.0
                        endif
                     endif
c     
c     form equations
c     
                     divi=1.0-s(i)
                     divpi=0.0
                     if(ieos(i).ne.2) then
                        divei=0.0
                     else
                        divei=-1.0
                     endif
c     
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     pxyh=t2(neighc)
                     pvxy=t4(neighc)
                     sx4h=t7(neighc)
                     vxyd=t8(neighc)
                     fid=t9(neighc)
                     fid1=1.0-fid
                     divkb=1.0-s(kb)
                     divpkb=0.0
                     if(ieos(kb).ne.2) then
                        divekb=0.0
                     else
                        divekb=-1.0
                     endif
                     dvpi=-pvxy+0.5*sx4h*dgvp(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvpkb=pvxy+0.5*sx4h*dgvp(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-
     &                    cord(iz,igrav))
                     dvekb=0.5*sx4h*dgve(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvci=0.5*sx4h*dgvc(i)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     dvckb=0.5*sx4h*dgvc(kb)*
     *                    (cord(kz,igrav)-cord(iz,igrav))
                     vxyf=(fid*divkb+fid1*divi)
                     vxy=vxyd*vxyf
                     dvapi=dvpi*vxyf+vxyd*fid1*divpi
                     dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
                     dvaei=dvei*vxyf+vxyd*fid1*divei
                     dvaekb=dvekb*vxyf+vxyd*fid*divekb
                     dvaci=dvci*vxyf
                     dvackb=dvckb*vxyf

                     a_vxy(iau+nmatavw)=vxy
                     a_vxy(ial+nmatavw)=-vxy

                     bp(iz+nrhs(1))=bp(iz+nrhs(1))+vxy
                     bp(kz+nrhs(1))=bp(kz+nrhs(1))-vxy
                     a(jmia+nmat(1))=a(jmia+nmat(1))+dvapi
                     a(jmia+nmat(2))=a(jmia+nmat(2))+dvaei
                     a(jmia+nmat(3))=a(jmia+nmat(3))+dvaci
                     a(ial+nmat(1))=a(ial+nmat(1))-dvapi
                     a(ial+nmat(2))=a(ial+nmat(2))-dvaei
                     a(ial+nmat(3))=a(ial+nmat(3))-dvaci
                     a(iau+nmat(1))=a(iau+nmat(1))+dvapkb
                     a(iau+nmat(2))=a(iau+nmat(2))+dvaekb
                     a(iau+nmat(3))=a(iau+nmat(3))+dvackb
                     a(jml+nmat(1))=a(jml+nmat(1))-dvapkb
                     a(jml+nmat(2))=a(jml+nmat(2))-dvaekb
                     a(jml+nmat(3))=a(jml+nmat(3))-dvackb
c     
c     
c     add air diffusion  
c     
                     cnvi=cnvf(i)+cnlf(i)
                     dcvi=dcvf(i)+dclf(i)
                     dcvei=dcvef(i)+dclef(i)
                     dcvci=dcvcf(i)+dclcf(i)
                     jm=1
                     kb=it8(jm)
                     kz=kb-icd
                     neighc=it9(jm)
                     heatc=t5(neighc)
                     iau=it11(jm)
                     ial=it12(jm)
                     jml=nelmdg(kz)-neqp1
                     heatc=t5(neighc)
                     cnvkb=cnvf(kb)+cnlf(kb)
                     dcvkb=dcvf(kb)+dclf(kb)
                     dcvekb=dcvef(kb)+dclef(kb)
                     dcvckb=dcvcf(kb)+dclcf(kb)
                     bp(iz+nrhs(3))=bp(iz+nrhs(3))+heatc*(cnvkb-cnvi)
                     bp(kz+nrhs(3))=bp(kz+nrhs(3))-heatc*(cnvkb-cnvi)
                     a(jmia+nmat(7))=a(jmia+nmat(7))-heatc*dcvi    
                     a(jmia+nmat(8))=a(jmia+nmat(8))-heatc*dcvei    
                     a(jmia+nmat(9))=a(jmia+nmat(9))-heatc*dcvci   
                     a(ial+nmat(7))=a(ial+nmat(7))+heatc*dcvi    
                     a(ial+nmat(8))=a(ial+nmat(8))+heatc*dcvei    
                     a(ial+nmat(9))=a(ial+nmat(9))+heatc*dcvci    
                     a(iau+nmat(7))=a(iau+nmat(7))+heatc*dcvkb    
                     a(iau+nmat(8))=a(iau+nmat(8))+heatc*dcvekb    
                     a(iau+nmat(9))=a(iau+nmat(9))+heatc*dcvckb      
                     a(jml+nmat(7))=a(jml+nmat(7))-heatc*dcvkb    
                     a(jml+nmat(8))=a(jml+nmat(8))-heatc*dcvekb    
                     a(jml+nmat(9))=a(jml+nmat(9))-heatc*dcvckb    
c     
c     add heat conduction
c     
                     ti=t(i)
                     jm=1
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
                     a(jmia+nmat(4))=a(jmia+nmat(4))-heatc*dtpa(i)
                     a(jmia+nmat(5))=a(jmia+nmat(5))-heatc*dtpae(i)
                     a(jmia+nmat(6))=a(jmia+nmat(6))-heatc*dtpac(i)
                     a(ial+nmat(4))=a(ial+nmat(4))+heatc*dtpa(i)
                     a(ial+nmat(5))=a(ial+nmat(5))+heatc*dtpae(i)
                     a(ial+nmat(6))=a(ial+nmat(6))+heatc*dtpac(i)
                     a(iau+nmat(4))=a(iau+nmat(4))+heatc*dtpa(kb)
                     a(iau+nmat(5))=a(iau+nmat(5))+heatc*dtpae(kb)
                     a(iau+nmat(6))=a(iau+nmat(6))+heatc*dtpac(kb)
                     a(jml+nmat(4))=a(jml+nmat(4))-heatc*dtpa(kb)
                     a(jml+nmat(5))=a(jml+nmat(5))-heatc*dtpae(kb)
                     a(jml+nmat(6))=a(jml+nmat(6))-heatc*dtpac(kb)
c     end ifblock for heat/water/air equation type
                  endif
               else
c     heat conduction only(idof=1)
c     start ifblock for heat/water/air equation type
c     
c     3-d and 2-d geometry
c     
c     
c     form equations
c     
                  ti=t(i)
                  jm=1
                  kb=it8(jm)
                  kz=kb-icd
                  t5(neighc)=sx_md*max(thx(i),thy(i),thz(i),
     &             thx(kb),thy(kb),thz(kb))
                  neighc=it9(jm)
                  heatc=t5(neighc)
                  iau=it11(jm)
                  ial=it12(jm)
                  jml=nelmdg(kz)-neqp1
                  heatc=t5(neighc)
                  bp(i)=bp(i)+heatc*(t(kb)-ti)
                  bp(kb)=bp(kb)-heatc*(t(kb)-ti)
                  a(jmia)=a(jmia)-heatc*dtpae(i)
                  a(ial)=a(ial)+heatc*dtpae(i)
                  a(iau)=a(iau)+heatc*dtpae(kb)
                  a(jml)=a(jml)-heatc*dtpae(kb)
c     end ifblock for heat conduction only
               endif
c     end ifblock for mdnodes
            endif
c     end do loop on neighbors of md nodes
         enddo
c     end do loop on nodes
      enddo
      
      if(idpdp.ne.0) then
         do i=1,idof*idof
            nmat(i)=(i-1)*ldna
         enddo
         do i=1,idof
            nrhs(i)=(i-1)*neq
         enddo
      endif

      return
      end


