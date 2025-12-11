      subroutine geneq2_uz_wt(i) 
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
CD1 To generate equations isothermal air-water solution at each node.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 22-May-02, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/geneq2_uz_wt.f_a  $
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.3.3 Noncondensible gas flow equations        
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************
c
c generate equations for 3-d schemes'isothermal air-water ,
c full derivatives'
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
      real*8 axyds
      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd,i_dir_gdkm,kb_pri    

c     changed by avw -- entered here by seh
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
c     storage for upwind
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
      swi=s(i)
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
      iz4m1 = 4*(iz-1)+1
c     
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
      a_axy(jmia+nmatavw)=sk(i)
c     Take care of source/sink term 
c     If this is an isothermal air-water simulation
      a_vxy(jmia+nmatavw)=qh(i)

      do 58 jm=jmi+1,ii2
         iq=iq+1
         kb=nelm(jm)+icd
         it8(iq)=kb
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
c     gaz 11-18-2001
c     if(imdnode.ne.0) then
c     imd = mdnodes(i) + mdnodes(kb)
c     if(imd.lt.2) it10(iq) = -abs(it10(iq))
c     endif
         it11(iq)=jm-neqp1
         ij1=nelm(nelm(jm))+1
         ij2=nelmdg(nelm(jm))-1
         do 68 ij=ij1,ij2
            if(nelm(ij).eq.iz) then
               it12(iq)=ij-neqp1
            endif
 68      continue
 58   continue
c     
c     3-d geometry
c     
      if(icnl.eq.0) then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
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
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(i_dir_gdkm.ge.0.and.reduction_factor.gt.2) then
               kb_pri = reduction_factor -2
               reduction_factor = 1.0 
c  gaz 050118 harmonic weighting to match hi res grid               
               if(i_dir_gdkm.eq.1) then
                 pxy = sx2c*perml(1)
               else if(i_dir_gdkm.eq.2) then
                 pxy = sx2c*perml(2) 
               else if(i_dir_gdkm.eq.3) then
                 pxy = sx2c*perml(3)
               else if(dis2.gt.dis_tol) then
                pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
               endif                                    
            elseif(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
            endif
            if(reduction_factor.gt.2) reduction_factor = 1.0
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue
         if(irdof.ne.11) then
c     
c     liquid phase calculations
c     
            do 60 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyi=t1(neighc)
               pxy = t3(neighc)
               sx4d=t6(neighc)
               axyds = 0.5*pxy*(head12(kb,1)+head12(i,1))
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
     &              + axyds*(s(kb)-s(i))         
               t8(neighc)=axyd
               t10(neighc)=axyds 
 60         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
            isl=1
            do 61 jm=1,iq
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
 61         continue
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
                  axyds = t10(neighc)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  pxyi=t1(neighc)
                  pxy=t3(neighc)
                  sx4d=t6(neighc)
                  dilkb=dil(kb)
                  dilpkb=dilp(kb)
                  dlpi=-pxy+0.5*sx4d*dglp(i)*
     2                 (cord(kz,igrav)-cord(iz,igrav))
                  dlpkb=pxy+0.5*sx4d*dglp(kb)*
     2                 (cord(kz,igrav)-cord(iz,igrav))
c     
                  axyf=(fid*dilkb+fid1*dili)
                  axy=axyd*axyf
                  if(irdof.ne.13) then
                     dilekb=dile(kb)
                     dlei=pxy*dpvti+0.5*sx4d*dgle(i)*
     2                    (cord(kz,igrav)-cord(iz,igrav))
     &                    - axyds
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
     &                    + axyds
                     dlaei=dlei*axyf+axyd*fid1*dilei
                     dlaekb=dlekb*axyf+axyd*fid*dilekb
                  endif
c     
                  dlapi=dlpi*axyf+axyd*fid1*dilpi
                  dlapkb=dlpkb*axyf+axyd*fid*dilpkb
c     
                  a_axy(iau+nmatavw)=axy
                  a_axy(ial+nmatavw)=-axy

                  bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
                  bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
                  a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
                  a(ial+nmat(1))=a(ial+nmat(1))-dlapi
                  a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
                  a(jml+nmat(1))=a(jml+nmat(1))-dlapkb
                  if(irdof.ne.13) then 
                     a(jmia+nmat(2))=a(jmia+nmat(2))+dlaei

                     a(ial+nmat(2))=a(ial+nmat(2))-dlaei

                     a(iau+nmat(2))=a(iau+nmat(2))+dlaekb

                     a(jml+nmat(2))=a(jml+nmat(2))-dlaekb

                  endif
 62            continue
            endif
         endif
         if(irdof.ne.13) then   
c     
c     vapour phase calculations
c     
            do 63 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
 63         continue
c     
c     determine upwind nodes and if vapour phase exists
c     
            isl=1
            do 64 jm=1,iq
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
 64         continue
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
                  dvpi=-pvxy+0.5*sx4h*dgvp(i)*
     2                 (cord(kz,igrav)-cord(iz,igrav))
                  dvpkb=pvxy+0.5*sx4h*dgvp(kb)*
     2                 (cord(kz,igrav)-cord(iz,igrav))
                  dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
                  dvekb=0.5*sx4h*dgve(kb)
     *                 *(cord(kz,igrav)-cord(iz,igrav))
                  vxyf=(fid*divkb+fid1*divi)
                  vxy=vxyd*vxyf
                  dvapi=dvpi*vxyf+vxyd*fid1*divpi
                  dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
                  dvaei=dvei*vxyf+vxyd*fid1*divei
                  dvaekb=dvekb*vxyf+vxyd*fid*divekb
c     
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

 65            continue
            endif
         endif
c     
c     
c     2-d geometry
c     
      elseif(icnl.ne.0) then
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
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
c gaz 051818            
            if(i_dir_gdkm.ge.0.and.reduction_factor.gt.2) then
               kb_pri = reduction_factor -2
               reduction_factor = 1.0 
               if(i_dir_gdkm.eq.1) then
                pxy = sx2c*perml(1)
               else if(i_dir_gdkm.eq.2) then
                pxy = sx2c*perml(2)
               else if(dis2.gt.dis_tol) then 
                pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2))
               endif                                
            elseif(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
            if(reduction_factor.gt.2) reduction_factor = 1.0
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 69      continue
         if(irdof.ne.11) then
c     
c     liquid phase calculations
c     
            do 70 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyi=t1(neighc)
               pxy = t3(neighc)
               sx4d=t6(neighc)
               axyds = 0.5*pxy*(head12(kb,1)+head12(i,1))
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
     &              + axyds*(s(kb)-s(i))         
               t8(neighc)=axyd
               t10(neighc)=axyds 
 70         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
            isl=1
            do 71 jm=1,iq
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
 71         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 72 jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  iau=it11(jm)
                  ial=it12(jm)
                  jml=nelmdg(kz)-neqp1
                  heatc=t5(neighc)
                  axyd=t8(neighc)
                  axyds = t10(neighc)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  pxyi=t1(neighc)
                  pxy=t3(neighc)
                  sx4d=t6(neighc)
                  dilkb=dil(kb)
                  dilpkb=dilp(kb)
                  dlpi=-pxy+0.5*sx4d*dglp(i)*(cord(kz,igrav)-
     2                 cord(iz,igrav))
                  dlpkb=pxy+0.5*sx4d*dglp(kb)*
     2                 (cord(kz,igrav)-cord(iz,igrav))
c     
                  axyf=(fid*dilkb+fid1*dili)
                  axy=axyd*axyf
                  if(irdof.ne.13) then
                     dilekb=dile(kb)
                     dlei=pxy*dpvti+0.5*sx4d*dgle(i)*
     2                    (cord(kz,igrav)-cord(iz,igrav))
     &                    - axyds
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
     &                    + axyds
                     dlaei=dlei*axyf+axyd*fid1*dilei
                     dlaekb=dlekb*axyf+axyd*fid*dilekb
                  endif
c     
                  dlapi=dlpi*axyf+axyd*fid1*dilpi
                  dlapkb=dlpkb*axyf+axyd*fid*dilpkb
c     
                  a_axy(iau+nmatavw)=axy
                  a_axy(ial+nmatavw)=-axy

                  bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
                  bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
                  a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
                  a(ial+nmat(1))=a(ial+nmat(1))-dlapi
                  a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
                  a(jml+nmat(1))=a(jml+nmat(1))-dlapkb

                  if(irdof.ne.13) then
                     a(jmia+nmat(2))=a(jmia+nmat(2))+dlaei
                     a(ial+nmat(2))=a(ial+nmat(2))-dlaei
                     a(iau+nmat(2))=a(iau+nmat(2))+dlaekb
                     a(jml+nmat(2))=a(jml+nmat(2))-dlaekb
                  endif
 72            continue
            endif
         endif
         if(irdof.ne.13) then
c     
c     vapour phase calculations
c     
            do 73 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
 73         continue
c     
c     determine upwind nodes and if vapour phase exists
c     
            isl=1
            do 74 jm=1,iq
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
 74         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 75 jm=1,iq
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
                  dvpi=-pvxy+0.5*sx4h*dgvp(i)*(cord(kz,igrav)-
     2                 cord(iz,igrav))
                  dvpkb=pvxy+0.5*sx4h*dgvp(kb)*(cord(kz,igrav)-
     2                 cord(iz,igrav))
                  dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
                  dvekb=0.5*sx4h*dgve(kb)
     *                 *(cord(kz,igrav)-cord(iz,igrav))
                  vxyf=(fid*divkb+fid1*divi)
                  vxy=vxyd*vxyf
                  dvapi=dvpi*vxyf+vxyd*fid1*divpi
                  dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
                  dvaei=dvei*vxyf+vxyd*fid1*divei
                  dvaekb=dvekb*vxyf+vxyd*fid*divekb
c     
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
 75            continue
            endif
         endif
c     
      endif
      
      return
      end


