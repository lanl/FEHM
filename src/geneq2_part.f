      subroutine geneq2_part(i,zone)
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

!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To generate equations isothermal air-water solution at each node.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/geneq2_part.f_a  $
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
!D4 SPECIAL COMMENTS
!D4
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
!D4    FEHM Application Version 2.20
!D4
!**********************************************************************
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
      use com_part
      implicit none

      integer i
      integer zone
      integer i1  
      integer i2
      integer iq_boun
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
      integer kb_new
      integer neighc
      integer iau
      integer ial
      integer kz
      integer nmatavw
      integer old_node_i
      integer old_node_iz
      integer old_node_kb
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
      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    
c gaz temporary array
      integer it14(326)
      integer it15(326)
      integer jmia_old
      integer neqp1_old
      integer neq_primary_new
      integer iq_old,iau_old,ial_old
      integer ij3,ij4,ij5,kk              

c changed by avw -- entered here by seh
      neqp1=neq_part(zone)+1
      if(i.gt.neq_part(zone)) then
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
c storage for upwind
c

      old_node_i = index_part(zone,i)
      sx1d=sx1(old_node_i)
      axi=pnx(old_node_i)
      ayi=pny(old_node_i)
      azi=pnz(old_node_i)
      alxi=axi
      avxi=axi
      alyi=ayi
      avyi=ayi
      alzi=azi
      avzi=azi
      pvii=phi(old_node_i)
      dili=dil(old_node_i)
      dilpi=dilp(old_node_i)
      if(irdof.ne.13) then
        phii=pvii-pcp(old_node_i)
        dpvti=dpcef(old_node_i)
        divi=div(old_node_i)
        divpi=divp(old_node_i)
        divei=dive(old_node_i)
        dilei=dile(old_node_i)
      else
        phii=pvii
      endif
      ti=t(old_node_i)
      swi=s(old_node_i)
c
c form constants for i>neq
c
      neq_primary_new=neq
      if(i.gt.neq_primary_new.and.idualp.eq.0) then
         icd=neq_primary_new
      else
         icd=0
      endif
      iz=i-icd
c
      iz4m1 = 4*(iz-1)+1
c
      ii1=nelm_part(zone,i-icd)+1
      ii2=nelm_part(zone,i-icd+1)
      idg=nelmdg_part(zone,i-icd)-ii1+1
      neqp1=neq_part(zone)+1
      iq=0
      jmi=nelmdg_part(zone,i-icd)
      jmia=jmi-neqp1

      old_node_iz = index_part(zone,iz)
c
c Still need some old number connections for filling in a_axy
c for transport calculations
c
      neqp1_old = neq_primary+1
      jmia_old =  nelmdg(old_node_iz)-neqp1_old
c Save source/sink term by putting it into empty slot (from 
c same node to same node). 
c If this is an isothermal air-water simulation
      a_axy(jmia_old+nmatavw)=sk(old_node_i)
c Take care of source/sink term 
c If this is an isothermal air-water simulation
      a_vxy(jmia_old+nmatavw)=qh(old_node_i)


      do 58 jm=jmi+1,ii2
         iq=iq+1
c this kb is in new number system
         kb=nelm_part(zone,jm)+icd
         it9(iq) = kb 
c this kb is in old number system
         kb=index_part(zone,nelm_part(zone,jm))+icd
         it8(iq)=kb
         it10(iq)=istrw_part(zone,jm-neqp1)
         it11(iq)=jm-neqp1
         it13(iq)=nelmdg_part(zone,nelm_part(zone,jm))-neqp1
         ij1=nelm_part(zone,nelm_part(zone,jm))+1
         ij2=nelmdg_part(zone,nelm_part(zone,jm))-1
         do 68 ij=ij1,ij2
            if(nelm_part(zone,ij).eq.iz) then
               it12(iq)=ij-neqp1
            endif
 68      continue
 58   continue
c
c add boundary cell information
c kb is in old number system
c
        iq_boun = iq
        i1 = bound_part(zone,iz) + 1
        i2 = bound_part(zone,iz + 1)
         do jm = i1,i2
           iq = iq + 1
           kb = bound_part(zone,jm)
           iw = bound_istrw(zone,jm-neqp1)
           it8(iq) = kb
           it10(iq) =iw
         enddo
c
c Need connections in old numbering for mass balance calculation
c
        iq_old = 0
        ij1 = nelm(old_node_i)+1
        ij2 = nelmdg(old_node_i)
        ij3 = nelm(old_node_i+1)
        do jm=1,iq                                        
         iq_old = iq_old+1
         kb = it8(jm)
         if(kb.gt.old_node_i) then
          do ij = ij2+1,ij3
            if(nelm(ij).eq.kb) then
               it14(iq_old)=ij-neqp1_old
               ij4 = nelm(kb)+1
               ij5 = nelmdg(kb) - 1
               do kk = ij4,ij5
                if(nelm(kk).eq.old_node_i) then
                   it15(iq_old)=kk-neqp1_old
                endif
               enddo
            endif
          enddo
         else if(kb.lt.old_node_i) then
          do ij = ij1,ij2-1
            if(nelm(ij).eq.kb) then
               it14(iq_old)=ij-neqp1_old
            endif
          enddo
         endif
        enddo

c
c 3-d geometry
c
      if(icnl.eq.0) then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iwd=it10(jm)
            iw =abs(iwd)
            axkb=pnx((kb))
            aykb=pny((kb))
            azkb=pnz((kb))
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
             reduction_factor = 1.0 
c     &              red_factor(istrw_itfc(it11(jm)))
             perml(1)=2.*alxkb*alxi/(alxkb+alxi)
             perml(2)=2.*alykb*alyi/(alykb+alyi)
             perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,isox)+
     &           sx(iw,isoy)+
     &           sx(iw,isoz)
            pvikb=phi((kb))
            phikb=pvikb-pcp((kb))
            delx2=(cord((kz),1)-
     &             cord((old_node_iz),1))**2
            dely2=(cord((kz),2)-
     &             cord((old_node_iz),2))**2
            delz2=(cord((kz),3)-
     &             cord((old_node_iz),3))**2
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
            t1(jm)=pxyi
            t2(jm)=pxyh
            t3(jm)=pxy
            t4(jm)=pxy
            t6(jm)=-grav*t3(jm)
            t7(jm)=-grav_air*t4(jm)
 59      continue
      elseif(icnl.ne.0) then
c     
c 2-d geometry
c
         radi=cord(iz,3)
         do 69 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iwd=it10(jm)
            iw = abs(iwd)
c           jm=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
             reduction_factor = 1.0 
c     &               red_factor(istrw_itfc(index_part(zone,it11(jm))))
             perml(1)=2.*alxkb*alxi/(alxkb+alxi)
             perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
c           sx2c=radkb*(sx(index_part(zone,iw),isox)+
c    &                  sx(index_part(zone,iw),isoy))
            sx2c=radkb*(sx(iw,isox)+
     &                  sx(iw,isoy))
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-
     &             cord(old_node_iz,1))**2
            dely2=(cord(kz,2)-
     &             cord(old_node_iz,2))**2
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
            t1(jm)=pxyi
            t2(jm)=pxyh
            t3(jm)=pxy
            t4(jm)=pxy
            t6(jm)=-grav*t3(jm)
            t7(jm)=-grav_air*t4(jm)
 69      continue
      endif
      if(irdof.ne.11) then
c     
c liquid phase calculations
c
         do 60 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
c           jm=it9(jm)
            pxyi=t1(jm)
            sx4d=t6(jm)
            axyd=pxyi+0.5*sx4d*(rolf(old_node_i)+
     *                          rolf(kb))
     *           *(cord(kz,igrav)-
     *             cord(old_node_iz,igrav))
            t8(jm)=axyd
 60      continue
c
c determine upwind nodes and if liquid phase exists
c
         isl=1
         do 61 jm=1,iq
c           jm=it9(jm)
c add coding to save upwind position
         if(iad.le.iad_up) then
            fid=0.5
            axyd=t8(jm)
            if(axyd.lt.0.0) fid=dnwgt
            if(axyd.gt.0.0) fid=upwgt
            t9(jm)=fid
c
            call setbit(nbits,jm,
     &                  upwind_l(iz4m1),fid)
c
         else
           if(bit(nbits,jm,
     &                  upwind_l(iz4m1))) then 
            t9(jm)=1.0
           else
            t9(jm)=0.0
           endif
         endif
 61      continue
c
c form equations
c
         if(isl.ne.0) then
            do 62 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               kb_new = it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=it13(jm)                          
               iau_old = it14(jm)
               ial_old = it15(jm)
               axyd=t8(jm)
               fid=t9(jm)
               fid1=1.0-fid
               pxyi=t1(jm)
               pxy=t3(jm)
               sx4d=t6(jm)
               dilkb=dil(kb)
               dilpkb=dilp(kb)
               dlpi=-pxy+0.5*sx4d*dglp(old_node_i)*
     2              (cord(kz,igrav)-
     3               cord(old_node_iz,igrav))
               dlpkb=pxy+0.5*sx4d*dglp(kb)*
     2              (cord(kz,igrav)-
     3               cord(old_node_iz,igrav))
c
               axyf=(fid*dilkb+fid1*dili)
               axy=axyd*axyf
               if(irdof.ne.13) then
                dilekb=dile(kb)
                dlei=pxy*dpvti+0.5*sx4d*dgle(old_node_i)*
     2               (cord(kz,igrav)-
     3                cord(old_node_iz,igrav))
                dlekb=-pxy*dpcef(kb)+0.5*sx4d*
     *                dgle(kb)
     *               *(cord(kz,igrav)-
     *                 cord(old_node_iz,igrav))
                dlaei=dlei*axyf+axyd*fid1*dilei
                dlaekb=dlekb*axyf+axyd*fid*dilekb
               endif
c
               dlapi=dlpi*axyf+axyd*fid1*dilpi
               dlapkb=dlpkb*axyf+axyd*fid*dilpkb
c     
c              a_axy(iau+nmatavw)=axy
c              a_axy(ial+nmatavw)=-axy
c  axy calcs in loops below

             if(jm.le.iq_boun) then
               a_axy(iau_old+nmatavw)=axy
               a_axy(ial_old+nmatavw)=-axy
               bp_part(zone,iz+nrhs_part(zone,1))=
     &                    bp_part(zone,iz+nrhs_part(zone,1))+axy
               bp_part(zone,kb_new+nrhs_part(zone,1))=
     &                    bp_part(zone,kb_new+nrhs_part(zone,1))-axy
               a(jmia+nmat_part(zone,1))=
     &                    a(jmia+nmat_part(zone,1))+dlapi
               a(ial+nmat_part(zone,1))=
     &                    a(ial+nmat_part(zone,1))-dlapi
               a(iau+nmat_part(zone,1))=
     &                    a(iau+nmat_part(zone,1))+dlapkb
               a(jml+nmat_part(zone,1))=
     &                    a(jml+nmat_part(zone,1))-dlapkb
               if(irdof.ne.13) then
                 a(jmia+nmat_part(zone,2))=
     &                    a(jmia+nmat_part(zone,2))+dlaei
                 a(ial+nmat_part(zone,2))=
     &                    a(ial+nmat_part(zone,2))-dlaei
                 a(iau+nmat_part(zone,2))=
     &                    a(iau+nmat_part(zone,2))+dlaekb
                 a(jml+nmat_part(zone,2))=
     &                    a(jml+nmat_part(zone,2))-dlaekb
               endif
             else
               a(jmia+nmat_part(zone,1))=
     &                    a(jmia+nmat_part(zone,1))+dlapi*overf
                 if(irdof.ne.13) then
                   a(jmia+nmat_part(zone,2))=
     &                    a(jmia+nmat_part(zone,2))+dlaei*overf
                 endif
               a_axy(iau_old+nmatavw)=axy
               bp_part(zone,iz+nrhs_part(zone,1))=
     &                    bp_part(zone,iz+nrhs_part(zone,1))+axy

             endif
 62         continue
         endif
         endif
       if(irdof.ne.13) then   
c     
c vapour phase calculations
c
         do 63 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
c           jm=it9(jm)
            pxyh=t2(jm)
            sx4h=t7(jm)
            vxyd=pxyh+0.5*sx4h*(rovf(old_node_i)+
     *                          rovf(kb))
     *           *(cord(kz,igrav)-
     *             cord(old_node_iz,igrav))
            t8(jm)=vxyd
 63      continue
c     
c     determine upwind nodes and if vapour phase exists
c     
         do 64 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
c           jm=it9(jm)
c add coding to save upwind position
         if(iad.le.iad_up) then
            fid=0.5
            vxyd=t8(jm)
            if(vxyd.lt.0.0) fid=dnwgt
            if(vxyd.gt.0.0) fid=upwgt
            t9(jm)=fid
c
            call setbit(nbits,jm,
     *                  upwind_v(iz4m1),fid)
c
         else
           if(bit(nbits,jm,
     *                  upwind_v(iz4m1))) then 
            t9(jm)=1.0
           else
            t9(jm)=0.0
           endif
         endif
 64      continue
c
c form equations
c
         isl=0
         if(isl.ne.1) then
            do 65 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               kb_new = it9(jm)
               iau=it11(jm)
               ial=it12(jm)           
               iau_old=it14(jm)
               ial_old=it15(jm)           
               jml=it13(jm)                          
               fid=t9(jm)
               fid1=1.0-fid
               pxyh=t2(jm)
               pvxy=t4(jm)
               sx4h=t7(jm)
               vxyd=t8(jm)
               divkb=div(kb)
               divpkb=divp(kb)
               divekb=dive(kb)
               dvpi=-pvxy+0.5*sx4h*dgvp(old_node_i)*
     2              (cord(kz,igrav)-
     3               cord(old_node_iz,igrav))
               dvpkb=pvxy+0.5*sx4h*dgvp(kb)*
     2              (cord(kz,igrav)-
     3               cord(old_node_iz,igrav))
               dvei=0.5*sx4h*dgve(old_node_i)*
     2              (cord(kz,igrav)-
     3               cord(old_node_iz,igrav))
               dvekb=0.5*sx4h*dgve(kb)
     *              *(cord(kz,igrav)-
     *                cord(old_node_iz,igrav))
               vxyf=(fid*divkb+fid1*divi)
               vxy=vxyd*vxyf
               dvapi=dvpi*vxyf+vxyd*fid1*divpi
               dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
               dvaei=dvei*vxyf+vxyd*fid1*divei
               dvaekb=dvekb*vxyf+vxyd*fid*divekb
c     
c              a_vxy(iau+nmatavw)=vxy
c              a_vxy(ial+nmatavw)=-vxy

             if(jm.le.iq_boun) then
               a_vxy(iau_old+nmatavw)=vxy
               a_vxy(ial_old+nmatavw)=-vxy
               bp_part(zone,iz+nrhs_part(zone,2))=
     &                   bp_part(zone,iz+nrhs_part(zone,2))+vxy
               bp_part(zone,kb_new+nrhs_part(zone,2))=
     &                   bp_part(zone,kb_new+nrhs_part(zone,2))-vxy
               a(jmia+nmat_part(zone,3))=
     &                   a(jmia+nmat_part(zone,3))+dvapi
               a(jmia+nmat_part(zone,4))=
     &                   a(jmia+nmat_part(zone,4))+dvaei
               a(ial+nmat_part(zone,3))=
     &                   a(ial+nmat_part(zone,3))-dvapi
               a(ial+nmat_part(zone,4))=
     &                   a(ial+nmat_part(zone,4))-dvaei
               a(iau+nmat_part(zone,3))=
     &                   a(iau+nmat_part(zone,3))+dvapkb
               a(iau+nmat_part(zone,4))=
     &                   a(iau+nmat_part(zone,4))+dvaekb
               a(jml+nmat_part(zone,3))=
     &                   a(jml+nmat_part(zone,3))-dvapkb
               a(jml+nmat_part(zone,4))=
     &                   a(jml+nmat_part(zone,4))-dvaekb
             else                    
               a(jmia+nmat_part(zone,3))=
     &                    a(jmia+nmat_part(zone,3))+dvapi
               a(jmia+nmat_part(zone,4))=
     &                    a(jmia+nmat_part(zone,4))+dvaei
               a_vxy(iau_old+nmatavw)=vxy
               bp_part(zone,iz+nrhs_part(zone,2))=
     &                    bp_part(zone,iz+nrhs_part(zone,2))+vxy
             endif                     
 65         continue
         endif
       endif
      return
      end

