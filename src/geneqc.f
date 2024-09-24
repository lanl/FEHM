      subroutine geneqc(i)
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
CD1  This subroutine generates the equations for 3-dimensional heat 
CD1  and mass transfer with non-condensible gas included (full
CD1  derivatives).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/geneqc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.13   Mon Mar 31 08:37:12 1997   gaz
CD2 minor changes for anisotropic properties
CD2 
CD2    Rev 1.12   Thu Feb 15 10:43:08 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.11   Wed Feb 07 11:03:24 1996   gaz
CD2 corrected radial term, storage for fluxes(a_lxy,a_vxy)
CD2 
CD2    Rev 1.10   Wed Jan 10 13:39:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.9   12/13/95 10:28:44   robinson
CD2 Incorporated new zeolite hydration module
CD2 
CD2    Rev 1.8   12/13/95 08:40:32   gaz
CD2 changed setbit counted to accomodate nbits=256
CD2 
CD2    Rev 1.7   08/07/95 11:51:14   awolf
CD2 Fixed for DPDP to account for all a_axy a_vxy
CD2 
CD2    Rev 1.6   06/21/95 10:36:36   llt
CD2 declared bit function as a logical for IBM
CD2 
CD2    Rev 1.5   06/01/95 16:53:20   gaz
CD2 upwind bitmap
CD2 
CD2    Rev 1.4   04/03/95 08:46:14   robinson
CD2 Corrected the source/sink term in mass flow rate array
CD2 
CD2    Rev 1.3   03/10/95 10:44:02   llt
CD2 modified to allow for Richard's Equation - gaz
CD2 
CD2    Rev 1.2   01/28/95 13:54:52   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.1   03/18/94 15:51:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:16   pvcs
CD2 original version in process of being certified
CD2
c 11/21/94 gaz
c water balance equation 
c 1/5/95 gaz added 1-xnl term here  l 284-289  
c 1/5/95 gaz added 1-xnv term here  l 419-430  
c 1/23/95 gaz harmonic weight of dva terms
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
C***********************************************************************
C---------- PHS 11 March 2004
c----    added comrxni so we have the mw_water and mw_air terms
c----    adding xnvmol as the mol fraction of the vapor phase
c----    adding in adif section two variables delmpv_wv delmpv_air
c-----    which are the delta moles per volume water vapor and air.
c------- adding varables
c          mpv_airi  mpv_airkb  mpv_wvi  mpv_wvkb
c          moles per volume air and watervapor for kb and i
c          
C -------- PHS  21 Aug 2003 
C ---  Removed coupling with gas diffusion in ADIF case

      use comzeoli
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
      use comrxni
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, icesd
      integer ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
      integer imd,iwd
      real*8 reduction_factor, reduction_factor_t
      real*8  acxy, acxyf
      real*8  aexy, aexyf, alxi, alxkb, alyi, alykb, alzi, alzkb
      real*8  avxi, avyi, avzi, axi, axkb, axy, axyd, axyf
      real*8  ayi, aykb, azi, azkb
      real*8  cnli, cnlkb, cnvi, cnvkb
      real*8  dclei, dclci, dclckb, dclekb, dcli, dclkb
      real*8  dcvei, dcvci, dcvckb, dcvekb, dcvi, dcvkb
      real*8  ddva_barci, ddva_barckb, ddva_barei, ddva_barekb 
      real*8  ddva_barpi, ddva_barpkb
      real*8  delci, delckb 
      real*8  delei, delekb, deli, delkb, delx2, dely2, delz2
      real*8  devei, devci, devckb, devekb, devi, devkb 
      real*8  devmci, devmei, devmpi
      real*8  dis2,dis_tol,sx_min 
      real*8  dilci, dilckb
      real*8  dilei, dilekb, dili, dilkb, dilpi, dilpkb
      real*8  divci, divckb
      real*8  divei, divekb, divi, divkb, divpi, divpkb
      real*8  dlaei, dlaekb, dlapi, dlapkb, dleei, dleekb 
      real*8  dlaci, dlackb, dlcci, dlcckb, dlcei, dlcekb 
      real*8  dlci, dlckb, dlcpi, dlcpkb
      real*8  dleci, dleckb 
      real*8  dlei, dlekb, dlepi, dlepkb, dlpi, dlpkb, dpdt
      real*8  dpvti, dqhzt, dqhzt1, dqhzt2, dskzt, dskzt1, dskzt2   
      real*8  dva_bar, dvaei, dvaekb, dvapi, dvapkb, dveei, dveekb
      real*8  dvaci, dvackb, dvai, dvakb, dvame, dvamet
      real*8  dvcci, dvcckb 
      real*8  dvcei, dvcekb, dvci, dvckb, dveci, dveckb, dvcpi, dvcpkb
      real*8  dvei, dvekb, dvepi, dvepkb 
      real*8  dvmckb, dvmekb, dvmpkb, dvpi, dvpkb
      real*8  fid, fid1, grav_air, heatc 
      real*8  enli, enlkb, envi, envkb, envmi, envmkb
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyi, pxyh
      real*8  psat, qhz, radi, radkb, skz
      real*8  swi, sx1d, sx2c, sx2t, sx3c, sx3t, sx4d, sx4h, sxzt
      real*8  thxi, thxkb, thyi, thykb, thzi, thzkb, ti
      real*8  vcxy, vcxyf, vexy, vexyf, vxy, vxyd, vxyf, water_pressure

      real*8  mpv_airi,mpv_airkb,mpv_wvi,mpv_wvkb,delmpv_air,delmpv_wv 

      real*8 heatt
      
      integer kb_pri, i_dir_gdkm

      parameter(dis_tol=1.d-12)

c gaz debug 100720     
c      kb = l+ pci(1)+pcio(1)
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
      phii=pvii-pcp(i)
      dpvti=dpcef(i)
      enli=enlf(i)
      dilpi=dilp(i)
      dilei=dile(i)
      divpi=divp(i)
      divei=dive(i)
      deli=delf(i)
      delei=delef(i)
      envi=envf(i)
      devi=devf(i)
      devei=devef(i)
      dili=dil(i)
      divi=div(i)
      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      ti=t(i)
      swi=s(i)
c
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
      dilci=dilc(i)
      divci=divc(i)
c gaz 020217      
c determine direction of model (define for both materials) in geneg2 and other geneq etc 
      if(gdkm_flag.eq.1) then
       if(i.le.neq_primary) then
        i_dir_gdkm = gdkm_dir(igdpm(i))
       else 
        i_dir_gdkm = gdkm_dir(igdpm(i-neq_primary))
       endif
      else
        i_dir_gdkm = -1
      endif      
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
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1
c changed by avw -- entered here by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif


c Take care of source/sink term by putting it into empty slot (from
c same node to same node).
      a_axy(jmia+nmatavw)=sk(i)
c Take care of source/sink term
      a_vxy(jmia+nmatavw)=qc(i)


c
c calculate variable for setbit
c
      iz4m1 = 4*(iz-1)+1

      do 58 jm=jmi+1,ii2
         iq=iq+1
         kb=nelm(jm)+icd
         it8(iq)=kb
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
c gaz 11-18-2001
c     if(imdnode.ne.0) then
c       imd = mdnodes(i) + mdnodes(kb)
c       if(imd.lt.2) it10(iq) = -abs(it10(iq))
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
c gaz 051616 
      allocate(grav_wgt(iq))
c
c 3-d geometry
c
      if ( icnl.eq.0 )  then
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
            reduction_factor_t = reduction_factor
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
            phikb=pvikb-pcp(kb)
c           pxy=sx2c*perml(1)+sx3c*perml(2)+sxzc*perml(3)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(i_dir_gdkm.ge.0.and.reduction_factor.gt.2) then
               kb_pri = reduction_factor -2
               reduction_factor = 1.0 
c gaz 051416 harmonic weighting in coordinate directions               
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
            if(reduction_factor.gt.2.) reduction_factor = 1.0
            pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
             if(i_dir_gdkm.ge.0.and.reduction_factor_t.gt.2) then
               kb_pri = reduction_factor_t -2
               reduction_factor_t = 1.0 
               if(i_dir_gdkm.eq.1) then
                 sx3c = sx2c*sx2t
               else if(i_dir_gdkm.eq.2) then
                 sx3c = sx2c*sx3t 
               else if(i_dir_gdkm.eq.3) then
                 sx3c = sx2c*sxzt
               else if(dis2.gt.dis_tol) then
                sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t+
     &           delz2/sxzt)
               endif                                    
            elseif(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t+
     &           delz2/sxzt)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
            if(reduction_factor_t.gt.2) reduction_factor_t = 1.0
            sx3c = reduction_factor_t*sx3c
c gaz 080118 added reduction factor to sx2c for air-water vapor diffusion       
            sx2c = reduction_factor_t*sx2c
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c                    
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
            t10(neighc)=sx2c
 59      continue
      else if(icnl.ne.0) then
c
c 2-d problem
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
            reduction_factor_t = reduction_factor
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
            dis2=delx2+dely2
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
            if(i_dir_gdkm.ge.0.and.reduction_factor_t.gt.2) then
               kb_pri = reduction_factor_t -2
               reduction_factor_t = 1.0 
               if(i_dir_gdkm.eq.1) then
                 sx3c = sx2c*sx2t
               else if(i_dir_gdkm.eq.2) then
                 sx3c = sx2c*sx3t  
               else if(dis2.gt.dis_tol) then
                 sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
               endif                                    
            elseif(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
            if(reduction_factor_t.gt.2) reduction_factor_t = 1.
             sx3c = reduction_factor_t*sx3c
             sx2c = reduction_factor_t*sx2c
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t5(neighc)=sx3c
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
            t10(neighc)=sx2c
 69      continue
      endif
c
c change coefficients if ice present
c
c gaz 10-18-2001 (set condition to unattainable value)
      if (ice.eq.-99)  then
         icesd=ices(i)
         do 49 jm=1,iq
            kb=it8(jm)
            neighc=it9(jm)
            if ( icesd.eq.2.or.ices(kb).eq.2 )  then
               t1(neighc)=0.0
               t3(neighc)=0.0
               t6(neighc)=0.0
            end if
 49      continue
      end if
c
c liquid phase calculations
c
c gaz 051616
c initialize grav_wgt
      grav_wgt = 0.5d0
      do 60 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         pxyi=t1(neighc)
         sx4d=t6(neighc)
         if(rolf(i).le.0.0.or.rolf(kb).le.0.0) grav_wgt(jm) = 1.0d0
         axyd=pxyi+grav_wgt(jm)*sx4d*(rolf(i)+rolf(kb))
     *        *(cord(kz,igrav)-cord(iz,igrav))
         t8(neighc)=axyd
 60   continue
c
c determine upwind nodes and if liquid phase exists
c
      isl=1
      do 61 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
c add coding to save upwind position
         if(iad.le.iad_up) then
            fid=0.5
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
c form equations
c
      if ( isl.ne.0 )  then
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
c
            dlpi=-pxy+grav_wgt(jm)*sx4d*dglp(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dlpkb=pxy+grav_wgt(jm)*sx4d*dglp(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
            dlei=pxy*dpvti+grav_wgt(jm)*sx4d*dgle(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+grav_wgt(jm)*sx4d*dgle(kb)
     &           *(cord(kz,igrav)-cord(iz,igrav))
            dlci=grav_wgt(jm)*sx4d*dglc(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dlckb=grav_wgt(jm)*sx4d*dglc(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
c
            aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            axyf=(fid*dilkb*(1.0-cnlkb)+fid1*dili*(1.0-cnli))
            acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
            aexy=axyd*aexyf
            axy=axyd*axyf
            acxy=axyd*acxyf
            dlapi=dlpi*axyf+axyd*fid1*(dilpi*(1.0-cnli)-dili*dcli) 
            dlapkb=dlpkb*axyf+axyd*fid*(dilpkb*(1.0-cnlkb)-dilkb*dclkb)
            dlaei=dlei*axyf+axyd*fid1*(dilei*(1.0-cnli)-dili*dclei)
            dlaekb=dlekb*axyf+axyd*fid*(dilekb*(1.0-cnlkb)-dilkb*dclekb)
            dlaci=dlci*axyf+axyd*fid1*(dilci*(1.0-cnli)-dili*dclci)
            dlackb=dlckb*axyf+axyd*fid*(dilckb*(1.0-cnlkb)-dilkb*dclckb)
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
c
            a_axy(iau+nmatavw)=axyd*(fid*dilkb+fid1*dili)
            a_axy(ial+nmatavw)=-a_axy(iau+nmatavw)

c s kelkar 3 July 2014, for calculating heat flow vectors
      if(flag_heat_out) then
         e_axy_adv(iau+nmatavw)=aexy
         e_axy_adv(ial+nmatavw)=-aexy
      endif

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
            bp(iz+nrhs(2))=bp(iz+nrhs(2))+aexy
            bp(kz+nrhs(2))=bp(kz+nrhs(2))-aexy
            a(jmia+nmat(4))=a(jmia+nmat(4))+dlepi
            a(jmia+nmat(5))=a(jmia+nmat(5))+dleei
            a(jmia+nmat(6))=a(jmia+nmat(6))+dleci
            a(ial+nmat(4))=a(ial+nmat(4))-dlepi
            a(ial+nmat(5))=a(ial+nmat(5))-dleei
            a(ial+nmat(6))=a(ial+nmat(6))-dleci
            a(iau+nmat(4))=a(iau+nmat(4))+dlepkb
            a(iau+nmat(5))=a(iau+nmat(5))+dleekb
            a(iau+nmat(6))=a(iau+nmat(6))+dleckb
            a(jml+nmat(4))=a(jml+nmat(4))-dlepkb
            a(jml+nmat(5))=a(jml+nmat(5))-dleekb
            a(jml+nmat(6))=a(jml+nmat(6))-dleckb
c
            bp(iz+nrhs(3))=bp(iz+nrhs(3))+acxy
            bp(kz+nrhs(3))=bp(kz+nrhs(3))-acxy
            a(jmia+nmat(7))=a(jmia+nmat(7))+dlcpi
            a(jmia+nmat(8))=a(jmia+nmat(8))+dlcei
            a(jmia+nmat(9))=a(jmia+nmat(9))+dlcci
            a(ial+nmat(7))=a(ial+nmat(7))-dlcpi
            a(ial+nmat(8))=a(ial+nmat(8))-dlcei
            a(ial+nmat(9))=a(ial+nmat(9))-dlcci
            a(iau+nmat(7))=a(iau+nmat(7))+dlcpkb
            a(iau+nmat(8))=a(iau+nmat(8))+dlcekb
            a(iau+nmat(9))=a(iau+nmat(9))+dlcckb
            a(jml+nmat(7))=a(jml+nmat(7))-dlcpkb
            a(jml+nmat(8))=a(jml+nmat(8))-dlcekb
            a(jml+nmat(9))=a(jml+nmat(9))-dlcckb
 62      continue
      end if
c
c vapour phase calculations
c
c gaz 051616
c initialize grav_wgt
      grav_wgt = 0.5d0
      do 63 jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         pxyh=t2(neighc)
         sx4h=t7(neighc)
         if(rovf(i).le.0.0.or.rovf(kb).le.0.0) grav_wgt(jm) = 1.0d0
         vxyd=pxyh+grav_wgt(jm)*sx4h*(rovf(i)+rovf(kb))
     &        *(cord(kz,igrav)-cord(iz,igrav))
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
c add coding to save upwind position
         if(iad.le.iad_up) then
            fid=0.5
            vxyd=t8(neighc)
            if(vxyd.lt.0.0) fid=dnwgt
            if(vxyd.gt.0.0) fid=upwgt
            if(t3(neighc).le.0.0) t9(neighc)=fid
            if(t3(neighc).gt.0.0) t9(neighc)=fid
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
c form equations
c
      if ( isl.ne.0 )  then
         do 65 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iau=it11(jm)
            ial=it12(jm)
            jml=nelmdg(kz)-neqp1
            vxyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            pxyh=t2(neighc)
            pvxy=t4(neighc)
            sx4h=t7(neighc)
c
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
c
            dvpi=-pvxy+grav_wgt(jm)*sx4h*dgvp(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dvpkb=pvxy+grav_wgt(jm)*sx4h*dgvp(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
            dvei=grav_wgt(jm)*sx4h*dgve(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dvekb=grav_wgt(jm)*sx4h*dgve(kb)
     &           *(cord(kz,igrav)-cord(iz,igrav))
            dvci=grav_wgt(jm)*sx4h*dgvc(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
            dvckb=grav_wgt(jm)*sx4h*dgvc(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
c
            vexyf=(fid*divkb*envkb+fid1*divi*envi)
            vxyf=(fid*divkb*(1.0-cnvkb)+fid1*divi*(1.0-cnvi))
            vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
            vexy=vxyd*vexyf
            vxy=vxyd*vxyf
            vcxy=vxyd*vcxyf
            dvapi=dvpi*vxyf+vxyd*fid1*(divpi*(1.0-cnvi)-divi*dcvi) 
            dvapkb=dvpkb*vxyf+vxyd*fid*(divpkb*(1.0-cnvkb)-divkb*dcvkb)
            dvaei=dvei*vxyf+vxyd*fid1*(divei*(1.0-cnvi)-divi*dcvei)
            dvaekb=dvekb*vxyf+vxyd*fid*(divekb*(1.0-cnvkb)-divkb*dcvekb)
            dvaci=dvci*vxyf+vxyd*fid1*(divci*(1.0-cnvi)-divi*dcvci)
            dvackb=dvckb*vxyf+vxyd*fid*(divckb*(1.0-cnvkb)-divkb*dcvckb)
c
c     vexyf=(fid*divkb*envkb+fid1*divi*envi)
c     vxyf=(fid*divkb+fid1*divi)
c     vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
c     vexy=vxyd*vexyf
c     vxy=vxyd*vxyf
c     vcxy=vxyd*vcxyf
c     dvapi=dvpi*vxyf+vxyd*fid1*divpi
c     dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
c     dvaei=dvei*vxyf+vxyd*fid1*divei
c     dvaekb=dvekb*vxyf+vxyd*fid*divekb
c     dvaci=dvci*vxyf+vxyd*fid1*divci
c     dvackb=dvckb*vxyf+vxyd*fid*divckb
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
c
            a_vxy(iau+nmatavw)=vxyd*(fid*divkb+fid1*divi)
            a_vxy(ial+nmatavw)=-a_vxy(iau+nmatavw)

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
            bp(iz+nrhs(2))=bp(iz+nrhs(2))+vexy
            bp(kz+nrhs(2))=bp(kz+nrhs(2))-vexy
            a(jmia+nmat(4))=a(jmia+nmat(4))+dvepi
            a(jmia+nmat(5))=a(jmia+nmat(5))+dveei
            a(jmia+nmat(6))=a(jmia+nmat(6))+dveci
            a(ial+nmat(4))=a(ial+nmat(4))-dvepi
            a(ial+nmat(5))=a(ial+nmat(5))-dveei
            a(ial+nmat(6))=a(ial+nmat(6))-dveci
            a(iau+nmat(4))=a(iau+nmat(4))+dvepkb
            a(iau+nmat(5))=a(iau+nmat(5))+dveekb
            a(iau+nmat(6))=a(iau+nmat(6))+dveckb
            a(jml+nmat(4))=a(jml+nmat(4))-dvepkb
            a(jml+nmat(5))=a(jml+nmat(5))-dveekb
            a(jml+nmat(6))=a(jml+nmat(6))-dveckb
c
            bp(iz+nrhs(3))=bp(iz+nrhs(3))+vcxy
            bp(kz+nrhs(3))=bp(kz+nrhs(3))-vcxy
            a(jmia+nmat(7))=a(jmia+nmat(7))+dvcpi
            a(jmia+nmat(8))=a(jmia+nmat(8))+dvcei
            a(jmia+nmat(9))=a(jmia+nmat(9))+dvcci
            a(ial+nmat(7))=a(ial+nmat(7))-dvcpi
            a(ial+nmat(8))=a(ial+nmat(8))-dvcei
            a(ial+nmat(9))=a(ial+nmat(9))-dvcci
            a(iau+nmat(7))=a(iau+nmat(7))+dvcpkb
            a(iau+nmat(8))=a(iau+nmat(8))+dvcekb
            a(iau+nmat(9))=a(iau+nmat(9))+dvcckb
            a(jml+nmat(7))=a(jml+nmat(7))-dvcpkb
            a(jml+nmat(8))=a(jml+nmat(8))-dvcekb
            a(jml+nmat(9))=a(jml+nmat(9))-dvcckb
 65      continue
         if(iadif.ne.0) then
c------------------------------------------------------------
c                                 MAJOR REVISIONS 3/25/05 PHS
c add air-water vapor diffusion
c
            if(s(i).lt.1.0) then
               dvai=dva(i)
               envmi=enva(i)
               devmpi=denvap(i)
               devmei=denvae(i)
               devmci=denvac(i)
               mpv_airi = rovf(i)*cnvf(i)
               mpv_wvi  = rovf(i)*(1-cnvf(i))
               do  jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  heatc=t10(neighc)
                  iau=it11(jm)
                  ial=it12(jm)
                  jml=nelmdg(kz)-neqp1
                  if(s(kb).lt.1.0) then

                     mpv_airkb = rovf(kb)*cnvf(kb)
                     mpv_wvkb  = rovf(kb)*(1-cnvf(kb))

c------  delta mass per  volume are
c-------  gives the same gradient as  PHS 3/17/04  
c-------   moles per volume in gas 

                     delmpv_air = (mpv_airkb - mpv_airi)
                     delmpv_wv  = (mpv_wvkb  - mpv_wvi)   
c
c only do for saturations lt 1.0
c
c upwinded solution         PHS need to upwind for wv and air also 
c                            separately.  
c                     but for now using wv in this eqtn
                     fid=.5
                     vxyd=t10(neighc)*delmpv_wv
                     if ( vxyd.lt.0.0) fid=dnwgt
                     if ( vxyd.gt.0.0) fid=upwgt
                     fid1=1.0-fid
c
c      Derivatives used in energy water air equations

                     dvakb=dva(kb)
                     envmkb=enva(kb)
                     dvmpkb=denvap(kb)
                     dvmekb=denvae(kb)
                     dvmckb=denvac(kb)
                     dva_bar=2.0*dvai*dvakb/(dvai+dvakb) 
                     ddva_barpi=2.0*(ddvap(i)*dvakb/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvap(i))
                     ddva_barpkb=2.0*(dvai*ddvap(kb)/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvap(kb))
                     ddva_barei=2.0*(ddvae(i)*dvakb/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvae(i))
                     ddva_barekb=2.0*(dvai*ddvae(kb)/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvae(kb))
                     ddva_barci=2.0*(ddvac(i)*dvakb/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvac(i))
                     ddva_barckb=2.0*(dvai*ddvac(kb)/(dvai+dvakb)
     &                    - dvai*dvakb/(dvai+dvakb)**2*ddvac(kb))
                     dvame=(fid*envmkb+fid1*envmi)
                     dvamet=dvame*dva_bar

c---------------------------------------------------
c   PHS   3/24/2004
c   Adding a_wvxy to store mass flux calculations 
c   for use in coneq1 
c----------------------------------------------------

                     a_wvxy(iau+nmatavw)=+heatc*dva_bar*delmpv_wv
                     a_wvxy(ial+nmatavw)=-heatc*dva_bar*delmpv_wv

c------------------------------------------------------------------
c energy equation          PHS  3/23/04 WV for energy 
c--------0---------2---------3---------4---------5---------6---------7
c                 PHS 3/19/04  Changing derivative in the a(blah)
c------------------------------------------------------------------

                     bp(iz+nrhs(2))=bp(iz+nrhs(2))+
     &                    heatc*dvamet*delmpv_wv
                     bp(kz+nrhs(2))=bp(kz+nrhs(2))-
     &                    heatc*dvamet*delmpv_wv

                     a(jmia+nmat(4))=a(jmia+nmat(4))+
     &                    heatc*((devmpi*fid1*dva_bar
     &                    +dvame*ddva_barpi)*delmpv_wv
     &                    +dvamet*(-dgvp(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvf(i) ) )
                     
                     a(jmia+nmat(5))=a(jmia+nmat(5))+
     &                    heatc*((devmei*fid1*dva_bar
     &                    +dvame*ddva_barei)*delmpv_wv
     &                    +dvamet*(-dgve(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvef(i) ) )
                     
                     a(jmia+nmat(6))=a(jmia+nmat(6))+
     &                    heatc*((devmci*fid1*dva_bar
     &                    +dvame*ddva_barci)*delmpv_wv
     &                    +dvamet*(-dgvc(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvcf(i) ) )
                     
c---- 
                     a(ial+nmat(4))=a(ial+nmat(4))-
     &                    heatc*((devmpi*fid1*dva_bar
     &                    +dvame*ddva_barpi)*delmpv_wv
     &                    +dvamet*(-dgvp(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvf(i) ) )

                     a(ial+nmat(5))=a(ial+nmat(5))-
     &                    heatc*((devmei*fid1*dva_bar
     &                    +dvame*ddva_barei)*delmpv_wv
     &                    +dvamet*(-dgve(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvef(i) ) )

                     a(ial+nmat(6))=a(ial+nmat(6))-
     &                    heatc*((devmci*fid1*dva_bar
     &                    +dvame*ddva_barci)*delmpv_wv
     &                    +dvamet*(-dgvc(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvcf(i) ) )

c=========================================

                     a(iau+nmat(4))=a(iau+nmat(4))+
     &                    heatc*((dvmpkb*fid*dva_bar
     &                    +dvame*ddva_barpkb)*delmpv_wv
     &                    +dvamet*(dgvp(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvf(kb) ) )
                     
                     a(iau+nmat(5))=a(iau+nmat(5))+
     &                    heatc*((dvmekb*fid*dva_bar
     &                    +dvame*ddva_barekb)*delmpv_wv
     &                    +dvamet*(dgve(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvef(kb) ) )

                     a(iau+nmat(6))=a(iau+nmat(6))+
     &                    heatc*((dvmckb*fid*dva_bar
     &                    +dvame*ddva_barckb)*delmpv_wv
     &                    +dvamet*(dgvc(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvcf(kb) ) )

c----
                     a(jml+nmat(4))=a(jml+nmat(4))-
     &                    heatc*((dvmpkb*fid*dva_bar
     &                    +dvame*ddva_barpkb)*delmpv_wv
     &                    +dvamet*(dgvp(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvf(kb) ) )
                     
                     a(jml+nmat(5))=a(jml+nmat(5))-
     &                    heatc*((dvmekb*fid*dva_bar
     &                    +dvame*ddva_barekb)*delmpv_wv
     &                    +dvamet*(dgve(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvef(kb) ) )
                     
                     a(jml+nmat(6))=a(jml+nmat(6))-
     &                    heatc*((dvmckb*fid*dva_bar
     &                    +dvame*ddva_barckb)*delmpv_wv
     &                    +dvamet*(dgvc(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvcf(kb) ) )

c-----------------------------------
c    PHS  21 AUG  2003
c-----------------------------------
c air equation
c        bp(iz+nrhs(3))=bp(iz+nrhs(3))+heatc*dva_bar*delmpv_air
c        bp(kz+nrhs(3))=bp(kz+nrhs(3))-heatc*dva_bar*delmpv_air
c        a(jmia+nmat(7))=a(jmia+nmat(7))+heatc*(ddva_barpi
c    &   *delmpv_air-dva_bar*dcvf(i))
c        a(jmia+nmat(8))=a(jmia+nmat(8))+heatc*(ddva_barei
c    &   *delmpv_air-dva_bar*dcvef(i))
c        a(jmia+nmat(9))=a(jmia+nmat(9))+heatc*(ddva_barci
c    &   *delmpv_air-dva_bar*dcvcf(i))
c        a(ial+nmat(7))=a(ial+nmat(7))-heatc*(ddva_barpi
c    &   *delmpv_air-dva_bar*dcvf(i))
c        a(ial+nmat(8))=a(ial+nmat(8))-heatc*(ddva_barei
c    &   *delmpv_air-dva_bar*dcvef(i))
c        a(ial+nmat(9))=a(ial+nmat(9))-heatc*(ddva_barci
c    &   *delmpv_air-dva_bar*dcvcf(i))
c        a(iau+nmat(7))=a(iau+nmat(7))+heatc*(ddva_barpkb
c    &   *delmpv_air+dva_bar*dcvf(kb))
c        a(iau+nmat(8))=a(iau+nmat(8))+heatc*(ddva_barekb
c    &   *delmpv_air+dva_bar*dcvef(kb))
c        a(iau+nmat(9))=a(iau+nmat(9))+heatc*(ddva_barckb
c    &   *delmpv_air+dva_bar*dcvcf(kb))
c        a(jml+nmat(7))=a(jml+nmat(7))-heatc*(ddva_barpkb
c    &   *delmpv_air+dva_bar*dcvf(kb))
c        a(jml+nmat(8))=a(jml+nmat(8))-heatc*(ddva_barekb
c    &   *delmpv_air+dva_bar*dcvef(kb))
c        a(jml+nmat(9))=a(jml+nmat(9))-heatc*(ddva_barckb
c    &   *delmpv_air+dva_bar*dcvcf(kb))
c
c-------------------------------------------------------------------
c water equation  PHS 3/11/04  Now a function of delta mol/m3 wv
c                  which means flipping the signs of +-heatc
c                  to match the air equations. . . 
c--------0---------2---------3---------4---------5---------6---------7
c                 PHS 3/17/04  Changing derivative in the a(blah)
c
                     bp(iz+nrhs(1))=bp(iz+nrhs(1))+
     &                    heatc*dva_bar*delmpv_wv
                     bp(kz+nrhs(1))=bp(kz+nrhs(1))-
     &                    heatc*dva_bar*delmpv_wv

c========================

                     a(jmia+nmat(1))=a(jmia+nmat(1))+
     &                    heatc*(ddva_barpi*delmpv_wv
     &                    + dva_bar*(-dgvp(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvf(i) )) 

                     a(jmia+nmat(2))=a(jmia+nmat(2))+
     &                    heatc*(ddva_barei*delmpv_wv
     &                    + dva_bar*(-dgve(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvef(i) ))

                     a(jmia+nmat(3))=a(jmia+nmat(3))+
     &                    heatc*(ddva_barci*delmpv_wv
     &                    + dva_bar*(-dgvc(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvcf(i) ))

c-------

                     a(ial+nmat(1))=a(ial+nmat(1))-
     &                    heatc*(ddva_barpi*delmpv_wv
     &                    + dva_bar*(-dgvp(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvf(i) ))

                     a(ial+nmat(2))=a(ial+nmat(2))-
     &                    heatc*(ddva_barei*delmpv_wv
     &                    + dva_bar*(-dgve(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvef(i) ))

                     a(ial+nmat(3))=a(ial+nmat(3))-
     &                    heatc*(ddva_barci*delmpv_wv
     &                    + dva_bar*(-dgvc(i)*(1-cnvf(i))+
     &                    rovf(i)*dcvcf(i) ))

c=======================

                     a(iau+nmat(1))=a(iau+nmat(1))+
     &                    heatc*(ddva_barpkb*delmpv_wv
     &                    + dva_bar*(dgvp(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvf(kb) ))

                     a(iau+nmat(2))=a(iau+nmat(2))+
     &                    heatc*(ddva_barekb*delmpv_wv
     &                    + dva_bar*(dgve(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvef(kb) ))

                     a(iau+nmat(3))=a(iau+nmat(3))+
     &                    heatc*(ddva_barckb*delmpv_wv
     &                    + dva_bar*(dgvc(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvcf(kb) ))

c-------

                     a(jml+nmat(1))=a(jml+nmat(1))-
     &                    heatc*(ddva_barpkb*delmpv_wv
     &                    + dva_bar*(dgvp(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvf(kb) ))

                     a(jml+nmat(2))=a(jml+nmat(2))-
     &                    heatc*(ddva_barekb*delmpv_wv
     &                    + dva_bar*(dgve(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvef(kb) ))

                     a(jml+nmat(3))=a(jml+nmat(3))-
     &                    heatc*(ddva_barckb*delmpv_wv
     &                    + dva_bar*(dgvc(kb)*(1-cnvf(kb))-
     &                    rovf(kb)*dcvcf(kb) ))

                  endif
               enddo
            endif
         endif
      endif
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
      
c s kelkar 3 July 2014, for calculating heat flow vectors
      heatt = heatc*(t(kb)-ti)
      if(flag_heat_out) then
         e_axy_cond(iau+nmatavw)=+heatt
         e_axy_cond(ial+nmatavw)=-heatt
      endif

         bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatt
         bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatt
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
 66   continue

      skz = 0.
      qhz = 0.
      dskzt1 = 0.
      dqhzt1 = 0.
      dskzt2 = 0.
      dqhzt2 = 0.

      if( izeolites .ne. 0) then
         if ( kzeol(i) .gt. 0. ) then
            if(ieos(i) .ne. 3) then
               water_pressure = psat(t(i), dpdt, 0)
            else
               water_pressure = psat(t(i), dpdt, 0)
            end if
            call zeolites(t(i), kzeol(i), fwater_old(i), fwater(i), 
     2           water_pressure, dpdt, skz, dskzt, qhz, dqhzt)
            skz = skz/dtot
            qhz = qhz/dtot
            if( ieos(i) .eq. 2 ) then
               dskzt1 = 0.
               dqhzt1 = 0.
               dskzt2 = dskzt/dtot
               dqhzt2 = dqhzt/dtot
            else
               dskzt2 = 0.
               dqhzt2 = 0.
               dskzt1 = dskzt/dtot
               dqhzt1 = dqhzt/dtot
            end if
         end if
      end if

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)+ skz
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)+ qhz
      bp(iz+nrhs(3))=bp(iz+nrhs(3))+sx1d*denpci(i)+qc(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i) + dskzt1
c gaz debug 121104 now have a derivative wrt air pressure;see thrmwc
      a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*dmc(i) +dqpc(i) + dskzt2
      a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(5))=a(jmia+nmat(5))+sx1d*deef(i)+deqh(i) + dqhzt1
      a(jmia+nmat(6))=a(jmia+nmat(6))+sx1d*dec(i)+dcqh(i) + dqhzt2
      a(jmia+nmat(7))=a(jmia+nmat(7))+sx1d*dcp(i)+dqc(i)
      a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*dce(i)+deqc(i)
      a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*dcc(i)+dcqc(i)

c gaz 051616
      deallocate(grav_wgt)

      return
      end
