      subroutine  geneq1  ( i )
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
CD1  and mass transfer without non-condensible gas.  
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
CD2 $Log:   /pvcs.config/fehm90/src/geneq1.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.12   Mon Mar 31 08:36:02 1997   gaz
CD2 minor changes for anisotropic properties
CD2 
CD2    Rev 1.11   Thu Feb 15 10:43:04 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.10   Wed Feb 07 10:56:36 1996   gaz
CD2 replaced if(nmat(1).eq.nmat(11)) then with if(i.gt.neq) then
CD2 
CD2    Rev 1.9   Wed Jan 10 13:26:54 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.8   12/13/95 08:37:58   gaz
CD2 changed setbit counted to accomodate nbits=256
CD2 
CD2    Rev 1.7   08/07/95 11:32:54   awolf
CD2 Sets a_axy and a_vxy terms for frac - matrix flow
CD2 
CD2    Rev 1.6   06/21/95 11:08:18   llt
CD2 declared bit function to be logical type for IBM
CD2 
CD2    Rev 1.5   06/01/95 16:51:08   gaz
CD2 upwind bitmap
CD2 
CD2    Rev 1.4   04/03/95 08:45:22   robinson
CD2 Minor change to mass flow rate array
CD2 
CD2    Rev 1.3   01/28/95 14:03:38   llt
CD2 modified for new particle tracking module
CD2 
CD2    Rev 1.2   05/11/94 16:10:40   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 15:51:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:08   pvcs
CD2 original version in process of being certified
CD2 
c version FEHM5.J changes
c 15-july-91
c in geneq1 checked for isl.ne.0 for vap phase calcs
c 16-july-91
c zerod out dfmp(i) etc in here rather than in optsj.f
c 17-july-91
c corrected dlei and dlekb terms in gensl1
c 5-august-91
c changed gravity term  using dx/dx
c 6-august-91
c more mods to gravity term and made changes fo co2 as well
c 8-august-91
c corrections made to gravity term
c 12-august-91
c corrected derivatives in geneqc(added divekb)
c 19-august-91
c removed reference to rlf and rvf and derivatives,also gl and gv
c 11-oct-91
c changed geneq1 so i>neq is calculated correctly(dpdp solution)
c 12-oct-91
c put i=iz in radi calc in geneq1
c 9-nov-91
c putting new sx system in
c 11-nov-91
c big changes for symmetry in assembly
c making like geneq2 in airsj.f
c 12-nov-91
c symmetry changes to geneqc
c got rid of fad,fbd,fcd in call sequence
c 13-nov-91
c modified geneq3 to put in symmetry
c complete rework
c changed geneq1 to simplify 2-d/3-d calcs
c 15-nov-91
c added coding to define icd
c 22-nov-91
c corrected if statement "if(nelm(ij).eq.i)"
c corrected "nelmdg(kb)"
c 23-nov-91
c more corrections  in geneq1:"bp(iz..."
c 31-mar-92
c changed cdir$ ivdep to cdir$ novector around line 703
c now runs on unicos,I think this is a compiler bug
c 19-may-92
c t3(neighc).lt to t3(neighc).le,also t4
c 26-may-92
c made t9=fid always
c 30-july-92
c in all routines set t6()=t3),t7()=t4()
c 15-jan-93
c added air-vapor diffusion 
c 19-jan-93
c stil working on air-vapor diffusion
c main coding routine geneqc
c 2-feb-93
c skipped air-water vapor diffusion if iadif=0
c 23-feb-93 llt
c tok out equivalences to put in memory management
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
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
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
      use comriv
      use comfem, only: permfactor
      use comsi, only: ihms
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
      integer imd,iwd
      integer edge
      integer i_pri
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor, reduction_factor_t
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
      real*8 heatt
c gaz 110822 debug
      real*8 red_tmp(20), sx_tmp(20)
      integer kb_pri,i_dir_gdkm
      integer igdkm_test
      parameter(dis_tol=1.d-12, igdkm_test = 0) 
c
c
! Guessed at value for grav_air
c      grav_air = 0.
       grav_air = grav
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
c gaz 020217      
c determine direction of model (define for both materials) in geneg2 and other geneq etc 
      if(gdkm_flag.eq.1) then
       if(i.le.neq_primary) then
        i_dir_gdkm = gdkm_dir(igdpm(i))
       else 
        i_pri = nelm(nelm(i)+1)
c gaz debug 113022 i_pri not i-neq_primary        
c        i_dir_gdkm = gdkm_dir(igdpm(i-neq_primary))
        i_dir_gdkm = gdkm_dir(igdpm(i_pri))
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
      a_axy(jmia+nmatavw)=sk(i)*s(i)
c Take care of souce/sink term
c If this is an isothermal air-water simulation
      a_vxy(jmia+nmatavw)=sk(i)*(1-s(i))

      do 58 jm=jmi+1,ii2
      kb=nelm(jm)+icd
      if(iriver.eq.2.and.kb.gt.neq_primary) go to 58
      iq=iq+1
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
68    continue
58    continue
c gaz 051616 
      allocate(grav_wgt(iq))
      if(icnl.eq.0) then
c
c 3-d geometry
c
         do 59 jm=1,iq
            kb=it8(jm)
            edge = it11(jm) + neqp1
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
            reduction_factor_t = reduction_factor 
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)

            if(allocated(permfactor)) then
              perml(1) = perml(1)*permfactor(edge,1)
              perml(2) = perml(2)*permfactor(edge,2)
              perml(3) = perml(3)*permfactor(edge,3)
            endif

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
            if(iriver.eq.2.and.kb.gt.neq_primary) then
c this connection broken, then added when well equations are generated            
              sx2c = 0.0
            endif
            dis2=delx2+dely2+delz2
            if(i_dir_gdkm.ge.0.and.reduction_factor.gt.2) then
               kb_pri = reduction_factor -2
               reduction_factor = 1.0 
c use directional harmonic weighting to match high res solution               
               if(i_dir_gdkm.eq.1) then
c                 pxy = sx2c*pnx(kb_pri)
                  pxy = sx2c*perml(1)
               else if(i_dir_gdkm.eq.2) then
c                 pxy = sx2c*pny(kb_pri) 
                  pxy = sx2c*perml(2) 
               else if(i_dir_gdkm.eq.3) then
c                 pxy = sx2c*pnz(kb_pri)
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
c gaz 110822            
            red_tmp(jm) = reduction_factor
            reduction_factor_t = reduction_factor
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
c  gaz 110822 
            sx_tmp(jm) = sx2c
            thxkb=thx(kb)
            thykb=thy(kb)
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(i_dir_gdkm.ge.0.and.reduction_factor.gt.2) then
               kb_pri = reduction_factor -2
               reduction_factor = 1.0 
c gaz 050118 harmonic weightging to match high resolution grid               
               if(i_dir_gdkm.eq.1) then
c                 pxy = sx2c*pnx(kb_pri)
                  pxy = sx2c*perml(1)
               else if(i_dir_gdkm.eq.2) then
c                 pxy = sx2c*pnx(kb_pri)
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
               reduction_factor = 1.0 
               if(i_dir_gdkm.eq.1) then
                 sx3c = sx2c*sx2t
               else if(i_dir_gdkm.eq.2) then
                 sx3c = sx2c*sx3t  
               else if(dis2.gt.dis_tol) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
               endif                                    
            elseif(dis2.gt.dis_tol.and.iwd.gt.0) then
c gaz 092922 had extra 3d term in sx3c (removed it)                
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
              if(reduction_factor_t.gt.2) reduction_factor_t = 1.
              sx3c = reduction_factor_t*sx3c
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
c gaz 110822 debug print out gdkm connections
      if(igdkm_test.ne.0) then
      if(i.eq.1) then
       write(ierr,399) 
       write(ierr,400) 
      endif
      if(i.le.neq_primary) then
        if(gdkm_dir(igdpm(i)).gt.0) then
         write(ierr,401) i,(it8(jm),jm = 1,iq)
         write(ierr,402) 'redfac',(red_tmp(jm),jm = 1,iq)
         write(ierr,403) gdkm_dir(igdpm(i)),(sx_tmp(jm),jm = 1,iq)
        endif
      else if(i.gt.neq_primary) then
        if(gdkm_dir(igdpm(nelm(nelm(i)+1))).gt.0) then
         write(ierr,401) i,(it8(jm),jm = 1,iq)
         write(ierr,402) 'redfac',(red_tmp(jm),jm = 1,iq)
         write(ierr,403) gdkm_dir(igdpm(nelm(nelm(i)+1))),
     %                  (sx_tmp(jm),jm = 1,iq)
        endif
      endif
399   format(/,'>>>>> Check Eq connections with GDKM   <<<<<<',/)
400   format(t1,'node',t15,'n1',t30,'n2',t45,'n3',t60,'n4',t75,'n5')
401   format(t1,i6,t15,i6,t30,i6,t45,i6,t60,i6,t75,i6)  
402   format(t1,a6,t15,g14.5,t30,g14.5,t45,g14.5,t60,g14.5,
     &       t75,g14.5) 
403   format(t1,'f dir ',i6,t15,g14.5,t30,g14.5,t45,g14.5,t60,g14.5,
     &       t75,g14.5)       
404   format('>>>>> End Check Eq connections with GDKM   <<<<<<',/)      
      if(i.eq.neq) then
       write(ierr,404) 
      endif
      endif
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
      enlkb=enlf(kb)
      delkb=delf(kb)
      delekb=delef(kb)
      dilekb=dile(kb)
      dlpi=-pxy+grav_wgt(jm)*sx4d*dglp(i)*
     &     (cord(kz,igrav)-cord(iz,igrav))
      dlpkb=pxy+grav_wgt(jm)*sx4d*dglp(kb)*
     &     (cord(kz,igrav)-cord(iz,igrav))
      dlei=pxy*dpvti+grav_wgt(jm)*sx4d*dgle(i)*
     &     (cord(kz,igrav)-cord(iz,igrav))
      dlekb=-pxy*dpcef(kb)+grav_wgt(jm)*sx4d*dgle(kb)*
     &     (cord(kz,igrav)-cord(iz,igrav))
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
c
c changed by avw -- entered here by seh
      a_axy(iau+nmatavw)=axy
      a_axy(ial+nmatavw)=-axy

c s kelkar 3 July 2014, for calculating heat flow vectors
      if(flag_heat_out) then
         e_axy_adv(iau+nmatavw)=aexy
         e_axy_adv(ial+nmatavw)=-aexy
      endif

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
      a(jmia+nmat(3))=a(jmia+nmat(3))+dlepi
      a(jmia+nmat(4))=a(jmia+nmat(4))+dleei
      a(ial+nmat(3))=a(ial+nmat(3))-dlepi
      a(ial+nmat(4))=a(ial+nmat(4))-dleei
      a(iau+nmat(3))=a(iau+nmat(3))+dlepkb
      a(iau+nmat(4))=a(iau+nmat(4))+dleekb
      a(jml+nmat(3))=a(jml+nmat(3))-dlepkb
      a(jml+nmat(4))=a(jml+nmat(4))-dleekb
   62 continue
      e n d i f
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
      envkb=envf(kb)
      devkb=devf(kb)
      devekb=devef(kb)
      divekb=dive(kb)
      dvpi=-pvxy+grav_wgt(jm)*sx4h*dgvp(i)*
     &     (cord(kz,igrav)-cord(iz,igrav))
      dvpkb=pvxy+grav_wgt(jm)*sx4h*dgvp(kb)*
     &     (cord(kz,igrav)-cord(iz,igrav))
      dvei=grav_wgt(jm)*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
      dvekb=grav_wgt(jm)*sx4h*dgve(kb)*
     &      (cord(kz,igrav)-cord(iz,igrav))
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
      a(jmia+nmat(3))=a(jmia+nmat(3))+dvepi
      a(jmia+nmat(4))=a(jmia+nmat(4))+dveei
      a(ial+nmat(3))=a(ial+nmat(3))-dvepi
      a(ial+nmat(4))=a(ial+nmat(4))-dveei
      a(iau+nmat(3))=a(iau+nmat(3))+dvepkb
      a(iau+nmat(4))=a(iau+nmat(4))+dveekb
      a(jml+nmat(3))=a(jml+nmat(3))-dvepkb
      a(jml+nmat(4))=a(jml+nmat(4))-dveekb
   65 continue
      e n d i f
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
      
      heatt = +heatc*(t(kb)-ti)
c s kelkar 3 July 2014, for calculating heat flow vectors
      if(flag_heat_out) then
         e_axy_cond(iau+nmatavw)=+heatt
         e_axy_cond(ial+nmatavw)=-heatt
      endif

      bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatt
      bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatt
      a(jmia+nmat(3))=a(jmia+nmat(3))-heatc*dtpa(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))-heatc*dtpae(i)
      a(ial+nmat(3))=a(ial+nmat(3))+heatc*dtpa(i)
      a(ial+nmat(4))=a(ial+nmat(4))+heatc*dtpae(i)
      a(iau+nmat(3))=a(iau+nmat(3))+heatc*dtpa(kb)
      a(iau+nmat(4))=a(iau+nmat(4))+heatc*dtpae(kb)
      a(jml+nmat(3))=a(jml+nmat(3))-heatc*dtpa(kb)
      a(jml+nmat(4))=a(jml+nmat(4))-heatc*dtpae(kb)
66    continue
c
c     bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*(aw*deni(i)+ay*denj(i))+sk(i)
c     bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*(aw*denei(i)+ay*denej(i))+qh(i)
c     a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*(aw*dmpf(i))+dq(i)
c     a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*(aw*dmef(i))+dqt(i)
c     a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*(depf(i)*aw)+dqh(i)
c     a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*(aw*deef(i))+deqh(i)
c gaz 092922 debug   
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i)
      a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*deef(i)+deqh(i)
c gaz 051616
      deallocate(grav_wgt)
      r e t u r n
      e    n    d
