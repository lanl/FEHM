      subroutine  geneq1_stress_coupl  ( i )
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
	use comsi
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
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

	real*8 rlsxij, rlsyij, rlszij, dpxyp1, dpxyp2, dpxyp3 
	integer kru1,kru2,krv1,krv2,krw1,krw2 
	integer iu1,iu2,iv1,iv2,iw1,iw2,k,kb1
      integer jm1,kbu1,kbu2,kbv1,kbv2,kbw1,kbw2

      real*8 diffp,dmep1,dmep2,dmep3
	real*8 flpart,flpartm,flparte
      real*8 deep1,deep2,deep3,alxavg,alyavg,alzavg
      real*8 dxiv1,dxiv2,dxiw1,dxiw2
      real*8 dxkbv1,dxkbv2,dxkbw1,dxkbw2
      real*8 dyiu1,dyiu2,dyiw1,dyiw2
      real*8 dykbu1,dykbu2,dykbw1,dykbw2   
      real*8 dziu1,dziu2,dziv1,dziv2
      real*8 dzkbu1,dzkbu2,dzkbv1,dzkbv2    
      integer idxiv1,idxiv2,idxiw1,idxiw2
      integer idxkbv1,idxkbv2,idxkbw1,idxkbw2
      integer idyiu1,idyiu2,idyiw1,idyiw2
      integer idykbu1,idykbu2,idykbw1,idykbw2   
      integer idziu1,idziu2,idziv1,idziv2
      integer idzkbu1,idzkbu2,idzkbv1,idzkbv2 

      parameter(dis_tol=1.d-12)
c
c
! Guessed at value for grav_air
c
      grav_air = 0.
c changed by avw 4/95 -- entered into new version by seh
      neqp1=neq+1
      ldna=nelm(neqp1)-neqp1
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

      do  jm=ii1,ii2
       kb=nelm(jm)+icd
	 if(kb.ne.i) then	 
        iq=iq+1
        it8(iq)=kb
        it9(iq)=jm-ii1+1
        it10(iq)=istrw(jm-neqp1)
        it11(iq)=jm-neqp1
       endif
      enddo
c
c check for inclusion in stencil of perm derivitives wrt displacements
c (zero out fixed displacements) set position to ii (jmia)
c determine position in a(jacobian) array
c
c  
      do jm1=1,iq
	  jm = it11(jm1)
	  kb1 = it8(jm1)
        kbu1 = ipermx(kb1,1)
	  kbu2 = ipermx(kb1,2)
	  kbv1 = ipermy(kb1,1)
	  kbv2 = ipermy(kb1,2)
	  kbw1 = ipermz(kb1,1)
	  kbw2 = ipermz(kb1,2)
	  kru1 = kr(kbu1,1)
	  kru2 = kr(kbu2,1)
	  krv1 = kr(kbv1,2)
	  krv2 = kr(kbv2,2)	  
	  krw1 = kr(kbw1,3)
	  krw2 = kr(kbw2,3)	 
	  if(kbu1.eq.0) then
	   kru1 = 1
	   kbu1 = i
	  endif
	  if(kbu2.eq.0) then
	   kru2 = 1
	   kbu2 = i
	  endif
	  if(kbv1.eq.0) then
	   krv1 = 2
	   kbv1 = i
	  endif
	  if(kbv2.eq.0) then
	   krv2 = 2
	   kbv2 = i
	  endif 
	  if(kbw1.eq.0) then
	   krw1 = 3
	   kbw1 = i
	  endif
	  if(kbw2.eq.0) then
	   krw2 = 3
	   kbw2 = i
	  endif 	  
	   its32(jm1,1) = 0
	   its32(jm1,2) = 0
	   its42(jm1,1) = 0
	   its42(jm1,2) = 0
	   its22(jm1,1) = 0
	   its22(jm1,2) = 0
	   its32(jm1,3) = 0
	   its32(jm1,4) = 0
	   its42(jm1,3) = 0
	   its42(jm1,4) = 0
	   its22(jm1,3) = 0
	   its22(jm1,4) = 0
	   do k=ii1,ii2     
	     kb= nelm(k)
	     if(kbu1.eq.kb) then
	      if(kru1.eq.1) then
	       its32(jm1,1) = jmia
	       its42(jm1,1) = jmia
	       ts32(jm1,1) = 0.0
	       ts42(jm1,1) = 0.0
	      else
	       its32(jm1,1) = k-neqp1
	       its42(jm1,1) = k-neqp1
	       ts32(jm1,1) = 0.5*drlys(kb,1)
	       ts42(jm1,1) = 0.5*drlzs(kb,1)
	      endif
	     endif
	     if(kbu2.eq.kb) then
	      if(kru2.eq.1) then	     
	       its32(jm1,2) = jmia
	       its42(jm1,2) = jmia
	       ts32(jm1,2) = 0.0
	       ts42(jm1,2) = 0.0
	      else
	       its32(jm1,2) = k-neqp1
	       its42(jm1,2) = k-neqp1	
	       ts32(jm1,2) = 0.5*drlys(kb,2)
	       ts42(jm1,2) = 0.5*drlzs(kb,2)
	      endif
	     endif
	     if(kbv1.eq.kb) then
	      if(krv1.eq.2) then
	       its22(jm1,1) = jmia
	       its42(jm1,3) = jmia
	       ts22(jm1,1) = 0.0
	       ts42(jm1,3) = 0.0
	      else
	       its22(jm1,1) = k-neqp1
	       its42(jm1,3) = k-neqp1
	       ts22(jm1,1) = 0.5*drlxs(kb,1)
	       ts42(jm1,3) = 0.5*drlzs(kb,3)
	      endif
	     endif
	     if(kbv2.eq.kb) then
	      if(krv2.eq.2) then
	       its22(jm1,2) = jmia
	       its42(jm1,4) = jmia
	       ts22(jm1,2) = 0.0
	       ts42(jm1,4) = 0.0
	      else
	       its22(jm1,2) = k-neqp1
	       its42(jm1,4) = k-neqp1	
	       ts22(jm1,2) = 0.5*drlxs(kb,2)
	       ts42(jm1,4) = 0.5*drlzs(kb,4)
	      endif
	     endif
	     if(kbw1.eq.kb) then
	      if(krw1.eq.3) then	     
	       its22(jm1,3) = jmia
	       its32(jm1,3) = jmia
	       ts22(jm1,3) = 0.0
	       ts32(jm1,3) = 0.0
	      else
	       its22(jm1,3) = k-neqp1
	       its32(jm1,3) = k-neqp1
	       ts22(jm1,3) = 0.5*drlxs(kb,3)
	       ts32(jm1,3) = 0.5*drlys(kb,3)
	      endif
	     endif
	     if(kbw2.eq.kb) then
	      if(krw2.eq.3) then		     
	       its22(jm1,4) = jmia
	       its32(jm1,4) = jmia
	       ts22(jm1,4) = 0.0
	       ts32(jm1,4) = 0.0
	      else
	       its22(jm1,4) = k-neqp1
	       its32(jm1,4) = k-neqp1
	       ts22(jm1,4) = 0.5*drlxs(kb,4)
	       ts32(jm1,4) = 0.5*drlys(kb,4)
	      endif
	     endif
	   enddo
	enddo
	
c
c now do i (zero out fixed displacements)
c
            iu1 = ipermx(i,1)
	    iu2 = ipermx(i,2)
	    iv1 = ipermy(i,1)
            iv2 = ipermy(i,2)
            iw1 = ipermz(i,1)
            iw2 = ipermz(i,2) 
	    kru1 = kr(i,1)
	    kru2 = kr(i,1)
	    krv1 = kr(i,2)
	    krv2 = kr(i,2)	  
	    krw1 = kr(i,3)
	    krw2 = kr(i,3)
	    if(iu1.eq.0) then
	     kru1 = 1
	     iu1 = i
	    endif
	    if(iu2.eq.0) then
	     kru2 = 1
	     iu2 = i
	    endif
	    if(iv1.eq.0) then
	     krv1 = 2
	     iv1 = i
	    endif
	    if(iv2.eq.0) then
	     krv2 = 2
	     iv2 = i
	    endif 
	    if(iw1.eq.0) then
	     krw1 = 3
	     iw1 = i
	    endif
	    if(iw2.eq.0) then
	     krw2 = 3
	     iw2 = i
	    endif 	  	        		   
	   its31(1,1) = 0
	   its31(1,2) = 0
	   its41(1,1) = 0
	   its41(1,2) = 0
	   its21(1,1) = 0
	   its21(1,2) = 0
	   its31(1,3) = 0
	   its31(1,4) = 0
	   its41(1,3) = 0
	   its41(1,4) = 0
	   its21(1,3) = 0
	   its21(1,4) = 0	 
	   do k=ii1,ii2
	     jm = k-ii1+1
	     kb= nelm(k)
	     if(iu1.eq.kb) then
	      if(kru1.eq.1) then
	       its31(1,1) = jmia
	       its41(1,1) = jmia
	       ts31(1,1) = 0.0
	       ts41(1,1) = 0.0
	      else
	       its31(1,1) = k-neqp1
	       its41(1,1) = k-neqp1
	       ts31(1,1) = 0.5*drlys(i,1)
	       ts41(1,1) = 0.5*drlzs(i,1)
	      endif
	     endif
	     if(iu2.eq.kb) then
	      if(kru2.eq.1) then
	       its31(1,2) = jmia
	       its41(1,2) = jmia
	       ts31(1,2) = 0.0
	       ts41(1,2) = 0.0
	      else
	       its31(1,2) = k-neqp1
	       its41(1,2) = k-neqp1
	       ts31(1,2) = 0.5*drlys(i,2)
	       ts41(1,2) = 0.5*drlzs(i,2)
	      endif
	     endif
	     if(iv1.eq.kb) then
	      if(krv1.eq.2) then	     
	       its21(1,1) = jmia
	       its41(1,3) = jmia
	       ts21(1,1) = 0.0	 
	       ts41(1,3) = 0.0
	      else
	       its21(1,1) = k-neqp1
	       its41(1,3) = k-neqp1
	       ts21(1,1) = 0.5*drlxs(i,1)
	       ts41(1,3) = 0.5*drlzs(i,3)
	      endif
	     endif
	     if(iv2.eq.kb) then
	      if(krv2.eq.2) then	     
	       its21(1,2) = jmia
	       its41(1,4) = jmia
	       ts21(1,2) = 0.0	 
	       ts41(1,4) = 0.0
	      else
	       its21(1,2) = k-neqp1
	       its41(1,4) = k-neqp1
	       ts21(1,2) = 0.5*drlxs(i,2)
	       ts41(1,4) = 0.5*drlzs(i,4)
	      endif
	     endif
	     if(iw1.eq.kb) then
	      if(krw1.eq.3) then
	       its21(1,3) = jmia
	       its31(1,3) = jmia
	       ts21(1,3) = 0.0
	       ts31(1,3) = 0.0
	      else
	       its21(1,3) = k-neqp1
	       its31(1,3) = k-neqp1
	       ts21(1,3) = 0.5*drlxs(i,3)
	       ts31(1,3) = 0.5*drlys(i,3)
	      endif
	     endif
	     if(iw2.eq.kb) then
	      if(krw2.eq.3) then	     
	       its21(1,4) = jmia
	       its31(1,4) = jmia
	       ts21(1,4) = 0.0
	       ts31(1,4) = 0.0
	      else
	       its21(1,4) = k-neqp1
	       its31(1,4) = k-neqp1
	       ts21(1,4) = 0.5*drlxs(i,4)
	       ts31(1,4) = 0.5*drlys(i,4)
	      endif
	     endif
	   enddo

c
c now zero out derivatives with fixed displacements
c

c  might have to forgo symmetric assembly 

      if(icnl.eq.0) then
c
c 3-d geometry
c
c iq is now all the neighbors of i
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
            perml(1)=2.*reduction_factor*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*reduction_factor*alykb*alyi/(alykb+alyi)
            perml(3)=2.*reduction_factor*alzkb*alzi/(alzkb+alzi)
c
c  new code for perm changes with stress 
c
	      rlsxij = 0.5*(rlxs(i) + rlxs(kb))
            rlsyij = 0.5*(rlys(i) + rlys(kb))
	      rlszij = 0.5*(rlzs(i) + rlzs(kb))
            alxavg = perml(1)
            alyavg = perml(2)
            alzavg = perml(3)
	      perml(1) = alxavg*rlsxij
	      perml(2) = alyavg*rlsyij
	      perml(3) = alzavg*rlszij

c            
c derivatives go here (position is defined in stress_perm)
c

c             ts22(jm,1) = 0.5*drlxs(kb,1)
c	      ts22(jm,2) = 0.5*drlxs(kb,2)
c             ts22(jm,3) = 0.5*drlxs(kb,3)
c	      ts22(jm,4) = 0.5*drlxs(kb,4)
	   
c             ts32(jm,1) = 0.5*drlys(kb,1)
c	      ts32(jm,2) = 0.5*drlys(kb,2)
c             ts32(jm,3) = 0.5*drlys(kb,3)
c	      ts32(jm,4) = 0.5*drlys(kb,4)

c             ts42(jm,1) = 0.5*drlzs(kb,1)
c	      ts42(jm,2) = 0.5*drlzs(kb,2)
c             ts42(jm,3) = 0.5*drlzs(kb,3)
c	      ts42(jm,4) = 0.5*drlzs(kb,4)

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
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
             pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2)+delz2/perml(3))
             dpxyp1=pxy/(delx2/perml(1)+dely2/perml(2)+delz2/perml(3))*
     &           delx2/perml(1)**2
             dpxyp2=pxy/(delx2/perml(1)+dely2/perml(2)+delz2/perml(3))*
     &           dely2/perml(2)**2
             dpxyp3=pxy/(delx2/perml(1)+dely2/perml(2)+delz2/perml(3))*
     &           delz2/perml(3)**2

	       ts51(jm) = dpxyp1*alxavg
	       ts52(jm) = dpxyp2*alyavg
	       ts53(jm) = dpxyp3*alzavg

            else	
c
c multiply defined nodes not allowed for stress
c
             pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
              write (*,*) 'MD nodes not allowed for coupled stress'
	         stop
            endif
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t+
     &           delz2/sxzt)
            else
              sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
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
c coupled stress not implemented for 2-D yet
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
            perml(1)=2.*reduction_factor*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*reduction_factor*alykb*alyi/(alykb+alyi)
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
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2))
            else
              pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t)
            else
              sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
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
c
c liquid phase calculations
c
      do 60 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      pxyi=t1(neighc)
      sx4d=t6(neighc)
      axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     & *(cord(kz,igrav)-cord(iz,igrav))
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
c      ial=it12(jm)
c      jml=nelmdg(kz)-neqp1
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
      dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-cord(iz,igrav))
      dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     **(cord(kz,igrav)-cord(iz,igrav))
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
c
       a_axy(iau+nmatavw)=axy
c      a_axy(ial+nmatavw)=-axy

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
c      bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
       a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
       a(jmia+nmat(2))=a(jmia+nmat(2))+dlaei
c      a(ial+nmat(1))=a(ial+nmat(1))-dlapi
c      a(ial+nmat(2))=a(ial+nmat(2))-dlaei
       a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
       a(iau+nmat(2))=a(iau+nmat(2))+dlaekb
c      a(jml+nmat(1))=a(jml+nmat(1))-dlapkb
c      a(jml+nmat(2))=a(jml+nmat(2))-dlaekb
c
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+aexy
c      bp(kz+nrhs(2))=bp(kz+nrhs(2))-aexy
       a(jmia+nmat(3))=a(jmia+nmat(3))+dlepi
       a(jmia+nmat(4))=a(jmia+nmat(4))+dleei
c      a(ial+nmat(3))=a(ial+nmat(3))-dlepi
c      a(ial+nmat(4))=a(ial+nmat(4))-dleei
       a(iau+nmat(3))=a(iau+nmat(3))+dlepkb
       a(iau+nmat(4))=a(iau+nmat(4))+dleekb
c      a(jml+nmat(3))=a(jml+nmat(3))-dlepkb
c      a(jml+nmat(4))=a(jml+nmat(4))-dleekb 

c
c add derivatives wrt displacements
c
c identify derivatives wrt internodal permeability
c      
	  dpxyp1	= ts51(jm) 
	  dpxyp2	= ts52(jm) 
	  dpxyp3	= ts53(jm) 
         pvikb=phi(kb)
         phikb=pvikb-pcp(kb)
         diffp = phikb-phii
         flpart = diffp-0.5*grav*(rolf(i)+rolf(kb))
     &   *(cord(kz,igrav)-cord(iz,igrav))
c         flpart = 0.0
         flpartm = flpart
         dmep1 = dpxyp1*flpartm*axyf
         dmep2 = dpxyp2*flpartm*axyf
         dmep3 = dpxyp3*flpartm*axyf
         flparte = flpart
         deep1 = dpxyp1*flparte*aexyf
         deep2 = dpxyp2*flparte*aexyf
         deep3 = dpxyp3*flparte*aexyf         
c
c  identify derivatives and position in a array
c
         dxiv1 = ts21(1,1)
         dxiv2 = ts21(1,2)
         dxiw1 = ts21(1,3)
         dxiw2 = ts21(1,4)
         idxiv1 = its21(1,1)
         idxiv2 = its21(1,2)
         idxiw1 = its21(1,3)
         idxiw2 = its21(1,4)  
         
         a(idxiv1+nmat(6)) = a(idxiv1+nmat(6)) + dmep1*dxiv1
         a(idxiv2+nmat(6)) = a(idxiv2+nmat(6)) + dmep1*dxiv2  
         a(idxiw1+nmat(7)) = a(idxiw1+nmat(7)) + dmep1*dxiw1
         a(idxiw2+nmat(7)) = a(idxiw2+nmat(7)) + dmep1*dxiw2 
         
         a(idxiv1+nmat(9)) = a(idxiv1+nmat(9)) + deep1*dxiv1
         a(idxiv2+nmat(9)) = a(idxiv2+nmat(9)) + deep1*dxiv2  
         a(idxiw1+nmat(10)) = a(idxiw1+nmat(10)) + deep1*dxiw1
         a(idxiw2+nmat(10)) = a(idxiw2+nmat(10)) + deep1*dxiw2            
               
         dxkbv1 = ts22(jm,1) 
         dxkbv2 = ts22(jm,2) 
         dxkbw1 = ts22(jm,3)
         dxkbw2 = ts22(jm,4) 
         idxkbv1 = its22(jm,1) 
         idxkbv2 = its22(jm,2) 
         idxkbw1 = its22(jm,3)
         idxkbw2 = its22(jm,4)  
	   
         a(idxkbv1+nmat(6)) = a(idxkbv1+nmat(6)) + dmep1*dxkbv1
         a(idxkbv2+nmat(6)) = a(idxkbv2+nmat(6)) + dmep1*dxkbv2  
         a(idxkbw1+nmat(7)) = a(idxkbw1+nmat(7)) + dmep1*dxkbw1
         a(idxkbw2+nmat(7)) = a(idxkbw2+nmat(7)) + dmep1*dxkbw2 
         
         a(idxkbv1+nmat(9)) = a(idxkbv1+nmat(9)) + deep1*dxkbv1
         a(idxkbv2+nmat(9)) = a(idxkbv2+nmat(9)) + deep1*dxkbv2  
         a(idxkbw1+nmat(10)) = a(idxkbw1+nmat(10)) + deep1*dxkbw1
         a(idxkbw2+nmat(10)) = a(idxkbw2+nmat(10)) + deep1*dxkbw2              

         dyiu1 = ts31(1,1)
         dyiu2 = ts31(1,2)
         dyiw1 = ts31(1,3)
         dyiw2 = ts31(1,4)
         idyiu1 = its31(1,1)
         idyiu2 = its31(1,2)
         idyiw1 = its31(1,3)
         idyiw2 = its31(1,4)   
	   
         a(idyiu1+nmat(5)) = a(idyiu1+nmat(5)) + dmep2*dyiu1
         a(idyiu2+nmat(5)) = a(idyiu2+nmat(5)) + dmep2*dyiu2  
         a(idyiw1+nmat(7)) = a(idyiw1+nmat(7)) + dmep2*dyiw1
         a(idyiw2+nmat(7)) = a(idyiw2+nmat(7)) + dmep2*dyiw2   	 
         
         a(idyiu1+nmat(8)) = a(idyiu1+nmat(8)) + deep2*dyiu1
         a(idyiu2+nmat(8)) = a(idyiu2+nmat(8)) + deep2*dyiu2  
         a(idyiw1+nmat(10)) = a(idyiw1+nmat(10)) + deep2*dyiw1
         a(idyiw2+nmat(10)) = a(idyiw2+nmat(10)) + deep2*dyiw2   	        
                
         dykbu1 = ts32(jm,1) 
         dykbu2 = ts32(jm,2) 
         dykbw1 = ts32(jm,3)
         dykbw2 = ts32(jm,4)
         idykbu1 = its32(jm,1) 
         idykbu2 = its32(jm,2) 
         idykbw1 = its32(jm,3)
         idykbw2 = its32(jm,4)          
         
         a(idykbu1+nmat(5)) = a(idykbu1+nmat(5)) + dmep2*dykbu1
         a(idykbu2+nmat(5)) = a(idykbu2+nmat(5)) + dmep2*dykbu2  
         a(idykbw1+nmat(7)) = a(idykbw1+nmat(7)) + dmep2*dykbw1
         a(idykbw2+nmat(7)) = a(idykbw2+nmat(7)) + dmep2*dykbw2  
         
         a(idykbu1+nmat(8)) = a(idykbu1+nmat(8)) + deep2*dykbu1
         a(idykbu2+nmat(8)) = a(idykbu2+nmat(8)) + deep2*dykbu2  
         a(idykbw1+nmat(10)) = a(idykbw1+nmat(10)) + deep2*dykbw1
         a(idykbw2+nmat(10)) = a(idykbw2+nmat(10)) + deep2*dykbw2          

         dziu1 = ts41(1,1)
         dziu2 = ts41(1,2)
         dziv1 = ts41(1,3)
         dziv2 = ts41(1,4)
         idziu1 = its41(1,1)
         idziu2 = its41(1,2)
         idziv1 = its41(1,3)
         idziv2 = its41(1,4)        
	   
         a(idziu1+nmat(5)) = a(idziu1+nmat(5)) + dmep3*dziu1
         a(idziu2+nmat(5)) = a(idziu2+nmat(5)) + dmep3*dziu2  
         a(idziv1+nmat(6)) = a(idziv1+nmat(6)) + dmep3*dziv1
         a(idziv2+nmat(6)) = a(idziv2+nmat(6)) + dmep3*dziv2  
         
         a(idziu1+nmat(8)) = a(idziu1+nmat(8)) + deep3*dziu1
         a(idziu2+nmat(8)) = a(idziu2+nmat(8)) + deep3*dziu2  
         a(idziv1+nmat(9)) = a(idziv1+nmat(9)) + deep3*dziv1
         a(idziv2+nmat(9)) = a(idziv2+nmat(9)) + deep3*dziv2   	           
                
         dzkbu1 = ts42(jm,1) 
         dzkbu2 = ts42(jm,2) 
         dzkbv1 = ts42(jm,3)
         dzkbv2 = ts42(jm,4) 
         idzkbu1 = its42(jm,1) 
         idzkbu2 = its42(jm,2) 
         idzkbv1 = its42(jm,3)
         idzkbv2 = its42(jm,4) 

         a(idzkbu1+nmat(5)) = a(idzkbu1+nmat(5)) + dmep3*dzkbu1
         a(idzkbu2+nmat(5)) = a(idzkbu2+nmat(5)) + dmep3*dzkbu2  
         a(idzkbv1+nmat(6)) = a(idzkbv1+nmat(6)) + dmep3*dzkbv1
         a(idzkbv2+nmat(6)) = a(idzkbv2+nmat(6)) + dmep3*dzkbv2  
         
         a(idzkbu1+nmat(8)) = a(idzkbu1+nmat(8)) + deep3*dzkbu1
         a(idzkbu2+nmat(8)) = a(idzkbu2+nmat(8)) + deep3*dzkbu2  
         a(idzkbv1+nmat(9)) = a(idzkbv1+nmat(9)) + deep3*dzkbv1
         a(idzkbv2+nmat(9)) = a(idzkbv2+nmat(9)) + deep3*dzkbv2   	         
                  
   62 continue
      e n d i f
c
c vapour phase calculations
c
      do 63 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      pxyh=t2(neighc)
      sx4h=t7(neighc)
      vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
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
c      ial=it12(jm)
c      jml=nelmdg(kz)-neqp1
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
     **(cord(kz,igrav)-cord(iz,igrav))
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
c      a_vxy(ial+nmatavw)=-vxy

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+vxy
c      bp(kz+nrhs(1))=bp(kz+nrhs(1))-vxy
      a(jmia+nmat(1))=a(jmia+nmat(1))+dvapi
      a(jmia+nmat(2))=a(jmia+nmat(2))+dvaei
c      a(ial+nmat(1))=a(ial+nmat(1))-dvapi
c      a(ial+nmat(2))=a(ial+nmat(2))-dvaei
      a(iau+nmat(1))=a(iau+nmat(1))+dvapkb
      a(iau+nmat(2))=a(iau+nmat(2))+dvaekb
c      a(jml+nmat(1))=a(jml+nmat(1))-dvapkb
c      a(jml+nmat(2))=a(jml+nmat(2))-dvaekb
c
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+vexy
c      bp(kz+nrhs(2))=bp(kz+nrhs(2))-vexy
      a(jmia+nmat(3))=a(jmia+nmat(3))+dvepi
      a(jmia+nmat(4))=a(jmia+nmat(4))+dveei
c      a(ial+nmat(3))=a(ial+nmat(3))-dvepi
c      a(ial+nmat(4))=a(ial+nmat(4))-dveei
      a(iau+nmat(3))=a(iau+nmat(3))+dvepkb
      a(iau+nmat(4))=a(iau+nmat(4))+dveekb
c      a(jml+nmat(3))=a(jml+nmat(3))-dvepkb
c      a(jml+nmat(4))=a(jml+nmat(4))-dveekb
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
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatc*(t(kb)-ti)
c      bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatc*(t(kb)-ti)
      a(jmia+nmat(3))=a(jmia+nmat(3))-heatc*dtpa(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))-heatc*dtpae(i)
c      a(ial+nmat(3))=a(ial+nmat(3))+heatc*dtpa(i)
c      a(ial+nmat(4))=a(ial+nmat(4))+heatc*dtpae(i)
      a(iau+nmat(3))=a(iau+nmat(3))+heatc*dtpa(kb)
      a(iau+nmat(4))=a(iau+nmat(4))+heatc*dtpae(kb)
c      a(jml+nmat(3))=a(jml+nmat(3))-heatc*dtpa(kb)
c      a(jml+nmat(4))=a(jml+nmat(4))-heatc*dtpae(kb)
66    continue
c
c     bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*(aw*deni(i)+ay*denj(i))+sk(i)
c     bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*(aw*denei(i)+ay*denej(i))+qh(i)
c     a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*(aw*dmpf(i))+dq(i)
c     a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*(aw*dmef(i))+dqt(i)
c     a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*(depf(i)*aw)+dqh(i)
c     a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*(aw*deef(i))+deqh(i)
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqt(i)
      a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*depf(i)+dqh(i)
      a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*deef(i)+deqh(i)

      r e t u r n
      e    n    d
