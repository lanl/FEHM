      subroutine geneq3(i)
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
CD1  This subroutine generates the equations for heat conduction.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/geneq3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Mon Mar 31 08:37:08 1997   gaz
CD2 minor changes for anisotropic properties
CD2 
CD2    Rev 1.3   Thu Feb 15 10:43:08 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.2   Wed Jan 10 13:31:20 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:51:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:14   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.1 Heat-conduction equations
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

      use comai
      use combi
      use comci
      use comdi 
      use comdti
      use comei
      use comflow, only : flag_heat_out, e_axy_cond, e_cond_nodal
      use comgi
      use comji
      use davidi
      implicit none

      integer i
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1
      integer imd,iwd
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8  dtpaei, heatc  
      real*8  radi, radkb
      real*8  sx1d, sx2c, sx2t, sx3c, sx3t, sxzt, sxzc
      real*8  thxi, thxkb, thyi, thykb, thzi, thzkb, ti
      real*8 heatt
      parameter(dis_tol=1.d-12)
c
c generate equations for  heat conduction
c
      sx1d=sx1(i)
      thxi=thx(i)
      thyi=thy(i)
      thzi=thz(i)
      ti=t(i)
      dtpaei=dtpae(i)
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
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1
      sx_min = 0.0d00
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
      if(icnl.eq.0) then
c     
c     3-d geometry
c     
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iwd=it10(jm)
            iw = abs(iwd)
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
            thxkb=thx(kb)
            thykb=thy(kb)
            thzkb=thz(kb)
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            sxzt=2.*thzi*thzkb/(thzi+thzkb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sxzc=sx2c*dis2/
     &              (delx2/sx2t+dely2/sx3t+
     &              delz2/sxzt)
            else
               sxzc=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
            t5(neighc)=sxzc               
 59      continue
      else
c     
c     2-d coding
c     
         radi=cord(iz,3)
         do 60 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iwd=it10(jm)
            iw = abs(iwd)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
            thxkb=thx(kb)
            thykb=thy(kb)
            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/
     &              (delx2/sx2t+dely2/sx3t)
            else
               sxzc=sx2c*sx_mult*max(sx2t,sx3t)
            endif
            t5(neighc)=sx3c
 60      continue
      endif
c     
c     form equations
c     
      do jm=1,iq
         kb=it8(jm)
         kz=kb-icd
         neighc=it9(jm)
         heatc=t5(neighc)
         iau=it11(jm)
         ial=it12(jm)
         jml=nelmdg(kz)-neqp1
         heatc=t5(neighc)

c     s kelkar 3 July 2014, for calculating heat flow vectors
         heatt = +heatc*(t(kb)-ti)
         if(flag_heat_out) then
            e_axy_cond(iau)=+heatt
            e_axy_cond(ial)=-heatt
         endif

         bp(i)=bp(i)+heatc*(t(kb)-ti)
         bp(kb)=bp(kb)-heatc*(t(kb)-ti)
         a(jmia)=a(jmia)-heatc*dtpae(i)
         a(ial)=a(ial)+heatc*dtpae(i)
         a(iau)=a(iau)+heatc*dtpae(kb)
         a(jml)=a(jml)-heatc*dtpae(kb)
      end do
c     
      bp(i)=bp(i)+sx1d*(aw*denei(i)+ay*denej(i))+qh(i)
      a(jmia)=a(jmia)+sx1d*aw*deef(i)+deqh(i)

      r e t u r n
      e    n    d
