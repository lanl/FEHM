      subroutine  dpdpta(ispecies,spec_numf,matnumf,
     2     spec_numm,matnumm)
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
CD1  solution matrix for a tracer simulation.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dpdpta.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Mon Jun 10 10:59:28 1996   hend
CD2 Updated to use variable diffusion (concadiff)
CD2 
CD2    Rev 1.6   Thu Feb 15 09:50:48 1996   zvd
CD2 Corrected purpose and requirements.
CD2 
CD2    Rev 1.5   Wed Jan 10 11:12:04 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.4   Wed Jan 10 09:42:24 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   08/07/95 11:13:54   awolf
CD2 Now is called with ispecies arguement. Fixed diskbl and diskbv terms too,
CD2 
CD2    Rev 1.2   03/23/94 14:41:04   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:26   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:16   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3  
CD3 2.3.4 Solute-transport equations
CD3 2.4.9 Double-porosity/double-permeability formulation
CD3 2.5.2 Solve nonlinear equation set at each time step
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

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comhi
      use comji
      use davidi
      use comrxni
      use comrxnb
      use comflow
      use comdti
      use comai
      implicit none

      integer ispecies
      integer matnumf
      integer spec_numf
      integer matnumm
      integer spec_numm
      integer i, id, idg, idl, idum, ireg, isl, itr, iz
      integer jm, kb, kz, mi, neighc, neq2, neqp1
      real*8  newdiff, concadiff, satr
      real*8  alen, al0, al1, area
      real*8  anli, anlkb, anlri, anvi, anvkb, anvri
      real*8  alxi,  alyi, alzi,  avxi, avyi, avzi, axi, ayi, azi
      real*8  athkb, axkb, axy, axyd, axyf
      real*8  danli, danlkb, danlri, danvi, danvkb, danvri, dfeei 
      real*8  dili, dilkb, divi, divkb, dlaei, dlaekb, dvaei, dvaekb 
      real*8  coef1, heatc, diskbl, diskbv, dist01
      real*8  fid, fid1, frac0, frac1
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb
      real*8  swi, sx1d, sx2c, sx2tl, sx2tv, sx4d, sx4h, tot
      real*8  vxy, vxyd, vxyf

      neq2   =  neq+neq
      neqp1  =  neq+1
c
c zero out transfer terms
c
c
c compute the fracture-matrix transfer terms
c
c loop on nodes
      do 20 id=1,neq
c identify diagonal member
         idg=nelmdg(id)
         idum=idg-neqp1
c identify matrix node
         i=id
         idl=id+neq
         tot=sx1(id)+sx1(idl)
         frac0=sx1(id)/tot
         frac1=sx1(idl)/tot
         alen=apuv1(id)
         al0=frac0*alen
         al1=frac1*alen
c
c determine contribution to node id
c
c
         dist01=0.5*(al0+al1)
         if(dist01.le.1.e-15) dist01=1.e-15
         area=tot/alen
c
         itr=i+npn
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
         swi=s(i)
         dili=dil(i)
         divi=div(i)
         anli=anl(itr)
         anvi=anv(itr)
         danli=danl(itr)
         danvi=danv(itr)
         anlri=anli*rolf(i)
         anvri=anvi*rovf(i)
         danlri=danli*rolf(i)
         danvri=danvi*rovf(i)
         dfeei=0.
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            satr = s(idl)
         else
            satr = 1.0d0
         end if
c
c form constants for i>neq
c
         iz=i
c
         coef1=-area/dist01
         axkb =  max( pnx(idl  ),pny(idl  ),
     *        pnz(idl  ),zero_t )
         athkb =  max( thx(idl  ),thy(idl  ),
     *        thz(idl  ),zero_t )
         mi = idl + npt(ispecies)
         ireg = itrc(mi)
c        diskbl = 1e-15 * ps(idl)
c        diskbv = 1e-30 * ps(idl)
         newdiff =  concadiff(1,mflagl(ispecies,ireg),
     &        diffmfl(ispecies,ireg),ps(idl),satr,
     &        phi(i),t(i))
c----- phs 9/26/01   added saturation term to diskbl and diskbv

         diskbl = newdiff * ps(idl) * satr
         newdiff = concadiff(2,mflagv(ispecies,ireg),
     &        diffmfv(ispecies,ireg),ps(idl),satr,
     &        phi(i),t(i))
         diskbv = newdiff * ps(idl) * (1-satr)
         diskbl =  max( diskbl,zero_t )
         diskbv =  max( diskbv,zero_t )
c
c 2-d geometry
c
         if ( icnl .ne. 0 )  then
            radi=cord(iz,3)
         elseif ( icnl .eq. 0 )  then
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
c      thxkb=athkb
         sx2tl=sx2c*diskbl
         sx2tv=sx2c*diskbv
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
         t5(neighc)=sx2tl
         t5v(neighc)=sx2tv
         t6(neighc)=0.0
         t7(neighc)=0.0
c
c choose vapour or liquid tracer
c
         if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2) then
c
c liquid phase calculations
c
            jm=1
            kb=idl
            kz=kb
            neighc=1
            pxyi=t1(neighc)
            sx4d=t6(neighc)
            axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=axyd
c
c determine upwind nodes and if liquid phase exists
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
            if(swi+satr.ne.0.0) isl=1
c
c form equations
c
c     if(isl.ne.0) then
            jm=1
            kb=idl
            kz=kb
            neighc=1
            axyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            dilkb=dil(kb)
            anlkb=anl(kb+npn)
            danlkb=danl(kb+npn)
            axyf=(fid*dilkb*anlkb+fid1*dili*anli)
            axy=axyd*axyf
            dlaei=axyd*fid1*dili*danli
            dlaekb=axyd*fid*dilkb*danlkb
c
            bp(i+nrhs(spec_numf))=bp(i+nrhs(spec_numf))+axy
            bp(i+nrhs(spec_numm))=bp(i+nrhs(spec_numm))-axy
            a(idum+nmat(matnumf))=a(idum+nmat(matnumf))+dlaei
            a(idum+nmat(matnumm-1))=a(idum+nmat(matnumm-1))-dlaei
            a(idum+nmat(matnumf+1))=a(idum+nmat(matnumf+1))+dlaekb
            a(idum+nmat(matnumm))=a(idum+nmat(matnumm))-dlaekb
c     endif
c
c add dispersion term for liquid
c
            jm=1
            kb=idl
            kz=kb
            neighc=1
            heatc=t5(neighc)
            axy=heatc*(anl(kb+npn)*rolf(kb)-anlri)
            dlaei=-heatc*danlri
            dlaekb=heatc*rolf(kb)*danl(kb+npn)
            bp(i+nrhs(spec_numf))=bp(i+nrhs(spec_numf))+axy
            bp(i+nrhs(spec_numm))=bp(i+nrhs(spec_numm))-axy
            a(idum+nmat(matnumf))=a(idum+nmat(matnumf))+dlaei
            a(idum+nmat(matnumm-1))=a(idum+nmat(matnumm-1))-dlaei
            a(idum+nmat(matnumf+1))=a(idum+nmat(matnumf+1))+dlaekb
            a(idum+nmat(matnumm))=a(idum+nmat(matnumm))-dlaekb
         end if
c
         if (icns(nsp).eq.-1.or.abs(icns(nsp)).eq.2) then
c
c vapour phase calculations
c
            jm=1
            kb=idl
            kz=kb
 
            neighc=1
            pxyh=t2(neighc)
            sx4h=t7(neighc)
            vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=vxyd
c
c determine upwind nodes and if vapour phase exists
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
c form equations
c
c     if(isl.ne.0) then
            jm=1
            kb=idl
            kz=kb
            neighc=1
            vxyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            divkb=div(kb)
            anvkb=anv(kb+npn)
            danvkb=danv(kb+npn)
            vxyf=(fid*divkb*anvkb+fid1*divi*anvi)
            vxy=vxyd*vxyf
            dvaei=vxyd*fid1*divi*danvi
            dvaekb=vxyd*fid*divkb*danvkb
c
            bp(i+nrhs(spec_numf))=bp(i+nrhs(spec_numf))+vxy
            bp(i+nrhs(spec_numm))=bp(i+nrhs(spec_numm))-vxy
            a(idum+nmat(matnumf))=a(idum+nmat(matnumf))+dvaei
            a(idum+nmat(matnumm-1))=a(idum+nmat(matnumm-1))-dvaei
            a(idum+nmat(matnumf+1))=a(idum+nmat(matnumf+1))+dvaekb
            a(idum+nmat(matnumm))=a(idum+nmat(matnumm))-dvaekb
c     endif
c
c add dispersion term
c
            jm=1
            kb=idl
            kz=kb
            neighc=1
            heatc=t5v(neighc)
            vxy=heatc*(anv(kb+npn)*rovf(kb)-anvri)
            dvaei=-heatc*danvri
            dvaekb=heatc*rovf(kb)*danv(kb+npn)
c
            bp(i+nrhs(spec_numf))=bp(i+nrhs(spec_numf))+vxy
            bp(i+nrhs(spec_numm))=bp(i+nrhs(spec_numm))-vxy
            a(idum+nmat(matnumf))=a(idum+nmat(matnumf))+dvaei
            a(idum+nmat(matnumm-1))=a(idum+nmat(matnumm-1))-dvaei
            a(idum+nmat(matnumf+1))=a(idum+nmat(matnumf+1))+dvaekb
            a(idum+nmat(matnumm))=a(idum+nmat(matnumm))-dvaekb
         endif
20    continue

      return
      end
