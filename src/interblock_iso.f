      subroutine interblock_iso(idum)
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
CD1 To generate interblock flow for air-water equations.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/interblock_iso.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Tue Jan 30 15:16:10 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   11/15/95 15:39:04   gaz
CD2 corrected bugs,still not fully operational
CD2 
CD2    Rev 1.1   08/18/95 10:14:46   llt
CD2 iw was already defined, removed for cray
CD2 
CD2    Rev 1.0   05/02/95 13:18:22   gaz
CD2 Initial revision.
CD2 
CD2    Rev 1.1   03/18/94 15:40:58   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:10   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name      Type        Description
CD3 
CD3 idum       I          Current node number
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 lenreal, nn, ipt1, ipt2, ipt3, ipt4, ipt5, ipt6, ipt7, ipt8, ipt9,
CD4 ipt10, ipit8, ipit9, ipit10, ipit11, ipit12, sx1, pnx, pny, pnz,
CD4 phi, pcp, dpvti, dil, div, dilp, dile, divp, dive, thx, thy, thz,
CD4 ti, t, s, neq, idualp, nelm, nelmdg, it8, it9, it10, istrw,
CD4 it11, it12, icnl, perml, permv, sx, t1, t2, t3, t4, t5, t6, grav,
CD4 t7, t8, dnwgt, upwgt, t9, dglp, dgle, a, a, rovf, igrav, cord, 
CD4 aw, ay, sk, deni, denj, denei, denej, qh, dmpf, drc, depf, deef,
CD4 deqh, dpcef, dtpa, dtpae, nmat
CD4 
CD4 
CD4 Global Suarograms
CD4
CD4 Name    Type     Description
CD4 
C**********************************************************************
CD5 
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 sx1d         real*8      Volume associated with the current node
CD5 axi          real*8      Permeability in x-direction
CD5 ayi          real*8      Permeability in y-direction
CD5 azi          real*8      Permeability in z-direction
CD5 alxi         real*8      Permeability in x-direction
CD5 alyi         real*8      Permeability in y-direction
CD5 alzi         real*8      Permeability in z-direction
CD5 avxi         real*8      Permeability in x-direction
CD5 avyi         real*8      Permeability in y-direction
CD5 avzi         real*8      Permeability in z-direction
CD5 pvii         real*8      Pressure of vapor
CD5 phii         real*8      Pressure of liquid
CD5 swi          real*8      Saturation
CD5 dili         real*8      Liquid density
CD5 dilkb        real*8      Liquid density, connected node
CD5 divi         real*8      Vapor density
CD5 divkb        real*8      Vapor density, connected node
CD5 icd          int         Integer index parameter
CD5 ii1          int         Integer index parameter
CD5 ii2          int         Integer index parameter
CD5 idg          int         Integer index parameter
CD5 iq           int         Integer index parameter
CD5 jmi          int         Integer index parameter
CD5 jml          int         Integer index parameter
CD5 jmia         int         Integer index parameter
CD5 jm           int         Do-loop index parameter
CD5 neqp1        int         Number of equations plus 1
CD5 ij           int         Do-loop index parameter
CD5 ij1          int         Integer index parameter
CD5 ij2          int         Integer index parameter
CD5 iz           int         Integer index parameter
CD5 kb           int         Integer index parameter
CD5 neighc       int         Integer index parameter
CD5 iw           int         Integer index parameter
CD5 axkb         real*8      Permeability in x-direction, connected node
CD5 aykb         real*8      Permeability in y-direction, connected node
CD5 azkb         real*8      Permeability in z-direction, connected node
CD5 alxkb        real*8      Permeability in x-direction, connected node
CD5 alykb        real*8      Permeability in y-direction, connected node
CD5 alzkb        real*8      Permeability in z-direction, connected node
CD5 sx2t         real*8      Real parameter used in calculation
CD5 sx3t         real*8      Real parameter used in calculation
CD5 sxzt         real*8      Real parameter used in calculation
CD5 sx2c         real*8      Real parameter used in calculation
CD5 sx4c         real*8      Real parameter used in calculation
CD5 sx4d         real*8      Real parameter used in calculation
CD5 sx4h         real*8      Real parameter used in calculation
CD5 sxzc         real*8      Real parameter used in calculation
CD5 pvikb        real*8      Vapor pressure at connected node
CD5 phikb        real*8      Liquid pressure at connected node
CD5 radi         real*8      Parameter used in calculation
CD5 radkb        real*8      Parameter used in calculation
CD5 fid          real*8      Parameter used in calculation
CD5 fid1         real*8      Parameter used in calculation
CD5 axyd         real*8      Parameter used in calculation
CD5 axy          real*8      Parameter used in calculation
CD5 axyf         real*8      Parameter used in calculation
CD5 vxy          real*8      Parameter used in calculation
CD5 vxyf         real*8      Parameter used in calculation
CD5 pvxy         real*8      Parameter used in calculation
CD5 pxy          real*8      Parameter used in calculation
CD5 pxyh         real*8      Parameter used in calculation
CD5 pxyi         real*8      Parameter used in calculation
CD5 heatc        real*8      Parameter used in calculation
CD5 vxyd         real*8      Parameter used in calculation
CD5 storage_term real*8      Parameter used in calculation
CD5 iau          int         Integer index parameter
CD5 ial          int         Integer index parameter
CD5 kz           int         Integer index parameter
CD5 dpvti        real*8      Derivative of capillary pressure with
CD5                              energy variable
CD5 dilpi        real*8      Derivative of liquid transmissibility
CD5                              with pressure
CD5 dilei        real*8      Derivative of liquid transmissibility
CD5                              with energy variable
CD5 divpi        real*8      Derivative of vapor transmissibility
CD5                              with pressure
CD5 divei        real*8      Derivative of vapor transmissibility
CD5                              with energy variable
CD5 dilpkb       real*8      Derivative of liquid transmissibility
CD5                              with pressure, connected node
CD5 dilekb       real*8      Derivative of liquid transmissibility
CD5                              with energy variable, connected node
CD5 divpkb       real*8      Derivative of vapor transmissibility
CD5                              with pressure, connected node
CD5 divekb       real*8      Derivative of vapor transmissibility
CD5                              with energy variable, connected node
CD5 dlaei        real*8      Derivative term used in calculation
CD5 dlaekb       real*8      Derivative term used in calculation
CD5 dvaei        real*8      Derivative term used in calculation
CD5 dvaekb       real*8      Derivative term used in calculation
CD5 dlpi         real*8      Derivative term used in calculation
CD5 dlpkb        real*8      Derivative term used in calculation
CD5 dvpi         real*8      Derivative term used in calculation
CD5 dvpkb        real*8      Derivative term used in calculation
CD5 dlei         real*8      Derivative term used in calculation
CD5 dvei         real*8      Derivative term used in calculation
CD5 dlekb        real*8      Derivative term used in calculation
CD5 dvekb        real*8      Derivative term used in calculation
CD5 dlapi        real*8      Derivative term used in calculation
CD5 dvapi        real*8      Derivative term used in calculation
CD5 dlapkb       real*8      Derivative term used in calculation
CD5 dvapkb       real*8      Derivative term used in calculation
CD5 icode        int         Return from call to mmgetblk
CD5 ti           real*8      Current temperature
CD5 isl          int         Flag denoting whether 1 or 2 phase
CD5 
CD5 Local Suarograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6     
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.2.3 Heat- and mass-transfer equations  
CD9 2.3.3 Noncondensible gas flow equations
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN geneq2
CPS 
CPS mmgetblk - allocate space in arrays used in this routine
CPS Compute terms used later in the calculation
CPS 
CPS IF this node is a matrix node
CPS   Set parameter used later to indicate matrix node
CPS ELSE
CPS   Set parameter used later to indicate primary node
CPS ENDIF
CPS 
CPS FOR each node connected to the current node
CPS 
CPS   Compute node numbers used later
CPS   
CPS   FOR each node connected to this node
CPS     IF this node is the current node
CPS       Set parameter accordingly
CPS     ENDIF
CPS     
CPS   ENDFOR
CPS   
CPS ENDFOR
CPS 
CPS IF this is a 3-D simulation
CPS 
CPS   FOR each node connected to the current node
CPS     Compute average permeability and other terms in equations
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute liquid advection terms
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute liquid upwinding term
CPS   ENDFOR
CPS   
CPS   IF there is liquid flow to or from this node
CPS   
CPS     FOR each node connected to the current node
CPS       Compute liquid advection terms and derivatives
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS   
CPS   FOR each node connected to the current node
CPS     Compute air advection terms
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute air upwinding term
CPS   ENDFOR
CPS   
CPS   IF there is air flow to or from this node
CPS   
CPS     FOR each node connected to the current node
CPS       Compute air advection terms and derivatives
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS   
CPS   FOR each node connected to the current node
CPS     Compute energy balance terms and derivatives
CPS   ENDFOR
CPS   
CPS ELSEIF this is a 2-D simulation
CPS 
CPS   FOR each node connected to the current node
CPS     Compute average permeability and other terms in equations
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute liquid advection terms
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute liquid upwinding term
CPS   ENDFOR
CPS   
CPS   IF there is liquid flow to or from this node
CPS   
CPS     FOR each node connected to the current node
CPS       Compute liquid advection terms and derivatives
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS   
CPS   FOR each node connected to the current node
CPS     Compute air advection terms
CPS   ENDFOR
CPS   
CPS   FOR each node connected to the current node
CPS     Compute air upwinding term
CPS   ENDFOR
CPS   
CPS   IF there is air flow to or from this node
CPS   
CPS     FOR each node connected to the current node
CPS       Compute air advection terms and derivatives
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS   
CPS   FOR each node connected to the current node
CPS     Compute energy balance terms and derivatives
CPS   ENDFOR
CPS   
CPS ENDIF
CPS 
CPS Add accumulation terms to the equations
CPS 
CPS mmrelblk - free space used in temporary arrays
CPS 
CPS END geneq2
CPS 
C**********************************************************************
c
c generate equations for 3-d schemes'isothermal air-water ,
c full derivatives'

      use comflow
      use davidi
      use comji
      use comfi
      use comii
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i,icd,ii1,ii2,idg,iq,jmi,jml,jmia,jm,neqp1,ij,ij1,ij2,iz
      integer kb,neighc,iau,ial,kz,isl,idum,idum1,idum2
      real*8 sx1d,axi,ayi,azi,alxi,alyi,alzi,avxi,avyi,avzi,pvii,phii
      real*8 swi,dili,dilkb,divi,divkb,axkb,aykb,azkb
      real*8 alxkb,alykb,alzkb,sx2c,sx4d,pvikb,phikb,radi
      real*8 radkb,fid,fid1,axyd,axy,axyf,pxy,pxyh,pxyi,sx4h
      real*8 vxyd,vxy,vxyf,ti


      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 dis2 
      real*8 fac1
      real*8 fac2
      real*8 deni_old
      real*8 denei_old
      real*8 fracm
      real*8 frace
      real*8 tol_m 
      real*8 frac_change
      real*8 wat_mass    
      real*8 air_mass    
      integer ifracm   
      integer ifrace   
      real*8 grav_air

      logical bit
      integer iz4m1

      parameter(fac1=1.0,fac2=1.0)
      parameter(tol_m=0.00001,frac_change=0.20)  
      neqp1=neq+1
      if(icons.le.maxit) then
       grav_air=0.0  
      else
       grav_air=grav
      endif

c ******************* do loop on nodes ***********************
      if(idum.eq.0) then
       idum1=1
       idum2=neq
      else if(idum.lt.0) then
       idum1=neq+1
       idum2=neq+neq
      else if(idum.gt.0) then
       idum1=idum
       idum2=idum
      endif
      do i=idum1,idum2
       a(i+nrhs(1))=0.0
       a(i+nrhs(2))=0.0
      enddo
      do i=idum1,idum2
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
      dili=dil(i)
      divi=div(i)
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
      iz4m1 = 4*(iz-1)+1
c
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1
      do 58 jm=jmi+1,ii2
         iq=iq+1
         it8(iq)=nelm(jm)+icd
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
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
c 3-d geometry
c
      if(icnl.eq.0) then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw=it10(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            azkb=pnz(kb)
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,1)+sx(iw,2)+sx(iw,3)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2)+delz2/perml(3))
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue

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
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=axyd
 60      continue
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
 61      continue

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
               dilkb=dil(kb)
               axyf=(fid*dilkb+fid1*dili)
               axy=axyd*axyf
               a(iz+nrhs(1))=a(iz+nrhs(1))+axy
               a(kz+nrhs(1))=a(kz+nrhs(1))-axy
 62         continue
         endif
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
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=vxyd
 63      continue
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
 64      continue

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
               vxyd=t8(neighc)
               divkb=div(kb)
               vxyf=(fid*divkb+fid1*divi)
               vxy=vxyd*vxyf
               a(iz+nrhs(2))=a(iz+nrhs(2))+vxy
               a(kz+nrhs(2))=a(kz+nrhs(2))-vxy
 65         continue
         endif
c     
c 2-d geometry
c
      elseif(icnl.ne.0) then
         radi=cord(iz,3)
         do 69 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iw=it10(jm)
            neighc=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,1)+sx(iw,2))
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            pxy=sx2c*dis2/
     &          (delx2/perml(1)+dely2/perml(2))
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 69      continue

c     
c     liquid phase calculations
c     
         do 70 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            pxyi=t1(neighc)
            sx4d=t6(neighc)
            axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=axyd
 70      continue
c    
c     determine upwind nodes and if liquid phase exists
c    
         isl=1
         do 71 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
c add coding to save upwind position
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
 71      continue
c    

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
               axyd=t8(neighc)
               fid=t9(neighc)
               fid1=1.0-fid
               dilkb=dil(kb)
               axyf=(fid*dilkb+fid1*dili)
               axy=axyd*axyf
               a(iz+nrhs(1))=a(iz+nrhs(1))+axy
               a(kz+nrhs(1))=a(kz+nrhs(1))-axy
 72         continue
         endif
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
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=vxyd
 73      continue
c     
c     determine upwind nodes and if vapour phase exists
c     
         isl=1
         do 74 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
c add coding to save upwind position
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
 74      continue

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
               vxyd=t8(neighc)
               divkb=div(kb)
               vxyf=(fid*divkb+fid1*divi)
               vxy=vxyd*vxyf
               a(iz+nrhs(2))=a(iz+nrhs(2))+vxy
               a(kz+nrhs(2))=a(kz+nrhs(2))-vxy
 75         continue
         endif
      endif
      enddo
c
      do i=idum1,idum2
          a(i+nrhs(1))=a(i+nrhs(1))+sk(i)
          a(i+nrhs(2))=a(i+nrhs(2))+qh(i)
      enddo
c
      fracm=0.0
      frace=0.0
      ifracm = 0 
      ifrace = 0 
      do i=idum1,idum2
c        s(i)      =(denh(i)*sx1(i)+dtot*(-a(i+nrhs(1))))
c    &    /(ps(i)*sx1(i)*rolf(i))
c        denj(i)   =(deneh(i)*sx1(i)+dtot*(-a(i+nrhs(2))))
c    &    /(ps(i)*sx1(i)*rovf(i))
         deni_old  = deni(i)
         denei_old = denei(i)
         wat_mass = ps(i)*rolf(i)/dtot 
         air_mass = ps(i)*rovf(i)/dtot 
         denj(i)   = -a(i+nrhs(1))/sx1(i)*fac1 +deni_old*(1.0-fac1)
         denej(i)  = -a(i+nrhs(2))/sx1(i)*fac2 +denei_old*(1.0-fac2)
         if(abs(deni_old).gt.tol_m) then
          if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
           fracm = (denj(i)-deni_old)/wat_mass
           ifracm = i
          endif
         endif
         if(abs(denei_old).gt.tol_m) then
          if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
           frace = (denej(i)-denei_old)/air_mass 
           ifrace = i
          endif
         endif

      enddo
      if (iout .ne. 0) then
         write(iout,*) 'node max mass correction occurs = ', ifracm 
         write(iout,*) 'max mass correction = ', fracm
         write(iout,*) 'node max air correction occurs = ', ifrace 
         write(iout,*) 'max air correction = ', frace
      end if
      if(iptty.gt.0) then
         write(iptty,*) 'node max mass correction occurs = ', ifracm 
         write(iptty,*) 'max mass correction = ', fracm
         write(iptty,*) 'node max air correction occurs = ', ifrace 
         write(iptty,*) 'max air correction = ', frace
      endif
      if(fracm.gt.frac_change) then
c        mlz=1
      endif
      return
      end
      subroutine interblock_ngas(idum)
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine generates the interblock fluxes for the heat 
CD1  and mass transfer with non-condensible gas included 
CD1
C***********************************************************************
CD2 
CD2    Rev 1.0   01/14/97    pvcs
CD2 original version in process of being certified
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.2.2     Heat and Mass Transfer Equations
CD3  3.2.3     Noncondensible Gas Flow Equations
CD3 
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
C***********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comji
      use davidi
      use comflow
      use comzeoli
      implicit none
c
      logical bit
      integer iz4m1
      integer nmatavw
      integer isl
      integer idum
      integer idum1    
      integer idum2
      integer i, iau, ial, icd, icesd, idg, ifracc, ifrace, ifracm 
      integer ii1, ii2, ij, ij1, ij2, iq, iz, jm, jmi, jmia, jml 
      integer kb, kz, neighc, neqp1
 
      real*8 fac1
      real*8 fac2
      real*8 fac3
      real*8 tol_m,frac_change
      real*8 grav_air
      real*8  dpvti, enli, dilpi, dilei, divpi, divei, deli
      real*8  delei, envi, devi, devei, cnli, dcli, dclei
      real*8  dclci, cnvi, dcvi, dcvei, dcvci, delci, devci
      real*8  dilci, divci, delx2, dely2, delz2, dis2, enlkb
      real*8  cnlkb, aexyf, acxyf, aexy, acxy, pvxy, envkb, cnvkb
      real*8  vexyf, vcxyf, vexy, vcxy, dvai, envmi, devmpi
      real*8  devmei, devmci, heatc, dvakb, envmkb, dvmpkb
      real*8  dvmekb, dvmckb, dva_bar
      real*8  dvame, dvamet, skz, qhz, dskzt1, dqhzt1, dskzt2
      real*8  dqhzt2, water_pressure, fracm, frace, fracc
      real*8  deni_old, denei_old, denpci_old, wat_mass, air_mass
      real*8  energy, sx1d, axi, ayi, azi, alxi, avxi, alyi, avyi
      real*8  alzi, avzi, pvii, phii, dili, divi, ti, swi
      real*8  axkb, aykb, azkb, alxkb, alykb, alzkb
      real*8  sx2c, pvikb, phikb, pxy, pxyi, pxyh, radi, radkb
      real*8  sx4d, axyd, fid, fid1, dilkb, axy, axyf, sx4h
      real*8  vxyd, divkb, vxyf, vxy, psat, dskzt, dqhzt, dpdt

      parameter(fac1=1.0,fac2=1.0,fac3=1.0)
      parameter(tol_m=0.00001,frac_change=0.20)  
      neqp1=neq+1
      if(icons.le.maxit) then
       grav_air=0.0  
      else
       grav_air=grav
      endif

c ******************* do loop on nodes ***********************
      if(idum.eq.0) then
       idum1=1
       idum2=neq
      else if(idum.lt.0) then
       idum1=neq+1
       idum2=neq+neq
      else
       idum1=idum
       idum2=idum
      endif
      do i=idum1,idum2
       a(i+nrhs(1))=0.0
       a(i+nrhs(2))=0.0
       a(i+nrhs(3))=0.0
      enddo
      do i=idum1,idum2
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
      it8(iq)=nelm(jm)+icd
      it9(iq)=jm-ii1+1
      it10(iq)=istrw(jm-neqp1)
      it11(iq)=jm-neqp1
      ij1=nelm(nelm(jm))+1
      ij2=nelmdg(nelm(jm))-1
      do 68 ij=ij1,ij2
      if(nelm(ij).eq.iz) then
      it12(iq)=ij-neqp1
      endif
68    continue
58    continue
c
c set flow factor for multiply defined nodes
c
c     call md_nodes(6,0,iz)
c
c 3-d geometry
c
      if ( icnl.eq.0 )  then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw=it10(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            azkb=pnz(kb)
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,1)+sx(iw,2)+sx(iw,3)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2)+delz2/perml(3))
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
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
            iw=it10(jm)
            neighc=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,1)+sx(iw,2))
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            pxy=sx2c*dis2/
     &          (delx2/perml(1)+dely2/perml(2))
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav*t4(neighc)
            t10(neighc)=sx2c
 69      continue
      endif
c
c change coefficients if ice present
c
c gaz 10-18-2001 ( set content to unattainable value)
!      if (ice.ne.0)  then
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
      do 60 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      pxyi=t1(neighc)
      sx4d=t6(neighc)
      axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
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
 61      continue
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
      dilkb=dil(kb)
c
      aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
      axyf=(fid*dilkb*(1.0-cnlkb)+fid1*dili*(1.0-cnli))
      acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
      aexy=axyd*aexyf
      axy=axyd*axyf
      acxy=axyd*acxyf
c
      a_axy(iau+nmatavw)=axyd*(fid*dilkb+fid1*dili)
      a_axy(ial+nmatavw)=-a_axy(iau+nmatavw)

      a(iz+nrhs(1))=a(iz+nrhs(1))+axy
      a(kz+nrhs(1))=a(kz+nrhs(1))-axy
c
      a(iz+nrhs(2))=a(iz+nrhs(2))+aexy
      a(kz+nrhs(2))=a(kz+nrhs(2))-aexy
c
      a(iz+nrhs(3))=a(iz+nrhs(3))+acxy
      a(kz+nrhs(3))=a(kz+nrhs(3))-acxy
   62 continue
      end if
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
 64      continue
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
      divkb=div(kb)
c
      vexyf=(fid*divkb*envkb+fid1*divi*envi)
      vxyf=(fid*divkb*(1.0-cnvkb)+fid1*divi*(1.0-cnvi))
      vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
      vexy=vxyd*vexyf
      vxy=vxyd*vxyf
      vcxy=vxyd*vcxyf
c
      a_vxy(iau+nmatavw)=vxyd*(fid*divkb+fid1*divi)
      a_vxy(ial+nmatavw)=-a_vxy(iau+nmatavw)

      a(iz+nrhs(1))=a(iz+nrhs(1))+vxy
      a(kz+nrhs(1))=a(kz+nrhs(1))-vxy
c
      a(iz+nrhs(2))=a(iz+nrhs(2))+vexy
      a(kz+nrhs(2))=a(kz+nrhs(2))-vexy
c
      a(iz+nrhs(3))=a(iz+nrhs(3))+vcxy
      a(kz+nrhs(3))=a(kz+nrhs(3))-vcxy
   65 continue
       if(iadif.ne.0) then
c
c add air-water vapor diffusion
c
        if(s(i).lt.1.0) then
         dvai=dva(i)
         envmi=enva(i)
         devmpi=denvap(i)
         devmei=denvae(i)
         devmci=denvac(i)
      do  jm=1,iq
       kb=it8(jm)
       kz=kb-icd
       neighc=it9(jm)
       heatc=t10(neighc)
       iau=it11(jm)
       ial=it12(jm)
       jml=nelmdg(kz)-neqp1
         if(s(kb).lt.1.0) then
c
c only do for saturations lt 1.0
c
c mid-point weight all terms
          fid=0.5            
          fid1=1.0-fid
c energy equation
         dvakb=dva(kb)
         envmkb=enva(kb)
         dvmpkb=denvap(kb)
         dvmekb=denvae(kb)
         dvmckb=denvac(kb)
         dva_bar=2.0*dvai*dvakb/(dvai+dvakb) 
         dvame=(fid*envmkb+fid1*envmi)
         dvamet=dvame*dva_bar
         a(iz+nrhs(2))=a(iz+nrhs(2))+heatc*dvamet*(cnvf(kb)-cnvi)
         a(kz+nrhs(2))=a(kz+nrhs(2))-heatc*dvamet*(cnvf(kb)-cnvi)
c air equation
         a(iz+nrhs(3))=a(iz+nrhs(3))+heatc*dva_bar*(cnvf(kb)-cnvi)
         a(kz+nrhs(3))=a(kz+nrhs(3))-heatc*dva_bar*(cnvf(kb)-cnvi)
c water equation (negative of air equation)
         a(iz+nrhs(1))=a(iz+nrhs(1))-heatc*dva_bar*(cnvf(kb)-cnvi)
         a(kz+nrhs(1))=a(kz+nrhs(1))+heatc*dva_bar*(cnvf(kb)-cnvi)
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
      a(iz+nrhs(2))=a(iz+nrhs(2))+heatc*(t(kb)-ti)
      a(kz+nrhs(2))=a(kz+nrhs(2))-heatc*(t(kb)-ti)
66    continue

      skz = 0.
      qhz = 0.
      dskzt1 = 0.
      dqhzt1 = 0.
      dskzt2 = 0.
      dqhzt2 = 0.
      if( izeolites .ne. 0 .and. kzeol(i) .gt. 0. ) then
         if(ieos(i) .ne. 3) then
            water_pressure = psat(t(i), dpdt, 0)
         else
            water_pressure = psat(t(i), dpdt, 0)
         end if
         call zeolites(t(i), kzeol(i), fwater_old(i),
     2        fwater(i), water_pressure, dpdt, skz, dskzt, qhz, dqhzt)
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
c
c add source terms
c
          a(iz+nrhs(1))=a(iz+nrhs(1))+sk(i)
          a(iz+nrhs(2))=a(iz+nrhs(2))+qh(i)
          a(iz+nrhs(3))=a(iz+nrhs(3))+qc(i)+qhz

      enddo
      if(idpdp.eq.0) then
c
      fracm=0.0
      frace=0.0
      fracc=0.0
      ifracm = 0 
      ifrace = 0 
      ifracc = 0 
      do i=idum1,idum2
         deni_old  = deni(i)
         denei_old = denei(i)
         denpci_old = denpci(i)
         wat_mass = 1.0                
         air_mass = 1.0                
         energy   = 1.0                
         denj(i)   = -a(i+nrhs(1))/sx1(i)*fac1 +deni_old*(1.0-fac1)
         denej(i)  = -a(i+nrhs(2))/sx1(i)*fac2 +denei_old*(1.0-fac2)
         denpcj(i) = -a(i+nrhs(3))/sx1(i)*fac3 +denpci_old*(1.0-fac3)
         if(abs(deni_old).gt.tol_m) then
          if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
           fracm = (denj(i)-deni_old)/wat_mass
           ifracm = i
          endif
         endif
         if(abs(denei_old).gt.tol_m) then
          if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
           frace = (denej(i)-denei_old)/air_mass 
           ifrace = i
          endif
         endif
         if(abs(denpci_old).gt.tol_m) then
          if(abs((denpcj(i)-denpci_old)/energy ).gt.abs(fracc)) then
           fracc = (denpcj(i)-denpci_old)/energy 
           ifracc = i
          endif
         endif

      enddo
      if (iout .ne. 0) then
         write(iout,*) 'node max mass correction occurs = ', ifracm 
         write(iout,*) 'max mass correction = ', fracm
         write(iout,*) 'node max energy correction occurs = ', ifrace 
         write(iout,*) 'max energy correction = ', frace
         write(iout,*) 'node max air correction occurs = ', ifracc 
         write(iout,*) 'max air correction = ', fracc
      end if
      if(iptty.gt.0) then
         write(iptty,*) 'node max mass correction occurs = ', ifracm 
         write(iptty,*) 'max mass correction = ', fracm
         write(iptty,*) 'node max energy correction occurs = ', ifrace 
         write(iptty,*) 'max energy correction = ', frace
         write(iptty,*) 'node max air correction occurs = ', ifracc 
         write(iptty,*) 'max air correction = ', fracc
      endif
      if(fracm.gt.frac_change) then
c        mlz=1
      endif
      endif
      return
      end
      subroutine  interblock_dpdp3
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
CD2 $Log:   /pvcs.config/fehmn/src/interblock_dpdp3.f_a  $
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
c 18-nov-92
c got rid of initialization of na and nb
c must be passed through
c 25-july-92
c started programming implementation of dpdp 3dof
c 1-13-92
c finished initial programming
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.2.2     Heat and Mass Transfer Equations
CD3  3.3.8     Double Porosity / Double Permeability Formulation
CD3  3.4.2     Solve Nonlinear Equation Set at Each Time Step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
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
      use comdti
      use comai
      use comflow
      implicit none

      integer i, ial, iau, id, idg, idl, idum, isl, iz, iz4m1
      integer jm, jml, kb, kz
      integer neighc, neq2, neqp1
      real*8  acxy, acxyf, aexy, aexyf, al0, al1, alen, area
      real*8  axkb, axy, axyd, axyf 
      real*8  cnli, cnlkb, cnvi, cnvkb, coef1, coefmx
      real*8  dfmlis, dfmlkbs, dfmvis, dfmvkbs 
      real*8  dili, dilkb 
      real*8  dist01
      real*8  divi, divkb
      real*8  dvai, dvakb, dvam, dvame
      real*8  enli, enlkb, envi, envkb, envmi, envmkb
      real*8  heatc
      real*8  fid, fid1, fmlkb, fmli, fmvi, fmvkb, frac0, frac1
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb
      real*8  swi, sx2c, sx4d, sx4h
      real*8  ti, tot
      real*8  vcxy, vcxyf, vexy, vxy, vxyf, vexyf, vxyd
      parameter (coefmx=1000000000.0)
      
      neq2   =  neq+neq
      neqp1  =  neq+1
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
         enli=enlf(i)
         envi=envf(i)
         cnli=cnlf(i)
         cnvi=cnvf(i)
         if(irlpt(irlp(i)).eq.7) then
c     add fracture-matrix factors(Sandia)
            fmvi=fmvf(i)
            fmvkb=fmvf(kb)
            fmli=fmlf(i)
            fmlkb=fmlf(kb)
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
            dilkb=dil(kb)/pny(kb)*fmlkb
            divkb=div(kb)/pny(kb)*fmvkb
         else
            dili=dil(i)*fmli
            divi=div(i)*fmvi
            dilkb=dil(kb)*fmlkb
            divkb=div(kb)*fmvkb
         endif
         ti=t(i)
         swi=s(i)
c     
c     form constants for i>neq
c     
         iz=i
         iz4m1 = 4*(iz-1)+1
c
         coef1=-area/dist01
         if(coef1.lt.-coefmx) coef1=coefmx
         if(iupk.ne.0) then
            axkb=max(pny(idl),zero_t)*pnx(idl)
         else
            axkb =  max( pnx(idl  ),pny(idl  ),
     *           pnz(idl  ),zero_t )
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
            dilkb=dil(kb)
            aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            axyf=(fid*dilkb+fid1*dili)
            acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
            aexy=axyd*aexyf
            axy=axyd*axyf
            acxy=axyd*acxyf
            
            a(i+nrhs(1))=a(i+nrhs(1))+axy
            a(i+nrhs(4))=a(i+nrhs(4))-axy
            
            a(i+nrhs(2))=a(i+nrhs(2))+aexy
            a(i+nrhs(5))=a(i+nrhs(5))-aexy
            
            a(i+nrhs(3))=a(i+nrhs(3))+acxy
            a(i+nrhs(6))=a(i+nrhs(6))-acxy
            
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
            divkb=div(kb)
            
            
            vexyf=(fid*divkb*envkb+fid1*divi*envi)
            vxyf=(fid*divkb+fid1*divi)
            vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
            vexy=vxyd*vexyf
            vxy=vxyd*vxyf
            vcxy=vxyd*vcxyf
            
            a(i+nrhs(1))=a(i+nrhs(1))+vxy
            a(i+nrhs(4))=a(i+nrhs(4))-vxy
            
            a(i+nrhs(2))=a(i+nrhs(2))+vexy
            a(i+nrhs(5))=a(i+nrhs(5))-vexy
            
            a(i+nrhs(3))=a(i+nrhs(3))+vcxy
            a(i+nrhs(6))=a(i+nrhs(6))-vcxy
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
               t9(neighc)=0.5
               dvai=dva(i) 
               envmi=enva(i)
               jm=1
               kb=idl    
               kz=kb
               neighc=1
               heatc=t10(neighc)
               fid=t9(neighc)
               fid1=1.0-fid
c     energy equation
               envmkb=enva(kb)
               dvakb=dva(kb)
               dvame=(fid*dvakb*envmkb+fid1*dvai*envmi)
               a(i+nrhs(2))=a(i+nrhs(2))+heatc*dvame*
     &              (cnvf(kb)-cnvi)
               a(i+nrhs(5))=a(i+nrhs(5))-heatc*dvame*
     &              (cnvf(kb)-cnvi)
c     air equation
               dvakb=dva(kb)
               dvam=(fid*dvakb+fid1*dvai)
               a(i+nrhs(3))=a(i+nrhs(3))+heatc*dvam*(cnvf(kb)-cnvi)
               a(i+nrhs(6))=a(i+nrhs(6))-heatc*dvam*(cnvf(kb)-cnvi)
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
         a(i+nrhs(2))=a(i+nrhs(2))+heatc*(t(kb)-ti)
         a(i+nrhs(5))=a(i+nrhs(5))-heatc*(t(kb)-ti)
      enddo
      
      return
      end
      subroutine mass_balance(iflg)
c
c this subroutine corrects the mass and energy balance
c
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comii
      use davidi
      implicit none

      real*8 fracm, frace, fracc
      real*8 deni_old, denei_old, denpci_old
      real*8 wat_mass, air_mass, energy
      real*8 fac1, fac2, fac3
      real*8 tol_m, frac_change
      parameter(fac1=1.0, fac2=1.0, fac3=1.0)
      integer ifracm, ifrace, ifracc
      integer iflg,i,ii
      integer nrhs1 
      integer nrhs2 
      integer nrhs3 
c
      if(iexrlp.ge.0) return
c
      if(iflg.eq.1) then
       if(ico2.lt.0) then
        if(idpdp.eq.0) then
         call airctr(3,0)
         call interblock_iso(0)
        else if(idpdp.ne.0) then
c *******************************
          tol_m=0.0
          frac_change=0.1

         do i=1,neq
            a(i+nrhs(1))=0.0
            a(i+nrhs(2))=0.0
            a(i+nrhs(3))=0.0
            a(i+nrhs(4))=0.0
         enddo
            nrhs1=nrhs(1)
            nrhs2=nrhs(2)
            nrhs(1)=nrhs1
            nrhs(2)=nrhs2
c     
            call interblock_iso(0)
c
            nrhs(1)=nrhs(3)
            nrhs(2)=nrhs(3)
c     
            call interblock_iso(-1)
c
            nrhs(1)=nrhs1
            nrhs(2)=nrhs2

c     call md_modes to complete equations for multiply defined nodes
c     
c        if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     now get coupling for influence terms
c     
        call interblock_dpdpfa
c
         fracm=0.0
         frace=0.0
         fracc=0.0
         ifracm = 0 
         ifrace = 0 
         ifracc = 0 
         do i=1,neq            
            deni_old  = deni(i)
            denei_old = denei(i)
            wat_mass = 1.0                
            air_mass = 1.0                
            denj(i)   = -a(i+nrhs(1))/sx1(i)*fac1 +deni_old*(1.0-fac1)
            denej(i)  = -a(i+nrhs(2))/sx1(i)*fac2 +denei_old*(1.0-fac2)
            denpcj(i) = -a(i+nrhs(3))/sx1(i)*fac3 +denpci_old*(1.0-fac3)
            if(abs(denj(i)).gt.tol_m) then
             if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
              fracm = (denj(i)-deni_old)/wat_mass
              ifracm = i
             endif
            endif
            if(abs(denej(i)).gt.tol_m) then
             if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
              fracc = (denej(i)-denei_old)/air_mass 
              ifracc = i
             endif
            endif
         enddo
         do ii=1,neq            
            i=ii+neq
            deni_old  = deni(i)
            denei_old = denei(i)
            wat_mass = 1.0                
            air_mass = 1.0                
            denj(i)   = -a(ii+nrhs(4))/sx1(i)*fac1+deni_old*(1.0-fac1)
            denej(i)  = -a(ii+nrhs(5))/sx1(i)*fac2+denei_old*(1.0-fac2)
            if(abs(denj(i)).gt.tol_m) then
             if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
              fracm = (denj(i)-deni_old)/wat_mass
              ifracm = i
             endif
            endif
            if(abs(denej(i)).gt.tol_m) then
             if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
              fracc = (denej(i)-denei_old)/air_mass 
              ifracc = i
             endif
            endif
         enddo
         if (iout .ne. 0) then
            write(iout,*) 'node max mass correction at  ', ifracm 
            write(iout,*) 'max mass correction = ', fracm
            write(iout,*) 'node max air correction at ', ifracc 
            write(iout,*) 'max air correction = ', fracc
         end if
         if(iptty.gt.0) then
            write(iptty,*) 'node max mass correction at  ', ifracm 
            write(iptty,*) 'max mass correction = ', fracm
            write(iptty,*) 'node max air correction at ', ifracc 
            write(iptty,*) 'max air correction = ', fracc
         endif
         if(fracm.gt.frac_change) then
c                    mlz=1
         endif
c ******************************
      endif
      else if(ico2.eq.0) then
c       if(idualp.eq.0) then
c        call thermw(0)
c        call interblock_h2o(0)
c       endif                    
      else if(ico2.gt.0) then
         if(idpdp.eq.0) then
            call interblock_ngas(0)
         else if(idpdp.ne.0) then
c *******************************
            tol_m=0.0
            frac_change=0.1
c         call thrmwc(neq)

            do i=1,neq
               a(i+nrhs(1))=0.0
               a(i+nrhs(2))=0.0
               a(i+nrhs(3))=0.0
               a(i+nrhs(4))=0.0
               a(i+nrhs(5))=0.0
               a(i+nrhs(6))=0.0
            enddo
            nrhs1=nrhs(1)
            nrhs2=nrhs(2)
            nrhs3=nrhs(3)
            nrhs(1)=nrhs1
            nrhs(2)=nrhs2
            nrhs(3)=nrhs3
c     
            call interblock_ngas(0)
c
            nrhs(1)=nrhs(4)
            nrhs(2)=nrhs(5)
            nrhs(3)=nrhs(6)
c     
            call interblock_ngas(-1)
c
            nrhs(1)=nrhs1
            nrhs(2)=nrhs2
            nrhs(3)=nrhs3

c     call md_modes to complete equations for multiply defined nodes
c     
c        if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     now get coupling for influence terms
c     
         call interblock_dpdp3
c
         fracm=0.0
         frace=0.0
         fracc=0.0
         ifracm = 0 
         ifrace = 0 
         ifracc = 0 
         do i=1,neq            
            deni_old  = deni(i)
            denei_old = denei(i)
            denpci_old = denpci(i)
            wat_mass = 1.0                
            air_mass = 1.0                
            energy   = 1.0                
            denj(i)   = -a(i+nrhs(1))/sx1(i)*fac1 +deni_old*(1.0-fac1)
            denej(i)  = -a(i+nrhs(2))/sx1(i)*fac2 +denei_old*(1.0-fac2)
            denpcj(i) = -a(i+nrhs(3))/sx1(i)*fac3 +denpci_old*(1.0-fac3)
            if(abs(denj(i)).gt.tol_m) then
             if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
              fracm = (denj(i)-deni_old)/wat_mass
              ifracm = i
             endif
            endif
            if(abs(denej(i)).gt.tol_m) then
             if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
              frace = (denej(i)-denei_old)/air_mass 
              ifrace = i
             endif
            endif
            if(abs(denpcj(i)).gt.tol_m) then
             if(abs((denpcj(i)-denpci_old)/energy ).gt.abs(fracc)) then
              fracc = (denpcj(i)-denpci_old)/energy 
              ifracc = i
             endif
            endif
         enddo
         do ii=1,neq            
            i=ii+neq
            deni_old  = deni(i)
            denei_old = denei(i)
            denpci_old = denpci(i)
            wat_mass = 1.0                
            air_mass = 1.0                
            energy   = 1.0                
            denj(i)   = -a(ii+nrhs(4))/sx1(i)*fac1+deni_old*(1.0-fac1)
            denej(i)  = -a(ii+nrhs(5))/sx1(i)*fac2+denei_old*(1.0-fac2)
            denpcj(i)= -a(ii+nrhs(6))/sx1(i)*fac3+denpci_old*(1.0-fac3)
            if(abs(denj(i)).gt.tol_m) then
             if(abs((denj(i)-deni_old)/wat_mass).gt.abs(fracm)) then
              fracm = (denj(i)-deni_old)/wat_mass
              ifracm = i
             endif
            endif
            if(abs(denej(i)).gt.tol_m) then
             if(abs((denej(i)-denei_old)/air_mass ).gt.abs(frace)) then
              frace = (denej(i)-denei_old)/air_mass 
              ifrace = i
             endif
            endif
            if(abs(denpcj(i)).gt.tol_m) then
             if(abs((denpcj(i)-denpci_old)/energy ).gt.abs(fracc)) then
              fracc = (denpcj(i)-denpci_old)/energy 
              ifracc = i
             endif
            endif
         enddo
         if (iout .ne. 0) then
            write(iout,*) 'node max mass correction at  ', ifracm 
            write(iout,*) 'max mass correction = ', fracm
            write(iout,*) 'node max energy correction at   ', ifrace 
            write(iout,*) 'max energy correction = ', frace
            write(iout,*) 'node max air correction at ', ifracc 
            write(iout,*) 'max air correction = ', fracc
         end if
         if(iptty.gt.0) then
            write(iptty,*) 'node max mass correction at  ', ifracm 
            write(iptty,*) 'max mass correction = ', fracm
            write(iptty,*) 'node max energy correction at   ', ifrace 
            write(iptty,*) 'max energy correction = ', frace
            write(iptty,*) 'node max air correction at ', ifracc 
            write(iptty,*) 'max air correction = ', fracc
         endif
         if(fracm.gt.frac_change) then
c                    mlz=1
         endif
c ******************************
      endif                    
      endif
      else if(iflg.eq.2) then
       if(ico2.lt.0) then
       do i=1,neq
        deni(i)=denj(i)
        denei(i)=denej(i)
       enddo
       else if(ico2.gt.0) then
        if(idpdp.eq.0) then
         do i=1,neq
          deni(i)=denj(i)
          denei(i)=denej(i)
          denpci(i)=denpcj(i)
         enddo
        else if(idpdp.ne.0) then
         do i=1,n
          deni(i)=denj(i)
          denei(i)=denej(i)
          denpci(i)=denpcj(i)
         enddo
        endif
       endif
      endif
      return
      end
      subroutine  interblock_dpdpfa
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
CD2 $Log:   /pvcs.config/fehmn/src/interblock_dpdpfa.f_a  $
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
CD3  3.2.2     Heat and Mass Transfer Equations
CD3  3.3.8     Double Porosity / Double Permeability Formulation
CD3  3.4.2     Solve Nonlinear Equation Set at Each Time Step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
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

      integer i, id, idg, idl, idum, isl, iz
      integer jm, kb, kz
      integer neighc, neq2, neqp1, nmat2avw
      real*8  al0, al1, alen, area
      real*8  axkb, axy, axyd, axyf 
      real*8  coef1, coefmx
      real*8  dfmlis, dfmlkbs, dfmvis, dfmvkbs 
      real*8  dili, dilkb 
      real*8  dist01
      real*8  divi, divkb
      real*8  dpvti
      real*8  fid, fid1, fmli, fmlkb, fmvi, fmvkb, frac0, frac1
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyh, pxyi 
      real*8  radi, radkb
      real*8  swi, sx2c, sx4d, sx4h
      real*8  ti, tot
      real*8  vxy, vxyf, vxyd
      parameter (coefmx=1000000000.0)

      neq2   =  neq+neq
      neqp1  =  neq+1
      nmat2avw = 2*(nelm(neqp1)-neqp1)
c     
c     compute the fracture-matrix transfer terms
c     
c     loop on nodes
c     
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
            dilkb=dil(kb)/pny(kb)*fmlkb
            divkb=div(kb)/pny(kb)*fmvkb
         else
            dili=dil(i)*fmli
            divi=div(i)*fmvi
            dilkb=dil(kb)*fmlkb
            divkb=div(kb)*fmvkb
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
            axkb =  max( pnx(idl  ),pny(idl  ),
     *           pnz(idl  ),zero_t )
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
            axyf=(fid*dilkb+fid1*dili)
            axy=axyd*axyf
            a_axy(nmat2avw+i) = axy
            
            a(i+nrhs(1))=a(i+nrhs(1))+axy
            a(i+nrhs(3))=a(i+nrhs(3))-axy
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
            vxyf=(fid*divkb+fid1*divi)
            vxy=vxyd*vxyf
            a_vxy(nmat2avw+i) = vxy
            
            a(i+nrhs(2))=a(i+nrhs(2))+vxy
            a(i+nrhs(4))=a(i+nrhs(4))-vxy
         endif
      enddo
      
      return
      end
