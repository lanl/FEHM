      subroutine geneq2(i)
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
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/geneq2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!!D2    Rev 2.4   29 Jan 2003 09:05:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:18   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.13   Mon Mar 31 08:37:04 1997   gaz
CD2 minor changes for anisotropic properties
CD2 
CD2    Rev 1.12   Fri Apr 26 15:15:42 1996   gaz
CD2 took out mdnode coding(put elsewhere)
CD2 
CD2    Rev 1.11   Thu Feb 15 10:43:06 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.10   Wed Feb 07 11:00:48 1996   gaz
CD2 replaced if(nmat(1).eq.nmat(11)) then with if(i.gt.neq) then
CD2 also modified for mdnodes
CD2 
CD2    Rev 1.9   Mon Jan 29 16:23:16 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   12/13/95 08:39:24   gaz
CD2 changed setbit counted to accomodate nbits=256
CD2 
CD2    Rev 1.7   08/18/95 10:13:46   llt
CD2 iw was already defined, removed for cray
CD2 
CD2    Rev 1.6   08/07/95 11:35:44   awolf
CD2 a_axy and a_vxy terms fixed for DPDP
CD2 
CD2    Rev 1.5   04/27/95 18:24:10   llt
CD2 added bit routines
CD2 
CD2    Rev 1.4   04/03/95 08:45:40   robinson
CD2 Minor change to mass flow rate array
CD2 
CD2    Rev 1.3   03/10/95 10:41:12   llt
CD2 modified to allow for Richard's Equation
CD2 
CD2    Rev 1.2   01/28/95 14:03:44   llt
CD2 modified for new particle tracking module
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
CD3 i          I          Current node number
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
CD4 t7, t8, dnwgt, upwgt, t9, dglp, dgle, bp, a, rovf, igrav, cord, 
CD4 aw, ay, sk, deni, denj, denei, denej, qh, dmpf, drc, depf, deef,
CD4 deqh, dpcef, dtpa, dtpae, nmat
CD4 
CD4 
CD4 Global Subprograms
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
CD5 thxi         real*8      Thermal conductivity in x-direction
CD5 thyi         real*8      Thermal conductivity in x-direction
CD5 thzi         real*8      Thermal conductivity in x-direction
CD5 thxkb        real*8      Thermal conductivity in x-direction,
CD5                              connected node
CD5 thykb        real*8      Thermal conductivity in x-direction,
CD5                              connected node
CD5 thzkb        real*8      Thermal conductivity in x-direction,
CD5                              connected node
CD5 icode        int         Return from call to mmgetblk
CD5 ti           real*8      Current temperature
CD5 isl          int         Flag denoting whether 1 or 2 phase
CD5 
CD5 Local Subprograms
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
CD9 2.3.2 Heat- and mass-transfer equations
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
      use comwellphys
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
      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd, iwd 
      integer kb_pri, i_dir_gdkm
c following variables are associated with the drift flux nodel    
      real*8 mdrifti,dmdriftpi,dmdriftei,mdrift_part
      real*8 mdriftkb,dmdriftpkb,dmdriftekb
      real*8 area_face,dmdrpkb,dmdrekb,dmdrpi,dmdrei
c gaz debug   
      alxi = ps(1)
c changed by avw -- entered here by seh
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
c storage for upwind
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
         swi=s(i)
      else
         phii=pvii
         swi = 1.0d0
      endif
      ti=t(i)
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
      iz4m1 = 4*(iz-1)+1
c
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
      a_axy(jmia+nmatavw)=sk(i)
c Take care of source/sink term 
c If this is an isothermal air-water simulation
      if (irdof .ne. 13) then
         a_vxy(jmia+nmatavw)=qh(i)
      end if

      do 58 jm=jmi+1,ii2
         iq=iq+1
         kb=nelm(jm)+icd
         it8(iq)=kb
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
c gaz 11-18-2001
c         if(imdnode.ne.0) then
c           imd = mdnodes(i) + mdnodes(kb)
c           if(imd.lt.2) it10(iq) = -abs(it10(iq))
c         endif
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
            if (irdof .ne. 13) then
               phikb=pvikb-pcp(kb)
            else
               phikb=pvikb
            end if
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
c added area term GAZ (11-12-08)      
            t5(neighc)=sx2c*sqrt(dis2)           
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue
         if(irdof.ne.11) then
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
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=axyd
 60         continue
c
c determine upwind nodes and if liquid phase exists
c
            isl=1
            do 61 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c add coding to save upwind position
cc               if(iad.le.iad_up) then
                  fid=0.5
                  axyd=t8(neighc)
                  if(axyd.lt.0.0) fid=dnwgt
                  if(axyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c
cc                  call setbit(nbits,neighc,upwind_l(iz4m1),fid)
c
cc               else
cc                  if(bit(nbits,neighc,upwind_l(iz4m1))) then 
cc                     t9(neighc)=1.0
cc                  else
cc                     t9(neighc)=0.0
cc                  endif
cc               endif
 61         continue
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
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
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
c vapour phase calculations
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
c  assumption- upwind direction for vdrift is the same as vtotal
c  can be different that liquid velocity   
c   
            isl=1
            do 64 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c add coding to save upwind position
cc               if(iad.le.iad_up) then
                  fid=0.5
                  vxyd=t8(neighc)
                  if(vxyd.lt.0.0) fid=dnwgt
                  if(vxyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c
cc                  call setbit(nbits,neighc,upwind_v(iz4m1),fid)
c
cc               else
cc                  if(bit(nbits,neighc,upwind_v(iz4m1))) then 
cc                     t9(neighc)=1.0
cc                  else
cc                     t9(neighc)=0.0
cc                  endif
cc               endif
 64         continue
c
c form equations
c
c
c identifty drift flux mass flows (and derivatives) at node i
c
            if(nwellphy.ne.0) then
               mdrifti = mdriftf(i)
               dmdriftpi = dmdriftp(i)
               dmdriftei = dmdrifte(i)
            else
               mdrift_part = 0.0
            endif  
c                       
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
c drift flux contribution (in terms of mass flux calculated in wellphysics.f)            
c 
                  if(nwellphy.ne.0) then
c face area of CV needed for drift mass flux                  
                     area_face = t5(neighc)          
                     mdriftkb = mdriftf(kb)
                     dmdriftpkb = dmdriftp(kb)
                     dmdriftekb = dmdrifte(kb)                 
                     mdrift_part = (fid*mdriftkb+fid1*mdrifti)*area_face
                     dmdrpkb = fid*dmdriftpkb*area_face
                     dmdrpi = fid1*dmdriftpi*area_face
                     dmdrekb = fid*dmdriftekb*area_face
                     dmdrei = fid1*dmdriftei*area_face
                  else
                     dmdrpi = 0.d0
                     dmdrei = 0.d0
                     dmdrpkb = 0.d0
                     dmdrekb = 0.d0
                  endif                           
                  a_vxy(iau+nmatavw)=vxy + mdrift_part
                  a_vxy(ial+nmatavw)=-(vxy + mdrift_part)

                  bp(iz+nrhs(2))=bp(iz+nrhs(2))+(vxy + mdrift_part)
                  bp(kz+nrhs(2))=bp(kz+nrhs(2))-(vxy + mdrift_part)
                  a(jmia+nmat(3))=a(jmia+nmat(3))+dvapi+dmdrpi
                  a(jmia+nmat(4))=a(jmia+nmat(4))+dvaei+dmdrei
                  a(ial+nmat(3))=a(ial+nmat(3))-dvapi-dmdrpi
                  a(ial+nmat(4))=a(ial+nmat(4))-dvaei-dmdrei
                  a(iau+nmat(3))=a(iau+nmat(3))+dvapkb+dmdrpkb
                  a(iau+nmat(4))=a(iau+nmat(4))+dvaekb+dmdrekb
                  a(jml+nmat(3))=a(jml+nmat(3))-dvapkb-dmdrpkb
                  a(jml+nmat(4))=a(jml+nmat(4))-dvaekb-dmdrekb
 65            continue
            endif
         endif
c     
c     
c 2-d geometry (including radial)
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
            if (irdof .ne. 13) then
               phikb=pvikb-pcp(kb)
            else
               phikb=pvikb
            end if
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
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
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
c added area term GAZ (11-12-08)      
            t5(neighc)=sx2c*sqrt(dis2)                  
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
               sx4d=t6(neighc)
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=axyd
 70         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
            isl=1
            do 71 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c add coding to save upwind position
cc               if(iad.le.iad_up) then
                  fid=0.5
                  axyd=t8(neighc)
                  if(axyd.lt.0.0) fid=dnwgt
                  if(axyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c
cc                  call setbit(nbits,neighc,upwind_l(iz4m1),fid)
c
cc               else
cc                  if(bit(nbits,neighc,upwind_l(iz4m1))) then 
cc                     t9(neighc)=1.0
cc                 else
cc                     t9(neighc)=0.0
cc                  endif
cc               endif
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
                     dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *                    *(cord(kz,igrav)-cord(iz,igrav))
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
c add coding to save upwind position
cc               if(iad.le.iad_up) then
                  fid=0.5
                  vxyd=t8(neighc)
                  if(vxyd.lt.0.0) fid=dnwgt
                  if(vxyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c
cc                  call setbit(nbits,neighc,upwind_v(iz4m1),fid)
c
cc               else
cc                  if(bit(nbits,neighc,upwind_v(iz4m1))) then 
cc                     t9(neighc)=1.0
cc                  else
cc                     t9(neighc)=0.0
cc                  endif
cc               endif
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


