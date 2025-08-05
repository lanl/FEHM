      subroutine geneq2_nondarcy(i)
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
CD1 To generate equations isothermal air-water solution at each node with
CD1 non darcy velocity model
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 04-28-24     G. Zyvoloski      Initial implementation.
CD2
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
      use com_nondarcy
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
c gaz 110424 nd_flow liq and gas phase
      real*8 vxy_nd, dvapi_nd, dvapkb_nd, vel_nd
      real*8 axy_nd, dlapi_nd, dlapkb_nd
      real*8 dlaei_nd, dlaekb_nd, dvaei_nd, dvaekb_nd
      real*8 velij_nd,dvelpi,dvelpj,aij,den_term
      real*8 axyd_nd,vxyd_nd,kij,kij_tol
      parameter(kij_tol=1.d-24)

c gaz debug   
      alxi = ps(1)+iad
c gaz debug 110724 
c      nd_test = 100   (turned off)
c changed by avw -- entered here by  seh
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
c gaz 102914 
            t15(neighc)= (pxy)/(sx2c+kij_tol)
 59      continue
c irdof.eq.0:vapor phase only
         if(irdof.ne.11) then
c     
c liquid phase calculations
c
c gaz 091324 if(nd_flow) then
            do jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxy = t3(neighc)
               pxyi=t1(neighc)
               sx4d=t6(neighc)
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
c axyd units m**2*(area/dis)*Mpa
               g_term = 0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
               dg_termpi = 0.5*sx4d*dglp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
               dg_termpkb = 0.5*sx4d*dglp(kb)*
     &                (cord(kz,igrav)-cord(iz,igrav))
               daxydpi = -pxy+dg_termpi
               daxydpkb = pxy+dg_termpkb
               
               daxydei=pxy*dpvti+0.5*sx4d*dgle(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
               daxydekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav))      
               t8(neighc)=axyd
               fid=0.5
c gaz 220125,250525 removed upwind for liq
                 if(axyd.lt.0.0) fid=dnwgt
                 if(axyd.gt.0.0) fid=upwgt
                 fid1= 1.d0-fid
                aij = abs(t5(neighc)) 
                xrl_nd = fid*rlf_nd(kb)+fid1*rlf_nd(i)
                kij = t15(neighc)*1.d-6*xrl_nd
                call nd_props(0,icd,1,iq,i,kb,jm,0.5d0)
                vel_nd = axyd/(aij*muij+kij_tol)

                call nd_flow_vel(1,icd,1,iq,axyd,vel_nd,aij,kij,
     &                 fid,dlapi_nd,dlapkb_nd,dlaei_nd,dlaekb_nd,i,kb,jm)
                aij = abs(t5(neighc)) 
c gaz 090125                
                axyd_nd = vel_nd*aij*muij
                call nd_flow_vel(2,icd,1,iq,axyd,vel_nd,aij,kij,
     &                 fid,dlapi_nd,dlapkb_nd,dlaei_nd,dlaekb_nd,i,kb,jm)  
    
c determine upwind direction
c find upwinding
                t8_nd(neighc)=axyd_nd 
                 fid=0.5
                 if(axyd_nd.lt.0.0) fid=dnwgt
                 if(axyd_nd.gt.0.0) fid=upwgt
                 t9_nd(neighc)=fid  
c gaz 103024 need derivatived muij
             aij = abs(t5(neighc))  
             t13(neighc) = dlapi_nd*aij*muij + vel_nd*aij*dmuijpi
             t14(neighc) = dlapkb_nd*aij*muij + vel_nd*aij*dmuijpj
             t18(neighc) = dlaei_nd*aij*muij + vel_nd*aij*dmuijpi
             t19(neighc) = dlaekb_nd*aij*muij + vel_nd*aij*dmuijpj
             continue   
           enddo       
c
c
            isl=1
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

c gaz 091624 aij = t5() area
c             if(nd_flow) then
c   gaz db 010125            aij=abs(t5(neighc)
                aij = abs(t5(neighc)) 
                aij = 1.d0  
                axyd_nd = t8_nd(neighc)
                fid=t9_nd(neighc)
                fid1=1.0-fid  
                axyf=(fid*dilkb+fid1*dili)
                axy_nd = axyd_nd*axyf
                
                dlapi_nd = t13(neighc)
                dlapkb_nd = t14(neighc)
                dlapi_nd = dlapi_nd*axyf+axyd_nd*fid1*dilpi
                dlapkb_nd = dlapkb_nd*axyf+axyd_nd*fid*dilpkb
                if(irdof.ne.13) then
                  dilekb=dile(kb)
                  dlaei_nd = t18(neighc)
                  dlaekb_nd = t19(neighc)
                  dlaei_nd=dlaei_nd*axyf+axyd_nd*fid1*dilei
                  dlaekb_nd=dlaekb_nd*axyf+axyd_nd*fid*dilekb          
                endif
                continue

                axy = axy_nd
                dlapi = dlapi_nd
                dlapkb = dlapkb_nd
                dlaei = dlaei_nd
                dlaekb = dlaekb_nd
c             endif     
   
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
c
c irdof.eq.11:vapor phase only
c irdof.eq.13:liquid phase only
         if(irdof.ne.13) then   
c     
c vapour phase calculations
c
            do jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               pvxy = t4(neighc)               
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
c vxyd units m**2*(area/dis)*Mpa

                  divkb=div(kb)
                  divpkb=divp(kb)
                  divekb=dive(kb)                
               g_term = 0.5*sx4d*(rovf(i)+rovf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
               dg_termpi = 0.5*sx4d*dgvp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
               dg_termpkb = 0.5*sx4d*dgvp(kb)*
     &                (cord(kz,igrav)-cord(iz,igrav))
               daxydpi = -pvxy+dg_termpi
               daxydpkb = pvxy+dg_termpkb

               daxydei= 0.5*sx4d*dgve(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
               daxydekb= 0.5*sx4d*dgve(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav))
               call nd_props(0,icd,1,iq,i,kb,jm,0.5d0)
c gaz 112424 
c vxyd is the darcy based velocity
c                fid = 0.5d0
c gaz 220525 added upwind calc
c                if(vxyd.lt.0.0) fid=dnwgt
c                if(vxyd.gt.0.0) fid=upwgt
c                fid1= 1.d0-fid
               aij=abs(t5(neighc))
               xrv_nd = fid*rvf_nd(kb)+fid1*rvf_nd(i)  
               kij = t15(neighc)*1.d-6
                vel_nd = vxyd/(aij*muvij+kij_tol)  
                call nd_flow_vel(3,icd,1,iq,vxyd,vel_nd,aij,kij,
     &               fid,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd,i,kb,jm)
                aij = abs(t5(neighc))
                vxyd_nd = vel_nd*aij*muvij
                call nd_flow_vel(4,icd,1,iq,vxyd,vel_nd,aij,kij,
     &               fid,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd,i,kb,jm)

c find upwinding
                t8_nd(neighc)=vxyd_nd 
                 fid=0.5
                 if(vxyd_nd.lt.0.0) fid=dnwgt
                 if(vxyd_nd.gt.0.0) fid=upwgt
                 t9_nd(neighc)=fid  
c gaz 103024 need derivatived muij
             t13(neighc) = dvapi_nd*aij*muvij + vel_nd*aij*dmuvijpi
             t14(neighc) = dvapkb_nd*aij*muvij + vel_nd*aij*dmuvijpj
c find 
             continue   
           enddo       
c           endif

c     
c     determine upwind nodes and if vapour phase exists
c  
c  assumption- upwind direction for vdrift is the same as vtotal
c  can be different that liquid velocity   
c   
            isl=1
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

c gaz 110424
c gaz 091624 aij = t5() area
c   gaz db 010125            aij=abs(t5(neighc)
                aij = abs(t5(neighc))
                aij = 1.d0   
                vxyd_nd = t8_nd(neighc)  
                fid=t9_nd(neighc)
                fid1=1.0-fid 
                sx4h=t7(neighc)
                divkb=div(kb)
                vxyf=(fid*divkb+fid1*divi)
                vxy_nd = vxyd_nd*vxyf
                dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
                dvekb=0.5*sx4h*dgve(kb)
     &                 *(cord(kz,igrav)-cord(iz,igrav))
                divkb=div(kb)
                divpkb=divp(kb)
                divekb=dive(kb)                
                dvapi_nd = t13(neighc)
                dvapkb_nd = t14(neighc)
                dvapi_nd = dvapi_nd*vxyf+vxyd_nd*fid1*divpi
                dvapkb_nd = dvapkb_nd*vxyf+vxyd_nd*fid*divpkb
c gaz 111124
                dvaei_nd=dvei*vxyf+vxyd_nd*fid1*divei
                dvaekb_nd=dvekb*vxyf+vxyd_nd*fid*divekb
                continue
c gaz_testing  121424
       if(nd_test_write.gt.0) then
c        if(i.eq.1.and.kb.eq.21.and.iad.eq.0) then
c         write(55,181)l,i,kb,vxy,dvapi,dvapkb,
c     &    dvaei,dvaekb,vxy_nd,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd
c       else if(i.eq.360.and.kb.eq.380.and.iad.eq.0) then
c         write(56,181)l,i,kb,vxy,dvapi,dvapkb,
c     &    dvaei,dvaekb,vxy_nd,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd
c       else if(i.eq.365.and.kb.eq.385.and.iad.eq.0) then
c         write(57,181)l,i,kb,vxy,dvapi,dvapkb,
c     &    dvaei,dvaekb,vxy_nd,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd
c       endif
181    format(t1,i5,t9,i6,t17,i6,1p,5(1x,g16.6),/,t23,5(1x,g16.6))
c gaz 230125 moved endif up (consistent with 2D ?)
       endif
                vxy = vxy_nd
                dvapi = dvapi_nd
                dvapkb = dvapkb_nd
                dvaei = dvaei_nd
                dvaekb = dvaekb_nd
c       endif                           
       
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
c gaz 102914 
            t15(neighc)= (pxy)/(sx2c+kij_tol)
 69      continue
         if(irdof.ne.11) then
c     
c liquid phase calculations
c
c gaz 091324 if(nd_flow) then
            do jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxy = t3(neighc)
               pxyi=t1(neighc)
               sx4d=t6(neighc)
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
c axyd units m**2*(area/dis)*Mpa
               g_term = 0.5*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
               dg_termpi = 0.5*sx4d*dglp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
               dg_termpkb = 0.5*sx4d*dglp(kb)*
     &                (cord(kz,igrav)-cord(iz,igrav))
               daxydpi = -pxy+dg_termpi
               daxydpkb = pxy+dg_termpkb
               if(irdof.ne.13) then
                daxydei=pxy*dpvti+0.5*sx4d*dgle(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
                daxydekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav))
               endif
               t8(neighc)=axyd
c   
                 fid=0.5
c gaz 220125 removed upwind for liq
                 if(axyd.lt.0.0) fid=dnwgt
                 if(axyd.gt.0.0) fid=upwgt
                 fid1= 1.d0-fid
                aij = abs(t5(neighc)) 
                xrl_nd = fid*rlf_nd(kb)+fid1*rlf_nd(i)
                kij = t15(neighc)*1.d-6*xrl_nd
                call nd_props(0,icd,1,iq,i,kb,jm,0.5d0)
                vel_nd = axyd/(aij*muij+kij_tol)
                call nd_flow_vel(1,icd,1,iq,axyd,vel_nd,aij,kij,
     &                fid,dlapi_nd,dlapkb_nd,dlaei_nd,dlaekb_nd,i,kb,jm) 
                 aij = abs(t5(neighc)) 
c gaz 090125
                 axyd_nd = vel_nd*muij*aij

                call nd_flow_vel(2,icd,1,iq,axyd,vel_nd,aij,kij,
     &                fid,dlapi_nd,dlapkb_nd,dlaei_nd,dlaekb_nd,i,kb,jm)  
    
c determine upwind direction
c find upwinding
                t8_nd(neighc)=axyd_nd 
                 fid=0.5
                 if(axyd_nd.lt.0.0) fid=dnwgt
                 if(axyd_nd.gt.0.0) fid=upwgt
                 t9_nd(neighc)=fid  
c gaz 103024 need derivatived muij
             aij = abs(t5(neighc))  
             t13(neighc) = dlapi_nd*aij*muij + vel_nd*aij*dmuijpi
             t14(neighc) = dlapkb_nd*aij*muij + vel_nd*aij*dmuijpj
             t18(neighc) = dlaei_nd*aij*muij + vel_nd*aij*dmuijpi
             t19(neighc) = dlaekb_nd*aij*muij + vel_nd*aij*dmuijpj
             continue   
           enddo       
   
c     liquid phase calculations
c     
c            do 70 jm=1,iq
c               kb=it8(jm)
c               kz=kb-icd
c               neighc=it9(jm)
c               pxyi=t1(neighc)
c               sx4d=t6(neighc)
c               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
c     *              *(cord(kz,igrav)-cord(iz,igrav))
c               t8(neighc)=axyd
c 70         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
c            isl=1
c            do 71 jm=1,iq
c               kb=it8(jm)
c               kz=kb-icd
c               neighc=it9(jm)
c add coding to save upwind position
cc               if(iad.le.iad_up) then
c                  fid=0.5
c                  axyd=t8(neighc)
c                  if(axyd.lt.0.0) fid=dnwgt
c                  if(axyd.gt.0.0) fid=upwgt
c                  t9(neighc)=fid
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
c 71         continue
c     
c     form equations
c  
c
            isl=1   
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
                  pxyi=t1(neighc)
                  pxy=t3(neighc)
                  sx4d=t6(neighc)
                  dilkb=dil(kb)
                  dilpkb=dilp(kb)

c gaz 091624 aij = t5() area
c             if(nd_flow) then
c   gaz db 010125            aij=abs(t5(neighc)
c                aij = abs(t5(neighc))  
                aij = 1.0d0 
                axyd_nd = t8_nd(neighc)
                fid=t9_nd(neighc)
                fid1=1.0-fid  
                axyf=(fid*dilkb+fid1*dili)
                axy_nd = axyd_nd*axyf
                
                dlapi_nd = t13(neighc)
                dlapkb_nd = t14(neighc)
                dlapi_nd = dlapi_nd*axyf+axyd_nd*fid1*dilpi
                dlapkb_nd = dlapkb_nd*axyf+axyd_nd*fid*dilpkb
                if(irdof.ne.13) then
                  dilekb=dile(kb)
                  dlaei_nd = t18(neighc)
                  dlaekb_nd = t19(neighc)
                  dlaei_nd=dlaei_nd*axyf+axyd_nd*fid1*dilei
                  dlaekb_nd=dlaekb_nd*axyf+axyd_nd*fid*dilekb    
                endif
                continue
c
                axy = axy_nd
                dlapi = dlapi_nd
                dlapkb = dlapkb_nd
                dlaei = dlaei_nd
                dlaekb = dlaekb_nd
 
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
            do jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               pvxy = t4(neighc)
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
                 fid=0.5
c                 if(vxyd.lt.0.0) fid=dnwgt
c                 if(vxyd.gt.0.0) fid=upwgt
                 fid1 = 1.0d0-fid
c vxyd units m**2*(area/dis)*Mpa
c   gaz db 010125            aij=abs(t5(neighc)
c               aij = abs(t5(neighc))
c               kij = t15(neighc)*1.d-6
                  divkb=div(kb)
                  divpkb=divp(kb)
                  divekb=dive(kb) 
                  g_term = 0.5*sx4h*(rovf(i)+rovf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))                
                  dg_termpi = 0.5*sx4h*dgvp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
                  dg_termpkb = 0.5*sx4h*dgvp(kb)*
     &                 (cord(kz,igrav)-cord(iz,igrav))

               daxydpi = -pvxy+dg_termpi
               daxydpkb = pvxy+dg_termpkb

               daxydei= 0.5*sx4d*dgve(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
               daxydekb= 0.5*sx4d*dgve(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav))
               call nd_props(0,icd,1,iq,i,kb,jm,0.5d0)
c gaz 112424                
c vxyd is the darcy based velocity
c                fid = 0.5d0
c gaz 220525 added upwind calc
c                if(vxyd.lt.0.0) fid=dnwgt
c                if(vxyd.gt.0.0) fid=upwgt
c                fid1= 1.d0-fid
                aij = abs(t5(neighc)) 
                xrv_nd = fid*rvf_nd(kb)+fid1*rvf_nd(i)                              
                kij = t15(neighc)*1.d-6*xrv_nd
                call nd_props(0,icd,1,iq,i,kb,jm,0.5d0)
                vel_nd = vxyd/(aij*muvij+kij_tol)           
                call nd_flow_vel(3,icd,1,iq,vxyd,vel_nd,aij,kij,
     &               fid,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd,i,kb,jm)
c gaz 050125
                aij = abs(t5(neighc))
                vxyd_nd = vel_nd*aij*muvij
                call nd_flow_vel(4,icd,1,iq,vxyd,vel_nd,aij,kij,
     &               fid,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd,i,kb,jm)

c find upwinding
                t8_nd(neighc)=vxyd_nd 
                 fid=0.5
                 if(vxyd_nd.lt.0.0) fid=dnwgt
                 if(vxyd_nd.gt.0.0) fid=upwgt
                 fid1 = 1.0d0-fid
                 t9_nd(neighc)=fid  
c gaz 103024 need derivatived muij
               aij = abs(t5(neighc))
             t13(neighc) = dvapi_nd*aij*muvij + vel_nd*aij*dmuvijpi
             t14(neighc) = dvapkb_nd*aij*muvij + vel_nd*aij*dmuvijpj
c find 
             continue   
           enddo       
c           endif
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
c   gaz db 010125            aij=abs(t5(neighc)
                aij = abs(t5(neighc))  
                aij = 1.d0 
                vxyd_nd = t8_nd(neighc)  
                fid=t9_nd(neighc)
                sx4h=t7(neighc)
                fid1=1.0-fid 
                divkb=div(kb)
                vxyf=(fid*divkb+fid1*divi)
                vxy_nd = vxyd_nd*vxyf
                dvei=0.5*sx4h*dgve(i)*(cord(kz,igrav)-cord(iz,igrav))
                dvekb=0.5*sx4h*dgve(kb)
     &                 *(cord(kz,igrav)-cord(iz,igrav))
                divkb=div(kb)
                divpkb=divp(kb)
                divekb=dive(kb)                
                dvapi_nd = t13(neighc)
                dvapkb_nd = t14(neighc)
                dvapi_nd = dvapi_nd*vxyf+vxyd_nd*fid1*divpi
                dvapkb_nd = dvapkb_nd*vxyf+vxyd_nd*fid*divpkb
c gaz 111124
                dvaei_nd=dvei*vxyf+vxyd_nd*fid1*divei
                dvaekb_nd=dvekb*vxyf+vxyd_nd*fid*divekb
                continue

                vxy = vxy_nd
                dvapi = dvapi_nd
                dvapkb = dvapkb_nd
                dvaei = dvaei_nd
                dvaekb = dvaekb_nd

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
 

 
      