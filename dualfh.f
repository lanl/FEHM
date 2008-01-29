      subroutine  dualfh
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
CD1 To compute the Jacobian and residual terms of the heat and mass
CD1 equations at each node for a dual porosity solution.
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
CD2 $Log:   /pvcs.config/fehm90/src/dualfh.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Mon Jan 29 15:37:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   08/18/95 10:19:38   llt
CD2 neigh was already defined, removed for cray
CD2 
CD2    Rev 1.2   12/09/94 16:02:22   llt
CD2 Changes to fix  radi (made by gaz).
CD2 
CD2    Rev 1.1   03/18/94 15:50:10   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:28   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 None
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
CD4 lenreal, nn, neq, nelmdg, sx1, apuv1, tmch, phi, pcp, 
CD4 pnx, pny, pnz, dil, div, s,  rolf, rovf,
CD4 igrav,
CD4 cord, dnwgt, upwgt, aw, ay, a, bp,
CD4 wb11, tb11, a21mpf, a32mpf, rb2mf, r3mf,
CD4  t1, t2, t3, t4, t5, t6, t7, t8, t9, zero_t
CD4 
CD4 
CD4 ipb, ipt1, ipt2, ipt3,ipt4,  ipt5, ipt6, ipt7, ipt8, ipt9, ipt10,
CD4 ipit8, ipit9, ipit10, ipit11, ipit12, istrw, nelm, dpcef, enlf,
CD4 delf, delef, envf, devf, devef, dilp, dile, divp, dive, thx, thy,
CD4 thz, t, s, perml, permv, dglp, dgle, devef, dgvp, dgve, dtpa,
CD4 dtpae, sk, deni, denj, qh, denei, denej, dq, dmpf, dqt, dmef, 
CD4 deef, deqh, wb12, wb21, wb22, tb12, tb21, tb22, a21mef, a21epf,
CD4 a21eef, a32mef, a32epf, a32eef, rb2ef, r3ef, nbd
CD4 
CD4 Global Subprograms
CD4
CD4 Name       Type     Description
CD4 
CD4 mmgetblk   N/A      Allocates memory for an array
CD4 mmrelblk   N/A      Frees array space
CD4 thermc     N/A      Computes solute storage and reaction terms
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
CD5 pvii         real*8      Pressure of vapor
CD5 phii         real*8      Pressure of liquid
CD5 swi          real*8      Saturation
CD5 dili         real*8      Liquid transmissibility
CD5 dilkb        real*8      Liquid transmissibility, connected node
CD5 divi         real*8      Vapor transmissibility
CD5 divkb        real*8      Vapor transmissibility, connected node
CD5 id           int         Integer index parameter
CD5 idg          int         Integer index parameter
CD5 jm           int        Integer index parameter
CD5 neqp1        int         Number of equations plus 1
CD5 neigh        int         Integer index parameter
CD5 neighc       int         Integer index parameter
CD5 axkb         real*8      Permeability in x-direction, connected node
CD5 sx2c         real*8      Real parameter used in calculation
CD5 sx4d         real*8      Real parameter used in calculation
CD5 sx4h         real*8      Real parameter used in calculation
CD5 pvikb        real*8      Vapor pressure at connected node
CD5 phikb        real*8      Liquid pressure at connected node
CD5 radi         real*8      Parameter used in calculation
CD5 radkb        real*8      Parameter used in calculation
CD5 fid          real*8      Parameter used in calculation
CD5 axyd         real*8      Parameter used in calculation
CD5 axyf         real*8      Parameter used in calculation
CD5 axy          real*8      Parameter used in calculation
CD5 aexyf        real*8      Parameter used in calculation
CD5 vxy          real*8      Parameter used in calculation
CD5 vexy         real*8      Parameter used in calculation
CD5 vxyf         real*8      Parameter used in calculation
CD5 vexyf        real*8      Parameter used in calculation
CD5 pxyh         real*8      Derivative term used in calculation
CD5 heatc        real*8      Derivative term used in calculation
CD5 vxyd         real*8      Derivative term used in calculation
CD5 icode        int         Returned flag from call to mmgetblk
CD5 neq2         int         Twice the number of equations
CD5 idum         int         Index for arrays
CD5 tot          real*8      Total volume at a node
CD5 idl          int         Unknown number of first matrix node
CD5 idll         int         Unknown number of second matrix node
CD5 frac0        real*8      Volume fraction of primary node
CD5 frac1        real*8      Volume fraction of first matrix node
CD5 frac2        real*8      Volume fraction of second matrix node
CD5 alen         real*8      Volume to area ratio
CD5 area         real*8      Area term used in calculation
CD5 al0          real*8      Length scale for primary node
CD5 al1          real*8      Length scale for first matrix node
CD5 al2          real*8      Length scale for second matrix node
CD5 dist01       real*8      Average length scale used in calculation
CD5 dist12       real*8      Average length scale of matrix nodes
CD5 sx1d0        real*8      Volume of primary node
CD5 sx1d1        real*8      Volume of first matrix node
CD5 sx1d2        real*8      Volume of second matrix node
CD5 coef1        real*8      Length scale used in calculation
CD5 coef2        real*8      Length scale used in calculation
CD5 i            int         Index for first matrix node number
CD5 kb           int         Index for second matrix node number
CD5 iz           int         Index for primary node number
CD5 kz           int         Index for primary node number
CD5 sx2t         real*8      Parameter used in calculation
CD5 pxy          real*8      Parameter used in calculation
CD5 pvxy         real*8      Parameter used in calculation
CD5 pxyi         real*8      Parameter used in calculation
CD5 a1           real*8      Argument to statement function definition
CD5 a2           real*8      Argument to statement function definition
CD5 a3           real*8      Argument to statement function definition
CD5 a4           real*8      Argument to statement function definition
CD5 b1           real*8      Argument to statement function definition
CD5 b2           real*8      Argument to statement function definition
CD5 iwa          int         Integer index used in calculation
CD5 nb0          int         Integer index used in calculation
CD5 nb1          int         Integer index used in calculation
CD5 nb2          int         Integer index used in calculation
CD5 r1m          real*8      Term used in correction of residual of
CD5                             mass balance in matrix nodes
CD5 r2m          real*8      Term used in correction of residual of
CD5                             mass balance in matrix nodes
CD5 r3m          real*8      Term used in correction of residual of
CD5                             mass balance in matrix nodes
CD5 r1e          real*8      Term used in correction of residual of
CD5                             energy balance in matrix nodes
CD5 r2e          real*8      Term used in correction of residual of
CD5                             energy balance in matrix nodes
CD5 r3e          real*8      Term used in correction of residual of
CD5                             energy balance in matrix nodes
CD5 a11mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a12mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a21mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a22mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a23mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a32mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a33mp        real*8      Derivative term used in mass balance
CD5                             (mass with respect to pressure)
CD5 a11me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a12me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a21me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a22me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a23me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a32me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a33me        real*8      Derivative term used in mass balance
CD5                             (mass with respect to energy)
CD5 a11ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a12ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a21ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a22ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a23ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a32ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a33ep        real*8      Derivative term used in energy balance
CD5                             (energy with respect to pressure)
CD5 a11ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a12ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a21ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a22ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a23ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a32ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 a33ee        real*8      Derivative term used in energy balance
CD5                             (energy with respect to energy)
CD5 dpvti        real*8      Cap. pressure derivative term
CD5 enli         real*8      Liquid enthalpy
CD5 enlkb        real*8      Liquid enthalpy, connected node
CD5 envkb        real*8      Vapor enthalpy, connected node
CD5 deli         real*8      Enthalpy derivative with pressure
CD5 delkb        real*8      Enthalpy derivative with pressure,
CD5                             connected node
CD5 devkb        real*8      Enthalpy derivative with pressure (vapor),
CD5                             connected node
CD5 delei        real*8      Enthalpy derivative with energy
CD5 delekb       real*8      Enthalpy derivative with energy,
CD5                             connected node
CD5 devekb       real*8      Enthalpy derivative with energy (vapor),
CD5                             connected node
CD5 envi         real*8      Vapor enthalpy
CD5 devi         real*8      Enthalpy derivative with pressure (vapor)
CD5 devei        real*8      Enthalpy derivative with energy (vapor)
CD5 dilpi        real*8      Trans. derivative with pressure
CD5 dilpkb       real*8      Trans. derivative with pressure,
CD5                             connected node
CD5 dilei       real*8      Trans. derivative with enthalpy
CD5 dilekb      real*8      Trans. derivative with enthalpy, connected
CD5                             node
CD5 divpi        real*8      Trans. derivative with pressure (vapor)
CD5 divpkb       real*8      Trans. derivative with pressure (vapor),
CD5                             connected node
CD5 divei        real*8      Trans. derivative with enthalpy (vapor)
CD5 divekb       real*8      Trans. derivative with enthalpy (vapor),
CD5                             connected node
CD5 thxi         real*8      Thermal conductivity, x-direction
CD5 thyi         real*8      Thermal conductivity, y-direction
CD5 thzi         real*8      Thermal conductivity, z-direction
CD5 ti           real*8      Temperature
CD5 athkb        real*8      Thermal conductivity, matrix nodes
CD5 thxkb        real*8      Thermal conductivity, matrix nodes
CD5 isl          int         Flag denoting if liquid or vapor exists
CD5 dlpi         real*8      Derivative term used in calculation
CD5 dvpi         real*8      Derivative term used in calculation
CD5 dvpkb        real*8      Derivative term used in calculation
CD5 dlpkb        real*8      Derivative term used in calculation
CD5 dlei         real*8      Derivative term used in calculation
CD5 dvei         real*8      Derivative term used in calculation
CD5 dlekb        real*8      Derivative term used in calculation
CD5 dvekb        real*8      Derivative term used in calculation
CD5 aexy         real*8      Term used in calculation
CD5 dlapi        real*8      Derivative term used in calculation
CD5 dlapkb       real*8      Derivative term used in calculation
CD5 dlaei        real*8      Derivative term used in calculation
CD5 dlaekb       real*8      Derivative term used in calculation
CD5 dlepi        real*8      Derivative term used in calculation
CD5 dlepkb       real*8      Derivative term used in calculation
CD5 dleei        real*8      Derivative term used in calculation
CD5 dleekb       real*8      Derivative term used in calculation
CD5 dvapi        real*8      Derivative term used in calculation
CD5 dvaei        real*8      Derivative term used in calculation
CD5 dvapkb       real*8      Derivative term used in calculation
CD5 dvepi        real*8      Derivative term used in calculation
CD5 dvepkb       real*8      Derivative term used in calculation
CD5 dveei        real*8      Derivative term used in calculation
CD5 dveekb       real*8      Derivative term used in calculation
CD5 dvaekb       real*8      Derivative term used in calculation
CD5 t11          real*8      Intermediate term used in calculation
CD5 t12          real*8      Intermediate term used in calculation
CD5 t21          real*8      Intermediate term used in calculation
CD5 t22          real*8      Intermediate term used in calculation
CD5 u11          real*8      Intermediate term used in calculation
CD5 u12          real*8      Intermediate term used in calculation
CD5 u21          real*8      Intermediate term used in calculation
CD5 u22          real*8      Intermediate term used in calculation
CD5 v11          real*8      Intermediate term used in calculation
CD5 v12          real*8      Intermediate term used in calculation
CD5 v21          real*8      Intermediate term used in calculation
CD5 v22          real*8      Intermediate term used in calculation
CD5 ab22mp       real*8      Intermediate term used in calculation
CD5 ab22me       real*8      Intermediate term used in calculation
CD5 ab22ep       real*8      Intermediate term used in calculation
CD5 ab22ee       real*8      Intermediate term used in calculation
CD5 w11          real*8      Intermediate term used in calculation
CD5 w12          real*8      Intermediate term used in calculation
CD5 w21          real*8      Intermediate term used in calculation
CD5 w22          real*8      Intermediate term used in calculation
CD5 x11          real*8      Intermediate term used in calculation
CD5 x12          real*8      Intermediate term used in calculation
CD5 x21          real*8      Intermediate term used in calculation
CD5 x22          real*8      Intermediate term used in calculation
CD5 y11          real*8      Intermediate term used in calculation
CD5 y12          real*8      Intermediate term used in calculation
CD5 y21          real*8      Intermediate term used in calculation
CD5 y22          real*8      Intermediate term used in calculation
CD5 rb2m         real*8      Term used in calculation
CD5 rb2e         real*8      Term used in calculation
CD5 
CD5 Local Subprograms
CD5 
CD5 Identifier    Type      Description
CD5 
CD5 ali           real*8    Function used in matrix calculation
CD5 alm           real*8    Function used in matrix calculation
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
CD9 2.4.7 Dual-porosity formulation
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
CPS BEGIN dualfh
CPS 
CPS mmrelblk - deallocate array space not needed here
CPS mmgetblk - allocate space in temporary arrays
CPS 
CPS varchk - determine equation of state parameters for first set of...
CPS ... matrix nodes
CPS varchk - determine equation of state parameters for second set of...
CPS ... matrix nodes
CPS 
CPS Compute integer parameters needed later
CPS 
CPS FOR each node
CPS 
CPS   Initialize parameters for computation of heat and mass...
CPS   ... calculations at the matrix nodes
CPS   
CPS   IF there is no volume in the second set of matrix nodes
CPS     Set parameters accordingly
CPS   ENDIF
CPS   
CPS   Compute terms used later in the heat and mass balance...
CPS   ... calculations for the matrix nodes
CPS   
CPS   IF there is liquid at this node
CPS     Compute liquid mass and energy balance contributions for...
CPS     ... liquid in the matrix nodes
CPS   ENDIF
CPS   
CPS   IF there is vapor at this node
CPS     Compute vapor mass and energy balance contributions for...
CPS     ... vapor in the matrix nodes
CPS   ENDIF
CPS 
CPS   Compute and add heat conduction terms to energy balance equation
CPS   Add accumulation term for first and second matrix nodes
CPS   
CPS   Compute terms used later in the heat and mass balance...
CPS   ... calculations for the primary nodes
CPS   
CPS   IF there is liquid at this node
CPS     Compute liquid mass and energy balance contributions for...
CPS     ... liquid in the primary nodes
CPS   ENDIF
CPS   
CPS   IF there is vapor at this node
CPS     Compute vapor mass and energy balance contributions for...
CPS     ... vapor in the primary nodes
CPS   ENDIF
CPS   
CPS   Compute and add heat conduction terms to energy balance equation
CPS   
CPS   Compute terms needed in matrix for computing heat and mass in...
CPS   ... the matrix nodes
CPS   
CPS   Add contributions of matrix nodes to derivative terms
CPS   Add contributions of matrix nodes to residual values
CPS   Save information needed to back-substitute to obtain values in...
CPS   ... matrix nodes
CPS   
CPS ENDFOR
CPS 
CPS mmrelblk - deallocate array space not needed here
CPS mmgetblk - allocate space in temporary arrays
CPS 
CPS END dualfh
CPS 
C**********************************************************************

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

      real*8 sx1d,pvii,phii,swi,dili,dilkb,divi,divkb,axkb,sx2c,sx4d
      real*8 sx4h,pvikb,phikb,radi,radkb,fid,axyd,axyf,aexyf,vxy
      real*8 vexy,vxyf,vexyf,pxyh,heatc,vxyd,tot,frac0,frac1,frac2
      real*8 alen,area,al0,al1,al2,dist01,dist12,sx1d0,sx1d1,sx1d2
      real*8 coef1,coef2,sx2t,pxy,pvxy,pxyi,a1,a2,a3,a4,b1,b2,r1m,r2m
      real*8 r3m,r1e,r2e,r3e,a11mp,a12mp,a21mp,a22mp,a23mp,a32mp,a33mp
      real*8 a11me,a12me,a21me,a22me,a23me,a32me,a33me,a11ep,a12ep
      real*8 a21ep,a22ep,a23ep,a32ep,a33ep,a11ee,a12ee,a21ee,a22ee
      real*8 a23ee,a32ee,a33ee,dpvti,enli,enlkb,envkb,deli,delkb,devkb
      real*8 delei,delekb,devekb,envi,devi,devei,dilpi,dilpkb,dilei
      real*8 dilekb,divpi,divpkb,divei,divekb,thxi,thyi,thzi,ti,athkb
      real*8 thxkb,dlpi,dvpi,dlpkb,dlei,dvei,dlekb,dvekb,aexy,dlapi
      real*8 dlapkb,dlaei,dlaekb,dlepi,dlepkb,dleei,dleekb,dvapi
      real*8 dvapkb,dvepi,dvepkb,dveei,dveekb,t11,t12,t21,t22,u11
      real*8 u12,u21,u22,v11,v12,v21,v22,ab22mp,ab22me,ab22ep,ab22ee
      real*8 w11,w12,w21,w22,x11,x12,x21,x22,y11,y12,y21,y22
      real*8 fid1,axy,dvpkb,dvaekb,dvaei,rb2m,rb2e,ali,alm
      integer id,idg,jm,neqp1,neighc,icode,neq2,idum,idl,idll,i
      integer kb,iz,kz,iwa,nb0,nb1,nb2,isl

c     
c     compute the fracture-matrix transfer terms
c     
c**** linear algebra ****
c     
      alm(a1,a2,b1,b2)    =  a1*b1+a2*b2
      ali(a1,a2,a3,a4,b1) =  b1/( a1*a4-a2*a3 )
      
      neq2   =  neq+neq
      neqp1  =  neq+1
c     
c**** call thermo properties of matrix material ****
c     
      call  varchk  ( 0,neq   )
      call  varchk  ( 0,neq*2 )
      
      iwa    =  istrw(neq)+neigh
      neq2   =  neq+neq
      neqp1  =  neq+1
      nb0     =  nelm(neqp1)-neqp1
      nb1    =  nb0 +nb0
      nb2    =  nb1+nb0
c     
c     modify geometric terms
c     save nelm terms for later use
c     
c     loop on nodes
      do id=1,neq
c     identify diagonal member
         idg =  nelmdg(id)
         idum =  idg-neqp1
c     
c     form equations at id+neq
c     
         idl =  id+neq
         idll =  id+neq2
         r1m=0.0
         r2m=0.0
         r3m=0.0
         r1e=0.0
         r2e=0.0
         r3e=0.0
         a11mp=0.0
         a12mp=0.0
         a21mp=0.0
         a22mp=0.0
         a23mp=0.0
         a32mp=0.0
         a33mp=0.0
         a11me=0.0
         a12me=0.0
         a21me=0.0
         a22me=0.0
         a23me=0.0
         a32me=0.0
         a33me=0.0
         a11ep=0.0
         a12ep=0.0
         a21ep=0.0
         a22ep=0.0
         a23ep=0.0
         a32ep=0.0
         a33ep=0.0
         a11ee=0.0
         a12ee=0.0
         a21ee=0.0
         a22ee=0.0
         a23ee=0.0
         a32ee=0.0
         a33ee=0.0
c     
c     generate transport terms
c     
         tot =  sx1(id)+sx1(idl)+sx1(idll)
         frac0 =  sx1(id)/tot
         frac1 =  sx1(idl)/tot
         frac2 =  sx1(idll)/tot
         alen =  apuv1(id)
         area=tot/alen
         al0 =  frac0*alen
         al1 =  frac1*alen
         al2 =  frac2*alen
         dist12=0.5*(al1+al2)
         coef2=-area/dist12
         sx1d0=sx1(id)
         sx1d1=sx1(idl)
         sx1d2=sx1(idll)
         if(frac2.le.tmch) then
            coef2=0.0
            sx1d2=1.0
         endif

         i=idl
         kb=idll
         iz=id
         kz=id

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
c     form constants for i>neq
c     
         axkb =  max( pnx(idll  ),pny(idll  ),
     *        pnz(idll  ),zero_t )
         athkb =  max( thx(idll  ),thy(idll  ),
     *        thz(idll  ),zero_t )
c     
c     2-d geometry
c     
         if(icnl.ne.0) then
            radi=cord(iz,3)
         else
            radi = 1.
         endif
         jm=1
         neighc=1
         perml(1)=axkb
         permv(1)=perml(1)
         radkb=radi
         sx2c=radkb*coef2
         thxkb=athkb
         sx2t=sx2c*thxkb
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

c     liquid phase calculations
c     
         jm=1
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
            neighc=1
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
            dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
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

            r2m=r2m+axy
            r3m=r3m-axy
            a22mp=a22mp+dlapi
            a22me=a22me+dlaei
            a32mp=a32mp-dlapi
            a32me=a32me-dlaei
            a23mp=a23mp+dlapkb
            a23me=a23me+dlaekb
            a33mp=a33mp-dlapkb
            a33me=a33me-dlaekb

            r2e=r2e+aexy
            r3e=r3e-aexy
            a22ep=a22ep+dlepi
            a22ee=a22ee+dleei
            a32ep=a32ep-dlepi
            a32ee=a32ee-dleei
            a23ep=a23ep+dlepkb
            a23ee=a23ee+dleekb
            a33ep=a33ep-dlepkb
            a33ee=a33ee-dleekb
         endif
c     
c     vapour phase calculations
         jm=1
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
            neighc=1
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
     *           *(cord(kz,igrav)-cord(iz,igrav))
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

            r2m=r2m+vxy
            r3m=r3m-vxy
            a22mp=a22mp+dvapi
            a22me=a22me+dvaei
            a32mp=a32mp-dvapi
            a32me=a32me-dvaei
            a23mp=a23mp+dvapkb
            a23me=a23me+dvaekb
            a33mp=a33mp-dvapkb
            a33me=a33me-dvaekb

            r2e=r2e+vexy
            r3e=r3e-vexy
            a22ep=a22ep+dvepi
            a22ee=a22ee+dveei
            a32ep=a32ep-dvepi
            a32ee=a32ee-dveei
            a23ep=a23ep+dvepkb
            a23ee=a23ee+dveekb
            a33ep=a33ep-dvepkb
            a33ee=a33ee-dveekb
         endif
c     
c     add heat conduction
c     
         jm=1
         
         neighc=1
         heatc=t5(neighc)
         vxy=heatc*(t(kb)-ti)
         dvapi=-heatc*dtpa(i)
         dvaei=-heatc*dtpae(i)
         dvapkb=heatc*dtpa(kb)
         dvaekb=heatc*dtpae(kb)
         r2e=r2e+vxy
         r3e=r3e-vxy
         a22ep=a22ep+dvapi
         a22ee=a22ee+dvaei
         a32ep=a32ep-dvapi
         a32ee=a32ee-dvaei
         a23ep=a23ep+dvapkb
         a23ee=a23ee+dvaekb
         a33ep=a33ep-dvapkb
         a33ee=a33ee-dvaekb
c     
c     add accumulation terms for second and first matrix levels
c     
         sx1d=sx1d1
c     r2m=r2m+sx1d*(aw*deni(i)+ay*denj(i))+sk(i)
c     r2e=r2e+sx1d*(aw*denei(i)+ay*denej(i))+qh(i)
      r2m=r2m+sx1d*deni(i)+sk(i)
      r2e=r2e+sx1d*denei(i)+qh(i)
c     a22mp=a22mp+sx1d*(aw*dmpf(i))+dq(i)
c     a22me=a22me+sx1d*(aw*dmef(i))+dqt(i)
c     a22ep=a22ep+sx1d*(depf(i)*aw)+dqh(i)
c     a22ee=a22ee+sx1d*(aw*deef(i))+deqh(i)
      a22mp=a22mp+sx1d*dmpf(i)+dq(i)
      a22me=a22me+sx1d*dmef(i)+dqt(i)
      a22ep=a22ep+sx1d*depf(i)+dqh(i)
      a22ee=a22ee+sx1d*deef(i)+deqh(i)

         sx1d=sx1d2
c     r3m=r3m+sx1d*(aw*deni(kb)+ay*denj(kb))+sk(kb)
c     r3e=r3e+sx1d*(aw*denei(kb)+ay*denej(kb))+qh(kb)
      r3m=r3m+sx1d*deni(kb)+sk(kb)
      r3e=r3e+sx1d*denei(kb)+qh(kb)
c     a33mp=a33mp+sx1d*(aw*dmpf(kb))+dq(kb)
c     a33me=a33me+sx1d*(aw*dmef(kb))+dqt(kb)
c     a33ep=a33ep+sx1d*(depf(kb)*aw)+dqh(kb)
c     a33ee=a33ee+sx1d*(aw*deef(kb))+deqh(kb)
      a33mp=a33mp+sx1d*dmpf(kb)+dq(kb)
      a33me=a33me+sx1d*dmef(kb)+dqt(kb)
      a33ep=a33ep+sx1d*depf(kb)+dqh(kb)
      a33ee=a33ee+sx1d*deef(kb)+deqh(kb)
c     
c     form equations at node id
c     
c     geometry
c     
         dist01=0.5*(al0+al1)
         coef1=-area/dist01
         if(frac1.le.0.0) then
            coef1=0.0
         endif

         i=id
         kb=idl
         iz=id
         kz=id

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
c     form constants for i>neq
c     
c     
c     coef1=-area/dist01
         axkb =  max( pnx(idl  ),pny(idl  ),
     *        pnz(idl  ),zero_t )
         athkb =  max( thx(idl  ),thy(idl  ),
     *        thz(idl  ),zero_t )
c     
c     2-d geometry
c     
         if(icnl.ne.0) then
            radi=cord(iz,3)
         else
            radi = 1.
         endif
         jm=1
         neighc=1
         perml(1)=axkb
         permv(1)=perml(1)
         radkb=radi
         sx2c=radkb*coef1
         thxkb=athkb
         sx2t=sx2c*thxkb
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
            neighc=1
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
            dlei=pxy*dpvti+0.5*sx4d*dgle(i)*(cord(kz,igrav)-
     &           cord(iz,igrav))
            dlekb=-pxy*dpcef(kb)+0.5*sx4d*dgle(kb)
     *           *(cord(kz,igrav)-cord(iz,igrav))
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

            r1m=r1m+axy
            r2m=r2m-axy
            a11mp=a11mp+dlapi
            a11me=a11me+dlaei
            a21mp=a21mp-dlapi
            a21me=a21me-dlaei
            a12mp=a12mp+dlapkb
            a12me=a12me+dlaekb
            a22mp=a22mp-dlapkb
            a22me=a22me-dlaekb

            r1e=r1e+aexy
            r2e=r2e-aexy
            a11ep=a11ep+dlepi
            a11ee=a11ee+dleei
            a21ep=a21ep-dlepi
            a21ee=a21ee-dleei
            a12ep=a12ep+dlepkb
            a12ee=a12ee+dleekb
            a22ep=a22ep-dlepkb
            a22ee=a22ee-dleekb
         endif
c     
c     vapour phase calculations
c     
         jm=1
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
            neighc=1
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
     *           *(cord(kz,igrav)-cord(iz,igrav))
            vxyf=(fid*divkb+fid1*divi)
            vxy=vxyd*vxyf
            dvapi=dvpi*vxyf+vxyd*fid1*divpi
            dvapkb=dvpkb*vxyf+vxyd*fid*divpkb
            dvaei=dvei*vxyf+vxyd*fid1*divei
            dvaekb=dvekb*vxyf+vxyd*fid*divekb
            r1e=r1e+vxy
            r2e=r2e-vxy
            a11ep=a11ep+dvapi
            a11ee=a11ee+dvaei
            a21ep=a21ep-dvapi
            a21ee=a21ee-dvaei
            a12ep=a12ep+dvapkb
            a12ee=a12ee+dvaekb
            a22ep=a22ep-dvapkb
            a22ee=a22ee-dvaekb
         endif
c     
c     add heat conduction
c     
         jm=1
         neighc=1
         heatc=t5(neighc)
         heatc=t5(neighc)
         vxy=heatc*(t(kb)-ti)
         dvapi=-heatc*dtpa(i)
         dvaei=-heatc*dtpae(i)
         dvapkb=heatc*dtpa(kb)
         dvaekb=heatc*dtpae(kb)
         r1e=r1e+vxy
         r2e=r2e-vxy
         a11ep=a11ep+dvapi
         a11ee=a11ee+dvaei
         a21ep=a21ep-dvapi
         a21ee=a21ee-dvaei
         a12ep=a12ep+dvapkb
         a12ee=a12ee+dvaekb
         a22ep=a22ep-dvapkb
         a22ee=a22ee-dvaekb
c     
c     contribution to id are already done
c     
c     
c     form matrix ab22=a22-a23*a33i*a32
c     
c     form a33i=t
         t11 =  ali(a33mp,a33me,a33ep,a33ee,a33ee)
         t12 =  -ali(a33mp,a33me,a33ep,a33ee,a33me)
         t21 =  -ali(a33mp,a33me,a33ep,a33ee,a33ep)
         t22 =  ali(a33mp,a33me,a33ep,a33ee,a33mp)
c     form a23*a33i =  u
         u11 =  alm(a23mp,a23me,t11,t21)
         u12 =  alm(a23mp,a23me,t12,t22)
         u21 =  alm(a23ep,a23ee,t11,t21)
         u22 =  alm(a23ep,a23ee,t12,t22)
c     form a23*a33i*a32=v
         v11 =  alm(u11,u12,a32mp,a32ep)
         v12 =  alm(u11,u12,a32me,a32ee)
         v21 =  alm(u21,u22,a32mp,a32ep)
         v22 =  alm(u21,u22,a32me,a32ee)
c     form a22-a23*a33i*a32=ab22
         ab22mp =  a22mp-v11
         ab22me =  a22me-v12
         ab22ep =  a22ep-v21
         ab22ee =  a22ee-v22
c     form rb2=r2-a23*a33i*r3
         rb2m =  r2m-alm(u11,u12,r3m,r3e)
         rb2e =  r2e-alm(u21,u22,r3m,r3e)
c     
c     form matrix ab11=a11-a12*ab22i*a21
c     
c     form ab22i =  w
         w11 =   ali(ab22mp,ab22me,ab22ep,ab22ee,ab22ee)
         w12 =  -ali(ab22mp,ab22me,ab22ep,ab22ee,ab22me)
         w21 =  -ali(ab22mp,ab22me,ab22ep,ab22ee,ab22ep)
         w22 =   ali(ab22mp,ab22me,ab22ep,ab22ee,ab22mp)
c     form a12*ab22i=x
         x11 =  alm(a12mp,a12me,w11,w21)
         x12 =  alm(a12mp,a12me,w12,w22)
         x21 =  alm(a12ep,a12ee,w11,w21)
         x22 =  alm(a12ep,a12ee,w12,w22)
c     form a12*ab22i*a21=y
         y11 =  alm(x11,x12,a21mp,a21ep)
         y12 =  alm(x11,x12,a21me,a21ee)
         y21 =  alm(x21,x22,a21mp,a21ep)
         y22 =  alm(x21,x22,a21me,a21ee)
c     form a11-a12*a22i*a12
         a(idum    ) =  a(idum    )-y11+a11mp
         a(idum+nb0 ) =  a(idum+nb0 )-y12+a11me
         a(idum+nb1) =  a(idum+nb1)-y21+a11ep
         a(idum+nb2) =  a(idum+nb2)-y22+a11ee
c     form rb1=r1-a12*ab22i*rb2
         bp(id) =  bp(id)+r1m-alm(x11,x12,rb2m,rb2e)
         bp(idl) =  bp(idl)+r1e-alm(x21,x22,rb2m,rb2e)
c     
c     save matrices to back out solution
c     ab22i=wb,a33i=tb,a21,a32
c     
         wb11(id) =  w11
         wb12(id) =  w12
         wb21(id) =  w21
         wb22(id) =  w22
         tb11(id) =  t11
         tb12(id) =  t12
         tb21(id) =  t21
         tb22(id) =  t22
         a21mpf(id) =  a21mp
         a21mef(id) =  a21me
         a21epf(id) =  a21ep
         a21eef(id) =  a21ee
         a32mpf(id) =  a32mp
         a32mef(id) =  a32me
         a32epf(id) =  a32ep
         a32eef(id) =  a32ee
c     
c     save vectors neccessary to get solution
c     
         rb2mf(id) =  rb2m
         rb2ef(id) =  rb2e
         r3mf(id) =  r3m
         r3ef(id) =  r3e
      enddo

      return
      end
