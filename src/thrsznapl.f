      subroutine thrsznapl(ndummy)
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
CD1 To calculate coefficients and derivatives for isothermal saturated
CD1 zone water-NAPL system.
CD1 The air phase properties have been replaced by napl properties. The
CD1 viscosity is a constant and the density is a simple linear function 
CD1 of pressure.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/thrsznapl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 ndummy                  int       I   Parameter to set the correct
CD3                                          node number for dual
CD3                                          porosity or dpdp
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
CD4 no, ipdeef, ipdepf, ipdmef, ipdmef, ipdq, ipdqh, dtot, iporos,
CD4 t, crl, neq, ieos, phi, s, pcp, dpcef, qc, ps, ka, esk, wellim,
CD4 pflow, avgmolwt, denh, deni, denei, deneh, dmef, dmpf, depf, deef,
CD4 dil, dilp, dile, dglp, dgle, rolf, div, divp, dive, dgvp, dgve,
CD4 rovf, sk, dq, dqh, qh, drc, deqh, dtpa, dtpae, den
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 None
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5 
CD5 pcl0         real*8      Reference pressure at which reference air
CD5                              density is given
CD5 roc0         real*8      Air density at reference pressure
CD5 
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 dtin         real*8      Reciprocal of time step
CD5 tempc        real*8      Dimensionless absolute temperature
CD5 drocp        real*8      Factor used in density calculation
CD5 rolref       real*8      Liquid density at reference conditions
CD5 xvisl        real*8      Liquid viscosity
CD5 comw         real*8      Factor in liquid density calculation
CD5 pref         real*8      Reference pressure
CD5 xvisv        real*8      Vapor viscosity
CD5 mid          int         Do loop over all nodes
CD5 xrl          real*8      Relative permeability
CD5 mi           int         Index that includes pointer to matrix nodes
CD5 drl          real*8      Derivative of relative permeability
CD5 drlp         real*8      Derivative of relative permeability
CD5 xrv          real*8      Vapor relative permeability
CD5 drv          real*8      Derivative of vapor relative permeability
CD5 drvp         real*8      Derivative of vapor relative permeability
CD5 ieosd        int         Equation of state parameter
CD5 pl           real*8      Pressure
CD5 sl           real*8      Saturation
CD5 pwl          real*8      Pressure difference
CD5 dpwlp        real*8      Derivative of pressure difference with
CD5                              pressure
CD5 dpwls        real*8      Derivative of pressure difference with
CD5                              saturation
CD5 qdis         real*8      Total mass flux source/sink
CD5 qwdis        real*8      Water mass flux source/sink
CD5 qadis        real*8      Air mass flux source/sink
CD5 dqws         real*8      Derivative of water flux with saturation
CD5 dqas         real*8      Derivative of air flux with saturation
CD5 dqwp         real*8      Derivative of water flux with pressure
CD5 dqap         real*8      Derivative of air flux with pressure
CD5 por          real*8      Porosity
CD5 kq           int         Flag denoting the type of source term used
CD5 sflux        real*8      Mass source term at this node
CD5 permsd       real*8      Impedance of source term
CD5 pflowd       real*8      Pressure of source stream
CD5 perml        real*8      Liquid permeability
CD5 dpls         real*8      Derivative of liquid permeability with
CD5                              saturation
CD5 dpas         real*8      Derivative of vapor permeability with
CD5                              saturation
CD5 perma        real*8      Air permeability
CD5 roc          real*8      Gas density
CD5 drocs        real*8      Derivative of gas density with saturation
CD5 rol          real*8      Liquid density
CD5 drolp        real*8      Derivative of liquid density with pressure
CD5 rcomd        real*8      Term used in calculation
CD5 drols        real*8      Derivative of liquid density with
CD5                              saturation
CD5 dena         real*8      Total air density (per unit total volume)
CD5 ddenp        real*8      Derivative of total liquid density with
CD5                                pressure
CD5 ddens        real*8      Derivative of total liquid density with
CD5                                saturation
CD5 ddenap       real*8      Derivative of total air density with
CD5                                pressure
CD5 ddenas       real*8      Derivative of total air density with
CD5                                saturation
CD5 dql          real*8      Liquid kinematic viscosity
CD5 dqv          real*8      Vapor kinematic viscosity
CD5 wimped       real*8      Permeability term used in calculation
CD5 
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
CD9 2.3.7 Sources and sinks
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
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
CPS BEGIN thrair
CPS 
CPS rlperm - compute relative permeabilities
CPS porosi - compute pressure dependent porosity contribution
CPS cappr - compute capillary pressures
CPS
CPS Set constants
CPS
CPS FOR each node
CPS
CPS   Set values used in calculation
CPS   IF the air has a relative permeability
CPS     Set terms
CPS   ELSE the air is fully mobile
CPS     Set terms
CPS   ENDIF
CPS   IF there is a two-phase source with inflow or outflow
CPS     Set source/sink terms
CPS     IF this is a liquid only source
CPS       Compute water source/sink term
CPS     ELSEIF it is a two-phase source
CPS       IF the boundary is held at partial saturation
CPS         Compute air source/sink
CPS       ENDIF
CPS       Compute the water source/sink
CPS     ELSEIF this is an air only source/sink
CPS       Compute air source/sink
CPS     ENDIF
CPS   ELSEIF this is an outflow only node
CPS     Compute air and water sink terms and derivatives
CPS   ELSEIF fluxes are specified directly
CPS     IF total flux and fraction is specified
CPS       Compute air and water source/sink terms and derivatives
CPS     ELSE the air pressure specification option is invoked
CPS       Compute air and water source/sink terms and derivatives
CPS     ENDIF
CPS   ENDIF
CPS   
CPS   Compute average molecular weight, air density, water density
CPS   Compute accumulation terms
CPS   Compute liquid and vapor transport terms
CPS   Set source/sink terms in arrays
CPS   Zero out heat terms
CPS   
CPS ENDFOR each node
CPS 
CPS END thrair
CPS 
C**********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comfi
      use comgi
      use comii
      use comrlp, only : rlpnew
      use comrxni
      use davidi
      implicit none

      integer ndummy,mid,mi,ieosd,kq
c gaz 110819 tref, pref (now global)       
      real*8 dtin,tempc,drocp,drocp0,rolref,xvisl,comw,xvisv
      real*8 xrl,drl,drlp,xrv,drv,drvp,pl,sl,svd,pwl,dpwlp,dpwls
      real*8 qdis,qwdis,qadis,dqws,dqas,dqwp,dqap,por,sflux,permsd
      real*8 pflowd,perml,dpls,dpas,perma,roc,drocs,rol,drolp
      real*8 rcomd,drols,dena,ddenp,ddens,ddenap,ddenas,dql,dqv
      real*8 pcl0,wimped,airmobile, rcomd_napl
	real*8 area, uperm, watterm, airterm, pdifa, pdifw
	real*8 qwmax, wat_ex, dwat_exs, prop_ratio
      parameter(pcl0 = 0.101325)
c gaz debug roc0 now in comai
c      parameter(roc0 = 1.292864)
      parameter(airmobile = 10.0)
	integer iadka 
	parameter(iadka=20) 
c
c NAPL properties(density in kg/m**3, viscosity in Pa*sec) 
c	We now read these in the szna macro
c      parameter(den_Soltrol=789.0) 
c      parameter(den_111TCA=1349.0)
c      parameter(vis_Soltrol=6.12e-3)
c      parameter(vis_111TCA=1.20e-3)
c gaz need to formalize this ratio(Soltrol)
c      
      prop_ratio = (dennapl/viscnapl)/(998./1.e-3)
c
      dtin=1.0/dtot
c     get relative perms
      if (rlpnew) then
         call rlp_cap(ndummy)
      else
         call rlperm(ndummy,1)
      end if
c     calculate variable porosity if enabled
      if(iporos.ne.0) call porosi(1)
c     get capillary pressures
      if (.not. rlpnew) call cappr(1,ndummy)

c     dependent variables vap p and sl
c     misc. constants
      tempc=(273.0)/(tref+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
      xvisl=crl(2,1)
      comw=crl(3,1)
c gaz 110819 pref, tref (global) read in scanin        
c      pref=crl(4,1)
      xvisv=crl(5,1)
      rcomd=comw*rolref
c
c napl is given water compressibility
c
      rcomd_napl=dennapl*comw
      xvisv=viscnapl

c     generate coef and derivatives
      do mid=1,neq
         mi=mid+ndummy
c     water relative perm
         xrl=rlf(mi)
         drl=drlef(mi)
         drlp=drlpf(mi)
c     air relative perm
         xrv=rvf(mi)
         drv=drvef(mi)
         drvp=drvpf(mi)
         ieosd=ieos(mi)
         pl=phi(mi)
         sl=s(mi)

c     add new codeing for sv
c     next line used to carry sv info for irdof=14
c     svd=denj(mi)   
         svd=1.0-s(mi)
         pwl=pl-pcp(mi)
         dpwlp=1.0
         dpwls=-dpcef(mi)
         qdis=qc(mi)
         qwdis=0.0
         qadis=0.0
         dqws=0.0
         dqas=0.0
         dqwp=0.0
         dqap=0.0
         por=ps(mi)
         kq=ka(mi)
c     form flow terms
         if(kq.eq.-1 .and. compute_flow) then
c     flow in or out
            sflux=esk(mi)
            permsd=abs(wellim(mi))
            pflowd=pflow(mi)
            if(sflux.eq.1.0.and.irdof.ne.13) then
               qwdis=permsd*(pl-pflowd) + permsd*(sl-sflux) 
               dqwp=permsd
               dqws=permsd
            else if(sflux.eq.1.0.and.irdof.eq.13) then
               qwdis=permsd*(pl-pflowd) 
               dqwp=permsd
            else if(sflux.gt.0.0.and.sflux.lt.1.0) then
               if(pflowd.ge.0.0) then
                  qadis=permsd*(pl-pflowd)
                  dqap=permsd
                  qwdis=permsd*(sl-sflux)
                  dqws=permsd
               else if(pflowd.lt.0.0) then
                  qwdis=permsd*(sl-sflux)
                  dqws=permsd
               endif
            else if(sflux.le.0.0) then
	         airterm= permsd*prop_ratio
               qadis=airterm*xrv*(pl-pflowd)
               dqap=airterm*xrv + drvp*airterm*(pl-pflowd)
               dqas= drv*airterm*(pl-pflowd)
               qwdis=permsd*xrl*(pl-pflowd)
               dqwp=permsd*xrl + drlp*permsd*(pl-pflowd)
               dqws= drl*permsd*(pl-pflowd)
            endif
            if(wellim(mi).lt.0.and.pl-pflowd.lt.0.0) then
               qwdis = 0.0d00
               dqws = 0.0d00
               dqwp = 0.0d00
            endif
         else if(kq.eq.-2) then
c ===Free drainage condition for water flow at the lower boundary
c ===(each node at the free drainage boundary needs to be specified
c ===explicitly, also the area or length associated with that node).
c           rol = rolref*(1.0+comw*(pl-pref))
            rol = rolref
            area=pflow(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl
            qwdis=uperm*xrl
            dqws=uperm*drl
            dqwp=uperm*drlp
            qadis=0.
            dqas=0.
            dqap=0.
            if(qwdis.lt.0.0) then
             qwdis=-qwdis
             dqws=-dqws
             dqwp=-dqwp
            endif
         else if(kq.eq.-3.and.irdof.ne.13) then
c ===Seepage face condition for air and water 
           if(iad.le.iadka) then
            permsd=abs(wellim(mi))
            watterm= xrl*permsd
            airterm= xrv*permsd*prop_ratio
            pdifa = pl-pflow(mi)
c            pdifw = pwl-pflow(mi)
            pdifw = pl-pflow(mi)
            if(pdifw.gt.0.0) then
             qwdis = watterm*pdifw
             dqws = drl*permsd*pdifw 
             dqwp = drlp*permsd*pdifw + watterm 
            else
             qwdis = 0.0
             dqws = 0.0
             dqwp = 0.0
            endif
            qadis = airterm*(pdifa) 
            dqas = drv*permsd*prop_ratio*pdifa
            dqap = drvp*permsd*prop_ratio*pdifa + airterm
           else
            qwdis = sk(mi)
            qadis = qh(mi) 
           endif
         else if(kq.eq.-3.and.irdof.eq.13) then
c Seepage face for simpliflied water table
c with maximum inflow value
            permsd=abs(wellim(mi))
            watterm= permsd
            pdifw = pwl-pflow(mi)
	      qwdis = watterm*pdifw
            qwmax = esk(mi)
	      if(qwdis.gt.-qwmax) then
             dqws = 0.0
             dqwp = watterm
            else
             qwdis = -qwmax
             dqws = 0.0
             dqwp = 0.0
            endif
         else if(kq.eq.-4.and.irdof.ne.13) then
c ===Ponding condition for air and water 
c ===area input in pflow term
            rol = rolref
            area=pflow(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl
	      sl = s(mi)
            wat_ex= (sl-1.0)/sl
            dwat_exs= 1.0/sl - (sl-1.0)/sl**2
            qwdis=uperm*wat_ex
            dqws=uperm*dwat_exs
            dqwp=0.0              
            qadis=0.
            dqas=0.
            dqap=0.
         else if(kq.gt.0) then
            if(esk(mi).ge.0.0) then
               if(sl.gt.0.0) then
                qwdis=qdis*esk(mi)
               else
                qwdis=0.0           
               endif
c gaz 051809                
c               if(svd.gt.0.0) then
               if(esk(mi).lt.1.0) then               
                qadis=qdis*(1.0-esk(mi))
               else
                qadis=0.0           
               endif
            else
               qwdis=qdis
               airterm = sx1(mi)*xrv
               qadis=airterm*(pl+esk(mi))
               dqap=airterm
               dqas=sx1(mi)*drv*(pl+esk(mi))
            endif
         endif
****  Average molecular weight is equal to the molecular
****  weight of air, since no water is present in the vapor phase
****  in this routine ****
         avgmolwt(mi) = mw_air
c
c incorporate napl changes
c
         if(ico2.eq.-3) then
          roc=dennapl*(1.0+comw*(pl-pref))
          drocp = rcomd_napl
         else
c     density of air
         drocp=drocp0                          
         roc=drocp*pl
         drocs=0.0
	end if
c     water density
         rol=rolref*(1.0+comw*(pl-pref))
         drolp=rcomd
         drols=0.0
c     accumulation terms
         den=por*rol*sl
         dena=por*roc*svd         
         ddenp=por*sl*drolp
         ddens=por*rol
         ddenap=por*drocp*svd         
         ddenas=por*(drocs*svd-roc)             
         deni(mi)=(den-denh(mi))*dtin
         denei(mi)=(dena-deneh(mi))*dtin
         dmef(mi)=ddens*dtin
         dmpf(mi)=ddenp*dtin
         depf(mi)=ddenap*dtin
         deef(mi)=ddenas*dtin
         if(icons.lt.maxit) then

c     make densities constant for transport terms
c     density of air
            roc=drocp0*pref
            drocs=0.0
            drocp=0.0
c     water density
            rol=rolref
            drolp=0.0
            drols=0.0		 
         endif
c     transport terms
c     liquid
         dql=rol/xvisl*xrl
         dil(mi)=dql
         dilp(mi)=drolp/xvisl*xrl
         dile(mi)=drols/xvisl*xrl+rol/xvisl*drl
         dglp(mi)=drolp
         dgle(mi)=drols
         rolf(mi)=rol
c     vapour
         dqv=roc/xvisv*xrv
         div(mi)=dqv
         divp(mi)=drocp/xvisv*xrv
         dive(mi)=drocs/xvisv*xrv+roc/xvisv*drv
         dgvp(mi)=drocp
         dgve(mi)=drocs
         rovf(mi)=roc
c     source terms
         if(ieosd.eq.2) then
            sk(mi)=qwdis
            dq(mi)=dqwp
            dqh(mi)=dqws
            qh(mi)=qadis
            drc(mi)=dqap
            deqh(mi)=dqas
         else if(ieosd.eq.1) then
            sk(mi)=qwdis
            dq(mi)=dqwp
            dqh(mi)=0.0  
            if(strd.ne.1.0.and.qadis.ge.0.0) then
               qh(mi)= sx1(mi)*(sl-1.0)            
               drc(mi)= 0.0 
               deqh(mi)=sx1(mi)
            else
               qh(mi)= qadis            
               drc(mi)= dqap
               deqh(mi)=dqas
            endif
         else if(ieosd.eq.3) then
            sk(mi)= sx1(mi)*(sl-0.0)            
            dq(mi)=sx1(mi)
            dqh(mi)= 0.0 
            qh(mi)=qadis
            drc(mi)=dqap
            deqh(mi)=dqas
         endif
c     zero heat terms
         dtpa(mi)=0.0
         dtpae(mi)=0.0
      enddo

      return
      end
