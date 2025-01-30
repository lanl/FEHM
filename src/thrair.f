       subroutine thrair(ndummy) 
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
CD1 To calculate coefficients and derivatives for isothermal air-water
CD1 system.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                        However, previous non-YMP
CD2                                        versions of FEHM exist, and
CD2                                        the current version differs
CD2                                        from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/thrair.f_a  $
CD2
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.15   Fri Nov 21 16:42:56 1997   gaz
CD2 added sub pcp_save to thrair
CD2 
CD2    Rev 1.14   Fri Sep 26 15:13:52 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.13   Mon Apr 14 12:44:22 1997   gaz
CD2 minor cosmetic change
CD2 
CD2    Rev 1.11   Fri Feb 16 13:01:56 1996   zvd
CD2 Added requirements.
CD2 
CD2    Rev 1.10   Wed Feb 07 11:59:14 1996   gaz
CD2 changes made for bous macro
CD2 
CD2    Rev 1.9   Fri Feb 02 12:50:04 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   11/15/95 15:27:54   gaz
CD2 slightly different formulation works a little better
CD2 
CD2    Rev 1.7   06/02/95 10:27:22   llt
CD2 gaz changes
CD2 
CD2    Rev 1.6   05/01/95 15:18:28   gaz
CD2 some changes to phase derivatives (gaz)
CD2 
CD2    Rev 1.3   03/23/95 19:14:00   gaz
CD2 gaz coding for single phase derivatives
CD2 
CD2    Rev 1.1   03/18/94 15:41:04   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:46   pvcs
CD2 original version in process of being certified
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
CD5 dporpl       real*8      Derivative of porosity with pressure
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
c
c calculate coefficients and derivatives for isothermal air-water
c system
c
      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comevap, only : evaporation_flag, evap_flag
      use com_prop_data, only : ihenryiso, ctest, xnl_ngas, xnl_max,
     &       den_h2o, visc_h2o
      use comfi
      use comgi
      use comii
      use comrlp, only : rlpnew
      use comrxni
      use comwt
      use davidi
      use com_nondarcy, only : nd_flow
      implicit none

      integer ndummy,mid,mi,ieosd,kq
c gaz 110819 removed tref, pref (now global)       
      real*8 dtin,tempc,drocp,drocp0,rolref,xvisl,comw,xvisv
      real*8 xrl,drl,drlp,xrv,drv,drvp,pl,sl,svd,pwl,dpwlp,dpwls
      real*8 qdis,qwdis,qadis,dqws,dqas,dqwp,dqap,por,sflux,permsd
      real*8 pflowd,area,uperm,roc,drocs,rol,drolp,dporpl
      real*8 rcomd,drols,dena,ddenp,ddens,ddenap,ddenas,dql,dqv
      real*8 pcl0,xvisl0,pdiff,plow,watfrac,dwfracp
      real*8 airterm,watterm,pdifa,pdifw   
      real*8 tl_last,viln,vild,vil 
      real*8 viln1,viln2,viln3,vild1,vild2,vild3
      real*8 vla0,vlpa1,vlpa2,vlpa3,vlta1,vlta2,vlta3,vlpta,vlp2ta
      real*8 vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3,vlptb,vlp2tb
      real*8 vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tlx,tl2x,tlx2
      real*8 dsatp, rlpmin, qwmax_fac, qwmax
      real*8 pld,dis_ex,wat_ex,dwat_exs, time_max
      real*8 seep_facv,seep_facl,permsdv,permsdl,plwt
      real*8 rol_b, rolref_b, tol_dis, area_term, dis_term, vap_ratio
      real*8 cden_correction, cden_cor

c gaz 081323
      real*8 xrvdum,drvdum,xrldum,drldum,air_fac
c gaz 111123
      real*8 xnl, dxnlp, dxnls
      integer i_mem_rlp
      integer iadka        
      save i_mem_rlp
      parameter(pcl0 = 0.101325)
c      parameter(roc0 = 1.292864)
      parameter(rlpmin = 0.01, tol_dis = 1.e-12)
      parameter(qwmax_fac = 1., air_fac = 0.05d0)

      parameter(iadka=1000)                
      
      real*8, allocatable :: s0(:)
      real*8, allocatable :: pcp0(:)
      real*8, allocatable :: dpcps0(:)
      real*8, allocatable :: rlf0(:)
      real*8, allocatable :: drlfs0(:)
      real*8, allocatable :: rvf0(:)
      real*8, allocatable :: drvfs0(:)
c gaz 110224
c      if(nd_flow) then
c         if(.not.allocated(den_h2o)) allocate(den_h2o(n,6))
c         if(.not.allocated(visc_h2o)) allocate(visc_h2o(n,6))
c      endif
      if(irdof.ne.13) then
         if(abs(iexrlp).ne.0.and.i_mem_rlp.eq.0) then
            i_mem_rlp=1
            allocate(s0(n),pcp0(n),dpcps0(n),rlf0(n))
            allocate(drlfs0(n),rvf0(n),drvfs0(n))       
         endif
         
c     get relative perms
         if(iad.lt.abs(iexrlp).or.abs(iexrlp).eq.0) then
            if (rlpnew) then
               call rlp_cap(ndummy)
            else
               call rlperm(ndummy,1)
            end if
         else if(iad.eq.abs(iexrlp)) then
            if (rlpnew) then
               call rlp_cap(ndummy)
            else
               call rlperm(ndummy,1)
            end if
            call pcp_save(0,neq,ndummy,0,pcp,dpcef,rlf,drlef,
     &           rvf,drvef,s,pcp0,dpcps0,rlf0,drlfs0,
     &           rvf0,drvfs0,s0)
         else if(iad.gt.abs(iexrlp)) then
            call pcp_save(2,neq,ndummy,0,pcp,dpcef,rlf,drlef,
     &           rvf,drvef,s,pcp0,dpcps0,rlf0,drlfs0,
     &           rvf0,drvfs0,s0)
            do mid=1,neq
               mi=mid+ndummy
               drlpf(mi)=0.0          
               drvpf(mi)=0.0
            enddo
         endif
c     get capillary pressures
         if (.not. rlpnew) call cappr(1,ndummy)

      endif

c     
c     calculate variable porosity if enabled
c     
      if(iporos.ne.0) call porosi(1)
c     
c     dependent variables vap p and sl
c     
c     misc. constants (roc is at tref)
c gaz 110819 pref, tref (global) read in scanin        
c      tref = crl(6,1)
      tempc=(273.0)/(tref+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
      xvisl0=crl(2,1)
      comw=crl(3,1)
c      pref=crl(4,1)
      xvisv=visc_gas
c      rcomd=comw*rolref

      seep_facv = roc0/xvisv
c     
c     liquid viscosity
c     numerator coefficients
      vla0=cvl(1,1)
      vlpa1=cvl(2,1)
      vlpa2=cvl(3,1)
      vlpa3=cvl(4,1)
      vlta1=cvl(5,1)
      vlta2=cvl(6,1)
      vlta3=cvl(7,1)
      vlpta=cvl(8,1)
      vlp2ta=cvl(9,1)
      vlpt2a=cvl(10,1)
c     denomenator coefficients
      vlb0=cvl(11,1)
      vlpb1=cvl(12,1)
      vlpb2=cvl(13,1)
      vlpb3=cvl(14,1)
      vltb1=cvl(15,1)
      vltb2=cvl(16,1)
      vltb3=cvl(17,1)
      vlptb=cvl(18,1)
      vlp2tb=cvl(19,1)
      vlpt2b=cvl(20,1)
      x = pref
      x2 = x*x
      x3 = x2*x
      tl_last = tref
c     
      dtin=1.0/dtot
c     
c     generate coef and derivatives
c gaz 071819  
c     save some source terms if iad = 0
c
c      if(iad.eq.0) then
c       call source_manage(0)
c      endif
      ifree1 = 0 
      do 100 mid=1,neq
         mi=mid+ndummy
         ieosd=2 
c     count non fully saturated cells	    
         pl =phi(mi) 
         if (ifree .ne. 0) then  
            sl = rlxyf(mi) - rlptol
            if(sl.lt.1.0 - rlptol) ifree1 = ifree1 + 1
         else if (irdof .ne. 13 ) then
            sl = s(mi) 
            if(s(mi).lt.1.0) ifree1 = ifree1 + 1
         else
            sl = 1.0d0
         end if
c     
c     liquid viscosity
c     
         if(t(mi).ne.tl_last) then
            tl = t(mi)
            tl2 = tl*tl
            tl3 = tl2*tl
            tlx = tl*x
            tl2x = tl2*x
            tlx2 = tl*x2
            viln1=vla0+vlpa1*x+vlpa2*x2+vlpa3*x3
            viln2=vlta1*tl+vlta2*tl2+vlta3*tl3
            viln3=vlpta*tlx+vlpt2a*tl2x+vlp2ta*tlx2
            viln=viln1+viln2+viln3
            vild1=vlb0+vlpb1*x+vlpb2*x2+vlpb3*x3
            vild2=vltb1*tl+vltb2*tl2+vltb3*tl3
            vild3=vlptb*tlx+vlpt2b*tl2x+vlp2tb*tlx2
            vild=vild1+vild2+vild3
            vil=viln/vild
            xvisl=vil
         else
            xvisl=xvisl0
         endif
c manage liquid phase density and viscosity    
       if(cden)then
        rolref_b= rolref+cden_correction(mi)
       else if(iden_vis.gt.0) then
c spatially variable density and viscosity
        rolref_b = den_spatial(mi)
        xvisl = vis_spatial(mi)
        xvisl0 = xvisl
        if(comp_spatial(mi).gt.0.0) then
          comw = comp_spatial(mi)
        endif
       else
        rolref_b= rolref
       endif
       rcomd=comw*rolref_b
       rolf(mi) = rolref_b
c manage gas phase density and viscosity 
       if(ideng_vis.gt.0) then
c spatially variable density and viscosity
        roc0 = deng_spatial(mi)
        drocp0=roc0*tempc/pcl0
        xvisv = visg_spatial(mi)
       endif
c
c
       seep_facl = rolref_b/xvisl0

c     
         if(irdof.ne.13) then
c     water relative perm
            xrl=rlf(mi)
            drl=drlef(mi)

            drlp=drlpf(mi)
            dsatp=0.0
c     air relative perm
            xrv=rvf(mi) 
            drv=drvef(mi)
            drvp=drvpf(mi)
c     cappilary pressure and derivatives
            pwl=pl-pcp(mi)
            dpwlp=1.0
            dpwls=-dpcef(mi)
         else
            xrl=1.0d00
            drl=0.0d00
            drlp=0.0d00
            pwl=pl
            dpwlp=1.0d00
            dpwls=0.0d00  
            xrv = 0.0   
            drvp = 0.0
            drv = 0.0 
            if(ifree.ne.0) then
               dsatp=drlxyf(mi)
            else
               dsatp = 0.0
            endif
         endif
c     
c     
c     add new coding for sv
c     
c     next line used to carry sv info for irdof=14
c     
         if(abs(irdof).eq.14) then
            svd=denj(mi)   
         else
            svd = 1.0 - sl
         endif
c     
         qdis=qc(mi)
         qwdis=0.0
         qadis=0.0
         dqws=0.0
         dqas=0.0
         dqwp=0.0
         dqap=0.0
         por=ps(mi)
         kq=ka(mi)
c gaz 111223 get liq mass frac dissolved gas
c gaz 010624   these calls duplicates not necessary
         if(l.lt.-1) then
c        if(l.gt.0) then
          call phase_change_mass_conv(1,mi,mi)
          call solubility_isothermal(1,mi)
         endif
c gaz 040624 invoke  only for  ihenryiso ne 0
          call solubility_isothermal(1,mi)
          if(sl.lt.1.0d0) then
           xnl    = xnl_ngas(mi,1)
           dxnlp  = xnl_ngas(mi,2)
           dxnls  = xnl_ngas(mi,3)
           cnlf(mi) = xnl
           dclf(mi) = dxnlp
           dclef(mi) = dxnls
          else
           xnl_max = xnl_ngas(mi,1)
           xnl = cnlf(mi)
           dxnlp =0.0d0
           dclf(mi) = 0.d0
c gaz 122623
           dclef(mi) = 1.0d0
          endif
c gaz 071224  check for head bc/elevation inconsistency
          if(ihead.ne.0) then
           call hdbc_mod(2,mi,mi,0)
           kq = ka(mi)
          endif
c 
c     form flow terms
         if(kq.eq.-1.or.kq.eq.-2.and.compute_flow) then
c     flow in or out
            sflux=esk(mi)
            permsd=abs(wellim(mi))
            if(iwelimd.ne.0)then
c 
c peaceman solution only available for models kq(-1,-2,1)
c           
             if(izonewel1(mi).ne.0) then
              permsd = wellim(mi)*rolref_b/xvisl
             endif
            endif
            pflowd=pflow(mi)
            if(sflux.eq.1.0.and.irdof.ne.13) then
               qwdis=permsd*(pl-pflowd) + permsd*(sl-sflux) 
               dqwp=permsd
               dqws=permsd
            else if(sflux.eq.1.0.and.irdof.eq.13) then
               qwdis=permsd*(pl-pflowd) 
               dqwp=permsd
            else if(sflux.eq.-1.000) then
                qadis=permsd*(pl-pflowd) 
                dqap=permsd
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
               permsdv = permsd*seep_facv
               qadis=permsdv*xrv*(pl-pflowd)
               dqap=permsdv*xrv + drvp*permsdv*(pl-pflowd)
               dqas= drv*permsdv*(pl-pflowd)
               permsdl = permsd*seep_facl
               qwdis=permsdl*xrl*(pl-pflowd)
               dqwp=permsdl*xrl + drlp*permsdl*(pl-pflowd)
               dqws= drl*permsdl*(pl-pflowd) 
            endif
            
            if(kq.eq.-2.and.sflux.lt.0.0) then 
               qadis = permsd*(pl-pflowd)
               dqap  = permsd 
               dqas = 0.0
               qwdis = permsd*(sl)
               dqws = permsd
               dqwp = 0.0d00     
            else if(kq.eq.-2.and.qwdis.lt.0.0) then
               qwdis = 0.0d00
               dqws = 0.0d00
               dqwp = 0.0d00
            endif
         else if(kq.eq.-23) then
c     ===Free drainage condition for water flow at the lower boundary
c     ===(each node at the free drainage boundary needs to be specified
c     ===explicitly, also the area or length associated with that node)
c     area is read in the flow rate or pressure slot and transferred to array wellim
c     rol = rolref*(1.0+comw*(pl-pref))
            rol = rolref_b
            area=wellim(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl
            qwdis=uperm*xrl
            dqws=uperm*drl
            dqwp=uperm*drlp
            qadis=0.
            dqas=0.
            dqap=0.
            if(qwdis.lt.0.0) then
c gaz debug 040316
               qwdis=0.0
               dqws=0.0
               dqwp=0.0
            endif
      else if(kq.eq.-24) then
c     same as -23 except  xrl = s and drl = 1    
c     ===Free drainage condition for water flow at the lower boundary
c     ===(each node at the free drainage boundary needs to be specified
c     ===explicitly, also the area or length associated with that node)
c     area is read in the flow rate or pressure slot and transferred to array wellim
c     rol = rolref*(1.0+comw*(pl-pref))
            rol = rolref_b
            area=wellim(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl
            qwdis=uperm*sl
            dqws=uperm
            dqwp=0.
            qadis=0.
            dqas=0.
            dqap=0.
            if(qwdis.lt.0.0) then
               qwdis=-qwdis
               dqws=-dqws
               dqwp=-dqwp
            endif
         else if(kq.eq.-3.and.ifree.eq.0) then
c     ===Seepage face condition for air and water 

            permsd=abs(wellim(mi))
            permsdv = permsd*seep_facv
            permsdl = permsd*seep_facl
c gaz debug 081119            
c            watterm= xrl*permsdl
c            airterm= xrv*permsdv
            watterm= sl*permsdl
            airterm= (1.0-sl)*permsdv   
c gaz 110819 pref (crl(4,1) now pref and inputed elsewhere            
            pflow(mi) = pref
            pdifa = pl-pflow(mi)
c     
            pdifw = pl-pflow(mi)
            if(pdifw.gt.0.0) then
               qwdis = watterm*pdifw
c  gaz debug  081119      
c               dqws = drl*permsdl*pdifw 
c               dqwp = drlp*permsdl*pdifw + watterm 
c  
               dqws = permsdl*pdifw 
               dqwp = watterm                
            else
               qwdis = 0.0
               dqws = 0.0
               dqwp = 0.0
            endif
            if(jswitch.eq.0) then
c gaz debug 121219                 
c             qadis = airterm*(pdifa) 
c             dqas = -1.0*permsdv*pdifa
             qadis = permsdv*(pdifa) 
             dqas = 0.
             dqap =  permsdv
            else
             qadis = 0.0d0
             dqas = 0.0d0
             dqap = 0.0d0                
            endif


         else if(kq.eq.-3.and.ifree.ne.0) then
c     Seepage face for simpliflied water table
c     outflow only based on highest cell pressure (sat = 1.)
c     
c     permsd= qc(mi)/wellim(mi)
            permsd=abs(wellim(mi))
            permsdl = permsd*seep_facl
            plow = head12(mi,2)
            watterm= permsdl*esk(mi)
            watfrac = rlxyf(mi)
            dwfracp = drlxyf(mi)
            pdifw = pwl-plow
            if(pdifw.ge.0.0d0) then
               qwdis = watterm*watfrac*pdifw
               dqws = 0.0
               dqwp = watterm*watfrac + watterm*dwfracp*pdifw
            else
	       qwdis = 0.0
	       dqws = 0.0
             dqwp = 0.0
            endif


         else if(kq.eq.-4.and.ifree.ne.0) then
c     Seepage face for simpliflied water table
c     
c     based on average cell pressure
c     permsd= qc(mi)/wellim(mi)
c            permsd=abs(wellim(mi))
            if(igrav.eq.2) then
             permsd = pnx(mi)*sx1(mi)**0.3333
            elseif(igrav.eq.3) then
             permsd = max(pnx(mi),pny(mi))*sx1(mi)**0.3333
            endif
            permsdl = permsd*seep_facl*esk(mi) 
            plow = 0.5*(head12(mi,2)+head12(mi,1))
            watterm= permsdl
            watfrac = rlxyf(mi) 
            dwfracp = drlxyf(mi)
            pdifw = pwl-plow
            if(pdifw.ge.0.0d0) then
               qwdis = watterm*watfrac*pdifw
               dqws = 0.0
               dqwp = watterm*watfrac + watterm*dwfracp*pdifw
            else
	       qwdis = 0.0
	       dqws = 0.0
             dqwp = 0.0
            endif
            
        else if(kq.eq.-4.and.ifree.eq.0) then
c gaz 081219
c     Seepage face uz and rich eq
c     fully saturated only
c     based on average cell pressure
c     permsd= qc(mi)/wellim(mi)
c            permsd=abs(wellim(mi))
            if(igrav.eq.2) then
             permsd = pnx(mi)*sx1(mi)**0.3333
            elseif(igrav.eq.3) then
             permsd = max(pnx(mi),pny(mi))*sx1(mi)**0.3333
            endif
c gaz 121219 see below permsdl = permsd*seep_facl*esk(mi) 
c gaz 081219             
c            plow = crl(4,1)
            plow = pflow(mi)            
            watterm = permsd*abs(wellim(mi))
            pdifw = pwl-plow
c gaz 121219 use  esk for sl            
            if(pdifw.ge.0.0d0.and.sl.ge.esk(mi)) then
               qwdis = watterm*pdifw
               dqws = 0.0
               dqwp = watterm
            else
	       qwdis = 0.0
	       dqws = 0.0
             dqwp = 0.0
            endif
            else if(kq.eq.-5.and.ifree.ne.0) then
c     Seepage face for simpliflied water table
c     when cell is just wet
c     
            permsd=abs(wellim(mi))
            permsdl = permsd*seep_facl
            watterm= permsdl*esk(mi)
            watfrac = max(rlxyf(mi)-rlptol,0.0d00)
            plow = head12(mi,1)    
c             plow = 0.1      
            qwdis = watterm*watfrac*(pwl-plow)
c
            dqwp = watterm*drlxyf(mi)*(pwl-plow)+watterm*watfrac
             if(qwdis.lt.0.0d0) then
               qwdis = 0.0
               dqws = 0.0
               dqwp = 0.0
            endif
         else if(kq.eq.-6.and.irdof.ne.13) then
c     ===Ponding condition for air and water 
c     ===area input in pflow term
            rol = rolref_b
            area=pflow(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl            
            wat_ex= (sl-1.0)/sl
            dwat_exs= 1.0/sl - (sl-1.0)/sl**2
            qwdis=uperm*wat_ex
            dqws=uperm*dwat_exs
            dqwp=0.0d00              
            qadis=0.0d00
            dqas=0.0d00
            dqap=0.0d00
         else if(kq.eq.-7.and.ifree.eq.0) then
c     ===specified air pressure no flow water
c     flow in or out
c gaz 081323 allow for outflow only
c gaz 082223 modied
           pflowd=pflow(mi)
           permsd=abs(wellim(mi))
           pdiff = pl - pflowd
           if(pdiff.le.0.0d0) then 
             permsd = 0.0d0                                        
             pdiff = 0.0d0                     
            endif
            xrvdum = 1.0d0-sl
            drvdum = -1.0d0
            qadis = xrvdum*air_fac*permsd*(pdiff)
            dqap = permsd*air_fac*xrvdum
            dqas = permsd*air_fac*(pdiff)*drvdum
            qwdis = 0.0d00
            dqws = 0.0d00
            dqwp = 0.0d00  
         else if(kq.eq.-18.and.ifree.eq.0) then
c gaz 082023 ===specified  pressure  flow water  and air
c  linear rel perm
c gaz 081323 allow for outflow only
           pflowd=pflow(mi)
           permsd=abs(wellim(mi))
           pdiff = pl - pflowd
           if(pdiff.le.0.0d0) then 
             permsd = 0.0d0                                        
             pdiff = 0.0d0                     
            endif
            xrldum = sl
            drldum = 1.0d0
            xrvdum = 1.0d0-sl
            drvdum = -1.0d0
            qadis = xrvdum*air_fac*permsd*(pdiff)
            dqap = permsd*air_fac*xrvdum
            dqas = permsd*air_fac*(pdiff)*drvdum
            qwdis = xrldum*permsd*(pdiff)
            dqwp = permsd*xrldum
            dqws = permsd*(pdiff)*drldum
c            if(days.ge.40.8d00) then
c             write(ierr,'(1x,1p,g16.8,1x,i6,7(1x,g16.8))') 
c     &       days, iad, sl, pdiff, xrldum, xrvdum,
c     &       permsd, qadis, qwdis
c            endif
         else if(kq.eq.-9.and.ifree.eq.0) then
c     ===specified air pressure 
c     ===specified saturation
c     no derivatives of rlps
c     flow in or out
            permsd=abs(wellim(mi))
            pflowd=pflow(mi)
            sflux = esk(mi)
            qadis =permsd*(pl-pflowd)
            dqap = permsd
            dqas = 0.0d0
            qwdis = permsd*(sl-sflux)
            dqws = permsd
            dqwp = 0.0d00
         else if(kq.eq.-10.and.ifree.eq.0) then
c     ===specified air pressure 
c     ===specified flux
c     no derivatives of rlps
            if(sl.lt.1.d0) then
               permsd=abs(wellim(mi))
               pflowd=pflow(mi)
               sflux = esk(mi)
               qadis =permsd*(pl-pflowd)
               dqap = permsd
               dqas = 0.0d0
               qwdis = sflux
               dqws = 0.0d00
               dqwp = 0.0d00
            else
               qadis = 0.0d0
               dqap = 0.0d0
               dqas = 0.0d0
               qwdis = sflux
               dqws = 0.0d00
               dqwp = 0.0d00
            endif
         else if(kq.eq.-13.and.ifree.eq.0) then
c     ===specified saturation (air source)
               permsd=abs(wellim(mi))
               pflowd=pflow(mi)
               sflux = esk(mi)             
               qadis =-permsd*(sl-sflux)
               dqap = 0.0
               dqas = -permsd
               qwdis = 0.0
               dqws = 0.0d00
               dqwp = 0.0d00  
         else if(kq.eq.-14.and.ifree.eq.0) then
c     ===specified saturation (air source), outflow only
               permsd=abs(wellim(mi))
               if(sl-sflux.ge.0.0d0) permsd = 0.0d0
               sflux = esk(mi)             
               qadis =-permsd*(sl-sflux)
               dqap = 0.0
               dqas = -permsd
               qwdis = 0.0
               dqws = 0.0d00
               dqwp = 0.0d00    
         else if(kq.eq.-15.and.ifree.eq.0) then
c     ===specified saturation (air source), outflow only
c     ===specified pressure (water source), outflow and inflow
c     ++++ assumes no capillary pressure
               permsd=abs(wellim(mi))
               sflux = esk(mi) 
               if(sl-sflux.le.0.0d0) then            
                qadis =-permsd*(sl-sflux)
                dqap = 0.0
                dqas = -permsd
               else
                qadis = 0.0
                dqap = 0.0
                dqas = 0.0
               endif
               pflowd = pflow(mi) 
               qwdis = permsd*(pl-pflowd)
               dqws = 0.0
               dqwp = permsd 
         else if(kq.eq.-16.and.ifree.eq.0) then
c     ===specified air source
c     ===specified water source
c     ++++ assumes no capillary pressure
          
                qadis =esk(mi)
                dqap = 0.0
                dqas = 0.0
                qwdis = sk(mi)
                dqws = 0.0
                dqwp = 0.0 
         else if(kq.eq.-17.and.ifree.eq.0) then
c     ===specified saturation (air source), outflow only
c     ===specified pressure (water source), outflow only
c     ++++ assumes no capillary pressure
               permsd=abs(wellim(mi))
               sflux = esk(mi) 
               if(sl-sflux.le.0.0d0) then            
                qadis =-permsd*(sl-sflux)
                dqap = 0.0
                dqas = -permsd
               else
                qadis = 0.0
                dqap = 0.0
                dqas = 0.0
               endif
               pflowd = pflow(mi) 
               if(pl-pflowd.ge.0.0) then
                qwdis = permsd*(pl-pflowd)
                dqws = 0.0
                dqwp = permsd  
               else
                qwdis = 0.0
                dqws = 0.0
                dqwp = 0.0 
               endif                           
         else if(kq.eq.-11.or.kq.eq.-12) then
c     ===x dir generalized head BC
c     ===wellim calculated in area_flow_bc
               permsd=abs(wellim(mi))
               pflowd=pflow(mi) 
               qwdis = permsd*(phi(mi)-pflowd)
               dqwp = permsd
               if(wellim(mi).le.0.0.and.qwdis.le.0.0) then
                 qwdis = 0.0
                 dqwp = 0.0
               endif
               qadis = 0.0d0
               dqap = 0.0d0
               dqas = 0.0d0                           
c     ===specified flux
c     no derivatives of rlps
         else if(kq.eq.-22.and.irdof.ne.13) then
c     ===manage outflow only conditions
            if(sflux.lt.0.0) then
               permsd=abs(wellim(mi))
               qadis = permsd*(pl-pflowd)
               dqap  = permsd 
               qwdis = permsd*(sl)
               dqws = permsd
               dqwp = 0.0d00     
            else if(pl-pflowd.lt.0.0) then
               qwdis = 0.0d00
               dqws = 0.0d00
               dqwp = 0.0d00
            endif
         else if(kq.eq.-19.and.irdof.ne.13) then  
c gaz 121519 conditional air source with water spec flowrate
            permsd=abs(wellim(mi))
            permsdv = permsd*seep_facv

             qadis = 0.0d0
             dqas = 0.0d0
             dqap = 0.0d0    
c pflow is flowrate             
               qwdis = pflow(mi)
               dqws  = 0.0d00
               dqwp  = 0.0d00 
          if(sl.lt.esk(mi)) then
               qadis = permsdv*(pl-pref)
               dqap  = permsdv 
          elseif(sl.lt.1.0) then
               qadis = permsdv*(pl-pref)
               dqap  = permsdv            
          endif                     
         else if(kq.eq.-20.and.irdof.ne.13) then  
c gaz 121519 conditional air source for boundaries
     
            if(igrav.le.2) then
             permsd = pnx(mi)*sx1(mi)**0.3333
            elseif(igrav.eq.3) then
             permsd = max(pnx(mi),pny(mi))*sx1(mi)**0.3333
            endif
            permsdv = permsd*abs(wellim(mi))
            permsdl = permsd*abs(wellim(mi))*esk(mi)
             qadis = 0.0d0
             dqas = 0.0d0
             dqap = 0.0d0        
               qwdis = 0.0d00
               dqws  = 0.0d00
               dqwp  = 0.0d00 
          if(sl.le.0.0) then
               if(jswitch.eq.0) then
                qadis = permsdv*(pl-pflow(mi))
                dqap  = permsdv 
              endif
          elseif(sl.lt.1.0) then
               if(jswitch.eq.0) then
                qadis = permsdv*(pl-pflow(mi))
                dqap  = permsdv 
              endif
               qwdis = permsdl*(sl)
               dqws  = permsdl
               dqwp = 0.0d00  
          else if(sl.ge.1.0) then
              if(pl-pflow(mi).ge.0.0) then
               qwdis = permsdl*(1.0 + (pl-pflow(mi)))
               dqwp  = permsdl  
              else
               qwdis = permsdl
               dqwp = 0.0
              endif 
          endif 
         else if(kq.eq.-25.and.irdof.ne.13) then  
c gaz 121519 conditional air source for boundaries
c gaz 122219
            
            if(igrav.le.2) then
             dis_term = dxrg(mi) + tol_dis
             area_term = sx1(mi)/(dis_term*dis_term)
             permsd = pnx(mi)*dil(mi)*area_term
            elseif(igrav.eq.3) then
             dis_term = max(dxrg(mi),dyrg(mi)) + tol_dis 
             area_term = sx1(mi)/(dis_term*dis_term)
             permsd = max(pnx(mi),pny(mi))*dil(mi)*area_term
            endif
            vap_ratio = max(0.,div(mi)/(dil(mi)+tol_dis))
     &        *abs(wellim(mi)*1.e-6)
            permsdv = permsd*vap_ratio/(esk(mi)+tol_dis)
            permsdl = permsd/(esk(mi)+tol_dis)
             qadis = 0.0d0
             dqas = 0.0d0
             dqap = 0.0d0        
               qwdis = 0.0d00
               dqws  = 0.0d00
               dqwp  = 0.0d00 
          if(sl.le.0.0) then
               if(jswitch.eq.0) then
                qadis = permsdv*(pl-pflow(mi))
                dqap  = permsdv 
              endif
          elseif(ieos(mi).eq.2) then
               if(jswitch.eq.0) then
                qadis = permsdv*(pl-pflow(mi))
                dqap  = permsdv 
              endif
               qwdis = permsdl*(sl)
               dqws  = permsdl
               dqwp = 0.0d00  
          else if(sl.gt.1.0) then
              if(pl-pflow(mi).ge.0.0) then
               qwdis = permsdl*(1.0 + (pl-pflow(mi)))
               dqwp  = permsdl  
              else
               qwdis = permsdl
               dqwp = 0.0
              endif 
          endif                   
         else if(kq.eq.-21.and.irdof.ne.13) then  
c gaz 111519 conditional air source for boundaries
c no water source allowed
c gaz 121219 added permsd = here
          if(sl.lt.esk(mi)) then
               permsd=abs(wellim(mi))
               pflowd = pflow(mi) 
               qadis = permsd*(pl-pflowd)
               dqap  = permsd 
          else
               qadis = 0.0
               dqap  = 0.0             
          endif
               qwdis = 0.0d00
               dqws = 0.0d00
               dqwp = 0.0d00
         else if(kq.eq.-22.and.irdof.eq.13) then
c     make sure there is a deivative wrt P for ifree ne 0
         else if(kq.eq.-101) then
c     change specified pressure to head condition
              permsd=abs(wellim(mi))
              pflowd = pflow(mi) 
              qwdis = permsd*(pl-pflowd)
              if(abs(qwdis).gt.esk(mi)) then
               pflow(mi) = qwdis/permsd + pl
               qwdis = 0.0
               dqwp = permsd
              else
               qwdis = permsd*(pl-pflowd)
               dqwp = permsd
              endif
         else if(kq.eq.-202) then
c gaz 010724 specified pressure and mass fraction
c pflow(mi) specified pressure(liq)
c esk(mi)   specified mass fraction
              permsd=abs(wellim(mi))
              pflowd = pflow(mi)
              sflux = esk(mi)
c gaz 041524 added sat = 1 constraint
              qwdis = permsd*(pl-pflowd) 
c removed 
c     &       + permsd*(sl-1.0)
              dqwp  = permsd 
              dqws = 0.0d0
              qadis = permsd*(cnlf(mi)-sflux)
              dqas = permsd*dclef(mi) 
              dqap = permsd*dclf(mi)
         else if(kq.eq.2) then
            qwdis = qc(mi)
            qadis = esk(mi)
            dqws  = 0.0
            dqwp = 0.0
            dqap  = 0.0
            dqas  = 0.0
         else if(kq.eq.-8) then
            pflowd = pflow(mi)
            sflux = esk(mi)
            permsd = wellim(mi) 
            if(sflux.ge.1.0) then
               qwdis=permsd*(pl-pflowd)
               dqwp=permsd
               dqws=0.0

            else if(sflux.gt.0.0.and.sflux.lt.1.0) then

               qadis=permsd*(pl-pflowd)
               dqap=permsd
               qwdis=permsd*(sl-sflux)
               dqws=permsd

            else if(sflux.lt.0.0) then
               qadis=permsd*(pl-pflowd)
               dqap=permsd
               qwdis = 0.0
               dqws = 0.0
            endif
         else if(kq.eq.1) then
            if(ifree.ne.0) then
c     outflow only if water present
               if(qdis.gt.0.0.and.esk(mi).lt.0.0) then
c     permsd = sx1(mi)
                  qwdis = qdis*rlxyf(mi)
                  dqwp  = qdis*drlxyf(mi)
c     qwdis = qdis
c     dqwp = 0.0
	       else
                  qwdis = qdis
                  dqwp = 0.0
	       endif
            else if(esk(mi).ge.0.0) then
               if(sl.gt.0.0) then
                  qwdis=qdis*esk(mi)
               else
                  qwdis=qdis*esk(mi)          
               endif
               if(svd.gt.0.0) then
                  qadis=qdis*(1.0-esk(mi))
               else
                  qadis=0.0   
c     GAZ 11-10-08                  
                  qadis=qdis*(1.0-esk(mi))        
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
c     density of air
c     GAZ had to move it here for constant density approx
c     added the derivative of porosity wrt pressure
         dporpl=dporp(mi)
c     drocp=drocp0                     
c     roc=drocp0*pl
         drocp=drocp0  
         
         roc=drocp0*(pl)
         drocs=0.0

c     water density
         if(ihead.ne.0) then
            pld = phi(mi) - phi_inc
         else
            pld=phi(mi)
         endif
c gaz 042720 error fix    
         rol=rolref*(1.0+comw*(pld-pref))
         if(cden) rol= rol+cden_correction(mi)
c         rol = rolf(mi)
         drolp=rcomd
         drols=0.0
c     accumulation terms
c gaz 111123 add Henry's law option
c        den=por*(sl*rol*(1.0-xnl)+sv*rov*(1.0-xnv))
c        denc=por*(sl*rol*xnl+sv*rov*xnv)
c pressure derivatives  (AWH)
c        rop=por*(sv*drovp*(1.0-xnv)+sl*drolp*(1.0-xnl))    
c        damp =rop*dtin + por*(sv*rov*(-dxnvp)
c     &  +sl*rol*(-dxnlp))*dtin    
c        dacp =(por*(sv*drovp*xnv+sv*rov*dxnvp+sl*drolp*xnl+sl*rol
c     &   *dxnlp))*dtin
        
         if(ihenryiso.eq.0) then
c no henry's law
          den=por*rol*sl
          dena=por*roc*svd
          ddenp=por*sl*drolp+dporpl*rol*sl 
          ddens=por*rol
          ddenap=por*drocp*svd + dporpl*roc*svd
          ddenas=por*(drocs*svd-roc) 
c gaz 040624
          cnlf(mi) = 0.0d0
          dclf(mi) = 0.0d0
          dclef(mi) = 0.0d0
         else      
c henry's law
         if(sl.lt.1.0d0.and.sl.gt.0.0d0) then
c 2 phase, Variables P and S, xnl= F(P,S)
c gaz 122623
c          den=por*rol*sl*(1.0-xnl)
          den=por*rol*sl
          dena=por*(roc*svd+rol*sl*xnl)
c          ddenp=por*(sl*drolp*(1.0-xnl)-rol*sl*dxnlp)+
c     &    dporpl*(rol*sl*(1-xnl))
          ddenp=por*(sl*drolp) + dporpl*(rol*sl)
c          ddens=por*(rol*(1.0-xnl) - rol*sl*dxnls)
          ddens=por*rol
          ddenap=por*(drocp*svd+drolp*sl*xnl+rol*sl*dxnlp)+ 
     &    dporpl*(roc*svd+rol*sl*xnl)
          ddenas=por*(drocs*svd-roc+rol*xnl+rol*sl*dxnls)  
c          cnlf(mi) = xnl
c          dclf(mi) = dxnlp
c          dclef(mi) = dxnls
         else
c one phase, variables P and xnl 
c          den=por*rol*sl*(1.0-xnl)
          den=por*rol*sl
          dena=por*(roc*svd+rol*sl*xnl)
          ddenp=por*sl*drolp+dporpl*rol*sl 

c ddens :derivative den wrt xnl
c          ddens=por*(-rol*sl)
          ddens=0.0d0
          ddenap=por*(drocp*svd+drolp*sl*xnl)+ 
     &    dporpl*(roc*svd+rol*sl*xnl)
          ddenas=por*rol*sl  

         endif
         endif  
c 110124 check for non darcy flow  (liquid part)
          if(nd_flow) then 
           den_h2o(mi,1)= rol 
           den_h2o(mi,2) = drolp 
           den_h2o(mi,4)= roc 
           den_h2o(mi,5) = drocp  
           visc_h2o(mi,1) = xvisl 
           visc_h2o(mi,2) = 0.0d0 
           visc_h2o(mi,4) = xvisv 
           visc_h2o(mi,5) = 0.0d0     
          endif  
         if(irdof.eq.13) then
            ddenp=ddenp + por*rol*dsatp
         endif
         deni(mi)=(den-denh(mi))*dtin
         dmpf(mi)=ddenp*dtin
         if(icons.lt.abs(maxit)) then 
c     
c     make densities constant for transport terms
c     
c     density of air
            roc=drocp0*pref  
            drocs=0.0
            drocp=0.0
c     water density
            rol=rolref_b
            drolp=0.0
            drols=0.0		 
         endif
c     transport terms
c     liquid
         dql=rol/xvisl*xrl
         dil(mi)=dql
         dilp(mi)=drolp/xvisl*xrl
         dglp(mi)=drolp
c         rolf(mi)=rol
c     vapour
         if(irdof.ne.13) then
            dqv=roc/xvisv*xrv
            div(mi)=dqv
            divp(mi)=drocp/xvisv*xrv
            dive(mi)=drocs/xvisv*xrv+roc/xvisv*drv
            dile(mi)=drols/xvisl*xrl+rol/xvisl*drl
            dgvp(mi)=drocp
            dgve(mi)=drocs
            dgle(mi)=drols
            rovf(mi)=roc
            denei(mi)=(dena-deneh(mi))*dtin
            dmef(mi)=ddens*dtin
            depf(mi)=ddenap*dtin
            deef(mi)=ddenas*dtin
         endif
c     source terms

         if(ieosd.eq.2) then
c     - - - - - PHS 2/08 changed to not
c     overwrite sk(mi) for evaporation nodes.
            if (evaporation_flag) then
               if(.not. evap_flag(mi)) sk(mi)=qwdis
            else
               sk(mi)=qwdis
            end if
            dq(mi)=dqwp
            if(irdof.ne.13) then
               dqh(mi)=dqws
               qh(mi)=qadis
               drc(mi)=dqap
               deqh(mi)=dqas
            endif
         else if(ieosd.eq.1) then
            sk(mi)=qwdis
            dq(mi)=dqwp
            if(irdof.ne.13) then
               dqh(mi)=0.0d00
               qh(mi)=qadis
               drc(mi)=dqap
               deqh(mi)=dqas
            endif
         else if(ieosd.eq.3) then
            sk(mi)= sx1(mi)*(sl-0.0)            
            dq(mi)=sx1(mi)
            if(irdof.ne.13) then
               dqh(mi)=0.0d00
               qh(mi)=qadis
               drc(mi)=dqap
               deqh(mi)=dqas
            endif
         endif
       if(kq.eq.1) then
c
c check for peaceman calculation
c       
            if(iwelimd.ne.0)then        
             if(izonewel1(mi).ne.0) then
              permsd = wellim(mi)*rolref_b/xvisl
             endif
             pflow(mi)=phi(mi) - sk(mi)/permsd
            endif
        endif      
 100        continue
c gaz 071819
c call routine to limit source change magitude
c  
c      call source_manage(1)
c      call source_manage(2)
c      
      if(allocated(s0)) then
         deallocate(s0,pcp0,dpcps0,rlf0)
         deallocate(drlfs0,rvf0,drvfs0)       
      endif
      do mi = 1,n
         if(nelm(mi)+1.eq.nelm(mi+1)) then
            sk(mi) = 0.0
            qh(mi) = 0.0
	 endif
      enddo
      return
      end
      subroutine source_manage(iflg)
      use comai
      use comdi
      implicit none
      integer iflg, i
      real*8 ratio, rat_tol      
      if(iflg.eq.0) then
c allocate space          
c save last step source    
       if(.not.allocated(sko)) allocate(sko(n))
       sko = sk
       elseif(iflg.eq.1) then
c limit max change in sk
        source_ratio_in = 5
        source_ratio_out = 2
        do i = 1, n
         if(sk(i).ne.0.0d0) then
          ratio = max(abs(sk(i))/abs(sko(i)+rat_tol),1.0)
          if(sk(i).gt.0.0) ratio = min(source_ratio_out,ratio)
          if(sk(i).lt.0.0) ratio = min(source_ratio_in,ratio)
          sk(i) = ratio*sk(i)
         endif
        enddo

       elseif(iflg.eq.2) then
c        deallocate(sko)
      endif
      end
      subroutine solubility_isothermal(iflg,mid)
c gaz 110923
c subroutine to calculate accurate saturations for gridblocks
c undergoing phase change
      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comevap, only : evaporation_flag, evap_flag
      use comfi, only : alpha0, q_gas, q_gas_tot, alpha_h2, alpha_co2,
     &         alpha_meth, alpha_air
      use comgi
      use comii
      use comrlp, only : rlpnew
      use com_prop_data, only : ctest, ihenryiso,xnl_max, xnl_ngas
      use comrxni
      use comwt
      use davidi
      implicit none
      integer iflg, mid, ii, i
      integer iphase_old, iphase_new
      real*8 pl, xnl,dxnlp,dxnlt,dxnls,alpha_tol 
c      real*8 ratio_mw_h2_air, ratio_mw_meth_air, ratio_mw_co2_air
      real*8 alpha_iso,p_air,dp_air_s,dp_air_p
c      parameter(ratio_mw_h2_air=0.07d0,ratio_mw_meth_air=0.55d0)
c      parameter(ratio_mw_co2_air= 1.52d0, alpha_h2 = 8.d-7)
      parameter(alpha_tol = 1.d-18)
       ctest = 1.0
       dp_air_p = 0.0d0
c      
      
      if(iflg.eq.0) then
c allocate memory
c xnl,xnv
       if(.not.allocated(xnl_ngas) ) then
         allocate(xnl_ngas(n0,4))
         xnl_ngas = 0.0d0
       endif
      else if(iflg.eq.1) then
c use Henry's law 
        if(ihenryiso.eq.0) then
         xnl_ngas(mid,1) = 0.0d0
         xnl_ngas(mid,2) = 0.0d0
         xnl_ngas(mid,3) = 0.0d0
         return
        endif
c gaz 111523 test 
         
        if(itype_air.ne.0) then
         alpha_air = alpha0
         alpha_iso = alpha_air
        else if(itype_co2.ne.0) then
          alpha_iso = alpha_co2
        else if(itype_meth.ne.0) then
          alpha_iso = alpha_meth
        else if(itype_h2.ne.0) then
          alpha_iso = alpha_h2
        endif 
c gaz 111723 gas pressure Pg = phi(mid) + pcap(mid) << depends on saturation
c        p_air = pcp(mid) + phi(mid)
c gaz 042024 gas pressure Pg = phi(mid)
c        p_air = pcp(mid) + phi(mid)
c        dp_air_s  = dpcef(mid) 
        p_air = phi(mid)
        dp_air_p =1.0d0
        dp_air_s  = 0.0d0
          if(s(mid).lt.1.d00) then
           xnl = alpha_iso*p_air + alpha_tol
           dxnlp =alpha_iso*dp_air_p  
           dxnlt =0.0d0
           dxnls= alpha_iso*dp_air_s
          else
c max xnl 
           xnl = alpha_iso*p_air + alpha_tol
           dxnlp=0.0
           dxnlt=0.0
           dxnls=0.0           
          endif
        xnl_ngas(mid,1) =  xnl*ctest
        xnl_ngas(mid,2) =  dxnlp*ctest
        xnl_ngas(mid,3) =  dxnls*ctest
      else if(iflg.eq.2) then
c gaz 122223 calculate Sl
c delp =pcp(s) calc S
c        call vgcap_inv_calc(2,mid) 
      else if(iflg.eq.3) then
c gaz 120423 gas diff boundary flow       
       q_gas_tot = 0.0
       do i = 1, n0
        if(ka(i).eq.-202) then
         q_gas_tot = q_gas_tot +q_gas(i)
        endif
       enddo
      endif
      end
      subroutine phase_change_mass_conv(iflg,mi1,mi2)
c subroutine to calculate mass fractions for gridblocks
c undergoing phase change
      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comevap, only : evaporation_flag, evap_flag
      use comfi
      use comgi
      use comii
      use comrlp, only : rlpnew
      use comrxni
      use comwt
      use davidi
      use com_prop_data, only : xnl_ngas, xnl_max
      implicit none
      integer iflg, mi1, mi2, ii, i, i1, i2, nr1, nr2
      integer iphase_old, iphase_new
      real*8 tol_frac, frac_gas, gas_mass, liq_mass, phi_old, sl
c gaz 041724
      real*8 den_vap,rolref,rol,pld,comw, sl_new, xnl
      parameter(tol_frac = 1.d-18)
c gaz notes 090719
c add arrays to comdi    
c pressure may decrease past 0.1 ; need starting pres for
c      mass (denity change)
      i =   l +iad
      if(iflg.eq.-1) then
c gaz 121823 this section not used yet
c allocate arrays if necessary
        allocate(s_prev(n))
        allocate(phi_prev(n))
        allocate(n_phase_nodes(n))
        allocate(ieos_prev(n))
        s_prev = 0.0
        phi_prev = 0.0
        n_phase_nodes = 0
        ieos_prev = 0
      else if(iflg.eq.0) then
c gaz 121823 allocate memory
       if(.not.allocated(frac_gas_iso)) then
        allocate(frac_gas_iso(n))
        frac_gas_iso = 0.0d0
       endif
      else if(iflg.eq.1) then  
c calculate mass fraction for noncondensible gas
       do i = mi1, mi2
        gas_mass = denei(i)*dtot+deneh(i)+tol_frac
        liq_mass = deni(i)*dtot+denh(i) 
        frac_gas = gas_mass/(gas_mass + liq_mass)
        frac_gas_iso(i) = frac_gas
       enddo
      else if(iflg.eq.2) then  
c NR update 2 phase isothermal  (P-S-xnl)
c comes through use module : strd
c
       nr1=nrhs(1)
       nr2=nrhs(2) 
       do i = mi1, mi2                           
         i1=i+nr1                 
         i2=i+nr2 
         sl = s(i)
         phi_old = phi(i)
c gaz 042124 
c         if(sl.ge.1.d0.or.ieos(i).ne.2) then 
         if(ieos(i).eq.1) then            
          phi(i)=phi_old-bp(i1)*strd
          cnlf(i) = cnlf(i)-bp(i2)*strd
         else
          phi(i)=phi_old-bp(i1)*strd
          s(i) = sl-bp(i2)*strd  
c          s(i) = min(1.d0,s(i))        
         endif
       enddo
       else if(iflg.eq.3) then  
c gaz 041724
c estimate saturation if phase change
c from solubilty gt sol max

       i = mi1
       if(s(i).ge.1.0) then
        xnl = cnlf(i)
        rolref=crl(1,1)
        comw=crl(3,1)
        pld  = phi(i)
        rol=rolref*(1.0+comw*(pld-pref))
        den_vap = roc0*(273.0/(tref+273.0))*pref/0.101325 
        sl_new = (den_vap-rol*xnl)/(den_vap-rol*xnl_max)
        s(i) = sl_new
        continue
       else
       endif 
      endif
      end                                    
