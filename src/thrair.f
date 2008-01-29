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
      use comrxni
      use davidi
      use comii
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comwt
      implicit none

      integer ndummy,mid,mi,ieosd,kq
      real*8 dtin,tempc,drocp,drocp0,rolref,xvisl,comw,pref,xvisv
      real*8 xrl,drl,drlp,xrv,drv,drvp,pl,sl,svd,pwl,dpwlp,dpwls
      real*8 qdis,qwdis,qadis,dqws,dqas,dqwp,dqap,por,sflux,permsd
      real*8 pflowd,area,uperm,roc,drocs,rol,drolp,dporpl
      real*8 rcomd,drols,dena,ddenp,ddens,ddenap,ddenas,dql,dqv
      real*8 pcl0,roc0,xvisl0,pdiff,plow,watfrac,dwfracp
      real*8 airterm,watterm,pdifa,pdifw
      real*8 tl_last,viln,vild,vil,tref
      real*8 viln1,viln2,viln3,vild1,vild2,vild3
      real*8 vla0,vlpa1,vlpa2,vlpa3,vlta1,vlta2,vlta3,vlpta,vlp2ta
      real*8 vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3,vlptb,vlp2tb
      real*8 vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tlx,tl2x,tlx2
      real*8 dsatp, rlpmin, qwmax_fac, qwmax
      real*8 pld,dis_ex,wat_ex,dwat_exs, time_max
      real*8 seep_facv,seep_facl,permsdv,permsdl,plwt
      integer i_mem_rlp
      integer iadka        
      save i_mem_rlp
      parameter(pcl0 = 0.101325)
      parameter(roc0 = 1.292864)
      parameter(rlpmin = 0.01)
      parameter(qwmax_fac = 1.)

      parameter(iadka=1000)                
      
      real*8, allocatable :: s0(:)
      real*8, allocatable :: pcp0(:)
      real*8, allocatable :: dpcps0(:)
      real*8, allocatable :: rlf0(:)
      real*8, allocatable :: drlfs0(:)
      real*8, allocatable :: rvf0(:)
      real*8, allocatable :: drvfs0(:)
      
      if(irdof.ne.13) then
         if(abs(iexrlp).ne.0.and.i_mem_rlp.eq.0) then
            i_mem_rlp=1
            allocate(s0(n),pcp0(n),dpcps0(n),rlf0(n))
            allocate(drlfs0(n),rvf0(n),drvfs0(n))       
         endif

c     get relative perms
         if(iad.lt.abs(iexrlp).or.abs(iexrlp).eq.0) then
            call rlperm(ndummy,1)
         else if(iad.eq.abs(iexrlp)) then
            call rlperm(ndummy,1)
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
         call cappr(1,ndummy)

      endif
c     
c     calculate variable porosity if enabled
c     
      if(iporos.ne.0) call porosi(1)
c     
c     dependent variables vap p and sl
c     
c     misc. constants
      tref = crl(6,1)
      tempc=273.0/(tref+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
      xvisl0=crl(2,1)
      comw=crl(3,1)
      pref=crl(4,1)
      xvisv=crl(5,1)
      rcomd=comw*rolref
      seep_facl = rolref/xvisl0
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
c     
      ifree1 = 0 
      do 100 mid=1,neq
         mi=mid+ndummy
         ieosd=2  
c count non fully saturated cells	    
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
	       
               permsdv = permsd*seep_facv
               qadis=permsdv*xrv*(pl-pflowd)
               dqap=permsdv*xrv + drvp*permsdv*(pl-pflowd)
               dqas= drv*permsdv*(pl-pflowd)
               if(pl-pflowd.gt.0.0) then
	          permsdl = permsd*seep_facl
                  qwdis=permsdl*xrl*(pl-pflowd)
                  dqwp=permsdl*xrl + drlp*permsdl*(pl-pflowd)
                  dqws= drl*permsdl*(pl-pflowd) 
               else
                  permsdl = permsd*seep_facl
                  qwdis=permsdl*xrl*(pl-pflowd)
                  dqwp=permsdl*xrl + drlp*permsdl*(pl-pflowd)
                  dqws= drl*permsdl*(pl-pflowd) 
               endif
               
            endif
            if(wellim(mi).lt.0.and.sflux.lt.0.0) then
               qadis = permsd*(pl-pflowd)
               dqap  = permsd 
               dqas = 0.0
               qwdis = permsd*(sl)
               dqws = permsd
               dqwp = 0.0d00     
            else if(wellim(mi).lt.0.and.pl-pflowd.lt.0.0) then
               qwdis = 0.0d00
               dqws = 0.0d00
               dqwp = 0.0d00
            endif
         else if(kq.eq.-2) then
c     ===Free drainage condition for water flow at the lower boundary
c     ===(each node at the free drainage boundary needs to be specified
c     ===explicitly, also the area or length associated with that node).
c     rol = rolref*(1.0+comw*(pl-pref))
            rol = rolref
            area=wellim(mi)
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
c     ===Seepage face condition for air and water 
            permsd=abs(wellim(mi))
            permsdv = permsd*seep_facv
            permsdl = permsd*seep_facl
            watterm= xrl*permsdl
            airterm= xrv*permsdv
            pflow(mi) = crl(4,1)
            pdifa = pl-pflow(mi)
            pdifw = pl-pflow(mi)
            if(pdifw.gt.0.0) then
               qwdis = watterm*pdifw
c     
               dqws = drl*permsdl*pdifw 
               dqwp = drlp*permsdl*pdifw + watterm 
            else
               qwdis = 0.0
               dqws = 0.0
               dqwp = 0.0
            endif
            qadis = airterm*(pdifa) 
            dqas = drv*permsdv*pdifa
            dqap = drvp*permsdv*pdifa + airterm


         else if(kq.eq.-3.and.ifree.ne.0) then
c     Seepage face for simpliflied water table
c     outflow only based on highest cell pressure (sat = 1.)
c     
            permsd= qc(mi)/wellim(mi)
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
c     outflow only based on highest cell pressure (sat = 1.)
c     based on average cell pressure
            permsd= qc(mi)/wellim(mi)
            permsdl = permsd*seep_facl
            plow = 0.5*(head12(mi,2)+head12(mi,1))
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
            
          else if(kq.eq.-5.and.ifree.ne.0) then
c Seepage face for simpliflied water table
c only flow when cell is full
c might need rlxyf here
            permsd=abs(wellim(mi))
            permsdl = permsd*seep_facl
            watterm= permsdl
            watfrac = max(rlxyf(mi)-rlptol,0.0d0)
            plow = head12(mi,2)          
            qwdis = watterm*watfrac*(pwl-plow)
            qwmax = esk(mi)
            dqwp = watterm*drlxyf(mi)*(pwl-plow)+watterm*watfrac
            if(qwdis.gt.qwmax) then
               qwdis = qwmax
               dqwp = 0.0
            else if(qwdis.lt.0.0d0) then
               qwdis = 0.
               dqwp = 0.0
            endif
         else if(kq.eq.-7.and.irdof.ne.13) then
c     ===Ponding condition for air and water 
c     ===area input in pflow term
            rol = rolref
            area=pflow(mi)
            uperm=-pnx(mi)*rol*rol*grav*area/xvisl
            
            wat_ex= (sl-1.0)/sl
            dwat_exs= 1.0/sl - (sl-1.0)/sl**2
            qwdis=uperm*wat_ex
            dqws=uperm*dwat_exs
            dqwp=0.0              
            qadis=0.
            dqas=0.
            dqap=0.
         else if(kq.eq.-22.and.irdof.ne.13) then
c     ===manage outflow only conditions
            if(sflux.lt.0.0) then
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
         else if(kq.eq.-22.and.irdof.eq.13) then
c     make sure there is a deivative wrt P for ifree ne 0
         else if(kq.eq.2) then
            qwdis = qc(mi)
            qadis = esk(mi)
            dqws  = 0.0
            dqwp = 0.0
            dqap  = 0.0
            dqas  = 0.0
         else if(kq.eq.-6) then
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
c outflow only if water present
               if(qdis.gt.0.0.and.esk(mi).lt.0.0) then
c	        permsd = sx1(mi)
                  qwdis = qdis*rlxyf(mi)
                  dqwp  = qdis*drlxyf(mi)
c		    qwdis = qdis
c	        dqwp = 0.0
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
         rol=rolref*(1.0+comw*(pld-pref))
         if(cden) rol = rol+factcden*anl((ispcden-1)*n0+mi)
         drolp=rcomd
         drols=0.0
c     accumulation terms
         den=por*rol*sl
         dena=por*roc*svd
         ddenp=por*sl*drolp+dporpl*rol*sl 
         ddens=por*rol
         ddenap=por*drocp*svd + dporpl*roc*svd
         ddenas=por*(drocs*svd-roc)         
c     
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
            rol=rolref
            drolp=0.0
            drols=0.0		 
         endif
c     transport terms
c     liquid
         dql=rol/xvisl*xrl
         dil(mi)=dql
         dilp(mi)=drolp/xvisl*xrl
         dglp(mi)=drolp
         rolf(mi)=rol
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
            sk(mi)=qwdis
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
 100  continue
      
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
