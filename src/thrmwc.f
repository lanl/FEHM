	subroutine thrmwc(ndummy)
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
CD1 To calculate the equation coeffients and derivatives for a
CD1 mixture of water and air, co2, or salt.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 NOV 1993     Z. Dash        22      Add prolog
CD2             G. Zyvoloski   N/A     Initial implementation
CD2
CD2  NOV 2014 G. Zyvoloski changed dqpc (gas source term ) to qng (q (n gas) ) for readability
CD2  Changed this everywhere in the code
CD2
CD2 $Log:   /pvcs.config/fehm90/src/thrmwc.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:44:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.15   Fri Nov 21 08:09:34 1997   gaz
CD2 corrected air flow term
CD2 
CD2    Rev 1.14   Fri Sep 26 15:13:56 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.13   Mon Mar 31 08:44:04 1997   gaz
CD2 major cleanup for source terms
CD2 
CD2    Rev 1.12   Fri Apr 26 16:36:28 1996   gaz
CD2 removed many lines (icons stuff now gone)
CD2 
CD2    Rev 1.11   Fri Mar 01 14:04:02 1996   gaz
CD2 added dpsatt=0.0, might help
CD2 
CD2    Rev 1.10   Fri Feb 16 13:01:56 1996   zvd
CD2 Added requirements.
CD2 
CD2    Rev 1.9   Wed Feb 07 12:04:28 1996   gaz
CD2 changes made for bous macro
CD2 
CD2    Rev 1.8   Fri Feb 02 12:52:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.7   08/07/95 13:48:00   awolf
CD2 NGAS fix - now can take prescribed saturation nodes.
CD2 
CD2    Rev 1.6   04/25/95 09:06:50   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.5   03/22/95 10:34:44   gaz
CD2 gaz : corrections made so humidity BC works
CD2 
CD2    Rev 1.4   03/10/95 11:08:08   llt
CD2 added humidity and reduced degree of air solution - gaz
CD2
CD2    Rev 1.3   01/28/95 13:56:08   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.2   08/22/94 16:22:24   robinson
CD2 BAR - Revised model for heat capacity of air
CD2 
CD2    Rev 1.1   03/18/94 16:12:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:50   pvcs
CD2 original version in process of being certified
CD2 
c  Note - this version combines thrmwc and thrmwcl (must check w/ GZ
c         and test before using)
c 11/21/94 gaz
c water balance equation mods
c 1/4/95 gaz workrd on source terms for water balance equations
c 1/5/95  gaz removed 1-xnl terms from transmissibilities
c 1/11/95 gaz cut off xnl at 0.0
c 1/11/95 gaz made derivatives ne 0 at cutoff for xnl,xnv
c 1/11/95 gaz made provision for just air source
c 1/12/95 gaz set enva(mi)=hcg, derivatives have similiar changes
c 1/17/95 gaz set enva(mi)=hcg-ens, derivatives have similiar changes
c 2/2/95  gaz set pointer to pnx for rlperms
c 2/7/95  gaz undid pointer above
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 ndummy       int     I       Parameter denoting whether primary of
CD3                                 matrix nodes are being treated.
CD3                                  
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 None
CD3 
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 n0, ipdeef, ipdepf, ipdmpf, ipdmef, ipdq, ipdqh, dtot, neq, 
CD4 qh, qc, dq, dqt, dqh, deqh, dcqh, dqc, deqc, dcqc, cel, crl, cvl,
CD4 cev, crv, t, pci, phi, pcp, dpcef, mw_air, mw_water, avgmolwt,
CD4 ice, ices, sii, tmelt, s, sk, denr, cpr, ps, ka, volume, pflow,
CD4 wellim, eflow, dmef, dmpf, dmc, depf, deef, dec, dcp, dce, dcc,
CD4 dtpa, dtpac, dtpae, a, dstm, dil, dilp, dile, enlf, dglp, dgle,
CD4 delf, delef, delcf, dclf, dclef, dclcf, dilc, dglc, dqv, div, divp, 
CD4 dive, envf, dgvp, dgve, dgvf, devf, devef, devcf, dcvf, dcvef, 
CD4 dcvcf, divc, dgvc, rolf, rovf, cnlf, cnvf, enva, denvap, denvae,
CD4 denvac, eskc, fimp, sx1, qflux, qflxm, ivfcal, denh, deni, denei,
CD4 denpch, denpci, ivapl, sv, den, dene
CD4 dqpc
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 psat            real*8   Computes saturation temperature or pressure
CD4 psatl           real*8   Computes saturation temperature or
CD4                             pressure with vapor pressure lowering
CD4 rlperm          N/A      generates relative permeabilities
CD4 cappr           N/A      generates capillary pressures
CD4 dvcalc          N/A      computes water vapor diffusion contribution
CD4 vfcal           real*8   Computes volume changes
CD4 air_cp          N/A      Computes air heat capacity and derivative
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 iieosl       int         Equation of state flag
CD5 mid          int         Do loop index over all nodes
CD5 mi           int         Current node number
CD5 ieosd        int         Phase state of fluid at current node
CD5 iieosd       int         Thermodynamics set at current node
CD5 xrl          real*8      Liquid relative permeability
CD5 xrv          real*8      Vapor relative permeability
CD5 drl          real*8      Derivative of liquid rel. perm.
CD5 drv          real*8      Derivative of vapor rel. perm.
CD5 drlp         real*8      Derivative of liquid rel. perm.
CD5 drvp         real*8      Derivative of vapor rel. perm.
CD5 rlf          real*8      Liquid relative permeability, each node
CD5 rvf          real*8      Vapor relative permeability, each node
CD5 drlef        real*8      Derivative of liquid rel. perm., each node
CD5 drvef        real*8      Derivative of vapor rel. perm., each node
CD5 drlpf        real*8      Derivative of liquid rel. perm., each node
CD5 drvpf        real*8      Derivative of vapor rel. perm., each node
CD5 ela0         real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpa1        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpa2        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpa3        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elta1        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elta2        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elta3        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpta        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elp2ta       real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpt2a       real*8      Polynomial coeff. for liquid water enthalpy
CD5 elb0         real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpb1        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpb2        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpb3        real*8      Polynomial coeff. for liquid water enthalpy
CD5 eltb1        real*8      Polynomial coeff. for liquid water enthalpy
CD5 eltb2        real*8      Polynomial coeff. for liquid water enthalpy
CD5 eltb3        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elptb        real*8      Polynomial coeff. for liquid water enthalpy
CD5 elp2tb       real*8      Polynomial coeff. for liquid water enthalpy
CD5 elpt2b       real*8      Polynomial coeff. for liquid water enthalpy
CD5 dla0         real*8      Polynomial coeff. for liquid water density
CD5 dlpa1        real*8      Polynomial coeff. for liquid water density
CD5 dlpa2        real*8      Polynomial coeff. for liquid water density
CD5 dlpa3        real*8      Polynomial coeff. for liquid water density
CD5 dlta1        real*8      Polynomial coeff. for liquid water density
CD5 dlta2        real*8      Polynomial coeff. for liquid water density
CD5 dlta3        real*8      Polynomial coeff. for liquid water density
CD5 dlpta        real*8      Polynomial coeff. for liquid water density
CD5 dlp2ta       real*8      Polynomial coeff. for liquid water density
CD5 dlpt2a       real*8      Polynomial coeff. for liquid water density
CD5 dlb0         real*8      Polynomial coeff. for liquid water density
CD5 dlpb1        real*8      Polynomial coeff. for liquid water density
CD5 dlpb2        real*8      Polynomial coeff. for liquid water density
CD5 dlpb3        real*8      Polynomial coeff. for liquid water density
CD5 dltb1        real*8      Polynomial coeff. for liquid water density
CD5 dltb2        real*8      Polynomial coeff. for liquid water density
CD5 dltb3        real*8      Polynomial coeff. for liquid water density
CD5 dlptb        real*8      Polynomial coeff. for liquid water density
CD5 dlp2tb       real*8      Polynomial coeff. for liquid water density
CD5 dlpt2b       real*8      Polynomial coeff. for liquid water density
CD5 vla0         real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpa1        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpa2        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpa3        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlta1        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlta2        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlta3        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpta        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlp2ta       real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpt2a       real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlb0         real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpb1        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpb2        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpb3        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vltb1        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vltb2        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vltb3        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlptb        real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlp2tb       real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 vlpt2b       real*8      Polynomial coeff. for liquid water
CD5                             viscosity
CD5 eva0         real*8      Polynomial coeff. for vapor enthalpy
CD5 evpa1        real*8      Polynomial coeff. for vapor enthalpy
CD5 evpa2        real*8      Polynomial coeff. for vapor enthalpy
CD5 evpa3        real*8      Polynomial coeff. for vapor enthalpy
CD5 evta1        real*8      Polynomial coeff. for vapor enthalpy
CD5 evta2        real*8      Polynomial coeff. for vapor enthalpy
CD5 evta3        real*8      Polynomial coeff. for vapor enthalpy
CD5 evpta        real*8      Polynomial coeff. for vapor enthalpy
CD5 evp2ta       real*8      Polynomial coeff. for vapor enthalpy
CD5 evpt2a       real*8      Polynomial coeff. for vapor enthalpy
CD5 evb0         real*8      Polynomial coeff. for vapor enthalpy
CD5 evpb1        real*8      Polynomial coeff. for vapor enthalpy
CD5 evpb2        real*8      Polynomial coeff. for vapor enthalpy
CD5 evpb3        real*8      Polynomial coeff. for vapor enthalpy
CD5 evtb1        real*8      Polynomial coeff. for vapor enthalpy
CD5 evtb2        real*8      Polynomial coeff. for vapor enthalpy
CD5 evtb3        real*8      Polynomial coeff. for vapor enthalpy
CD5 evptb        real*8      Polynomial coeff. for vapor enthalpy
CD5 evp2tb       real*8      Polynomial coeff. for vapor enthalpy
CD5 evpt2b       real*8      Polynomial coeff. for vapor enthalpy
CD5 dva0         real*8      Polynomial coeff. for vapor density
CD5 dvpa1        real*8      Polynomial coeff. for vapor density
CD5 dvpa2        real*8      Polynomial coeff. for vapor density
CD5 dvpa3        real*8      Polynomial coeff. for vapor density
CD5 dvta1        real*8      Polynomial coeff. for vapor density
CD5 dvta2        real*8      Polynomial coeff. for vapor density
CD5 dvta3        real*8      Polynomial coeff. for vapor density
CD5 dvpta        real*8      Polynomial coeff. for vapor density
CD5 dvp2ta       real*8      Polynomial coeff. for vapor density
CD5 dvpt2a       real*8      Polynomial coeff. for vapor density
CD5 dvb0         real*8      Polynomial coeff. for vapor density
CD5 dvpb1        real*8      Polynomial coeff. for vapor density
CD5 dvpb2        real*8      Polynomial coeff. for vapor density
CD5 dvpb3        real*8      Polynomial coeff. for vapor density
CD5 dvtb1        real*8      Polynomial coeff. for vapor density
CD5 dvtb2        real*8      Polynomial coeff. for vapor density
CD5 dvtb3        real*8      Polynomial coeff. for vapor density
CD5 dvptb        real*8      Polynomial coeff. for vapor density
CD5 dvp2tb       real*8      Polynomial coeff. for vapor density
CD5 dvpt2b       real*8      Polynomial coeff. for vapor density
CD5 vva0         real*8      Polynomial coeff. for vapor viscosity
CD5 vvpa1        real*8      Polynomial coeff. for vapor viscosity
CD5 vvpa2        real*8      Polynomial coeff. for vapor viscosity
CD5 vvpa3        real*8      Polynomial coeff. for vapor viscosity
CD5 vvta1        real*8      Polynomial coeff. for vapor viscosity
CD5 vvta2        real*8      Polynomial coeff. for vapor viscosity
CD5 vvta3        real*8      Polynomial coeff. for vapor viscosity
CD5 vvpta        real*8      Polynomial coeff. for vapor viscosity
CD5 vvp2ta       real*8      Polynomial coeff. for vapor viscosity
CD5 vvpt2a       real*8      Polynomial coeff. for vapor viscosity
CD5 vvb0         real*8      Polynomial coeff. for vapor viscosity
CD5 vvpb1        real*8      Polynomial coeff. for vapor viscosity
CD5 vvpb2        real*8      Polynomial coeff. for vapor viscosity
CD5 vvpb3        real*8      Polynomial coeff. for vapor viscosity
CD5 vvtb1        real*8      Polynomial coeff. for vapor viscosity
CD5 vvtb2        real*8      Polynomial coeff. for vapor viscosity
CD5 vvtb3        real*8      Polynomial coeff. for vapor viscosity
CD5 vvptb        real*8      Polynomial coeff. for vapor viscosity
CD5 vvp2tb       real*8      Polynomial coeff. for vapor viscosity
CD5 vvpt2b       real*8      Polynomial coeff. for vapor viscosity
CD5 tl           real*8      Temperature at this node
CD5 pcl          real*8      Gas pressure at this node
CD5 pl           real*8      Pressure at this node
CD5 x            real*8      Term used in polynomial calculation
CD5 x2           real*8      x squared
CD5 x3           real*8      x cubed
CD5 dpsats       real*8      Derivative of saturation pressure with
CD5                             saturation
CD5 dpsatt       real*8      Derivative of saturation pressure with
CD5                             temperature
CD5 dpct         real*8      Negative of dpsatt
CD5 pv           real*8      Vapor pressure
CD5 xv           real*8      Vapor pressure
CD5 xv2          real*8      xv squared
CD5 xv3          real*8      xv cubed
CD5 dtps         real*8      Term used in calculation
CD5 tl2          real*8      tl squared
CD5 tl3          real*8      tl cubed
CD5 tlx          real*8      x times tl
CD5 tl2x         real*8      tl2 times x
CD5 tlx2         real*8      tl times x2
CD5 tlxv         real*8      xv times tl
CD5 tl2xv        real*8      tl squared time xv
CD5 tlxv2        real**      tl times xv squared
CD5 enwn1        real*8      Term in liquid enthalpy calculation
CD5 enwn2        real*8      Term in liquid enthalpy calculation
CD5 enwn3        real*8      Term in liquid enthalpy calculation
CD5 enwn         real*8      Term in liquid enthalpy calculation
CD5 enwd1        real*8      Term in liquid enthalpy calculation
CD5 enwd2        real*8      Term in liquid enthalpy calculation
CD5 enwd3        real*8      Term in liquid enthalpy calculation
CD5 enwd         real*8      Term in liquid enthalpy calculation
CD5 enw          real*8      Liquid water enthalpy
CD5 enl          real*8      Liquid enthalpy
CD5 env          real*8      Vapor enthalpy
CD5 dhwpn1       real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwpn2       real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwpn        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwpd        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwp         real*8      Water enthalpy derivative with pressure
CD5 dhwtn1       real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwtn2       real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwtn        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwtd        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhwt         real*8      Water enthalpy derivative with temperature
CD5 dhlt         real*8      Liquid enthalpy derivative with temperature
CD5 dhlp         real*8      Liquid enthalpy derivative with pressure
CD5 dhlpc        real*8      Liquid enthalpy derivative with capillary
CD5                             pressure
CD5 rnwn1        real*8      Term used in liquid density calculation
CD5 rnwn2        real*8      Term used in liquid density calculation
CD5 rnwn3        real*8      Term used in liquid density calculation
CD5 rnwd1        real*8      Term used in liquid density calculation
CD5 rnwd2        real*8      Term used in liquid density calculation
CD5 rnwd3        real*8      Term used in liquid density calculation
CD5 rnwn         real*8      Term used in liquid density calculation
CD5 rnwd         real*8      Term used in liquid density calculation
CD5 rnw          real*8      Liquid water density
CD5 rol          real*8      Liquid density
CD5 drlpn1       real*8      Term used in derivative of density
CD5                              calculation
CD5 drlpn2       real*8      Term used in derivative of density
CD5                              calculation
CD5 drlpn        real*8      Term used in derivative of density
CD5                              calculation
CD5 drolpd       real*8      Term used in derivative of density
CD5                              calculation
CD5 drolp        real*8      Derivative of liquid density with pressure
CD5 drlen1       real*8      Term used in derivative of density
CD5                              calculation
CD5 drlen2       real*8      Term used in derivative of density
CD5                              calculation
CD5 drlen        real*8      Term used in derivative of density
CD5                              calculation
CD5 droled       real*8      Term used in derivative of density
CD5                              calculation
CD5 drolt        real*8      Derivative of liquid density with
CD5                              temperature
CD5 drolpc       real*8      Derivative of liquid density with
CD5                              capillary pressure
CD5 viln1       real*8       Term used in liquid viscosity calculation
CD5 viln2       real*8       Term used in liquid viscosity calculation
CD5 viln3       real*8       Term used in liquid viscosity calculation
CD5 viln        real*8       Term used in liquid viscosity calculation
CD5 vild1       real*8       Term used in liquid viscosity calculation
CD5 vild2       real*8       Term used in liquid viscosity calculation
CD5 vild3       real*8       Term used in liquid viscosity calculation
CD5 vild        real*8       Term used in liquid viscosity calculation
CD5 vil         real*8       Water viscosity
CD5 xvisl       real*8       Liquid viscosity
CD5 dvlpn1      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvlpn2      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvlpn       real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvilpd      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvislp      real*8       Derivative of liquid viscosity with
CD5                              pressure
CD5 dvlen1      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvlen2      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvlen       real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dviled      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvislt      real*8       Derivative of liquid viscosity with
CD5                              temperature
CD5 dvlpc       real*8       Derivative of liquid viscosity with
CD5                              capillary pressure
CD5 alpha       real*8       Inverse of Henry's law constant for air
CD5                              in water
CD5 dalpca      real*8       Derivative of Henry's constant with
CD5                              pressure
CD5 dxnlp       real*8       Derivative of air mass fraction with
CD5                              pressure
CD5 dxnlpc      real*8       Derivative of air mass fraction with
CD5                              capillary pressure
CD5 dxnlt       real*8       Derivative of air mass fraction with
CD5                              temperature
CD5 xnl         real*8       Mass fraction of air in water
CD5 xnv         real*8       Mass fraction of air in vapor
CD5 hsol        real*8       Heat of solution of air in water
CD5 cpa         real*8       Air heat capacity
CD5 dcpat       real*8       Derivative of air heat capacity with
CD5                             temperature
CD5 hcg         real*8       Air enthalpy
CD5 hcl         real*8       Heat capacity of air
CD5 dhsolt      real*8       Derivative of heat of solution with
CD5                              temperature
CD5 dhsolp      real*8       Derivative of heat of solution with
CD5                              pressure
CD5 dhcgt       real*8       Derivative of air enthalpy with
CD5                              temperature
CD5 dhcgp       real*8       Derivative of air enthalpy with
CD5                              pressure
CD5 dhslpc      real*8       Derivative of heat of solution with
CD5                              capillary pressure
CD5 ensn1        real*8      Term in vapor enthalpy calculation
CD5 ensn2        real*8      Term in vapor enthalpy calculation
CD5 ensn3        real*8      Term in vapor enthalpy calculation
CD5 ensn         real*8      Term in vapor enthalpy calculation
CD5 ensd1        real*8      Term in vapor enthalpy calculation
CD5 ensd2        real*8      Term in vapor enthalpy calculation
CD5 ensd3        real*8      Term in vapor enthalpy calculation
CD5 ensd         real*8      Term in vapor enthalpy calculation
CD5 ens          real*8      Vapor water enthalpy
CD5 dhvp1        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvp2        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvpn        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvpd        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvt1        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvt2        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvtn        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvtd        real*8      Term used in derivative of enthalpy
CD5                              calculation
CD5 dhvt         real*8      Water enthalpy derivative with temperature
CD5 dhvp         real*8      Liquid enthalpy derivative with pressure
CD5 dhvpc        real*8      Liquid enthalpy derivative with capillary
CD5                             pressure
CD5 rnsn1        real*8      Term used in vapor density calculation
CD5 rnsn2        real*8      Term used in vapor density calculation
CD5 rnsn3        real*8      Term used in vapor density calculation
CD5 rnsd1        real*8      Term used in vapor density calculation
CD5 rnsd2        real*8      Term used in vapor density calculation
CD5 rnsd3        real*8      Term used in vapor density calculation
CD5 rnsn         real*8      Term used in vapor density calculation
CD5 rnsd         real*8      Term used in vapor density calculation
CD5 ros          real*8      Density of water vapor
CD5 rns          real*8      Density of steam
CD5 rov          real*8      Density of vapor phase
CD5 drspn1       real*8      Term used in derivative of density
CD5                              calculation
CD5 drspn2       real*8      Term used in derivative of density
CD5                              calculation
CD5 drspn        real*8      Term used in derivative of density
CD5                              calculation
CD5 drospd       real*8      Term used in derivative of density
CD5                              calculation
CD5 drsen1       real*8      Term used in derivative of density
CD5                              calculation
CD5 drsen2       real*8      Term used in derivative of density
CD5                              calculation
CD5 drsen        real*8      Term used in derivative of density
CD5                              calculation
CD5 drostd       real*8      Term used in derivative of density
CD5                              calculation
CD5 visn1       real*8       Term used in vapor viscosity calculation
CD5 visn2       real*8       Term used in vapor viscosity calculation
CD5 visn3       real*8       Term used in vapor viscosity calculation
CD5 visn        real*8       Term used in vapor viscosity calculation
CD5 visd1       real*8       Term used in vapor viscosity calculation
CD5 visd2       real*8       Term used in vapor viscosity calculation
CD5 visd3       real*8       Term used in vapor viscosity calculation
CD5 visd        real*8       Term used in vapor viscosity calculation
CD5 vis         real*8       Water vapor viscosity
CD5 xvisv       real*8       Vapor viscosity
CD5 dvspn1      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvspn2      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvspn       real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvispd      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvisvp      real*8       Derivative of vapor viscosity with
CD5                              pressure
CD5 dvsen1      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvsen2      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvsen       real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvised      real*8       Term used in derivative of viscosity
CD5                              calculation
CD5 dvisvt      real*8       Derivative of vapor viscosity with
CD5                              temperature
CD5 dvvpc       real*8       Derivative of vapor viscosity with
CD5                              capillary pressure
CD5 pcl0        real*8       1 atm in MPa
CD5 roc0        real*8       Density of air at std. conditions
CD5 drocpc      real*8       Derivative of air with pressure
CD5 roc         real*8       Density of air
CD5 rocmol      real*8       Molar density of air
CD5 rosmol      real*8       Molar density of water
CD5 rovmol      real*8       Molar density of vapor
CD5 xnvmol      real*8       Mole fraction of air in vapor
CD5 drocp       real*8       Derivative of air density with pressure
CD5 droct       real*8       Derivative of air density with temperature
CD5 drovp       real*8       Derivative of vapor density with pressure
CD5 drovt       real*8       Derivative of vapor density with
CD5                              temperature
CD5 drovpc      real*8       Derivative of vapor density with
CD5                              capillary pressure
CD5 dxnvt       real*8       Derivative of air mass fraction with
CD5                              temperature
CD5 dxnvp       real*8       Derivative of air mass fraction with
CD5                              pressure
CD5 dxnvpc      real*8       Derivative of air mass fraction with
CD5                              capillary pressure
CD5 dhvsp       real*8       Derivative of air enthalpy with pressure
CD5 dhvst       real*8       Derivative of air enthalpy with temperature
CD5 dhvspc      real*8       Derivative of air enthalpy with capillary
CD5                              pressure
CD5 dhcgpc      real*8       Derivative of air enthalpy with capillary
CD5                              pressure
CD5 xvisc       real*8       Air viscosity
CD5 dxvsct      real*8       Derivative of air viscosity with
CD5                              temperature
CD5 dtpcs       real*8       Derivatice od saturation with temperature
CD5 siid        real*8       Ice saturation at this node
CD5 siie        real*8       Water saturation (when ice is present) at
CD5                              this node
CD5 roll        real*8       Liquid density
CD5 roli        real*8       constant of 950
CD5 enll        real*8       Liquid enthalpy
CD5 visll       real*8       Liquid viscosity
CD5 visli       real*8       Constant of 1.e20
CD5 xnll        real*8       Mass fraction of air in water
CD5 cph         real*8       Water heat capacity, std. conditions
CD5 drolpc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 drovpc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dhlpc0      real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dhvpc0      real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dvlpc0      real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dvvpc0      real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dxnlpc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dxnvpc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dhcgpc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dhvspc0     real*8       Term used in derivative calculations for
CD5                              ice solution
CD5 dtd         real*8       Factor used in derivative calculations
CD5 sl          real*8       Saturation at this node
CD5 qdis        real*8       Water source flow rate at this node
CD5 qdist       real*8       Total source flow rate at this node
CD5 dqdistp     real*8       Derivative of total flow wrt pressure    
CD5 dqdiste     real*8       Derivative of total flow wrt temperture  
CD5 dqdistc     real*8       Derivative of total flow wrt air variable
CD5 cp          real*8       Density time heat capacity
CD5 por         real*8       Porosity at this node
CD5 kq          int          Flag used in flow rate calculation
CD5 vol         real*8       Volume of this node
CD5 pldif       real*8       Delta P driving force for sink or source
CD5 permsd      real*8       Well impedance at current node
CD5 eskd        real*8       Energy flow from source at this node
CD5 eqdum       real*8       Energy stored at this node
CD5 denc        real*8       Total density at this node
CD5 rag         real*8       Ratio of vapor kinematic viscosities
CD5 sig         real*8       Ratio used in calculation
CD5 rop         real*8       Derivative of den with pressure
CD5 damp        real*8       rop per unit time
CD5 daep        real*8       Parameter used in calculation
CD5 dacp        real*8       Parameter used in calculation
CD5 damh        real*8       Parameter used in calculation
CD5 daeh        real*8       Parameter used in calculation
CD5 ropc        real*8       Parameter used in calculation
CD5 daepc       real*8       Parameter used in calculation
CD5 dacpc       real*8       Parameter used in calculation
CD5 dragp       real*8       Parameter used in calculation
CD5 drage       real*8       Parameter used in calculation
CD5 dragpc      real*8       Parameter used in calculation
CD5 dsigp       real*8       Derivative of sig with pressure
CD5 dsige       real*8       Derivative of sig with enthalpy
CD5 dsigpc      real*8       Derivative of sig with capillary pressure
CD5 den1        real*8       Liquid density
CD5 roe         real*8       Term used in calculation
CD5 dampc       real*8       Term used in calculation
CD5 dql         real*8       Term used in calculation
CD5 hprod       real*8       Term used in calculation
CD5 cprod       real*8       Term used in calculation
CD5 dhprdp      real*8       Derivative term used in calculation
CD5 dhprdc      real*8       Derivative term used in calculation
CD5 dhprde      real*8       Derivative term used in calculation
CD5 dcprdp      real*8       Derivative term used in calculation
CD5 dcprdc      real*8       Derivative term used in calculation
CD5 dcprde      real*8       Derivative term used in calculation
CD5 eskcd       real*8       Pressure  term for gas at this node
CD5 dqpc        real*8       Source term for gas at this node
CD5 pbounc      real*8       Negative of source term for gas at this
CD5                             node
CD5 htc         real*8       Heat flux impedance of source at this node
CD5 tbound      real*8       Heat flux of source at this node
CD5 hflux       real*8       Heat flux at this node
CD5 dhflxp      real*8       Derivative of heat flux
CD5 dhflxe      real*8       Derivative of heat flux
CD5 dhflxc      real*8       Derivative of heat flux
CD5 cprd        real*8       Heat capacity at this node
CD5 edif        real*8       Term used in calculation
CD5 denrd       real*8       Rock density heat capacity product
CD5 vfd         real*8       Term used in calculation
CD5 rl          real*8       Parameter used in calculation
CD5 dvfp        real*8       Parameter used in calculation
CD5 dtin        real*8       Reciprocal of time step
CD5 dach        real*8       Term used in calculation
CD5 dqv         real*8       Term used in calculation
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
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN thrmwcl
CPS 
CPS rlperm - generate relative permeabilities
CPS cappr - generate capillary pressures
CPS 
CPS FOR each node
CPS 
CPS   Initialize variables
CPS   Initialize flow terms
CPS   
CPS   IF thermodynamic coefficients need to be assigned
CPS   
CPS     IF coefficients are needed
CPS       Assign liquid enthalpy coefficients
CPS       Assign liquid density coefficients
CPS       Assign liquid viscosity coefficients
CPS       Assign vapor enthalpy coefficients
CPS       Assign vapor density coefficients
CPS       Assign vapor viscosity coefficients
CPS     ENDIF
CPS     
CPS     IF we are not in the superheated region
CPS       Assign parameter values
CPS     ENDIF
CPS     
CPS     IF this is a two-phase mixture
CPS       IF there is vapor pressure lowering
CPS         psatl - Compute partial pressure of gas in vapor and
CPS             derivative
CPS       ELSE
CPS         psat - Compute partial pressure of gas in vapor and
CPS             derivative
CPS       ENDIF
CPS     ENDIF
CPS     
CPS     IF this is not a one-phase liquid
CPS       Compute vapor pressure and terms used in derivative...
CPS       ... calculations
CPS     ENDIF
CPS     
CPS     IF this is not a superheated vapor
CPS       Compute terms used in derivative calculations
CPS     ENDIF
CPS   
CPS     IF this is not a compressed liquid
CPS       Compute terms used in derivative calculations
CPS     ENDIF
CPS     
CPS     IF we are not superheated
CPS       Compute liquid enthalpy
CPS       Compute derivatives of liquid enthalpy
CPS       Compute liquid density
CPS       Compute derivatives of liquid density
CPS       Compute liquid viscosity
CPS       Compute derivatives of liquid viscosity
CPS       
CPS       IF we are a compressed liquid
CPS         Set derivatives with respect to capillary pressure to 0
CPS       ENDIF
CPS       
CPS       Calculate mass fraction of air in liquid and derivatives
CPS       
CPS       IF this is a vapor pressure lowering simulation
CPS         IF computed value of mass fraction is less than 0
CPS           Set to 0
CPS         ENDIF
CPS       ENDIF
CPS       
CPS       air_cp - calculate air enthalpy and derivatives
CPS       
CPS     ENDIF
CPS   
CPS   ENDIF thermdynamic coefficients need to be assigned
CPS   
CPS   IF this is not a compressed liquid
CPS   
CPS     Compute vapor enthalpy
CPS     Compute derivatives of vapor enthalpy
CPS     Compute vapor density
CPS     Compute derivatives of vapor density
CPS     Compute vapor viscosity
CPS     Compute derivatives of vapor viscosity
CPS     
CPS     Compute air density, mass fraction, average molecular weight
CPS     Compute derivatives of air and mixture density, mass fraction
CPS     air_cp - compute enthalpy of air, mixture, and derivatives
CPS     Compute vapor viscosity (air and water vapor), derivatives
CPS   
CPS   ELSE we must still set the average molecular weight of vapor
CPS     Set average molecular weight to the value for air
CPS   ENDIF this is not a compressed liquid
CPS 
CPS   IF this is a two-phase mixture
CPS     IF this simulation involves ice
CPS     
CPS       IF liquid is not fully frozen
CPS         Modify values and derivatives
CPS       ENDIF
CPS     
CPS       IF liquid is fully frozen
CPS         Modify values and derivatives
CPS       ENDIF
CPS       
CPS     ENDIF
CPS     
CPS     Compute additional derivative terms
CPS     IF vapor pressure lowering is included
CPS       Compute derivatives for this special case
CPS     ELSE
CPS       Zero out these derivative terms
CPS     ENDIF
CPS     
CPS   ENDIF this is a two-phase mixture
CPS     
CPS   IF a pressure-dependent source/sink is included
CPS     Compute term for mass and energy sink
CPS     IF this is a fluid source
CPS       Recompute source of mass and energy
CPS     ENDIF
CPS   ENDIF
CPS   
CPS   IF this is a two-phase mixture
CPS   
CPS     Compute accumulation terms
CPS     Compute Steam production
CPS     Compute derivatives of accumulation terms
CPS     
CPS     IF there is a fluid sink here
CPS       Compute derivative terms
CPS       IF the vapor contains both steam and gas
CPS         Compute derivative terms
CPS       ELSE
CPS         Zero out these derivative terms
CPS       ENDIF
CPS     ENDIF
CPS   
CPS   ENDIF this is a two-phase mixture
CPS   
CPS   IF this is a compressed liquid
CPS   
CPS     Zero out vapor derivative terms
CPS     Compute accumulation terms
CPS     Compute derivatives of accumulation terms
CPS   
CPS   ENDIF this is a compressed liquid
CPS     
CPS   IF this is a superheated vapor
CPS   
CPS     Compute accumulation terms
CPS     Compute derivatives of accumulation terms
CPS   
CPS   ENDIF this is a superheated vapor
CPS   
CPS   IF this is not a heat flow only problem
CPS     Store derivatives of accumulation terms in arrays
CPS
CPS     IF this is not a compressed liquid
CPS       Store derivatives of terms involving vapor in arrays
CPS     ENDIF this is not a compressed liquid
CPS   
CPS     Store density terms in arrays
CPS   
CPS     IF there is a sink at this node
CPS   
CPS       IF this is a two-phase air-water system
CPS         Compute derivative terms
CPS       ENDIF
CPS   
CPS       IF this is a compressed liquid system
CPS         Compute derivative terms
CPS       ENDIF
CPS   
CPS       IF this is a two-phase system
CPS         Compute derivative terms
CPS       ENDIF
CPS     
CPS       Store sink terms in arrays
CPS   
CPS     ENDIF there is a sink at this node
CPS     
CPS     IF there is a source at this node
CPS       Store source terms in arrays
CPS     ENDIF there is a source at this node
CPS     
CPS   ENDIF this is not a heat flow only problem
CPS   
CPS   IF there is a specified noncondensible pressure boundary
CPS   
CPS     IF this is not an air-water isothermal simulation
CPS       Compute and store source/sink terms in arrays
CPS     ELSE
CPS       Compute and store source/sink terms in arrays
CPS     ENDIF
CPS   
CPS   ENDIF there is a specified noncondensible pressure boundary
CPS   
CPS   IF there is a heat source/sink term at this node
CPS   
CPS     Compute heat flux
CPS     IF the term is a sink
CPS       IF this is not an isothermal air-water system
CPS         Set derivative terms
CPS       ELSE
CPS         Set derivative terms
CPS       ENDIF
CPS       
CPS       Set values in arrays
CPS       
CPS     ELSEIF the term is a sink
CPS       IF this is not an isothermal air-water system
CPS         Set derivative term
CPS       ELSE
CPS         Set derivative term
CPS       ENDIF
CPS       
CPS       Set values in arrays
CPS       
CPS     ELSE
CPS       Set value in array
CPS     ENDIF
CPS   
CPS   ENDIF there is a heat source/sink term at this node
CPS     
CPS   IF this is a heat conduction only simulation
CPS     IF there is a heat source at this node
CPS       Set values in arrays
CPS     ENDIF
CPS     Set derivative terms
CPS   ENDIF
CPS   
CPS ENDFOR
CPS 
CPS IF there are volume changes in this simulation
CPS   FOR each node
CPS     Adjust derivative terms
CPS   ENDFOR
CPS ELSE
CPS   FOR each node
CPS     Set derivative values in arrays
CPS   ENDFOR
CPS ENDIF
CPS 
CPS dvcalc - compute water vapor diffusion contribution
CPS 
CPS END thrmwcl
CPS
C**********************************************************************
c
c this subroutine calculates equation coeffients and derivatives
c mixture of water and air,co2

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

      implicit none

      integer ndummy
      integer iieosl
      integer mid
      integer mi
      integer ieosd
      integer iieosd
      integer kq
      real*8 dtin
      real*8 psatl
      real*8 dach
      real*8 dqv
      real*8 vfcal
      real*8 xrl
      real*8 xrv
      real*8 drl
      real*8 drv
      real*8 drlp
      real*8 drvp
      real*8 ela0
      real*8 elpa1
      real*8 elpa2
      real*8 elpa3
      real*8 elta1
      real*8 elta2
      real*8 elta3
      real*8 elpta
      real*8 elp2ta
      real*8 elpt2a
      real*8 elb0
      real*8 elpb1
      real*8 elpb2
      real*8 elpb3
      real*8 eltb1
      real*8 eltb2
      real*8 eltb3
      real*8 elptb
      real*8 elp2tb
      real*8 elpt2b
      real*8 dla0
      real*8 dlpa1
      real*8 dlpa2
      real*8 dlpa3
      real*8 dlta1
      real*8 dlta2
      real*8 dlta3
      real*8 dlpta
      real*8 dlp2ta
      real*8 dlpt2a
      real*8 dlb0
      real*8 dlpb1
      real*8 dlpb2
      real*8 dlpb3
      real*8 dltb1
      real*8 dltb2
      real*8 dltb3
      real*8 dlptb
      real*8 dlp2tb
      real*8 dlpt2b
      real*8 vla0
      real*8 vlpa1
      real*8 vlpa2
      real*8 vlpa3
      real*8 vlta1
      real*8 vlta2
      real*8 vlta3
      real*8 vlpta
      real*8 vlp2ta
      real*8 vlpt2a
      real*8 vlb0
      real*8 vlpb1
      real*8 vlpb2
      real*8 vlpb3
      real*8 vltb1
      real*8 vltb2
      real*8 vltb3
      real*8 vlptb
      real*8 vlp2tb
      real*8 vlpt2b
      real*8 eva0
      real*8 evpa1
      real*8 evpa2
      real*8 evpa3
      real*8 evta1
      real*8 evta2
      real*8 evta3
      real*8 evpta
      real*8 evp2ta
      real*8 evpt2a
      real*8 evb0
      real*8 evpb1
      real*8 evpb2
      real*8 evpb3
      real*8 evtb1
      real*8 evtb2
      real*8 evtb3
      real*8 evptb
      real*8 evp2tb
      real*8 evpt2b
      real*8 dva0
      real*8 dvpa1
      real*8 dvpa2
      real*8 dvpa3
      real*8 dvta1
      real*8 dvta2
      real*8 dvta3
      real*8 dvpta
      real*8 dvp2ta
      real*8 dvpt2a
      real*8 dvb0
      real*8 dvpb1
      real*8 dvpb2
      real*8 dvpb3
      real*8 dvtb1
      real*8 dvtb2
      real*8 dvtb3
      real*8 dvptb
      real*8 dvp2tb
      real*8 dvpt2b
      real*8 vva0
      real*8 vvpa1
      real*8 vvpa2
      real*8 vvpa3
      real*8 vvta1
      real*8 vvta2
      real*8 vvta3
      real*8 vvpta
      real*8 vvp2ta
      real*8 vvpt2a
      real*8 vvb0
      real*8 vvpb1
      real*8 vvpb2
      real*8 vvpb3
      real*8 vvtb1
      real*8 vvtb2
      real*8 vvtb3
      real*8 vvptb
      real*8 vvp2tb
      real*8 vvpt2b
      real*8 tl
      real*8 pcl
      real*8 pl
      real*8 x
      real*8 x2
      real*8 x3
      real*8 dpsats
      real*8 dpsatt
      real*8 dpct
      real*8 pv
      real*8 xv
      real*8 xv2
      real*8 xv3
      real*8 dtps
      real*8 tl2
      real*8 tl3
      real*8 tlx
      real*8 tl2x
      real*8 tlx2
      real*8 tlxv
      real*8 tl2xv
      real*8 tlxv2
      real*8 enwn1
      real*8 enwn2
      real*8 enwn3
      real*8 enwn
      real*8 enwd1
      real*8 enwd2
      real*8 enwd3
      real*8 enwd
      real*8 enw
      real*8 enl
      real*8 dhwpn1
      real*8 dhwpn2
      real*8 dhwpn
      real*8 dhwpd
      real*8 dhwp
      real*8 dhwtn1
      real*8 dhwtn2
      real*8 dhwtn
      real*8 dhwtd
      real*8 dhwt
      real*8 dhlt
      real*8 dhlp
      real*8 dhlpc
      real*8 rnwn1
      real*8 rnwn2
      real*8 rnwn3
      real*8 rnwd1
      real*8 rnwd2
      real*8 rnwd3
      real*8 rnwn
      real*8 rnwd
      real*8 rnw
      real*8 rol
      real*8 drlpn1
      real*8 drlpn2
      real*8 drlpn
      real*8 drolpd
      real*8 drolp
      real*8 drlen1
      real*8 drlen2
      real*8 drlen
      real*8 droled
      real*8 drolt
      real*8 drolpc
      real*8 viln1
      real*8 viln2
      real*8 viln3
      real*8 viln
      real*8 vild1
      real*8 vild2
      real*8 vild3
      real*8 vild
      real*8 vil
      real*8 xvisl
      real*8 dvlpn1
      real*8 dvlpn2
      real*8 dvlpn
      real*8 dvilpd
      real*8 dvislp
      real*8 dvlen1
      real*8 dvlen2
      real*8 dvlen
      real*8 dviled
      real*8 dvislt
      real*8 dvlpc
      real*8 alpha
      real*8 dalpca
      real*8 dxnlp
      real*8 dxnlpc
      real*8 dxnlt
      real*8 xnl
      real*8 xnv
      real*8 hsol
      real*8 cpa
      real*8 dcpat
      real*8 hcg
      real*8 hcl
      real*8 dhsolt
      real*8 dhsolp
      real*8 dhcgt
      real*8 dhcgp
      real*8 dhslpc
      real*8 ensn1
      real*8 ensn2
      real*8 ensn3
      real*8 ensn
      real*8 ensd1
      real*8 ensd2
      real*8 ensd3
      real*8 ensd
      real*8 ens
      real*8 env
      real*8 dhvp1
      real*8 dhvp2
      real*8 dhvpn
      real*8 dhvpd
      real*8 dhvt1
      real*8 dhvt2
      real*8 dhvtn
      real*8 dhvtd
      real*8 dhvt
      real*8 dhvp
      real*8 dhvpc
      real*8 rnsn1
      real*8 rnsn2
      real*8 rnsn3
      real*8 rnsd1
      real*8 rnsd2
      real*8 rnsd3
      real*8 rnsn
      real*8 rnsd
      real*8 ros
      real*8 rns
      real*8 rov
      real*8 drspn1
      real*8 drspn2
      real*8 drspn
      real*8 drospd
      real*8 drsen1
      real*8 drsen2
      real*8 drsen
      real*8 drostd
      real*8 visn1
      real*8 visn2
      real*8 visn3
      real*8 visn
      real*8 visd1
      real*8 visd2
      real*8 visd3
      real*8 visd
      real*8 vis
      real*8 xvisv
      real*8 dvspn1
      real*8 dvspn2
      real*8 dvspn
      real*8 dvispd
      real*8 dvisvp
      real*8 dvsen1
      real*8 dvsen2
      real*8 dvsen
      real*8 dvised
      real*8 dvisvt
      real*8 dvvpc
      real*8 pcl0
c      real*8 roc0
      real*8 drocpc
      real*8 roc
      real*8 rocmol
      real*8 rosmol
      real*8 rovmol
      real*8 xnvmol
      real*8 drocp
      real*8 droct
      real*8 drovp
      real*8 drovt
      real*8 drovpc
      real*8 dxnvt
      real*8 dxnvp
      real*8 dxnvpc
      real*8 dhvsp
      real*8 dhvst
      real*8 dhvspc
      real*8 dhcgpc
      real*8 xvisc
      real*8 dxvsct
      real*8 dtpcs
      real*8 siid
      real*8 siie
      real*8 roll
      real*8 roli
      real*8 enll
      real*8 visll
      real*8 visli
      real*8 xnll
      real*8 cph
      real*8 drolpc0
      real*8 drovpc0
      real*8 dhlpc0
      real*8 dhvpc0
      real*8 dvlpc0
      real*8 dvvpc0
      real*8 dxnlpc0
      real*8 dxnvpc0
      real*8 dhcgpc0
      real*8 dhvspc0
      real*8 dtd
      real*8 sl
      real*8 qdis
      real*8 cp
      real*8 por
      real*8 vol
      real*8 pldif
      real*8 permsd, permsda
      real*8 eskd
      real*8 eqdum
      real*8 denc
      real*8 rag
      real*8 sig
      real*8 rop
      real*8 damp
      real*8 daep
      real*8 dacp
      real*8 damh
      real*8 daeh
      real*8 ropc
      real*8 daepc
      real*8 dacpc
      real*8 dragp
      real*8 drage
      real*8 dragpc
      real*8 dsigp
      real*8 dsige
      real*8 dsigpc
      real*8 roe
      real*8 dampc
      real*8 dql
      real*8 hprod
      real*8 cprod
      real*8 dhprdp
      real*8 dhprdc
      real*8 dhprde
      real*8 dcprdp
      real*8 dcprdc
      real*8 dcprde
      real*8 eskcd
      real*8 pbounc
      real*8 htc
      real*8 tbound
      real*8 hflux
      real*8 dhflxp
      real*8 dhflxe
      real*8 dhflxc
      real*8 cprd
      real*8 edif
      real*8 denrd
      real*8 vfd
      real*8 rl
      real*8 dvfp
      real*8 xtol
      real*8 flow_tol
c      real*8 fimped 
      real*8 p_energy
      real*8 cden_correction, cden_cor
      real*8 spec1, fracc,ms,xf,af,bf,cf,df,ef,dumf,tltemp
      real*8 pv_tol, satdif
c gaz debug 121714
      real*8 dum_gaz
      real*8 xair, pflowa_tol
      real*8 sk_hum, dqt_hum, sk_air, dq_air, sk_h2o
      real*8 qng_old,eh_air,qng_mod
      real*8 diff_w,diff_e,diff_a,permd_air_mult, permd_hum_mult
      real*8 dsk_hump,dsk_humpc,dsk_humt,psatld,pdiff,huma_fixed
c gaz 062916      
      real*8 psatl_100, dpsatt_100, dhumidp, dhumidpc, dhumidt
      real*8 dsk_hums, dsk_airp, dsk_h2op, deh_airp, deh_airt, deh_airpc
      real*8 dsk_airt, dsk_airpc, dsk_h2ot, dsk_h2opc
      real*8  sk_humf, dsk_humfp, dsk_humfpc, dsk_humft
      real*8  sk_airhf, dsk_airhfp, huma_tol
      real*8 pv_hum, t_hum, p_hum   
c gaz  081317,082917  
      real*8 ur, dur_dt
      integer kang, ipv_tol
c gaz 110715
      real*8 dum1,dumb,dumc,value(9)
      integer istate, ifail        
      parameter (kang = 1, pflowa_tol= 1.d-12, huma_tol = 1.d-12)
      parameter (permd_air_mult = 1.d-2, permd_hum_mult = 1.d-3)
c gaz 010519      
       real*8 pcrit_h2o, tcrit_h2o
       parameter(pcrit_h2o=22.00d0, tcrit_h2o=373.95)      
c gaz  081917  
      real*8 dcp_dt, dcprt_dum
c gaz 010719
      real*8 dtdp
c gaz 112718
      integer isk_key
      save isk_key
c
c     rol  -  density liquid
c     ros  -  density steam
c     roc  -  density of air
c     rov  -  density vapour
c     enw  -  enthalpy water
c     hcl  -  enthalpy dissolved air
c     ens  -  enthalpy steam
c     hcg  -  enthalpy air gas
c     env  -  enthalpy vapour
c     xnl  -  mass fraction of air in liquid
c     xnv  -  mass fraction of air in vapour
c     visl -  viscosity of liquid
c     visv -  viscosity of vapour
c     rl   -  relative permeability of liquid phase
c          else if(ieosd.eq.3) then
c calculate water vapor pressure at 100 % humidity
c needed for rel perm calc
c     tl -  temperature
c     sl   -  saturation liquid
c
      parameter(xtol=1.d-16)
      parameter(flow_tol=1.d-31)
      parameter(pv_tol=1.e-9)
c note flow_tol=1.d-31 because of default dqpc=1.d-30
c and default dqpc should look like 0.0 but specified flow
c
      integer i_mem_rlp
      save i_mem_rlp
     
      real*8, allocatable :: s0(:)
      real*8, allocatable :: pcp0(:)
      real*8, allocatable :: dpcps0(:)
      real*8, allocatable :: rlf0(:)
      real*8, allocatable :: drlfs0(:)
      real*8, allocatable :: rvf0(:)
      real*8, allocatable :: drvfs0(:)
c gaz 112818 
      dcp_dt = 0.0

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
     &        rvf,drvef,s,pcp0,dpcps0,rlf0,drlfs0,
     &        rvf0,drvfs0,s0)
      else if(iad.gt.abs(iexrlp)) then
         call pcp_save(2,neq,ndummy,0,pcp,dpcef,rlf,drlef,
     &        rvf,drvef,s,pcp0,dpcps0,rlf0,drlfs0,
     &        rvf0,drvfs0,s0)
         do mid=1,neq
            mi=mid+ndummy
            drlpf(mi)=0.0          
            drvpf(mi)=0.0          
         enddo
      endif
      iieosl=0
      dtin=1.0/dtot
c
c generate relative permeabilities
c
      if (.not. rlpnew) call rlperm(ndummy,1)
c
c call capillary pressure routine
c
c      call cappr(1,ndummy)
      if (.not. rlpnew) call cappr(1,ndummy)
c
      if(rlpnew) call rlp_cap(ndummy)
c
c  call variable rock properties if appropriate  
c  rock state started at last time step temperatures     
      if(iad.eq.0) call vrock_ctr(3,0)
      call vrock_ctr(1,0)
      call vrock_ctr(2,0)
c gaz 111118
c allocate temp storage for sources/sinks
c      
      if(.not.allocated(sk_temp)) allocate(sk_temp(neq))   
      
       do 100 mid=1,neq
         mi=mid+ndummy
         if(igrav.ne.0) then
          p_energy = -grav*cord(mi,igrav)
         else
          p_energy = 0.0d0
         endif
c gaz 081317         
         if(ivrock.ne.0) then
          ur = urock(mi)   
          dur_dt = durockt(mi)
         else
          urock(mi) = denr(mi)*cpr(mi)*t(mi)
          durockt(mi) = denr(mi)*cpr(mi)    
          ur = urock(mi)   
          dur_dt  = durockt(mi)
         endif
         ieosd=ieos(mi)
         iieosd=iieos(mi)
         if(ieosd.ge.2) ifree1 = ifree1 + 1
c
c undo equivalence relations for relative perms
c
         if (rlp_flag .eq. 1) then
            xrl=rlf(mi)
            drl=drlef(mi)
            drlp=drlpf(mi)
            xrv=rvf(mi)
            drv=drvef(mi)
            drvp=drvpf(mi)
         else
            if (ieosd .eq. 1) then
               xrl = 1.0
               xrv = 0.0
            else if (ieosd.eq.3) then
               xrl = 0.0
               xrv = 1.0
            else if (ieosd.eq.4) then
               xrl = 0.0
               xrv = 0.0
            end if
            drl = 0.0
            drlp = 0.0
            drv = 0.0
            drvp = 0.0
         end if
c
c zero flow terms
c
         qh(mi)=0.0
         qc(mi)=0.0
         dq(mi)=0.0
         dqt(mi)=0.0
         dqh(mi)=0.0
         deqh(mi)=0.0
         dcqh(mi)=0.0
         dqc(mi)=0.0
         deqc(mi)=0.0
         dcqc(mi)=0.0
c gaz debug 082115
         dqpc(mi) = 0.0
c
c adjust coefficients for thermo fits
c
         if(iieosd.ge.0) then
c iieos(mi) refers to the thermo set
            if(iieosd.ne.iieosl) then
c
c  liquid phase coefficients
c
c liquid enthalpy
c numerator coefficients
      ela0=cel(1,iieosd)
      elpa1=cel(2,iieosd)
      elpa2=cel(3,iieosd)
      elpa3=cel(4,iieosd)
      elta1=cel(5,iieosd)
      elta2=cel(6,iieosd)
      elta3=cel(7,iieosd)
      elpta=cel(8,iieosd)
      elp2ta=cel(9,iieosd)
      elpt2a=cel(10,iieosd)
c denomenator coefficients
      elb0=cel(11,iieosd)
      elpb1=cel(12,iieosd)
      elpb2=cel(13,iieosd)
      elpb3=cel(14,iieosd)
      eltb1=cel(15,iieosd)
      eltb2=cel(16,iieosd)
      eltb3=cel(17,iieosd)
      elptb=cel(18,iieosd)
      elp2tb=cel(19,iieosd)
      elpt2b=cel(20,iieosd)
c
c liquid density
c numerator coefficients
      dla0=crl(1,iieosd)
      dlpa1=crl(2,iieosd)
      dlpa2=crl(3,iieosd)
      dlpa3=crl(4,iieosd)
      dlta1=crl(5,iieosd)
      dlta2=crl(6,iieosd)
      dlta3=crl(7,iieosd)
      dlpta=crl(8,iieosd)
      dlp2ta=crl(9,iieosd)
      dlpt2a=crl(10,iieosd)
c denomenator coefficients
      dlb0=crl(11,iieosd)
      dlpb1=crl(12,iieosd)
      dlpb2=crl(13,iieosd)
      dlpb3=crl(14,iieosd)
      dltb1=crl(15,iieosd)
      dltb2=crl(16,iieosd)
      dltb3=crl(17,iieosd)
      dlptb=crl(18,iieosd)
      dlp2tb=crl(19,iieosd)
      dlpt2b=crl(20,iieosd)
c
c
c liquid viscosity
c numerator coefficients
      vla0=cvl(1,iieosd)
      vlpa1=cvl(2,iieosd)
      vlpa2=cvl(3,iieosd)
      vlpa3=cvl(4,iieosd)
      vlta1=cvl(5,iieosd)
      vlta2=cvl(6,iieosd)
      vlta3=cvl(7,iieosd)
      vlpta=cvl(8,iieosd)
      vlp2ta=cvl(9,iieosd)
      vlpt2a=cvl(10,iieosd)
c denomenator coefficients
      vlb0=cvl(11,iieosd)
      vlpb1=cvl(12,iieosd)
      vlpb2=cvl(13,iieosd)
      vlpb3=cvl(14,iieosd)
      vltb1=cvl(15,iieosd)
      vltb2=cvl(16,iieosd)
      vltb3=cvl(17,iieosd)
      vlptb=cvl(18,iieosd)
      vlp2tb=cvl(19,iieosd)
      vlpt2b=cvl(20,iieosd)
c
c
c  vapor phase coefficients
c
c vapor enthalpy
c numerator coefficients
      eva0=cev(1,iieosd)
      evpa1=cev(2,iieosd)
      evpa2=cev(3,iieosd)
      evpa3=cev(4,iieosd)
      evta1=cev(5,iieosd)
      evta2=cev(6,iieosd)
      evta3=cev(7,iieosd)
      evpta=cev(8,iieosd)
      evp2ta=cev(9,iieosd)
      evpt2a=cev(10,iieosd)
c denomenator coefficients
      evb0=cev(11,iieosd)
      evpb1=cev(12,iieosd)
      evpb2=cev(13,iieosd)
      evpb3=cev(14,iieosd)
      evtb1=cev(15,iieosd)
      evtb2=cev(16,iieosd)
      evtb3=cev(17,iieosd)
      evptb=cev(18,iieosd)
      evp2tb=cev(19,iieosd)
      evpt2b=cev(20,iieosd)
c
c vapor density
c numerator coefficients
      dva0=crv(1,iieosd)
      dvpa1=crv(2,iieosd)
      dvpa2=crv(3,iieosd)
      dvpa3=crv(4,iieosd)
      dvta1=crv(5,iieosd)
      dvta2=crv(6,iieosd)
      dvta3=crv(7,iieosd)
      dvpta=crv(8,iieosd)
      dvp2ta=crv(9,iieosd)
      dvpt2a=crv(10,iieosd)
c denomenator coefficients
      dvb0=crv(11,iieosd)
      dvpb1=crv(12,iieosd)
      dvpb2=crv(13,iieosd)
      dvpb3=crv(14,iieosd)
      dvtb1=crv(15,iieosd)
      dvtb2=crv(16,iieosd)
      dvtb3=crv(17,iieosd)
      dvptb=crv(18,iieosd)
      dvp2tb=crv(19,iieosd)
      dvpt2b=crv(20,iieosd)
c
c vapor viscosity
c numerator coefficients
c
      vva0=cvv(1,iieosd)
      vvpa1=cvv(2,iieosd)
      vvpa2=cvv(3,iieosd)
      vvpa3=cvv(4,iieosd)
      vvta1=cvv(5,iieosd)
      vvta2=cvv(6,iieosd)
      vvta3=cvv(7,iieosd)
      vvpta=cvv(8,iieosd)
      vvp2ta=cvv(9,iieosd)
      vvpt2a=cvv(10,iieosd)
c denomenator coefficients
      vvb0=cvv(11,iieosd)
      vvpb1=cvv(12,iieosd)
      vvpb2=cvv(13,iieosd)
      vvpb3=cvv(14,iieosd)
      vvtb1=cvv(15,iieosd)
      vvtb2=cvv(16,iieosd)
      vvtb3=cvv(17,iieosd)
      vvptb=cvv(18,iieosd)
      vvp2tb=cvv(19,iieosd)
      vvpt2b=cvv(20,iieosd)
      iieosl=iieosd
      endif
c evaluate variables
      tl=t(mi)
      pcl=pci(mi)
      pl=phi(mi)
c
c evaluate thermo functions and derivatives
c
c use total pressure in compressed liquid equations
c GAZ 5/8/97
c     if(ieosd.ne.3) then
        x=pl
        x2=x*x
        x3=x2*x
c     endif
c
c two phase conditions
c
      dpsatt=0.0
      dpct = 0.0
      psatl_100 = 0.0
      if(ieosd.eq.2) then
c        if(isalt.ne.0) then
c DRH 12/03/12
c gaz 070813 call added saltctr
c            call saltctr(1,mi,dpsatt,dpsats)
c            pcl = pci(mi)
c            pv=pl-pcl
c            dpct=-dpsatt
c         else
            pcl=pl-psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,
     &                   0,an(mi))
            pv=pl-pcl
            pci(mi)=pcl
            dpct=-dpsatt
            dtdp = 1./dpsatt
      else if(ieosd.eq.3) then
c calculate water vapor pressure at 100 % humidity
c needed for rel perm calc
          if(tl.ge.tcrit_h2o) then
              psatl_100 = pcrit_h2o
          else
            psatl_100 = psatl(tl,pcp(mi),dpcef(mi),dpsatt_100,dpsats,
     &                   0,an(mi))
          endif
      endif
c 
c calculate vapor pressure
c
c GAZ 5/8/97
c     if(ieosd.ne.1) then
      pv=pl-pcl
      if(pv.lt.pv_tol) then
       ipv_tol = 0
       xv= pv_tol
      else
       xv= pv
       ipv_tol = 0
      endif
      if(ieosd.eq.3) then          
       humida(mi) = xv/max(psatl_100,1.d-15)
       dhumidp = 1./max(psatl_100,1.d-15)
       dhumidpc = -dhumidp
       dhumidt = (-xv/max(psatl_100,1.d-15)**2)*dpsatt_100
      else 
       humida(mi) = 1.  
       dhumidp = 0.0
       dhumidpc = 0.0
       dhumidt = 0.0       
      endif    
      xv2=xv*xv
      xv3=xv2*xv
c     endif
      dtps=0.0
c
c calculate liquid terms
c
c GAZ 5/8/97
c     if(ieosd.ne.3) then
      tl2=tl*tl
      tl3=tl2*tl
      tlx=x*tl
      tl2x=tl2*x
      tlx2=tl*x2
c     endif
c
c calculate vapor terms
c
c GAZ 5/8/97
c     if(ieosd.ne.1) then
      tl2=tl*tl
      tl3=tl2*tl
      tlxv=xv*tl
      tl2xv=tl2*xv
      tlxv2=tl*xv2
c     endif
c GAZ 5/8/97
c     if(ieosd.ne.3) then
c
c liquid enthalpy
c
            if(iwater_table.ne.1) then
c*****
c     liquid enthalpy
      enwn1=ela0+elpa1*x+elpa2*x2+elpa3*x3
      enwn2=elta1*tl+elta2*tl2+elta3*tl3
      enwn3=elpta*tlx+elpt2a*tl2x+elp2ta*tlx2
      enwn=enwn1+enwn2+enwn3
      enwd1=elb0+elpb1*x+elpb2*x2+elpb3*x3
      enwd2=eltb1*tl+eltb2*tl2+eltb3*tl3
      enwd3=elptb*tlx+elpt2b*tl2x+elp2tb*tlx2
      enwd=enwd1+enwd2+enwd3
      enw=enwn/enwd
      enl=enw + p_energy
c
c derivatives of enthalpy
c
      dhwpn1=elpa1+2*elpa2*x+3*elpa3*x2+elpta*tl
      dhwpn1=enwd*(dhwpn1+elpt2a*tl2+elp2ta*2*tlx)
      dhwpn2=elpb1+2*elpb2*x+3*elpb3*x2+elptb*tl
      dhwpn2=enwn*(dhwpn2+elpt2b*tl2+elp2tb*2*tlx)
      dhwpn=dhwpn1-dhwpn2
      dhwpd=enwd**2
      dhwp=dhwpn/dhwpd
      dhwtn1=elta1+2*elta2*tl+3*elta3*tl2+elpta*x
      dhwtn1=enwd*(dhwtn1+elpt2a*2*tlx+elp2ta*x2)
      dhwtn2=eltb1+2*eltb2*tl+3*eltb3*tl2+elptb*x
      dhwtn2=enwn*(dhwtn2+elpt2b*2*tlx+elp2tb*x2)
      dhwtn=dhwtn1-dhwtn2
      dhwtd=enwd**2
      dhwt=dhwtn/dhwtd
      dhlt=dhwt
      dhlp=dhwp
      dhlpc=0.0
c
c liquid density
c
      rnwn1=dla0+dlpa1*x+dlpa2*x2+dlpa3*x3
      rnwn2=dlta1*tl+dlta2*tl2+dlta3*tl3
      rnwn3=dlpta*tlx+dlpt2a*tl2x+dlp2ta*tlx2
      rnwn=rnwn1+rnwn2+rnwn3
      rnwd1=dlb0+dlpb1*x+dlpb2*x2+dlpb3*x3
      rnwd2=dltb1*tl+dltb2*tl2+dltb3*tl3
      rnwd3=dlptb*tlx+dlpt2b*tl2x+dlp2tb*tlx2
      rnwd=rnwd1+rnwd2+rnwd3
      rnw=rnwn/rnwd
      rol=rnw
      if(cden) then
c     Add correction for liquid species
         cden_cor = cden_correction(mi)
         rol = rol + cden_cor
      end if

c
c       derivatives of density
c
      drlpn1=dlpa1+2*dlpa2*x+3*dlpa3*x2+dlpta*tl
      drlpn1=rnwd*(drlpn1+2*dlp2ta*tlx+dlpt2a*tl2)
      drlpn2=dlpb1+2*dlpb2*x+3*dlpb3*x2+dlptb*tl
      drlpn2=rnwn*(drlpn2+2*dlp2tb*tlx+dlpt2b*tl2)
      drlpn=drlpn1-drlpn2
      drolpd=rnwd**2
      drolp=drlpn/drolpd
      drlen1=dlta1+2*dlta2*tl+3*dlta3*tl2+dlpta*x
      drlen1=rnwd*(drlen1+dlp2ta*x2+2*dlpt2a*tlx)
      drlen2=dltb1+2*dltb2*tl+3*dltb3*tl2+dlptb*x
      drlen2=rnwn*(drlen2+dlp2tb*x2+2*dlpt2b*tlx)
      drlen=drlen1-drlen2
      droled=rnwd**2
      drolt=drlen/droled
      drolpc=0.0
c
c liquid viscosity
c
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
c
c derivatives of liquid viscosity
c
      dvlpn1=vlpa1+2*vlpa2*x+3*vlpa3*x2+vlpta*tl
      dvlpn1=vild*(dvlpn1+2*vlp2ta*tlx+vlpt2a*tl2)
      dvlpn2=vlpb1+2*vlpb2*x+3*vlpb3*x2+vlptb*tl
      dvlpn2=viln*(dvlpn2+2*vlp2tb*tlx+vlpt2b*tl2)
      dvlpn=dvlpn1-dvlpn2
      dvilpd=vild**2
      dvislp=dvlpn/dvilpd
      dvlen1=vlta1+2*vlta2*tl+3*vlta3*tl2+vlpta*x
      dvlen1=vild*(dvlen1+vlp2ta*x2+2*vlpt2a*tlx)
      dvlen2=vltb1+2*vltb2*tl+3*vltb3*tl2+vlptb*x
      dvlen2=viln*(dvlen2+vlp2tb*x2+2*vlpt2b*tlx)
      dvlen=dvlen1-dvlen2
      dviled=vild**2
      dvislt=dvlen/dviled
      dvlpc=0.0
      
c gaz 110715
c new supercritical table similiar to doherty fast table as modified by rajesh pawar
            else if(iwater_table.ne.0) then
c subroutine h2o_properties_new(iflg,iphase,var1,var2,var3,istate,var4,var5,var6)  
c    real*8 var1,var2,var3,var4,var5(9),var6
              if(ieosd.eq.1.or.ieosd.eq.4)then
                  call h2o_properties_new(4,ieosd,pl,tl,dum1,istate,
     &                 dumb,value,dumc)

              elseif (ieosd.eq.2) then                  
c                  call h2o_properties_new(5,ieosd,pl,tl,dum1,istate,
c     &                 dumb,value,dumc)
c gaz 120618 call with ieosd = 1 and use chain rule 
c might need to do something else                  
                  call h2o_properties_new(4,1,pl,tl,dum1,istate,
     &                 dumb,value,dumc)
             endif    
                   rol = value(1)
                   drolt = value(2)
                   drolp = value(3)
                   enl = value(4) + p_energy
c gaz 120518 enw has slight (grav) correction                   
                   enw = enl
                   dhlt = value(5)
                   dhlp = value(6)
                   xvisl = value(7)
                   dvislt = value(8)
                   dvislp = value(9)  
c gaz 120518                   
                   dhlpc=0.0
                   drolpc=0.0 
                   dvlpc=0.0
c chain rule not needed for table (d/dp already includes)
c gaz 120618 need d/dp now                   
c                 if(ieosd.eq.2) then
c                 drolt = 0.0
c                  dhlt = 0.0
c                  dvislt = 0.0  
c                 endif
c gaz 110117
c might need the following un commented
c     if(ieosd.ne.3) then
c      drolpc=0.
c      dhlpc=0.
c      dvlpc=0.
c     endif
c-----------------------------------------------------------------------
C*****
         endif      
c
c set derivatives wrt pc equal to 0 in compressed liquid
c
c GAZ 5/8/97
c     if(ieosd.ne.3) then
c     drolpc=0.
c     dhlpc=0.
c     dvlpc=0.
c     endif
c
c       calculate mass fraction of air in liquid
c
c gaz 010519
        call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
        go to 130
        alpha=1.6111e-04
        dalpca=0.0
        xnl=alpha*pcl + xtol
c
c       derivatives of liquid mass fraction
c
         dxnlp=0.0
         dxnlpc=alpha
         dxnlt=0.0
130     continue         
c           if(xnl.le.0.0) then
c              xnl=0.0
c next line is new 4/13/98
c              dxnlpc=0.0   
c           endif
c
c       enthalpies of mixture
c
c enthalpy of air
c
         hsol=0.0
c
c  Added routine to compute air heat capacity and derivative, old 
c  correlation was for water! (BAR - 8-15-94)
c
      call air_cp(tl, cpa, dcpat)
        hcg=1.e-6*cpa*tl
        hcl=hcg+hsol
        enl=xnl*hcl+(1.d0-xnl)*enw
c
c derivatives of enthalpy of air
c
      dhsolt=0.0
      dhsolp=0.0
c
c  Use new correlation for air heat capacity (BAR - 8-15-94)
c
      dhcgt = 1.e-6 * (cpa + dcpat * tl)
        dhcgp=0.0
        dhlp=dxnlp*hcl+xnl*(dhsolp+dhcgp)+(1.d0-xnl)*dhlp-dxnlp*enw
        dhlt=dxnlt*hcl+xnl*(dhsolt+dhcgt)+(1.d0-xnl)*dhlt-dxnlt*enw
c
        dhslpc=-dhsolp
        dhcgpc=-dhcgp
        dhlpc=dxnlpc*hcl+xnl*(dhslpc+dhcgpc)+(1.d0-xnl)*dhlpc-dxnlpc*enw
c     endif
c
c GAZ 5/8/97
c     if(ieosd.ne.1) then
c
c water vapor enthalpy
c
          if(iwater_table.eq.0) then
c     vapor enthalpy
      ensn1=eva0+evpa1*xv+evpa2*xv2+evpa3*xv3
      ensn2=evta1*tl+evta2*tl2+evta3*tl3
      ensn3=evpta*tlxv+evpt2a*tl2xv+evp2ta*tlxv2
      ensn=ensn1+ensn2+ensn3
      ensd1=evb0+evpb1*xv+evpb2*xv2+evpb3*xv3
      ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
      ensd3=evptb*tlxv+evpt2b*tl2xv+evp2tb*tlxv2
      ensd=ensd1+ensd2+ensd3
      ens=ensn/ensd
      env=ens + p_energy
c
c       derivatives of water vapor enthalpy
c
      dhvp1=evpa1+2*evpa2*xv+3*evpa3*xv2+evpta*tl
      dhvp1=ensd*(dhvp1+2*evp2ta*tlxv+evpt2a*tl2)
      dhvp2=evpb1+2*evpb2*xv+3*evpb3*xv2+evptb*tl
      dhvp2=ensn*(dhvp2+2*evp2tb*tlxv+evpt2b*tl2)
      dhvpn=dhvp1-dhvp2
      dhvpd=ensd**2
      dhvp=dhvpn/dhvpd
      dhvt1=evta1+2*evta2*tl+3*evta3*tl2+evpta*xv
      dhvt1=ensd*(dhvt1+evp2ta*xv2+2*evpt2a*tlxv)
      dhvt2=evtb1+2*evtb2*tl+3*evtb3*tl2+evptb*xv
      dhvt2=ensn*(dhvt2+evp2tb*xv2+2*evpt2b*tlxv)
      dhvtn=dhvt1-dhvt2
      dhvtd=ensd**2
      dhvt=dhvtn/dhvtd
      dhvpc=-dhvp
c
c water vapor density
c
      rnsn1=dva0+dvpa1*xv+dvpa2*xv2+dvpa3*xv3
      rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
      rnsn3=dvpta*tlxv+dvpt2a*tl2xv+dvp2ta*tlxv2
      rnsn=rnsn1+rnsn2+rnsn3
      rnsd1=dvb0+dvpb1*xv+dvpb2*xv2+dvpb3*xv3
      rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
      rnsd3=dvptb*tlxv+dvpt2b*tl2xv+dvp2tb*tlxv2
      rnsd=rnsd1+rnsd2+rnsd3
      rns=rnsn/rnsd
      ros=rns
      rov=ros
c
c       water derivatives of vapor density
c
      drspn1=dvpa1+2*dvpa2*xv+3*dvpa3*xv2+dvpta*tl
      drspn1=rnsd*(drspn1+2*dvp2ta*tlxv+dvpt2a*tl2)
      drspn2=dvpb1+2*dvpb2*xv+3*dvpb3*xv2+dvptb*tl
      drspn2=rnsn*(drspn2+2*dvp2tb*tlxv+dvpt2b*tl2)
      drspn=drspn1-drspn2
      drospd=rnsd**2
      drovp=drspn/drospd
      drsen1=dvta1+2*dvta2*tl+3*dvta3*tl2+dvpta*xv
      drsen1=rnsd*(drsen1+dvp2ta*xv2+2*dvpt2a*tlxv)
      drsen2=dvtb1+2*dvtb2*tl+3*dvtb3*tl2+dvptb*xv
      drsen2=rnsn*(drsen2+dvp2tb*xv2+2*dvpt2b*tlxv)
      drsen=drsen1-drsen2
      drostd=rnsd**2
      drovt=drsen/drostd
      drovpc=-drovp
c
c water vapor viscosity
c
      visn1=vva0+vvpa1*xv+vvpa2*xv2+vvpa3*xv3
      visn2=vvta1*tl+vvta2*tl2+vvta3*tl3
      visn3=vvpta*tlxv+vvpt2a*tl2xv+vvp2ta*tlxv2
      visn=visn1+visn2+visn3
      visd1=vvb0+vvpb1*xv+vvpb2*xv2+vvpb3*xv3
      visd2=vvtb1*tl+vvtb2*tl2+vvtb3*tl3
      visd3=vvptb*tlxv+vvpt2b*tl2xv+vvp2tb*tlxv2
      visd=visd1+visd2+visd3
      vis=visn/visd
      xvisv=vis
c
c derivatives of water vapor viscosity
c
      dvspn1=vvpa1+2*vvpa2*xv+3*vvpa3*xv2+vvpta*tl
      dvspn1=visd*(dvspn1+2*vvp2ta*tlxv+vvpt2a*tl2)
      dvspn2=vvpb1+2*vvpb2*xv+3*vvpb3*xv2+vvptb*tl
      dvspn2=visn*(dvspn2+2*vvp2tb*tlxv+vvpt2b*tl2)
      dvspn=dvspn1-dvspn2
      dvispd=visd**2
      dvisvp=dvspn/dvispd
      dvsen1=vvta1+2*vvta2*tl+3*vvta3*tl2+vvpta*xv
      dvsen1=visd*(dvsen1+vvp2ta*xv2+2*vvpt2a*tlxv)
      dvsen2=vvtb1+2*vvtb2*tl+3*vvtb3*tl2+vvptb*xv
      dvsen2=visn*(dvsen2+vvp2tb*xv2+2*vvpt2b*tlxv)
      dvsen=dvsen1-dvsen2
      dvised=visd**2
      dvisvt=dvsen/dvised
      dvvpc=-dvisvp
           else if(iwater_table.ne.0) then
c subroutine h2o_properties_new(iflg,iphase,var1,var2,var3,istate,var4,var5,var6)  
c    real*8 var1,var2,var3,var4,var5(9),var6
                if(ieosd.eq.3)then
                  call h2o_properties_new(4,ieosd,pv,tl,dum1,istate,
     &                 dumb,value,dumc)

                elseif (ieosd.eq.2) then
c gaz 120718  derivatives match rat polynomials                
c                  call h2o_properties_new(6,ieosd,pl,tl,dum1,istate,
c     &                 dumb,value,dumc)
                  call h2o_properties_new(4,3,pv,tl,dum1,istate,
     &                 dumb,value,dumc)
               endif    
                   rov = value(1)
c gaz  120418                  
                   ros = rov
                   drovt = value(2)
                   drovp = value(3)
                   env = value(4) + p_energy
c gaz  120418                   
                   ens = env
                   dhvt = value(5)
                   dhvp = value(6)
                   xvisv = value(7)
                   vis = xvisv
                   dvisvt = value(8)
                   dvisvp = value(9)  
                   dhvpc=-dhvp
                   drovpc=-drovp
                   dvvpc=-dvisvp
                   
c chain rule not needed for table (d/dp already includes)
c gaz 120718 need to have d/dp                   
c               if(ieosd.eq.2) then
c                drovt = 0.0
c                dhvt = 0.0
c                dvisvt = 0.0  
c               endif
          endif           
c gaz 111915
c also changed pv_tol to 1.e-3
c
      if(ipv_tol.eq.1) then
       dhvpc = 0.0
       dhvp = 0.0
       drovpc = 0.0
       drovp = 0.0
       dvvpc = 0.0
       dvisvp = 0.0
      endif
c
c density of air
c
      pcl0=0.101325
      roc0=1.292864
**** density of air eqn (31)
      drocpc=roc0*(273./(tl+273.))/pcl0
      roc=drocpc*pcl
      rov = roc + ros
      xnv = roc / rov + xtol
**** molar densities used for avg. mol wt. calculation
      rocmol = roc/mw_air
      rosmol = ros/mw_water
      rovmol = rocmol + rosmol
**** mole fraction calculated ****
      xnvmol = rocmol/rovmol
**** mass fraction calculated ****
      avgmolwt(mi) = xnvmol*mw_air+(1-xnvmol)*mw_water
c
c derivatives of the density of air and mixture
c
      droct=-roc/(tl+273.)
      drocp=0.0
      drovp=drovp+drocp
      drovt=drovt+droct
      drovpc=drovpc+drocpc
c
c vapor mass fraction
c
c
c       derivatives of mass fraction
c
        dxnvt=(rov*droct-roc*drovt)/(rov*rov)
        dxnvp=-roc*drovp/(rov*rov)
        dxnvpc=(rov*drocpc-roc*drovpc)/(rov*rov)
c
c       enthalpy of air and mixture
c
c
c  Added routine to compute air heat capacity and derivative, old 
c  correlation was for water! (BAR - 8-15-94)
c
        call air_cp(tl, cpa, dcpat)
        hcg=1.e-6*cpa*tl
        env=xnv*hcg+ens*(1.d0-xnv)
c
c derivatives of enthalpy of air
c      
        dhvsp=dhvp
        dhvst=dhvt
        dhvspc=dhvpc
c
        dhcgpc=0.0
c
c  Use new correlation for air heat capacity (BAR - 8-15-94)
c
        dhcgt = 1.e-6 * (cpa + dcpat * tl)
c        dhcgt=(-11.76+0.0354*(tl+273.))*1.e-6*tl+cpa
        dhcgp=0.0
c
        dhvt=dxnvt*hcg+xnv*dhcgt+(1.d0-xnv)*dhvt-dxnvt*ens
        dhvp=dxnvp*hcg+(1.d0-xnv)*dhvp-dxnvp*ens
        dhvpc=dxnvpc*hcg+(1.d0-xnv)*dhvpc-dxnvpc*ens
c
c vapor viscosity
c
      xvisc=182.7e-07+0.41e-07*(tl-18.0)
      dxvsct=0.41e-07
      xvisv=(1.00-xnv)*vis+xnv*xvisc
c
c dervatives of vapor viscosity
c
      dvisvt=(1.d0-xnv)*dvisvt-dxnvt*(vis-xvisc)+xnv*dxvsct
      dvisvp=(1.d0-xnv)*dvisvp-dxnvp*(vis-xvisc)
      dvvpc=(1.d0-xnv)*dvvpc-dxnvpc*(vis-xvisc)
c
c GAZ 5/8/97
c     else
c        avgmolwt(mi) = mw_air   PHS 4/2/04  removed this line
c GAZ 5/8/97
c     endif
c
c modify derivatives for 2-phase
c
      if(ieosd.eq.2) then
c gaz 010519    dtpcs important for t as variable when ieosd = 2      
      dtpcs=1.0
c gaz 10-18-2001 set condition to unattainable value
      if(ice .eq. -99) then
         if(ices(mi).eq.1) then
c
c modify values and derivatives when non gas phase is partially frozen
c
c liquid values
c
            siid=sii(mi)
            siie=1.0-siid
            roll=rol
            roli=950.0
            rol=siid*roli+siie*roll
            drolp=siie*drolp
            drolpc=siie*drolpc
            drolt=(roli-roll)
            enll=enl
            enl=siie*enll
            dhlp=siie*dhlp
            dhlt=-enll
            dhlpc=siie*dhlpc
            visll=xvisl
            visli=1.e20
            xvisl=siie*visll+siid*visli
            dvislp=siie*dvislp
            dvlpc=siie*dvlpc
            dvislt=visli-visll
            xnll=xnl
            xnl=siie*xnl
            dxnlp=siie*dxnlp
            dxnlpc=siie*dxnlpc
            dxnlt=-xnl
c vapor terms
c zero out temperature dependence
            drovt=0.0
            dhvt=0.0
            dvisvt=0.0
            dxnvt=0.0
            dtpcs=0.0
         endif
         if(ices(mi).eq.2) then
c
c non gas phase is fully frozen
c
            rol=950.
            drolp=0.0
            drolt=0.0
            drolpc=0.0
            cph=4.18e-3
            enl=cph*(tl-tmelt)
            dhlp=0.
            dhlpc=0.0
            dhlt=cph
            xvisl=1.e20
            dvislp=0.0
            dvislt=0.0
            dvlpc=0.0
c vapor terms leave alone
         endif
c
      endif
c
c derivative wrt pc contain wrt t terms (remember ieos(mi)=2)
c drolpc (for example) is now the temperature derivative
c
      drolp=drolp+drolpc
      drovp=drovp+drovpc
      dhlp=dhlp+dhlpc
      dhvp=dhvp+dhvpc
      dvisvp=dvisvp+dvvpc
      dvislp=dvislp+dvlpc
      dxnlp=dxnlp+dxnlpc
      dxnvp=dxnvp+dxnvpc
      dhcgp=dhcgp+dhcgpc
      dhvsp=dhvsp+dhvspc
      drolpc0=drolpc
      drovpc0=drovpc
      dhlpc0=dhlpc
      dhvpc0=dhvpc
      dvlpc0=dvlpc
      dvvpc0=dvvpc
      dxnlpc0=dxnlpc
      dxnvpc0=dxnvpc
      dhcgpc0=dhcgpc
      dhvspc0=dhvspc
c
      drolpc=drolt+drolpc0*dpct
      drovpc=drovt+drovpc0*dpct
      dhlpc=dhlt+dhlpc0*dpct
      dhvpc=dhvt+dhvpc0*dpct
      dvvpc=dvisvt+dvvpc0*dpct
      dvlpc=dvislt+dvlpc0*dpct
      dxnlpc=dxnlt+dxnlpc0*dpct
      dxnvpc=dxnvt+dxnvpc0*dpct
      dhcgpc=dhcgt+dhcgpc0*dpct
      dhvspc=dhvst+dhvspc0*dpct
c
      if( ivapl .ne. 0 ) then
         drolt=drolpc0*dpsats
         drovt=drovpc0*dpsats
         dhlt=dhlpc0*dpsats
         dhvt=dhvpc0*dpsats
         dvisvt=dvvpc0*dpsats
         dvislt=dvlpc0*dpsats
         dxnlt=dxnlpc0*dpsats
         dxnvt=dxnvpc0*dpsats
         dhcgt=dhcgpc0*dpsats
         dhvst=dhvspc0*dpsats
      else
         drolt=0.0
         drovt=0.0
         dhlt=0.0
         dhvt=0.0
         dvisvt=0.0
         dvislt=0.0
         dxnlt=0.0
         dxnvt=0.0
         dhcgt=0.0
         dhvst=0.0
      end if
c
      dtd=0.0
      endif
      sl=s(mi)
      cp=denr(mi)*cpr(mi)
      por=ps(mi)
      vol=volume(mi)
c      
      endif
c      
c     two phase conditions
c
      if(ieosd.eq.2) then
c
c       saturation and relative permeabilities
c
        sv=1.d0-sl
c
c accumulation terms
c water balance
      den=por*(sl*rol*(1.0-xnl)+sv*rov*(1.0-xnv))
c energy balance      
      eqdum=sl*rol*enl+sv*rov*env-pl
      dene=((1.-por)*cp*tl+por*eqdum)
c air balance      
      denc=por*(sl*rol*xnl+sv*rov*xnv)
c
c       production of steam
c
        rag=rol*xvisv/rov/xvisl
        sig=xrv/(xrv+rag*xrl)
c
c       derivatives of accumulation terms
c
c pressure derivatives
        rop=por*(sv*drovp*(1.0-xnv)+sl*drolp*(1.0-xnl))    
        damp =rop*dtin + por*(sv*rov*(-dxnvp)
     &  +sl*rol*(-dxnlp))*dtin     
c 011219     replace dtdp with 0.   
        daep =((1.d0-por)*(dur_dt*0.0)+
     &   por*(sv*drovp*env+sv*rov*dhvp+
     &   sl*drolp*enl+sl*rol*dhlp)-por)*dtin
        dacp =(por*(sv*drovp*xnv+sv*rov*dxnvp+sl*drolp*xnl+sl*rol
     &   *dxnlp))*dtin
c saturation derivatives    
        damh =por*dtin*((rol*(1.0-xnl)-rov*(1.0-xnv))
     &  +(sl*drolt*(1.0-xnl)+sv*drovt*(1.0-xnv))
     &  +(-sl*rol*dxnlt-sv*rov*dxnvt))
        daeh =por*dtin*((rol*enl-rov*env)
     &  +(sl*drolt*enl+sv*drovt*env)
     &  +(sl*rol*dhlt+sv*rov*dhvt))
        dach =por*dtin*((rol*xnl-rov*xnv)
     &  +(sl*drolt*xnl+sv*drovt*xnv)
     &  +(sl*rol*dxnlt+sv*rov*dxnvt))
c temperature derivatives          
        ropc=por*(sv*drovpc*(1.0-xnv)+sl*drolpc*(1.0-xnl))
        dampc =ropc*dtin+por*(sv*rov*(-dxnvpc)+sl*rol*(-dxnlpc))*dtin
        daepc =((1.d0-por)*(dur_dt*dtpcs)+
     &   por*(sv*drovpc*env+sv*rov*dhvpc+
     &   sl*drolpc*enl+sl*rol*dhlpc))*dtin
        dacpc =por*(sv*drovpc*xnv+sv*rov*dxnvpc
     &   +sl*drolpc*xnl+sl*rol*dxnlpc)*dtin
c
      end if
c
c     compressed liquid
c
      if(ieosd.eq.1) then
      dtps=0.0
      dtpcs=0.0
        sl=1.d0
        sv=0.0
        sig=0.d0
        rov=0.d0
        env=0.d0
      drovp=0.0
      drovt=0.0
      dhvp=0.0
      dhvt=0.0
      dvisvp=0.0
      dvisvt=0.0
      drovpc=0.0
      dhvpc=0.0
      dvvpc=0.0
      dxnvpc=0.0
      dtd=1.0
      xvisv=1.d20
c
c       accumulation terms
c
c water balance        
c
        den=rol*por*(1.0-xnl)
        dene=(1.d0-por)*ur+rol*por*enl-por*pl
        denc=rol*por*xnl
c
c       derivatives of accumulation terms
c
        rop=por*drolp
c water balance        
        damp = (rop*(1.0-xnl)-por*rol*dxnlp)*dtin        
        daep =(rop*enl+por*rol*dhlp-por)*dtin
        dacp =(rop*xnl+por*rol*dxnlp)*dtin
        roe=por*drolt
c water balance        
        damh =(roe*(1.0-xnl)-por*rol*dxnlt)*dtin              
        daeh =(por*rol*dhlt+roe*enl+(1.d0-por)*dur_dt)*dtin
        dach =(por*rol*dxnlt+roe*xnl)*dtin
c water balance        
c assumes drolpc=0.0
        dampc = -por*rol*dxnlpc*dtin              
        daepc = por*rol*dhlpc*dtin
        dacpc = por*rol*dxnlpc*dtin
c
      end if
c
c     superheated vapour
c
      if(ieosd.eq.3) then
        dtps=0.0
        dtpcs=0.0
        sl=0.d0
        sv=1.d0
        sig=1.d0
        rol=0.d0
        enl=0.d0
        drolp=0.0
        drolt=0.0
        dhlp=0.0
        dhlt=0.0
        dvislp=0.0
        dvislt=0.0
        drolpc=0.0
        dhlpc=0.0
        dvlpc=0.0
        dxnlpc=0.0
        dtd=1.0
        xvisl=1.d20
c
c       accumulation terms
c
c water balance        
        den=rov*por*(1.0-xnv)
c enthalpy balance      
        dene=(1.d0-por)*ur+rov*por*env-por*pl
c air balance      
        denc=rov*por*xnv
c
c       derivatives of accumulation terms
c
        rop=por*drovp
        damp =rop*dtin
c water balance        
        damp =(rop*(1.0-xnv)-por*rov*dxnvp)*dtin
        daep =(rop*env+por*rov*dhvp-por)*dtin
        dacp =(por*rov*dxnvp+xnv*rop)*dtin
        roe=por*drovt
c enthalpy balance      
        damh =(roe*(1.0-xnv)-por*rov*dxnvt)*dtin
c gaz 090117        
        daeh =((1.d0-por)*dur_dt+por*rov*dhvt+env*roe)*dtin    
        dach =(por*rov*dxnvt+xnv*roe)*dtin
        ropc=por*(drovpc)
c air balance      
        dampc =(ropc*(1.0-xnv)-por*rov*dxnvpc)*dtin        
        daepc =(ropc*env+por*rov*dhvpc)*dtin
        dacpc =(por*rov*dxnvpc+xnv*ropc)*dtin
      end if
c
      if(ieosd.ne.4) then
c
c store derivatives of accumulation terms
c
      dmef(mi)=damh
      dmpf(mi)=damp
      dmc(mi)=dampc
      depf(mi)=daep
      deef(mi)=daeh
      dec(mi)=daepc
      dcp(mi)=dacp
      dce(mi)=dach
      dcc(mi)=dacpc
      dtpa(mi)=dtps
      dtpac(mi)=dtpcs
      dtpae(mi)=dtd
      t(mi)=tl
      s(mi)=sl
      a(mi)=den
      a(mi+neq)=dene
      a(mi+neq+neq)=denc
      dstm(mi)=por*rov*sv*vol*(1.0-xnv)
c
c GAZ 5/8/97
c     if(ieosd.ne.3) then
        dql=rol*xrl/xvisl
        dil(mi)=dql
        dilp(mi)=drolp*xrl/xvisl-rol*xrl/xvisl**2*dvislp
     *  +rol*drlp/xvisl
        dile(mi)=drolt*xrl/xvisl+rol*drl/xvisl
     *  -rol*xrl/xvisl**2*dvislt
        enlf(mi)=enl
        dglp(mi)=drolp
        dgle(mi)=drolt
        delf(mi)=dhlp
        delef(mi)=dhlt
        delcf(mi)=dhlpc
        dclf(mi)=dxnlp
        dclef(mi)=dxnlt
        dclcf(mi)=dxnlpc
        dilc(mi)=(drolpc/xvisl-rol/xvisl**2*dvlpc)*xrl
        dglc(mi)=drolpc
c     endif
c
c GAZ 5/8/97
c     if(ieosd.ne.1) then
        dqv=rov*xrv/xvisv
        div(mi)=dqv
        divp(mi)=drovp*xrv/xvisv-rov*xrv/xvisv**2*dvisvp
     *  +rov*drvp/xvisv
        dive(mi)=drovt*xrv/xvisv+rov*drv/xvisv
     *  -rov*xrv/xvisv**2*dvisvt
        envf(mi)=env
        dgvp(mi)=drovp
        dgve(mi)=drovt
        devf(mi)=dhvp
        devef(mi)=dhvt
        devcf(mi)=dhvpc
        dcvf(mi)=dxnvp
        dcvef(mi)=dxnvt
        dcvcf(mi)=dxnvpc
        divc(mi)=(drovpc/xvisv-rov/xvisv**2*dvvpc)*xrv
        dgvc(mi)=drovpc
c     endif
c
      rolf(mi)=rol
      rovf(mi)=rov
      cnlf(mi)=xnl
      cnvf(mi)=xnv
c
      enva(mi)=ens
      denvap(mi)=dhvsp
      denvae(mi)=dhvst
      denvac(mi)=dhvspc
c
c ********** source term code start *****************
c
c gaz debug terms (so intel debugger recognizes fdum and l from use module)
      dum_gaz=fdum+l 
c
c gaz  110218 fix for not zeroing out qng(mi)     
       if(abs(pflowa(mi)).gt.pflowa_tol) qng(mi) = 0.0
       if(xairfl(mi).gt.0.0) qng(mi) = 0.0
c gaz 111418 can now set qc for specified air flowrate(qng ne 0)
       qc(mi)=qng(mi)
       dqc(mi)=0.0
       deqc(mi) = 0.0
       dcqc(mi) = 0.0
c
c identify water source term
      if(ka(mi).gt.0) then
c specified water source/sink
       qdis=sk(mi)
      else
c calculated water source/sink (via specified pressure or saturation or humidity)
       qdis = 0.0
       qh(mi) = 0.0
       sk(mi) = 0.0
      endif
      pldif = 0.0
      satdif = 0.0
c      
c calculate capillary pressure terms
c no capp pres terms because liquid evaluated at total pressure
c form pressure dependent flow term
c
      kq=ka(mi)
      sk_hum = 0.0  
      dsk_hums = 0.0
      dsk_humpc = 0.0
      dsk_humt = 0.0   
      sk_air = 0.0
      sk_h2o = 0.0
c start loop on specified pressure      
      if(kq.lt.0 .and. compute_flow) then
c pflow >= 0 , specified pressure
c pflow <  0 , specified saturation (sometimes resulting from humidity )
         if(pflow(mi) .lt. 0.) then 
            permsd=wellim(mi)
            if(ieosd.eq.2) then
c gaz debug 081815 only allow outflow
c basically for two phase conditions, saturation is managed (fixed BC) with water flow
c pressure is managed (fixed BC) with airflow
c -888 disables humidity:dqt is d/ds for 2 phase, d/dt for 1 phase
c gaz 112418  try allowing flow in or out             
c              if(sl+pflow(mi).gt.0.0.and.pflow(mi).ne.-888.) then
               satdif = sl+pflow(mi)
               sk_hum = permsd*(satdif)
               dsk_hums = permsd
               dsk_humpc = 0.0
               dsk_humt = 0.0
c gaz 112418 send to different loops               
c               qdis = sk_hum
c               sk(mi) = qdis
c               dqt(mi) = dsk_hums
c              else
c               sk_hum = 0.0
c               dsk_hums = 0.0
c               dsk_humpc = 0.0
c               dsk_humt = 0.0               
c              endif
c gaz debug 081815
c only outflow  
c gaz 112518 do nothing here pflowa addressed later             
c              if(pflowa(mi).gt.pflowa_tol) then
c                if(phi(mi)-pflowa(mi).gt.0.0) then
c                 permsda = permsd*permd_air_mult
c                 sk_air =   permsda*(phi(mi)-pflowa(mi))
c                 dsk_airp=  permsda
c                 qc(mi) =   permsda*(phi(mi)-pflowa(mi))
c                 dqc(mi) =  permsda                 
c                else
c do nothing      
c                 sk_air = 0.0
c                 dsk_airp = 0.0          
c                endif
c              endif
            else if(ieosd.eq.3) then
c do nothing here
            endif
c            qc(mi)=0.
c            dqc(mi)=0.
c dqpc (now qng) should be set in co2ctr and input
         else if(abs(wellim(mi)).ne.0.0) then
c now for positive pflow()             
            pldif=pl-pflow(mi)
            permsd=wellim(mi)
c kq=-2 means inflow not allowed
            if(pldif.le.0.0d00.and.kq.eq.-2) permsd=0.0
c gaz 092418 - can't hold fixed pressure if air partial pressure (pcl) is larger  than pl   
c gaz 110518 - changed definition 
c gaz 111818 - changed definition again            
c if air pressure gets larger than the specified  total pressure; turn off specified
c pressure   
c gaz 112718 changed    pflow to only inflow
             if(ieosd.eq.2) then
              if(phi(mi).gt.pflow(mi).and.qc(mi).lt.0.0) then
               sk(mi)=0.0
               dq(mi)=0.0 
               qdis = 0.0
              else
               qdis=permsd*(pldif)
               sk(mi)=qdis
               dq(mi)=permsd  
              endif
c not 2-phase              
             else
              if(phi(mi).gt.pflow(mi).and.qc(mi).lt.0.0) then
               sk(mi)=0.0
               dq(mi)=0.0 
               qdis = 0.0
              else
               qdis=permsd*(pldif)
               sk(mi)=qdis
               dq(mi)=permsd  
              endif 
            endif
c gaz debug 122814
c check for air fraction of water source sk();(only for inflow)
c if we assume dry air, then we could calculate humidity   
c xairfl splits a water  inflow source (including fixed pressure) into 
c a water and air source (avialable from ngas or boun input)            
           if(xairfl(mi).gt.0.0.and.qdis.lt.0.0) then
            xair = xairfl(mi)
            sk(mi) = (1.0-xair)*qdis
            qc(mi) = xair*qdis
            qng(mi) = qc(mi)
            dqc(mi) = dq(mi)*xair
            dq(mi) = dq(mi)*(1.0-xair)
            qdis = sk(mi)
c gaz debug 070416: could add enthalpy here  with call to flow_humidity_bc(2...
c need to add conditional 2 for enthalpy calc (no d/dp,d/dy etc)            
           endif
         endif
c end loop on specified saturation or total pressure
      else
c
c gaz 111118 add "xfa" capability to specified water flow
c water flow, however is split into water and air flow
c 
       if(xairfl(mi).gt.0.0.and.sk(mi).lt.0.0) then
          if(iad.eq.0) then
c save source info   
            sk_temp(mi) = sk(mi)
          endif
          sk(mi) = sk_temp(mi)
          qdis = sk(mi)
          xair = xairfl(mi)
          sk(mi) = (1.0-xair)*qdis
          qc(mi) = xair*qdis
          qng(mi) = qc(mi)
          dqc(mi) = dq(mi)*xair
          dq(mi) = dq(mi)*(1.0-xair)
          qdis = sk(mi) 
       endif
       continue
      end if
c end loop specified pressure
c
c start section on air inflow with specified humidity (inflow only)
c no derivatives because of inflow only condition
c xnva air mass fraction in air with relative humidity huma()
c enva is the enthalpy of the air (Mw/kg)
c xnva() and entha() have been caluculated in flow_hmuidity_bc
c not perfect because of dissolved air in water is assumed to be zero
c
c can't use huma and xairfl at same node
c note; we assume air inflow for this section to 
c be humidified air
c by evaluating qng we omit the qc(mi) calculated above
        qng_old = qng(mi)
        dq_air = 0.0
c check if pflowa exists and  inflow exists
c gaz 092518 changed  so specified pa uses pci(mi) significant
c need derivatives of pci  when  ieosd = 2    
c with pa the partial pressure pa means dry 
c gaz 111318 need to zero permsda        
        permsda = 0.0
        if(pflowa(mi).gt.pflowa_tol) then
            if(pci(mi)-pflowa(mi).lt.0.0) then
                 permsd=abs(wellim(mi))
                 permsda = permsd*permd_air_mult
                 qc(mi) = permsda*(pci(mi)-pflowa(mi))
                 if(ieosd.eq.2) then
c gaz 111318 pci = phi - pv(t)                     
                  dqc(mi) = permsda
                  dcqc(mi) = 0.0
                  deqc(mi) = permsda*dpct 
                 else
                  dqc(mi) = 0.0
                  dcqc(mi) = permsda
                  deqc(mi) = 0.0
                 endif
             else        
                 permsd=abs(wellim(mi))
                 permsda = permsd*permd_air_mult
c                 sk_air should be non zero only if air is not dry
c                  pflowa represents the partial pressure of dry air                 
c                 sk_air = permsda*(pci(mi)-pflowa(mi))
                 qc(mi) = permsda*(pci(mi)-pflowa(mi))
                 if(ieosd.eq.2) then
                  dqc(mi) = permsda
                  deqc(mi) = permsda*dpct 
                 else
                  dqc(mi) = 0.0
                  dcqc(mi) = permsda
                  deqc(mi) = 0.0
                 endif 
c gaz  092518    might be not needed           
c                 dq_air = permsda
c                 dsk_airp =  0.0                 
c                 qng_old = 0.0
             endif
        endif
       if(iha.ne.0) then
c gaz 111318 allow humidity flow (h20)to be calculated 
c if no humidity or air fraction then dry air assumed             
       if(huma(mi).gt.0.0.and.xairfl(mi).eq.0.0) then
c only for inflow of gas phase             
        if(qc(mi).lt.0.0) then
c water flow in rate sk()
c air flow in is sk_air
c entha() is correct mixture enthalpy
c if qc() is from a specified air pressure then permsda ne 0
c if qc() is from a specified air flowrate then permsda eq 0  
c at this point sk is calculated and if qc <0 , then qc is split  
         if(pflowa(mi).gt.pflowa_tol) then
c specified air pressure             
          if(ieosd.eq.2) then
           dqc(mi) = permsda
           deqc(mi) = permsda*dpct 
          else
           dqc(mi) = 0.0
           dcqc(mi) = permsda
           deqc(mi) = 0.0
          endif
c specified air flowrate (no derivatives)          
         else
           dqc(mi) = 0.0
           dcqc(mi) = 0.0
           dcqh(mi) = 0.0          
         endif
                
         sk_h2o = (1.0-xnva(mi))*qc(mi)
         dsk_h2op = (1.0-xnva(mi))*dqc(mi)
         dsk_h2ot = (1.0-xnva(mi))*deqc(mi)
         dsk_h2opc = (1.0-xnva(mi))*dcqc(mi)
         sk_air = xnva(mi)*qc(mi)
         dsk_airp = xnva(mi)*dqc(mi)
         dsk_airt = xnva(mi)*deqc(mi)
         dsk_airpc = xnva(mi)*dcqc(mi)
         eh_air = entha(mi)*qc(mi)
         deh_airp = entha(mi)*dqc(mi)
         deh_airt = entha(mi)*deqc(mi)
         deh_airpc = entha(mi)*dcqc(mi)
         qng_old = 0.0
c gaz 111318
c new air source and derivatives here
c they will be added later so zero out
c air flow and derivatives                 
         qc(mi) = 0.0
         dqc(mi) = 0.0
         deqc(mi) = 0.0
         dcqc(mi) = 0.0
        endif
      endif
      endif
c
c end section on air inflow with specified humidity
c
c gaz 070516
c start section on specified humidity (inflow only and outflow)
c water flows in or out
c       
c
              t_hum = 0.0
              p_hum = 0.0
              huma_fixed = 0.0
              if(iha.ne.0) then
               if(huma(mi).lt.0) then
                 huma_fixed = abs(huma(mi))                         
                 permsd=abs(wellim(mi))*permd_hum_mult
                 t_hum = thuma(mi)
                 p_hum = phuma(mi) 
                 pv_hum = psatl(t_hum,pcp(mi),dpcef(mi),dpsatt_100,
     &                   dpsats,0,an(mi))
                 pdiff = p_hum-pci(mi)
                 sk_humf = permsd*(pdiff-huma_fixed*pv_hum)  
                 if(ieosd.eq.2) then
                  dsk_humfp = -permsd
                  dsk_humfpc = 0.0
c  dsk_humft uses pci = p - pv(t)  (pci is not a variable for ieosd = 2)              
                  dsk_humft = -permsd*dpct
                 else
                  dsk_humfp = 0.0
                  dsk_humfpc = -permsd
                  dsk_humft = 0.0                   
                 endif
                if(p_hum.gt.pflowa_tol) then
                 permsda = permsd*permd_air_mult
                 sk_airhf = permsda*(phi(mi)-p_hum)
                 dsk_airhfp =  permsda 
                else
                 sk_airhf = 0.0
                 dsk_airhfp = 0.0
                endif                
               endif
              endif   
555           continue 

c
c end section on air inflow with specified humidity
c
      hprod=0.0
      cprod=0.0

c
c       derivatives of sink terms for outflow
c
        if(ieosd.eq.2.and.qdis.gt.0.0) then
          dragp=(rov*xvisl*(drolp*xvisv+rol*dvisvp)-rol*xvisv*(drovp
     &     *xvisl+rov*dvislp))/(xvisl*rov)**2
          drage=(rov*xvisl*(drolt*xvisv+rol*dvisvt)-rol*xvisv*(drovt
     &     *xvisl+rov*dvislt))/(xvisl*rov)**2
          dragpc=(rov*xvisl*(drolpc*xvisv+rol*dvvpc)-rol*xvisv*(drovpc
     &     *xvisl+rov*dvlpc))/(xvisl*rov)**2
c
          if(sig.ne.0.d0.and.sig.ne.1.d0) then
            dsigp=-xrv*dragp*xrl/(xrv+rag*xrl)**2
            dsige=((xrv+rag*xrl)*drv-xrv*(drv+rag*drl))
     &       /(xrv+rag*xrl)**2
     &      -xrv*drage*xrl/(xrv+rag*xrl)**2
            dsigpc=-xrv*dragpc*xrl/(xrv+rag*xrl)**2
          else
            dsigp=0.d0
            dsige=0.d0
            dsigpc=0.d0
          end if
        end if
c
c organize source terms and derivatives
c
c qng(mi) is the specified air flow
c water flows out - vapor phase calculated
         if(qdis.gt.0.0.and.abs(qng_old).lt.flow_tol) then
c flow out of reservoir
          if(ieosd.eq.2) then
c calculate energy flow out 
            hprod=sig*env+(1.0-sig)*enl
            dhprdp=dsigp*env+sig*dhvp-dsigp*enl+(1.0-sig)*dhlp
            dhprdc=dsigpc*env+sig*dhvpc-dsigpc*enl+(1.0-sig)*dhlpc
            dhprde=dsige*env+sig*dhvt-dsige*enl+(1.0-sig)*dhlt 
c            if(eskc(mi).ge.0.0) then
            if(ka(mi).lt.0) then
c this is the humidity condition (specified saturation) and specified pressure
c calculate dry air flow out
              cprod=sig*xnv+(1.0-sig)*xnl
              dcprdp=dsigp*xnv+sig*dxnvp-dsigp*xnl+(1.0-sig)*dxnlp
              dcprdc=dsigpc*xnv+sig*dxnvpc-dsigpc*xnl+(1.0-sig)*dxnlpc
              dcprde=dsige*xnv+sig*dxnvt-dsige*xnl+(1.0-sig)*dxnlt
            else   
c this is for specified water flow only            
              cprod=0.0
              dcprdp=0.0                                           
              dcprdc=0.0
              dcprde=0.0
            endif                    
          else if(ieosd.eq.1) then
c single phase - just water
            hprod=enl
            dhprdp=dhlp
            dhprdc=dhlpc
            dhprde=dhlt
c gaz debug 120814
c            if(eskc(mi).gt.0.0) then
             if(ka(mi).lt.0) then
              cprod=xnl
              dcprdp=dxnlp
              dcprdc=dxnlpc
              dcprde=dxnlt
            else                     
              cprod=0.0
              dcprdp=0.0                                           
              dcprdc=0.0 
              dcprde=0.0
            endif                    
          else if(ieosd.eq.3) then
            hprod=env
            dhprdp=dhvp
            dhprdc=dhvpc
            dhprde=dhvt
c gaz debug 120814
c            if(eskc(mi).ge.0.0) then
            if(ka(mi).lt.0) then
              cprod=xnv
              dcprdp=dxnvp
              dcprdc=dxnvpc
              dcprde=dxnvt
            else                     
              cprod=0.0
              dcprdp=0.0                                           
              dcprdc=0.0
              dcprde=0.0
            endif                    
          endif
          qh(mi)=qdis*hprod       
          dqh(mi)=dhprdp*qdis+hprod*dq(mi)
          deqh(mi)=dhprde*qdis+hprod*dqt(mi)
          dcqh(mi)=dhprdc*qdis 
c gaz debug 120814
c            if(eskc(mi).ge.0.0) then
c gaz debug 121114 (Big change)
c water fraction in outlet stream
c            qdis=permsd*(pldif)
           if(pflow(mi).gt.0.0.and.pldif.gt.0.0) then
c specified pressure
            sk(mi) = (1.0-cprod)*permsd*pldif
            dq(mi) = (1.0-cprod)*permsd - permsd*pldif*dcprdp
            dqt(mi) = - permsd*pldif*dcprde
            dqpc(mi) = - permsd*pldif*dcprdc
            qc(mi)=qc(mi) + cprod*permsd*pldif 
            dqc(mi)= dqc(mi)+cprod*permsd + permsd*pldif*dcprdp 
            deqc(mi)= deqc(mi)+ permsd*pldif*dcprde 
            dcqc(mi)= dcqc(mi)+ permsd*pldif*dcprdc
           else if(pflow(mi).lt.0.0.and.satdif.gt.0.0) then
C (-pflow) is the specified saturation
c specified saturation: qdis = fimped*(sl+pflow(mi))
            satdif = (sl+pflow(mi))
            sk(mi) = (1.0-cprod)*permsd*satdif
            dq(mi) =  - permsd*satdif*dcprdp
            dqt(mi) = (1.0-cprod)*permsd - permsd*satdif*dcprde
            dqpc(mi) = - permsd*satdif*dcprdc 
            if(qc(mi).ge.0.0) then
             qc(mi)=qc(mi) + cprod*permsd*satdif
             dqc(mi)= dqc(mi)+ permsd*satdif*dcprdp 
             deqc(mi)= deqc(mi)+ permsd*satdif*dcprde 
             dcqc(mi)= dcqc(mi)+ permsd*satdif*dcprdc
            else 
             deqc(mi) = 0.0
             dcqc(mi) = 0.0
           endif
           endif
        else if(qdis.gt.0.0.and.abs(qng_old).ge.flow_tol) then
c variable pressure(outflow)
c no specified air pressure (but specified airflow in or out)
          if(ieosd.eq.2) then
            hprod=sig*env+(1.0-sig)*enl
            dhprdp=dsigp*env+sig*dhvp-dsigp*enl+(1.0-sig)*dhlp
            dhprdc=dsigpc*env+sig*dhvpc-dsigpc*enl+(1.0-sig)*dhlpc
            dhprde=dsige*env+sig*dhvt-dsige*enl+(1.0-sig)*dhlt 
          else if(ieosd.eq.1) then
            hprod=enl
            dhprdp=dhlp
            dhprdc=dhlpc
            dhprde=dhlt
          else if(ieosd.eq.3) then
            hprod=env
            dhprdp=dhvp
            dhprdc=dhvpc
            dhprde=dhvt
          endif
          qh(mi)=qdis*hprod       
          dqh(mi)=dhprdp*qdis+hprod*dq(mi) 
          deqh(mi)=dhprde*qdis+hprod*dqt(mi) 
          dcqh(mi)=dhprdc*qdis
c 5-4-2001 air can now be in or out
c constant and no derivatives
          qc(mi) = qng_old
c derivatives zero out before 
c          dqc(mi)=0.0
c          deqc(mi)=0.0
c          dcqc(mi)=0.0
        else if(qdis.le.0.0.and.abs(qng_old).lt.flow_tol) then
c variable flow (pressure dependent)
c specified air pressure allowed(already accounted for)
c water flow into reservoir
          qh(mi)=qdis*eflow(mi)
          dqh(mi)=dq(mi)*eflow(mi)
        else if(qdis.le.0.0.and.abs(qng_old).ge.flow_tol) then
c variable flow (pressure dependent)
c specified air flow allowed(already accounted for)
c water flow into reservoir
          qh(mi)=qdis*eflow(mi)
          dqh(mi)=dq(mi)*eflow(mi)
          qc(mi)=qng_old  
c derivatives zero out before, may ne non zero (xairfl(mi))
c          dqc(mi)=0.0
c          deqc(mi)=0.0     
c          dcqc(mi)=0.0 
         endif     
c
c now add humidity terms (outflow)
c 
c specified saturation just acts on water in liquid phase
c and only outflow
c need energy term here
c
        if(sk_hum.gt.0.0) then
         sk(mi) = sk(mi) + sk_hum
         dqt(mi) = dqt(mi) + dsk_hums
         qh(mi) = qh(mi) + sk_hum*enl
         dqh(mi) = dqh(mi) + dhlp*sk_hum 
         deqh(mi)=deqh(mi) + enl*dsk_hums
         dcqh(mi)=dcqh(mi) + dhlpc*sk_hum  
c gaz 112518  inflow or outflow       
        else if(sk_hum.le.0.0) then
         sk(mi) = sk(mi) + sk_hum
         dqt(mi) = dqt(mi) + dsk_hums
         qh(mi) = qh(mi) + sk_hum*enl
         dqh(mi) = dqh(mi) + dhlp*sk_hum 
         deqh(mi)=deqh(mi) + enl*dsk_hums
         dcqh(mi)=dcqh(mi) + dhlpc*sk_hum            
        endif 
c air and water leave in the gas phase (outflow)
c need energy term here
        if(sk_air.gt.0.0) then
         sk(mi) = sk(mi) + (1.0-xnv)*sk_air
         dq(mi) = dq(mi) + (1.0-xnv)*dq_air - dxnvp*sk_air
         dqpc(mi) = dqpc(mi) - dxnvpc*sk_air
         qc(mi) = qc(mi) + sk_air*xnv
         dqc(mi) = dqc(mi) + xnv*dq_air + dxnvp*sk_air
         dcqc(mi) = dcqc(mi) + dxnvpc*sk_air
         qh(mi) = qh(mi) + sk_air*env
         dqh(mi)=dqh(mi) + env*dq_air + dhvp*sk_air
         dcqh(mi)=dcqh(mi) + dhvpc*sk_air      
        elseif (sk_air.lt.0.0) then
c air and water enter in the gas phase (inflow)
c includes specified humidity in inflow            
c need energy term here
c no derivatives wrt T or pci because inflow
c gaz corrected 111318 pci = p-pv(t) need more derivatives for ieosd =2            
         sk(mi) = sk(mi) + sk_h2o
         dq(mi) = dq(mi) + dsk_h2op
         dqt(mi) = dqt(mi) + dsk_h2ot
         dqpc(mi) = dqpc(mi) + dsk_h2opc
         qc(mi) = qc(mi) + sk_air
         dqc(mi) = dqc(mi) + dsk_airp
         deqc(mi) = deqc(mi) + dsk_airt         
         dcqc(mi) = dcqc(mi) + dsk_airpc
         qh(mi) = qh(mi)+ eh_air
         dqh(mi) = dqh(mi) + deh_airp
         deqh(mi) = deqh(mi) + deh_airt
         dcqh(mi) = dcqh(mi) + deh_airpc
       endif
c
c section on fixed relative humidity
c temperature should be constrained via hflux below
c
         if(huma_fixed.gt.huma_tol) then
          if(sk_humf.lt.0) then
             sk(mi) = sk(mi) + sk_humf
             dq(mi) = dq(mi) + dsk_humfp
             dqpc(mi) = dqpc(mi) + dsk_humfpc
             dqt(mi) = dqt(mi) + dsk_humft
             qc(mi) = qc(mi) + sk_airhf
             dqc(mi) = dqc(mi) + dsk_airhfp    
          else
             sk(mi) = sk(mi) + sk_humf
             dqpc(mi) = dqpc(mi) + dsk_humfpc  
             qc(mi) = qc(mi) + sk_airhf
             dqc(mi) = dqc(mi) + dsk_airhfp  
          endif
         endif
c
c
c
c ********** fluid source term code end *****************
c
c add heat source term
c
      htc = 0.0
      if(qflux(mi).ne.0.0.or.t_hum.ne.0.0) then
      if(qflxm(mi).gt.0.0) then
       tbound = qflux(mi)
       htc = qflxm(mi)
      elseif(t_hum.gt.0.0) then
       tbound = t_hum
       htc = abs(wellim(mi))
      endif
      if(htc.gt.0.0) then
       hflux=htc*(tl-tbound)
       if(ieosd.ne.2) then
        dhflxp=0.0
        dhflxe=htc
        dhflxc=0.0
       else
        dhflxp=0.0
        dhflxe=0.0
        dhflxc=htc
       endif
       qh(mi)=qh(mi)+hflux
       dqh(mi)=dqh(mi)+dhflxp
       deqh(mi)=deqh(mi)+dhflxe
       dcqh(mi)=dcqh(mi)+dhflxc
      else if(qflxm(mi).lt.0.0) then
c do nothing
      else
       qh(mi)=qh(mi)+qflux(mi)
      endif
      endif
c
c gaz 2-26-2002
      if(ieosd.eq.0) then
c heat conduction only
c
      cprd=cpr(mi)
      if(kq.lt.0) then
      eskd=eflow(mi)
      edif=cprd*tl-eskd
      permsd=wellim(mi)
      qh(mi)=permsd*edif
      deqh(mi)=permsd*cprd
      endif
      if(kq.ge.0.and.qflux(mi).eq.0.0) deqh(mi)=0.
      dtpae(mi)=1.
      denrd=denr(mi)*cprd
      a(mi+neq)=denrd*tl
      a(mi)=0.0
      a(mi+neq+neq)=0.0
      deef(mi)=denrd*dtin
      endif
      endif
  100 continue
c
c modify accumulation terms for volume changes
c
      if(ivfcal.eq.0) then
      do 105 mid=1,neq
      mi=mid+ndummy
      deni(mi)=(a(mi)   -denh(mi))*dtin
      denei(mi)=(a(mi+neq)-deneh(mi))*dtin
      denpci(mi)=(a(mi+neq+neq)-denpch(mi))*dtin
105   continue
      else
      do 110 mid=1,neq
      mi=mid+ndummy
      den=a(mi)
      dene=a(mi+neq)
      denc=a(mi+neq+neq)
      vfd=vfcal(mi,rl,dvfp,drlp)
      dmpf(mi)=dmpf(mi)*vfd+den*dvfp*dtin
      depf(mi)=depf(mi)*vfd+dene*dvfp*dtin
      dmef(mi)=dmef(mi)*vfd
      deef(mi)=deef(mi)*vfd
      dmc(mi)=dmc(mi)*vfd
      dec(mi)=dec(mi)*vfd
      dcp(mi)=dcp(mi)*vfd+denc*dvfp*dtin
      dce(mi)=dce(mi)*vfd
      dcc(mi)=dcc(mi)*vfd
      deni(mi)=(den*vfd-denh(mi))*dtin
      denei(mi)=(dene*vfd-deneh(mi))*dtin
      denpci(mi)=(denc*vfd-denpch(mi))*dtin
c must put in new relations for dil,dilp,dile
110   continue
      endif
c
c check for mass and energy balance decrepencies
c 
c     call check_balance_eq()
      do mi = 1, neq
        if(ieos_bal(mi).ne.0) then
        diff_w = deni(mi)-deni_ch(mi)
        diff_e = denei(mi)-denei_ch(mi)
        diff_a = denpci(mi)-denpci_ch(mi)
       endif
      enddo
c
c call air water-vapor diffusion
c
      call dvacalc

      if(allocated(s0)) then
         deallocate(s0,pcp0,dpcps0,rlf0)
         deallocate(drlfs0,rvf0,drvfs0)       
      endif

      return
      end
