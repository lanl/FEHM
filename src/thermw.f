      subroutine thermw(ndummy)
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
CD1 pure water simulation.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 NOV 1993     Z. Dash        22      Add prolog
CD2              G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/thermw.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Feb 16 13:01:54 1996   zvd
CD2 Added requirements.
CD2 
CD2    Rev 1.3   Fri Feb 02 12:48:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   Fri Jan 12 18:01:12 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.1   03/18/94 16:12:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:42   pvcs
CD2 original version in process of being certified
CD2 
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
CD4 n0, ipdeef, ipdepf, ipdmpf, ipdmef, ipdq, ipdqh, !dtot, !neq, 
CD4 qh, dq, dqt, dqh, deqh, cel, crl, cvl,
CD4 cev, crv, t, phi, s, sk, denr, cpr, ps, ka, volume, pflow,
CD4 wellim, eflow, dmef, dmpf, depf, deef,
CD4 dtpa, dtpac, dtpae, dstm, dil, dilp, dile, enlf, dglp, dgle,
CD4 delf, delef, div, divp, 
CD4 dive, envf, dgvp, dgve, devf, devef, rolf, rovf, qflux,
CD4 qflxm, iporos, ieos, iieos, cvv, dporp, dport,
CD4 b
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
CD4 rlperm          N/A      Generates relative permeabilities
CD4 porosi          N/A      Computes pressure dependent porosity terms
CD4 cappr           N/A      Generates capillary pressures
CD4 vfcal           real*8   Computes volume changes
CD4 wellbore        N/A      Compute wellbore terms when this option is
CD4                              specified
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
CD5 dtin         real*8      Reciprocal of time step
CD5 tfun         real*8      Term used in temperature calculation
CD5 tfunn        real*8      Term used in temperature calculation
CD5 tfund        real*8      Term used in temperature calculation
CD5 dtpsn        real*8      Term used in temperature calculation
CD5 dtpsd        real*8      Term used in temperature calculation
CD5 dporpl       real*8      Derivative of porosity with pressure
CD5 dportl       real*8      Derivative of porosity with temperature
CD5 dpldt        real*8      Derivative term
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
CD5 x4           real*8      x to the fourth
CD5 dpsats       real*8      Derivative of saturation pressure with
CD5                             saturation
CD5 dpsatt       real*8      Derivative of saturation pressure with
CD5                             temperature
CD5 dpct         real*8      Negative of dpsatt
CD5 dtps         real*8      Term used in calculation
CD5 tl2          real*8      tl squared
CD5 tl3          real*8      tl cubed
CD5 tlx          real*8      x times tl
CD5 tl2x         real*8      tl2 times x
CD5 tlx2         real*8      tl times x2
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
CD5 rnsn1        real*8      Term used in vapor density calculation
CD5 rnsn2        real*8      Term used in vapor density calculation
CD5 rnsn3        real*8      Term used in vapor density calculation
CD5 rnsd1        real*8      Term used in vapor density calculation
CD5 rnsd2        real*8      Term used in vapor density calculation
CD5 rnsd3        real*8      Term used in vapor density calculation
CD5 rnsn         real*8      Term used in vapor density calculation
CD5 rnsd         real*8      Term used in vapor density calculation
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
CD5 drovp       real*8       Derivative of vapor density with pressure
CD5 drovt       real*8       Derivative of vapor density with
CD5                              temperature
CD5 dhvst       real*8       Derivative of air enthalpy with temperature
CD5 dtd         real*8       Factor used in derivative calculations
CD5 sl          real*8       Saturation at this node
CD5 qdis        real*8       Source flow rate at this node
CD5 cp          real*8       Density time heat capacity
CD5 por         real*8       Porosity at this node
CD5 kq          int          Flag used in flow rate calculation
CD5 vol         real*8       Volume of this node
CD5 pldif       real*8       Delta P driving force for sink or source
CD5 permsd      real*8       Well impedance at current node
CD5 eskd        real*8       Energy flow from source at this node
CD5 eqdum       real*8       Energy stored at this node
CD5 rag         real*8       Ratio of vapor kinematic viscosities
CD5 sig         real*8       Ratio used in calculation
CD5 rop         real*8       Derivative of den with pressure
CD5 rop1        real*8       Derivative of den with pressure
CD5 damp        real*8       rop per unit time
CD5 daep        real*8       Parameter used in calculation
CD5 daep1       real*8       Parameter used in calculation
CD5 damh        real*8       Parameter used in calculation
CD5 daeh        real*8       Parameter used in calculation
CD5 dragp       real*8       Parameter used in calculation
CD5 dsigp       real*8       Derivative of sig with pressure
CD5 dsige       real*8       Derivative of sig with enthalpy
CD5 den1        real*8       Liquid density
CD5 roe         real*8       Term used in calculation
CD5 dql         real*8       Term used in calculation
CD5 hprod       real*8       Term used in calculation
CD5 dhprdp      real*8       Derivative term used in calculation
CD5 dhprde      real*8       Derivative term used in calculation
CD5 dcprdp      real*8       Derivative term used in calculation
CD5 htc         real*8       Heat flux impedance of source at this node
CD5 tbound      real*8       Heat flux of source at this node
CD5 hflux       real*8       Heat flux at this node
CD5 dhflxp      real*8       Derivative of heat flux
CD5 dhflxe      real*8       Derivative of heat flux
CD5 sbound      real*8       Saturation of boundary node
CD5 cprd        real*8       Heat capacity at this node
CD5 edif        real*8       Term used in calculation
CD5 denrd       real*8       Rock density heat capacity product
CD5 vfd         real*8       Term used in calculation
CD5 rl          real*8       Parameter used in calculation
CD5 dvfp        real*8       Parameter used in calculation
CD5 rop2        real*8       Parameter used in calculation
CD5 daep2       real*8       Parameter used in calculation
CD5 dqv         real*8       Parameter used in calculation
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.7 Sources and sinks
CD9 2.4.1 Pressure- and temperature-dependent water properties
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
CPS BEGIN thermw
CPS 
CPS rlperm - generate relative permeabilities
CPS porosi - compute pressure dependent porosity terms
CPS cappr - generate capillary pressures
CPS mmgetblk - allocate space in temporary arrays
CPS 
CPS FOR each node
CPS 
CPS   Initialize variables
CPS   Initialize flow terms
CPS   
CPS   IF coefficients are needed
CPS     Assign liquid enthalpy coefficients
CPS     Assign liquid density coefficients
CPS     Assign liquid viscosity coefficients
CPS     Assign vapor enthalpy coefficients
CPS     Assign vapor density coefficients
CPS     Assign vapor viscosity coefficients
CPS   ENDIF
CPS     
CPS   IF this is a two-phase mixture
CPS     Compute tmerpature and dt/dp
CPS   ELSE
CPS     Set temperature
CPS   ENDIF
CPS     
CPS   IF we are not superheated
CPS     Compute liquid enthalpy
CPS     Compute derivatives of liquid enthalpy
CPS     Compute liquid density
CPS     Compute derivatives of liquid density
CPS     Compute liquid viscosity
CPS     Compute derivatives of liquid viscosity
CPS   ENDIF
CPS       
CPS   IF we are not a compressed liquid
CPS     Compute vapor enthalpy
CPS     Compute derivatives of vapor enthalpy
CPS     Compute vapor density
CPS     Compute derivatives of vapor density
CPS     Compute vapor viscosity
CPS     Compute derivatives of vapor viscosity
CPS   ENDIF
CPS       
CPS   IF this is a two-phase mixture
CPS     Modify derivatives where needed
CPS   ENDIF
CPS       
CPS   IF a pressure-dependent source/sink is included
CPS     IF there is a flowing source term
CPS       Compute pressure difference and derivative
CPS     ELSE there is a flowing sink
CPS       psat - compute pressure difference and derivative
CPS     ENDIF
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
CPS     IF this is not a superheated vapor
CPS       Store derivatives of terms involving liquid in arrays
CPS     ENDIF this is not a compressed liquid
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
CPS       IF this is a superheated vapor system
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
CPS   IF there is a heat source/sink term at this node
CPS   
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
CPS mmrelblk - free space in temporary storage arrays
CPS 
CPS welbor - set volume factors if wellbore simulation is specified
CPS 
CPS END thermw
CPS
C**********************************************************************
C*****
C***** AF 11/15/10 updated for lookup table
c----------------------------------------------------------
cphs  adding lookup table capabilities      4/23/99
cphs  adding comtable.h 
cphs  tableFLAG = 1 means use the lookup table for 
cphs                water properties
c----------------------------------------------------------
catf
catf  modified to give more complete error messages 
catf  when lookup table errors out
catf
c-----------------------------------------------------------

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
      use comtable

      implicit none
C*****
C***** AF 11/15/10
c      include 'comtable.h'                   ! phs 4/23/99
C*****
      integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq
      real*8 dtin,dporpl,dportl,xrl,xrv,drl,drv,drlp,drvp,ela0,elpa1
      real*8 elpa2,elpa3,elta1,elta2,elta3,elpta,elp2ta,elpt2a
      real*8 elb0,elpb1,elpb2,elpb3,eltb1,eltb2,eltb3,elptb
      real*8 elp2tb,elpt2b,dla0,dlpa1,dlpa2,dlpa3,dlta1,dlta2,dlta3
      real*8 dlpta,dlp2ta,dlpt2a,dlb0,dlpb1,dlpb2,dlpb3,dltb1,dltb2
      real*8 dltb3,dlptb,dlp2tb,dlpt2b,vla0,vlpa1,vlpa2,vlpa3,vlta1
      real*8 vlta2,vlta3,vlpta,vlp2ta,vlpt2a,vlb0,vlpb1,vlpb2,vlpb3
      real*8 vltb1,vltb2,vltb3,vlptb,vlp2tb,vlpt2b,eva0,evpa1,evpa2
      real*8 evpa3,evta1,evta2,evta3,evpta,evp2ta,evpt2a,evb0,evpb1
      real*8 evpb2,evpb3,evtb1,evtb2,evtb3,evptb,evp2tb,evpt2b,dva0
      real*8 dvpa1,dvpa2,dvpa3,dvta1,dvta2,dvta3,dvpta,dvp2ta,dvpt2a
      real*8 dvb0,dvpb1,dvpb2,dvpb3,dvtb1,dvtb2,dvtb3,dvptb,dvp2tb
      real*8 dvpt2b,vva0,vvpa1,vvpa2,vvpa3,vvta1,vvta2,vvta3,vvpta
      real*8 vvp2ta,vvpt2a,vvb0,vvpb1,vvpb2,vvpb3,vvtb1,vvtb2,vvtb3
      real*8 vvptb,vvp2tb,vvpt2b,tl,pl,x,x2,x3,x4,dpsatt,tl2,tl3
      real*8 tlx,tl2x,tlx2,enwn1,enwn2,enwn3,enwn,enwd1,enwd2,enwd3
      real*8 enwd,enw,enl,dhwpn1,dhwpn2,dhwpn,dhwpd,dhwp,dhwtn1
      real*8 dhwtn2,dhwtn,dhwtd,dhwt,dhlt,dhlp,rnwn1,rnwn2,rnwn3
      real*8 rnwd1,rnwd2,rnwd3,rnwn,rnwd,rnw,rol,drlpn1,drlpn2
      real*8 drlpn,drolpd,drolp,drlen1,drlen2,drlen,droled,drolt,viln1
      real*8 viln2,viln3,viln,vild1,vild2,vild3,vild,vil,xvisl
      real*8 dvlpn1,dvlpn2,dvlpn,dvilpd,dvislp,dvlen1,dvlen2,dvlen
      real*8 dviled,dvislt,ensn1,ensn2,ensn3,ensn,ensd1,ensd2,ensd3
      real*8 ensd,ens,env,dhvp1,dhvp2,dhvpn,dhvpd,dhvt1,dhvt2
      real*8 dhvtn,dhvtd,dhvt,dhvp,rnsn1,rnsn2,rnsn3,rnsd1,rnsd2
      real*8 rnsd3,rnsn,rnsd,rns,rov,drspn1,drspn2,drspn,drospd
      real*8 drsen1,drsen2,drsen,drostd,visn1,visn2,visn3,visn,visd1
      real*8 visd2,visd3,visd,vis,xvisv,dvspn1,dvspn2,dvspn,dvispd
      real*8 dvisvp,dvsen1,dvsen2,dvsen,dvised,dvisvt,sl,qdis
      real*8 cp,por,vol,pldif,permsd,eskd,eqdum,rag,sig
      real*8 rop,rop1,daep1,dragp,dsigp,dsige
      real*8 den1,roe,dql,hprod,dhprdp,dhprde,htc,tbound,hflux
      real*8 dhflxp,sbound,cprd,edif,denrd,vfd,rl,dvfp,tfun
      real*8 tfunn,tfund,dtpsn,dtpsd,dpldt,psat,vfcal,rop2,daep2
      real*8 dqv,dhflxe,drovp,drovt
C*****
C***** AF 11/15/10
      real*8 zwp, zwt                                ! phs 4/23/99
      integer indexp, indext, point(4),izerrFLAG     ! phs 4/23/99
C*****
      real*8 xa,xa2,xa3,xa4
      real*8 p_energy,prop,dpropt,dpropp
      real*8, allocatable :: sto1(:)
      real(8) :: damh = 0., damp = 0., daep = 0., daeh = 0.
      real(8) :: dtps = 0., dtd = 0.

      integer i_mem_rlp
      save i_mem_rlp
      
      real*8, allocatable :: s0(:)
      real*8, allocatable :: pcp0(:)
      real*8, allocatable :: dpcps0(:)
      real*8, allocatable :: rlf0(:)
      real*8, allocatable :: drlfs0(:)
      real*8, allocatable :: rvf0(:)
      real*8, allocatable :: drvfs0(:)
      
      allocate(sto1(n0*2))
c     
c     rol  -  density liquid
c     rov  -  density vapour
c     enl  -  enthalpy liquid
c     env  -  enthalpy vapour
c     visl -  viscosity of liquid
c     visv -  viscosity of vapour
c     rl   -  relative permeability of liquid phase
c     rv   -  relative permeability of vapour phase
c     tfun -  temperature
c     sw   -  saturation liquid
c     
****  Avg molecular weight is set to molecular weight of water ****
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

c     generate relative permeabilities
c     calculate pressure dependant porosity and derivatives
      if(iporos.ne.0) call porosi(1)
c     call capillary pressure models
      if (.not. rlpnew) call cappr(1,ndummy)
      ifree1 = 0
      do mid=1,neq
         mi=mid+ndummy
         avgmolwt(mi) = mw_water
         ieosd=ieos(mi)
         iieosd=iieos(mi)
         if(ieosd.ge.2) ifree1 = ifree1 + 1

c     undo equivalence relations for relative perms
         xrl=rlf(mi)
         drl=drlef(mi)
         drlp=drlpf(mi)
         xrv=rvf(mi)
         drv=drvef(mi)
         drvp=drvpf(mi)
         if(igrav.ne.0) then
            p_energy = -grav*cord(mi,igrav)
         else
            p_energy = 0.0d0
         endif
c check for aux eos(iieosd.gt.10)
         if(itsat.le.10) then
   
c     adjust coefficients for thermo fits
         if(iieosd.ne.iieosl) then

c     liquid phase coefficients
c     
c     liquid enthalpy
c     numerator coefficients
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
c     denomenator coefficients
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
c     liquid density
c     numerator coefficients
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
c     denomenator coefficients
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
c     liquid viscosity
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor phase coefficients
c     
c     vapor enthalpy
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor density
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor viscosity
c     numerator coefficients
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
c     denomenator coefficients
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
         pl=phi(mi)

c     evaluate thermo functions and derivatives
         x=pl - phi_inc
         x2=x*x
         x3=x2*x
         x4=x3*x
         xa=pl
         xa2=xa*xa
         xa3=xa2*xa
         xa4=xa3*xa 
         if(ieosd.eq.2) then

c     two phase conditions

c     calculate temperature and dt/dp
            tfunn=tsa0+tspa1*xa+tspa2*xa2+tspa3*xa3+tspa4*xa4
            tfund=tsb0+tspb1*xa+tspb2*xa2+tspb3*xa3+tspb4*xa4
            tfun=tfunn/tfund
            tl=tfun
            dtpsn=((tspa1+2.*tspa2*xa+3.*tspa3*xa2+4.*tspa4*xa3)*tfund)-
     &           (tfunn*(tspb1+2.*tspb2*xa+3.*tspb3*xa2+4.*tspb4*xa3))
            dtpsd=tfund**2
            dtps=dtpsn/dtpsd
         else
            tl=t(mi)
         endif
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         if(ieosd.ne.3) then
C*****
C*****AF 11/15/10
c-----------------------------------------
c     phs      Lookup table if tableFLAG = 1
c-----------------------------------------
c     
            if(tableFLAG.NE.1) then
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

c     derivatives of enthalpy
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

c     liquid density
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
               if(cden) rol = rol+factcden*anl((ispcden-1)*n0+mi)


c     derivatives of density
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

c     liquid viscosity
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

c     derivatives of liquid viscosity
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
C*****
C*****AF 11/15/10
c--------------------------------------------------------------------------
            else    ! USE LOOKUP TABLE       phs 4/23/99         LOOKUP

               izerrFLAG = 0.

               indexp = pmin(1) + incp*dint((pl-pmin(1))/incp)
               indext = tmin(1) + inct*dint((tl-tmin(1))/inct)
c---  find 4 points            LOOKUP

               point(1) = 1 + (((indexp-pmin(1))/incp)*numt)
     x              +  ((indext-tmin(1))/inct)
               point(2) = point(1) + numt
               point(3) = point(2) + 1
               point(4) = point(1) + 1
c---  
               if(PP(point(1),3).LT.0) izerrFLAG = 1.
               if(PP(point(2),3).LT.0) izerrFLAG = 1.
               if(PP(point(3),3).LT.0) izerrFLAG = 1.
               if(PP(point(4),3).LT.0) izerrFLAG = 1.
               if(izerrFLAG.EQ.1.) then
                  write (ierr, 10)
                  write (ierr, 20) pl, tl
                  write (ierr, 30) pmin(1), tmin(1)
                  write (ierr, 40) point(1),point(2),point(3),point(4)
                  write (ierr, 50) PP(point(1),3), PP(point(4),3)
                  write (ierr, 60) PP(point(2),3), PP(point(3),3)
                  if (iptty .ne. 0) then
c STOP execution if any lookuppoints are out of bounds.
c     *****AF
                     write (iptty, 10)
                     write (iptty, 20) pl, tl
                     write (iptty, 30) pmin(1), tmin(1)
                     write (iptty, 40) point(1),point(2),point(3),
     $                 point(4)
                     write (iptty, 50) PP(point(1),3), PP(point(4),3)
                     write (iptty, 60) PP(point(2),3), PP(point(3),3)
                  end if
c     *****AF
c     *****AF  It would be ideal to just set values to the limits when the limits
c     *****AF  are exceeded, then tell the user about the problem but let the run
c     *****AF  continue. Perhaps add this later...
c     *****AF
                  stop                                      
               endif
 10            format ('Out of bounds in Thermw.f')
 20            format ('Target values of P and T:', 2(1x, g16.9))
 30            format ('Min values pmin(1) and tmin(1):', 2(1x, g16.9))
 40            format ('point(1), (2), (3), (4):', 4(1x, i6))
 50            format ('Rho at T before/after:', 2(1x,g16.9))
 60            format ('Rho at P before/after:', 2(1x,g16.9))

c---  compute weights for the P and T direction
               zwp = (pl-indexp)/incp
               zwt = (tl-indext)/inct
c---  find values as function of P+T                 LOOKUP

               rol   = (1-zwp)*(1-zwt)*PP(point(1),3) + (1-zwt)*zwp*
     &              PP(point(2),3) + zwt*zwp*PP(point(3),3) +
     &              (1-zwp)*zwt*PP(point(4),3)

               drolp = (1-zwp)*(1-zwt)*PP(point(1),4) + (1-zwt)*zwp*
     &              PP(point(2),4) + zwt*zwp*PP(point(3),4) +
     &              (1-zwp)*zwt*PP(point(4),4)

               drolt = (1-zwp)*(1-zwt)*PP(point(1),5) + (1-zwt)*zwp*
     &              PP(point(2),5) + zwt*zwp*PP(point(3),5) +
     &              (1-zwp)*zwt*PP(point(4),5)

               enl   = (1-zwp)*(1-zwt)*PP(point(1),6) + (1-zwt)*zwp*
     &              PP(point(2),6) + zwt*zwp*PP(point(3),6) +
     &              (1-zwp)*zwt*PP(point(4),6)

               dhlp  = (1-zwp)*(1-zwt)*PP(point(1),7) + (1-zwt)*zwp*
     &              PP(point(2),7) + zwt*zwp*PP(point(3),7) +
     &              (1-zwp)*zwt*PP(point(4),7)

               dhlt  = (1-zwp)*(1-zwt)*PP(point(1),8) + (1-zwt)*zwp*
     &              PP(point(2),8) + zwt*zwp*PP(point(3),8) +
     &              (1-zwp)*zwt*PP(point(4),8)

               xvisl = (1-zwp)*(1-zwt)*PP(point(1),9) + (1-zwt)*zwp*
     &              PP(point(2),9) + zwt*zwp*PP(point(3),9) +
     &              (1-zwp)*zwt*PP(point(4),9)

               dvislp = (1-zwp)*(1-zwt)*PP(point(1),10) + (1-zwt)*zwp*
     &              PP(point(2),10) + zwt*zwp*PP(point(3),10) +
     &              (1-zwp)*zwt*PP(point(4),10)

               dvislt = (1-zwp)*(1-zwt)*PP(point(1),11) + (1-zwt)*zwp*
     &              PP(point(2),11) + zwt*zwp*PP(point(3),11) +
     x              (1-zwp)*zwt*PP(point(4),11)  ! LOOKUP


            end if              !  tableFLAG.NE.1           phs 4/23/99
c-----------------------------------------------------------------------
C*****
         endif

         if(ieosd.ne.1) then

c     vapor enthalpy
            ensn1=eva0+evpa1*x+evpa2*x2+evpa3*x3
            ensn2=evta1*tl+evta2*tl2+evta3*tl3
            ensn3=evpta*tlx+evpt2a*tl2x+evp2ta*tlx2
            ensn=ensn1+ensn2+ensn3
            ensd1=evb0+evpb1*x+evpb2*x2+evpb3*x3
            ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
            ensd3=evptb*tlx+evpt2b*tl2x+evp2tb*tlx2
            ensd=ensd1+ensd2+ensd3
            ens=ensn/ensd
            env=ens + p_energy

c     derivatives of vapor enthalpy
            dhvp1=evpa1+2*evpa2*x+3*evpa3*x2+evpta*tl
            dhvp1=ensd*(dhvp1+2*evp2ta*tlx+evpt2a*tl2)
            dhvp2=evpb1+2*evpb2*x+3*evpb3*x2+evptb*tl
            dhvp2=ensn*(dhvp2+2*evp2tb*tlx+evpt2b*tl2)
            dhvpn=dhvp1-dhvp2
            dhvpd=ensd**2
            dhvp=dhvpn/dhvpd
            dhvt1=evta1+2*evta2*tl+3*evta3*tl2+evpta*x
            dhvt1=ensd*(dhvt1+evp2ta*x2+2*evpt2a*tlx)
            dhvt2=evtb1+2*evtb2*tl+3*evtb3*tl2+evptb*x
            dhvt2=ensn*(dhvt2+evp2tb*x2+2*evpt2b*tlx)
            dhvtn=dhvt1-dhvt2
            dhvtd=ensd**2
            dhvt=dhvtn/dhvtd

c     vapor density
            rnsn1=dva0+dvpa1*x+dvpa2*x2+dvpa3*x3
            rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
            rnsn3=dvpta*tlx+dvpt2a*tl2x+dvp2ta*tlx2
            rnsn=rnsn1+rnsn2+rnsn3
            rnsd1=dvb0+dvpb1*x+dvpb2*x2+dvpb3*x3
            rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
            rnsd3=dvptb*tlx+dvpt2b*tl2x+dvp2tb*tlx2
            rnsd=rnsd1+rnsd2+rnsd3
            rns=rnsn/rnsd
            rov=rns

c     derivatives of vapor density
            drspn1=dvpa1+2*dvpa2*x+3*dvpa3*x2+dvpta*tl
            drspn1=rnsd*(drspn1+2*dvp2ta*tlx+dvpt2a*tl2)
            drspn2=dvpb1+2*dvpb2*x+3*dvpb3*x2+dvptb*tl
            drspn2=rnsn*(drspn2+2*dvp2tb*tlx+dvpt2b*tl2)
            drspn=drspn1-drspn2
            drospd=rnsd**2
            drovp=drspn/drospd
            drsen1=dvta1+2*dvta2*tl+3*dvta3*tl2+dvpta*x
            drsen1=rnsd*(drsen1+dvp2ta*x2+2*dvpt2a*tlx)
            drsen2=dvtb1+2*dvtb2*tl+3*dvtb3*tl2+dvptb*x
            drsen2=rnsn*(drsen2+dvp2tb*x2+2*dvpt2b*tlx)
            drsen=drsen1-drsen2
            drostd=rnsd**2
            drovt=drsen/drostd

c     vapor viscosity
            visn1=vva0+vvpa1*x+vvpa2*x2+vvpa3*x3
            visn2=vvta1*tl+vvta2*tl2+vvta3*tl3
            visn3=vvpta*tlx+vvpt2a*tl2x+vvp2ta*tlx2
            visn=visn1+visn2+visn3
            visd1=vvb0+vvpb1*x+vvpb2*x2+vvpb3*x3
            visd2=vvtb1*tl+vvtb2*tl2+vvtb3*tl3
            visd3=vvptb*tlx+vvpt2b*tl2x+vvp2tb*tlx2
            visd=visd1+visd2+visd3
            vis=visn/visd
            xvisv=vis

c     derivatives of vapor viscosity
            dvspn1=vvpa1+2*vvpa2*x+3*vvpa3*x2+vvpta*tl
            dvspn1=visd*(dvspn1+2*vvp2ta*tlx+vvpt2a*tl2)
            dvspn2=vvpb1+2*vvpb2*x+3*vvpb3*x2+vvptb*tl
            dvspn2=visn*(dvspn2+2*vvp2tb*tlx+vvpt2b*tl2)
            dvspn=dvspn1-dvspn2
            dvispd=visd**2
            dvisvp=dvspn/dvispd
            dvsen1=vvta1+2*vvta2*tl+3*vvta3*tl2+vvpta*x
            dvsen1=visd*(dvsen1+vvp2ta*x2+2*vvpt2a*tlx)
            dvsen2=vvtb1+2*vvtb2*tl+3*vvtb3*tl2+vvptb*x
            dvsen2=visn*(dvsen2+vvp2tb*x2+2*vvpt2b*tlx)
            dvsen=dvsen1-dvsen2
            dvised=visd**2
            dvisvt=dvsen/dvised
         endif

c     modify derivatives for 2-phase
         if(ieosd.eq.2) then
            drolp=drolp+drolt*dtps
            drovp=drovp+drovt*dtps
            dhlp=dhlp+dhlt*dtps
            dhvp=dhvp+dhvt*dtps
            dvisvp=dvisvp+dvisvt*dtps
            dvislp=dvislp+dvislt*dtps
            drolt=0.0
            drovt=0.0
            dhlt=0.0
            dhvt=0.0
            dvisvt=0.0
            dvislt=0.0
            dtd=0.0
         endif
c exclude 2phase in aux eos because of smeared phase lines         
         else
          ieosd = 1
          tl = t(mi)
          pl = phi(mi)
c density and derivatives          
          call eos_aux(itsat,tl,pl,0,0,prop,dpropt,dpropp)
          rol = prop
          drolt = dpropt
          drolp = dpropp
c enthalpy and derivatives (add potential energy term)          
          call eos_aux(itsat,tl,pl,1,1,prop,dpropt,dpropp)
          enl = prop + p_energy
          dhlt = dpropt
          dhlp = dpropp               
c viscosity and derivatives           
          call eos_aux(itsat,tl,pl,2,1,prop,dpropt,dpropp)
          xvisl = prop
          dvislt = dpropt
          dvislp = dpropp         
         endif

         sl=s(mi)
         dq(mi)=0.0
         qh(mi)=0.0
         dqh(mi)=0.0
         deqh(mi)=0.0
         dqt(mi)=0.0
         qdis=sk(mi)
         cp=denr(mi)*cpr(mi)
         por=ps(mi)
         dporpl=dporp(mi)
         dportl=dport(mi)
         kq=ka(mi)
         vol=volume(mi)

c     form pressure dependent flow term
         if(kq.lt.0 .and. compute_flow) then
            qh(mi)=0.
            sk(mi)=0.0
            if(pflow(mi).gt.0.0) then
               pldif=pl-pflow(mi)
               dpldt=0.0
            else
               pldif=pl-(psat(tl,dpsatt,0)-pflow(mi))
               dpldt=-dpsatt
            endif
            permsd=abs(wellim(mi))
            if(iwelimd.ne.0)then
c     
c     peaceman solution only available for models kq(-1,-2,1)
c     (need) to put derivative terms in later (gaz)  
               if(izonewel1(mi).ne.0) then
                  permsd = wellim(mi)*rol/xvisl
               endif
            endif            
            if(pldif.le.0.0d00.and.kq.eq.-2) permsd=0.0d0
            qdis=permsd*(pldif)
            sk(mi)=qdis
            dq(mi)=permsd
            if(qdis.le.0.) then
               eskd=eflow(mi)
               qh(mi)=eskd*qdis
               sk(mi)=qdis
               dqh(mi)=eskd*dq(mi)
            endif
            dqt(mi)=dpldt*permsd
         endif

c     two phase conditions
         if(ieosd.eq.2) then
c     saturation and relative permeabilities
            sv=1.d0-sl
c     accumulation terms
            den=por*(sl*rol+sv*rov)
            eqdum=sl*rol*enl+sv*rov*env-pl
            dene=((1.-por)*cp*tl+por*eqdum)
c     production of steam
            rag=rol*xvisv/rov/xvisl
            sig=xrv/(xrv+rag*xrl)
c     derivatives of accumulation terms
            rop=por*(sv*drovp+sl*drolp)
            damp =rop*dtin
            daep =((1.d0-por)*cp*dtps+por*(sv*drovp*env+sv*rov*dhvp+
     &           sl*drolp*enl+sl*rol*dhlp)-por)*dtin
            damh =por*(rol-rov)*dtin
            daeh =por*(rol*enl-rov*env)*dtin
c     derivatives of sink terms
            if(qdis.gt.0.0) then
               dragp=(rov*xvisl*(drolp*xvisv+rol*dvisvp)-rol*
     2              xvisv*(drovp
     3              *xvisl+rov*dvislp))/(xvisl*rov)**2
               if(sig.ne.0.d0.and.sig.ne.1.d0) then
                  dsigp=-xrv*dragp*xrl/(xrv+rag*xrl)**2
                  dsige=((xrv+rag*xrl)*drv-xrv*(drv+rag*drl))
     3                 /(xrv+rag*xrl)**2
               else
                  dsigp=0.d0
                  dsige=0.d0
               end if
            end if
         end if

c     compressed liquid
         if(ieosd.eq.1) then
            dtps=0.0
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
            dtd=1.0
            xvisv=1.d20

c     accumulation terms
            den1=rol
            den=den1*por
            dene=(1.d0-por)*cp*tl+den*enl-por*pl

c     derivatives of accumulation terms
            rop=por*drolp
            rop2=rol*dporpl
            rop1=rop
            damp=(rop1+rop2)*dtin
            daep1=(rop*enl+por*rol*dhlp-por)
            daep2=dporpl*(-cp*tl+rol*enl-pl)
            daep=(daep1+daep2)*dtin
            roe=por*drolt
            damh=(rol*dportl+drolt*por)*dtin
            daeh=(-dportl*cp*tl+(1.-por)*cp+
     &           dportl*(enl*rol-pl)+por*(rol*dhlt+enl*drolt))*dtin

         end if

c     superheated vapour
         if(ieosd.eq.3) then
            dtps=0.0
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
            dtd=1.0
            xvisl=1.d20

c     accumulation terms
            den1=rov
            den=den1*por
            dene=(1.d0-por)*cp*tl+den*env-por*pl

c     derivatives of accumulation terms
            rop=por*drovp
            damp =rop*dtin
            daep=(env*rop-por)*dtin
            roe=por*drovt
            damh =roe*dtin
            daeh =((1.d0-por)*cp+por*rov*dhvt+env*roe)*dtin
         end if
         if(ieosd.ne.-1) then

c     store derivatives of accumulation terms
            dmef(mi)=damh
            dmpf(mi)=damp
            depf(mi)=daep
            deef(mi)=daeh
            dtpa(mi)=dtps
            dtpae(mi)=dtd
            t(mi)=tl
            s(mi)=sl

c     save accumulation terms for possible volume changes
            sto1(mi)=den
            sto1(mi+neq)=dene
            dstm(mi)=por*rov*sv*vol
            dil(mi)=0.0
            dilp(mi)=0.0
            dile(mi)=0.0
            div(mi)=0.0
            divp(mi)=0.0
            dive(mi)=0.0
            if(ieosd.ne.3) then

c     modify flow terms for new upwind scheme
               dql=rol*xrl/xvisl
               dil(mi)=dql
               dilp(mi)=drolp*xrl/xvisl-rol*xrl/xvisl**2*dvislp
     *              +rol*drlp/xvisl
               dile(mi)=drolt*xrl/xvisl+rol*drl/xvisl
     *              -rol*xrl/xvisl**2*dvislt
            endif
            dglp(mi)=drolp
            dgle(mi)=drolt
            enlf(mi)=enl
            delf(mi)=dhlp
            delef(mi)=dhlt

c     modify flow terms for new upwind scheme
            if(ieosd.ne.1) then
               dqv=rov*xrv/xvisv
               div(mi)=dqv
               divp(mi)=drovp*xrv/xvisv-rov*xrv/xvisv**2*dvisvp
     *              +rov*drvp/xvisv
               dive(mi)=drovt*xrv/xvisv+rov*drv/xvisv
     *              -rov*xrv/xvisv**2*dvisvt
            endif
            dgvp(mi)=drovp
            dgve(mi)=drovt
            envf(mi)=env
            devf(mi)=dhvp
            devef(mi)=dhvt
            rolf(mi)=rol
            rovf(mi)=rov

            if(qdis.gt.0.0) then

c     organize source terms and derivatives
               if(ieosd.eq.2) then
                  hprod=sig*env+(1.0-sig)*enl
                  dhprdp=dsigp*env+sig*dhvp-dsigp*enl+(1.0-sig)*dhlp
                  dhprde=dsige*(env-enl)
               endif
               if(ieosd.eq.1) then
                  hprod=enl
                  dhprdp=dhlp
                  dhprde=dhlt
               endif
               if(ieosd.eq.3) then
                  hprod=env
                  dhprdp=dhvp
                  dhprde=dhvt
               endif
               qh(mi)=hprod*qdis
               dqh(mi)=dhprdp*qdis+hprod*dq(mi)
               deqh(mi)=dhprde*qdis+hprod*dqt(mi)
            endif
            if(qdis.lt.0.) then
               qh(mi)=qdis*eflow(mi)
               dqh(mi)=dq(mi)*eflow(mi)
            endif
         endif

c     add heat source term
         if(qflux(mi).ne.0.0) then
            if(qflxm(mi).gt.0.0) then
               htc=qflxm(mi)
               tbound=qflux(mi)
               hflux=htc*(tl-tbound)
               if(ieosd.ne.2) then
                  dhflxp=0.0
                  dhflxe=htc
               else
                  dhflxp=htc*dtps
                  dhflxe=0.0
               endif
               qh(mi)=qh(mi)+hflux
               dqh(mi)=dqh(mi)+dhflxp
               deqh(mi)=deqh(mi)+dhflxe
            else if(qflxm(mi).lt.0.0) then
               htc=abs(qflxm(mi))
               sbound=qflux(mi)
               hflux=htc*(tl-sbound)
               if(ieosd.ne.2) then
                  dhflxp=0.0
                  dhflxe=0.0
               else
                  dhflxp=0.0
                  dhflxe=htc
               endif
               qh(mi)=qh(mi)+hflux
               dqh(mi)=dqh(mi)+dhflxp
               deqh(mi)=deqh(mi)+dhflxe
            else
               qh(mi)=qh(mi)+qflux(mi)
            endif
         endif

c     gaz 2-26-2002
         if(ieosd.eq.0) then

c     heat conduction only
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
            sto1(mi)=0.0
            sto1(mi+neq)=denrd*tl
            deef(mi)=denrd*dtin
         endif
         if(kq.eq.1) then
c     
c     check for peaceman calculation
c     
            if(iwelimd.ne.0)then        
               if(izonewel1(mi).ne.0) then
                  permsd = wellim(mi)*rolf(mi)/xvisl
               endif
               pflow(mi)=phi(mi) - sk(mi)/permsd
            endif
         endif               
      enddo

c     modify accumulation terms for volume changes
      if(ivfcal.eq.0) then
         do 105 mid=1,neq
            mi=mid+ndummy
            deni(mi)=(sto1(mi)  -denh(mi))*dtin
            denei(mi)=(sto1(mi+neq)-deneh(mi))*dtin
 105     continue
      else
         do 110 mid=1,neq
            mi=mid+ndummy
            den=sto1(mi)
            dene=sto1(mi+neq)
            vfd=vfcal(mi,rl,dvfp,drlp)
            dmpf(mi)=dmpf(mi)*vfd+den*dvfp*dtin
            depf(mi)=depf(mi)*vfd+dene*dvfp*dtin
            dmef(mi)=dmef(mi)*vfd
            deef(mi)=deef(mi)*vfd
            deni(mi)=(den*vfd-denh(mi))*dtin
            denei(mi)=(dene*vfd-deneh(mi))*dtin
c     need new relations for dil,dilp,dile
 110     continue
      endif
c     call welbore to set volume factors if wellbore simulation is
c     specified
c     call welbor(3)
      
      deallocate(sto1)
      if(allocated(s0)) then
         deallocate(s0,pcp0,dpcps0,rlf0)
         deallocate(drlfs0,rvf0,drvfs0)       
      endif

      return
      end
