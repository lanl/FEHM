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
      use com_nondarcy
c gaz 101221
      use com_exphase
c gaz 101321
      use com_prop_data, only : den_h2o, enth_h2o, visc_h2o, humid_h2o,
     & psat_h2o, den_ngas, enth_ngas, visc_ngas, xnl_ngas, ieval_flag  
     
      use property_interpolate_1
      use comsi, only : ihms, density, internal_energy
      use comtable

      implicit none
C*****
C***** AF 11/15/10
c      include 'comtable.h'                   ! phs 4/23/99
C*****
      integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq
      real*8 psatl,dtsatp,dpsats
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
      real*8 cden_correction, cden_cor
c gaz  081317,082917  
      real*8 ur, dur_dt  
c gaz 101321  
      integer ipv_tpl
      real*8 ros,pv,dtdp,dpvt,dpct
      real*8 dum_gaz
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
c gaz 110715
      real*8 dum1,dumb,dumc,value(9)
      integer istate, ifail
c gaz 101221      
      integer loop_start, loop_end
      
      integer i_mem_rlp
      save i_mem_rlp
      
      real*8, allocatable :: s0(:)
      real*8, allocatable :: pcp0(:)
      real*8, allocatable :: dpcps0(:)
      real*8, allocatable :: rlf0(:)
      real*8, allocatable :: drlfs0(:)
      real*8, allocatable :: rvf0(:)
      real*8, allocatable :: drvfs0(:)
c gaz 101221      
cDEC$ FIXEDFORMLINESIZE:132      
      
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
c gaz 090623 rlpm requires s(:) to be allocated
c try to skip if heat conduction only (idoff = -1)
      if(idoff.ne.-1) then
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
c endif for heat conduction
      endif
      iieosl=0
      dtin=1.0/dtot
c gaz 020425 modified for non darcy flow
      if(nd_flow) then
       rlf_nd(1:n0) = rlf(1:n0)
       rvf_nd(1:n0) = rvf(1:n0)
       drlef_nd(1:n0) = drlef(1:n0)
       drvef_nd(1:n0) = drvef(1:n0)
      endif
c     generate relative permeabilities
c     calculate pressure dependant porosity and derivatives
      if(iporos.ne.0) call porosi(1)
c     call capillary pressure models
      if (.not. rlpnew) call cappr(1,ndummy)
c gaz 103017      
c  call variable rock properties if appropriate  
c  rock state started at last time step temperatures     
      if(iad.eq.0) call vrock_ctr(3,0)
      call vrock_ctr(1,0)
      call vrock_ctr(2,0)
c gaz 101221 added fluid control module
c initialize and allocate memory   
       call fluid_props_control(0, 0, 0, 
     &   fluid(1), 'all      ', '         ')   
c gaz 123020 manage explicit update 
       if(i_ex_update.ne.0.and.ieq_ex.gt.0) then
         loop_start = ieq_ex
         loop_end = ieq_ex
       else
         loop_start = 1
         loop_end = neq          
       endif       
      ifree1 = 0
c       do 100 mid=loop_start,loop_end
c gaz 101221
       do mid=loop_start,loop_end
         mi=mid+ndummy
         if(igrav.ne.0) then
          p_energy = -grav*cord(mi,igrav)
         else
          p_energy = 0.0d0
         endif
         avgmolwt(mi) = mw_water
         ieosd=ieos(mi)
c gaz 111415 modification to include supercritical
         if(ieosd.eq.4) ieosd = 1
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
            if (ieosd.eq.2) then
             xrl=rlf(mi)
             drl=drlef(mi)
             drlp=drlpf(mi)
             xrv=rvf(mi)
             drv=drvef(mi)
             drvp=drvpf(mi)
            else if (ieosd .eq. 1) then
               xrl = 1.0 
               xrv = 0.0
            else if (ieosd.eq.3) then
               xrl = 0.0
               xrv = 1.0
            else if (ieosd.eq.4) then
c
c gaz 102621 sc phase associated with liquid phase                 
c              
                xrl = 1.0
                xrv =0.0
            end if
            drl = 0.0
            drlp = 0.0
            drv = 0.0
            drvp = 0.0
         end if

c check for aux eos(iieosd.gt.10)
         if(itsat.le.10) then
          pl = phi(mi)
          tl = t(mi)

c       go to 699
       kq = l
       call fluid_props_control(1, mi, mi,
     &   fluid(1), 'all      ', '         ')
c     
c gaz 101321 set variables to new code values from fluid_props_control
c
        rol   =  den_h2o(mi,1) 
        drolp	=  den_h2o(mi,2) 
        drolt	=  den_h2o(mi,3) 
        rov	=  den_h2o(mi,4) 
        ros   =  rov
        drovp	=  den_h2o(mi,5) 
        drovt =  den_h2o(mi,6) 
        enl	=  enth_h2o(mi,1)
        dhlp	=  enth_h2o(mi,2)
        dhlt	=  enth_h2o(mi,3)
        env	=  enth_h2o(mi,4)
        ens   =  env
        dhvp	=  enth_h2o(mi,5)
        dhvt	=  enth_h2o(mi,6)
        xvisl   =  visc_h2o(mi,1)
        dvislp  =  visc_h2o(mi,2)
        dvislt  =  visc_h2o(mi,3)
        xvisv   =  visc_h2o(mi,4)
        vis     =  xvisv
        dvisvp  =  visc_h2o(mi,5)  
        dvisvt  =  visc_h2o(mi,6)
        
      
       pv         = psat_h2o(mi,1)  
       dtdp       = psat_h2o(mi,2)
       dpct       = psat_h2o(mi,3)
c gaz 100721 (dpsats needed for vapor pressure lowering)       
       dpsats     = psat_h2o(mi,4)      
       dtps       = dtdp       

c gaz 072520 (different from ngas)    
      if(ieosd.eq.3)then
       ieosd = 3
       xrv = 1
c gaz 102621 possible error   xrl = 1     
       xrl = 0
       drl = 0
       drv = 0
       drlp = 0
       drvp = 0 
      endif                  
c      
c
c gaz moved cden correction here
c
699    continue
      if(cden) then
c     Add correction for liquid species
         cden_cor = cden_correction(mi)
         rol = rol + cden_cor
      end if

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
c gaz 081317      (gaz moved 102917)   
         if(ivrock.ne.0) then
          ur = urock(mi)   
          dur_dt = durockt(mi)
         else
          urock(mi) = denr(mi)*cpr(mi)*tl
          durockt(mi) = denr(mi)*cpr(mi)    
          ur = urock(mi)   
          dur_dt  = durockt(mi)
         endif
c gaz 110123 modified for heat conduction
         if(idoff.eq.-1) then
          sl = 1.0
         else
          sl=s(mi)
         endif
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
            dene=((1.-por)*ur+por*eqdum)
c     production of steam
            rag=rol*xvisv/rov/xvisl
            sig=xrv/(xrv+rag*xrl)
c     derivatives of accumulation terms
            rop=por*(sv*drovp+sl*drolp)
            damp =rop*dtin
c gaz 081317            
            daep =((1.d0-por)*(dur_dt*dtps)+
     &           por*(sv*drovp*env+sv*rov*dhvp+       
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
            dene=(1.d0-por)*ur+den*enl-por*pl

c......................................................
c s kelkar, 28 feb 2011, for derivatives pore volume wrt displacements
            if(ihms.eq.15.or.ihms.eq.17) then
               density(mi) = rol*dtin
               internal_energy(mi) = (rol*enl-pl)*dtin
            endif
c.....................................................

c     derivatives of accumulation terms
            rop=por*drolp
            rop2=rol*dporpl
            rop1=rop
            damp=(rop1+rop2)*dtin
            daep1=(rop*enl+por*rol*dhlp-por)
            daep2=dporpl*(-ur+rol*enl-pl)
            daep=(daep1+daep2)*dtin
            roe=por*drolt
            damh=(rol*dportl+drolt*por)*dtin
c gaz 081317            
            daeh=(-dportl*ur+(1.-por)*dur_dt +
     &       dportl*(enl*rol-pl)+por*(rol*dhlt+enl*drolt))*dtin
            
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
            dene=(1.d0-por)*ur+den*env-por*pl

c     derivatives of accumulation terms
            rop=por*drovp
            damp =rop*dtin
            daep=(env*rop-por)*dtin
            roe=por*drovt
            damh =roe*dtin
c gaz 081317            
            daeh=((1.d0-por)*dur_dt+por*rov*dhvt+env*roe)*dtin          
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
c gaz 110123
            if(idoff.ne.-1) then
             s(mi)=sl
            endif

c     save accumulation terms for possible volume changes
            sto1(mi)=den
            sto1(mi+neq)=dene
            dstm(mi)=por*rov*sv*vol
            dil(mi)=0.0
            dilp(mi)=0.0
            dile(mi)=0.0
c gaz 110123
            if(idoff.ne.-1) then
             div(mi)=0.0
             divp(mi)=0.0
             dive(mi)=0.0
            endif
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
c gaz 110123
            if(idoff.ne.-1) then
             rovf(mi)=rov
            endif

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
c eflow should have fixed temperature (no cprd)
            if(kq.lt.0) then
               eskd=eflow(mi)
               edif=tl-eskd
               permsd=wellim(mi)
               qh(mi)=permsd*edif
c gaz   081317            
               deqh(mi)=permsd
            endif
            if(kq.ge.0.and.qflux(mi).eq.0.0) deqh(mi)=0.
            dtpae(mi)=1.
c gaz 082917            
c            denrd=denr(mi)*cprd
            sto1(mi)=0.0
            sto1(mi+neq)=ur
            deef(mi)=dur_dt*dtin
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
