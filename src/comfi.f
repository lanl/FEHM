      module comfi 
!    comfi
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Include file containing passed parameters and pointers related to
!D1 memory management.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 09-13-93     B. Robinson    00022   Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comfi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   01/28/95 14:07:10   llt
!D2 water balance equation was modified
!D2 
!D2    Rev 1.1   03/18/94 16:23:08   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:06   pvcs
!D2 original version in process of being certified
!D2 10/28/94 GAZ  added dqpc , the derivative of total source wrt air pressure
!D2
!**********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!**********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   None
!D4 
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block co2 pointers and associated variables *****
!D4   ipcnlf          POINTER  co2    Pointer for variable cnlf
!D4   ipcnvf          POINTER  co2    Pointer for variable cnvf
!D4   ipdcc           POINTER  co2    Pointer for variable dcc
!D4   ipdce           POINTER  co2    Pointer for variable dce
!D4   ipdclcf         POINTER  co2    Pointer for variable dclcf
!D4   ipdclef         POINTER  co2    Pointer for variable dclef
!D4   ipdclf          POINTER  co2    Pointer for variable dclf
!D4   ipdcp           POINTER  co2    Pointer for variable dcp
!D4   ipdcqc          POINTER  co2    Pointer for variable dcqc
!D4   ipdcqh          POINTER  co2    Pointer for variable dcqh
!D4   ipdcvcf         POINTER  co2    Pointer for variable dcvcf
!D4   ipdcvef         POINTER  co2    Pointer for variable dcvef
!D4   ipdcvf          POINTER  co2    Pointer for variable dcvf
!D4   ipdec           POINTER  co2    Pointer for variable dec
!D4   ipdelcf         POINTER  co2    Pointer for variable delcf
!D4   ipdenpch        POINTER  co2    Pointer for variable denpch
!D4   ipdenpci        POINTER  co2    Pointer for variable denpci
!D4   ipdenpcj        POINTER  co2    Pointer for variable denpcj
!D4   ipdeqc          POINTER  co2    Pointer for variable deqc
!D4   ipdevcf         POINTER  co2    Pointer for variable devcf
!D4   ipdglc          POINTER  co2    Pointer for variable dglc
!D4   ipdgvc          POINTER  co2    Pointer for variable dgvc
!D4   ipdilc          POINTER  co2    Pointer for variable dilc
!D4   ipdivc          POINTER  co2    Pointer for variable divc
!D4   ipdmc           POINTER  co2    Pointer for variable dmc
!D4   ipdqc           POINTER  co2    Pointer for variable dqc
!D4   ipdtpac         POINTER  co2    Pointer for variable dtpac
!D4   ipeskc          POINTER  co2    Pointer for variable eskc
!D4   ippci           POINTER  co2    Pointer for variable pci
!D4   ippcio          POINTER  co2    Pointer for variable pcio
!D4   ipqc            POINTER  co2    Pointer for variable qc
!D4   ipdqpc          POINTER  co2    Pointer for variable dqpc
!D4   ipsici          POINTER  co2    Pointer for variable sici
!D4   
!D4   cnlf            REAL*8   co2    Gas concentration in the liquid phase
!D4   cnvf            REAL*8   co2    !oncentration of gas in the vapor phase
!D4   dcc             REAL*8   co2    Derivative of gas accumulation term with
!D4                                     respect to gas
!D4   dce             REAL*8   co2    Derivative of gas accumulation term with
!D4                                     respect to energy at each node
!D4   dclcf           REAL*8   co2    Derivative of liquid concentration with
!D4                                     respect to gas
!D4   dclef           REAL*8   co2    Derivative of liquid concentration with
!D4                                     respect to energy
!D4   dclf            REAL*8   co2    Derivative of liquid concentration with
!D4                                     respect to pressure
!D4   dcp             REAL*8   co2    Derivative of gas accumulation term with
!D4                                     respect to pressure at each node
!D4   dcqc            REAL*8   co2    Derivative of gas source term with
!D4                                     respect to gas
!D4   dcqh            REAL*8   co2    Derivative of gas source term with
!D4                                     respect to energy
!D4   dcvcf           REAL*8   co2    Derivative of gas concentration with
!D4                                     respect to gas
!D4   dcvef           REAL*8   co2    Derivative of gas concentration with
!D4                                     respect to energy
!D4   dcvf            REAL*8   co2    Derivative of gas concentration with
!D4                                     respect to pressure
!D4   dec             REAL*8   co2    Derivative of energy accumulation term
!D4                                     with respect to gas at each node
!D4   delcf           REAL*8   co2    Derivative of liquid energy with respect
!D4                                     to gas
!D4   denpch          REAL*8   co2    Last time step value of the mass
!D4                                     accumulation term at each node
!D4   denpci          REAL*8   co2    Gas accumulation term at each node
!D4   denpcj          REAL*8   co2    Last time step accumulation term of gas
!D4                                     equation
!D4   deqc            REAL*8   co2    Derivative of energy source term with
!D4                                     respect to gas
!D4   devcf           REAL*8   co2    Derivative of energy with respect to gas
!D4   dglc            REAL*8   co2    Derivative of liquid gravity term with
!D4                                     respect to gas
!D4   dgvc            REAL*8   co2    Derivative of vapor gravity term with
!D4                                     respect to gas
!D4   dilc            REAL*8   co2    Derivative of liquid transmissibility
!D4                                     with respect to gas
!D4   divc            REAL*8   co2    Derivative of vapor transmissibility
!D4                                     with respect to gas
!D4   dmc             REAL*8   co2    Derivative of mass accumulation term
!D4                                     with respect to gas at each node
!D4   dqc             REAL*8   co2    Derivative of mass source term with
!D4                                     respect to gas
!D4   dtpac           REAL*8   co2    Derivative of temperature with respect
!D4                                     to gas
!D4   eskc            REAL*8   co2    Source term for gas equation
!D4   pci             REAL*8   co2    Gas pressure
!D4   pcio            REAL*8   co2    Last time step gas pressure
!D4   qc              REAL*8   co2    Source term for the gas equation
!D4   dqpc            REAL*8   co2    derivative of total source wrt air pressure
!D4   sici            REAL*8   co2    Array in ice solution
!D4   
!D4   ***** COMMON Block co2r variables *****
!D4   acner           REAL*8   co2r   Total mass accumulation at a time step
!D4   amc             REAL*8   co2r   Total mass accumulation at a time step
!D4   difc            REAL*8   co2r   Noncondensible gas mass balance error
!D4   qtc             REAL*8   co2r   Total mass injected in source term
!D4   qtotc           REAL*8   co2r   Total mass injected in source term
!D4   
!**********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!**********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!**********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!**********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!**********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 None
!D9
!**********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!**********************************************************************
!PS
!PS PSEUDOCODE
!PS 
!PS None
!PS
!**********************************************************************

      real*8, allocatable ::  cnlf(:)
      real*8, allocatable ::  cnlof(:)
      real*8, allocatable ::  cnvf(:) 
      real*8, allocatable ::  dcc(:) 
      real*8, allocatable ::  dce(:) 
      real*8, allocatable ::  dclcf(:) 
      real*8, allocatable ::  dclef(:) 
      real*8, allocatable ::  dclf(:)
      real*8, allocatable ::  dcp(:) 
      real*8, allocatable ::  dcqc(:) 
      real*8, allocatable ::  dcqh(:) 
      real*8, allocatable ::  dcvcf(:) 
      real*8, allocatable ::  dcvef(:)
      real*8, allocatable ::  dcvf(:) 
      real*8, allocatable ::  dec(:) 
      real*8, allocatable ::  delcf(:) 
      real*8, allocatable ::  denpch(:)
      real*8, allocatable ::  denpci(:) 
      real*8, allocatable ::  denpcj(:) 
      real*8, allocatable ::  deqc(:) 
      real*8, allocatable ::  devcf(:) 
      real*8, allocatable ::  dglc(:)
      real*8, allocatable ::  dgvc(:) 
      real*8, allocatable ::  dilc(:) 
      real*8, allocatable ::  divc(:) 
      real*8, allocatable ::  dmc(:) 
      real*8, allocatable ::  dqc(:)
      real*8, allocatable ::  dtpac(:)
      real*8, allocatable ::  eskc(:) 
      real*8, allocatable ::  pci(:) 
      real*8, allocatable ::  pcio(:)
      real*8, allocatable ::  qc(:) 
      real*8, allocatable ::  qng(:)
      real*8, allocatable ::  dqpc(:)  
      real*8, allocatable ::  sici(:) 
      real*8, allocatable ::  flux_ts(:) 
c gaz 100920 added last iteration variable values (air-water-heat)   
      real*8, allocatable :: phi_last_iter(:)
      real*8, allocatable :: t_last_iter(:) 
      real*8, allocatable :: pci_last_iter(:)
      real*8, allocatable :: s_last_iter(:) 

!     ***** Variables in COMMON Block co2r *****
      real*8 acner, amc, difc, qtc, qtotc, qtotin
c gaz 101020 alpha0 (Henry's law constant made global)
c gaz 121223 added Henry's law constant, isothermal, two phase
      real*8, allocatable :: frac_gas_iso(:)
      real*8 alpha0, alpha_air, alpha_h2, alpha_meth, alpha_co2
      integer imod_sol
c gaz 051722 added some variables for mass calculations
      real*8, allocatable :: xntot(:)
      real*8, allocatable :: zntotf(:)
      real*8, allocatable :: zntot(:)
c gaz 080722 more global mass variables   
      real*8, allocatable :: mass_h2o(:)
      real*8, allocatable :: mass_ngas(:)
      real*8, allocatable :: energy_tot(:)
      real*8, allocatable :: mass0_h2o(:)
      real*8, allocatable :: mass0_ngas(:)
      real*8, allocatable :: energy0_tot(:)
c gaz 090222  
      character*4, allocatable :: dum_calc_eos(:)
c gaz 090622
      real*8, allocatable :: dtot_test(:)
c gaz 090722 make reset variables global    
      real*8 phi_reset, pci_reset, t_reset, s_reset, dtot_min
c gaz 120223 added idiff_iso
      integer ieos_reset, iphase_chng, idiff_iso
c gaz 040723, 050123  
      character*22, allocatable :: descrip(:,:)
      real*8 phi_max, phi_min
      integer node_awh, inr_count, idelp, imass_unit, nawh, imass_chk
      integer imout
      integer, allocatable ::  inode_chng(:)
      integer, allocatable ::  nawh_nodes(:)
      integer, allocatable ::  ieos0(:)
c gaz 030223, 120223
      real*8, allocatable ::  sk_awh(:), dsk_awhdp(:), q_gas(:)
      real*8 q_gas_tot
      real*8 bpt(3),bpt0(3),a33t(3,3),a33ti(3,3),del(3)
c gaz032623
      logical  ts_chng
      logical, allocatable :: read_mfrac_iso(:) 
      end module comfi
 