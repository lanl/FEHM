      module comrxni
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
!D1 solute module.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 09-13-93     B. Robinson    000??    Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comrxni.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:59:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:24   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:20   pvcs
!D2 original version in process of being certified
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
!D4   Identifier    Type     Description
!D4 
!D4   conc_min      REAL*8   Minimum concentration that values less than this
!D4                            get set to
!D4   gas_constant  REAL*8   Universal gas constant
!D4   max_rxns      INT      Maximum number of reactions allowed
!D4   max_species   INT      Maximum number of species allowed
!D4   mw_air        REAL*8   Molecular weight of air (kg/mol)
!D4   mw_water      REAL*8   Molecular weight of water (kg/mol)
!D4   rtol          REAL*8   Minimum value of certain parameters used in
!D4                            calculations
!D4   temp_conv     REAL*8   Conversion from !elcius to Kelvin temperature
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block     Description
!D4
!D4   ***** COMMON Block readrxn variables *****
!D4   ar_for          REAL*8   readrxn   Pre-exponential factors for forward
!D4                                        reactions
!D4   ar_rev          REAL*8   readrxn   Pre-exponential factors for reverse
!D4                                        reactions
!D4   ea_for          REAL*8   readrxn   Activation energies for forward
!D4                                        reactions
!D4   ea_rev          REAL*8   readrxn   Activation energies for reverse
!D4                                        reactions
!D4   fl_mult         REAL*8   readrxn   Flag to denote if fluid species
!D4                                        participates in reaction
!D4   nrxns           INT      readrxn   Number of reactions
!D4   rate_power      REAL*8   readrxn   Exponents for each species in the
!D4                                        reactions
!D4   sb_mult         REAL*8   readrxn   Flag to denote if sorbed-phase
!D4                                        species participates in reaction
!D4   stoic           REAL*8   readrxn   Stoiciometric coefficients for
!D4                                        reactions
!D4 
!D4   ***** COMMON Block stor pointers and associated variables *****
!D4    ipavgmolwt     POINTER  stor      Pointer for variable avgmolwt
!D4    ipdsccl        POINTER  stor      Pointer for variable dsccl
!D4    ipdsccv        POINTER  stor      Pointer for variable dsccv
!D4    ipscl          POINTER  stor      Pointer for variable scl
!D4    ipscv          POINTER  stor      Pointer for variable scv
!D4 
!D4    avgmolwt       REAL*8   stor      Average molecular weight of vapor
!D4    dsccl          REAL*8   stor      Derivative of sorbed phase
!D4                                        concentration with respect to
!D4                                        concentration
!D4    dsccv          REAL*8   stor      Derivative of sorbed phase
!D4                                        concentration with respect to
!D4                                        concentration
!D4    scl            REAL*8   stor      Sorbed phase concentration from
!D4                                        liquid to rock
!D4    scv            REAL*8   stor      Sorbed phase concentration from vapor
!D4                                        to rock
!D4 
!D4   ***** COMMON Block setup variables *****
!D4   forward_solids  INT      setup     Number of solid species in each
!D4                                        forward reaction
!D4   forward_species INT      setup     Number of species in each forward
!D4                                        reaction
!D4   rxn_array       INT      setup     Species number for a given reaction
!D4                                        and species number in reaction
!D4   rxn_interval    INT      setup     Parameter for determining when to
!D4                                        perform full iteration of
!D4                                        concentrations 
!D4   solid_flag      INT      setup     Array no longer used (can be omitted)
!D4   solid_rxns      INT      setup     Number of reactions involving a solid
!D4                                        that each species participates in
!D4   species_array   INT      setup     Reaction number for each species in
!D4                                        each reaction the current species
!D4                                        is involved in
!D4   total_rxns      INT      setup     Total number of reactions each
!D4                                        species participates in
!D4   total_solids    INT      setup     Number of solid species in each
!D4                                        reaction
!D4   total_species   INT      setup     Total number of species in each
!D4                                        reaction 
!D4   
!D4   ***** COMMON Block henry variables *****
!D4   a_henry         REAL*8   henry     Pre-exponential term in Henry's Law
!D4                                        coefficient
!D4   density_flag    INT      henry     Flag denoting the fluid density term
!D4                                        to use in the kinetic expression
!D4   dh_henry        REAL*8   henry     Exponential term in Henry's Law
!D4                                        coefficient
!D4   h_mult          REAL*8   henry     Flag denoting whether reaction takes
!D4                                        place in liquid or vapor for Henry's
!D4                                        Law species
!D4   rxn_flag        INT      henry     Flag denoting whether to call
!D4                                        chemical reaction routine
!D4   
!D4   ***** COMMON Block dispersion pointers and associated variables *****
!D4   ipdisplx        POINTER  dispersion Pointer for variable displx
!D4   ipdisply        POINTER  dispersion Pointer for variable disply
!D4   ipdisplz        POINTER  dispersion Pointer for variable displz
!D4   ipdispvx        POINTER  dispersion Pointer for variable dispvx
!D4   ipdispvy        POINTER  dispersion Pointer for variable dispvy
!D4   ipdispvz        POINTER  dispersion Pointer for variable dispvz
!D4   displx          REAL*8   dispersion Dispersion coefficient in 
!D4                                         x-direction in liquid
!D4   disply          REAL*8   dispersion Dispersion coefficient in 
!D4                                         y-direction in liquid
!D4   displz          REAL*8   dispersion Dispersion coefficient in 
!D4                                         z-direction in liquid
!D4   dispvx          REAL*8   dispersion Dispersion coefficient in 
!D4                                         x-direction in vapor
!D4   dispvy          REAL*8   dispersion Dispersion coefficient in 
!D4                                         y-direction in vapor
!D4   dispvz          REAL*8   dispersion Dispersion coefficient in 
!D4                                         z-direction in vapor
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

      real*8    conc_min
      parameter (conc_min = 1.e-20)
      real*8    gas_const
      parameter (gas_const = 8.314)
      real*8    strac_max

!      integer   max_rxns
!      parameter (max_rxns = 10)
! Above too!
! The following should remain but conflicts because in 
! comcouple too!!!!!!!!!!!!!!!!!!
!      integer   max_species
!      parameter (max_species = 10)

      real*8    mw_air 
      parameter (mw_air = 28.8e-3)
      real*8    mw_water
      parameter (mw_water = 18.e-3)
      real*8    rtol
      parameter (rtol = 1.e-7)
      real*8    temp_conv
      parameter (temp_conv = 273.16)

!     ***** Variables in COMMON Block readrxn *****
c      integer nrxns
c      real*8, allocatable ::  ar_for(:) 
c      real*8, allocatable ::  ar_rev(:)
c      real*8, allocatable ::  ea_for(:)
c      real*8, allocatable ::  ea_rev(:)
c      real*8, allocatable ::  fl_mult(:, :)
c      real*8, allocatable ::  rate_power(:, :)
c      real*8, allocatable ::  sb_mult(:, :)
c      real*8, allocatable ::  stoic(:, :)

!     ***** Pointers in COMMON Block stor *****
      real*8, allocatable ::  avgmolwt(:)
      real*8, allocatable ::  dsccl(:)
      real*8, allocatable ::  dsccv(:)
      real*8, allocatable ::  scl(:)
      real*8, allocatable ::  scv(:)

!     ***** Variables in COMMON Block setup *****
c      integer, allocatable ::  forward_solids(10)
c      integer, allocatable ::  forward_species(10)
c      integer, allocatable ::  rxn_array(10, 10)
c      integer rxn_interval
c      integer, allocatable ::  solid_flag(10)
c      integer, allocatable ::  solid_rxns(10)
c      integer, allocatable ::  species_array(10, 10)
c      integer, allocatable ::  total_rxns(10)
c      integer, allocatable ::  total_solids(10)
c      integer, allocatable ::  total_species(10)

!     ***** Variables in COMMON Block henry *****
c zvd 02-Sep-08 Make henry model variables allocatable
      integer density_flag(10)
      integer rxn_flag
      real*8, allocatable ::  a_henry(:)
      real*8, allocatable ::  dh_henry(:)
c      real*8, allocatable ::  h_mult(:,:)

!     ***** Pointers in COMMON Block dispersion *****
      real*8, allocatable ::  displx(:)
      real*8, allocatable ::  disply(:)
      real*8, allocatable ::  displz(:)
      real*8, allocatable ::  dispvx(:)
      real*8, allocatable ::  dispvy(:)
      real*8, allocatable ::  dispvz(:)

c zvd 07-Mar-12 Variables for concentration-dependent density
      integer :: cden_flag = 0, cden_sp = 0
      real, allocatable :: mw_speci(:)
      character(20) :: cden_spnam

      end module comrxni

