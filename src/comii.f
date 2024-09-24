      module comii
!    comii
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
!D2 $Log:   /pvcs.config/fehm90/src/comii.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:57:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:16   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:12   pvcs
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
!D4   ***** COMMON Block coeff pointers and associated variables *****
!D4   ipcel           POINTER  coeff  Pointer for variable cel
!D4   ipcev           POINTER  coeff  Pointer for variable cev
!D4   ipcrl           POINTER  coeff  Pointer for variable crl
!D4   ipcrv           POINTER  coeff  Pointer for variable crv
!D4   ipcvl           POINTER  coeff  Pointer for variable cvl
!D4   ipcvv           POINTER  coeff  Pointer for variable cvv
!D4
!D4   cel             REAL*8   coeff  Polynomial coefficients for liquid water
!D4                                     enthalpy equations
!D4   cev             REAL*8   coeff  Polynomial coefficients for vapor water
!D4                                     enthalpy equations
!D4   crl             REAL*8   coeff  Polynomial coefficients for liquid water
!D4                                     density equations
!D4   crv             REAL*8   coeff  Polynomial coefficients for vapor water
!D4                                     density equations
!D4   cvl             REAL*8   coeff  Polynomial coefficients for liquid water
!D4                                     viscosity equations
!D4   cvv             REAL*8   coeff  Polynomial coefficients for vapor water
!D4                                     viscosity equations
!D4   
!D4   ***** COMMON Block coeffr variables *****
!D4   psa0            REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   psb0            REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   psta1           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   psta2           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   psta3           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   psta4           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   pstb1           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   pstb2           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   pstb3           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   pstb4           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     pressure equation
!D4   tsa0            REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tsb0            REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspa1           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspa2           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspa3           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspa4           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspb1           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspb2           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspb3           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   tspb4           REAL*8   coeffr Polynomial coefficient for saturation
!D4                                     temperature equation
!D4   
!D4   ***** COMMON Block coeff1 pointers and associated variables *****
!D4   ippmax          POINTER  coeff1 Pointer for variable pmax
!D4   ippmin          POINTER  coeff1 Pointer for variable pmin
!D4   iptmax          POINTER  coeff1 Pointer for variable tmax
!D4   iptmin          POINTER  coeff1 Pointer for variable tmin

!D4   pmax            REAL*8   coeff1 Maximum pressure allowed for each
!D4                                     coefficient set
!D4   pmin            REAL*8   coeff1 Minimum pressure allowed for each
!D4                                     coefficient set
!D4   tmax            REAL*8   coeff1 Maximum temperature allowed for each
!D4                                     coefficient set
!D4   tmin            REAL*8   coeff1 Minimum temperature allowed for each
!D4                                     coefficient set
!D4   
!D4   ***** COMMON Block coeff2 variables *****
!D4   ew1             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew2             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew3             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew4             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew5             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew6             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew7             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew8             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew9             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew10            REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ew11            REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev1             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev2             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev3             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev4             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev5             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev6             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev7             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev8             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev9             REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev10            REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
!D4   ev11            REAL*8   coeff2 Coefficient used in simplifying
!D4                                     thermodynamics relations
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

!     ***** Pointers in COMMON Block coeff *****
      real*8, allocatable ::  cel(:,:) 
      real*8, allocatable ::  crl(:,:) 
      real*8, allocatable ::  cev(:,:) 
      real*8, allocatable ::  crv(:,:)
      real*8, allocatable ::  cvl(:,:) 
      real*8, allocatable ::  cvv(:,:)
 
!     ***** Variables in COMMON Block coeffr *****
      real*8  psa0, psb0, psta1, psta2, psta3, psta4, pstb1, pstb2
      real*8  pstb3, pstb4, tsa0, tsb0, tspa1, tspa2, tspa3, tspa4
      real*8  tspb1, tspb2, tspb3, tspb4

!     ***** Pointers in COMMON Block coeff2 *****
      real*8, allocatable ::  pmax(:) 
      real*8, allocatable ::  pmin(:) 
      real*8, allocatable ::  tmax(:) 
      real*8, allocatable ::  tmin(:)
! gaz 070720 add (p,t)info for air_tables      
      real*8, allocatable ::  pmax_air_tabl(:) 
      real*8, allocatable ::  pmin_air_tabl(:) 
      real*8, allocatable ::  tmax_air_tabl(:) 
      real*8, allocatable ::  tmin_air_tabl(:)
! gaz 0081421 add (p,t)info for co2_tables      
      real*8, allocatable ::  pmax_co2wh_tabl(:) 
      real*8, allocatable ::  pmin_co2wh_tabl(:) 
      real*8, allocatable ::  tmax_co2wh_tabl(:) 
      real*8, allocatable ::  tmin_co2wh_tabl(:)      

!     ***** Variables in COMMON Block coeff2 *****
      real*8  ew1, ew2, ew3, ew4, ew5, ew6, ew7, ew8, ew9, ew10, ew11
      real*8  ev1, ev2, ev3, ev4, ev5, ev6, ev7, ev8, ev9, ev10, ev11
      real*8 dennapl, viscnapl
! gaz 082621 fluids for ngas
      character*9 fluid(3)

      end module comii
