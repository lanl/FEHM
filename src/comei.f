      module comei
!     comei
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global include file for array variables and !pointers (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 30-SEP-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comei.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   03/28/94 16:41:08   robinson
!D2 Removed unneeded array.
!D2 
!D2    Rev 1.1   03/18/94 16:23:06   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:02   pvcs
!D2 original version in process of being certified
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
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
!D4   ***** COMMON Block fee !pointer and associated variable *****
!D4   ipa             !pointer  fee    !pointer to variable array a     
!D4
!D4   a               REAL*8          array containing the jacobian matrix 
!D4
!D4   ***** COMMON Block fff !pointers and associated variables *****
!D4   ipb             !pointer  fff    !pointer to variable array b     
!D4   ipbcoef         !pointer  fff    !pointer to variable array bcoef     
!D4
!D4   b               REAL*8          array containing the incomplete lu 
!D4                                     decomposition of the jacobian matrix
!D4   bcoef           REAL*8          
!D4
!D4   ***** COMMON Block fhh !pointers and associated variables *****
!D4   ipiirb          !pointer  fhh    !pointer to variable array iirb  
!D4   ipirb           !pointer  fhh    !pointer to variable array irb   
!D4   ipnopt          !pointer  fhh    !pointer to variable array nopt  
!D4   ipnpvt          !pointer  fhh    !pointer to variable array npvt  
!D4   ippiv           !pointer  fhh    !pointer to variable array piv   
!D4
!D4   iirb            INT      fhh    Inverse of irb 
!D4   irb             INT      fhh    Array containing the reordered node  
!D4                                     numbers 
!D4   nopt            INT      fhh    Array indicating active variables 
!D4   npvt            INT      fhh    Pivot information for the lu  
!D4                                     decomposition matrix 
!D4   piv             REAL*8   fhh    The pivot elements of the lu  
!D4                                     decomposition 
!D4
!D4   ***** COMMON Block fii !pointers and associated variables *****
!D4   ipc             !pointer  fii    !pointer to variable array c
!D4   ipg             !pointer  fii    !pointer to variable array g 
!D4   iph             !pointer  fii    !pointer to variable array h
!D4   ipss            !pointer  fii    !pointer to variable array ss 
!D4   ipy             !pointer  fii    !pointer to variable array y
!D4
!D4   c               REAL*8   fii    ?          
!D4   g               REAL*8   fii    ?          
!D4   h               REAL*8   fii    ?          
!D4   ss              REAL*8   fii    ?           
!D4   y               REAL*8   fii    ?          
!D4
!D4 Global Subprograms
!D4
!D4   None
!D4
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 N/A
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS None
!PS
!***********************************************************************

!     ***** !pointer in COMMON Block fee *****
      real*8, allocatable ::  a(:)
 
!     ***** !pointers in COMMON Block fff *****
      real*8, allocatable ::  bcoef(:,:)

!     ***** !pointers in COMMON Block fhh *****
      integer, allocatable :: iirb(:)
      integer, allocatable :: irb(:)
      integer, allocatable :: nopt(:)
      integer, allocatable :: npvt(:)

!     ***** !pointers in COMMON Block fii *****
      real*8, allocatable ::  c(:)
      real*8, allocatable ::  g(:)
      real*8, allocatable ::  h(:,:)
      real*8, allocatable ::  ss(:)
      real*8, allocatable ::  y(:)

!     ***** active-passive variables *****

      real*8, allocatable ::  a_active(:) 
      real*8, allocatable ::  phi_base(:) 
      real*8, allocatable ::  t_base(:) 
      real*8, allocatable ::  phi_base_grad(:)  
      real*8, allocatable ::  t_base_grad(:) 
      real*8, allocatable ::  time_phi_base(:) 
      real*8, allocatable ::  time_t_base(:) 
      real*8, allocatable ::  bp_base(:)
      integer, allocatable ::  nelm_active(:) 
      integer, allocatable ::  node_active(:) 
      integer, allocatable ::  ncount(:) 
      integer, allocatable ::  nmat_active(:)
      integer, allocatable ::  nrhs_active(:)
      integer n_active, av_stop, av_write
      integer iav_chkp,iav_chkt,iav_chkbp,ij_change
      integer iacfile
      real*8  day_tol
      real*8 p_act_tol, t_act_tol, h_act_tol, bp_act_tol, resid_tol
      parameter(day_tol = 1.d-12)

      end module comei
