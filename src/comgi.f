      module comgi
!    comgi
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
!D2 $Log:   /pvcs.config/fehm90/src/comgi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:10   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:08   pvcs
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
!D4   ***** COMMON Block fgg pointers and associated variables *****
!D4   ipbp            POINTER  fgg    Pointer for variable bp
!D4   ipbpc           POINTER  fgg    Pointer for variable bpc
!D4   ipdfee          POINTER  fgg    Pointer for variable dfee
!D4   ipdfep          POINTER  fgg    Pointer for variable dfep
!D4   ipdfepc         POINTER  fgg    Pointer for variable dfepc
!D4   ipdfme          POINTER  fgg    Pointer for variable dfme
!D4   ipdfmp          POINTER  fgg    Pointer for variable dfmp
!D4   ipdfmpc         POINTER  fgg    Pointer for variable dfmpc
!D4   ipdfpce         POINTER  fgg    Pointer for variable dfpce
!D4   ipdfpcp         POINTER  fgg    Pointer for variable dfpcp
!D4   ipdfpcpc        POINTER  fgg    Pointer for variable dfpcpc
!D4   
!D4   bp              REAL*8   fgg    Array of Newton-Raphson residuals, after
!D4                                     solution an array of Newton-Raphson
!D4                                     corrections
!D4   bpc             REAL*8   fgg    Array of Newton-Raphson residuals for
!D4                                     gas, after solution, an array of
!D4                                     Newton-Raphson corrections
!D4   dfee            REAL*8   fgg    Derivative of energy accumulation term
!D4                                     with respect to energy for neighbor
!D4                                     nodes
!D4   dfep            REAL*8   fgg    Derivative of energy accumulation term
!D4                                     with respect to pressure for neighbor
!D4                                     nodes
!D4   dfepc           REAL*8   fgg    Derivative of energy accumulation term
!D4                                     with respect to gas for neighbor nodes
!D4   dfme            REAL*8   fgg    Derivative of mass accumulation term
!D4                                     with respect to energy for neighbor
!D4                                     nodes
!D4   dfmp            REAL*8   fgg    Derivative of mass accumulation term
!D4                                     with respect to pressure for neighbor
!D4                                     nodes
!D4   dfmpc           REAL*8   fgg    Derivative of mass accumulation term
!D4                                     with respect to gas for neighbor nodes
!D4   dfpce           REAL*8   fgg    Derivative of gas accumulation term with
!D4                                     respect to energy for neighbor nodes
!D4   dfpcp           REAL*8   fgg    Derivative of gas accumulation term with
!D4                                     respect to pressure for neighbor nodes
!D4   dfpcpc          REAL*8   fgg    Derivative of gas accumulation term with
!D4                                     respect to gas for neighbor nodes
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
!     ***** Pointers in COMMON Block fgg *****
      real*8, allocatable :: bp(:)
      real*8, allocatable :: bpc(:)
      real*8, allocatable :: dfee(:)
      real*8, allocatable :: dfep(:)
      real*8, allocatable :: dfepc(:)
      real*8, allocatable :: dfme(:)
      real*8, allocatable :: dfmp(:)
      real*8, allocatable :: dfmpc(:)
      real*8, allocatable :: dfpce(:)
      real*8, allocatable :: dfpcp(:)
      real*8, allocatable :: dfpcpc(:)

      end module comgi
