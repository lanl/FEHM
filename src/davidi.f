      module davidi
!     davidi           
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
!D1 Global include file for array variables and pointers (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 19-OCT-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/davidi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:30   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:54   pvcs
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
!D4   ***** COMMON Block david1 variables *****
!D4   iback           INT      david1 LU factorization save parameter 
!D4   icoupl          INT      david1 Number of SOR iterations 
!D4   irdof           INT      david1 Reduced degree of freedom model used 
!D4   islord          INT      david1 Parameter used in the reduced degree of
!D4                                     freedom model 
!D4   itest           INT      david1 Parameter used in the reduced degree of
!D4                                     freedom model 
!D4
!D4   ***** COMMON Block david2 pointers and associated variables *****
!D4   ipnb            POINTER  david2 Pointer to variable array nb
!D4   ipnmat          POINTER  david2 Pointer to variable array nmat
!D4   ipnmatb         POINTER  david2 Pointer to variable array nmatb
!D4   ipnrhs          POINTER  david2 Pointer to variable array nrhs
!D4   ipnrhsb         POINTER  david2 Pointer to variable array nrhsb
!D4   
!D4   nb              INT      david2 Array used in the reduced degree of
!D4                                      freedom method 
!D4   nmat            INT      david2 Array used in the reduced degree of
!D4                                      freedom method 
!D4   nmatb           INT      david2 Array used in the reduced degree of
!D4                                      freedom method 
!D4   nrhs            INT      david2 Array used in the reduced degree of
!D4                                      freedom method 
!D4   nrhsb           INT      david2 Array used in the reduced degree of
!D4                                      freedom method 
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

      integer iback, icoupl, irdof, islord, itest

!     ***** Pointers in COMMON Block david2 *****
      integer, allocatable :: nb(:)
      integer, allocatable :: nmat(:)
      integer, allocatable :: nmatb(:)
      integer, allocatable :: nrhs(:)
      integer, allocatable :: nrhsb(:)
      integer, allocatable :: nvar_select(:)
      integer max_variables
      parameter(max_variables = 12)

      end module davidi
