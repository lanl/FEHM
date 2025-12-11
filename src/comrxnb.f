      module comrxnb
!     comrxnb
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
!D2 09-13-93     B. Robinson    00022   Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comrxnb.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:59:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:22   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:18   pvcs
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
!D4                                 COMMON
!D4   Identifier           Type     Block    Description
!D4   
!D4   ***** COMMON Block comrxnb pointers and associated variables *****
!D4   ipstored_derivative  POINTER  comrxnb  Pointer for array stored_derivative
!D4   ipstored_residual    POINTER  comrxnb  Pointer for array stored_residual
!D4   stored_derivatives   REAL*8   comrxnb  Array of derivatives for each
!D4                                            solute at each node
!D4   stored_residual      REAL*8   comrxnb  Array of residuals for each
!D4                                            solute at each node
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

      real*8, allocatable ::  stored_derivative(:)
      real*8, allocatable ::  stored_residual(:)

      end module comrxnb

