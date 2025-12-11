      module comsplitts
!           comsplitts
!***********************************************************************
!  Copyright, 1995, 2004,  The  Regents of the University of California.
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
!D1 the option for split timesteps for phase changes
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 Initial implementation: 24-OCT-95, Programmer: B. Robinson 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comsplitts.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:30:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
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

      integer isplitts, ivf_sv, ipr_tets
      real*8 dtot_split, dtot_next, explicit_factor
      logical sv_hex_tet 

!     ***** Variables in COMMON Block splitts_pointers *****
      real*8, allocatable :: bp_last(:)

      end module comsplitts 

