      module comhi
!    comhi
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
!D2 $Log:   /pvcs.config/fehm90/src/comhi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:57:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:46 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:14   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:10   pvcs
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
!D4   ***** COMMON Block dualp pointers and associated variables *****
!D4   ipa21epf        POINTER  dualp  Pointer for variable a21epf
!D4   ipa21mef        POINTER  dualp  Pointer for variable a21mef
!D4   ipa21mpf        POINTER  dualp  Pointer for variable a21mpf
!D4   ipa32eef        POINTER  dualp  Pointer for variable a32eef
!D4   ipa32epf        POINTER  dualp  Pointer for variable a32epf
!D4   ipa32mef        POINTER  dualp  Pointer for variable a32mef
!D4   ipa32mpf        POINTER  dualp  Pointer for variable a32mpf
!D4   ipapuv1         POINTER  dualp  Pointer for variable apuv1
!D4   iprb2ef         POINTER  dualp  Pointer for variable rb2ef
!D4   iprb2mf         POINTER  dualp  Pointer for variable rb2mf
!D4   iprb3ef         POINTER  dualp  Pointer for variable rb3ef
!D4   iprb3mf         POINTER  dualp  Pointer for variable rb3mf
!D4   iptb11          POINTER  dualp  Pointer for variable tb11
!D4   iptb12          POINTER  dualp  Pointer for variable tb12
!D4   iptb21          POINTER  dualp  Pointer for variable tb21
!D4   iptb22          POINTER  dualp  Pointer for variable tb22
!D4   ipvolf1         POINTER  dualp  Pointer for variable volf1
!D4   ipvolf2         POINTER  dualp  Pointer for variable volf2
!D4   ipwb11          POINTER  dualp  Pointer for variable wb11
!D4   ipwb12          POINTER  dualp  Pointer for variable wb12
!D4   ipwb21          POINTER  dualp  Pointer for variable wb21
!D4   ipwb22          POINTER  dualp  Pointer for variable wb22
!D4   
!D4   a21epf          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a21mef          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a21mpf          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a32eef          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a32epf          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a32mef          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   a32mpf          REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   apuv1           REAL*8   dualp  Area per unit volume for the first
!D4                                     matrix layer
!D4   rb2ef           REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   rb2mf           REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   rb3ef           REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   rb3mf           REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   tb11            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   tb12            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   tb21            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   tb22            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   volf1           REAL*8   dualp  Volume fraction at each node for the
!D4                                     first matrix layer
!D4   volf2           REAL*8   dualp  Volume fraction at each node for the
!D4                                     second matrix layer
!D4   wb11            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   wb12            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   wb21            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
!D4   wb22            REAL*8   dualp  Array needed to store intermediate dual
!D4                                     porosity results
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

!     ***** Pointers in COMMON Block dualp *****
      real*8, allocatable ::  a21eef(:)
      real*8, allocatable ::  a21epf(:)
      real*8, allocatable ::  a21mef(:)
      real*8, allocatable ::  a21mpf(:)
      real*8, allocatable ::  a32eef(:)
      real*8, allocatable ::  a32epf(:)
      real*8, allocatable ::  a32mef(:)
      real*8, allocatable ::  a32mpf(:)
      real*8, allocatable ::  apuv1(:)
      real*8, allocatable ::  r3ef(:)
      real*8, allocatable ::  r3mf(:)
      real*8, allocatable ::  rb2ef(:)
      real*8, allocatable ::  rb2mf(:)
      real*8, allocatable ::  tb11(:)
      real*8, allocatable ::  tb12(:)
      real*8, allocatable ::  tb21(:)
      real*8, allocatable ::  tb22(:)
      real*8, allocatable ::  volf1(:)
      real*8, allocatable ::  volf2(:)
      real*8, allocatable ::  wb11(:)
      real*8, allocatable ::  wb12(:)
      real*8, allocatable ::  wb21(:)
      real*8, allocatable ::  wb22(:) 

      end module comhi
