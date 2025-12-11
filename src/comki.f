      module comki
!    comki
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
!D2 12-21-93     B. Robinson    00022   Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comki.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:57:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   03/18/94 16:23:20   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:16   pvcs
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
!D4   Identifier      Type     Description
!D4   
!D4   max_inputs      int      Maximum number of variables being read
!D4                               on an input line
!D4 
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4   None
!D4   
!**********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5   Identifier      Type     Description
!D5   
!D5   default         REAL*8   Array of default values for input arrays
!D5   igroup          INT      Current group number in this macro
!D5   ireturn         INT      Returned error flag from input subroutine
!D5   itype           INT      Array of variable types being input
!D5   macro           CHAR     Current macro being read
!D5   narrays         INT      Number of arrays being read in
!D5   pointer         INT      Integer array of pointer values for
!D5                              variables being read
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

      integer narrays
      integer max_inputs
      parameter(max_inputs = 10)
      integer itype(max_inputs)
      real*8 default(max_inputs)
      character*4 macro
      integer igroup
      integer ireturn

      interface

         subroutine initdata2(in_number, out_number, npoints,
     2        narrays, itype, default, readflag, macro, igroup,ireturn,
     3        r8_1,r8_2,r8_3,r8_4,r8_5,i4_1,i4_2,i4_3,i4_4,i4_5)
         integer in_number,out_number,npoints,narrays,ireturn
         integer itype(*)
         integer igroup
         integer max_arrays
         parameter(max_arrays=10)
         real*8 default(max_arrays),values(max_arrays)
         logical readflag
         character*4 macro
         
         real*8, optional :: r8_1(:)
         real*8, optional :: r8_2(:)
         real*8, optional :: r8_3(:)
         real*8, optional :: r8_4(:)
         real*8, optional :: r8_5(:)
         integer, optional :: i4_1(:)
         integer, optional :: i4_2(:)
         integer, optional :: i4_3(:)
         integer, optional :: i4_4(:)
         integer, optional :: i4_5(:)
         end subroutine

      end interface

      end module comki
