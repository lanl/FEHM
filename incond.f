      subroutine incond
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Read thermal conductivity data.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 21-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/incond.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jun 12 15:21:02 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2
CD2    Rev 1.4   Mon Jun 03 11:17:50 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.3   Fri May 31 10:57:36 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.2   Tue Jan 30 09:48:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:03:20   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:48   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3   
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   inpt            INT      faai   Unit number for input file
CD4   ipthx           POINTER  fdd    Pointer to array variable thx
CD4   ipthy           POINTER  fdd    Pointer to array variable thy
CD4   ipthz           POINTER  fdd    Pointer to array variable thz
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(10)   LOGICAL  macro  Flag denoting if macro cond has been read 
CD4   thx             REAL*8   fdd    Thermal conductivity x-direction
CD4   thy             REAL*8   fdd    Thermal conductivity y-direction
CD4   thz             REAL*8   fdd    Thermal conductivity z-direction
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   initdata                 Read data values and set parameter values at
CD4                              given nodes
CD4 
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   i               INT      Loop index
CD5   igroup          INT      Current group number in this macro
CD5   ireturn         INT      Returned error flag from input subroutine
CD5   itype           INT      Array of variable types being input
CD5   default         REAL*8   Array of default values for input arrays
CD5   macro           CHAR     Current macro being read
CD5   narrays         INT      Number of arrays being read in
CD5   pointer         INT      Integer array of pointer values for
CD5                              variables being read
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN incond
CPS 
CPS   set arrays used, pointers, data type and default values
CPS   call initdata to read data values and set conductivity values at . . .
CPS   . . . given nodes
CPS   
CPS   FOR each node
CPS       set conductivity value to maximum of value read or zero_t
CPS   END FOR
CPS   
CPS   set macroread to true
CPS
CPS END incond
CPS
C***********************************************************************

      use comdi
      use comdti
      use comai
      use comki
      implicit none
      
      integer i
      
      real*8, allocatable :: dumx(:)
      real*8, allocatable :: dumy(:)
      real*8, allocatable :: dumz(:)
      macro = 'cond'
      
c**** read thermal conductivity data ****
      
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = zero_t
      default(2) = zero_t
      default(3) = zero_t
      allocate(dumx(n0),dumy(n0),dumz(n0))
      igroup = 1
      if(ico2.ge.0.or.ice.ne.0) then
         dumx = thx
         dumy = thy
         dumz = thz
      else
         dumx = thx(1)
         dumy = thy(1)
         dumz = thz(1)
      endif

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     &     default, macroread(10), macro, igroup, ireturn,
     &     r8_1=dumx(1:n0),r8_2=dumy(1:n0),r8_3=dumz(1:n0)) 

      if(ico2.ge.0.or.ice.ne.0) then
         do i = 1, n0
            thx(i) = max (zero_t,dumx(i))
            thy(i) = max (zero_t,dumy(i))
            thz(i) = max (zero_t,dumz(i))
         end do
      endif

      macroread(10) = .TRUE.

      deallocate(dumx,dumy,dumz)
      end
