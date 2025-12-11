      subroutine inperm
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
CD1 Read permeability data.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 22-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/inperm.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jun 12 15:21:06 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.4   Mon Jun 03 11:18:02 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.3   Fri May 31 10:56:16 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.2   Tue Jan 30 13:11:56 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:03:42   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:12   pvcs
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
CD4   ippnx           POINTER  fdd    Pointer to array variable pnx
CD4   ippny           POINTER  fdd    Pointer to array variable pny
CD4   ippnz           POINTER  fdd    Pointer to array variable pnz
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(14)   LOGICAL  macro  Flag denoting if macro perm has been read 
CD4   pnx             REAL*8   fdd    Permeability in the x-direction, liquid
CD4   pny             REAL*8   fdd    Permeability in the y-direction liquid
CD4   pnz             REAL*8   fdd    Permeability in the z-direction, liquid
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
CPS BEGIN inperm
CPS 
CPS   set arrays used, pointers, data type and default values
CPS   call initdata to read data values and set permeability values at . . .
CPS   . . . given nodes
CPS   
CPS   FOR each node
CPS       set permeability value to maximum of value read or zero_t
CPS   END FOR
CPS   
CPS   set macroread to true
CPS
CPS END inperm
CPS
C***********************************************************************

      use comdi
      use comdti
      use comai
      use comki
      implicit none

      integer i

      macro = 'perm'
      read(inpt,'(a8)') wdd1(1:8)
      if(wdd1(1:8).eq.'modflow') then
c**** read hydraulic conductivity(and storitivity from modflow file
       call structured(2)
      else
      backspace inpt
c**** read permeability data ****
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = zero_t
      default(2) = zero_t
      default(3) = zero_t
      igroup = 1

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(14), macro, igroup, ireturn,
     2     r8_1=pnx(1:n0),r8_2=pny(1:n0),r8_3=pnz(1:n0)) 

      do i=1,n0
         if(pnx(i).ge.0.0d00) then
            pnx(i) = max (zero_t, pnx(i))
         else
            pnx(i) = 1.0d00/10.0d00**(abs(pnx(i)))
         endif
         if(pny(i).ge.0.0d00) then
            pny(i) = max (zero_t, pny(i))
         else
            pny(i) = 1.0d00/10.0d00**(abs(pny(i)))
         endif
         if(pnz(i).ge.0.0d00) then
            pnz(i) = max (zero_t, pnz(i))
         else
            pnz(i) = 1.0d00/10.0d00**(abs(pnz(i)))
         endif
      end do

      macroread(14) = .TRUE.
   
      endif

      end
