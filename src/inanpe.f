      subroutine inanpe
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
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
CD1
CD1 PURPOSE
CD1
CD1 Read anisotropic permeability data.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/inanpe.f_a  $
!D2
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
CD4   ipanxy          POINTER  fdd    Pointer to array variable anxy
CD4   ipanxz          POINTER  fdd    Pointer to array variable anxz
CD4   ipanyz          POINTER  fdd    Pointer to array variable anyz
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(14)   LOGICAL  macro  Flag denoting if macro perm has been read 
CD4   anxy            REAL*8   fdd    Permeability in the xy-direction, liquid
CD4   anxz             REAL*8   fdd    Permeability in the xz-direction liquid
CD4   anyz             REAL*8   fdd    Permeability in the yz-direction, liquid
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
      character*4 duma1
      macro = 'anpe'
      ianpe = 2
      backspace inpt
      read (inpt,'(a80)') wdd(1:80)
      read (wdd,*) duma1
      do i = 5,80
       if(wdd(i:i).eq.'o'.or.wdd(i:i).eq.'O') then
        ianpe = 3
        go to 100
       endif
       if(wdd(i:i).eq.'a'.or.wdd(i:i).eq.'A') then
        ianpe = 2
        go to 100
       endif 
      enddo
100   continue      
c**** read anisotropic permeability
      if(.not.allocated(anxy)) then
       allocate(anxy(n0),anxz(n0),anyz(n0))
      endif
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = 0.0d00
      default(2) = 0.0d00
      default(3) = 0.0d00
      igroup = 1

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     &     default, macroread(14), macro, igroup, ireturn,
     &     r8_1=anxy(1:n0),r8_2=anxz(1:n0),r8_3=anyz(1:n0)) 

      macroread(14) = .TRUE.

      end
