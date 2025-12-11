       subroutine inpres
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
CD1 Read nonuniform pressure and temperature or saturation data.
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
CD2 $Log:   /pvcs.config/fehm90/src/inpres.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:08 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:18:02 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:58:00 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Tue Jan 30 13:13:48 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   Fri Jan 12 17:50:52 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.1   03/18/94 16:03:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:14   pvcs
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
CD4   ieos            INT      fddi1  Phase state of fluid at each node
CD4   ipieos          POINTER  fddi1  Pointer to array variable ieos
CD4   ippho           POINTER  fdd    Pointer to array variable pho
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(15)   LOGICAL  macro  Flag denoting if macro pres has been read 
CD4   pein            REAL*8   faar   Initial pressure of problem
CD4   pho             REAL*8   fdd    Last time step pressure at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   to              REAL*8   fdd    Last time step temperature at each node
CD4   
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   initdata                 Read data values and set parameter values at
CD4                              given nodes
CD4   mmgetblk                 Allocate memory to an array
CD4   mmrelblk                 Deallocate array memory
CD4   psat            REAL*8   Calculate the saturation temperature or pressure
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
CD5   default         REAL*8   Array of default values for input arrays
CD5   dtsatp          REAL*8   Derivative of saturation pressure with
CD5                              temperature
CD5   i               INT      Loop index
CD5   icode           INT      Return value from mmgetblk, mmrelblk routines
CD5   igroup          INT      Current group number in this macro
CD5   iptmp           POINTER  Pointer to variable array tmp
CD5   ireturn         INT      Returned error flag from input subroutine
CD5   itype           INT      Array of variable types being input
CD5   macro           CHAR     Current macro being read
CD5   narrays         INT      Number of arrays being read in
CD5   pointer         INT      Integer array of pointer values for
CD5                              variables being read
CD5   tmp           REAL*8     Temporary storage for temperatures/saturations
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
CPS BEGIN inpres
CPS 
CPS   define temporary storage arrays and allocate space
CPS   
CPS   set arrays used, pointers, data type and default values
CPS   
CPS   FOR each node
CPS       set temporary storage to default value
CPS   END FOR
CPS   
CPS   set macroread to true (don't change any values previously set)
CPS   
CPS   call initdata to read data values and set parameter values at given . . .
CPS   . . . nodes
CPS   
CPS   FOR each node
CPS       IF the temporary variable is not the default value
CPS          IF the thermodynamic region is 1
CPS             set the temperature to temporary value
CPS             set the saturation to 1
CPS          ELSE IF the thermodynamic region is 2
CPS             set the saturation to temporary value
CPS             call psat to calculate the temperature
CPS          ELSE IF the thermodynamic region is 3
CPS             set the temperature to temporary value
CPS             set the saturation to 0
CPS          END IF
CPS       END IF
CPS   END FOR
CPS   
CPS   deallocate temporary storage
CPS
CPS END inpres
CPS
C***********************************************************************

      use comdi
      use comdti
      use comai
      use comki
      use davidi
      use comco2, only: icarb
      implicit none

      integer i,icode
      real*8 dtsatp,psat
      real*8, allocatable ::  tmp(:)

      macro = 'pres'

      allocate(tmp(n0))

c**** read flow data ****
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 4
      default(1) = pein
      default(2) = -99.
      default(3) = 1
      igroup = 1

C**** Only want to set specified nodes ?
      do i = 1, n0
         tmp(i) = default(2)
      end do

      macroread(15) = .TRUE.

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(15), macro, igroup, ireturn,
     2     r8_1=pho(1:n0),r8_2=tmp(1:n0),i4_1=ieos(1:n0)) 

      do i = 1, n0
         if (tmp(i) .ne. default(2)) then
c gaz 110715 added sc phase
            if (abs (ieos(i)) .eq. 1.or.abs (ieos(i)).eq.4) then
               to(i) =  tmp(i)
c single phase set saturation for wtsi
               if(pho(i).le.0.) then
c                 s(i)=-pho(i)
c                 pho(i)=head0
               else
                  if (irdof .ne. 13 .or. ifree .ne. 0) s (i) =  1.0
!                 s (i) =  1.0
               endif
            else if (abs (ieos(i)) .eq. 2) then
c two phase
               if (irdof .ne. 13 .or. ifree .ne. 0) s (i) =  tmp(i)
               if(icarb.eq.0) then
c gaz 120512 don't change T for co2-h20 problems               
                if (ico2.ge.0) to(i) =  psat(pho(i), dtsatp, 1)
               endif
            else if (abs (ieos(i)) .eq. 3) then
c gas only
               to(i) =  tmp(i)
               s (i) =  0.0
           else if(ieos(i).eq.0) then
c gaz 091023
c heat conduction only 
               to(i) =  tmp(i)
               s (i) =  1.0          
           end if
         endif
      end do

      deallocate(tmp)

      end      
      
      subroutine inpres_solubility
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
CD1 gaz 030924 initial coding      
CD1 Read nonuniform pressure and temperature or saturation data.
CD1 Inpres_w_solubility read nonuniform pressure,saturation,fraction data.      
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
CD2 $Log:   /pvcs.config/fehm90/src/inpres.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:08 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:18:02 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:58:00 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Tue Jan 30 13:13:48 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   Fri Jan 12 17:50:52 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.1   03/18/94 16:03:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:14   pvcs
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
CD4   ieos            INT      fddi1  Phase state of fluid at each node
CD4   ipieos          POINTER  fddi1  Pointer to array variable ieos
CD4   ippho           POINTER  fdd    Pointer to array variable pho
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(15)   LOGICAL  macro  Flag denoting if macro pres has been read 
CD4   pein            REAL*8   faar   Initial pressure of problem
CD4   pho             REAL*8   fdd    Last time step pressure at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   to              REAL*8   fdd    Last time step temperature at each node
CD4   
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   initdata                 Read data values and set parameter values at
CD4                              given nodes
CD4   mmgetblk                 Allocate memory to an array
CD4   mmrelblk                 Deallocate array memory
CD4   psat            REAL*8   Calculate the saturation temperature or pressure
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
CD5   default         REAL*8   Array of default values for input arrays
CD5   dtsatp          REAL*8   Derivative of saturation pressure with
CD5                              temperature
CD5   i               INT      Loop index
CD5   icode           INT      Return value from mmgetblk, mmrelblk routines
CD5   igroup          INT      Current group number in this macro
CD5   iptmp           POINTER  Pointer to variable array tmp
CD5   ireturn         INT      Returned error flag from input subroutine
CD5   itype           INT      Array of variable types being input
CD5   macro           CHAR     Current macro being read
CD5   narrays         INT      Number of arrays being read in
CD5   pointer         INT      Integer array of pointer values for
CD5                              variables being read
CD5   tmp           REAL*8     Temporary storage for temperatures/saturations
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
CPS BEGIN inpres
CPS 
CPS   define temporary storage arrays and allocate space
CPS   
CPS   set arrays used, pointers, data type and default values
CPS   
CPS   FOR each node
CPS       set temporary storage to default value
CPS   END FOR
CPS   
CPS   set macroread to true (don't change any values previously set)
CPS   
CPS   call initdata to read data values and set parameter values at given . . .
CPS   . . . nodes
CPS   
CPS   FOR each node
CPS       IF the temporary variable is not the default value
CPS          IF the thermodynamic region is 1
CPS             set the temperature to temporary value
CPS             set the saturation to 1
CPS          ELSE IF the thermodynamic region is 2
CPS             set the saturation to temporary value
CPS             call psat to calculate the temperature
CPS          ELSE IF the thermodynamic region is 3
CPS             set the temperature to temporary value
CPS             set the saturation to 0
CPS          END IF
CPS       END IF
CPS   END FOR
CPS   
CPS   deallocate temporary storage
CPS
CPS END inpres
CPS
C***********************************************************************

      use comdi
      use comdti
      use comai
      use comki
      use comfi, only: frac_gas_iso, cnlf, cnlof
      use davidi
      use comco2, only: icarb
      implicit none

      integer i,icode
      real*8 dtsatp,psat
      real*8, allocatable ::  tmp(:)

      macro = 'pres'
      allocate(tmp(n0))
      if(.not.allocated(frac_gas_iso)) then
       allocate(frac_gas_iso(n))
       frac_gas_iso = 0.0d0
      endif

c**** read pres data ****
      narrays = 4
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      itype(4) = 4
      default(1) = pein
      default(2) = -99.
      default(3) = -99.
      default(4) = 1
      igroup = 1

C**** Only want to set specified nodes ?
      do i = 1, n0
         tmp(i) = default(2)
      end do

      macroread(15) = .TRUE.

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(15), macro, igroup, ireturn,
     2     r8_1=pho(1:n0),r8_2=tmp(1:n0),r8_3=frac_gas_iso(1:n0),
     3     i4_1=ieos(1:n0)) 

      do i = 1, n0
       if(ieos(i).eq.0) then
c gaz 091023
c heat conduction only 
        to(i) =  tmp(i)
        s(i) =  1.0          
       else
        s(i) = tmp(i)
       end if
      end do

      deallocate(tmp)
c gaz 042124  added old mass fraction  cnlof()   
      do i = 1, n0
        cnlf(i) = frac_gas_iso(i)  
        cnlof(i) =  cnlf(i)
      enddo

      end
