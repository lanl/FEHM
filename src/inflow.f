      subroutine inflow
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
CD1 Read flow data.
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
CD2 $Log:   /pvcs.config/fehm90/src/inflow.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Wed Jun 12 15:21:04 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.6   Mon Jun 03 11:17:56 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.5   Fri May 31 10:44:04 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.4   Fri Feb 16 08:59:44 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.3   Tue Jan 30 11:50:36 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   Fri Jan 12 17:50:20 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.1   03/18/94 16:03:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:00   pvcs
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
CD4   ipesk           POINTER  fdd    Pointer to array variable esk
CD4   ippflow         POINTER  fdd    Pointer to array variable pflow
CD4   ischk           INT      faai   Unit number for input data check file
CD4   ka              INT      fbb    Contains boundary type information for
CD4                                     each node
CD4   macroread(12)   LOGICAL  macro  Flag denoting if macro flow has been read 
CD4   pflow           REAL*8   fdd    Flowing pressure at each source node
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   wellim          REAL*8   fdd    Well impedance at each source nodeCD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   initdata                 Read data values and set parameter values at
CD4                              given nodes
CD4   mmgetblk                 Allocate memory to an array
CD4   mmrelblk                 Deallocate array memory
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
CD5   aiped           REAL*8   Temporary space for impedance parameter
CD5   default         REAL*8   Array of default values for input arrays
CD5   esktmp          REAL*8   Temporary storage for enthalpy/temperature
CD5   i               INT      Loop index
CD5   icode           INT      Return value from mmgetblk, mmrelblk routines
CD5   igroup          INT      Current group number in this macro
CD5   ipaiped         POINTER  Pointer to variable array aiped
CD5   ipesktmp        POINTER  Pointer to variable array esktmp
CD5   ipsktmp         POINTER  Pointer to variable array sktmp
CD5   ireturn         INT      Returned error flag from input subroutine
CD5   itype           INT      Array of variable types being input
CD5   macro           CHAR     Current macro being read
CD5   narrays         INT      Number of arrays being read in
CD5   pointer         INT      Integer array of pointer values for
CD5                              variables being read
CD5   sktmp           REAL*8   Temporary storage for source strength
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
CD9 2.3.7 Sources and sinks
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN inflow
CPS 
CPS   define temporary storage arrays and allocate space
CPS   
CPS   set arrays used, pointers, data type and default values
CPS   call initdata to read data values and set parameter values at 
CPS   given nodes
CPS   
CPS   FOR each node
CPS       IF the temporary variables are not the default values
CPS          set the enthalpy or temperature
CPS          IF the temporary impedance variable is zero
CPS             set the flow and boundary type to 1
CPS          ELSE
CPS             IF the temporary impedance variable is negative
CPS                set the boundary type to -2
CPS             ELSE
CPS                set the boundary type to -1
CPS             END IF
CPS             set the impedance and flowing pressure
CPS          END IF
CPS       END IF
CPS   END FOR
CPS   
CPS   deallocate temporary storage
CPS
CPS END inflow
CPS
C***********************************************************************

      use comai
      use combi
      use comdi
      use comdti
      use comki
      use comii
      implicit none

      integer i,icode
      real*8, allocatable :: aiped(:)
      real*8, allocatable::  sktmp(:)
      real*8, allocatable ::  esktmp(:)
      real*8  sat_mult
      parameter (sat_mult = 1.d-5)
      macro = 'flow'
c gaz 112920 added sk0 to keep original water source    
      allocate(aiped(n0),esktmp(n0),sktmp(n0))
	      
c**** read flow data ****
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = 0.
      default(2) = 0.
      default(3) = 0.
      igroup = 1
      
      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(12), macro, igroup, ireturn,
     *     r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0)) 
      
      do i = 1, n0
         if (sktmp(i) .ne. default(1) .or. esktmp(i) .ne. default(2)
     *        .or. aiped(i) .ne. default(3)) then
            esk(i) = esktmp(i)
            if(esk(i).gt.0.and.igrav.ne.0.and.ico2.ge.0) then
             esk(i) = esktmp(i)-grav*cord(i,igrav)
            endif
            if (abs(aiped(i)) .lt. zero_t) then
               sk(i) = sktmp(i) 
c gaz 112920 keep copy of water source               
               sk0(i) = sk(i)
               ka(i) = 1
            else
               if (aiped(i) .lt. 0.) then
                  ka(i) = -2
               else
                  ka(i) = -1
               end if
               wellim(i) = abs(aiped(i)) * 1.0e+06
               pflow(i) = sktmp(i)
c gaz 110719 if pflow< pref then turn off node (sort of)
c rich only               
               if(pflow(i).lt.pref.and.jswitch.ne.0) then
                wellim(i) = wellim(i)*sat_mult
               endif
            end if
         end if
      end do
      
      deallocate(aiped,esktmp,sktmp)
           
      end
