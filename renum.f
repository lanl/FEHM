      subroutine renum(iflg)
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
CD1 Read node renumbering data.
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
CD2 $Log:   /pvcs.config/fehm90/src/renum.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Mon Mar 31 08:41:56 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:18 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:18:28 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:54:56 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Thu Feb 01 15:48:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   03/28/94 16:40:38   robinson
CD2 Removed unneeded array.
CD2 
CD2    Rev 1.1   03/18/94 16:03:52   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:12   pvcs
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
CD4   ig              INT      fhh    Variable order lu decomposition
CD4   inpt            INT      faai   Unit number for input file
CD4   ipig            POINTER  fhh    Pointer to variable array ig
CD4   ischk           INT      faai   Unit number for input data check file
CD4   macroread(16)   LOGICAL  macro  Flag denoting if macro renm has been read 
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
CPS BEGIN renum
CPS 
CPS   set array used, pointer, data type and default value
CPS   call initdata to read data values and set node number for given . . .
CPS   . . . nodes
CPS   
CPS   set macroread to true
CPS
CPS END renum
CPS
C***********************************************************************

      use comei
      use comdti
      use comai
      use combi
      use comki
      implicit none

    
      integer i,iflg

      macro = 'renm'

c  Note - still need to compute irb based on iirb, set noppar to 1,
c         put noppar in common block (probably comai.h), and change
c         the call to slvesu (startup.f) to pass in noppar instead
c         of 0, which is currently hardwired in the calling statement.
c GAZ- parameter ireord is the global variable that 
c      noppar represents in slvesu
c      pointer should be ipirb
c
      if(ireord.eq.0) return
      if(iflg.eq.0) then
c read data
       read(inpt,*) ireord
       if(ireord.lt.0) then
        narrays = 1
        itype(1) = 4
        default(1) = 0
        igroup = 1
 
        call initdata2 (inpt, ischk, neq, narrays,  itype, 
     *       default, macroread(16), macro, igroup, ireturn,
     2       i4_1=iirb(1:neq)) 

        macroread(16) = .TRUE.
       endif
      else if(iflg.eq.1) then
c compute ordering
       if(ireord.eq.1) then
        call genrcm(neq,nelm,irb,nop(1),nop(neq+1))
       else if(ireord.ge.10) then
         do i=1,neq
          irb(i)=i
          iirb(i)=i
         enddo
       endif
      endif
      end
