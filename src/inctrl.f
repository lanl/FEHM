      subroutine inctrl
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
CD1 Read control variables.
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
CD2 $Log:   /pvcs.config/fehm90/src/inctrl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:02 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:17:52 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:59:56 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Tue Jan 30 11:34:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   01/28/95 13:55:10   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 16:03:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:52   pvcs
CD2 original version in process of being certified
CD2 
c 12/22/94 gaz read fill-in in directly to nar
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
CD3   Name                     Use  Description
CD3
CD3   inpt                     I    Input data file
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
CD4   aiaa            REAL*8   faar   Time step multiplication factor
CD4   aw              REAL*8   faar   Time step weighting parameter for heat and
CD4                                     mass solution
CD4   awt             REAL*8   faar   Value of implicitness factor
CD4   ay              REAL*8   faar   Time step weighting parameter for tracer
CD4   daymax          REAL*8   faar   Maximum time step allowed
CD4   daymin          REAL*8   faar   Minimum time step allowed
CD4   epe             REAL*8   faar   Tolerance for newton-raphson iteration
CD4   grav            REAL*8   faar   Value of gravity
CD4   iamm            INT      faai   Maximum iterations allowed for time step
CD4                                     increase (heat and mass solution)
CD4   iamx            INT      faai   Iteration count after which the time step
CD4                                     will be halved
CD4   icnl            INT      faai   Problem dimension
CD4   igrav           INT      faai   Direction of gravity in problem
CD4   inpt            INT      faai   Unit number for input file
CD4   ipnar           POINTER  fbb    Pointer to variable array nar    
CD4   ischk           INT      faai   Unit number for input data check file
CD4   lda             INT      faai   Parameter which specifies if the geometric
CD4                                     coefficients are saved
CD4   macroread(11)   LOGICAL  macro  Flag denoting if macro ctrl has been read 
CD4   maxit           INT      faai   Maximum number of iterations allowed
CD4   north           INT      faai   Maximum number of orthogonalizatins
CD4   upwgt           REAL*8   faar   Upwind weighting parameter
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
CD5   aaw             REAL*8   Implicitness factor flag
CD5   agrav           REAL*8   Gravity direction flag
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
CPS BEGIN inctrl
CPS 
CPS   read control variables in group 1
CPS   
CPS   set array used, pointer, data type and default value
CPS   call initdata to read data values and set parameter values at 
CPS   given nodes
CPS   
CPS   set macroread to true
CPS   
CPS   read control variables in groups 3, 4, 5 and set appropriate 
CPS   values based on data read
CPS   
CPS END inctrl
CPS
C***********************************************************************

      use comai
      use combi
      use comdti
      use comki
      implicit none

      real*8 aaw, agrav
      character*6 chdum

      macro = 'ctrl'

c**** read control variables ****

c zvd Add option to suppress output of nodal equation residuals
      resid_out = .true.
      backspace (inpt)
      read (inpt,'(a80)') wdd1
      read (wdd1, *, end = 5) chdum, chdum
      if (chdum(1:2) .eq. 'no' .or. chdum(1:2) .eq. 'NO')
     &     resid_out = .false.
c gaz read (inpt,*) maxit, epe, north
 5    read (inpt,'(a80)') wdd1         
      read (wdd1,*,end=10) maxit, epe, north, maxsolve, accm
      go to 20
 10   continue
      if(north.le.0) then
         accm = 'bcgs'
         maxsolve = max(3*abs(north),60)
         north = 1
      else if(north.gt.0) then
         accm = 'gmre'
         maxsolve = 3*north
      else
         write(ierr, 15)
         if (iout .ne. 0) write(iout, 15)
         if (iatty .ne. 0) write(iatty, 15)
         stop
      endif
 15    format ('ERROR : north = 0, stopping')
 20    continue
      north = abs(north)

      narrays = 1
      itype(1) = 4
      default(1) = 1
      igroup = 2

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(11), macro, igroup, ireturn,
     2     i4_1=nar(1:n0)) 

      macroread(11) = .TRUE.

      read (inpt,*) aaw, agrav, upwgt
      if (aaw .le. 1.0) then
         aw = 1.0
      else
         aw = 1.5
      end if

      awt = aw
      aw = 1.0
      ay = 1.0 - aw

      if (agrav .gt. 3.0) agrav = 3.0
      if (abs(agrav) .gt. zero_t) then
         igrav = agrav
         grav = 9.81
      else
         igrav = 3
         grav = 0.0
      end if
      if (upwgt .gt. 1.0) then
         upwgt = 1.0
      else if (upwgt .lt. 0.5) then
         upwgt = 0.5
      end if

      aiar = 0.
      read (inpt, '(a80)') wdd1
      read (wdd1,*,end=30) iamm, aiaa, daymin, daymax, aiar
      aiar = max (aiar, 1.5d0)
      goto 40
 30   continue
      aiar = 1.5
 40   continue
      iamx = min0(iamm, maxit)
      if(aiaa .lt. 1.0) aiaa = 1.0
      read (inpt, *) icnl,lda
      if(ivf.eq.-1) then
         if(icnl.ne.0) then
            if (iout .ne. 0) write(iout, 100) 
            if (iptty.ne.0) write(iptty, 100) 
            icnl = 0
         endif
         if(lda.ne.0) then
            if (iout .ne. 0) write(iout, 110) 
            if (iptty.ne.0) write(iptty, 110) 
            lda = 0
         endif
      endif
 100  format('>>> Warning: problem dimension = 3 (icnl=0) for FDM <<<')
 110  format ('>>> Warning: No stor file needed for FDM <<<')
      end
