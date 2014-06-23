      subroutine  storag  ( neq,ncon,nop,north,idof,nsza,nszb,nszom,
     *                      iout,iatty,irdof,nelmd,nnop,accm)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To calculate storage requirements and terminate execution if any
CD1 array has insufficient space.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-13-92     G. Zyvoloski   00097   Initial Implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/storag.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:25:40 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:08:42 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   Thu Jan 11 12:51:22 1996   gaz
CD2 fixed requirements for reduced degree of freedom
CD2 
CD2    Rev 1.1   03/18/94 16:07:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:22   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier      Type  Use  Description
CD3
CD3   SOME are I and O if the change works!!!!!!!!!!1
CD3
CD3 neq             int    I   Number of equations
CD3 ncon            int    I   Connectivity matrix for the equation set
CD3 nop             int    O   Connectivity matrix for factorization
CD3                            matrix
CD3 north           int    I   Number of orthogonalizations for the
CD3                            GMRES acceleration algorithm
CD3 idof            int    I   Number of degress of freedom
CD3 nsza            int    I/O Size of storage available for solution
CD3                            matrix (used to check if storage supplied
CD3                            is adequate)
CD3 nszb            int    I/O Size of matrix available for
CD3                            factorization matrix (used to check if
CD3                            storage supplied is adequate)
CD3 nszom           int    I/O Size of matrix available for GMRES
CD3                            algorithm calculations (used to check if
CD3                            storage supplied is adequate and correct
CD3                            if it is not)
CD3 iout            int    I   Unit number for writing storage
CD3                            information
CD3 irdof           int    I   Flag to specify whether reduced degree of
CD3                            freedom algorithm is to be used
CD3 nelmd           int    I/O Storage space alloted for ncon
CD3 nnop            int    I/O Storage space alloted for nop
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Unit number iout   O   FORTRAN I/O unit number iout specifies the
CD3                        device that error, warning, and informational
CD3                        messages are written to.
CD3
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 None
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                     Common
CD4 Identifier  Type    Block         Description
CD4
CD4
CD4 Global Subprograms
CD4
CD4 None
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier  Type        Description
CD5
CD5 error_flag  int         Flag used to determine whether to return to
CD5                         calling routine or exit
CD5 mcs         int         Parameter used in computation of storage
CD5                         space required
CD5 mct         int         Parameter used in computation of storage
CD5                         space required
CD5 north_old   int         Value of north upon input (north is lowered
CD5                         if necessary)
CD4 ist1        int         Maximum allowable size of
CD4                            ncon array
CD4 ist2        int         Maximum allowable size of
CD4                            nop array
CD4 ist3        int         Maximum allowable size of
CD4                            the array a
CD4 ist4        int         Maximum allowable size of
CD4                            the array b
CD4 ist5        int         Maximum allowable size of
CD4                            acceleration arrays
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 The code computes the storage requirements for the arrays, then
CD6 checks to see if this is a reduced degree of freedom solution and
CD6 corrects the storage requirements accordingly.  Then, the code
CD6 writes the values just computed to the specified file (or screen,
CD6 depending on the value of iout), and performs a series of error 
CD6 checks, terminating execution if the space for ncon, nop, a, or b 
CD6 exceeds the maximum.  Then, the code reduces the number of
CD6 orthogonalizations if the space for the GMRES arrays is too small.
CD6 In this case, the code issues a warning message.  Then, execution
CD6 returns to the calling routine.
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 N/A
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This routine stops execution of the program within the routine if
CD8 memory has not been allocated properly, rather than passing a flag
CD8 to the calling routine.  Although the latter is probably better
CD8 design, the practice of stopping execution within this routine has
CD8 been used in past versions of this code.  Changing it now would
CD8 require that those applications would have to carry out their own
CD8 error check in order to terminate execution.  It is better to
CD8 preserve the current behavior of this reuse component so that these
CD8 other applications need not be changed.  Furthermore, the user would
CD8 never wish to continue execution, since insufficient allocation
CD8 of array space would always be a fatal error.
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 3.1.1. Compute Storage Requirements
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA GZSOLVE SRS.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN storag
CPS
CPS Compute storage required for all arrays
CPS
CPS IF a reduced degree of freedom problem is specified
CPS   Recompute storage for b matrix and GMRES arrays
CPS ENDIF
CPS
CPS Write storage requirements to specified output device
CPS
CPS Set error flag to indicate successful execution
CPS
CPS IF space is insufficient for ncon array
CPS   Write error message to specified output device
CPS   Set error flag to indicate termination of execution
CPS   ERROREXIT
CPS ENDIF
CPS
CPS IF space is insufficient for nop array
CPS   Write error message to specified output device
CPS   Set error flag to indicate termination of execution
CPS   ERROREXIT
CPS ENDIF
CPS
CPS IF space is insufficient for the array a
CPS   Write error message to specified output device
CPS   Set error flag to indicate termination of execution
CPS   ERROREXIT
CPS ENDIF
CPS
CPS IF space is insufficient for the array b
CPS   Write error message to specified output device
CPS   Set error flag to indicate termination of execution
CPS   ERROREXIT
CPS ENDIF
CPS
CPS IF space is insufficient for acceleration arrays
CPS   Compute maximum number of orthogonalizations allowed
CPS   Set acceleration storage parameter
CPS   Write warning message to specified output device
CPS ELSE
CPS   Set acceleration storage parameter
CPS ENDIF
CPS
CPS ERRORSEGMENT
CPS
CPS   IF a fatal error was detected
CPS     exit - terminate execution
CPS   ELSE
CPS     Set array length parameters to their exact values
CPS   ENDIF
CPS   
CPS ENDSEGMENT
CPS
CPS END storag
CPS 
C**********************************************************************
c****--------------------------------------------------------------****c
c**** calculates storage requirements this routine called          ****c
c****                                                  from slvesu ****c

      implicit none

      integer error_flag
      integer iatty
      integer idof
      integer iout
      integer irdof
      integer ist1
      integer ist2
      integer ist3
      integer ist4
      integer ist5
      integer mcs
      integer mct
      integer ncon(*)
      integer nelmd
      integer neq
      integer nnop
      integer nop(*)
      integer north
      integer north_old
      integer nsza
      integer nszb
      integer nszom
      character*4 accm


c
      mct    =  idof
      mcs    =  idof**2
      ist1   =  ncon(neq+1)
      ist2   =  nop(neq+1)
      ist3   =  ( ist1-neq-1 )*mcs
      ist4   =  ( ist2-neq-1 )*mcs
      if(accm.eq.'bcgs') then
       ist5   =  ( 3 + 2 * ( abs(north)+1 ) )*mct*neq
      elseif(accm.eq.'gmre') then
       ist5   =  ( north+1 )*mct*neq
      elseif(accm.eq.'bcgn') then
      endif
c
c**** modify storage for reduced dof systems ****
c
      if ( irdof .ne. 0 )  then
       ist4   = ( ist2-neq-1 )*irdof**2
       ist5   = ( 3 + 2 * ( abs(north)+1 ) )*irdof*neq
       if(accm.eq.'bcgs') then
        ist5   =  ( 3 + 2 * ( abs(north)+1 ) )*irdof*neq
       elseif(accm.eq.'gmre') then
        ist5   =  ( north+1 )*irdof*neq
       elseif(accm.eq.'bcgn') then
       endif
      endif
c
      if (iout .ne. 0)     write(iout  ,6000)  ist1 , nelmd
      if ( iatty .gt. 0 )  write(iatty ,6000)  ist1 , nelmd
 6000 format(/,1x,'storage needed for ncon     ',i10,' available ',i10)
      if (iout .ne. 0)     write(iout  ,6001)  ist2 , nnop
      if ( iatty .gt. 0 )  write(iatty ,6001)  ist2 , nnop
 6001 format(  1x,'storage needed for nop      ',i10,' available ',i10)
      if (iout .ne. 0)     write(iout  ,6002)  ist3 , nsza
      if ( iatty .gt. 0 )  write(iatty ,6002)  ist3 , nsza
 6002 format(  1x,'storage needed for a matrix ',i10,' available ',i10)
      if (iout .ne. 0)     write(iout  ,6003)  ist4 , nszb
      if ( iatty .gt. 0 )  write(iatty ,6003)  ist4 , nszb
 6003 format(  1x,'storage needed for b matrix ',i10,' available ',i10)
      if (iout .ne. 0)     write(iout  ,6004)  ist5 , nszom
      if ( iatty .gt. 0 )  write(iatty ,6004)  ist5 , nszom
 6004 format(  1x,'storage for acceleration',i10,' available ',i10)
c
      error_flag = 0
      if(ist1.gt.nelmd) then
         if (iout .ne. 0) then
            write(iout, 100) 
            write(iout, 200) 'ncon'
         end if
         if ( iatty .gt. 0 ) then
           write(iatty, 100) 
           write(iatty, 200) 'ncon'
         end if
         error_flag = -1
         goto 9000
      end if
      if(ist2.gt.nnop) then
         if (iout .ne. 0) then
            write(iout, 100) 
            write(iout, 200) 'nop'
         end if
         if ( iatty .gt. 0 ) then
            write(iatty, 100) 
            write(iatty, 200) 'nop'
         end if
         error_flag = -1
         goto 9000
      end if
      if(ist3.gt.nsza) then
         if (iout .ne. 0) then
            write(iout, 100)
            write(iout, 200) 'a'
         end if
         if ( iatty .gt. 0 ) then
            write(iatty, 100) 
            write(iatty, 200) 'a'
         end if
         error_flag = -1
         goto 9000
      end if
      if(ist4.gt.nszb) then
         if (iout .ne. 0) then
            write(iout, 100) 
            write(iout, 200) 'b'
         end if
         if ( iatty .gt. 0 ) then
            write(iatty, 100)
            write(iatty, 200) 'b'
         end if
         error_flag = -1
         goto 9000
      end if
      if(ist5.gt.nszom) then
         nszom  =  ist5                              
         if (iout .ne. 0) then
            write(iout,*) 'Warning (Issued from GZSOLVE-storag)'
            write(iout,*) 'Storage for acceleration method changed '
         end if
         if ( iatty .gt. 0 ) then
            write(iatty,*) 'Warning (Issued from GZSOLVE-storag)'
            write(iatty,*) 'Storage for acceleration method changed '
         end if
      else
         nszom = ist5
      end if
 100  format ('Fatal error in GZSOLVE-storag')
 200  format ('Not enough space allocated in the ', 4a, ' array')
 9000 continue
      if( error_flag .ne. 0 ) then
         stop
      else
         nelmd = ist1
         nnop = ist2
         nsza = ist3
         nszb = ist4
      end if
      return
      end
