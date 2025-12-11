      subroutine  slvesu ( neq,ib,ncon,nop,npvc,npvt,igaus,irb,iirb,
     *                     nopt,idum,nsza,nszb,nszom,nelmd,nnop,
     *                     idof,north,noppar,iout,iatty,irdof,nbnd,accm)
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
CD1 To perform the symbolic factorization and calculate computer
CD1 storage necessary for the solvers.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-12-92     G. Zyvoloski   00097   Initial Implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/slvesu.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Thu Sep 12 08:25:28 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.5   Fri Feb 02 11:42:12 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   05/12/94 10:07:10   llt
CD2 corrected pvcs log info
CD2
CD2    Rev 1.3   05/11/94 16:22:12   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/28/94 16:42:24   robinson
CD2 Removed unneeded array.
CD2 
CD2    Rev 1.1   03/18/94 16:07:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:52   pvcs
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
CD3 neq             int    I   Number of equations
CD3 ib              int    I   Scratch space used to perform symbolic
CD3                               factorization
CD3 ncon            int    I   Connectivity matrix for the equation set
CD3 nop             int    O   Connectivity matrix for factorization
CD3                               matrix
CD3 npvc            int    I   Array of positions in ncon array where
CD3                               ncon(ncon(npvc(i))) = i
CD3 npvt            int    O   Array of pivot positions in nop matrix
CD3 igaus           int    I   Factorization level for each node
CD3 irb             int    I   Renumber array used when renumbering is
CD3                               chosen
CD3 iirb            int    I   Renumber array used when renumbering is
CD3                               chosen
CD3 nopt            int    I   Storage space needed for operations
CD3 idum            int    I   Storage space needed for operations
CD3 nsza            int    I   Size of storage available for solution
CD3                               matrix (used to check if storage
CD3                               supplied is adequate)
CD3 nszb            int    I   Size of matrix available for
CD3                                factorization matrix (used to check
CD3                                if storage supplied is adequate)
CD3 nszom           int    I/O Size of matrix available for GMRES
CD3                                algorithm calculations (used to
CD3                                check if storage supplied is
CD3                                adequate and correct if it is not)
CD3 nelmd           int    I   Storage space alloted for ncon
CD3 nnop            int    I   Storage space alloted for nop
CD3 idof            int    I   Number of degress of freedom
CD3 north           int    I   Number of orthogonalizations for the
CD3                                GMRES acceleration algorithm
CD3 noppar          int    I   Flag to specify whether renumbering is
CD3                                to be performed
CD3 iout            int    I   Unit number for writing storage
CD3                                information
CD3 iatty           int    I   Tape number for writing messages to
CD3                               output
CD3 irdof           int    I   Flag to specify whether reduced degree
CD3                                of freedom algorithm is to be used
CD3 nbnd            int    O   Maximum number of connections in ncon
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Identifier        Use  Description
CD3 
CD3 Unit number iout   O   FORTRAN I/O unit number iout specifies the
CD3                            device that error, warning, and
CD3                            informational messages are written to.
CD3 Unit number iatty  O   FORTRAN I/O unit number iatty specifies the
CD3                            device that error, warning, and
CD3                            informational messages are written to.
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
CD4 None
CD4
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4
CD4 nopcnv          N/A      Performs symbolic factorization for nodes
CD4                             in natural order
CD4 nopcnr          N/A      Performs symbolic factorization for nodes
CD4                             in permuted order
CD4 storag          N/A      Calculates the storag requirements for
CD4                             the solver routines
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
CD5 i           int         Do loop index
CD5 j           int         Do loop index
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 First, the code decides if renumbering is being employed (noppar=0)
CD6 and sets the values in the renumbering arrays if it is.  Then, if
CD6 renumbering is employed, nopcnr is called to perform symbolic
CD6 factorization; otherwise, nopcnv is called for this pupose.  Then,
CD6 storag is called to compute the storage requirements.  Finally,
CD6 the maximum number of connections in ncon is found by looping
CD6 through the pointer section of ncon and replacing the value
CD6 each time a larger number of connections is found.  The code then
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
CD8 The requirements listed in Section D9 are actually performed in
CD8 routines called by slvesu - slvesu is the routine that controls
CD8 these operations.
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 3.1.1. Compute Storage Requirements
CD9 3.1.2. Perform Symbolic Factorization
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
CPS BEGIN slvesu
CPS
CPS IF renumbering is being used
CPS
CPS   FOR all equations
CPS     Set the value of each renumber array to the current equation...
CPS...    number
CPS   ENDFOR
CPS
CPS ENDIF
CPS
CPS IF renumbering is not being used
CPS   nopcnv - perform symbolic factorization in natural order
CPS ELSE renumbering is being used
CPS   nopcnr - perform symbolic factorization in permuted order
CPS ENDIF
CPS
CPS storag - compute storage requirements
CPS
CPS Initialize the maximum number of connections in ncon to 0
CPS
CPS FOR each equation
CPS   Update maximum number of connections in ncon
CPS ENDFOR
CPS
CPS END slvesu
C**********************************************************************

! comxi added to allow use of input name for nop.temp file
      use comxi
      implicit none

      integer i
      integer iatty
      integer ib(*)
      integer idof
      integer idum(*)
      integer igaus(*)
      integer iirb(*)
      integer iout
      integer irb(*)
      integer irdof
      integer j
      integer k
      integer nbnd
      integer ncon(*)
      integer neq
      integer nelmd
      integer nnop
      integer nop(*)
      integer noppar
      integer nopt(*)
      integer north
      integer npvc(*)
      integer npvt(*)
      integer nsza
      integer nszb
      integer nszom

      logical used
      integer neq_temp, nsize_ncon, nsize_nop, igaus_temp
      integer io_status, io_nop
      integer, allocatable ::  npvc_temp(:)
      character*4 accm
      used = .false.
c
c**** if noppar eq 0 set default renumbering(to grid input numbers ****
c
      if ( noppar .eq. 0 )  then
        do   10   j=1,neq
           iirb(j) =  j
           irb (j) =  j
   10   continue
      endif
      
      io_nop = nufilb(28)
      inquire(file=nmfil(28), exist=used)
      open(io_nop,file=nmfil(28),status=cstats(28),form=cform(28))
      if(used) then
         if (iout .ne. 0) 
     &        write(iout,*) '>>>reading nop from file ', 
     &        trim(nmfil(28)), ' .....'
         if(iatty.ne.0) 
     &        write(iatty,*) '>>>reading nop from file ', 
     &        trim(nmfil(28)), ' .....'
         read(io_nop,iostat=io_status)
     &        neq_temp, nsize_ncon, nsize_nop, igaus_temp
         if(neq_temp.eq.neq.and.nsize_ncon.eq.ncon(neq+1)
     &        .and.igaus_temp.eq.igaus(1)) then
            nnop = nsize_nop
            allocate (npvc_temp(neq))
            do i = 1,neq
             npvc_temp(i) = npvc(i)
            enddo
            read(io_nop,iostat=io_status) (nop(i),i=1,nsize_nop)
            read(io_nop,iostat=io_status) (npvt(i),i=1,neq)            
            read(io_nop,iostat=io_status) (npvc(i),i=1,neq)            
            if (iout .ne. 0) write(iout,*) 
     &           '>>>reading nop was succesful .....'
            if(iatty.ne.0) 
     &           write(iatty,*) '>>>reading nop was succesful .....'
             do i = 1,neq
              if(npvc(i).ne.npvc_temp(i)) then
               if (iout .ne. 0) write(iout,*) 
     &           '>>>diagonal check failed .....'
               if(iatty.ne.0) write(iatty,*) 
     &           '>>>diagonal check failed .....'
               if (iout .ne. 0) 
     &          write(iout,*) '>>>generating nop .....'
               if(iatty.ne.0) 
     &          write(iatty,*) '>>>generating nop .....'
               used = .false.  
               do k = 1,neq
                npvc(k) = npvc_temp(k)
               enddo
               go to 100            
              endif
             enddo  
         else
            if (iout .ne. 0) write(iout,*) 
     &           '>>>stopped reading, nop did not match .....'
            if(iatty.ne.0) write(iatty,*) 
     &           '>>>stopped reading, nop did not match .....'
            used = .false.
         endif
      else
         if (iout .ne. 0) 
     &        write(iout,*) '>>>generating nop .....'
         if(iatty.ne.0) 
     &        write(iatty,*) '>>>generating nop .....'
      endif
100   close(io_nop)
      if(allocated(npvc_temp))deallocate (npvc_temp)
c
      if ( noppar .le. 0 .and. .not.used)  then
c     
         call  nopcnv  (ib,igaus,neq,ncon,nop,npvc,npvt,nopt,idum,nszb)
         if(nszb.lt.0) return
c
      else if(.not.used) then
c
         call  nopcnr  ( ib,igaus,neq,ncon,nop,npvc,npvt,nopt,irb,iirb,
     *        idum,nszb )
c
      endif
      if(.not.used) then
         open(io_nop,file=nmfil(28),status=cstats(28),form=cform(28))
         if (iout .ne. 0) 
     &        write(iout,*) '>>>writing nop to file ', trim(nmfil(28)),
     &        ' .....'
         if(iatty.ne.0) 
     &        write(iatty,*) '>>>writing nop to file ', trim(nmfil(28)),
     &        ' .....'
         nsize_nop = nop(neq+1)
         nnop = nsize_nop
         write(io_nop) neq,ncon(neq+1),nsize_nop,igaus(1)
         write(io_nop) (nop(i),i=1,nsize_nop)
         write(io_nop) (npvt(i),i=1,neq)
         write(io_nop) (npvc(i),i=1,neq)
         write(io_nop) 
         close(io_nop)
      endif
c
      call  storag  ( neq,ncon,nop,north,idof,nsza,nszb,nszom,iout,
     *     iatty,irdof,nelmd,nnop,accm)
c
      nbnd   =  0
c
      do   20   i=1,neq
         nbnd =  max0( ncon(i+1)-ncon(i),nbnd )
 20   continue
c
      return
      end
