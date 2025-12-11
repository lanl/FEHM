      subroutine namefile2(icall,lund,iformat,tail,idz)
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
CD1 Generate name for AVS format output file.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JUN-94    Z. Dash        22      Add partial prolog and fix file name usage.
CD2              C. Gable               Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/namefile2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:06   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:46 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Jan 18 13:19:42 1996   zvd
CD2 Modified prolog
CD2 
CD2    Rev 1.4   04/25/95 09:55:26   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.3   06/20/94 10:57:16   zvd
CD2 Added ierr unit number for error output.
CD2 Modified to use longer file prefix for naming AVS output files.
CD2 
CD2    Rev 1.2   04/06/94 10:51:12   tam
CD2 added function c_open() to open binary files for write
CD2
CD2    Rev 1.1   03/18/94 15:43:10   gaz
CD2 Added solve_new and cleaned up memory management.
CD2
CD2    Rev 1.0   01/20/94 10:25:44   pvcs
CD2 original version in process of being certified
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   icall           INT      I/O  Number of times this routine has been called
CD3   lund            INT      I    Logical unit number of file to be opened
CD3   iformat         INT      I    Flag to indicate if file is to be formatted
CD3   root            CHAR     I    Filename root to which suffix will be added
CD3   tail            CHAR     O    Filename suffix
CD3   nchar_root      INT      I    Number of characters in input
CD3                                   filename prefix
CD3   ierr            INT      I    Unit number for error output
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
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
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4   
CD4   c_open          INT      Open AVS binary format output file
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   Identifier      Type     Description
CD5
CD5
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
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
CPS
C***********************************************************************
c     icall= the number of times this routine has been called.
c            This value is used in putting a suffix on the output
c            file. The file names will have sequential suffixes. 
c            Suffix will be between 10001 and 99999 if icall 
c            .gt. 0 and .lt. 89999.
c            You are on your own if you want more than 90000 files!
c
c     lund = logical unit number of file to be opened
c
c     root = raster file's name root to which suffix will be added.
c            must be declared as character*94 in calling program.
c
c     iformat = 1 unformatted file is opened
c             = 2 formatted file is opened
c
c     nchar_root = number of characters in input file name prefix
C Modified 12/03/93 by ZVD
C Remodified 06/10/94 by ZVD
C Added nchar_root passed from calling routine 
C (max = 94 + 1 same as other routines using the file name prefix)
C Allowed file names to be 125 characters (rest of code uses 100)
C***********************************************************************

      use avsio
      use comai, only : altc, contim, days, ierr
      implicit none

      integer i, icall, lund, iformat, idz, open_file
      integer nchar, nchar_tail, c_open, len_char, k1, k2, kl
      character*165 fname
      character*94 root
      character*40 daychar
      character*(*) tail
      character*5 ch5
      character*4 ch4

      write(ch4, '(i4.4)') idz

      if (contim .ge. 0) then
         write(ch5, '(i5.5)') icall
      else
         daychar = ''
c gaz 102623 g16.9 tp g16.8
         if(days.eq.0.0) then
          days = 0.0
          write(daychar,'(1p,g16.9)') days
         else
          write(daychar,'(1p,g16.8)') days
         endif
         k1 = 0
         do 
            k1 = k1 + 1
            if (daychar(k1:k1) .ne. ' ') exit
         end do
         k2 = k1
         do 
            k2 = k2 + 1
            if (daychar(k2:k2) .eq. ' ') exit
         end do
         len_char = k2 - k1
         k2 = k2 - 1
      end if

c Blank out the file name
      fname = ''
      nchar_tail = len_trim(tail)

c     name output file
      fname(1:iaroot) = avs_root(1:iaroot)
      nchar = iaroot
      if (altc(1:3) .eq. 'sur' .and. idz .ne. 0) then
c Add zone number to name
         fname(nchar+1:nchar+4) = ch4(1:4)
         fname(nchar+5:nchar+5) = '_'
         nchar = len_trim(fname)
      end if
      if (contim .ge. 0) then
         fname(nchar+1:nchar+5) = ch5(1:5)
      else
         fname(nchar+1:nchar+len_char) = daychar(k1:k2)
         nchar = len_trim(fname)
         fname(nchar+1:nchar+5) = '_days'
      end if
      nchar = len_trim(fname)
      fname(nchar+1:nchar+nchar_tail)=tail(1:nchar_tail)
      nchar = len_trim(fname)
      if (altc(1:3) .eq. 'tec') then
         fname(nchar+1:nchar+4) = '.dat'
      else if (altc(1:3) .eq. 'sur') then
         fname(nchar+1:nchar+4) = '.csv'
      else if (altc(1:4) .eq. 'avsx') then
         fname(nchar+1:nchar+5) = '.avsx'
      else if (altc(1:3) .eq. 'avs') then
         fname(nchar+1:nchar+4) = '.avs'
      end if
      nchar = len_trim(fname)
      
c     open to output file
      if(iformat .eq. 1)then
         continue
      elseif(iformat .eq. 2)then
!         open(unit=lund,file=fname(1:nchar),form='formatted')
         lund = open_file(fname(1:nchar),'unknown')
      endif         

      return
      end
