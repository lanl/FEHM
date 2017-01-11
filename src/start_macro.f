      subroutine start_macro(inputnum, locunitnum, macro)
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
CD1 Allows input to be read from any file.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/start_macro.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Wed Jun 26 12:40:50 1996   hend
CD2 Updated Prolog
CD2 
CD2    Rev 1.2   Wed Jun 12 15:21:22 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.1   Mon Jun 03 11:18:38 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.0   Fri May 31 10:41:40 1996   hend
CD2 Added optional input from specified file
CD2
C**********************************************************************
CD3
CD3 REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3
C**********************************************************************
CD4
CD4 SPECIAL COMMENTS AND/OR REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      use comai
      implicit none

      integer inputnum, locunitnum
      character*4 lockeyword, macro, macroname1, macroname2
      character*200 locfilename
      logical ex
c check is user wants to read macro data from alternative data file

      read(inputnum,'(a4)') lockeyword
      if (lockeyword(1:4) .eq. 'file') then

         read(inputnum, '(a200)', err = 10) locfilename
         
         if (iout .ne. 0) write(iout, 100) macro, locfilename      
         if (iptty .gt. 0) write(iptty, 100) macro, locfilename
 100     format (1x, a4, ' read from optional input file: ', a200)

* Make sure file exists
         inquire (file = locfilename, exist = ex)
         if (.not. ex) then
            write (ierr, 200) locfilename
            if (iptty .gt. 0) write (iptty, 200) locfilename
            goto 30
 200        format ('ERROR nonexistant file', a200)
         end if

         open(iocntl, file = locfilename, err = 20)
         locunitnum = iocntl

* Skip over comments at start of file
 1       read(locunitnum,'(a4)') lockeyword
         if (lockeyword(1:1).eq.'#') goto 1

* Verify reading correct macro
         if (macro(1:3) .eq. 'air') then
            macroname1 = 'air '
            macroname2 = 'airw'
         else if (macro .eq. 'ngas' .or. macro .eq. 'co2i') then
            macroname1 = 'ngas'
            macroname2 = 'co2i'
         else
            macroname1 = macro
            macroname2 = macro
         end if
            
         if (lockeyword .ne. macroname1 .and. 
     .        lockeyword .ne. macroname2) then
            write(ierr, 300) macro, lockeyword
            if (iptty .gt. 0) write (iptty, 300) macro, lockeyword
            goto 30
 300        format ('ERROR  --> Macro name in file for macro ', a4,
     &           ' is ', a4)
         endif

      else

         backspace inputnum
         locunitnum = inputnum

      endif

      return

* ERROR EXIT
 10   write (ierr, *) 'ERROR reading optional input file name'
      if (iptty .gt. 0) write (iptty, *) 'ERROR reading optional input',
     &     ' file name'
      goto 30
 20   write (ierr, *) 'ERROR opening ', locfilename
      if (iptty .gt. 0) write (iptty, *) 'ERROR opening ', locfilename
      goto 30
 30   write (ierr, *) 'STOPPED trying to use optional input file'
      if (iptty .gt. 0) write (iptty, *) 'STOPPED trying to use ',
     &     'optional input file'
      stop

      end
