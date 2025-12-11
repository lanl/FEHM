      subroutine cntlin
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
CD1 Read I/O file names from control file, set unit numbers and open files.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 16-NOV-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/cntlin.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Mon Jan 29 13:33:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   08/18/95 10:24:44   llt
CD2 on the cray, since the end do has a statement label of 5, the
CD2 do also needs a label of 5
CD2 
CD2    Rev 1.3   06/20/94 11:09:14   zvd
CD2  
CD2    Rev 1.2   03/18/94 15:53:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:35:08   zvd
CD2 Modified so input is compatible with terminal I/O.
CD2 
CD2    Rev 1.0   01/20/94 10:21:42   pvcs
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
CD3   See name descriptors under Global Variables
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
CD4   Identifier      Type      Key   Description
CD4
CD4   cform           CHAR            File formats
CD4   cstats          CHAR            File status
CD4   iowork          CHAR            File use/type designator
CD4   nmfil           CHAR ARY        I/O file names:
CD4                                   1 - Control file name
CD4                             input 2 - Main input file name
CD4                             grid  3 - Coordinate input file name
CD4                             zone  4 - Zone input file name
CD4                             outpt 5 - Main output file name
CD4                             rstin 6 - Restart input file name
CD4                             rstot 7 - Restart output file name
CD4                             hist  8 - Simulation history output file name
CD4                             trac  9 - Solute history output file name
CD4                             cont  10 - Contour plot output file name
CD4                     dual or dpdp  11 - Dual porosity or dpdp contour plot
CD4                                          output file name
CD4                             stor  12 - Coefficient storage output file name
CD4                             check 13 - Input check output file name
CD4                             error 14 - Error output filename
CD4                             nopf  92 - Symbolic factorization filename
CD4   nmmax           INT             Number of I/O files
CD4   nufilb          CHAR            File numbers
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   null1            LOGICAL  Check for null lines or 0's.
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
CD5   numfil          INT ARY  Unit numbers for I/O files
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
CPS BEGIN cntlin
CPS 
CPS   FOR each I/O file
CPS   
CPS       set name to blank
CPS       read the file name from the control file (for read errors force . . .
CPS       . . . loop completion)
CPS       
CPS   END FOR
CPS   
CPS   FOR each I/O file (except error file)
CPS       IF the file name exists (not blank or zero)
CPS       
CPS          set the file unit number
CPS          
CPS          IF file status is 'old'
CPS             determine if file exists
CPS             IF file doesn't exist
CPS                IF this is the coefficient storage file
CPS                   set file status to new
CPS                ELSE   
CPS                   set error flag (to 2) and go to file assignments
CPS                END IF
CPS             END IF
CPS          END IF
CPS          
CPS          IF this isn't the save file or the coordinate and zone files . . .
CPS          . . . and they aren't the same as the input file
CPS             open and rewind the file
CPS          ELSE IF this is the coordinate or zone file and it is the same . . .
CPS          . . . as the input file
CPS             set the file unit number to the input file number
CPS          END IF
CPS          
CPS       ELSE
CPS       
CPS         set the file unit number to 0 (no file)
CPS         set the file name to 'not using'
CPS         
CPS       END IF
CPS       
CPS       IF file_i [inpt, iout] and the file number is zero 
CPS          set unit_i to standard in [inpt] or standard out [iout]
CPS          set file_i to terminal input / output
CPS       ELSE IF [incoor, inzone] and the file number is zero 
CPS          set unit_i to inpt
CPS          set file_i to input file
CPS       END IF
CPS   
CPS   END FOR
CPS   
CPS   call setunits to set I/O file unit numbers
CPS   
CPS END cntlin
CPS
C***********************************************************************

      use comai
      use comtable, only : lookup_file
      use comxi
      implicit none

      logical null1, opnd
      integer i, j, flen, numfil(32), itest, istat
      character*1 dummy1
      character*4 macro
      character*6 dummy
      character*100 filename
      character*120 input_string
c gaz 070720 need to zero out water table here because scanin called later      
      iwater_table = 0
      iair_table = 0
      ico2WH_table = 0
      itable_files = 0
      num_eos_table = 0
      numfil(1) = iocntl 

      do i = 2,nmmax-1
         nmfil(i) = ''
c         nmfil(i) = blank
      end do
c gaz 082521  think
      nmfil(33) = ''
      root_name = ''
c      root_name = blank
      read (iocntl, '(a6)') dummy
      backspace (iocntl)
      if (dummy(5:5) .eq. ":" .or. dummy(6:6) .eq. ":") then
!Read I/O file names using new format
         do
            read (iocntl, 1001, end = 101, err = 101) input_string
 1001       format (a120)
            if (null1(input_string)) exit
            j = 0
            do
               j = j + 1
! Exit if keyword delimiter ":" found
               if (input_string(j:j) .eq. ":" .or. j .gt. 6) exit
               dummy(j:j) = input_string(j:j)
            end do
            if (j .gt. 6) then
               write (ierr, 2000) nmfil(1), input_string
! Use unit 6 here as iptty not yet defined
               write (6, 2000) nmfil(1), input_string
               stop
            end if
            do
               j = j + 1
               if (input_string(j:j) .ne. " " .and.
     &              input_string(j:j) .ne. achar(9)) exit
            end do  
            flen = len_trim (input_string)
            filename = input_string(j:flen)
            
            read_names: select case (dummy(1:4))
            case ('inpu')
               nmfil(2) = filename
            case ('grid')
               nmfil(3) = filename
               if (dummy(5:5) .eq. 'a' .or. dummy(5:5) .eq. 'f') then
                  cform(3)='formatted'
               else if (dummy(5:5) .eq. 'u' .or. 
     &                 dummy(5:5) .eq. 'b') then
                  cform(3)='unformatted'
               end if
            case ('zone')
               nmfil(4) = filename
            case ('outp')
               nmfil(5) = filename
            case ('rsti')
               nmfil(6) = filename
            case ('rsto')
               nmfil(7) = filename
            case ('hist')
               nmfil(8) = filename
            case ('trac')
               nmfil(9) = filename
            case ('cont')
               nmfil(10) = filename
            case ('dual', 'dpdp')
               nmfil(11) = filename
            case ('stor')
               nmfil(12) = filename
            case ('chec')
               nmfil(13) = filename
            case ('erro')
               nmfil(14) = filename
            case ('root')
               root_name = filename
            case ('colu')
               nmfil(27) = filename
            case ('nopf')
               nmfil(28) = filename
            case ('co2i')
               nmfil(29) = filename
            case ('h2oi')
               nmfil(31) = filename
               iwater_table = 1
            case ('airi')
               nmfil(32) = filename  
               iair_table = 1
            case ('co2w')
               nmfil(33) = filename  
               ico2wh_table = 1               
            case ('look')
               lookup_file = filename
            case default
               write (ierr, 2000) nmfil(1), dummy
! Use unit 6 here as iptty not yet defined
               write (6, 2000) nmfil(1), dummy
               stop
            end select read_names
         end do
 2000    format ('Error in file input:', /, 100a, /, 'Check ', 100a)
 2001    format ('Error in file input:', /, 100a, /, 'Check ', a6,
     &        ' invalid file keyword')
      else
!Read I/O file names using old format
         do i=2,nmmax-1
            read (iocntl, 1000, end = 5, err = 5) nmfil(i)
 1000       format(a100)
 5       enddo
      end if
c gaz 111621 if no error file read in create 'fehmn.err'
      if(nmfil(14)(1:20).eq.'                    ') then
       ierr = 23
       nmfil(14)(1:9) = 'fehmn.err'
      endif    
 101  do i=2,nmmax-1
         if (.not. null1(nmfil(i))) then

            numfil(i) = nufilb(i)
            if (cstats(i) .eq. 'old    ') then
               inquire  (file = trim(nmfil(i)), exist = ex)
               if (.not. ex) then
                  if (iowork(i) .eq. 'isstor') then
                     cstats(i) = 'new    '
                  else
                     isw(i) = 2
                     go to 20
                  end if
               end if
            end if

            if (iowork(i) .ne. 'isave ' .and.
     *           .not. ((iowork(i) .eq. 'incoor' .or.
     *           iowork(i) .eq. 'inzone') .and.
     *           nmfil(i) .eq. nmfil(2))) then
               if (iowork(i) .eq. 'incoor') then
! Determine type of coordinate file (if not already determined for msim)
                  inquire (numfil(i), opened = opnd)
                  if (.not. opnd) then
                     open (numfil(i), file = nmfil(i), status=cstats(i),
     *                    form =  cform(i), err = 12)
                     if (cform(i) .eq. 'formatted') then
                        read (numfil(i), '(a4)', err=8, iostat=istat) 
     &                       macro
                        if (macro(1:3) .ne. 'fdm') then
                           read (numfil(i), *, err=8, iostat=istat) 
     &                          itest
                        end if
 8                      if (istat .ne. 0) then
! Check to see if file is unformatted
                           close (numfil(i))
                           cform(i) = 'unformatted'
                           open (numfil(i), file=nmfil(i), 
     *                          status = cstats(i), form = cform(i), 
     *                          err = 12)
                           read (numfil(i), err=12) macro
                           read (numfil(i), err=12) itest
                        end if
                     else if (cform(i) .eq. 'unformatted') then
                        read (numfil(i), err=10, iostat=istat) macro
                        read (numfil(i), err=10, iostat=istat) itest
 10                     if (istat .ne. 0) then
! Check to see if file is formatted
                           close (numfil(i))
                           cform(i) = 'formatted'
                           open (numfil(i), file=nmfil(i), 
     *                          status = cstats(i), form = cform(i), 
     *                          err = 12)
                           read (numfil(i), '(a4)', err=12, 
     &                       iostat=istat) macro
                           read (numfil(i), *, err=12, iostat=istat) 
     &                          itest
                        end if
                     end if
                     goto 23
                  end if
               else
                 open (numfil(i), file = nmfil(i), status = cstats(i),
     *                 form = cform(i), err = 12)
               end if
               goto 23
 12            isw(i)=2
               goto 20
 23            rewind  numfil(i)
            else if ((iowork(i) .eq. 'incoor' .or.
     *              iowork(i) .eq. 'inzone') .and.
     *              nmfil(i) .eq. nmfil(2)) then
               numfil(i) = numfil(2)
               cform(i) = cform(2)
            end if
            
         else
            numfil(i) = 0
            nmfil(i) = nmfily(3)
         end if

 20      if (iowork(i) .eq. 'inpt  ' .and. numfil(i) .eq. 0) then
               numfil(i) = 5
               nmfil(i) = nmfily(1)
         else if (iowork(i) .eq. 'incoor' .and. numfil(i) .eq. 0) then
               nmfil(i) = nmfil(2)
               numfil(i) = numfil(2)
               cform(i) = cform(2)
         else if (iowork(i) .eq. 'inzone' .and. numfil(i) .eq. 0) then
               nmfil(i) = nmfil(2)
               numfil(i) = numfil(2)
         end if

      end do

      call setunits (numfil)
         
      end
