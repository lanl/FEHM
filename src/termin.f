      subroutine termin
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
CD1 Read I/O file names from terminal, set unit numbers and open files.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 22-NOV-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/termin.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Fri Feb 02 12:45:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   06/20/94 11:09:10   zvd
CD2  
CD2 
CD2    Rev 1.5   03/18/94 15:55:58   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.4   02/28/94 11:59:44   zvd
CD2 Edited Rev 1.3 comment to reflect actual update.
CD2 
CD2    Rev 1.3   02/28/94 11:52:38   zvd
CD2 Corrected problem of unit assignment for coordinate and zone 
CD2 input files.
CD2 
CD2    Rev 1.2   02/17/94 10:17:20   zvd
CD2 Fixed bug when opening alternate zone or coordinate data files.
CD2 
CD2    Rev 1.1   02/14/94 11:37:20   zvd
CD2 Modified so terminal I/O is compatible with control file input. 
CD2 Added single character response acceptance for terminal queries.
CD2 
CD2    Rev 1.0   01/20/94 10:28:32   pvcs
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
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   cform           CHAR ARY faax   I/O file formats
CD4   cstats          CHAR ARY faax   I/O file status
CD4   ex              LOGICAL  faay   Logical for existence check (T/F)
CD4   iowork          CHAR ARY faax   File use/type designator
CD4   itert           INT      faai   Intermediate iteration counter
CD4   nmfil           CHAR ARY faax   I/O file names:
CD4                                   1 - Control file name
CD4                                   2 - Main input file name
CD4                                   3 - Coordinate input file name
CD4                                   4 - Zone input file name
CD4                                   5 - Main output file name
CD4                                   6 - Restart input file name
CD4                                   7 - Restart output file name
CD4                                   8 - Simulation history output file name
CD4                                   9 - Solute history output file name
CD4                                   10 - Contour plot output file name
CD4                                   11 - Dual porosity or dpdp contour plot
CD4                                          output file name
CD4                                   12 - Coefficient storage output file name
CD4                                   13 - Input check output file name
CD4                                   14 - Error output file name
CD4   nmfily          CHAR ARY faax   File messages
CD4   nmmax           INT      faay   Number of I/O files
CD4   nufilb          CHAR ARY faax   File unit numbers
CD4   suffix          CHAR ARY faax   Default file suffixes
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   prefix                   Determine prefix of input file name
CD4   setunits                 Set file unit numbers
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
CD5   namfil          CHAR     Input file name
CD5   numfil          INT ARY  Unit numbers for I/O files
CD5   yesno           CHAR     String used for input yes or no response
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
CPS BEGIN termin
CPS 
CPS   initialize yesno flag to 'n'
CPS   
CPS   FOR each I/O file (except error file)
CPS   
CPS       set file unit number
CPS       
CPS       IF this isn't the control or input file and the file name 
CPS        prefix is different from the input file name prefix
CPS          IF this is the coordinate or zone file
CPS             set default file name to input file name
CPS          ELSE
CPS             set default file name using on input file name prefix and 
CPS              default suffix
CPS          END IF
CPS       END IF
CPS       
CPS       IF the file name should be set (flag is 'n')
CPS       
CPS          initialize input file string to blank
CPS          determine file name or if using default or if not using 
CPS          
CPS          IF default input is to be used
CPS             set input file string to default
CPS          END IF
CPS          
CPS          IF not using file
CPS          
CPS             set input file string to 'not using'
CPS             set file unit number to 0
CPS         
CPS          ELSE IF file status is old 
CPS          
CPS             determine if file exists
CPS             IF the file doesn't exist
CPS                IF this is the coefficient storage file
CPS                   set file status to new
CPS                ELSE
CPS                   write message to screen and request another name
CPS                   return to file name query
CPS                END IF
CPS             END IF
CPS             
CPS             IF this is the control file
CPS             
CPS                set control file unit number and name
CPS                return to calling routine to invoke control file I/O
CPS                
CPS             ELSE IF this is the input file
CPS             
CPS                set input unit number and name
CPS                
CPS                call file_prefix to determine the input file name prefix
CPS                IF prefix is longer than 94 characters               
CPS                   warn user that prefix was truncated to 94 
CPS                    charcaters and set prefix length tho 94
CPS                END IF
CPS                
CPS                determine if all file names should use the input file prefix
CPS                IF flag isn't set to 'n'
CPS                   set flag to 'y'
CPS                END IF
CPS                
CPS             END IF
CPS             
CPS          END IF

CPS       ELSE
CPS       
CPS          IF this is the coordinate or zone input file
CPS             set file unit and name equal to the input file
CPS          ELSE   
CPS             use the default file name 
CPS          END IF
CPS          
CPS          IF this is the initialization file
CPS             determine if file existss
CPS             IF the file does not exist
CPS                set file unit to 0
CPS             END IF
CPS          ELSE IF this is the coefficient storage file
CPS             determine if file existss
CPS             IF the file does not exist
CPS                set file status to new
CPS             END IF
CPS          END IF
CPS          
CPS       END IF
CPS       
CPS       IF the file unit number is 0
CPS       
CPS          IF this is the input file
CPS             set unit to tty (5) and name to 'terminal'
CPS          ELSE if this is the coordiante or zone input file
CPS             set file unit and name equal to the input file
CPS          ELSE if this is the output file
CPS             set unit to tty (6) and name to 'terminal'
CPS          ELSE
CPS             set input file name string to 'not using'
CPS          END IF
CPS          
CPS       ELSE   
CPS       
CPS          IF this isn't the save file or the coordinate and zone files 
CPS           and they aren't the same as the input file
CPS             open and rewind the file
CPS          END IF
CPS       
CPS       END IF
CPS       
CPS       set file name to input file string
CPS   
CPS   END FOR
CPS   
CPS   call setunits to assign file unit numbers
CPS 
CPS END termin
CPS
C***********************************************************************

      use comai
      use comxi
      implicit none
c gaz 052222 changed dim of numfil from 13 to 14
      integer i, numfil(14)
      character*1 yesno
      character*100 namfil
      
      yesno = 'n'
C     Make sure coefficient storage file starts out with status of old
      cstats(12) = 'old    ' 
      
      do i = 1, nmmax - 1
         
C     **** Append file suffix to prefix from input file name ****
         if (i .gt. 2 .and. 
     *        nmfil(2)(1:itert) .ne. nmfil(i)(1:itert)) then 
            nmfil(i) = blank
            if (i .eq. 3 .or. i .eq. 4) then
               nmfil(i) = nmfil(2)
            else
               nmfil(i)(1:itert) = nmfil(2)(1:itert)
               nmfil(i)(itert+1:itert+5) = suffix(i)
            end if
         endif
         
         numfil(i) = nufilb(i)
         
         if (yesno .eq. 'n') then 
            namfil = blank
 10         write(6, 6000)  iowork(i), nmfil(i)
 6000       format(/, 1x, 'Enter name for ', a6, 1x, 
     *           '-- default file name: ', a100, /, 1x,
     *           '[(name/na or not using), RETURN = DEFAULT]')
            
            read (5, '(a100)')  namfil
            if (namfil(1:4) .eq. '    ') then
               namfil = nmfil(i)
            end if
            if (namfil(1:4) .eq. 'na  ' .or. 
     *           namfil(1:4) .eq. 'not ') then
               namfil = nmfily(3)
               numfil(i) = 0
            else if (cstats(i) .eq. 'old    ')  then
               inquire (file = namfil, exist = ex)
               if (.not. ex) then
                  if (iowork(i) .eq. 'isstor') then
                     cstats(i) = 'new    '
                  else
                     write(6,*) iowork(i), ' Input file does not exist,'
     *                    , ' enter another'
                     go to 10
                  endif
               endif
               if (iowork(i) .eq. 'iocntl') then
                  iocntl = numfil(i)
                  nmfil(i) = namfil
                  return
               else if (iowork(i) .eq. 'inpt  ') then  
                  inpt = numfil(i)
                  nmfil(i) = namfil
                  call file_prefix(namfil, itert)
                  if (itert .gt. 94) then
                     itert = 94
                     write (6, *) ' **** Input file name prefix was',
     *                    ' truncated to 94 characters ****'
                  end if
                  write (6, *)
                  write(6, *) ' Do you want all file names of the form '
     *                 , namfil(1:itert), '.*',' ? [(y/n), RETURN = y]',
     *                 ' *** Note: If "y" incoor and inzone will equal',
     *                 ' inpt ***'
                  
                  read(5, '(a1)') yesno(1:1)
                  if (yesno(1:1) .ne. 'n') yesno = 'y'
               end if
            end if
         else
C     **** Set incoor and inzone to inpt ****
            if (iowork(i) .eq. 'incoor' .or. 
     *           iowork(i) .eq. 'inzone') then
               namfil = nmfil(2)
               nmfil(i) = nmfil(2)
               numfil(i) = numfil(2)
            else 
               namfil = nmfil(i)
            end if
C     **** Set iread to 'not using' if file not found ****
            if (iowork(i) .eq. 'iread ') then
               inquire (file = namfil, exist = ex)
               if (.not. ex) then
                  numfil(i) = 0
               end if
C     **** Set isstor to status new if file not found ****
            else if (iowork(i) .eq. 'isstor') then
               inquire (file = namfil, exist = ex)
               if (.not. ex) then
                  cstats(i) = 'new    '
               end if
            end if
         end if
         
         if (numfil(i) .eq. 0)  then
            if (iowork(i) .eq. 'inpt  ')  then
               numfil(i) = 5
               namfil = nmfily(1)
            else if (iowork(i) .eq. 'incoor' .or. 
     *              iowork(i) .eq. 'inzone')  then
               numfil(i) = numfil(2)
               namfil = nmfil(2)
            else if (iowork(i) .eq. 'iout  ')  then
               numfil(i) = 6
               namfil = nmfily(2)
            else
               namfil = nmfily(3)
            end if
         else  if (numfil(i) .eq. 5)  then
            namfil = nmfily(1)
         else if (numfil(i) .eq. 6)  then
            namfil = nmfily(2)
         else 
            if (iowork(i) .ne. 'isave ' .and. 
     *           .not. ((iowork(i) .eq. 'incoor' .or. 
     *           iowork(i) .eq. 'inzone') .and.
     *           namfil .eq. nmfil(2))) then
               open (numfil(i), file = namfil, status = cstats(i),
     *              form = cform (i))
               rewind  numfil(i)
            else if ((iowork(i) .eq. 'incoor' .or.
     *              iowork(i) .eq. 'inzone') .and.
     *              nmfil(i) .eq. nmfil(2)) then
               numfil(i) = numfil(2)
            end if
         end if
         nmfil(i) = namfil
      end do   
      
      call setunits (numfil)

      return
      end
