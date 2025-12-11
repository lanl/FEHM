      subroutine cntlio (usub_num)
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
CD1 Manage the opening and closing of files using control file input.
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
CD2 $Log:   /pvcs.config/fehm90/src/cntlio.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Mon Jan 29 13:36:08 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   06/20/94 11:08:58   zvd
CD2   
CD2    Rev 1.3   03/18/94 15:53:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.2   02/17/94 10:19:32   zvd
CD2 Fixed bug when setting user subroutine call number.
CD2 
CD2    Rev 1.1   02/14/94 11:35:10   zvd
CD2 Modified so input is compatible with terminal I/O.
CD2 
CD2    Rev 1.0   01/20/94 10:21:44   pvcs
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
CD3   usub_num        INT      O    User subroutine number
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
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
CD4   cform           CHAR     faay   File formats
CD4   cstats          CHAR     faay   File status
CD4   cuser           CHAR     faax   String written when not using a
CD4                                     particular input/output file
CD4   ex              LOGICAL  faay   Logical for existence check (T/F)
CD4   iatty           INT      faai   Unit number for all tty output
CD4   iocntl          INT      faai   Unit number for control file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   iowork          CHAR     faay   File use/type designator
CD4   isw             INT      faay   Error flag for opening files,
CD4                                     isw = 1 -- No error
CD4                                     isw = 2 -- Error
CD4   nmfils          CHAR     faax   Control file name
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   cntlin                   Read names of files, set unit numbers and open
CD4                              files
CD4   writeio                  Write assigned file names, unit numbers, and
CD4                              file purpose to specified output unit                           
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
CD5   chk             CHAR     String form of user subroutine number 
CD5   i               INT      Loop index
CD5   nmfilc          CHAR     Terminal usage on/off 
CD5   nmfild          CHAR     Terminal usage on/off 
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
CPS BEGIN cntlio
CPS 
CPS   open the control file
CPS   call cntlin to read file names, set unit numbers and open files
CPS   
CPS   read tty usage
CPS   
CPS   IF all reference nodes should be printed
CPS      set unit numbers to 6 (iatty, iptty), identify all tty output
CPS   ELSE IF only first node should be printed
CPS      set unit numbers to 6 (iptty) and 0 (iatty), identify single . . .
CPS      . . . node tty output
CPS   ELSE 
CPS      identify as not using
CPS   END IF
CPS   
CPS   determine if user subroutine is to be called and subroutine . . .
CPS   . . . number to use in call
CPS   IF user subroutine number is not input
CPS      set number to 0 and set message about not using
CPS   ELSE
CPS      set number to input value
CPS      IF number is 0
CPS         set message about not using
CPS      ELSE
CPS         set subroutine call message
CPS      END IF
CPS   END IF
CPS   
CPS   IF tty output is being used
CPS      write tty output comment string
CPS      call writeio to write assigned file names, unit numbers, and . . .
CPS      . . . file purpose to tty output unit 
CPS   END IF
CPS   
CPS   initialize error existence flag to false
CPS   FOR each I/O file
CPS       IF there was an error opening the file
CPS          set error existence flag to true
CPS          write message about file causing error to error output
CPS          IF tty output is being used
CPS             write message about file causing error to tty
CPS          END IF
CPS       END IF
CPS   END FOR
CPS   
CPS   IF an error was detected
CPS      call writeio to write assigned file names, unit numbers, and . . .
CPS      . . . file purpose to error output unit 
CPS      write message about stopping job to error output
CPS      IF tty output is being used
CPS         write message about stopping job to error output
CPS      END IF
CPS      terminate program
CPS   END IF     
CPS 
CPS   close the control file
CPS 
CPS END cntlio
CPS
C***********************************************************************

      use comai
      use comxi
      implicit none

      integer i, usub_num

      character*4  chk
      character*60  comment

      open(iocntl, file = nmfil(1), status = cstats(1), form = cform(1))

      call cntlin

      chk = '    '
      read (iocntl, '(a4)', end = 10, err = 10) chk

 10   if (chk .eq. 'all ' .or. chk .eq. 'ALL '
     *        .or. chk .eq. 'a   ' .or. chk .eq. 'A   ') then
         iatty = 6
         iptty = 6
         comment = 'All reference output nodes will be written to tty'
      else if (chk .eq. 'some' .or. chk .eq. 'SOME'
     *           .or. chk .eq. 's   ' .or. chk .eq. 'S   ') then
         iatty = 0
         iptty = 6
         comment = 
     .        'First reference output node will be written to tty'
      else
         comment = 'Not using tty output'
      end if

      chk = '    '
      read (iocntl, '(a4)', end = 20, err = 20) chk

 20   if (chk .eq. '    ' .or. chk .eq. 'none' .or. chk .eq. 'NONE' 
     *        .or. chk .eq. 'n   ' .or. chk .eq. 'N   ') then
         usub_num = 0
         cuser =  'not using'
      else
         read(chk, *) usub_num 
         if (usub_num .ne. 0)  then
            write(cuser , '(a4)') chk
         else
            cuser =  'not using'
         end if
      end if

      if (iptty .ne. 0) then
         write (iptty, '(/, 1x, 60a)') comment
         call writeio (iptty)
      end if

      ex = .FALSE.
      do i = 2, 12
         if (isw(i) .eq. 2) then
            ex = .TRUE.
            write(ierr, 6105) iowork(i)
            if (iptty .ne. 0) write(iptty, 6105) iowork(i) 
 6105       format(/, ' **** Error opening file ', a6, ' ****', /)
         end if
      end do

      if (ex) then
         call  writeio (ierr)
         write(ierr, 6110)
         if (iptty .ne. 0) write(iptty, 6110)

 6110    format( ' ****---------------------------****', /,
     *           ' ****       JOB STOPPED         ****', /,
     *           ' ****---------------------------****', /)
         stop
      end if

      close (iocntl)

      end
