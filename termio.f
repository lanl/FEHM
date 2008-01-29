      subroutine termio (usub_num)
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
CD1 Manage the opening and closing of files using terminal input.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 29-OCT-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/termio.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:16:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Fri Feb 02 12:46:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   12/01/94 10:45:08   llt
CD2 Changed spacing on format, for ibm
CD2 
CD2    Rev 1.5   11/29/94 18:22:42   llt
CD2 Changed length of jdate to 11 characters for ibm
CD2 
CD2    Rev 1.4   06/20/94 11:09:12   zvd
CD2  
CD2 
CD2    Rev 1.3   03/18/94 15:56:00   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.2   02/17/94 10:18:42   zvd
CD2 Fixed bug when setting user subroutine call number.
CD2 
CD2    Rev 1.1   02/14/94 11:37:22   zvd
CD2 Modified so terminal I/O is compatible with control file input. Added single
CD2 character response acceptance for terminal queries.
CD2 
CD2    Rev 1.0   01/20/94 10:28:34   pvcs
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
CD4   cuser           CHAR     faax   String written for user subroutinenumber
CD4                                     use
CD4   iatty           INT      faai   Unit number for selected tty output
CD4   iocntl          INT      faai   Unit number for control file
CD4   iptty           INT      faai   Unit number for all tty output
CD4   verno           CHAR     faac   Contains version number of FEHMN code
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   close_files              Close all open files
CD4   cntlio                   Manage the opening and closing of files using
CD4                              control file input
CD4   termin                   Read I/O file names from terminal, set unit
CD4                              numbers and open files
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
CD5   comment         CHAR     Terminal output on/off comment
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
CPS BEGIN termio
CPS 
CPS   write terminal input initiation messages
CPS   call termin to read file names from the terminal set unit numbers 
CPS   and open files
CPS   
CPS   IF iocntl has been set to 1 (an input control file is present)
CPS      call cntlio to manage opening and closing of files specified in 
CPS      a control file
CPS   ELSE
CPS   
CPS      determine if tty output is wanted
CPS      IF all reference nodes will be written
CPS         set tty unit values to 6 (iatty, iptty) and message about tty use
CPS      ELSE IF selected reference nodes will be written
CPS         set tty unit values (iatty to 0, iptty to 6) and message 
CPS         about selected tty output
CPS      ELSE none
CPS         set tty unit values to 0 (iatty, iptty) and message about not 
CPS         using tty 
CPS      END IF
CPS      
CPS      determine if user subroutine is to be called and subroutine 
CPS      number to use in call
CPS      IF user subroutine number is not input
CPS         set number to 0 and set message about not using
CPS      ELSE
CPS         set number to input value
CPS         IF number is 0
CPS            set message about not using
CPS         ELSE
CPS            set subroutine call message
CPS         END IF
CPS      END IF
CPS      
CPS      write tty usage message to terminal
CPS      call writeio to display values set on terminal
CPS      determine if values are ok and program execution should continue
CPS      IF values are not ok
CPS         call close_files to close all files that were opened and 
CPS         return to beginning of routine to reset selections
CPS      ELSE IF the program should be stopped
CPS         write message about stopping program and terminate program
CPS      END IF
CPS      
CPS   END IF
CPS 
CPS END termio
CPS
C***********************************************************************

      use comai
      use comxi
      implicit none

      integer usub_num
      character*4  chk
      character*60 comment

c**** read and write data ****
      write (6, 6000) verno, jdate, jtime
 6000 format(/, 1x, 'version  ', a30, 5x, a11, 1x, a8,
     *     //, 1x, '**** Default names for I/O files ****',
     *     //, 5x,'control file                          : fehmn.files',
     *     /,  5x,'input file                            : filen.*',
     *     /,  5x,'geometry data file                    : filen.*',
     *     /,  5x,'zone data file                        : filen.*',
     *     /,  5x,'output file                           : filen.out',
     *     /,  5x,'read file (if it exists)              : filen.ini',
     *     /,  5x,'write file (if it exists)             : filen.fin',
     *     /,  5x,'history plot file                     : filen.his',
     *     /,  5x,'tracer history plot file              : filen.trc',
     *     /,  5x,'contour plot file                     : filen.con',
     *     /,  5x,'dual or dpdp contour plot file        : filen.dp',
     *     /,  5x,'stiffness matrix data read/write file : filen.stor',
     *     /,  5x,'input check file                      : filen.chk',
     *     /)

      write (6, 6001)
 6001 format(1x, '**** where ****',
     *     /, 5x, '"filen.*" may be 100 characters maximum.  ',
     *     'If a name is not entered', /, 5x, 'when prompted for, ',
     *     'a default file name is used. "fehmn.dat" is the',
     *     /, 5x, 'default used for the input file name.',
     *     //, 1x, '**** note ****',
     *     /, 5x, 'A save file and input check file are always',
     *     ' written, if you do not', /, 5x, 'provide a name for these',
     *     ' files, the following defaults will be used: ', 
     *     /, 5x, 'fehmn.fin, fehmn.chk', //)

 10   continue

      call termin

      if (iocntl .eq. 1)  then
         call cntlio (usub_num)
      else

c**** Check if tty output wanted ****
         write(6, 6010)
 6010    format(/,1x,'tty output -- show all reference nodes,',
     *        ' selected reference nodes, or none:', /, 1x,
     *        '[(all/some/none), RETURN = none]')

         chk = '    '
         read (5, '(a4)') chk
    
         if (chk .eq. 'all ' .or. chk .eq. 'ALL '
     *        .or. chk .eq. 'a   ' .or. chk .eq. 'A   ') then
            iatty = 6
            iptty = 6
            comment = 
     *           'All reference output nodes will be written to tty'
         else  if (chk .eq. 'some' .or. chk .eq. 'SOME'
     *           .or. chk .eq. 's   ' .or. chk .eq. 'S   ')  then
            iatty = 0
            iptty = 6
            comment = 
     *           'First reference output node will be written to tty'
         else
            iatty = 0
            iptty = 0
            comment = 'Not using tty output'
         end if
     
c**** user's subroution user call step of number ****

         write(6, 6012)
 6012    format(/, 1x, 'user subroutine number (provided to ',
     *        'subroutine USER before every time step):', /, 1x,
     *        '[RETURN = none]')
     
         chk = '    '
         read (5, '(a4)') chk

         if (chk .eq. '    ' .or. chk .eq. 'none' .or. chk .eq. 'NONE' 
     *        .or. chk .eq. 'n   ' .or. chk .eq. 'N   ') then
            usub_num = 0
            cuser =  'not using'
         else
            read(chk, '(i4)') usub_num 
            if (usub_num .ne. 0)  then
               write(cuser , '(a4)') chk
            else
               cuser =  'not using'
            end if
         end if

         write (6, '(/, 1x, 60a)') comment
         call writeio (6)

         write(6, 6110)
 6110    format(//, 1x, 'If data is OK enter yes to continue, no to',
     *        ' restart terminal input,', /, 1x, 'or stop to end ',
     *        'program: [(yes/no/stop), RETURN = yes]', /)
         
         read (5,'(a4)')  chk

         if (chk .eq. 'no  ' .or. chk(1:1) .eq. 'NO  '
     *        .or. chk .eq. 'n   ' .or. chk(1:1) .eq. 'N   ')  then

            call close_files
            write(6, 6002)
 6002      format(////,1x,'**** restart terminal input, ',
     *           'reenter file names ****',//)
            go to 10

         else if (chk .eq. 'stop' .or. chk .eq. 'STOP'
     *           .or. chk .eq. 's   ' .or. chk .eq. 'S   ')  then
     
            write(6, 6120)
 6120       format(/,1x,'**** job stopped * set data again ****',
     *           /,1x,'****------------------------------****',/)

            stop
         end if
      end if

      return
      end
