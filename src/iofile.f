      subroutine iofile (usub_num)
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
CD1 Manage the opening of input and output files.
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
CD2 $Log:   /pvcs.config/fehm90/src/iofile.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Wed Jun 12 15:09:52 1996   zvd
CD2 Modified order of write to iout to compensate for writes from 
CD2 optional input file routine
CD2 
CD2    Rev 1.7   Tue Jan 30 15:24:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   09/08/95 11:44:46   zvd
CD2 
CD2    Rev 1.5   08/09/95 12:07:12   llt
CD2 changed ierr from unit 0 to unit 23 (for IBM)
CD2 
CD2    Rev 1.4   11/29/94 18:22:14   llt
CD2 Changed length of jdate to 11 characters for ibm
CD2 
CD2    Rev 1.3   06/20/94 11:09:02   zvd
CD2  
CD2    Rev 1.2   03/18/94 15:55:28   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:32:38   zvd
CD2 Correct misnumbered file designator.
CD2 
CD2    Rev 1.0   01/20/94 10:25:24   pvcs
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
CD4   iout            INT      faai   Unit number for output file
CD4   isave           INT      faai   ?Unit number for restart file (to write)
CD4   ischk           INT      faai   Unit number for input data check file
CD4   iscon           INT      faai   Unit number for contour data file
CD4   iscon1          INT      faai   Unit number for dual porosity or dpdp
CD4                                     contour data file
CD4   ishis           INT      faai   Unit number for history data file
CD4   istrc           INT      faai   Unit number for tracer history data file
CD4   jdate           CHAR     faac1  Contains the date (mm/dd/yr)
CD4   jtime           CHAR     faac1  Contains the time (hr:mn:sc)
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
!D4                                   14 - Error output file name
!D4                                   15 - Pest output file name
!D4                                   16 - Pest auxilliary output file name
!D4                                   17 - Streamline particle tracking
!D4                                          output file 1 name
!D4                                   18 - Streamline particle tracking
!D4                                          output file 2 name
!D4                                   19 - Streamline particle tracking
!D4                                          output file 3 name
!D4                                   20 - Streamline particle tracking
!D4                                          output file 4 name
!D4                                   21 - Streamline particle tracking
!D4                                          output file 5 name
!D4                                   22 - Streamline particle tracking
!D4                                          output file 6 name
!D4                                   23 - Streamline particle tracking
!D4                                          output file 7 name
!D4                                   24 - Extracted flow model (flow 
!D4                                          macro) file name 
!D4                                   25 - Streamline particle tracking
!D4                                          output file 8 name
!D4                                   26 - Streamline particle tracking
!D4                                          output file 9 name
CD4   verno           CHAR     faac   Contains version number of FEHMN code
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   cntlin                   Read control file for file name input
CD4   termin                   Invoke terminal I/O for file name input
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
CD5   ex              LOGICAL  Boolean for existence of file
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
CPS BEGIN iofile
CPS 
CPS   initialize I/O file default values (names, status, format, etc.)
CPS   open output error file and write code version, date and time to 
CPS    file
CPS   call inquire to determine if input control file exists
CPS   
CPS   IF the control file does not exist
CPS      call termin to invoke terminal input for files used
CPS   ELSE
CPS      call cntlin to read control file for files used
CPS   END IF
CPS   
CPS   IF save file is not being used
CPS      set up and open miscellaneous file 
CPS   END IF   
CPS   
CPS   [for each output file (iout, ishis, istrc, iscon, iscon1, ischk)]
CPS   IF the output file is open 
CPS      write code version, date and time to file
CPS   END IF
CPS 
CPS END iofile
CPS
C***********************************************************************

      use comai
      use comxi
      use comuserc, only : in
      implicit none

      logical opnd
      character*200 cmdline
      integer iargc, len, nargc, usub_num

      isw=1
      nufilb(1)=1
      nufilb(2)=11
      nufilb(3)=12
      nufilb(4)=13
      nufilb(5)=14
      nufilb(6)=15
      nufilb(7)=16
      nufilb(8)=17
      nufilb(9)=18
      nufilb(10)=19
      nufilb(11)=20
      nufilb(12)=21
      nufilb(13)=22
      nufilb(14)=23
      nufilb(15)=24
      nufilb(16)=25
      nufilb(17)=26
      nufilb(18)=27
      nufilb(19)=28
      nufilb(20)=29
      nufilb(21)=30
      nufilb(22)=31
      nufilb(23)=32
      nufilb(24)=33
      nufilb(25)=34
      nufilb(26)=35
      nufilb(27)=36
      nufilb(28)=37
      nufilb(29)=38
      nufilb(30)=39 
      nufilb(31)=40 
c gaz 072220 added unit nunber for air EOS table      
      nufilb(32)=41
      suffix(1)='.files'
      suffix(2)='.dat'
      suffix(3)='.dat'
      suffix(4)='.dat'
      suffix(5)='.out'
      suffix(6)='.ini'
      suffix(7)='.fin'
      suffix(8)='.his'
      suffix(9)='.trc'
      suffix(10)='.con'
      suffix(11)='.dp'
      suffix(12)='.stor'
      suffix(13)='.chk'
      suffix(14)='.err'
      suffix(15)='.pest'
      suffix(16)='.pest1'
      suffix(17)='.sptr1'
      suffix(18)='.sptr2'
      suffix(19)='.sptr3'
      suffix(20)='.sptr4'
      suffix(21)='.sptr5'
      suffix(22)='.sptr6'
      suffix(23)='.sptr7'
      suffix(25)='.sptrs'
      suffix(26)='.sptrx'
      suffix(24)='.subbc'
      suffix(27)='.col'
      suffix(28)='.nop'
      suffix(29)='.txt'
      suffix(30)='.well2'
      suffix(31)='.txt'
      suffix(32)='.txt'
      iowork(1)='iocntl'
      iowork(2)='inpt  '
      iowork(3)='incoor'
      iowork(4)='inzone'
      iowork(5)='iout  '
      iowork(6)='iread '
      iowork(7)='isave '
      iowork(8)='ishis '
      iowork(9)='istrc '
      iowork(10)='iscon '
      iowork(11)='iscon1'
      iowork(12)='isstor'
      iowork(13)='ischk '
      iowork(14)='ierr  '
      iowork(15)='ispest'
      iowork(16)='ispst1'
      iowork(17)='isptr1'
      iowork(18)='isptr2'
      iowork(19)='isptr3'
      iowork(20)='isptr4'
      iowork(21)='isptr5'
      iowork(22)='isptr6'
      iowork(23)='isptr7'
      iowork(25)='isptr8'
      iowork(26)='isptr9'
      iowork(24)='isubm '
      iowork(27)='iswt  '
      iowork(28)='ionop '
      iowork(29)='ioco2 '
      iowork(30)='well2 ' 
      iowork(31)='ioh2o '   
      iowork(32)='ioair '   
      cstats(1)='old    '
      cstats(2)='old    '
      cstats(3)='old    '
      cstats(4)='old    '
      cstats(5)='unknown'
      cstats(6)='old    '
      cstats(7)='unknown'
      cstats(8)='unknown'
      cstats(9)='unknown'
      cstats(10)='unknown'
      cstats(11)='unknown'
      cstats(12)='old    '
      cstats(13)='unknown'
      cstats(14)='unknown'
      cstats(15)='unknown'
      cstats(16)='unknown'
      cstats(17)='unknown'
      cstats(18)='unknown'
      cstats(19)='unknown'
      cstats(20)='unknown'
      cstats(21)='unknown'
      cstats(22)='unknown'
      cstats(23)='unknown'
      cstats(24)='unknown'
      cstats(25)='unknown'
      cstats(26)='unknown'
      cstats(27)='unknown'
      cstats(28)='unknown'
      cstats(29)='old    '
      cstats(30)='unknown'  
      cstats(31)='old    '   
c gaz 072220 air table modification      
      cstats(32)='old    '   
      cform(1)='formatted'
      cform(2)='formatted'
! Coordinate file can be formatted or unformatted
      cform(3)='formatted'
      cform(4)='formatted'
      cform(5)='formatted'
      cform(6)='formatted'
      cform(7)='formatted'
      cform(8)='formatted'
      cform(9)='formatted'
      cform(10)='formatted'
      cform(11)='formatted'
      cform(12)='formatted'
      cform(13)='formatted'
      cform(14)='formatted'
      cform(15)='formatted'
      cform(16)='formatted'
      cform(17)='formatted'
      cform(18)='formatted'
      cform(19)='formatted'
      cform(20)='formatted'
      cform(21)='formatted'
      cform(22)='formatted'
      cform(23)='formatted'
      cform(24)='formatted'
      cform(25)='formatted'
      cform(26)='formatted'
      cform(27)='formatted'
      cform(28)='unformatted'
      cform(29)='formatted'
      cform(31)='formatted'  
      cform(32)='formatted'   
      blank=' '
      if(in(4).NE.666) nmfil( 1) = 'fehmn.files'
      nmfil( 2) = 'fehmn.dat'
      nmfil( 3) = 'fehmn.dat'
      nmfil( 4) = 'fehmn.dat'
      nmfil( 5) = 'fehmn.out'
      nmfil( 6) = 'fehmn.ini'
      nmfil( 7) = 'fehmn.fin'
      nmfil( 8) = 'fehmn.his'
      nmfil( 9) = 'fehmn.trc'
      nmfil(10) = 'fehmn.con'
      nmfil(11) = 'fehmn.dp'
      nmfil(12) = 'fehmn.stor'
      nmfil(13) = 'fehmn.chk'
      nmfil(14) = 'fehmn.err'
      nmfil(15) = 'fehmn.pest'
      nmfil(16) = 'fehmn.pest1'
      nmfil(17) = 'fehmn.sptr1'
      nmfil(18) = 'fehmn.sptr2'
      nmfil(19) = 'fehmn.sptr3'
      nmfil(20) = 'fehmn.sptr4'
      nmfil(21) = 'fehmn.sptr5'
      nmfil(22) = 'fehmn.sptr6'
      nmfil(23) = 'fehmn.sptr7'
      nmfil(25) = 'fehmn.sptrs'
      nmfil(26) = 'fehmn.sptrx'
      nmfil(24) = 'fehmn.subbc'
      nmfil(27) = ''
      nmfil(28) = 'nop.temp'
      nmfil(29) = 'co2_interp_table.txt'
      nmfil(30) = 'fehmn.well2' 
c gaz 122020   changed 31 and 32 (EOS tables) to blank   
      nmfil(31) = '' 
      nmfil(32) = '' 
      nmfily( 1) = 'terminal console input'
      nmfily( 2) = 'terminal console output'
      nmfily( 3) = 'not using'

      ex = .false.
      cmdline = ''
! Has the control file name been provided on the command line

      if(in(4).NE.666) then 
         nargc = iargc()
         if (nargc .ge. 1) then
            call getarg(1, cmdline)
            len = len_trim(cmdline)
            inquire (file = cmdline(1:len), exist = ex)
            if (ex) nmfil(1) =  cmdline(1:len)
         end if
      end if

! If the control file name was not entered on the command line,
! Is there an fehmn control file named fehmn.files?
      if (.not. ex) inquire (file = nmfil(1), exist = ex)

      if (.not. ex) then
         nmfil( 1) = nmfily(3)
         call termio (usub_num)
      else
         iocntl = nufilb(1)
         call cntlio (usub_num)
      end if

! Open error output file [if not already open for msim] 
      if(in(4).NE.666) ierr = nufilb(14)
      if(nmfil(14)(1:9).eq.'not using') go to 444
       inquire (ierr, opened = opnd)
       if (.not. opnd) then
         open (ierr, file = nmfil(14), status = cstats(14),
     *        form = cform(14))
         write (ierr, 1000)  verno, jdate, jtime
       end if
444    continue       

! No longer create default save or fin files

! Write code version, date and time to open output files  
      if (iout   .gt. 0) call write_copyright (iout)
      if (ishis  .gt. 0) write(ishis , 1000)  verno, jdate, jtime
      if (istrc  .gt. 0) write(istrc , 1000)  verno, jdate, jtime
      if (iscon  .gt. 0) write(iscon , 1000)  verno, jdate, jtime
      if (iscon1 .gt. 0) write(iscon1, 1000)  verno, jdate, jtime
      if (ischk  .gt. 0) write(ischk , 1000)  verno, jdate, jtime
 1000 format(a30, 3x, a11, 3x, a8)
      
      end
