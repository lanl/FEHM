      subroutine intime
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
CD1 Read time step input.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 22-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/intime.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Wed Jun 12 15:21:10 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.8   Mon Jun 03 11:18:10 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.7   Fri May 31 10:49:32 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.6   Mon May 20 15:43:40 1996   hend
CD2 Fixed parser usage to work for real input as integer
CD2 
CD2    Rev 1.5   Fri May 17 13:00:06 1996   hend
CD2 Added optional parameter to time macro for initial time
CD2 
CD2    Rev 1.4   Fri Feb 16 09:45:46 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.3   Tue Jan 30 15:19:06 1996   hend
CD2 Changed Purpose and added Requirements Traceability
CD2 
CD2    Rev 1.2   08/16/95 16:28:14   robinson
CD2 Changed name of variable to set print out interval
CD2 
CD2    Rev 1.1   03/18/94 16:03:50   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:22   pvcs
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
CD3   macro           CHAR     I    Macro data being read
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
CD4   day             REAL*8   faar   Current time step size in days
CD4   deef            REAL*8   fcc    Used for temporary storage of time step
CD4                                     size
CD4   depf            REAL*8   fcc    Used for temporary storage of time step
CD4                                     multiplier
CD4   dit             REAL*8   fdd1   Array containing time step changes
CD4   icgts           INT      faai   Parameter controlling the time of solution
CD4                                     parameter changes
CD4   iprtout         INT      faai   Print-out interval, number of time steps
CD4   inpt            INT      faai   Unit number for input file
CD4   itc             INT      fdd1i  New print interval for time step changes
CD4   iyear           INT      faai   Current year in simulation
CD4   month           INT      faai   Current month in simulation
CD4   nstep           INT      faai   Maximum number of time steps
CD4   tims            REAL*8   faar   Ending simulation time
CD4   wdd1            CHAR     faac   Alternate character input string
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   null1            LOGICAL  Check for null lines or 0's in lines
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
CD5   itci            INT      New print interval
CD5   diti            REAL*8   Time to change time step
CD5   ditip           REAL*8   Time step size
CD5   ditip2          REAL*8   Time step multiplier
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
CD9 2.5.1 Implement time-step mechanism
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9 2.8   Provide Multiple Realization Option
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
CPS BEGIN intime
CPS 
CPS   read initial time step information  
CPS   
CPS   REPEAT
CPS     read input line
CPS   EXIT IF line is null
CPS     reread input line using time step change format
CPS     IF time change is zero
CPS        exit the loop
CPS     END IF
CPS     set time step values
CPS   UNTIL all time step change data is read
CPS   
CPS   If time step changes were input
CPS      FOR each time step change entered
CPS          set 
CPS      END FOR
CPS   END IF
CPS   
CPS END intime
CPS
C***********************************************************************

      use comci
      use comdi
      use comdti
      use comai
      implicit none

      logical null1
      integer i,itci
      real*8 diti,ditip,ditip2
      character*80 input_msg
      integer msg(10)
      integer nwds
      real*8 xmsg(10)
      integer imsg(10),i2,index
      character*32 cmsg(10)

      read (inpt,'(a)') input_msg
      call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
      if (msg(1).eq.1) then
         day=imsg(1)
      else
         day=xmsg(1)
      endif
      if (msg(2).eq.1) then
         tims=imsg(2)
      else
         tims=xmsg(2)
      endif
      enday=tims
      nstep=imsg(3)
      iprtout=imsg(4)
      iyear=imsg(5)
      month=imsg(6)
      if (nwds.ge.7) then
         irsttime=1
         if (msg(7).eq.1) then
            rsttime=imsg(7)
         else
            rsttime=xmsg(7)
         endif
      else
         irsttime=0
      endif
      
! icgts dtermined in scanin
!      icgts=0

c Read first time to see how many "dits" if any exist
C icgts determined in scanin now
!      ditnumlines=0
! 100  continue
!      read (inpt,'(a80)') wdd1
!      if (null1(wdd1)) go to 105
!      ditnumlines=ditnumlines+1
!      goto 100
      
 105  if (icgts.eq.0) goto 99
!      do i2=1,ditnumlines+2
!         backspace inpt
!      enddo
!      icgts=ditnumlines
      if(irun.ne.1) deallocate(dit,itc)
c     Add max time step values to the end of dit array
c     BAR 8-26-98
c      allocate(dit(icgts*3),itc(icgts))
      if(allocated(dit)) then
       deallocate(dit)
       allocate(dit(icgts*4))
      else 
       allocate(dit(icgts*4))
      endif
      if(allocated(itc)) then
       deallocate(itc)
       allocate(itc(icgts))
      else
       allocate(itc(icgts))
      endif
      dit=0
      itc=0
      index=1
      do i2=1,icgts
c     Added max time step to the end of this line
c     Parse the line to determine the number of input
c     numbers - this makes the code backward compatible
c     if max time step numbers are not input, the value
c     is left at the valu used in ctrl macro
         read (inpt,'(a)') input_msg
         call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
c     Needs to be made backward compatible!
c         read(inpt,*) diti, ditip, ditip2, itci
         if (msg(1).eq.1) then
            diti=imsg(1)
         else
            diti=xmsg(1)
         endif
         if (msg(2).eq.1) then
            ditip=imsg(2)
         else
            ditip=xmsg(2)
         endif
         if (msg(3).eq.1) then
            ditip2=imsg(3)
         else
            ditip2=xmsg(3)
         endif
         dit(index)=diti
         dit(index+icgts)=ditip
         dit(index+icgts*2)=ditip2
c gaz 062920 added coding to read a real or integer  itc()       
         if (msg(4).eq.1) then
            itc(index)=imsg(4)
         else
            itc(index)=xmsg(4)
         endif         
c     Store max time step value
         if(nwds.gt.4) then
            if (msg(5).eq.1) then
               dit(index+icgts*3)=imsg(5)
            else
               dit(index+icgts*3)=xmsg(5)
            endif
         end if
         index=index+1
      enddo
 99     read(inpt,'(a80)') wdd1

      end
      
