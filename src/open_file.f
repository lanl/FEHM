      integer function open_file(filename, filestat)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To find an unused unit number and open the file specified.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: 25-JAN-00, Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/open_file.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:34   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comai, only : ierr, iptty
      implicit none

      character(*) filename
      character(*) filestat
      logical ex, used
      integer irfile

c     Find an unused unit number, open the file if specified
      used = .true.
c zvd 04/06/2007 Make sure we don't conflict with reserved file numbers
c units, 1,2,3,11-35,41,42,79,80,83,92,121
      irfile = 50
      do while(used)
         inquire(unit = irfile, opened = used)
         if(.not.used) then
            
c     Open file to read
            used = .false.
            open_file = irfile
	    open(irfile, file = trim(filename), status=filestat, 
     &           err = 100)
         else
c     Unit number used, try the next one
            used = .true.
            irfile = irfile + 1
         end if
      end do

      return

 100  write (ierr, *) 'ERROR opening ', trim(filename)
      write (ierr, *) 'STOPPING execution'
      if (iptty .ne. 0) then
         write (iptty, *) 'ERROR opening ', trim(filename)
         write (iptty, *) 'STOPPING execution'
      end if
      inquire (file = filename, exist = ex)
      if (.not. ex)  then
         write (ierr, *) 'FILE does not exist'
         if (iptty .ne. 0) write (iptty, *) 'FILE does not exist'
      end if
      stop

      end

