      subroutine read_rtd(iocntl,fractional_area,
     2     pivottime,fact0,factfinal)
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read rtd curve for the mixing model.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 17-SEP-04,  Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/read_rtd.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************


c     4. read f(t)

      use comrtd, only : rtdcount, time_rtd, rtd
      implicit none
      logical done
      integer i, iocntl,isize
      real*8 tdummy
      real*8 rtddummy
      real*8 fractional_area,pivottime,fact0,factfinal, factor
c     first read through entire file and count the number of points

      done = .false.
      rtdcount = 0
      do while(.not.done)
         read(iocntl,*,end=1000) tdummy, rtddummy
         rtdcount = rtdcount + 1
      end do
 1000 continue
      rewind iocntl
c     Add 25 more spaces for extrapolated tail of curve

      rtdcount = rtdcount + 25

c     Allocate space in rtd arrays
      if (.not. allocated(time_rtd)) then
         allocate(time_rtd(rtdcount))
         allocate(rtd(rtdcount))
      else
         isize=size(time_rtd,1)
         if(isize.ne.rtdcount) then
            deallocate(time_rtd,rtd)
            allocate(time_rtd(rtdcount),rtd(rtdcount))
         end if       
      end if

      time_rtd = 0.
      rtd = 0.

c     Read in rtd curve

      do i = 1, rtdcount-25
         read(iocntl,*) time_rtd(i), rtd(i)
      end do

c     Apply factors to the curve
      do i = 1, rtdcount-25
         if(time_rtd(i).lt.pivottime) then
            factor=fact0+(1.-fact0)*time_rtd(i)/pivottime
         else
            factor=1.+(factfinal-1.)*(time_rtd(i)-pivottime)/
     2           (time_rtd(rtdcount)-pivottime)
            end if
         rtd(i)=rtd(i)*factor
      end do

      end subroutine read_rtd
