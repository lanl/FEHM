    
       subroutine indiff
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
CD1 Read diffusivity data for two phase isothermal problem (liquid phase water)
CD1 initial coding gaz 112023 
CD1
******************************************************************

      use comdi
      use comdti
      use comfi, only : idiff_iso
      use comai
      use comki
      implicit none
      
      integer i
      
      real*8, allocatable :: dumx(:)
      real*8, allocatable :: dumy(:)
      real*8, allocatable :: dumz(:)
      macro = 'diff'

      
c**** read diffusivity data ****
      
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = zero_t
      default(2) = zero_t
      default(3) = zero_t
      allocate(dumx(n0),dumy(n0),dumz(n0))
      igroup = 1
      if(ico2.lt.0.and.l.eq.0) then
       if(allocated(thx)) deallocate (thx)
       if(allocated(thy)) deallocate (thy)
       if(allocated(thz)) deallocate (thz)
       allocate(thx(n0),thy(n0),thz(n0))
      endif

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     &     default, macroread(10), macro, igroup, ireturn,
     &     r8_1=dumx(1:n0),r8_2=dumy(1:n0),r8_3=dumz(1:n0)) 

      if(ico2.lt.0) then
         do i = 1, n0
            thx(i) = max (zero_t,dumx(i))
            thy(i) = max (zero_t,dumy(i))
            thz(i) = max (zero_t,dumz(i))
         end do
      endif

      macroread(10) = .TRUE.

      deallocate(dumx,dumy,dumz)
      end
