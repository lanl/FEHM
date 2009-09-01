      subroutine infloa
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
!D1 To read flow data in additive fashion.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 8-Jan-02, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/infloa.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
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

      use combi
      use comdi
      use comdti
      use comai
      use comki
      implicit none

      integer i,icode
      real*8 aiptmp 
      real*8, allocatable ::  aiped(:)
      real*8, allocatable ::  sktmp(:)
      real*8, allocatable ::  esktmp(:)
      real*8, allocatable ::  weltmp(:)

      macro = 'floa'
      
      allocate(aiped(n0),esktmp(n0),sktmp(n0),weltmp(n0))
      
c**** read flow data ****
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = 0.
      default(2) = 0.
      default(3) = 0.
      igroup = 1
      
      call initdata2 (inpt, ischk, n0, narrays, itype, 
     &     default, macroread(12), macro, igroup, ireturn,
     &     r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0)) 
      
      do i = 1, n0
         if (sktmp(i) .ne. default(1) .or. esktmp(i) .ne. default(2)
     &        .or. aiped(i) .ne. default(3)) then
            esk(i) = esktmp(i)
            if(esk(i).gt.0.and.igrav.ne.0.and.ico2.ge.0) then
             esk(i) = esktmp(i)-grav*cord(i,igrav)
            endif
            if (abs(aiped(i)) .lt. zero_t) then
               sk(i) = sk(i)+sktmp(i)
               ka(i) = 1
            else
               if (aiped(i) .lt. 0.) then
                  ka(i) = -2
               else
                  ka(i) = -1
               end if
               aiptmp = abs(aiped(i)) * 1.0e+06
               
               pflow(i) = sktmp(i)*aiptmp + pflow(i)*wellim(i)
               wellim(i) = aiptmp + wellim(i)
               pflow(i) = pflow(i)/wellim(i)
            end if
         end if
      end do
      
      deallocate(aiped,esktmp,sktmp,weltmp)
      
      end

