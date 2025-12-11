      subroutine inflow_wt
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
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
!D1 To read flow data. in additive fashion.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Date 2005, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inflow_wt.f_a  $
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
      use comwt
      implicit none

      integer i,icode
      real*8, allocatable :: aiped(:)
      real*8, allocatable::  sktmp(:)
      real*8, allocatable ::  esktmp(:)

      macro = 'flow'
      allocate(aiped(n0),esktmp(n0),sktmp(n0))
      if(.not.allocated(move_type)) then
         allocate(move_type(n0))
         move_type=0
      endif

      
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
     *     default, macroread(12), macro, igroup, ireturn,
     *     r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0)) 
      
      do i = 1, n0
         if (sktmp(i) .ne. default(1) .or. esktmp(i) .ne. default(2)
     *        .or. aiped(i) .ne. default(3)) then
            esk(i) = esktmp(i)
            if (abs(aiped(i)) .lt. zero_t) then
c     this is a specified flux
               move_wt=1
               move_type(i)= 1
               sk(i) = sktmp(i)
               ka(i) = 1
            else
               if (aiped(i) .lt. 0.) then
                  ka(i) = -22
               else
c     this is a specified head
                  ka(i) = -1
               end if
               wellim(i) = abs(aiped(i)) * 1.0e+06
               pflow(i) = sktmp(i)
            end if
         end if
      end do
      
      deallocate(aiped,esktmp,sktmp)
      
      end
