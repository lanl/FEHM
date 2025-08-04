      subroutine innondarcy
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
CD1 Read non darcy  data.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 24-JAN-2025    G Zyvoloski        NA      Initial implementation.
CD2
CD2 
C***********************************************************************

      use comdi
      use comci
      use comdti
      use comai
      use comki
      use com_nondarcy
      implicit none
      
      macro = 'ndar'
c      
c**** read non darcy parameter beta ****
c
      if(.not.allocated(nd_beta))  then
        allocate(nd_beta(n0))
        allocate(rlf_nd(n0),rvf_nd(n0),drlef_nd(n0),drvef_nd(n0))
      endif
      narrays = 1
      itype(1) = 8
      default(1) = 1.d-15
      igroup = 1
      
      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(17), macro, igroup, ireturn,
     2     r8_1=nd_beta(1:n0)) 
      
      macroread(17) = .TRUE.
      
      end
 