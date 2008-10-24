      subroutine namefile1(lund,iformat,root,tail,nchar_root,ierr)
!***********************************************************************
!  Copyright, 2006,  The  Regents of the University of California.
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
!D1 Generate name for AVS format output header file.
!D1
C***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 23-Feb-06, Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/namefile1.f_a  $
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
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

      implicit none

      integer i, lund, iformat, ierr, open_file
      integer nchar_root, nchar_tail, c_open
      character*125 fname
      character*94 root
      character*(*) tail
      character*5 ch5

c Blank out the file name 
      do i = 1, 125
         fname(i:i)=' '
      end do

      nchar_tail=len_trim(tail)

      fname(1:nchar_root) = root(1:nchar_root)
      fname(nchar_root+1:nchar_root+nchar_tail) = tail(1:nchar_tail)

c     open to output file
      if(iformat .eq. 1)then
         continue
      elseif(iformat .eq. 2)then
         lund = open_file(fname,'unknown')
!         open(unit=lund,file=fname,form='formatted')
      endif

      return
      end
