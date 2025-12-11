      logical function null_new(wd)
!***********************************************************************
!  Copyright, 2006,  The  Regents  of the  University of California.
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
!D1 Check for null lines.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 06-May-02, Programmer: Z. Dash        
!D2
!D2 $Log:   /pvcs.config/fehm90/src/null_new.f_a  $
!D2 
!***********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 N/A
!D3
!***********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 None
!D4
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS BEGIN null_new
!PS 
!PS   determine the length of the input string
!PS 
!PS   FOR each character in the string
!PS   
!PS       IF the character is not a " " or tab
!PS          set the value to .FALSE. and return
!PS       ENDIF   
!PS       
!PS   ENDFOR
!PS   
!PS   set the value to .TRUE. and return
!PS       
!PS END null_new
!PS
!***********************************************************************

      implicit none

      integer i,length
      character*(*) wd

      length = len(wd)
      do i = 1, length
         if (wd(i:i) .ne. ' ' .and. 
     *       wd(i:i) .ne. char(9)) then
            null_new = .FALSE.
            return
         endif
      end do

      null_new = .TRUE.

      return
      end
