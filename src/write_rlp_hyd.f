      subroutine write_rlp_hyd(neq,lu,ifdual)
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
!D1 To write hydrate relative permeabilities.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Date 2005, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/write_rlp_hyd.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
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

      use comai, only : altc, days
      use commeth

      implicit none

      integer maxscalar 
      parameter (maxscalar =4) 

      integer i,ifdual,istep,nout,lu,neq
      character*80 title(maxscalar)

      data title( 1) /'Rel Perm Water (no dim)'/,
     1     title( 2) /'Rel Perm Meth (no dim)'/,
     2     title(3) /'Dual Rel Perm Water (no dim)'/
     3     title(4) /'Dual Rel Perm Meth (no dim)'/

      if(ifdual .ne. 0)then
         istep = 2
      else
         istep = 0
      endif
c     allocate rlp memory and calculate relative perms
      call ther_meth_h2o(0,0)

C---  Max number of scalars is 5, output all 5 variables
C---  Output all of these
      nout=2
      if(nout .eq. 2)then

         if(altc(1:4).NE.'avsx') then
            write(lu,90)nout,1,1,1,1,1
            write(lu,'(a56)')title(1+istep)
            write(lu,'(a56)')title(2+istep)
            
            do i = 1,neq
               write(lu,100)i,
     &              max(rl_h2o(i),1d-20), max(rl_meth(i),1d-20)
            enddo
         else
            write(lu,667) days,
     &           (title(i+istep),i=1,6)
            do i = 1,neq
               write(lu,666)i,
     &              max(rl_h2o(i),1d-20), max(rl_meth(i),1d-20)
            enddo
         end if
      end if
c     release rlp memory         
      call ther_meth_h2o(-1,0)

 90   format(i1,2x,5(i5,2x))
 100  format(i10.10,2x,2(e16.9,2x))
 666  format(i10.10,2(' : ',e16.9))
 667  format('nodes at ',e10.4,' days ',2(' : ',a40))

      return
      end

