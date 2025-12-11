      subroutine find_v_3lsq(i,xxx,b3)
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
!***********************************************************************
!D1 
!D1 PURPOSE
!D1 
!D1 Perform velocity calculation for corner node. 
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/find_v_3lsq.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
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
c s kelkar july 9 2004 temp fix- this is a corner node, not many 
c such nodes in the model, so take ggg to be the average of the 
c existing neighbours

      use comsptr
      implicit none

      integer i,j,k,kb

      real*8 xxx,b3(3)

      do j=1,3
         ggg(i,+j)=0.
         ggg(i,-j)=0.
      enddo

      do k=1,3
         kb=irray(i,k)
         if(kb.gt.0) then
            do j=1,3
               ggg(i,+j)=ggg(i,+j)+0.333*ggg(kb,+j)
               ggg(i,-j)=ggg(i,-j)+0.333*ggg(kb,-j)
            enddo
         endif
      enddo
      do k=-3,-1
         kb=irray(i,k)
         if(kb.gt.0) then
            do j=1,3
               ggg(i,+j)=ggg(i,+j)+0.333*ggg(kb,+j)
               ggg(i,-j)=ggg(i,-j)+0.333*ggg(kb,-j)
            enddo
         endif
      enddo

      return

      end
