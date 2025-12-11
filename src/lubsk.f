      subroutine lubsk(a,n,np,indx,b)
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
!D1 To .
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Mar-04, Programmer: S Kelkar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/lubsk.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.6 Streamline particle-tracking module
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      implicit none

      integer n,np,indx(n)
      integer i,ii,j,ll

      real*8 a(np,np),b(np),sum

      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if(ii.ne.0) then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            enddo
         elseif(sum.ne.0.) then
            ii=i
         endif
         b(i)=sum
      enddo
      do i=n,1,-1
         sum=b(i)
         do j=i+1,n
            sum=sum-a(i,j)*b(j)
         enddo
         b(i)=sum/a(i,i)
      enddo
      return
      end
