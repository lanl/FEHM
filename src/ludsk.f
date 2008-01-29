      subroutine ludsk(a,n,np,indx,d,ludflag)
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
!D2 $Log:   /pvcs.config/fehm90/src/ludsk.f_a  $
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

      integer n,np,indx(n),nmax,ludflag
      integer i,imax,j,k

      real*8 d,a(np,np),tiny

      parameter (nmax=500,tiny=1.e-20)

      real*8 aamax,dum,sum,vv(nmax)

      ludflag = 0
      d=1.

      do i=1,n
         aamax=0.
         do j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
         if(aamax.eq.0.) then
            ludflag = i
            return
         endif
         vv(i)=1./aamax
      enddo

      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
         enddo
         aamax=0.
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if(dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if(j.ne.imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.) a(j,j)=tiny
         if(j.ne.n) then
            dum=1./a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo
         endif
      enddo
      return
      end
