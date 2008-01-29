      subroutine lstsq(ndata,y,t)
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
!D1 This program calculates the least squares fit to a quartic 
!D1 polynomial using QR factorization and solve, and a cholesky
!D1 solution to the normal equations of the system.  It also 
!D1 compares the solution vectors.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/lstsq.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3  2.4.6 Multiple, interacting solutes
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!D4 y(ndata,y2(ndata)) = dependent variable array
!D4 t(ndata) = independent variable array
!D4 A(ndata,ncoeff),a2(ndata,ncoeff) = resulting matrix
!D4 AT(ndata,ncoeff) = transpose of resulting matrix
!D4 ch(ncoeff,ncoeff) = AT*A (by matrix multiplication)
!D4 p(ncoeff) = vector required by cholesky routine
!D4 x(ncoeff) = solution vector for cholesky
!D4 b(ncoeff) = altered LHS in cholesky/normal equations
!D4 ans1(ndata) = predicted points according to cholesky
!D4 n,np,m = matrix size parameters
!D4 
!***********************************************************************

      implicit none

      integer ndata, ncoeff
      real*8 y(*)
      real*8 t(*)
      real*8, allocatable :: A(:,:)
      real*8, allocatable :: AT(:,:)
      real*8, allocatable :: y2(:)
      real*8, allocatable :: a2(:,:)
      real*8, allocatable :: ch(:,:)
      real*8, allocatable :: b(:)
      real*8, allocatable :: pred(:)
      real*8, allocatable :: x(:)
      real*8, allocatable :: p(:)
      integer n,np,m,i,j,k
      ncoeff = 3
      allocate (A(ndata,ncoeff),AT(ncoeff,ndata))
      allocate(y2(ndata),a2(ndata,ncoeff),ch(ncoeff,ncoeff),b(ncoeff))
      allocate (pred(ndata),x(ndata))
      allocate(p(ncoeff))

      n=ncoeff
      np=n
      m=ndata
      do i=1,ncoeff
         do j=1,ncoeff
            ch(i,j)=0
         enddo
         b(i)=0
      enddo

      
      do i=1,ndata

         y2(i)=y(i)
         A(i,1)=1
         A(i,2)=t(i)
         A(i,3)=t(i)**2

         do j=1,ncoeff
            a2(i,j)=A(i,j)
            AT(j,i)=A(i,j)
         enddo

      enddo

*     CHOLESKY FACTORIZATION/NORMAL EQUATION METHOD
      
      do j=1,ncoeff
         do k=1,ncoeff
            do i=1,ndata
               ch(k,j)=AT(j,i)*a2(i,k)+ch(k,j)
            enddo
         enddo
      enddo

      do j=1,ncoeff
         do i=1,ndata
            b(j)=AT(j,i)*y2(i)+b(j)
         enddo
      enddo

      np=n
       
      call choldc(ch,n,np,p)
      call cholsl(ch,n,np,p,b,x)
      do i=1,ndata
         pred(i)=x(1)+x(2)*t(i)+x(3)*t(i)**2
      enddo
      do i = 1,ncoeff
         y(i)=x(i)
      enddo
      deallocate(A,AT)
      deallocate(y2,a2,ch,b)
      deallocate (pred,x)
      deallocate(p)

      end

*****************************************************************
*
* 	ALTERED SUBROUTINES
*
*****************************************************************
  
      SUBROUTINE cholsl(a,n,np,p,b,x)
      implicit none
      INTEGER n,np
      REAL*8 a(np,np),b(n),p(n),x(n)
      INTEGER i,k
      REAL*8 sum
      do 12 i=1,n
        sum=b(i)
        do 11 k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
11      continue
        x(i)=sum/p(i)
12    continue
      do 14 i=n,1,-1
        sum=x(i)
        do 13 k=i+1,n
          sum=sum-a(k,i)*x(k)
13      continue
        x(i)=sum/p(i)
14    continue
      return
      END

      SUBROUTINE choldc(a,n,np,p)
      implicit none
      INTEGER n,np
      REAL*8 a(np,np),p(n)
      INTEGER i,j,k
      REAL*8 sum
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.)pause 'choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END

      SUBROUTINE rsolv(a,n,np,d,b)
      implicit none
      INTEGER n,np
      REAL*8 a(np,np),b(n),d(n)
      INTEGER i,j
      REAL*8 sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      END

      SUBROUTINE qrsolv(a,n,np,c,d,b,m)
      implicit none
      INTEGER n,np,m
      REAL*8 a(m,np),b(n),c(n),d(n)
CU    USES rsolv
      INTEGER i,j
      REAL*8 sum,tau
      do 13 j=1,n
        sum=0.
        do 11 i=j,n
          sum=sum+a(i,j)*b(i)
11      continue
        tau=sum/c(j)
        do 12 i=j,n
          b(i)=b(i)-tau*a(i,j)
12      continue
13    continue
      call rsolv(a,n,np,d,b)
      return
      END
  
      SUBROUTINE qrdcmp(a,n,np,c,d,m)
      implicit none
      INTEGER n,np,m
      REAL*8 a(m,np),c(n),d(n)
      INTEGER i,j,k
      REAL*8 scale,sigma,sum,tau
      logical sing
      sing=.false.
      scale=0.
      do 17 k=1,n
        do 11 i=k,m
          scale=max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.0.)then
          c(k)=0.
          d(k)=0.
        else
          do 12 i=k,m
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.
          do 13 i=k,m
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k,n
            sum=0.
            do 14 i=k,m
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,m
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      return
      END


