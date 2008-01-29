      subroutine lubksb(a,ia,iabase,n,indx,indbas,b,ib,ibbase,m)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1  
CD1  To solve the set of n linear equations a*x=b. 
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/lubksb.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 16 10:04:26 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.2   Wed Jan 10 14:12:34 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 16:07:08   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:26   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable. This is a general utility routine used in the code as 
CD3  needed for solving a matrix.
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4  Solves the set of n linear equations a*x=b. Here a is input as
CD4  its lu decomposition (determined by the routine ludcmp).  indx is
CD4  input as the permutation vector returned by ludcmp.  b is input
CD4  as the right-hand side vector and returns with the solution
CD4  vector.  a, n, np and indx are not modified by this routine and
CD4  can be left in place for successive calls with different right-
CD4  hand sides.  this routine takes into account the possibility that
CD4  b will begin with many zero elements, so it is efficient for use
CD4  in matrix inversion.
CD4
C***********************************************************************
c     8 august 1988.  modified to reflect the way that the coefficient
c     matrix is stored for george's routines.
c
c     9 august 1988.  modified to reflect the way that the rhs is stored
c     for george's routines and to allow for multiple rhss.  the rhs
c     is not stored in contiguous positions of the array b so that a
c     reference position (ibbase) and an indexing array (ib) are
c     needed.  m is the number of rhss.
C***********************************************************************

      implicit none
      integer i, ii, j, k, ll
      integer ia(*), iabase, ib(*), ibbase, indx(*), indbas, m, n      
      real*8  a(*), b(*)
      real*8  sum
      
      do k=1,m
         ii=0
         do i=1,n
            ll=indx(indbas+i-1)
            sum=b(ibbase+ib((ll-1)*m+k))
            b(ibbase+ib((ll-1)*m+k))=b(ibbase+ib((i-1)*m+k))
            if (ii.ne.0) then
               do j=ii,i-1
                  sum=sum-a(iabase+ia((i-1)*n+j))*
     &                 b(ibbase+ib((j-1)*m+k))
               enddo
            else if (sum.ne.0.) then
               ii=i
            endif
            b(ibbase+ib((i-1)*m+k))=sum
         enddo
         do i=n,1,-1
            sum=b(ibbase+ib((i-1)*m+k))
            if (i.lt.n) then
               do j=i+1,n
                  sum=sum-a(iabase+ia((i-1)*n+j))*
     &                 b(ibbase+ib((j-1)*m+k))
               enddo
            endif
            b(ibbase+ib((i-1)*m+k))=sum/a(iabase+ia((i-1)*n+i))
         enddo
      enddo
      
      return
      end
