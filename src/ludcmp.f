      subroutine ludcmp(a,ia,iabase,n,indx,indbas,d,vv,iptty,ierr)
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
CD1  Given an n*n matrix a, with physical dimension np, this routine
CD1  replaces it by the lu decomposition of a rowwise decomposition
CD1  of itself.  a and n are input.  a is output with l below the
CD1  diagonal and u on and above the diagonal.  indx is an output
CD1  vector which records the row permutation effected by the
CD1  partial pivotting.  d is output as +/- 1 depending on whether
CD1  the number of row interchanges was even or odd, respectively.
CD1  this routine is used in conjunction with lubksb to solve
CD1  linear equations or invert a matrix.  vv is used for working
CD1  space.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/ludcmp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Feb 16 10:04:28 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.3   Wed Jan 10 14:28:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   06/20/94 11:13:22   zvd
CD2 Added ierr unit number for error output.
CD2 
CD2    Rev 1.1   03/18/94 16:07:10   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:30   pvcs
CD2 original version in process of being certified
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable. This is a general utility routine used in the code as 
CD3  needed for solving a matrix.
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  these two routines for lu decomposition and back-substitution,
CD4  with row pivotting are from
CD4  'numerical recipes the art of scientific computing' by
CD4  w.h.press, b.p.flannery, s.a.teukolsky and w.t.vetterling
CD4  they will be used in george's routines, for inverting the
CD4  small 'block' matricies.
CD4
CD4  Given an n*n matrix a, with physical dimension np, this routine
CD4  replaces it by the lu decomposition of a rowwise decomposition
CD4  of itself.  a and n are input.  a is output with l below the
CD4  diagonal and u on and above the diagonal.  indx is an output
CD4  vector which records the row permutation effected by the
CD4  partial pivotting.  d is output as +/- 1 depending on whether
CD4  the number of row interchanges was even or odd, respectively.
CD4  this routine is used in conjunction with lubksb to solve
CD4  linear equations or invert a matrix.  vv is used for working
CD4  space.
CD4
C***********************************************************************
c       8 august 1988.  modified to reflect the way that the coefficient
c       matrix is stored for george's routines.  it is stored in
c       form in a 1-d array.  iabase is a base position in a for the
c       sub-matrix of interest.  ia contains positions relative to this
c       base.  n is the number of equations in the sub-matrix.  so that
c       coefficient (i,j) is in a(iabase+ia((i-1)*n+j)).
C***********************************************************************

      implicit none
      integer i, imax, j, k
      integer ia(*), iabase, indx(*), indbas, iptty, ierr, n
      real*8  a(*), d, vv(*)
      real*8  aamax, dum, sum

 999  format ('singular matrix in ludcmp')
      d=1.
      do i=1,n
         aamax=0.
         do j=1,n
            if (abs(a(iabase+ia((i-1)*n+j))).gt.aamax)
     &           aamax=abs(a(iabase+ia((i-1)*n+j)))
         enddo
         if (aamax.eq.0.) then
            write (ierr, 999) 
            if (iptty .ne. 0) write(iptty, 999)  
            stop
         endif
         vv(i)=1./aamax
      enddo
      do j=1,n
         if (j.gt.1) then
            do i=1,j-1
               sum=a(iabase+ia((i-1)*n+j))
               if (i.gt.1) then
                  do k=1,i-1
                     sum=sum-a(iabase+ia((i-1)*n+k))
     &                    *a(iabase+ia((k-1)*n+j))
                  enddo
                  a(iabase+ia((i-1)*n+j))=sum
               endif
            enddo
         endif
         aamax=0.
         do i=j,n
            sum=a(iabase+ia((i-1)*n+j))
            if (j.gt.1) then
               do k=1,j-1
                  sum=sum-a(iabase+ia((i-1)*n+k))*
     &                 a(iabase+ia((k-1)*n+j))
               enddo
               a(iabase+ia((i-1)*n+j))=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax) then
            do k=1,n
               dum=a(iabase+ia((imax-1)*n+k))
               a(iabase+ia((imax-1)*n+k))=a(iabase+ia((j-1)*n+k))
               a(iabase+ia((j-1)*n+k))=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(indbas+j-1)=imax
         if (j.ne.n) then
            if (a(iabase+ia((j-1)*n+j)).eq.0.) then
               write (ierr, 999)
               if (iptty .ne. 0) write(iptty, 999) 
               stop
            endif
            dum=1./a(iabase+ia((j-1)*n+j))
            do i=j+1,n
               a(iabase+ia((i-1)*n+j))=a(iabase+ia((i-1)*n+j))*dum
            enddo
         endif
      enddo
      if (a(iabase+ia(n*n)).eq.0.) then
         write (ierr, 999)
         if (iptty .ne. 0) write(iptty, 999)
         stop
      endif

      return
      end
