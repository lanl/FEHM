      subroutine nrmlz4(nel,a,na,r,nr,ncon,rnmfin,rw)
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
CD1  This routine normalizes the matrix equations.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/nrmlz4.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:34   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:18   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Wed Jan 10 15:05:52 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   11/15/95 15:16:00   gaz
CD2 added nr1 where appropriate(minor correction)
CD2 
CD2    Rev 1.1   03/18/94 16:07:20   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:02   pvcs
CD2 original version in process of being certified
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.4.2     Solve Nonlinear Equation at Each Time Step
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
C**********************************************************************

      implicit none

      integer nel,ncon(*),na(16),nr(4)
      integer i,ifin,ii,iiind,ij,ijind,ista,nelp1,im1neq,na11,na12,
     &     na13,na14,na21,na22,na23,na24,na31,na32,na33,na34,na41,
     &     na42,na43,na44,nr1,nr2,nr3,nr4,inddd1,inddd2,inddd3,inddd4
      real*8 epn,a(*),r(nel*4),rw(nel*4)
      real*8 rnmfin,rnmsta,dd11,dd12,dd13,dd14,dd21,dd22,dd23,dd24
     &  ,dd31,dd32,dd33,dd34,dd41,dd42,dd43,dd44,swap,vv1,vv2,vv3,vv4
     &  ,pivot1,pivot2,pivot3,pivot4,pivmax,bksb1,bksb2,bksb3,bksb4

      nelp1=nel+1
      nr1=nr(1)
      nr2=nr(2)
      nr3=nr(3)
      nr4=nr(4)
      na11=na( 1)
      na12=na( 2)
      na13=na( 3)
      na14=na( 4)
      na21=na( 5)
      na22=na( 6)
      na23=na( 7)
      na24=na( 8)
      na31=na( 9)
      na32=na(10)
      na33=na(11)
      na34=na(12)
      na41=na(13)
      na42=na(14)
      na43=na(15)
      na44=na(16)
      
c     normalize a matrix and r.  this involves multplying ax=r by the
c     inverse of the matrix which has the diagonal blocks of a on its
c     diagonal and zeros elsewhere.
      im1neq=0
      do i=1,nel
         ista=ncon(i)+1
         ifin=ncon(i+1)
c     find the diagonal (can't use npvt because it is for b)
         ii=ista
 30      if ((ncon(ii).eq.i).or.(ii.gt.ifin)) goto 40
         ii=ii+1
         goto 30
 40      continue
c     copy the block on the diagonal into dd
         iiind=ii-nelp1
         dd11=a(iiind+na11)
         dd12=a(iiind+na12)
         dd13=a(iiind+na13)
         dd14=a(iiind+na14)
         dd21=a(iiind+na21)
         dd22=a(iiind+na22)
         dd23=a(iiind+na23)
         dd24=a(iiind+na24)
         dd31=a(iiind+na31)
         dd32=a(iiind+na32)
         dd33=a(iiind+na33)
         dd34=a(iiind+na34)
         dd41=a(iiind+na41)
         dd42=a(iiind+na42)
         dd43=a(iiind+na43)
         dd44=a(iiind+na44)
c     calculate the inverse of dd (the block on the diagonal)
         vv1=1.0/max(abs(dd11),abs(dd12),abs(dd13),abs(dd14))
         vv2=1.0/max(abs(dd21),abs(dd22),abs(dd23),abs(dd24))
         vv3=1.0/max(abs(dd31),abs(dd32),abs(dd33),abs(dd34))
         vv4=1.0/max(abs(dd41),abs(dd42),abs(dd43),abs(dd44))
         pivot1=vv1*abs(dd11)
         pivot2=vv2*abs(dd21)
         pivot3=vv3*abs(dd31)
         pivot4=vv4*abs(dd41)
         pivmax=max(pivot1,pivot2,pivot3,pivot4)
         inddd1=1
         if (pivot2.eq.pivmax) then
            swap=dd11
            dd11=dd21
            dd21=swap
            swap=dd12
            dd12=dd22
            dd22=swap
            swap=dd13
            dd13=dd23
            dd23=swap
            swap=dd14
            dd14=dd24
            dd24=swap
            inddd1=2
            vv2=vv1
         else
            if (pivot3.eq.pivmax) then
               swap=dd11
               dd11=dd31
               dd31=swap
               swap=dd12
               dd12=dd32
               dd32=swap
               swap=dd13
               dd13=dd33
               dd33=swap
               swap=dd14
               dd14=dd34
               dd34=swap
               inddd1=3
               vv3=vv1
            else
               if (pivot4.eq.pivmax) then
                  swap=dd11
                  dd11=dd41
                  dd41=swap
                  swap=dd12
                  dd12=dd42
                  dd42=swap
                  swap=dd13
                  dd13=dd43
                  dd43=swap
                  swap=dd14
                  dd14=dd44
                  dd44=swap
                  inddd1=4
                  vv4=vv1
               endif
            endif
         endif
         dd21=dd21/dd11
         dd31=dd31/dd11
         dd41=dd41/dd11
         dd22=dd22-dd21*dd12
         dd32=dd32-dd31*dd12
         dd42=dd42-dd41*dd12
         pivot2=vv2*abs(dd22)
         pivot3=vv3*abs(dd32)
         pivot4=vv4*abs(dd42)
         pivmax=max(pivot2,pivot3,pivot4)
         inddd2=2
         if (pivot3.eq.pivmax) then
            swap=dd21
            dd21=dd31
            dd31=swap
            swap=dd22
            dd22=dd32
            dd32=swap
            swap=dd23
            dd23=dd33
            dd33=swap
            swap=dd24
            dd24=dd34
            dd34=swap
            inddd2=3
            vv3=vv2
         else
            if (pivot4.eq.pivmax) then
               swap=dd21
               dd21=dd41
               dd41=swap
               swap=dd22
               dd22=dd42
               dd42=swap
               swap=dd23
               dd23=dd43
               dd43=swap
               swap=dd24
               dd24=dd44
               dd44=swap
               inddd2=4
               vv4=vv2
            endif
         endif
         dd32=dd32/dd22
         dd42=dd42/dd22
         dd23=dd23-dd21*dd13
         dd33=dd33-dd31*dd13-dd32*dd23
         dd43=dd43-dd41*dd13-dd42*dd23
         pivot3=vv3*abs(dd33)
         pivot4=vv4*abs(dd43)
         pivmax=max(pivot3,pivot4)
         inddd3=3
         if (pivot4.eq.pivmax) then
            swap=dd31
            dd31=dd41
            dd41=swap
            swap=dd32
            dd32=dd42
            dd42=swap
            swap=dd33
            dd33=dd43
            dd43=swap
            swap=dd34
            dd34=dd44
            dd44=swap
            inddd3=4
            vv4=vv3
         endif
         dd43=dd43/dd33
         dd24=dd24-dd21*dd14
         dd34=dd34-dd31*dd14-dd32*dd24
         dd44=dd44-dd41*dd14-dd42*dd24-dd43*dd34
c     calculate dd(-1)*r and put into the solution vector, rw
         bksb1=r(i+nr1)
         bksb2=r(i+nr2)
         bksb3=r(i+nr3)
         bksb4=r(i+nr4)
         if (inddd1.eq.2) then
            swap=bksb1
            bksb1=bksb2
            bksb2=swap
         else
            if (inddd1.eq.3) then
               swap=bksb1
               bksb1=bksb3
               bksb3=swap
            else
               if (inddd1.eq.4) then
                  swap=bksb1
                  bksb1=bksb4
                  bksb4=swap
               endif
            endif
         endif
         if (inddd2.eq.3) then
            swap=bksb2
            bksb2=bksb3
            bksb3=swap
         else
            if (inddd2.eq.4) then
               swap=bksb2
               bksb2=bksb4
               bksb4=swap
            endif
         endif
         bksb2=bksb2-dd21*bksb1
         if (inddd3.eq.4) then
            swap=bksb3
            bksb3=bksb4
            bksb4=swap
         endif
         bksb3=bksb3-dd31*bksb1-dd32*bksb2
         bksb4=bksb4-dd41*bksb1-dd42*bksb2-dd43*bksb3
         bksb4=bksb4/dd44
         bksb3=bksb3-dd34*bksb4
         bksb3=bksb3/dd33
         bksb2=bksb2-dd23*bksb3-dd24*bksb4
         bksb2=bksb2/dd22
         bksb1=bksb1-dd12*bksb2-dd13*bksb3-dd14*bksb4
         bksb1=bksb1/dd11
         rw(im1neq+1)=bksb1
         rw(im1neq+2)=bksb2
         rw(im1neq+3)=bksb3
         rw(im1neq+4)=bksb4
c     calculate dd(-1) times the other blocks in the row
         ista=ista-nelp1
         ifin=ifin-nelp1
         do ijind=ista,ifin
            if (ijind.ne.iiind) then
               bksb1=a(ijind+na11)
               bksb2=a(ijind+na21)
               bksb3=a(ijind+na31)
               bksb4=a(ijind+na41)
               if (inddd1.eq.2) then
                  swap=bksb1
                  bksb1=bksb2
                  bksb2=swap
               else
                  if (inddd1.eq.3) then
                     swap=bksb1
                     bksb1=bksb3
                     bksb3=swap
                  else
                     if (inddd1.eq.4) then
                        swap=bksb1
                        bksb1=bksb4
                        bksb4=swap
                     endif
                  endif
               endif
               if (inddd2.eq.3) then
                  swap=bksb2
                  bksb2=bksb3
                  bksb3=swap
               else
                  if (inddd2.eq.4) then
                     swap=bksb2
                     bksb2=bksb4
                     bksb4=swap
                  endif
               endif
               bksb2=bksb2
     &              -dd21*bksb1
               if (inddd3.eq.4) then
                  swap=bksb3
                  bksb3=bksb4
                  bksb4=swap
               endif
               bksb3=bksb3-dd31*bksb1-dd32*bksb2
               bksb4=bksb4-dd41*bksb1-dd42*bksb2-dd43*bksb3
               bksb4=bksb4/dd44
               bksb3=bksb3-dd34*bksb4
               bksb3=bksb3/dd33
               bksb2=bksb2-dd23*bksb3-dd24*bksb4
               bksb2=bksb2/dd22
               bksb1=bksb1-dd12*bksb2-dd13*bksb3-dd14*bksb4
               bksb1=bksb1/dd11
               a(ijind+na11)=bksb1
               a(ijind+na21)=bksb2
               a(ijind+na31)=bksb3
               a(ijind+na41)=bksb4
               bksb1=a(ijind+na12)
               bksb2=a(ijind+na22)
               bksb3=a(ijind+na32)
               bksb4=a(ijind+na42)
               if (inddd1.eq.2) then
                  swap=bksb1
                  bksb1=bksb2
                  bksb2=swap
               else
                  if (inddd1.eq.3) then
                     swap=bksb1
                     bksb1=bksb3
                     bksb3=swap
                  else
                     if (inddd1.eq.4) then
                        swap=bksb1
                        bksb1=bksb4
                        bksb4=swap
                     endif
                  endif
               endif
               if (inddd2.eq.3) then
                  swap=bksb2
                  bksb2=bksb3
                  bksb3=swap
               else
                  if (inddd2.eq.4) then
                     swap=bksb2
                     bksb2=bksb4
                     bksb4=swap
                  endif
               endif
               bksb2=bksb2-dd21*bksb1
               if (inddd3.eq.4) then
                  swap=bksb3
                  bksb3=bksb4
                  bksb4=swap
               endif
               bksb3=bksb3-dd31*bksb1-dd32*bksb2
               bksb4=bksb4-dd41*bksb1-dd42*bksb2-dd43*bksb3
               bksb4=bksb4/dd44
               bksb3=bksb3-dd34*bksb4
               bksb3=bksb3/dd33
               bksb2=bksb2-dd23*bksb3-dd24*bksb4
               bksb2=bksb2/dd22
               bksb1=bksb1-dd12*bksb2-dd13*bksb3-dd14*bksb4
               bksb1=bksb1/dd11
               a(ijind+na12)=bksb1
               a(ijind+na22)=bksb2
               a(ijind+na32)=bksb3
               a(ijind+na42)=bksb4
               bksb1=a(ijind+na13)
               bksb2=a(ijind+na23)
               bksb3=a(ijind+na33)
               bksb4=a(ijind+na43)
               if (inddd1.eq.2) then
                  swap=bksb1
                  bksb1=bksb2
                  bksb2=swap
               else
                  if (inddd1.eq.3) then
                     swap=bksb1
                     bksb1=bksb3
                     bksb3=swap
                  else
                     if (inddd1.eq.4) then
                        swap=bksb1
                        bksb1=bksb4
                        bksb4=swap
                     endif
                  endif
               endif
               if (inddd2.eq.3) then
                  swap=bksb2
                  bksb2=bksb3
                  bksb3=swap
               else
                  if (inddd2.eq.4) then
                     swap=bksb2
                     bksb2=bksb4
                     bksb4=swap
                  endif
               endif
               bksb2=bksb2-dd21*bksb1
               if (inddd3.eq.4) then
                  swap=bksb3
                  bksb3=bksb4
                  bksb4=swap
               endif
               bksb3=bksb3-dd31*bksb1-dd32*bksb2
               bksb4=bksb4-dd41*bksb1-dd42*bksb2-dd43*bksb3
               bksb4=bksb4/dd44
               bksb3=bksb3-dd34*bksb4
               bksb3=bksb3/dd33
               bksb2=bksb2-dd23*bksb3-dd24*bksb4
               bksb2=bksb2/dd22
               bksb1=bksb1-dd12*bksb2-dd13*bksb3-dd14*bksb4
               bksb1=bksb1/dd11
               a(ijind+na13)=bksb1
               a(ijind+na23)=bksb2
               a(ijind+na33)=bksb3
               a(ijind+na43)=bksb4
               bksb1=a(ijind+na14)
               bksb2=a(ijind+na24)
               bksb3=a(ijind+na34)
               bksb4=a(ijind+na44)
               if (inddd1.eq.2) then
                  swap=bksb1
                  bksb1=bksb2
                  bksb2=swap
               else
                  if (inddd1.eq.3) then
                     swap=bksb1
                     bksb1=bksb3
                     bksb3=swap
                  else
                     if (inddd1.eq.4) then
                        swap=bksb1
                        bksb1=bksb4
                        bksb4=swap
                     endif
                  endif
               endif
               if (inddd2.eq.3) then
                  swap=bksb2
                  bksb2=bksb3
                  bksb3=swap
               else
                  if (inddd2.eq.4) then
                     swap=bksb2
                     bksb2=bksb4
                     bksb4=swap
                  endif
               endif
               bksb2=bksb2-dd21*bksb1
               if (inddd3.eq.4) then
                  swap=bksb3
                  bksb3=bksb4
                  bksb4=swap
               endif
               bksb3=bksb3-dd31*bksb1-dd32*bksb2
               bksb4=bksb4-dd41*bksb1-dd42*bksb2-dd43*bksb3
               bksb4=bksb4/dd44
               bksb3=bksb3-dd34*bksb4
               bksb3=bksb3/dd33
               bksb2=bksb2-dd23*bksb3-dd24*bksb4
               bksb2=bksb2/dd22
               bksb1=bksb1-dd12*bksb2-dd13*bksb3-dd14*bksb4
               bksb1=bksb1/dd11
               a(ijind+na14)=bksb1
               a(ijind+na24)=bksb2
               a(ijind+na34)=bksb3
               a(ijind+na44)=bksb4
            endif
         enddo
c     put identity blocks onto the diagonal of a
         a(iiind+na11)=1.
         a(iiind+na12)=0.
         a(iiind+na13)=0.
         a(iiind+na14)=0.
         a(iiind+na21)=0.
         a(iiind+na22)=1.
         a(iiind+na23)=0.
         a(iiind+na24)=0.
         a(iiind+na31)=0.
         a(iiind+na32)=0.
         a(iiind+na33)=1.
         a(iiind+na34)=0.
         a(iiind+na41)=0.
         a(iiind+na42)=0.
         a(iiind+na43)=0.
         a(iiind+na44)=1.
         im1neq=im1neq+4
      enddo
c     calculate the norm of the normalized right-hand side (rnmfin)
      rnmfin=0.
      im1neq=0
      do i=1,nel
         r(i+nr1)=rw(im1neq+1)
         r(i+nr2)=rw(im1neq+2)
         r(i+nr3)=rw(im1neq+3)
         r(i+nr4)=rw(im1neq+4)
         rnmfin=rnmfin+r(i+nr1)*r(i+nr1)
     &        +r(i+nr2)*r(i+nr2)
     &        +r(i+nr3)*r(i+nr3)
     &        +r(i+nr4)*r(i+nr4)
         im1neq=im1neq+4
      enddo

      return
      end
