      subroutine normal_uc(neq,idegree,a,nmat,r,nrhs,nelmdg,nelm,iptty,
     & fdum2)
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
CD1  This subroutine normalizes the equations in an uncoupled manner.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/normal_uc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:34   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Wed Jan 10 14:55:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:59:42   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:00   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.4.2     Solve Nonlinear Equation at Each Time Step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
C***********************************************************************

      implicit none

      integer idegree, iptty, neq 
      integer nelm(*), nelmdg(*), nrhs(*), nmat(*)
      integer i, i1, i2, idegp1, idiag, ii, j, jj1, jj2
      integer k, kk, mm, neqp1, nmatd, nmatmm
      real*8  a(*), r(*), aii, fdum2, tolaii
      parameter(tolaii=0.000000000000000000001)
      
      idegp1=idegree+1
      neqp1=neq+1
      fdum2=0.0
      do j = 1,idegree
         idiag=(j-1)*idegp1 + 1
         nmatd=nmat(idiag)
         do i=1,neq
            ii=nelmdg(i)-neqp1
            aii=a(nmatd+ii)
c     check for zero on diagonal
            if(abs(aii).le.tolaii) then
               write(*,*) '$$$$$$$$$$$$$$$$$    $$$$$$$$$$$$$$$$$$$'
               write(*,*)'zero diagonal found in normal_uc'
               write(*,*) 'dof= ',j,' node = ',i,' aii= ',aii
               write(iptty,*) '$$$$$$$$$$$$$$$$$    $$$$$$$$$$$$$$$$$$$'
               write(iptty,*)'zero diagonal found in normal_uc'
               write(iptty,*) 'dof= ',j,' node = ',i,' aii= ',aii
               stop
            endif
            r(i+nrhs(j))=r(i+nrhs(j))/aii
            fdum2=fdum2+r(i+nrhs(j))*r(i+nrhs(j))
            jj1=(j-1)*idegree+1
            jj2=j*idegree
            i1=nelm(i)+1
            i2=nelm(i+1)
            do mm=jj1,jj2
               nmatmm=nmat(mm)
               do k=i1,i2
                  kk=k-neqp1
                  a(kk+nmatmm)=a(kk+nmatmm)/aii
               enddo
            enddo
         enddo
      enddo

      return
      end
