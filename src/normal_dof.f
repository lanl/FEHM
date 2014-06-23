      subroutine normal_dof(neq,a,r,nelm,nmat,nrhs,nelmdg,indx
     &                      ,idof,diag,diagi,rdum,iflg,fdum2)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To normalize equations for coupled problems.
CD1 
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/normal_dof.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Fri Apr 26 15:55:28 1996   gaz
CD2 small change when singular matrix occurs
CD2 
CD2    Rev 1.1   Fri Feb 16 10:37:04 1996   zvd
CD2 Added prolog
CD2 
C**********************************************************************
CD3
CD3 REQUIREMENTS TRACEABILITY
CD3
CD3 2.5.2 Solve nonlinear equation set at each time step
CD3 
C**********************************************************************
CD4
CD4 SPECIAL COMMENTS
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4 
C**********************************************************************
c equation to be normalized
c Ax=r
c ******************** define arguments **************************
c neq    - number of rows and columns
c a      - matrix to be normalized
c r      - right hand side of Ax=r linear equation
c nelm   - connectivity of a matrix
c nrhs   - array to partition r for multiple degrees of freedom
c nelmdg - position on nelm of diagonal element
c indx   - storage(idof) needed for algorithm(integer array)
c idof   - degree of freedom for problem
c diag   - storage(idof,idof) needed for diagonal element of a
c diagi  - storage(idof,idof) needed for diagonal inverse element of a
c rdum   - storage(idof) needed for algorithm
c iflg   - flag used for type of normalization
c fdum2  - sum-squared of the residuals
C***********************************************************************

      implicit none

      integer i,j,k,l,i1,i2,ii,ipos,neq,idof
      integer idiag,neqp1,iflg,idofdiag
      integer nelm(*),nmat(*),nrhs(*),nelmdg(*),indx(*)
      real*8 diag(idof,*),diagi(idof,*),rdum(*),a(*),r(*),det,fdum2
      real*8 diag_tol
      parameter (diag_tol=1.d-20)
      
      neqp1=neq+1
      fdum2 = 0.0
c     ****** start loop on rows *************************************** 
      if(iflg.eq.0) then     
         do i=1,neq
            idiag=nelmdg(i)-neqp1
c     ****************** store diagonal ****************************** 
            do j=1,idof
               do k=1,idof
                  diag(j,k)=a(nmat((j-1)*idof+k)+idiag)
               enddo
c     if any diagonal is small don't normalize     
            enddo
c     ****************** form inverse of diagonal ********************* 
c     diagonal is in diag(i,j),diagi(i,j) contains its inverse
            do j=1,idof
               do k=1,idof
                  diagi(j,k)=0.0
               enddo
               diagi(j,j)=1.0
            enddo
            call ludcmp0(diag,idof,idof,indx,det)
            if(indx(1).lt.0) then
               indx(1)=-i
               return
            endif
            do j=1,idof
               call lubksb0(diag,idof,idof,indx,diagi(1,j))
            enddo
c     ****************** normalize equation ************************** 
            i1=nelm(i)+1
            i2=nelm(i+1)
            do ii=i1,i2
c     skip over diagonal position
               if(ii.ne.idiag+neqp1)then
c     identify row elements
                  ipos=ii-neqp1
                  do j=1,idof
                     do k=1,idof
                        diag(j,k)=a(nmat((j-1)*idof+k)+ipos)
                        a(nmat((j-1)*idof+k)+ipos)=0.0
                     enddo
                  enddo
c     multiply row elements by inverse of diagonal
                  do j=1,idof
                     do k=1,idof
                        do l=1,idof
                           a(nmat((j-1)*idof+k)+ipos)=
     $                          + a(nmat((j-1)*idof+k)+ipos)+
     &                          diagi(j,l)*diag(l,k)
                        enddo
                     enddo
                  enddo
               else 
c     set diagonal to identity matrix
                  do j=1,idof
                     do k=1,idof
                        a(nmat((j-1)*idof+k)+idiag)=0.0
                     enddo
                     a(nmat((j-1)*idof+j)+idiag)=1.0
                  enddo
               endif
            enddo
            do j=1,idof
               rdum(j)=r(i+nrhs(j))
               r(i+nrhs(j))=0.0
            enddo
            do j=1,idof
               do k=1,idof
                  r(i+nrhs(j))=r(i+nrhs(j))+diagi(j,k)*rdum(k)         
               enddo
               fdum2 =fdum2 + r(i+nrhs(j))*r(i+nrhs(j))
            enddo
c     ****** end loop on rows ***************************************   
         enddo
c     ************* uncoupled normalization *************************
      elseif(iflg.ne.0) then
         do i=1,neq
            idiag=nelmdg(i)-neqp1
            i1=nelm(i)+1
            i2=nelm(i+1)
            do j=1,idof
               idofdiag=(j-1)*idof+j
               diagi(1,1)=1.0/a(nmat(idofdiag)+idiag)
               do k=(j-1)*idof+1,j*idof
                  do ii=i1,i2
                     a(nmat(k)+ii-neqp1) = a(nmat(k)+ii-neqp1)/
     &                    diagi(1,1)
                  enddo
               enddo
               r(i+nrhs(j))=r(i+nrhs(j))/diagi(1,1)
               fdum2 =fdum2 + r(i+nrhs(j))*r(i+nrhs(j))
            enddo
         enddo
      endif

      return 
      end
