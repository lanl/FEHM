
      subroutine combine_a(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *     ,tollr,irb,iirb,npvt,gmres,dum,piv
     *     ,h,c,ss,g,y,iter,iback,itype,iptty,maxor
     *     ,overf,irdof,icoupl,mcount,sto5,a_save,nelmdg,accm,mdof_sol)
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
CD1  This subroutine reduces the 3n*3n matrix into an n*n matrix or
CD1  into a 2n*2n matrix using the RDOF or IRDOF algorithms.  
CD1  Used when a mtrix is not normalized.                       
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 FEHM Version 2.20
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/combine_a.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 08:56:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2 
CD2 Split subroutine out of rdof_new 08-Feb-02
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  These algorithms do not assume normalized equations (ie aii=1.0)
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      implicit none

      integer nmat_save(16)
      integer neq,nmat(*), nb(*), nrhs(*), nelm(*), nop(*), nelmdg(*)
      integer north, irb(*), iirb(*), npvt(*), iter,iback, itype 
      integer iptty, maxor, irdof, icoupl, mcount, nmat_calc
      integer i, i1, i2, idiag, iterg, j, k, kb, kk
      integer neqp1, mdof_sol
      real*8  a(*),b(*),bp(*),tollr,gmres(*),a_save(*)
      real*8  dum(*),sto5(neq,*),piv(neq,*)
      real*8  h(maxor,*),c(*),ss(*),g(*),y(*),overf
      real*8  sum, sum1, sum2, sum3, sum4, a_piv
      real*8  alm, ali
      real*8  a1, a11, a11i, a12, a12a22i, a12i, a2, a21, a21i 
      real*8  a22, a22i, a3, a4, b1, b2
      real*8  fdum2
      character*4 accm

c neq       - number of equations (in block form)
c nmat      - position array for matrix a (for multiple degrees of freedom)
c nb        - position array for matrix b (for multiple degrees of freedom)
c nrhs      - position array for right hand side vector
c nelm      - connectivity array for solution matrix
c nop       - connectivity array for solution matrix
c north     - number of othogonalizations for gmres acceleration
c irb       - array of renumbered nodes
c iirb      - inverse of irb
c npvt      - array of positions of diagonal for the solution array
c iter      - maximum number of iterations
c iback     - iback=0 means caluclate the LU decomposition
c           - iback ne 0 means use existing LU decomposition
c itype     - process type(3-air-water-heat,2-water-heat)
c iptty     - unit number for printout
c maxor     - maximum number of orthogonalizations
c irdof     - parameter defining reduced degree of freedom procedures
c           - irdof=0, no action
c           - irdof=1, reduce 3dof to 2dof,reduce 2dof to 1dof
c           - reduce 4dof to 2dof
c           - irdof=2, reduce 3dof to 1dof,reduce 2dof to 1dof
c           - reduce 4dof to 2dof
c icoupl    - maximum number of SOR iterations 
c mcount    - parameter for SOR iterations mcount=0, no SOR iterations before solve_new
c           - mcount ne 0,
c a         - solution matrix
c b         - incomplete LU factorization matrix
c bp        - right hand side
c tollr     - tolerance for solution
c gmres     - scratch storage for GMRES acceleration
c sto1      - scratch storage
c sto2      - scratch storage
c sto3      - scratch storage
c sto4      - scratch storage
c piv       - array of pivot values for the LU factor matrix
c h         - array array used in GMRES acceleration
c c         - array array used in GMRES acceleration
c ss        - array array used in GMRES acceleration
c g         - array array used in GMRES acceleration
c y         - array array used in GMRES acceleration
c overf     - overrelaxation factor for SOR iteration 
c
c local variables
c 
c neqp1     - number of equations(neq) +1
c iterg     - number of iteratins
c i         - do loop index
c i1        - do loop starting value
c i2        - do loop ending value 
c j         - position in connectivity matrix
c k         - position in solution matrix
c kb        - neighbor node of node i
c sum       - intermediate sum            
c nmat_save - positions of nmat array on input
c
c save locations of submatrices of a
c no reorder here,independent of order in gensl4
c

c
c***     linear algebra     ***
c
      alm(a1,a2,b1,b2)=a1*b1+a2*b2
      ali(a1,a2,a3,a4,b1)=b1/(a1*a4-a2*a3)
c
      neqp1=neq+1
      iterg=0
      if(mcount.eq.0) then
      if(irdof.ge.1.and.itype.eq.2) then
      if (icoupl.ge.1) then
        do i=1,neq
          sto5(i,1)=bp(nrhs(1)+i)
          sto5(i,2)=bp(nrhs(2)+i)
        enddo
      endif
c
c       first reduce degree of freedom
c
        do i=1,neq
          i1=nelm(i)+1
          i2=nelm(i+1)
          idiag = nelmdg(i)
          sum=bp(nrhs(1)+i)
          a12=a(nmat(2)-neqp1+idiag)
          if(a12.ne.0.0) then
           a22=a(nmat(4)-neqp1+idiag)
           a12a22i=a12/a22
          else
           a12a22i=0.0      
          endif
          do j=i1,i2
            a_save(j-neqp1)=a(nmat(1)+j-neqp1)-
     &      a12a22i*a(nmat(3)+j-neqp1)
          enddo
          sum=sum-a12a22i*bp(nrhs(2)+i)
          bp(i+nrhs(1))=sum
c normalize reduced degree of freedom
          a11=a_save(idiag-neqp1)
          do j=i1,i2
            a_save(j-neqp1)=a_save(j-neqp1)/a11
          enddo
          bp(i+nrhs(1))=bp(i+nrhs(1))/a11
        enddo
c
c
c assumes diagonal is an identity matrix
c
c       reduce degree of freedom from 2 to 1
c
        call solve_new(neq,a_save,b,bp,nmat,nb,nrhs,nelm,nop
     &   ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     &   ,h,c,ss,g,y,iter,iback,4,iptty,maxor,accm)
        iterg=iterg+iter
c
c       first back substitute
c
        do i=1,neq
          i1=nelm(i)+1
          i2=nelm(i+1)
          sum=bp(nrhs(2)+i)
          idiag = nelmdg(i)
          a22=a(idiag+nmat(4)-neqp1)
          a12=a(idiag+nmat(2)-neqp1)
           if(a22.gt.a12) then
            nmat_calc=nmat(3)
            a_piv=a22
           else if(a12.ge.a22) then
            nmat_calc=nmat(1)
            a_piv=a12
           endif
           do j=i1,i2
             sum=sum-a(nmat_calc+j-neqp1)*bp(nelm(j)+nrhs(1))
           enddo
           bp(nrhs(2)+i)=sum/a_piv
        enddo
c
      endif
      if(itype.eq.2) then
c
c
c
      do  j=1,icoupl
        fdum2=0.0
        do  i=1,neq
          i1=nelm(i)+1
          i2=nelm(i+1)
          sum1=sto5(i,1)
          sum2=sto5(i,2)
          do  k=i1,i2
            kk=k-neqp1
            kb=nelm(k)
            sum1=sum1-a(nmat(1)+kk)*bp(nrhs(1)+kb)
     &        -a(nmat(2)+kk)*bp(nrhs(2)+kb)
            sum2=sum2-a(nmat(3)+kk)*bp(nrhs(1)+kb)
     &        -a(nmat(4)+kk)*bp(nrhs(2)+kb)
          enddo
            fdum2=fdum2 + abs(sum1)+abs(sum2)
c remove diagonal term
            idiag = nelmdg(i)
            sum1 = sum1 + a(idiag-neqp1+nmat(1))*bp(nrhs(1)+i)
     &        + a(idiag+nmat(2)-neqp1)*bp(nrhs(2)+i)
            sum2 = sum2 + a(idiag-neqp1+nmat(3))*bp(nrhs(1)+i)
     &        + a(idiag+nmat(4)-neqp1)*bp(nrhs(2)+i)
             a11 = a(idiag-neqp1+nmat(1))
             a22 = a(idiag-neqp1+nmat(4))
             a12 = a(idiag-neqp1+nmat(2))
             a21 = a(idiag-neqp1+nmat(3))
             a11i=ali(a11,a12,a21,a22,a22)
             a12i=-ali(a11,a12,a21,a22,a12)
             a21i=-ali(a11,a12,a21,a22,a21)
             a22i=ali(a11,a12,a21,a22,a11)
             bp(nrhs(1)+i)=alm(a11i,a12i,sum1,sum2)*overf +
     &                   (1.0-overf)*bp(nrhs(1)+i)
             bp(nrhs(2)+i)=alm(a21i,a22i,sum1,sum2)*overf +
     &                   (1.0-overf)*bp(nrhs(2)+i)
       enddo

        if(fdum2.le.tollr) go to 500
      enddo
      endif
c
500   continue
      return
      endif
      return
      end
