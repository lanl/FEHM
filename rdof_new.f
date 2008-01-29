      subroutine rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *     ,tollr,irb,iirb,npvt,gmres,dum,piv
     *     ,h,c,ss,g,y,iter,iback,itype,iptty,maxor
     *     ,overf,irdof,icoupl,mcount,sto5,accm)       
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
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rdof_new.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!!D2    Rev 2.3   14 Nov 2001 13:12:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Fri Sep 26 15:13:34 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.8   Mon Apr 14 12:44:10 1997   gaz
CD2 minor correction
CD2 
CD2    Rev 1.7   Fri May 31 15:09:24 1996   gaz
CD2 added different rdof algorithm for 6dof
CD2 
CD2    Rev 1.6   Tue Jan 16 14:31:18 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.5   Thu Jan 11 13:17:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   12/13/95 08:44:12   gaz
CD2 added 4 to 2 and 6 to 2 capabilities
CD2 
CD2    Rev 1.3   11/15/95 15:17:46   gaz
CD2 minor corrections plus added 4> 2 dof
CD2 
CD2    Rev 1.2   04/25/95 09:06:40   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.1   03/23/95 19:10:56   gaz
CD2 gaz update reduced degree of freedon solver
CD2 
CD2    Rev 1.0   01/28/95 13:55:20   llt
CD2 water balance equation was modified
CD2
c
c 11/22/94 gaz put in the solve new calls 
c check sto4
c combined two and three dof
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
CD4  These algorithms assume normalized equations (ie aii=1.0)
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      implicit none

      integer nmat_save(16)
      integer neq,nmat(*),nb(*),nrhs(*),nelm(*),nop(*)
      integer north,irb(*),iirb(*),npvt(*),iter,iback,itype 
      integer iptty,maxor,irdof,icoupl,mcount
      integer i, i1, i2, iterg, j, k, kb, kk
      integer neqp1
      real*8  a(*),b(*),bp(*),tollr,gmres(*)
      real*8  dum(*),sto5(neq,*),piv(neq,*)
      real*8  h(maxor,*),c(*),ss(*),g(*),y(*),overf
      real*8  sum, sum1, sum2, sum3, sum4, sum5, sum6 
      real*8  fdum2,toll_test
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
      neqp1=neq+1
      iterg=0
      if(mcount.eq.0) then
         if(irdof.eq.2.and.itype.eq.3) then
            if (icoupl.ge.1) then
               do i=1,neq
                  sto5(i,1)=bp(nrhs(1)+i)
                  sto5(i,2)=bp(nrhs(2)+i)
                  sto5(i,3)=bp(nrhs(3)+i)
               enddo
            endif
            nmat_save(3)=nmat(3)
            nmat_save(4)=nmat(4)
            nmat(3)=nmat(4)
            nmat(4)=nmat(5)
c     
c     reduce degree of freedom from 3 to 2
c     
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)     
            iterg=iterg+iter
c     
c     back substitute to find bp(nrhs(3)+i)
c     
            do i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum=bp(nrhs(3)+i)
               do j=i1,i2
                  k=j-neqp1
                  kb=nelm(j)
                  sum=sum-a(nmat(7)+k)*bp(nrhs(1)+kb)
     &                 -a(nmat(8)+k)*bp(nrhs(2)+kb)
               enddo
               bp(nrhs(3)+i)=sum
            enddo
            nmat(3)=nmat_save(3)
            nmat(4)=nmat_save(4)
c     
         endif
         if(irdof.eq.1.and.itype.eq.3) then     
c     
c     reduce degree of freedom from 3 to 1
c     
c     *** may not need this ******
            if (icoupl.ge.1) then
               do i=1,neq
                  sto5(i,1)=bp(nrhs(1)+i)
                  sto5(i,2)=bp(nrhs(2)+i)
                  sto5(i,3)=bp(nrhs(3)+i)
               enddo
            endif
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)     
            iterg=iterg+iter
c     
c     back substitute to find bp(nrhs(2)+i)
c     
            do i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum=bp(nrhs(2)+i)
               do j=i1,i2
                  k=j-neqp1
                  kb=nelm(j)
                  sum=sum-a(nmat(4)+k)*bp(nrhs(1)+kb)
               enddo
               bp(nrhs(2)+i)=sum
            enddo        
c     
c     back substitute to find bp(nrhs(3)+i)
c     
            do i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum=bp(nrhs(3)+i)
               doj=i1,i2
               k=j-neqp1
               kb=nelm(j)
               sum=sum-a(nmat(7)+k)*bp(nrhs(1)+kb)
c     &        -a(nmat(8)+k)*bp(nrhs(2)+kb)
            enddo
            bp(nrhs(3)+i)=sum
         enddo
c     
      endif      
      if(irdof.ge.1.and.itype.eq.2) then
         if (icoupl.ge.1) then
            do i=1,neq
               sto5(i,1)=bp(nrhs(1)+i)
               sto5(i,2)=bp(nrhs(2)+i)
            enddo
         endif
c     
c     assumes diagonal is an identity matrix
c     
c     reduce degree of freedom from 2 to 1
c     
         call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *        ,tollr,irb,iirb,npvt,gmres,dum,piv
     *        ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)     
         iterg=iterg+iter
c     
c     first back substitute
c     
         do i=1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            sum=bp(nrhs(2)+i)
            do j=i1,i2
               sum=sum-a(nmat(3)+j-neqp1)*bp(nrhs(1)+nelm(j))
            enddo
            bp(nrhs(2)+i)=sum
         enddo
c     
      endif
      if(irdof.ne.0.and.itype.eq.4) then     
         if (icoupl.ge.1) then
            do i=1,neq
               sto5(i,1)=bp(nrhs(1)+i)
               sto5(i,2)=bp(nrhs(2)+i)
               sto5(i,3)=bp(nrhs(3)+i)
               sto5(i,4)=bp(nrhs(4)+i)
            enddo
         endif
         nmat_save(3)=nmat(3)
         nmat_save(4)=nmat(4)
         nmat(3)=nmat(5)
         nmat(4)=nmat(6)
c     
c     reduce degree of freedom from 4 to 2
c     
         call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *        ,tollr,irb,iirb,npvt,gmres,dum,piv
     *        ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)     
         iterg=iterg+iter
c     
c     back substitute to find bp(nrhs(3)+i) and bp(nrhs(4)+i)
c     
         do i=1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            sum1=bp(nrhs(3)+i)
            sum2=bp(nrhs(4)+i)
            do j=i1,i2
               k=j-neqp1
               kb=nelm(j)
               sum1=sum1-a(nmat(9)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(10)+k)*bp(nrhs(2)+kb)
               sum2=sum2-a(nmat(13)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(14)+k)*bp(nrhs(2)+kb)
            enddo
            bp(nrhs(3)+i)=sum1
            bp(nrhs(4)+i)=sum2
         enddo        
         nmat(3)=nmat_save(3)
         nmat(4)=nmat_save(4)
c     
      endif      
c     
      if(irdof.eq.1.and.itype.eq.6) then
         if (icoupl.ge.1) then
            do i=1,neq
               sto5(i,1)=bp(nrhs(1)+i)
               sto5(i,2)=bp(nrhs(2)+i)
               sto5(i,3)=bp(nrhs(3)+i)
               sto5(i,4)=bp(nrhs(4)+i)
               sto5(i,5)=bp(nrhs(5)+i)
               sto5(i,6)=bp(nrhs(6)+i)
            enddo
         endif
         nmat_save(1)=nmat(3)
         nmat_save(2)=nmat(4)
         nmat(3)=nmat(7)
         nmat(4)=nmat(8)
c     
c     reduce degree of freedom from 6 to 2
c     
         call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *        ,tollr,irb,iirb,npvt,gmres,dum,piv
     *        ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
         iterg=iterg+iter
c     
c     back substitute to find bp(nrhs(3)+i) and bp(nrhs(4)+i)
c     back substitute to find bp(nrhs(5)+i) and bp(nrhs(6)+i)
c     
         do i=1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            sum1=bp(nrhs(3)+i)
            sum2=bp(nrhs(4)+i)
            sum3=bp(nrhs(5)+i)
            sum4=bp(nrhs(6)+i)
            do j=i1,i2
               k=j-neqp1
               kb=nelm(j)
               sum1=sum1-a(nmat(13)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(14)+k)*bp(nrhs(2)+kb)
               sum2=sum2-a(nmat(19)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(20)+k)*bp(nrhs(2)+kb)
               sum3=sum3-a(nmat(25)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(26)+k)*bp(nrhs(2)+kb)
               sum4=sum4-a(nmat(31)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(32)+k)*bp(nrhs(2)+kb)
            enddo
            bp(nrhs(3)+i)=sum1
            bp(nrhs(4)+i)=sum2
            bp(nrhs(5)+i)=sum3
            bp(nrhs(6)+i)=sum4
         enddo
         nmat(3)=nmat_save(1)
         nmat(4)=nmat_save(2)
c     
      endif
      if(irdof.eq.2.and.itype.eq.6) then
         if (icoupl.ge.1) then
            do i=1,neq
               sto5(i,1)=bp(nrhs(1)+i)
               sto5(i,2)=bp(nrhs(2)+i)
               sto5(i,3)=bp(nrhs(3)+i)
               sto5(i,4)=bp(nrhs(4)+i)
               sto5(i,5)=bp(nrhs(5)+i)
               sto5(i,6)=bp(nrhs(6)+i)
            enddo
         endif
         nmat_save(1)=nmat(5)
         nmat_save(2)=nmat(6)
         nmat_save(3)=nmat(7)
         nmat_save(4)=nmat(8)
         nmat_save(5)=nmat(9)
         nmat_save(6)=nmat(10)
         nmat_save(7)=nmat(11)
         nmat_save(8)=nmat(12)
         nmat_save(9)=nmat(13)
         nmat_save(10)=nmat(14)
         nmat_save(11)=nmat(15)
         nmat_save(12)=nmat(16)
         nmat(5)=nmat_save(3)
         nmat(6)=nmat_save(4)
         nmat(7)=nmat_save(5)
         nmat(8)=nmat_save(6)
         nmat(9)=nmat_save(9)
         nmat(10)=nmat_save(10)
         nmat(11)=nmat_save(11)
         nmat(12)=nmat_save(12)
         nmat(13)=nmat(19)
         nmat(14)=nmat(20)
         nmat(15)=nmat(21)
         nmat(16)=nmat(22)
c     
c     reduce degree of freedom from 6 to 4
c     
         call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *        ,tollr,irb,iirb,npvt,gmres,dum,piv
     *        ,h,c,ss,g,y,iter,iback,4,iptty,maxor,accm)
         iterg=iterg+iter
c     
c     back substitute to find bp(nrhs(5)+i) and bp(nrhs(6)+i)
c     
         do i=1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            sum3=bp(nrhs(5)+i)
            sum4=bp(nrhs(6)+i)
            do j=i1,i2
               k=j-neqp1
               kb=nelm(j)
               sum3=sum3-a(nmat(25)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(26)+k)*bp(nrhs(2)+kb)
     &              -a(nmat(27)+k)*bp(nrhs(3)+kb)
     &              -a(nmat(28)+k)*bp(nrhs(4)+kb)
               sum4=sum4-a(nmat(31)+k)*bp(nrhs(1)+kb)
     &              -a(nmat(32)+k)*bp(nrhs(2)+kb)
     &              -a(nmat(33)+k)*bp(nrhs(3)+kb)
     &              -a(nmat(34)+k)*bp(nrhs(4)+kb)
            enddo
            bp(nrhs(5)+i)=sum3
            bp(nrhs(6)+i)=sum4
         enddo
         nmat(5)=nmat_save(1)
         nmat(6)=nmat_save(2)
         nmat(7)=nmat_save(3)
         nmat(8)=nmat_save(4)
         nmat(9)=nmat_save(5)
         nmat(10)=nmat_save(6)
         nmat(11)=nmat_save(7)
         nmat(12)=nmat_save(8)
         nmat(13)=nmat_save(9)
         nmat(14)=nmat_save(10)
         nmat(15)=nmat_save(11)
         nmat(16)=nmat_save(12)
c     
      endif
      if(itype.eq.2) then
c     
c     successive over-relaxation
c     
         toll_test = 1.d-12 

         do  j=1,icoupl
            fdum2=0.0
            do  i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum1=sto5(i,1)
               sum2=sto5(i,2)
               if(abs(sum1).ge.toll_test) then
                  do  k=i1,i2
                     kk=k-neqp1
                     kb=nelm(k)
                     sum1=sum1-a(nmat(1)+kk)*bp(nrhs(1)+kb)
     &                    -a(nmat(2)+kk)*bp(nrhs(2)+kb)
                  enddo
               endif
               if(abs(sum2).ge.toll_test) then	   	 
                  do  k=i1,i2
                     kk=k-neqp1
                     kb=nelm(k)
                     sum2=sum2-a(nmat(3)+kk)*bp(nrhs(1)+kb)
     &                    -a(nmat(4)+kk)*bp(nrhs(2)+kb)
                  enddo
               endif
               bp(nrhs(1)+i)=bp(nrhs(1)+i)+overf*sum1
               bp(nrhs(2)+i)=bp(nrhs(2)+i)+overf*sum2
               fdum2=fdum2 + abs(sum1)+abs(sum2)
            enddo
            if(fdum2.le.tollr) go to 500
         enddo
      else if(itype.eq.3) then
         do  j=1,icoupl
            fdum2=0.0d00
            do  i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum1=sto5(i,1)
               sum2=sto5(i,2)
               sum3=sto5(i,3)
               do  k=i1,i2
                  kk=k-neqp1
                  kb=nelm(k)
                  sum1=sum1-a(nmat(1)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(2)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(3)+kk)*bp(nrhs(3)+kb)
                  sum2=sum2-a(nmat(4)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(5)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(6)+kk)*bp(nrhs(3)+kb)
                  sum3=sum3-a(nmat(7)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(8)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(9)+kk)*bp(nrhs(3)+kb)
               enddo
               bp(nrhs(1)+i)=bp(nrhs(1)+i)+overf*sum1
               bp(nrhs(2)+i)=bp(nrhs(2)+i)+overf*sum2
               bp(nrhs(3)+i)=bp(nrhs(3)+i)+overf*sum3
               fdum2=fdum2 + abs(sum1)+abs(sum2)+abs(sum3)
            enddo
            if(fdum2.le.tollr) go to 500
         enddo
      else if(itype.eq.4) then
         do  j=1,icoupl
            fdum2=0.0d00
            do  i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum1=sto5(i,1)
               sum2=sto5(i,2)
               sum3=sto5(i,3)
               sum4=sto5(i,4)
               do  k=i1,i2
                  kk=k-neqp1
                  kb=nelm(k)
                  sum1=sum1-a(nmat(1)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(2)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(3)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(4)+kk)*bp(nrhs(4)+kb)
                  sum2=sum2-a(nmat(5)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(6)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(7)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(8)+kk)*bp(nrhs(4)+kb)
                  sum3=sum3-a(nmat(9)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(10)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(11)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(12)+kk)*bp(nrhs(4)+kb)
                  sum4=sum4-a(nmat(13)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(14)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(15)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(16)+kk)*bp(nrhs(4)+kb)
               enddo
               bp(nrhs(1)+i)=bp(nrhs(1)+i)+overf*sum1
               bp(nrhs(2)+i)=bp(nrhs(2)+i)+overf*sum2
               bp(nrhs(3)+i)=bp(nrhs(3)+i)+overf*sum3
               bp(nrhs(4)+i)=bp(nrhs(4)+i)+overf*sum4
               fdum2=fdum2+abs(sum1)+abs(sum2)+abs(sum3)+abs(sum4)
            enddo
            if(fdum2.le.tollr) go to 500
         enddo
      else if(itype.eq.6) then
         do  j=1,icoupl
            fdum2=0.0d00
            do  i=1,neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               sum1=sto5(i,1)
               sum2=sto5(i,2)
               sum3=sto5(i,3)
               sum4=sto5(i,4)
               sum5=sto5(i,5)
               sum6=sto5(i,6)
               do  k=i1,i2
                  kk=k-neqp1
                  kb=nelm(k)
                  sum1=sum1-a(nmat(1)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(2)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(3)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(4)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(5)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(6)+kk)*bp(nrhs(6)+kb)
                  sum2=sum2-a(nmat(7)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(8)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(9)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(10)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(11)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(12)+kk)*bp(nrhs(6)+kb)
                  sum3=sum3-a(nmat(13)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(14)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(15)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(16)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(17)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(18)+kk)*bp(nrhs(6)+kb)
                  sum4=sum4-a(nmat(19)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(20)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(21)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(22)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(23)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(24)+kk)*bp(nrhs(6)+kb)
                  sum5=sum5-a(nmat(25)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(26)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(27)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(28)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(29)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(30)+kk)*bp(nrhs(6)+kb)
                  sum6=sum6-a(nmat(31)+kk)*bp(nrhs(1)+kb)
     &                 -a(nmat(32)+kk)*bp(nrhs(2)+kb)
     &                 -a(nmat(33)+kk)*bp(nrhs(3)+kb)
     &                 -a(nmat(34)+kk)*bp(nrhs(4)+kb)
     &                 -a(nmat(35)+kk)*bp(nrhs(5)+kb)
     &                 -a(nmat(36)+kk)*bp(nrhs(6)+kb)
               enddo
               bp(nrhs(1)+i)=bp(nrhs(1)+i)+overf*sum1
               bp(nrhs(2)+i)=bp(nrhs(2)+i)+overf*sum2
               bp(nrhs(3)+i)=bp(nrhs(3)+i)+overf*sum3
               bp(nrhs(4)+i)=bp(nrhs(4)+i)+overf*sum4
               bp(nrhs(5)+i)=bp(nrhs(5)+i)+overf*sum5
               bp(nrhs(6)+i)=bp(nrhs(6)+i)+overf*sum6
               fdum2=fdum2+abs(sum1)+abs(sum2)+abs(sum3)
     &              +abs(sum4)+abs(sum5)+abs(sum6)
            enddo
            if(fdum2.le.tollr) go to 500
         enddo
      endif
c     
 500  continue
      return
      endif
      return
      end
