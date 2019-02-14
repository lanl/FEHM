      subroutine solve_dual(neq_primary,neq,a,b,bp,na,nb
     &  ,nrhs,ncon,ncon_primary,nop,inorth,epn,irb,iirb
     &  ,npvt,stor1,dum,piv,h,c,s,g,y,iter,iback,idof,iptty
     &  ,maxor,igdpm,max_layers,ngdpm_layers,nelmdg,accm,
     &   mdof_sol)
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
!D1 To perform pre-conditioned conjugate gradient solution of a set of
!D1 linear, algebraic equations. Generalized Dual Porosity.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: ?, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/solve_dual.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:16:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 
!D4 
!***********************************************************************

      implicit none
      integer idof, neq_primary, neq, maxor, mdof_sol, iarr_size
      real*8 a(*),b(*),bp(*)
      real*8 stor1(*),piv(neq,*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*), epn, anorm
      real*8 dum(*)
      integer ncon(*),nop(*)
      integer irb(*),iirb(*),nelmdg(*),npvt(*)
      integer igdpm(*),ngdpm_layers(0:*)
      integer na(*),nb(*),nrhs(*)
      integer ncon_primary(*)
      integer inorth, iter, iback, iptty
      integer isolve
      integer ibcgs
c     variables for general dual porosity
      integer max_layers
      integer, allocatable :: iiconn(:),iim1conn(:)
      integer, allocatable :: iip1conn(:), iarr(:)
      real*8, allocatable :: a_primary(:)
      integer, allocatable :: idum(:),na_primary(:)
      real*8, allocatable :: adum(:)
      save idum,adum,iarr
      integer i,ic,ii,i1,i2
      integer ncon_size, ncon_primary_size, na_ratio 
      character*4 accm
c     
c     calculate maximum number of layers
c     
      ncon_primary_size =ncon_primary(neq_primary+1)
     &     -(neq_primary+1)
      ncon_size = ncon(neq+1) - (neq+1)
      allocate(iiconn(max_layers))
      allocate(iim1conn(max_layers))
      allocate(iip1conn(max_layers))

c     first reduce solution to neq by neq system

      call general_dual(1,neq_primary,neq,a,bp,na,nrhs,ncon,
     &     ngdpm_layers,nelmdg,igdpm,idof,
     &     iiconn,iim1conn,iip1conn)
c     
      allocate(na_primary(idof*idof))
      do i =1,idof*idof
         na_ratio = na(i)/ncon_size
         na_primary(i) = na_ratio*ncon_primary_size
      enddo
      iarr_size = ncon_primary_size*idof*idof
      if(.not.allocated(iarr)) then
c         allocate(iarr(ncon_primary_size*idof*idof))
         allocate(iarr(iarr_size))
         allocate(a_primary(1), adum(1), idum(1))
         call simplify_gdpm(4,neq,neq_primary,ncon
     &        ,idum,ncon_primary,a,adum,a_primary,na
     &        ,na_primary,ngdpm_layers,igdpm,iarr,idof)
         deallocate(adum)
         allocate (adum(iarr_size))
      endif

c     call neq by neq solution                     

c      call solve_new(neq_primary,a(iarr),b,bp,na_primary, gaz 102416
      adum(1:iarr_size) = a(iarr)
      call solve_new(neq_primary,adum,b,bp,na_primary,
     &     nb,nrhs,ncon_primary,nop,inorth,
     &     epn,irb,iirb,npvt,stor1,dum,piv,
     &     h,c,s,g,y,iter,iback,idof,iptty,maxor,accm)     

      deallocate(na_primary)

c     back out full solution                      

      call general_dual(2,neq_primary,neq,a,bp,na,nrhs,ncon,
     &     ngdpm_layers,nelmdg,igdpm,idof,
     &     iiconn,iim1conn,iip1conn)

      deallocate(iiconn)
      deallocate(iim1conn)
      deallocate(iip1conn)

      return
      end
      subroutine general_dual(iflg,neq_primary,neq,a,bp,
     &     na,nrhs,nelm_gen,
     &     ngdpm_layers,nelmdg_gen,igdpm,idof,
     &     iiconn,iim1conn,iip1conn)
c     
c     perform condensation and back substitution for  
c     generalized dual porosity
c     
      implicit none
c     
      integer iflg,i,j,ii,neq_primary,neq,neqp1
      integer ip1,ip2,ip3,ip1ip1,ip1i,iip1, ip1p1
      integer ip1m1,ip2m1, iip2p1, iip1p1, iim1, idof
      integer na(*),nrhs(*),nelm_gen(*),igdpm(*)
      integer ngdpm_layers(0:*)                  
      integer nelmdg_gen(*)
      integer iiconn(*)
      integer iim1conn(*)
      integer iip1conn(*)
      integer jj,l,k
      real*8 a(*),bp(*)
      real*8 ad11,ad12,ad21,ad22,det
      real*8, allocatable :: aii(:,:),aiip1(:,:),aiim1(:,:)
      real*8, allocatable :: aip1i(:,:),aip1ip1(:,:)
      real*8, allocatable :: bi(:),bip1(:),bim1(:)
c     
      neqp1 = neq + 1
c     
      if(iflg.eq.1) then
         if(idof.eq.1) then
c     condense solution 1 dof
            do i=1,neq_primary
               if(ngdpm_layers(igdpm(i)).gt.0) then
c     set up 1-d connectivity
c     ip1 = i-1, ip2 =i, ip3 = i+1
                  ip2 = i
                  ip3= nelm_gen(nelm_gen(ip2+1))
                  iiconn(1)=ip2
                  iip1conn(1)=ip3
                  do j=1,ngdpm_layers(igdpm(i))-1
                     ip2= nelm_gen(nelm_gen(ip2+1))
                     ip3= nelm_gen(nelm_gen(ip2)+3)
                     iiconn(j+1)=ip2
                     iip1conn(j+1)=ip3
                  enddo
c     perform forward substitution
c     a(ii) = a(ii) - a(i,i+1)/a(i+1,i+1)*a(i+1,i)
c     b(i)  =  b(i) - a(i,i+1)/a(i+1,i+1)*b(i+1)
                  do j = ngdpm_layers(igdpm(i)),1,-1
                     ip3 = iip1conn(j)
                     ip2 = iiconn(j)
                     ii = nelmdg_gen(ip2) - neqp1 
                     ip1ip1 = nelmdg_gen(ip3) - neqp1 
c     iip1 = ii + 1
                     iip1 = nelm_gen(ip2+1) - neqp1 
                     ip1i  = ip1ip1 - 1
                     a(ii) = a(ii) - a(iip1)*a(ip1i)/a(ip1ip1)
                     bp(ip2) = bp(ip2) - a(iip1)*bp(ip3)/a(ip1ip1)
                  enddo
               endif
            enddo
         else if (idof.eq.2) then
c     condense solution 2 dof
            allocate(aii(idof,idof),aiip1(idof,idof),aiim1(idof,idof))
            allocate(aip1i(idof,idof),aip1ip1(idof,idof))
            allocate(bi(idof),bip1(idof),bim1(idof))
            do i=1,neq_primary
               if(ngdpm_layers(igdpm(i)).gt.0) then
c     set up 1-d connectivity
c     ip1 = i-1, ip2 =i, ip3 = i+1
                  ip2 = i
                  ip3= nelm_gen(nelm_gen(ip2+1))
                  iiconn(1)=ip2
                  iip1conn(1)=ip3
                  do j=1,ngdpm_layers(igdpm(i))-1
                     ip2= nelm_gen(nelm_gen(ip2+1))
                     ip3= nelm_gen(nelm_gen(ip2)+3)
                     iiconn(j+1)=ip2
                     iip1conn(j+1)=ip3
                  enddo
c     perform forward substitution
c     a(ii) = a(ii) - a(i,i+1)/a(i+1,i+1)*a(i+1,i)
c     b(i)  =  b(i) - a(i,i+1)/a(i+1,i+1)*b(i+1)
                  do j = ngdpm_layers(igdpm(i)),1,-1
                     ip3 = iip1conn(j)
                     ip2 = iiconn(j)
                     ii = nelmdg_gen(ip2) - neqp1 
                     ip1ip1 = nelmdg_gen(ip3) - neqp1 
c     iip1 = ii + 1
                     iip1 = nelm_gen(ip2+1) - neqp1 
                     ip1i  = ip1ip1 - 1
                     jj = 0
                     do l = 1,idof
                        do k = 1,idof
                           jj=jj+1
                           aii(l,k) = a(ii+na(jj))
                           aiip1(l,k) = a(iip1+na(jj))
                           aip1i(l,k) = a(ip1i+na(jj))
                           aip1ip1(l,k) = a(ip1ip1+na(jj))
                        enddo
                        bi(l) = bp(ip2+nrhs(l))
                        bip1(l) = bp(ip3+nrhs(l))
                     enddo
c     find inverse of aip1ip1
                     ad11 = aip1ip1(1,1)
                     ad12 = aip1ip1(1,2)
                     ad21 = aip1ip1(2,1)
                     ad22 = aip1ip1(2,2)
                     det = ad11*ad22-ad12*ad21
c     put inverse back into aip1ip1
                     aip1ip1(1,1) = ad22/det
                     aip1ip1(1,2) = -ad12/det
                     aip1ip1(2,1) = -ad21/det
                     aip1ip1(2,2) = ad11/det
                     aii = aii - 
     &                    matmul(matmul(aiip1,aip1ip1),aip1i)
                     bi = bi - 
     &                    matmul(matmul(aiip1,aip1ip1),bip1)
                     jj = 0
                     do l = 1,idof
                        do k = 1,idof
                           jj=jj+1
                           a(ii+na(jj)) = aii(l,k) 
                        enddo
                        bp(ip2+nrhs(l)) = bi(l)
                     enddo
                  enddo
               endif
            enddo
            deallocate(aii,aiip1,aiim1)
            deallocate(aip1i,aip1ip1)
            deallocate(bi,bip1,bim1)
         endif    
      else if(iflg.eq.2) then
c     back out full solution
         if(idof.eq.1) then
            do i=1,neq_primary
               if(ngdpm_layers(igdpm(i)).gt.0) then
c     set up 1-d connectivity
c     ip1 = i-1, ip2 =i, ip3 = i+1
                  ip2 = i
                  ip2= nelm_gen(nelm_gen(ip2+1))
                  iiconn(1)=ip2
                  iim1conn(1) = i
                  do j=1,ngdpm_layers(igdpm(i))-1
                     ip2= nelm_gen(nelm_gen(ip2+1))
                     ip1= nelm_gen(nelm_gen(ip2)+1)
                     iiconn(j+1)=ip2
                     iim1conn(j+1) = ip1
                  enddo
c     perform back substitution
c     b(i)  =  (b(i) - a(i-1,i-1)*b(i-1))/a(i,i-1)
                  do j = 1,ngdpm_layers(igdpm(i))
                     ip2 = iiconn(j)
                     ip1 = iim1conn(j)
                     ii = nelmdg_gen(ip2) - neqp1 
                     iim1 = nelmdg_gen(ip2) -1 - neqp1 
                     bp(ip2) = (bp(ip2) - bp(ip1)*a(iim1))/a(ii)
                  enddo
               endif
            enddo
         else if (idof.eq.2) then
            allocate(aii(idof,idof),aiip1(idof,idof),aiim1(idof,idof))
            allocate(aip1i(idof,idof),aip1ip1(idof,idof))
            allocate(bi(idof),bip1(idof),bim1(idof))
            do i=1,neq_primary
               if(ngdpm_layers(igdpm(i)).gt.0) then
c     set up 1-d connectivity
c     ip1 = i-1, ip2 =i, ip3 = i+1
                  ip2 = i
                  ip2= nelm_gen(nelm_gen(ip2+1))
                  iiconn(1)=ip2
                  iim1conn(1) = i
                  do j=1,ngdpm_layers(igdpm(i))-1
                     ip2= nelm_gen(nelm_gen(ip2+1))
                     ip1= nelm_gen(nelm_gen(ip2)+1)
                     iiconn(j+1)=ip2
                     iim1conn(j+1) = ip1
                  enddo
c     perform back substitution
c     b(i)  =  (b(i) - a(i,i-1)*b(i-1))/a(i,i)
                  do j = 1,ngdpm_layers(igdpm(i))
                     ip2 = iiconn(j)
                     ip1 = iim1conn(j)
                     ii = nelmdg_gen(ip2) - neqp1 
                     iim1 = nelmdg_gen(ip2) -1 - neqp1 
                     jj = 0
                     do l = 1,idof
                        do k = 1,idof
                           jj=jj+1
                           aii(l,k) = a(ii+na(jj))
                           aiim1(l,k) = a(iim1+na(jj))
                        enddo
                        bi(l) = bp(ip2+nrhs(l))
                        bim1(l) = bp(ip1+nrhs(l))
                     enddo
c     find inverse of aii
                     ad11 = aii(1,1)
                     ad12 = aii(1,2)
                     ad21 = aii(2,1)
                     ad22 = aii(2,2)
                     det = ad11*ad22-ad12*ad21
c     put inverse back into aii
                     aii(1,1) = ad22/det
                     aii(1,2) = -ad12/det
                     aii(2,1) = -ad21/det
                     aii(2,2) = ad11/det
                     bim1 = matmul(aiim1,bim1)
                     bi = matmul(aii,bi - bim1)
                     do l = 1,idof
                        bp(ip2+nrhs(l)) = bi(l)
                     enddo
                  enddo
               endif
            enddo
            deallocate(aii,aiip1,aiim1)
            deallocate(aip1i,aip1ip1)
            deallocate(bi,bip1,bim1)
         endif
      endif
      end
      subroutine simplify_gdpm(iflg,neq,neq_primary,nelm_gen 
     &     ,idum,ncon_primary,a,adum,a_primary,na
     &     ,na_primary,ngdpm_layers,igdpm,iarr,idof)
c     
c     subroutine to simplify connectivity in gdpm simulation
c     
      implicit none
      real*8 a(*),a_primary(*),adum(*)
      integer iflg, neq, neq_primary, neqp1, i, ii, ic, i1, i2
      integer l, idof
      integer nelm_gen(*), igdpm(*), idum(*)
      integer ngdpm_layers(0:*)
      integer ncon_primary(*), na(*), na_primary(*), iarr(*)
c     
      if(iflg.eq.1) then
      else if(iflg.eq.3) then
      else if(iflg.eq.4) then
c     create iarr array for a matrix (jacobian)
         ic = 0
         do i = 1,neq_primary
            if(ngdpm_layers(igdpm(i)).ge.1) then
               do ii= nelm_gen(i)+1, nelm_gen(i+1) - 1
                  ic=ic+1
                  iarr(ic) = ii - (neq+1)
               enddo
            else
               do ii= nelm_gen(i)+1, nelm_gen(i+1) 
                  ic=ic+1
                  iarr(ic) = ii - (neq+1)
               enddo
            endif
         enddo
c     fill in degrees of freedom
         do i = 2,idof*idof
            do ii = 1,ic
               iarr(ic*(i-1)+ii) = iarr(ii)+na(i) 
            enddo
         enddo
      endif
      return
      end
