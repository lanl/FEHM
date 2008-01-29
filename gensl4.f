      subroutine gensl4
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
CD1  This subroutine sets up for and calls the necessary routines to
CD1  solve the heat and mass equations with non-condensible gas
CD1  (full jacobian, unsymmetric 3n by 3n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gensl4.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.20   Fri Apr 26 15:28:52 1996   gaz
CD2 changes for mdnodes and normalization( many)
CD2 
CD2    Rev 1.19   Wed Feb 07 11:15:30 1996   gaz
CD2 commented out storage for new normalize routine
CD2 
CD2    Rev 1.18   Wed Feb 07 11:13:14 1996   gaz
CD2 added(but commented out new normalize routine
CD2 
CD2    Rev 1.17   Wed Jan 17 08:59:24 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.16   Tue Jan 16 14:31:16 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.15   Fri Jan 12 17:49:12 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.14   Wed Jan 10 13:33:36 1996   hend
CD2 Ammended Requirements Traceability
CD2 
CD2    Rev 1.12   Wed Jan 10 12:30:20 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.11   Wed Jan 10 12:16:24 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.10   12/13/95 08:41:44   gaz
CD2 clean up :position of memory calls
CD2 
CD2    Rev 1.9   08/02/95 17:31:44   llt
CD2 moved where sto5 was getting released
CD2 
CD2    Rev 1.8   06/01/95 16:49:34   gaz
CD2 minor corrections
CD2 
CD2    Rev 1.7   04/25/95 09:06:36   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.6   03/23/95 19:17:28   gaz
CD2 gaz changed bp to abs(bp in one if block
CD2 
CD2    Rev 1.5   03/10/95 10:49:10   llt
CD2 removed unneeded pointers - gaz
CD2
CD2    Rev 1.4   05/11/94 16:14:54   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.3   03/24/94 11:42:38   llt
CD2 Fixed it so the arrays in storage_derivate won't be deallocated, 
CD2 unless they were allocated.
CD2 
CD2    Rev 1.2   03/23/94 14:43:14   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:59:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:28   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3  
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.3.3 Noncondensible gas flow equations
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use davidi
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none
      
      integer ndex(3)
      integer i, id, neqp1, nsizea, nsizea1
      real*8, allocatable :: sto5(:,:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: dumn(:)
      real*8  facr, fdm, fdum2, tollr, tolls
      parameter(fdm=20.0)

      neqp1=neq+1
c     zero out arrays
      do i=1,neq
         bp(i+nrhs(1))=0.0d00
         bp(i+nrhs(2))=0.0d00
         bp(i+nrhs(3))=0.0d00
      enddo
      nsizea1=nelm(neqp1)-neqp1
      nsizea=9*nsizea1
      do i=1,nsizea
         a(i)=0.0d00
      enddo
      fdum2=0.d00
      do id=1,neq
c     
c     decide on equation type
c     
         if(ps(id).gt.0.0) then
            call geneqc(id)
         else
            call geneqc(id)
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0d00
            a(nelmdg(id)-neqp1+nmat(9))=sx1(id)
            bp(id+nrhs(3))=0.0d00
         endif
      enddo

c     call md_modes to complete equations for multiply defined nodes
c     
c don't need to call geneqmdnode (through md_nodes)
c     if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     check for pci=0 and modify equations
c     
c     call mod_eqs_ngas(neq,nelm,nelmdg,a,nmat,bp,nrhs,
c    &     phi,pci,pcp,dpcef,t,ieos)

      call dual(10)
c     
c     normalize and obtain sum square of residuals
c     
      if(islord.ne.0) call switch(nmat,nmatb,islord,3,1,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,3,1,neq)
c
c
c     first call air_rdof to rearrange equations
c     
        if (irdof.ne.0) then
           call air_rdof(0,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
        endif
c     
      allocate(dumn(100))
      call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &     ,ndex,3,dumn(1),dumn(37),dumn(73),0,fdum2)
      deallocate(dumn)
c     
c     return and stop if singular
c     
      if(ndex(1).lt.0 ) then
         if (iout .ne. 0) then
            write(iout,*) '   '
            write(iout, 100)
            write(iout,*) '   '
            write(iout,'(i8,3g12.3)')
     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         end if
         if(iptty.gt.0)then
            write(iptty,*) '   '
            write(iptty, 100) 
            write(iptty,*) '   '
            write(iptty,'(i8,3g12.3)') 
     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         endif
         iad=maxit
         return
      endif
 100  format ('* singular matrix found during normalization *')
      
      fdum=sqrt(fdum2)
      if(fdum.eq.0.0d00) go to 999
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
         do i=1,neq
            if(abs(bp(i+nrhs(1))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(2))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(3))).gt.tmch) go to 99
         enddo
         fdum=-1.d000
         go to 999
 99      continue
         f0=-1.d000
      endif
      if(f0.gt.0.0d00) then
         if(fdum.le.f0.and.iad.ne.0) goto 999
         facr=1.0d00
         tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
         tollr=tolls*g3
         if(tmch*g3.gt.tollr) tollr=tmch*g3 
      else if(f0.le.0.0d00) then
         tollr=tmch *g3      
      endif

      iter=maxsolve
      call storage_derivatives(1,1)
      
      if(irdof.eq.0) then

c     full solution
         allocate(dum(neq*3*4))
         if (igauss .gt. 1) then
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,3,iptty,maxor,accm)
         else
            call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,3,iptty,maxor,accm)
         end if
         itert=itert+iter
         itotals=itotals+iter
         minkt=minkt+mink
         deallocate(dum)
      else if(irdof.eq.-3) then
c     reduced degree of freedom with full grmres
         allocate(dum(neq*3*4))
         if (igauss .gt. 1) then
            call solve_rdof(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,-3,iptty
     *           ,maxor,icoupl,tollr,overf,accm)
         else
            call solve_rdof(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,-3,iptty
     *           ,maxor,icoupl,tollr,overf,accm)
         end if
         
         itert=itert+iter
         itotals=itotals+iter
         minkt=minkt+mink
         deallocate(dum)
      else
         allocate(dum(neq*3*4))
         allocate(sto5(neq,3))
         if (igauss .gt. 1) then
            call rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,3,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm) 
         else
            call rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,3,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm) 
         endif
         itert=itert+iter
         itotals=itotals+iter
         minkt=minkt+mink
         deallocate(dum)
         deallocate(sto5)
      endif

      call storage_derivatives(0,1)
c
c     extract proper solution
c
      if(irdof.ne.0) then
         call air_rdof(1,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
      endif
c


      if(islord.ne.0) call switch(nmat,nmatb,islord,3,2,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,3,2,neq)

      call dual(11)
      
 999  continue
      
      return
      end
      subroutine a_block(iflg,a,bp,ncon,nelmdg,nmat,nrhs,neq,idof)
c
c debug routine to test jacobian array
c     
      implicit none
      integer iflg,neq,idof
      integer ncon(*),nelmdg(*),nmat(*),nrhs(*)
      real*8 a(*),bp(*)

      integer neqp1,i,nsizea

      neqp1=neq+1
      nsizea = (ncon(neqp1)-neqp1)*9
      if(idof.eq.3) then
         do i=nmat(2)+1,nsizea
            a(i)=0.0
         enddo
         do i=1,neq
            bp(i+nrhs(2))=0.0
            bp(i+nrhs(3))=0.0
            a(nelmdg(i)-neqp1+nmat(5))=1.0
            a(nelmdg(i)-neqp1+nmat(9))=1.0
         enddo
      endif
      return
      end
