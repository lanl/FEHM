      subroutine gensdp
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
CD1  This subroutine solves the isothermal air-water equations with 
CD1  full jacobian (unsymmetric, 2n by 2n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gensdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.16   Mon Jul 01 14:17:22 1996   hend
CD2 Changed to use normal_dof
CD2 
CD2    Rev 1.15   Fri May 31 15:10:44 1996   gaz
CD2 cleaned up some comments
CD2 
CD2    Rev 1.14   Tue May 21 13:50:52 1996   gaz
CD2 added call to md_nodes
CD2 
CD2    Rev 1.13   Thu Feb 15 10:48:24 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.12   Wed Jan 17 08:56:42 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.11   Tue Jan 16 14:31:04 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.10   Fri Jan 12 17:48:40 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.9   Wed Jan 10 12:27:52 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.8   Wed Jan 10 11:47:44 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.7   Wed Jan 10 10:28:18 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.6   11/15/95 10:29:58   gaz
CD2 changes to get rdof working properly
CD2 
CD2    Rev 1.5   03/23/95 23:55:40   gaz
CD2 gaz put abs(bp instead of bp ...
CD2 
CD2    Rev 1.4   05/11/94 16:07:10   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.3   03/24/94 12:53:12   llt
CD2 Changed mmgetblk for array b at bottom to mmrelblk. Moved 
CD2 deallocation of array sto1 to after continue, so can't be missed.
CD2 
CD2    Rev 1.2   03/23/94 14:41:06   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:18   pvcs
CD2 original version in process of being certified
CD2 
CD2    c 17-mar-94 gaz
CD2    c getting and releasing b matrix
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

      integer nc(16)
      integer nmatd,ndex(6)
      integer i, id, idl, neqp1, nrhs1, nrhs2, nsizea 
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto1(:)
      real*8, allocatable :: sto5(:)
      real*8, allocatable :: dum(:)
      real*8  facr, fdm, fdum2, tollr, tolls
      parameter(fdm=20.0)

c call thermo for next level (i+neq)
c
      call dpdp(2)
c     
c     adjust storage
c     
      neqp1=neq+1
      fdum2=0.
      nmatd=nelm(neqp1)-neqp1
      nsizea=16*nmatd
      do i=1,nsizea
         a(i)=0.0
      enddo
      do i=1,neq
         bp(i+nrhs(1))=0.0
         bp(i+nrhs(2))=0.0
         bp(i+nrhs(3))=0.0
         bp(i+nrhs(4))=0.0
      enddo
      nrhs1=nrhs(1)
      nrhs2=nrhs(2)

      do id=1,neq
         nmat(1)=0         
         nmat(2)=nmatd    
         nmat(3)=nmatd*(5-1)
         nmat(4)=nmatd*(6-1)
         nrhs(1)=nrhs1
         nrhs(2)=nrhs2
c     
c     decide on equation type
c     
         if(ico2.lt.0) then
            if(ieos(id).ne.4) then
               call geneq2( id)
               call add_accumulation(id)
            else
               call geneq3( id)
               a(nelmdg(id)+nmat(1))=sx1(id)
               bp(id+nrhs(1))=0.0
            endif
         else if (ico2.eq.0) then
            if(ieos(id).ne.4) then
               call geneq1( id)
            else
               call geneq3( id)
               a(nelmdg(id)+nmat(1))=sx1(id)
               bp(id+nrhs(1))=0.0
            endif
         endif
c     
c     decide on equation type(for node idl=id+neq)
c     
         idl=id+neq
         nmat(1)=nmatd*(11-1)
         nmat(2)=nmatd*(12-1)
         nmat(3)=nmatd*(15-1)
         nmat(4)=nmatd*(16-1)
         nrhs(1)=nrhs(3)
         nrhs(2)=nrhs(4)
         if(ico2.lt.0) then
            if(ieos(idl).ne.4) then
               call geneq2(idl)
               call add_accumulation(idl)
            else
               call geneq3(idl)
               a(nelmdg(id)+nmat(1))=sx1(id)
               bp(id+nrhs(1))=0.0
            endif
         else if (ico2.eq.0) then
            if(ieos(idl).ne.4) then
               call geneq1(idl)
            else
               call geneq3(idl)
               a(nelmdg(id)+nmat(1))=sx1(id)
               bp(id+nrhs(1))=0.0
            endif
         endif
      enddo
      do i=1,16
         nmat(i)=nmatd*(i-1)
      enddo
      nrhs(1)=nrhs1
      nrhs(2)=nrhs2
c     
c     call md_modes to complete equations for multiply defined nodes
c     
c      if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     now get coupling for influence terms
c     
      if(ico2.eq.0) then
         call dpdpfh
      else if(ico2.lt.0) then
         call dpdpfa
      endif

c     normalize the equations
c     
c     some changes to call sequence here
c     
      if(islord.ne.0) then
         call switch(nmat,nmatb,islord,4,1,nmatd)
         call switchb(bp,nrhs,nrhsb,islord,4,1,neq)
      end if
      
      allocate(dumn(36))
c     call nrmlz4(neq,a,nmat,bp,nrhs,nelm,fdum2,sto1)
      call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &    ,ndex,4,dumn(1),dumn(17),dumn(33),0,fdum2)

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
         if(iptty.gt.0) then
            write(iptty,*) '   '
            write(iptty,100) 
            write(iptty,*) '   '
            write(iptty,'(i8,3g12.3)') 
     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         endif
         iad=maxit
         return
      endif
 100  format ('* singular matrix found during normalization *')

      fdum=sqrt(fdum2)
      if(fdum.eq.0.0) go to 999
      mink=n
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
         do i=1,neq
            if(abs(bp(i+nrhs(1))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(2))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(3))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(4))).gt.tmch) go to 99
         enddo
         fdum=-1.0
         go to 999
 99      continue
         f0=-1.0       
      endif
      if(fdum.le.f0.and.iad.ne.0) goto 999
      facr=1.0
      tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
      tollr=tolls*g3
      if(g3*tmch.gt.tollr) tollr=g3*tmch

c     set maximum iterations in solve_new
c     
      iter=maxsolve 
      call storage_derivatives(1,1)

      allocate(dum(neq*4*4))
      if(irdof.eq.0) then
c     
c     full solution
c  
         if (igauss .gt. 1) then
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,4,iptty,maxor,accm)
         else
            call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,4,iptty,maxor,accm)
         end if
      else
         allocate(sto5(neq*4))
         if (igauss .gt. 1) then
            call  rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,4,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm)
         else
            call  rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,4,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm)
         end if
         deallocate(sto5)
      endif
      
      deallocate(dum)
      call storage_derivatives(0,1)
      itert=itert+iter
      itotals=itotals+iter
      minkt=minkt+mink
 999  continue
      if(islord.ne.0) then
         call switch(nmat,nmatb,islord,4,2,nmatd)     
         call switchb(bp,nrhs,nrhsb,islord,4,2,neq)
      end if

c     extract dpdp solution
c     
      call dpdp(4)
c     
c     times=times-second(0.0)
c     write(59,*) 'time to solve jacobian(full)', times
      
      return
      end
