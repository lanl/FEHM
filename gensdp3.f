      subroutine gensdp3
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
CD1  This subroutine solves the non-isothermal air-water equations with
CD1  full jacobian (unsymmetric, 3n by 3).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gensdp3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.14   Fri May 31 15:08:34 1996   gaz
CD2 correction to call to rdof_new
CD2 
CD2    Rev 1.13   Wed May 29 08:21:06 1996   gaz
CD2 corrected call to switch and setting of nmatb
CD2 
CD2    Rev 1.12   Tue May 21 13:52:24 1996   gaz
CD2 added call md_nodes
CD2 
CD2    Rev 1.11   Wed Feb 07 11:09:18 1996   gaz
CD2 corrections to make dpdp air/water/heat work better
CD2 
CD2    Rev 1.10   Tue Jan 16 14:31:06 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.9   Fri Jan 12 17:48:46 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.8   Wed Jan 10 11:49:50 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.7   Wed Jan 10 10:32:32 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.6   Tue Jan 09 14:30:08 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.5   03/23/95 23:59:52   gaz
CD2 gaz replaced bp with abs(bp
CD2 
CD2    Rev 1.4   11/28/94 14:21:58   llt
CD2 Took out extra comma in argurment list, so would run on ibm.
CD2 
CD2    Rev 1.3   06/20/94 11:07:32   zvd
CD2  
CD2 
CD2    Rev 1.2   03/23/94 14:41:08   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:20   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
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

      use comai
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comhi
      use davidi
      use comdti
      implicit none

      integer nmatd, ndex(6)
      integer i, icoupl_iter, id, idl 
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls
      
c     call thermo for next level (i+neq)
c     
      call dpdp(2)
c     
c     adjust storage
c     
      neqp1=neq+1
      nmatd=nelm(neqp1)-neqp1
      nsizea=36*nmatd
      do i=1,nsizea
         a(i)=0.0
      enddo

      do i=1,neq
         bp(i+nrhs(1))=0.0
         bp(i+nrhs(2))=0.0
         bp(i+nrhs(3))=0.0
         bp(i+nrhs(4))=0.0
         bp(i+nrhs(5))=0.0
         bp(i+nrhs(6))=0.0
      enddo
      nrhs1=nrhs(1)
      nrhs2=nrhs(2)
      nrhs3=nrhs(3)

      do id=1,neq
         nmat(1)=0         
         nmat(2)=nmatd  
         nmat(3)=nmatd*2  
         nmat(4)=nmatd*(7-1)
         nmat(5)=nmatd*(8-1)
         nmat(6)=nmatd*(9-1)
         nmat(7)=nmatd*(13-1)
         nmat(8)=nmatd*(14-1)
         nmat(9)=nmatd*(15-1)
         nrhs(1)=nrhs1
         nrhs(2)=nrhs2
         nrhs(3)=nrhs3
c     
c     decide on equation type
c     
         if(ico2.gt.0) then
            if(ieos(id).ne.4) then
               call geneqc( id)
            else
               call geneq3( id)
               a(nelmdg(id)+nmat(1))=sx1(id)
               a(nelmdg(id)+nmat(9))=sx1(idl)
               bp(id+nrhs(1))=0.0
               bp(id+nrhs(3))=0.0
            endif
         else if (ico2.eq.0) then
            if(ieos(id).ne.4) then
               call geneq1( id)
            else
               call geneq3( id)
               a(nelmdg(id)+nmat(1))=sx1(id)
               a(nelmdg(id)+nmat(9))=sx1(idl)
               bp(id+nrhs(1))=0.0
               bp(id+nrhs(3))=0.0
            endif
         endif
         idl=id+neq

         nmat(1)=nmatd*(22-1)
         nmat(2)=nmatd*(23-1)
         nmat(3)=nmatd*(24-1)
         nmat(4)=nmatd*(28-1)
         nmat(5)=nmatd*(29-1)
         nmat(6)=nmatd*(30-1)
         nmat(7)=nmatd*(34-1)
         nmat(8)=nmatd*(35-1)
         nmat(9)=nmatd*(36-1)
         nrhs(1)=nrhs(4)
         nrhs(2)=nrhs(5)
         nrhs(3)=nrhs(6)
c     
c     decide on equation type(for node idl=id+neq)
c     
         if(ico2.gt.0) then
            if(ieos(idl).ne.4) then
               call geneqc(idl)
            else
               call geneq3(idl)
               a(nelmdg(id)+nmat(1))=sx1(id)
               a(nelmdg(id)+nmat(9))=sx1(idl)
               bp(id+nrhs(1))=0.0
               bp(id+nrhs(3))=0.0
            endif
         else if (ico2.eq.0) then
            if(ieos(idl).ne.4) then
               call geneq1(idl)
            else
               call geneq3(idl)
               a(nelmdg(id)+nmat(1))=sx1(idl)
               a(nelmdg(id)+nmat(9))=sx1(idl)
               bp(id+nrhs(1))=0.0
               bp(id+nrhs(3))=0.0
            endif
         endif
      enddo
      do i=1,36
         nmat(i)=nmatd*(i-1)
      enddo

      nrhs(1)=nrhs1
      nrhs(2)=nrhs2
      nrhs(3)=nrhs3

c     call md_modes to complete equations for multiply defined nodes
c     
c      if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     now get coupling for influence terms
c     
      call dpdp3

      if(islord.ne.0) call switch(nmat,nmatb,islord,6,1,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,6,1,neq)

c     normalize the equations
c     
c     call mmgetblk ("dumn", "dpdp", ipdumn, 108  , 2, icode)	
      allocate(dumn(108))
      
      call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     $     ,ndex,6,dumn(1),dumn(37),dumn(73),0,fdum2)
      
c     call mmrelblk ("dumn", "dpdp", ipdumn, icode)	
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
      
c     find residual
c     
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
            if(abs(bp(i+nrhs(5))).gt.tmch) go to 99
            if(abs(bp(i+nrhs(6))).gt.tmch) go to 99
         enddo
         fdum=-1.0
         go to 999
 99      continue
         f0=-1.0      
      endif
      if(f0.gt.0.0) then
         if(fdum.le.f0.and.iad.ne.0) goto 999
         facr=1.0
         tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
         tollr=tolls*g3
         if(g3*tmch.gt.tollr) tollr=g3*tmch
       else if(f0.le.0.0) then
         tollr=tmch *g3  
       endif

c
c     set maximum iterations in solve_new
c     
      iter=maxsolve
      call storage_derivatives(1,1)
      allocate(dumz(neq*5*6))
c
c     set maximum recoupling iterations if necesary
c     
      if(icoupl.gt.0) then
        icoupl_iter = icoupl
      else if(icoupl.lt.0) then
        icoupl_iter =
     &   max(float(abs(icoupl)),
     &   (float(iad)/float(abs(maxit)))*float(abs(2*icoupl)))
      endif
      if(irdof.eq.0) then
c     
c     full solution
c
         if (igauss .gt. 1) then
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north,tollr
     *           ,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,6,iptty,maxor,accm)
         else
            call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,6,iptty,maxor,accm)
         end if
      else if(irdof.eq.-6) then
c     reduced degree of freedom with full gmres
         if (igauss .gt. 1) then
            call solve_rdof(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,-6,iptty
     *           ,maxor,icoupl_iter,tollr,overf,accm)
         else 
            call solve_rdof(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,-6,iptty
     *           ,maxor,icoupl_iter,tollr,overf,accm)
         end if
      else
c     
c     approximate solution
c     
         allocate(sto5(neq*6))
         if (igauss .gt. 1) then
            call rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,6,iptty,maxor
     *           ,overf,irdof,icoupl_iter,0,sto5,accm)
         else
            call rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dumz,piv
     *           ,h,c,ss,g,y,iter,iback,6,iptty,maxor
     *           ,overf,irdof,icoupl_iter,0,sto5,accm)
         end if
         deallocate(sto5)
      endif
      itert=itert+iter
      itotals=itotals+iter
      minkt=minkt+mink
      deallocate(dumz)
      call storage_derivatives(0,1)

 999  continue
      
      if(islord.ne.0) call switch(nmat,nmatb,islord,6,2,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,6,2,neq)

c     extract dpdp solution
c     
      call dpdp(4)

      return
      end
