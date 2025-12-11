       subroutine gensdp_switch
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!D1
!D1  PURPOSE
!D1
!D1  This subroutine solves the isothermal air-water equations with 
!D1  full jacobian (unsymmetric, 2n by 2n). Richards Equation
!D1
!***********************************************************************
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: Oct-2008, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/gensdp_switch.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.3.2 Heat- and mass-transfer equations
!D3  2.3.3 Noncondensible gas flow equations
!D3  2.5.2 Solve nonlinear equation set at each time step
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

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

      integer nc(16),inorm
      integer nmatd,ndex(6)
      integer i, id, idl, neqp1, nrhs1, nrhs2, nsizea
      integer j, i1, i2 
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto1(:)
      real*8, allocatable :: sto5(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: piv_gaz(:,:)
      real*8  facr, fdm, fdum2, tollr, tolls
      real*8 a_total
      parameter(fdm=20.0, inorm = 1)
c     
c     call thermo for next level (i+neq)
c     
      call dpdp(2)
c     
c     adjust storage
c     
      nstorepiv = neq*4
      neqp1=neq+1
      fdum2=0.
      nmatd=nelm(neqp1)-neqp1
      nsizea=8*nmatd
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
      if(.not.allocated(piv_gaz)) then
         allocate(piv_gaz(neq,4))
      endif
      do id=1,neq
         nmat(1)=0         
         nmat(2)=nmatd    
c         nmat(3)=nmatd*(5-1)
c         nmat(4)=nmatd*(6-1)
         nrhs(1)=nrhs1
         nrhs(2)=nrhs2
c     
c     decide on equation type
c     
         if(ico2.lt.0) then
            if(ieos(id).ne.4) then
               call geneq2_switch( id)
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
         nmat(1)=nmatd*(7-1)
         nmat(2)=nmatd*(8-1)
c         nmat(3)=nmatd*(15-1)
c         nmat(4)=nmatd*(16-1)
         nrhs(1)=nrhs(3)
         nrhs(2)=nrhs(4)
         if(ico2.lt.0) then
            if(ieos(idl).ne.4) then
               call geneq2_switch(idl)
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
         continue
      enddo
      do i=1,8
         nmat(i)=nmatd*(i-1)
      enddo
      nrhs(1)=nrhs1
      nrhs(2)=nrhs2
c     
c     call md_modes to complete equations for multiply defined nodes
c     
c     if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     now get coupling for influence terms
c     
      if(ico2.eq.0) then
         call dpdpfh
      else if(ico2.lt.0) then
         call dpdpfa
      endif

c     
c     arrange arrays for 2 dof solve (ie variable switching)
c     
      call air_combine(1)
c     
c     normalize the equations
c     
c     some changes to call sequence here
c     
      if(islord.ne.0) then
         call switch(nmat,nmatb,islord,2,1,nmatd)
         call switchb(bp,nrhs,nrhsb,islord,2,1,neq)
      end if

      ndex(1) = 1
      if(inorm.ne.0) then      
         allocate(dumn(36))
c     call nrmlz4(neq,a,nmat,bp,nrhs,nelm,fdum2,sto1)
         call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &        ,ndex,2,dumn(1),dumn(17),dumn(33),0,fdum2)
         deallocate(dumn)
      else
	 fdum2 = 0.0
	 do i = 1,neq
            fdum2 = fdum2 + bp(i+nrhs(1))*bp(i+nrhs(1)) +
     &           bp(i+nrhs(2))*bp(i+nrhs(2))
	 enddo
      endif
      
c     
c     return and stop if singular
c     
      if(ndex(1).lt.0 ) then
         if (iout .ne. 0) then
            write(iout,*) '   '
            write(iout, 100)
            write(iout,*) '   '
c            write(iout,'(i8,3g12.3)')
c     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
            write(iout,'(i8,3g12.3)')
     &           abs(ndex(1)),(cord(abs(ndex(1)),i),i=1,3)
         end if
         if(iptty.gt.0) then
            write(iptty,*) '   '
            write(iptty,100) 
            write(iptty,*) '   '
c            write(iptty,'(i8,3g12.3)') 
c     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
            write(iptty,'(i8,3g12.3)') 
     &           abs(ndex(1)),(cord(abs(ndex(1)),i),i=1,3)
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
c     
c     set maximum iterations in solve_new
c     

      iter=maxsolve 
      call storage_derivatives(1,1)
c     allocate(dum(neq*4*4))
      allocate(dum(neq*8))      
      if(irdof.eq.0) then
c     
c     full solution
c     
         
         if (igauss .gt. 1) then
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv_gaz
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
         else
            call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv_gaz
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
         end if
      else
         allocate(sto5(neq*4))
         if (igauss .gt. 1) then
            call  rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv_gaz
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm)
         else
            call  rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv_gaz
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor
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
         call switch(nmat,nmatb,islord,2,2,nmatd)     
         call switchb(bp,nrhs,nrhsb,islord,2,2,neq)
      end if

c     sort out variable updates
      call air_combine(2)

c     extract dpdp solution
      call dpdp(4)
c     
c     
c     times=times-second(0.0)
c     write(59,*) 'time to solve jacobian(full)', times
      
      return
      end
