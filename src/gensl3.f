      subroutine gensl3
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
CD1  solve the heat conduction equation.
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
CD2 $Log:   /pvcs.config/fehm90/src/gensl3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.10   Thu May 09 15:11:16 1996   gaz
CD2 added call to md_node
CD2 
CD2    Rev 1.9   Wed Jan 17 08:58:40 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.8   Tue Jan 16 14:31:14 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.7   Fri Jan 12 17:49:06 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.6   Wed Jan 10 12:08:16 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.5   Tue Jan 09 14:30:52 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.4   03/23/95 19:03:52   gaz
CD2 gaz put abs(bp instead of bp in 1 if block
CD2 
CD2    Rev 1.3   05/11/94 16:14:52   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:43:14   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:59:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:26   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.1 Heat-conduction equations
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
CC**********************************************************************

      use combi
      use comdi
      use comei
      use comgi
      use comci
      use davidi
      use comdti
      use comai
      implicit none
      
      integer i, i1, i2, id, idg, ikd, neqp1, nsizea
      real*8, allocatable :: anz(:)
      real*8, allocatable :: bn(:)
      real*8, allocatable :: sto5(:)
      real*8, allocatable :: sto6(:)
      real*8, allocatable :: sto7(:)
      real*8, allocatable :: sto8(:)
      real*8, allocatable :: sto9(:)
      real*8, allocatable :: dum(:)
      real*8  bp1, facr, fdm, fdum2, pive, tollr, tolls
      
      parameter(fdm=20.0)
      
      neqp1=neq+1
c     zero out arrays
      do i=1,neq
         bp(i)=0.0
      enddo
      nsizea=nelm(neqp1)-neqp1
      do i=1,nsizea
         a(i)=0.0
      enddo
      fdum2=0.
      do id=1,neq
         call geneq3(id)
      enddo
c
c simplify a matrix if activated
c
      call active_nodes_ctr(5)

c     call md_modes to complete equations for multiply defined nodes
c     
c don't need to call geneqmdnode (through md_nodes)
c     if(imdnode.ne.0) call md_nodes(3,0,0)
c     
c     normalize wrt diagonal
c     
      do id=1,neq
         i1=nelm(id)+1
         i2=nelm(id+1)
         idg=nelmdg(id)
         pive=1./a(idg-neqp1)
         bp(id)=bp(id)*pive
         do ikd=i1,i2
            a(ikd-neqp1)=a(ikd-neqp1)*pive
         enddo
         fdum2=fdum2+bp(id)*bp(id)
      enddo
      fdum=sqrt(fdum2)
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
         do i=1,neq
            if(abs(bp(i)).gt.tmch) go to 99
         enddo
         fdum=-1.0
         go to 999
 99      continue
         f0=-1.0
      endif
      if(fdum.le.f0.and.iad.ne.0) goto 999
      facr=1.0
      tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
c     tollr=g3*tolls
c     if ((tmch*0.01).gt.tollr) tollr=tmch*0.01
      tollr=g3*max(tolls,tmch)
      
      call storage_derivatives(1,1)
      iter=maxsolve     
      if(irdof.eq.0) then
         allocate(dum(neq*1*4))
         if (igauss .gt. 1) then
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
           else if (igauss .eq. -1) then
c in the active nodes call the equations are renumbered and the 
c jacobean matrix simplified (call active_nodes_ctr(5)) the connectivity
c for a and its incomplete factorization is in nop
            call solve_new(neq,a,b,bp,nb,nb,nrhs,nop,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
           else
            call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
         end if
         deallocate(dum)
      else
         allocate(anz(nbd/2),bn(nbd/2))
         allocate(dum(neq*1*4))
c zvd 12-Jan-12 modified call to rd1dof to use correct number of calling parameters -- used nmat, nb, and nrhs based on other routines and their calls to solve_new
         if (igauss .gt. 1) then
            call rd1dof(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north,mink,
     *           gmres,tollr,irb,
     *           iirb,nopt,npvt,dum,dum(neq+1),dum(neq*2+1),
     *           dum(neq*3+1),piv,iter,irdof,icoupl,overf,anz,bn,sto5,
     *           sto6,sto7,sto8,sto9,nbnd,iback,iout,iptty,
     *           maxor,h,c,ss,g,y,maxsolve,accm)
         else
            call rd1dof(neq,a,b,bp,nmat,nb,nrhs,nelm,nelm,north,mink,
     *           gmres,tollr,irb,
     *           iirb,nopt,npvt,dum,dum(neq+1),dum(neq*2+1),
     *           dum(neq*3+1),piv,iter,irdof,icoupl,overf,anz,bn,sto5,
     *           sto6,sto7,sto8,sto9,nbnd,iback,iout,iptty,
     *           maxor,h,c,ss,g,y,maxsolve,accm)
         end if
         deallocate(anz,bn)
         deallocate(dum)
      endif
      call storage_derivatives(0,1)
      
      itert=itert+iter
      itotals=itotals+iter
      minkt=minkt+mink
      do id=1,neq
         bp1=bp(id)
         bp(id)=0.
         bp(id+neq)=bp1      
      enddo

 999  continue
      
      return
      end
