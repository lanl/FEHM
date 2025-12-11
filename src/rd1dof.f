      subroutine rd1dof(neq,a,b,r,na,nb,nrhs,ncon,nop,north,mink,
     *sorthm,epn,irb,iirb,nopt,npvt,xtemp,rw,xx,dum,piv,iter,irdof,niter
     *,omega,an,bn,resid,delx,dum1,irbd,npvtr,nbnd,mcount,iout,iptty
     * ,maxor,h,c,ss,g,y,maxsolve,accm)
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
CD1  This subroutine reorders the nodes based on the residuals.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rd1dof.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Mon Apr 14 12:44:04 1997   gaz
CD2 minor correction
CD2 
CD2    Rev 1.2   Thu Jan 11 09:47:30 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:59:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:46   pvcs
CD2 original version in process of being certified
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
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      implicit none

      integer iout, iptty, irdof, iter, maxor, mcount 
      integer nbnd, neq, niter, north, mink
      integer it1(20),it2(20),npvtr(1),irbd(1)
      integer irb(*),iirb(*),nopt(*),npvt(*)
      integer i, i1, i2, i3, i4, iback, ic, icount, id, ig, ik , im
      integer in, inn, io, ioo, iterg, iterh, itr
      integer j, jd, jdd, je, kb, kbn, kbns
      integer na(*),nb(*),nrhs(*) 
      integer nact, nactp1, neqp1, nopc, nrest
      integer maxsolve
      integer ncon(*),nop(*)
      real*8  epn, omega
      real*8  a(*),b(*),r(*),sorthm(*)
      real*8  xtemp(*),dum(*),rw(*),piv(*),xx(*)
      real*8  h(*),c(*),ss(*),g(*),y(*)
      real*8  an(*),bn(*),resid(*),delx(*),dum1(*)
      real*8  dum1d
      real*8  fdum, fdum2, epn2
      real*8  residi, rsum, sum
      real*8  tolfac, tollrd 
      character*4 accm
c
      icount=0
      iterg=0
      iback=0
! Guessed at value for iterh
      iterh = 0
      epn2=epn*epn
c
c
c find residuals after SOR iterations
c delx is the solution vector
c resid is the residual of the matrix equation
c nact is the number of active varibles
c
c
      neqp1=neq+1
      if(irdof.lt.0) then
      tolfac=10**(-irdof)
      else
      tolfac=1./10**(irdof)
      endif
      fdum2=0.0
      iback=0
c
      do 49 i=1,neq
      delx(i)=0.0
      resid(i)=r(i)
      r(i)=0.0
      irbd(i)=i
      fdum2=fdum2+resid(i)*resid(i)
49    continue
      fdum=sqrt(fdum2)
c
c
c SOR loop
c
      do 1 icount=1,mcount
c
      do 51 i=1,neq
      delx(i)=0.0
51    continue
c
c
c reorder equations as necessary
c
      if(icount.eq.1.or.iback.eq.0) then
c
c find difficult equations (ie the ones that are gt fdum/neq*tolfac)
      tollrd=fdum/neq*tolfac
c
c reset reorder arrays
c
      do 350 i=1,neq
      irb(i)=0
      iirb(i)=0
350   continue
c
      nact=0
c
      do 400 i=1,neq
      residi=resid(i)
      if(residi.gt.tollrd) then
      if(irb(i).eq.0) then
      nact=nact+1
      irb(i)=nact
      iirb(nact)=i
      endif
      endif
c
400   continue
c
c now label the rest
c
      nrest=nact
      do 550 i=1,neq
      if(irb(i).eq.0) then
      nrest=nrest+1
      irb(i)=nrest
      iirb(nrest)=i
      endif
550   continue
c
c if all active,then don't reorder
      if(nact.eq.neq) then
c
c note changes to calling sequence(irbd)
      iter=maxsolve
       call solve_new(neq,a,b,resid,na,nb,nrhs,
     2ncon,nop,north,epn,irb,iirb,npvt,xtemp,
     3dum,piv,h,c,ss,g,y,iter,0,1,iptty,maxor,accm)
c      call solve_new(neq,a,b,resid,ncon,ncon,north,sorthm,tolld
c     *,irbd,irbd,nopt,npvt,xtemp,rw,xx,dum,piv,iter
c     *,nbnd,iback,iptty,maxor,h,c,ss,g,y,accm)
c
      do 560 i=1,neq
      r(i)=resid(i)
560   continue
c
      return
      endif
c if all passive,then go to SOR| iteration
      if(nact.eq.0) then
      go to 306
      endif
c
c
c coding to arrange the a matrix for the PCG solution of the hard part
c
c
      nactp1=nact+1
      nopc=nactp1
      nop(1)=nopc
c
c loop on active varibles
c
      do 101 in=1,nact
      io=iirb(in)
      i1=ncon(io)+1
      i2=ncon(io+1)
c
c find connections gt nact
c
      ic=0
      id=0
      do 110 jd=i1,i2
      kb=ncon(jd)
      kbn=irb(kb)
      if(kbn.gt.nact) then
c connections gt nact
      ic=ic+1
      it1(ic)=jd
      else
c connection le nact
      id=id+1
      it2(id)=jd
      endif
110   continue
c
c zero and load dummy array
c
      do 115 ik=1,neq
      sorthm(ik)=0.0
      sorthm(neq+ik)=0.0
115   continue
c
c dum will contain row of aa
      do 120 jd=1,id
      je=it2(jd)
      kb=ncon(je)
      kbn=irb(kb)
      sorthm(kbn)=a(je-neqp1)
120   continue
c
c load vector of row a(in,kbn)
c
      do 130 jd=1,ic
      jdd=it1(jd)
      kb=ncon(jdd)
      kbns=irb(kb)
      sorthm(neq+kbns)=a(jdd-neqp1)
130   continue
c
c add contribution from [ap][pa]
c
      do 140 inn=nact+1,neq
      ioo=iirb(inn)
      dum1d=sorthm(neq+inn)
      if(dum1d.ne.0.0) then
      i3=ncon(ioo)+1
      i4=ncon(ioo+1)
      do 150 ig=i3,i4
      kb=ncon(ig)
      kbn=irb(kb)
      sorthm(kbn)=sorthm(kbn)-dum1d*a(ig-neqp1)
150   continue
      endif
c
140   continue
c
c row of [a]-[ap][pa] is now complete
c
      do 160 im=1,nact
      if(sorthm(im).ne.0.0) then
      nopc=nopc+1
      nop(nopc)=im
      an(nopc-nactp1)=sorthm(im)
      endif
      if(im.eq.in) npvtr(in)=nopc
160   continue
c
      nop(in+1)=nopc
101   continue
      endif
c
c calulate residual for active varibles
c
      do 102 in=1,nact
      io=iirb(in)
      rsum=resid(io)
      i1=ncon(io)+1
      i2=ncon(io+1)
      do 170 jdd=i1,i2
      kb=ncon(jdd)
      kbn=irb(kb)
      if(kbn.gt.nact) rsum=rsum-a(jdd-neqp1)*resid(kb)
170   continue
      dum(in)=rsum
102   continue
c
c
      iter=maxsolve
c call routine solve with the reduced matrix
c
c
      iter=maxsolve
       call solve_new(neq,a,b,resid,na,nb,nrhs,
     2ncon,nop,north,epn,irb,iirb,npvt,xtemp,
     3dum,piv,h,c,ss,g,y,iter,0,1,iptty,maxor,accm)
c      call solve_new(nact,an,bn,dum,nop,nop,north,sorthm,epn
c     *,irbd,irbd,nopt,npvtr,xtemp,rw,xx,dum1,piv,iterh
c     *,nbnd,iback,iptty,maxor,h,c,ss,g,y,accm)
      iterg=iterg+iterh
c
c
c coding to extract solution
c
c
c
c first load delx for active varibles
c
      do 195 in=1,nact
      io=iirb(in)
      delx(io)=dum(in)
195   continue
c
      do 190 i=nact+1,neq
c find old node number
      io=iirb(i)
      i1=ncon(io)+1
      i2=ncon(io+1)
      sum=resid(io)
      do 180 j=i1,i2
      kb=ncon(j)
      sum=sum-a(j-neqp1)*delx(kb)
180   continue
c adjust diagonal and extract solution
      delx(io)=sum*omega+delx(io)
190   continue
c
      iback=0
c
c  SOR coding to check residuals
c
c     loop on neighbors of i
306   do  300 itr=1,niter
      do 100 i=1,neq
      i1=ncon(i)+1
      i2=ncon(i+1)
      sum=resid(i)
      do 200 j=i1,i2
      kb=ncon(j)
      sum=sum-a(j-neqp1)*delx(kb)
200   continue
      sorthm(i)=sum
c
c adjust diagonal term(rem aii=1.0!)
c
      delx(i)=sum*omega+delx(i)
c
c
100   continue
300   continue
c
c check if residuals are small enough(ie fdum2 le epn**2)
c
      fdum2=0.0
      do 304 i=1,neq
      fdum2=fdum2+sorthm(i)*sorthm(i)
      r(i)=r(i)+delx(i)
      resid(i)=sorthm(i)
304   continue
      if(fdum2.le.epn2) then
      iter=iterg
      return
      endif
      fdum=sqrt(fdum2)
1     continue
      if (iout .ne. 0) then
         write(iout,*) 'overall iterations exceeded in rd1dof'
         write(iout,*) 'fdum2 ',fdum2,' epn**2 ',epn2
      end if
      if (iptty .ne. 0) then
         write(iptty,*) 'overall iterations exceeded in rd1dof'
         write(iptty,*) 'fdum2 ',fdum2,' epn**2 ',epn2
      end if
      return
      end
