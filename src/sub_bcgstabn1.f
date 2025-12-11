      subroutine sub_bcgstabn1(n,a,bb,b,ncon,nop,l,work,tol
     *     ,irb,iirb,npvt,rw,dum,x,xtemp,piv
     *     ,h,c,s,g,yy,iter,idof,iptty,maxor)
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
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/sub_bcgstabn1.f_a  $ 
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2
!***********************************************************************c
c      subroutine bistbl (l, n, x, b, mv, solve, tol,
c     $     mxmv, work, ldw, rwork, ldrw, iwork, info)
c
c subroutine bistbl v1.0 1995
c
c Copyright (c) 1995 by D.R. Fokkema.
c Permission to copy all or part of this work is granted,
c provided that the copies are not made or distributed
c for resale, and that the copyright notice and this
c notice are retained.
c
c THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE AUTHOR
c PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESSED OR IMPLIED,
c REGARDING THE WORK, INCLUDING WARRANTIES WITH RESPECT TO ITS
c MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE.
c
C Parameter list (by M.Botchev, October, 1997):
C
C l      == (input) INTEGER BiCGstab's dimension
C n      == (input) INTEGER size of the system to solve 
C x      == (input/output) DOUBLE PRECISION array dimension n
C           initial guess on input, solution on output
C b      == (input) DOUBLE PRECISION array dimension n
C           right-hand side (rhs) vector
C mv     == (input) EXTERNAL name of matrix vector subroutine
C           to deliver y:=A*x by CALL mv(n,x,y)
C solve  == (input) EXTERNAL name of subroutine to perform
C           a preconditioner solve by CALL solve(n,x). This call 
C           delivers solution to M*x=x, M is the left preconditioning 
C           matrix, x is rhs on input and solution on output.
C tol    == (input) DOUBLE PRECISION tolerance for the relative
C           residual norm stopping criterion: ||r_k|| <= tol*||r_0||
C mxmv   == (input/output) INTEGER.  On input: maximum number of matrix 
C           vector multiplications allowed to be done.  On input: 
C           actual number of matrix vector multiplications done
C work   == (workspace) DOUBLE PRECISION array dimension (n,3+2*(l+1))
C ldw    == (input) INTEGER leading dimension of work (??)
C rwork  == (workspace) DOUBLE PRECISION array dimension (l+1,3+2*(l+1))
C ldrw   == (input) INTEGER leading dimension of rwork (??)
C iwork  == (workspace) INTEGER array dimension (l+1)
C info   == (output) INTEGER.  Info = 0 in case of normal computations
C           and 
C           info = m < 0 - means paramater number m has an illegal value
C           info = 1 - means no convergence achieved (stopping criterion
C           is not fulfilled)
C           info = 2 - means breakdown of the algorithm (try to enlarge
C           parameter l to get rid of this)
C Three parameters were added to link the code to VAC:
C oktest == (input) LOGICAL.  Set to .true. if intermediate residuals should 
C           be printed
C nonzer == (input) LOGICAL, indicates whether the initial guess is zero
C           vector.  If (nonzer.eq..false.) one matrix-vector 
C           multiplication will be saved
C stc    == (input) CHARACTER*3, stopping criterion parameter:  
C           stc='abs' -- use absolute stopping criterion ||r_k|| <= tol,
C           otherwise use relative stopping criterion, as indicated above

      implicit none
c
c     .. Parameters ..
c
      integer l, n, ldw, ldrw, iwork(l+1), info
      integer :: mxmv = 0
      double precision x(*), b(*), tol
      double precision work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1))
c gaz 1-18-2002 ldw and ldrw not defined
      parameter (ldw=0,ldrw=0)
c
c     .. Matrix ..
c
c
c     .. Local ..
c
      logical rcmp, xpdt
      integer i, j, k, nmv
      double precision alpha, beta, omega, rho0, rho1, sigma
      double precision varrho, hatgamma
      double precision rnrm0, rnrm, rnrm2
      double precision mxnrmx, mxnrmr, kappa0, kappal
c
c     .. Work Aliases ..
c
      integer z, zz, y0, yl, y
      integer rr, r, u, xp, bp
c
c     .. Constants ..
c
      double precision  zero, one, delta
      parameter (zero = 0d0, one = 1d0, delta = 1d-2)
c
c     .. BLAS ..
c
c     subroutine daxpy
c     subroutine dcopy
c     subroutine dgemv
c     subroutine dsymv
c     function ddot
c     function dnrm2
c
c     .. LAPACK ..
c
c     subroutine dgetrf
c     subroutine dgetrs
c     subroutine dlacpy
cM    subroutine dlaset
c
      double precision dnrm2, ddot
c
c     .. Intrinsic ..
c
      intrinsic abs, max, sign, sqrt
c
      real*8 a(*),bb(*),rw(*),dum(*),xtemp(*)
      integer ncon(*),nop(*)
      real*8 piv(*)
      integer irb(*),iirb(*),npvt(*)
c gaz 1-18-2002 rearranged so would complile on the PC
      integer maxit, iter, idof, iptty, maxor
      real*8 h(maxor,*),c(*),s(*),g(*),yy(*),epn
c     integer maxit, iter, idof, iptty, maxor
      integer nrhs(1),nrhs_dum(1)

      maxit=iter
      iter=0

      nrhs(1)=0
      nrhs_dum(1)=0

c     ===========================
c     .. Executable Statements ..
c     ===========================
c
      info = 0

      if (l.lt.1) info = -1
      if (n.lt.1) info = -2
      if (tol.le.zero) info = -7
      if (mxmv.lt.0) info = -8

      rr = 1
      r = rr+1
      u = r+(l+1)
      xp = u+(l+1)
      bp = xp+1
      if (bp*n.gt.ldw) info = -10

      z = 1
      zz = z+(l+1)
      y0 = zz+(l+1)
      yl = y0+1
      y = yl+1
      if (y*(l+1).gt.ldrw) info = -12

c     if (info.ne.0) return
c
c     --- Initialize first residual
c
cM      do i=1,n
cM	 work(i,r) = b(i)
cM      enddo
      call equal_array(work(1,r),b,n,idof,nrhs_dum,nrhs)
      call constant_value(x,0.0d00,n,idof,nrhs)
      call solve(n,bb,nop,irb,iirb,npvt,piv,nrhs_dum,nrhs,idof,
     *           work(1,r),rw)
      call equal_array(dum,work(1,r),n,idof,nrhs,nrhs_dum)
cM      do i=1,n
cM         dum(i)=work(i,r+1)
cM      enddo
c
c     --- Initialize iteration loop
c
      nmv = 0

	call equal_array(work(1,rr),work(1,r),n,idof,nrhs_dum,nrhs_dum)
	call equal_array(work(1,bp),work(1,r),n,idof,nrhs_dum,nrhs_dum)
	call equal_array(work(1,xp),x,n,idof,nrhs_dum,nrhs)
	call constant_value(x,zero,n,idof,nrhs)
	call residual(n,work(1,r),1,rnrm0,idof,nrhs_dum)
cM      call dcopy (n, work(1,r), 1, work(1,rr), 1)
cM      call dcopy (n, work(1,r), 1, work(1,bp), 1)
cM      call dcopy (n, x, 1, work(1,xp), 1)
cM      call dlaset ('n', n, 1, zero, zero, x, 1)
cM      rnrm0 = dnrm2 (n, work(1,r), 1)
      rnrm = rnrm0

      mxnrmx = rnrm0
      mxnrmr = rnrm0
      rcmp = .false.
      xpdt = .false.

      alpha = zero
      omega = one
      sigma = one
      rho0 = one
c
c     --- Iterate
c
      do while (rnrm.gt.tol .and. iter.lt.maxit)
	iter=iter+1

c
c     =====================
c     --- The BiCG part ---
c     =====================
c
	 rho0 = -omega*rho0
	 do k=1,l
	    rho1 = ddot (n, work(1,rr), 1, work(1,r+k-1), 1)
	    if (rho0.eq.zero) then
	       info = 2
		goto 5000
	    endif
	    beta = alpha*(rho1/rho0)
	    rho0 = rho1
	    do j=0,k-1
	       do i=1,n
		  work(i,u+j) = work(i,r+j) - beta*work(i,u+j)
	       enddo
	    enddo
	    call mv (n, a, ncon, work(1,u+k-1), work(1,u+k))
	    call solve(n,bb,nop,irb,iirb,npvt,piv,nrhs_dum,nrhs,idof,
     *                 work(1,u+k),rw)
	    nmv = nmv+1
	    sigma = ddot (n, work(1,rr), 1, work(1,u+k), 1)
	    if (sigma.eq.zero) then
	       info = 2
		goto 5000
	    endif
	    alpha = rho1/sigma
	    call daxpy (n, alpha, work(1,u), 1, x, 1)
	    do j=0,k-1
	       call daxpy (n, (-alpha), work(1,u+j+1), 1,
     $              work(1,r+j), 1)
	    enddo
	    call mv (n, a, ncon, work(1,r+k-1), work(1,r+k))
	    call solve(n,bb,nop,irb,iirb,npvt,piv,nrhs_dum,nrhs,idof,
     *                 work(1,r+k),rw)
	    nmv = nmv+1
cM	    rnrm = dnrm2 (n, work(1,r), 1)
		call residual(n,work(1,r),1,rnrm,idof,nrhs_dum)
	    mxnrmx = max (mxnrmx, rnrm)
	    mxnrmr = max (mxnrmr, rnrm)
	 enddo

*         print *,'2: work(1,r) = ',work(1,r),' ',work(2,r),' ',work(3,r)
*         print *,'3:         x = ',x(1),' ',x(2),' ',x(3)
*         pause

c
c     ==================================
c     --- The convex polynomial part ---
c     ==================================
c
c        --- Z = R'R
c
	 do i=1,l+1
	    call dgemv ('t', n, l+1-(i-1), one, work(1,r+i-1),
     $           n, work(1,r+i-1), 1, zero, rwork(i,z+i-1), 1)
	    if(l-(i-1).ne.0) then
                 call dcopy (l-(i-1), rwork(i+1,z+i-1), 1,
     $           rwork(i,z+i), l+1)
            endif
	 enddo
*         do i=1,l+1
*            write(*,'(20f9.2)') (rwork(i,j),j=1,3+2*(l+1))
*         enddo

	 call dlacpy ('a', l+1, l+1, rwork(1,z), l+1,
     $        rwork(1,zz), l+1)
	 call dgetrf (l-1, l-1, rwork(2,zz+1), l+1,
     $        iwork, info)
c
c        --- tilde r0 and tilde rl (small vectors)
c
	 rwork(1,y0) = -one
	 call dcopy (l-1, rwork(2,z), 1, rwork(2,y0), 1)
	 call dgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork,
     $        rwork(2,y0), l+1, info)
	 rwork(l+1,y0) = zero

	 rwork(1,yl) = zero
	 call dcopy (l-1, rwork(2,z+l), 1, rwork(2,yl), 1)
	 call dgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork,
     $        rwork(2,yl), l+1, info)
	 rwork(l+1,yl) = -one
c
c        --- Convex combination
c
	 call dsymv ('u', l+1, one, rwork(1,z), l+1,
     $        rwork(1,y0), 1, zero, rwork(1,y), 1)
         kappa0 = sqrt(ddot (l+1, rwork(1,y0), 1,
     $        rwork(1,y), 1))

	 call dsymv ('u', l+1, one, rwork(1,z), l+1,
     $        rwork(1,yl), 1, zero, rwork(1,y), 1)
	 kappal = sqrt(ddot (l+1, rwork(1,yl), 1,
     $        rwork(1,y), 1))

	 call dsymv ('u', l+1, one, rwork(1,z), l+1,
     $        rwork(1,y0), 1, zero, rwork(1,y), 1)
	 varrho = ddot (l+1, rwork(1,yl), 1, rwork(1,y), 1)
     $            / (kappa0*kappal)
         hatgamma =
     $        sign(1d0,varrho)*max(abs(varrho),7d-1)
     $        * (kappa0/kappal)
	 call daxpy (l+1, (-hatgamma), rwork(1,yl), 1,
     $        rwork(1,y0), 1)
c
c        --- Update
c
	 omega = rwork(l+1,y0)

	 call dgemv ('n', n, l, (-one), work(1,u+1), n,
     $        rwork(2,y0), 1, one, work(1,u), 1)
	 call dgemv ('n', n, l, one, work(1,r), n,
     $        rwork(2,y0), 1, one, x, 1)
	 call dgemv ('n', n, l, (-one), work(1,r+1), n,
     $        rwork(2,y0), 1, one, work(1,r), 1)

	 call dsymv ('u', l+1, one, rwork(1,z), l+1,
     $        rwork(1,y0), 1, zero, rwork(1,y), 1)
	 rnrm2 = ddot (l+1, rwork(1,y0), 1,
     $        rwork(1,y), 1)
         if(rnrm2.lt.zero) then 
		rnrm = zero
         else 
		rnrm = sqrt (rnrm2)
         endif
c
c     ================================
c     --- The reliable update part ---
c     ================================
c
	 mxnrmx = max (mxnrmx, rnrm)
	 mxnrmr = max (mxnrmr, rnrm)
	 xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
	 rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr)
     $        .or.xpdt)
	 if (rcmp) then
	    call mv (n, a, ncon, x, work(1,r))
	    call solve(n,bb,nop,irb,iirb,npvt,piv,nrhs_dum,nrhs,idof,
     *                 work(1,r),rw)
            nmv = nmv + 1
	    do i=1,n
	       work(i,r) =  work(i,bp) - work(i,r)
	    enddo
	    mxnrmr = rnrm
	    if (xpdt) then
	       call daxpy (n, one, x, 1, work(1,xp), 1)
cM	       call dlaset ('n', n, 1, zero, zero, x, 1)
	call constant_value(x,zero,n,idof,nrhs)
cM	       call dcopy (n, work(1,r), 1, work(1,bp), 1)
	call equal_array(work(1,bp),work(1,r),n,idof,nrhs_dum,nrhs_dum)
	       mxnrmx = rnrm
	    endif
	 endif
c         print *,'rnrm =',rnrm
      enddo
c
c     =========================
c     --- End of iterations ---
c     =========================
c
cM      call daxpy (n, one, work(1,xp), 1, x, 1)
c
c     --- Check stopping criterion
c
cM      call mv (n, a, ncon, x, work(1,r))
cM      do i=1,n
cM	 work(i,r) = b(i) - work(i,r)
cM      enddo
cM      call solve(n,bb,nop,irb,iirb,npvt,piv,nrhs_dum,nrhs,idof,
cM     *  work(1,r),rw)
cM      rnrm = dnrm2 (n, work(1,r), 1)
cM	call residual(n,work(1,r),1,rnrm,idof,nrhs_dum)
      if (rnrm.gt.tol) info = 1
c
c     --- Return
c
c     tol = rnrm/rnrm0
      mxmv = nmv

 5000 call equal_array(b,x,n,idof,nrhs,nrhs)
      call equal_array(rw,work(1,r),n,idof,nrhs,nrhs)

      if( info.eq.1 ) then
            if(iptty .ne. 0 ) then
               write(iptty,34)
 34            format(/,1x,'Warning issued by 1 degree of freedom '
     &              ,'solver:')
               if (iter .gt. maxit) write(iptty,35) maxit
 35            format(2x,'Maximum number of iterations (',i4,') '
     &              ,'exceeded.')
               write(iptty,36) rnrm,tol
 36            format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &              e14.6,/)
            endif
      endif

      return
      end
