      subroutine stone_rlp(rlc, drlc, ii, rl, drlw, drlg, drlk,
     &		rv, drvw, drvg, drvk, r, drw, drg, drk)

!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
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
!D1 PURPOSE
!D1
!D1 Calculate relative permeabilities and capillary pressures
!D1 and derivatives.
!D1
!***********************************************************************
!     rlc = krowc oil relative perm in presence of connate water only
!     rw  = krw   water relative perm
!     rl  = krow  oil/water relative perm
!     rg  = krg   gas relative perm
!     rv  = krog  oil/gas relative perm
!     r   = kro   oil relative perm

      use comcomp, only: rl_g
	use comrlp
	use comai, only : neq
	use comco2, only: rl_w
      implicit none

	integer ii
      real*8 rlc, drlc, rw, rg, rl, drlw, drlg, drlk
	real*8 rv, drvw, drvg, drvk, r, drw, drg, drk
  
c Oil relative perm    
      r = rlc * ((rl/rlc + rl_w(ii))*(rv/rlc + rl_g(ii)) - rl_w(ii) 
     &	-rl_g(ii))
	drw = rlc*((drlw/rlc + rl_w(ii+neq))*(rv/rlc + rl_g(ii)) +
     &	(rl/rlc + rl_w(ii))*(rv/rlc + rl_g(ii+neq)) -
     &	rl_w(ii+neq) - rl_g(ii+neq))
	drg = rlc*((drlg/rlc + rl_w(ii+2*neq))*(rv/rlc + rl_g(ii)) +
     &	(rl/rlc + rl_w(ii))*(rv/rlc + rl_g(ii+2*neq)) - 
     &	rl_w(ii+2*neq) - rl_g(ii+2*neq))
	drg = rlc*((drlg/rlc + rl_w(ii+2*neq))*(rv/rlc + rl_g(ii)) +
     &	(rl/rlc + rl_w(ii))*(rv/rlc + rl_g(ii+2*neq)) -
     &	rl_w(ii+2*neq) - rl_g(ii+2*neq))
	if(r.le.0.d0) then
		r = 0.d0
		drw = 0.d0
		drg = 0.d0
		drk = 0.d0
	endif      
      end subroutine stone_rlp
