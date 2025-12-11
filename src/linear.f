	subroutine linear(s,r1,r2,r,dr)
!***********************************************************************
! Copyright 2009 Los Alamos National Security, LLC  All rights reserved
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
!D1 Linear relative permeability and capillary pressure model.
!D1
!***********************************************************************

	implicit none

	real*8 s,r1,r2,r,dr,d
! r1 is smin; r2 is smax
	r=0.;dr=0.
	if(s.ge.r1) then
	   r=1.0
	   dr=0.
	   if(s.le.r2) then
	      d=r2-r1
	      r=(s-r1)/d
	      dr=1.0/d
	   endif
	endif

	end subroutine linear
