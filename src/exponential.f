      subroutine exponential(s,p1,p2,p3,p4,r,dr,rn,drn)
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
!D1 Exponential relative permeability and capillary pressure model.
!D1
!***********************************************************************
c     input:  s saturation of wetting phase
c     p1 residual sat of wetting phase
c     p2 max sat of wetting phase
c     p3 = exponent
c     p4 = multiplier
c     return:  r (wetting phase), rn (non-wetting phase), and derivatives 
c     with respect to wetting phase 
c gaz 112623 cosmetic only       
      implicit none
      real*8 s,p1,p2,p3,r,dr,star,p4,rn,drn
      r=0;dr=0
      if(s.ge.p1.and.s.le.p2) then
         star=(s-p1)/(p2-p1)
         r= p4*star**p3
         dr = p4*p3*star**(p3-1.0)/(p2-p1)
	   rn=p4*(1.-star)**p3
         drn= -1.*p4*p3*(1.-star)**(p3-1.0)/(p2-p1)
      elseif (s.gt.p2) then
         r=1.
         dr=0
         rn=0.
         drn=0.
      elseif (s.lt.p1) then
         r=0
         dr=0
         rn=1.
         drn=0.
      endif
      
      end subroutine exponential

