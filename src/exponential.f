      subroutine exponential(s,p1,p2,p3,r,dr)
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
!***********************************************************************c     input:  s saturation
c         p1 residual sat
c         p2 max sat
c         p3 = exponent
c     return:  r and its derivatives (dr1,dr2)
       
      implicit none
      real*8 s,p1,p2,p3,r,dr,star
      star=(s-p1)/(p2-p1)
      r= star**p3
      dr = p3*star**(p3-1.0)/(p2-p1)
c     if((s-p1).eq.0.d0) dr=0.d0
      if((s-p1).le.0.d0) dr=0.d0
      if((s-p1).le.0.d0) r=0.d0
      
      end subroutine exponential

