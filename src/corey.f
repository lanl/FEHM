      subroutine corey(sl,rp1,rp2,rl,drls,rv,drvs)
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
!D1 Corey relative permeability model.
!D1
!***********************************************************************

      implicit none
      real*8 sl,rp1,rp2,rl,drls,rv,drvs,star,denom,ds
c     this is the formulation in Gaz's original rlperm.f (geothermal)
c     and in the models and methods document
c     and in Preuss(2002)
c     rp1 and rp2 are both residual saturations

      rl = 0.0
      rv = 0.0
      drls = 0.0
      drvs = 0.0
      denom = 1.0-rp1-rp2
      star=(sl-rp1)/denom
      ds = 1.0/denom
      if(star.lt.0.d0) then
         rv = 1.0
      else if(star.gt.1.d0) then
         rl = 1.0
      else
         rl = star**4
         drls = 4.0*star**3*ds
         rv=(1.0d0-star)**2*(1.0d0-star**2)
         drvs=(2.0*(1.0-star)*(1.0-star**2)+
     $        (1.0-star)**2*(-2.0*star))*(-ds)
      endif

      end subroutine corey
