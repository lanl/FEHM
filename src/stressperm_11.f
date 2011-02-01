      subroutine stressperm_11(i)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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

c     s kelkar June 9 2009, simple gangi (1978, Int.J.Rock Mech)
c     model in 2-D
c     assumed model in x-y plane and fracture in x-z plane, so only 
c     x-perm is modified. This formulation is for fracture faces
c     in contact, so str_y = effective stress in y-dir has to be
c     compressive

      use comai
      use combi
      use comdi
      use comsi

      implicit none
      integer i, iispmd
      real*8 gk0, gpmod, gmexp, gn, sigy_eff

      iispmd = ispm(i)    

      if(icnl.ne.0) then
         gk0  =spm1f(iispmd)
         gpmod=spm2f(iispmd)
         gmexp =spm3f(iispmd)
         gn=1./gmexp
         sigy_eff=str_y(i)
         if(sigy_eff.gt.0.d0) then
            pnx(i)=gk0*(1.-(sigy_eff/gpmod)**gmexp)**3.
            e2(i) =  gn*sigy_eff*(gpmod/sigy_eff)**gn
         endif
      endif
      
      return
      
      end
c.....................................................................
