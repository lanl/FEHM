      subroutine brooks_corey_rlp(s, p1, lambda, prop1, dprop11, 
     &     dprop12, prop2, dprop21, dprop22)
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
!D1 Calculate Brooks-Corey relative permeability
! 2-phase 
! input:  s is wetting phase saturation
!         p1 is residual wetting phase saturation       
!         lambda (Li and Horne, 2006) 
! output: prop1 is rel perm (wetting phase saturation)
!	  dprop11 is dr/ds
!	  dprop12 is dr/ds(gas)
!	  prop2 is rel perm (non-wetting phase saturation)
!	  dprop21 is dr/ds
!	  dprop22 is dr/ds(gas)
!***********************************************************************

c      use comdi
      use comrlp

      implicit none
      real*8 s, p1, p2, lambda, prop1, dprop11, dprop12
      real*8 rp31, sw_star, ds, prop2, dprop21, dprop22

      rp31 = (2.d0+3.d0*lambda)/lambda
      sw_star = (s-p1)/(1.d0-p1)
      ds = 1.d0/(1.d0-p1)
c     calculate rel perms
      if(sw_star.le.0.d0) then
         prop1 = 0.d0
         dprop11 = 0.d0
         dprop12 = 0.d0
         prop2 = 1.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
      elseif(sw_star.ge.1.d0) then
         prop1 = 1.d0
         dprop11 = 0.d0
         dprop12 = 0.d0
         prop2 = 0.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
      else
         prop1 = sw_star**rp31
         dprop11=rp31*((s-p1)**(rp31-1))/((1-p1)**rp31)
c     check this wtih gaz, raj
         dprop12= -1.* dprop11
         prop2 = (1.d0-sw_star)**2.*(1.0-sw_star**rp31)
         dprop21=ds*(2.0*(1.0-sw_star)*(1.0-sw_star**rp31)+
     $        (1.0-sw_star)**2.0*(-rp31*sw_star**(rp31-1.0)))
         
      endif
      
      end subroutine brooks_corey_rlp
