      subroutine brooks_corey_rlp3(sw, sl, rp1, rp7, rp3, prop2,
     &     dprop21, dprop22, prop3, dprop31, dprop32)
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
!D1 Brooks-Corey 3-phase relative permeability model.
!D1
!***********************************************************************
! 3-phase rel perm model  (presumes you've called 2-phase model first)
! input:  sw is wetting phase saturation
!         sl is liquid phase saturation
!         rp1 is residual wetting phase saturation
!         rp3 is lambda (Li and Horne, 2006) 
!	  rp7 is residual CO2 liquid sat 
! output: prop2 is rel perm (wetting phase saturation)
!	  dprop21 is dr/ds
!	  dprop22 is dr/ds(gas)
!	  prop3 is rel perm (non-wetting phase saturation)
!	  dprop31 is dr/ds
!	  dprop32 is dr/ds(gas)
!***********************************************************************

      implicit none
      real*8 s,p1,p2,prop2,dprop21,dprop22,prop3,dprop31,dprop32
      real*8 rp31,sw_star,sl_star,sw,sl,rp1,rp3,rp7,rp32,dsl_stardw,
     &     dsl_stardg
      real*8 sc_star,dsc_stardw,dsc_stardg,dsw_stardg
c     for 3 phase, called after brooks-corey_rlp  
      sl_star = (sw+sl-rp1)/(1.d0-rp1)  
      sw_star = (sw-rp1)/(1.d0-rp1)
      rp31=((2.d0+3.d0*rp3)/rp3)
      rp32=(2.d0+rp3)/rp3
      dsl_stardw = 0.d0
      dsl_stardg = -1.d0/(1.d0-rp1)
      sc_star = (sl-rp7)/(1-rp1-rp7)
      dsc_stardw = -1.d0/(1-rp1-rp7)
      dsc_stardg = -1.d0/(1-rp1-rp7)
      dsw_stardg = 0.d0
      if(sl_star.le.0.d0) then
         prop3 = 1.d0
         dprop31 = 0.d0
         dprop32 = 0.d0
      elseif(sl_star.ge.1.d0) then
         prop3 = 0.d0
         dprop31 = 0.d0
         dprop32 = 0.d0
      else
         prop3=((1.d0-sl_star)**2.d0)*(1.d0-(sl_star**rp32))
         dprop31 = (-2.d0*(1.d0-sl_star)*
     &        (1.d0-(sl_star**rp32))+((1.d0-sl_star)**2.d0)*
     &        (-rp32*sl_star**(rp32-1.d0)))*(dsl_stardw)
         dprop32 = (-2.d0*(1.d0-sl_star)*
     &        (1.d0-(sl_star**rp32))+((1.d0-sl_star)**2.d0)*
     &        (-rp32*sl_star**(rp32-1.d0)))*(dsl_stardg)
      endif
      if(sc_star.le.0.d0) then
         prop2 = 1.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
      elseif(sc_star.ge.1.d0) then
         prop2 = 0.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
      else
         prop2=(sc_star**2.d0)*((sl_star**rp32)-
     &        (sw_star**rp32))
         dprop21=((sl_star**rp32)-(sw_star**rp32))*
     &        2.d0*sc_star*dsc_stardw+(sc_star**2.d0)*rp32*
     &        (sl_star**(rp32-1.d0)*dsl_stardw-
     &        sc_star**(rp32-1.d0)*dsc_stardw)
         dprop22=((sl_star**rp32)-(sw_star**rp32))*
     &        2.d0*sc_star*dsc_stardg+(sc_star**2.d0)*rp32*
     &        (sl_star**(rp32-1.d0)*dsl_stardg-
     &        sc_star**(rp32-1.d0)*dsc_stardg)
      endif

      end subroutine brooks_corey_rlp3
