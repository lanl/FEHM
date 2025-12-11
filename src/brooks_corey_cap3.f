      subroutine brooks_corey_cap3(sg, sr_g, smax_g, cmax, lambda, 
     &     prop2, dprop21, dprop22)
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
!D1 Brooks-Corey 3-phase capillary pressure model.
!D1
!***********************************************************************
! 3-phase cap pressure model 
! input:  sg is gas saturation
!         sr_g is residual gas saturation
!         smax_g is max gas saturation
!         cmax is maximum cap pressure
!         lambda is exponent
!***********************************************************************
	
      implicit none
      real*8 sg,sr_g,smax_g,cmax,lambda,prop2,dprop21,dprop22
      if(sg.le.sr_g) then
         prop2 = 0.d0
         dprop21 = 0.d0
         dprop22 = 0.d0
      elseif(sg.le.lambda) then
         if(lambda.eq.1) then
            prop2 = cmax*(sg-sr_g)/(smax_g-sr_g)
            dprop21 = 0.d0
            dprop22 = cmax/(smax_g-sr_g)
         else
            prop2 = cmax*(sg-sr_g)**lambda/(smax_g-sr_g)
            dprop21 = 0.0
            dprop22 = cmax*lambda*((sg-sr_g)**(lambda-1))/(smax_g-sr_g)
         endif
      else
         prop2 = cmax
         dprop21 = 0.d0
         dprop22 = 0.d0
      endif

      end subroutine brooks_corey_cap3
