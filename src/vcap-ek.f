      subroutine vgcap_ek	( sl, slr, smr,alpha, beta,smcut, 
     2     sucut,hp, dhp)
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

! this calculates capillary pressure (hp) as a function of liquid saturation (sl)
! slr = minimum residual saturation
! smr = maximum residual saturation
! alpha is 1/m (inverse air entry pressure) = Pa/(density*g)
! beta is 'm' in VG 1980 paper and 'n' in FEHM MMS description
! smcut, sucut are limits of the function -  if sl< smcut or sl > sucut,
! linear interpolation is used.  
! subroutine returns hp (capillary pressure) and dhp (dhp/dsl)
! 
      implicit none
      integer ireg
      real*8 sl
      real*8 slr
      real*8 smr
      real*8 alpha
      real*8 beta
      real*8 smcut
      real*8 sucut
      real*8 hp
      real*8 dhp,star2,hp_cut
      real*8 lambda
      real*8 ds
      real*8 denom
      real*8 star
     
      lambda = 1.0-1.0/beta
      denom = smr-slr
      ds = 1.0/denom
      if(sl.gt.smcut.and.sl.lt.sucut) then
      	 star=(sl-slr)/denom
      	 call cc(star,lambda,alpha,denom,hp,dhp)
      elseif(sl.le.smcut)  then
         star=(smcut-slr)/denom

       	 call cc(star,lambda,alpha,denom,hp_cut,dhp)
         dhp=hp_cut/smcut
         hp=hp_cut+(smcut-sl)*dhp
         dhp=-1.*dhp
      else
      
         star =(sucut-slr)/denom
         call cc(star,lambda,alpha,denom,hp_cut,dhp)
         dhp=hp_cut/(1.d0-sucut)
         hp=hp_cut-(sl-sucut)*dhp
         dhp=-1.*dhp
      endif
!      write(*,*) 'line 63 ',sl,slr,smr,alpha,beta,hp
      
      return
      end
      subroutine cc(star,lambda,alpha,denom,hp,dhp)
      implicit none
      real*8 star,lambda,alpha,denom,hp,chp,termstar1,termstar2,termb1
      real*8 termb2,dhp
         termstar1 = star**(1./lambda)
         termstar2 = star*termstar1
         termb1 = 1./termstar1 - 1.
         termb2 = termb1**(-lambda)  
         hp = 1./alpha * termb1*termb2  
         dhp = (1.-lambda)/alpha * termb2 / (-1.*termstar2*lambda*denom)
         return
         end
      
      
      
      
      
      
      
      
      
