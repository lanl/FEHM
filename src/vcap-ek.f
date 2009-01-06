      subroutine vgcap_ek	( sl, slr, smr,alpha, beta,smcut, 
     2	sucut,hp, dhp)
! this calculates capillary pressure (hp) as a function of liquid saturation (sl)
! slr = minimum residual saturation
! smr = maximum residual saturation
! alpha, beta are defined in other vgcap subroutines
! smcut, sucut are limits of the function -  if sl< smcut or sl > sucut,
! linear interpolation is used.  
! cp(sl=0) = 2 X cp(sl=smcut)
! subroutine returns hp (capillary pressure) and dhp (dhp/dsl)
! 
      implicit none
      integer ireg
      real*8 sl
      real*8 slr
      real*8 smr
      real*8 alpha
      real*8 beta
      real*8 ac1
      real*8 ac2
      real*8 ac3
      real*8 ac4
      real*8 bc3
      real*8 bc4
      real*8 smcut
      real*8 sucut
      real*8 hp
      real*8 dhp,star2,hp_cut


      real*8 alamda
      real*8 ds
      real*8 denom
      real*8 star
      real*8 alpi
      real*8 termstar1
      real*8 termstar2
      real*8 termb1
      real*8 termb2
      alamda = 1.0-1.0/beta
      alpi = 1.0/alamda
      denom = smr-slr
      star=(sl-slr)/denom
      ds = 1.0/denom
      if(sl.gt.smcut.and.sl.lt.sucut) then
            termstar1 = star**alpi
            termstar2 = star*termstar1
            termb1 = 1./termstar1 - 1.
            termb2 = termb1**(-alamda)
            hp = 1./alpha * termb1*termb2
            dhp = (1.-alamda)/alpha * termb2 * (-alpi/termstar2) * ds
         elseif(sl.le.smcut)  then
	star2=(smcut-slr)/denom
            termstar1 = star2**alpi
            termstar2 = star2*termstar1
            termb1 = 1./termstar1 - 1.
            termb2 = termb1**(-alamda)
            hp_cut = 1./alpha * termb1*termb2
	dhp=hp_cut/smcut
	hp=hp_cut+(smcut-sl)*dhp
	dhp=-1.*dhp
	else
	star2=(sucut-slr)/denom
            termstar1 = star2**alpi
            termstar2 = star2*termstar1
            termb1 = 1./termstar1 - 1.
            termb2 = termb1**(-alamda)
            hp_cut = 1./alpha * termb1*termb2
	dhp=hp_cut/(1.d0-sucut)
	hp=hp_cut-(sl-sucut)*dhp
	dhp=-1.*dhp
	endif
       if(sl.lt.0.0.or.sl.gt.1.0) then
c     lower residual cutoff
         hp = 0.
         dhp= 0.0
      endif
      return
      end




























































































































































































































































































































