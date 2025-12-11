      subroutine brooks_corey_cap(s, p1, p2, lambda, cmax, cut1, cut2,
     &     prop1, dprop11)
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
!D1 calculate Brooks-Corey capillary pressure
! 2-phase 
! input:  s is wetting phase saturation
!         p1 is residual wetting phase saturation        
!         p2 is max wetting phase saturation 
!         cut1,cut2 are cutoff parameters 
!  cut2 is saturation value; if s > cut2; brooks-corey, otherwise, linear
!  cut1 is only used for s < cut2 calculations
!  cp at s=0 will be cut1*cp(at s=cut2)
!
!         cmax = entry pressure 
!         lambda (Li and Horne, 2006) 
! output: prop1 is cap pressure
!		dprop11 is dcp/ds
!***********************************************************************

c      use comdi
      use comrlp

      implicit none
      real*8 s, p1, p2, lambda, prop1, dprop11
      real*8 rp31, sw_star, ds, fac_use, fac_min
      real*8 smcut, smcutm, cut1, cut2, slcut, hcut, cmax
      real*8 slope, dhp, hmax, hmin
      integer iflag, fac
      parameter (fac_min=2.d0)
      parameter (hmin=0.d0)

      rp31 = (2.d0+3.d0*lambda)/lambda
      sw_star = (s-p1)/(1.d0-p1)
      ds = 1.d0/(1.d0-p1)
c     calculate cap pressures
c     rp12 is entry pressure
c     sucut = 1.d0
c     rp15, rp16 are cutoff fitting parameters
      smcut=(cut2-p1)/(p2-p1)
      smcutm=max(smcut,1.d-3)
      slcut=smcutm*(p2-p1)+p1
      if(s.gt.0.d0) then
         if(s.le.cut2) then
            hcut=cmax*(smcutm**(-1.d0/lambda))
            if(cut1.gt.0.0) then
c     a multiple of value of
               fac_use=max(cut1,fac_min)
               hmax=hcut*cut1
               slope=-(hmax-hcut)/slcut
            else
c     linear fit from saturation curve
               ds=1.d0/(p2-p1)
               dhp=cmax*(-1.d0/lambda)*
     $              smcutm**((-lambda-1.d0)/lambda)*ds
               hmax=hcut-dhp*smcutm
               slope=-(hmax-hcut)/smcutm
            endif
            prop1=slope*s+hmax
            dprop11=slope
         else
            prop1 = cmax*((s-p1)/(p2-p1))**(-1.d0/lambda)
            dprop11 = cmax*(-1.d0/lambda)*
     $           ((s-p1)**((-lambda-1.d0)/lambda))/
     $           ((p2-p1)**(-1.d0/lambda))
         endif
      elseif(s.le.0.d0) then
         hcut=cmax*(smcutm**(-1.d0/lambda))
         if(cut1.gt.0.0) then
c     a multiple of value of
            fac_use=max(cut1,fac_min)
            hmax=hcut*cut1
c           write(45,*) 'fac_use',hcut,cut1
         else
c     linear fit from saturation curve
            ds=1.d0/(lambda-p1)
            dhp=cmax*(-1.d0/lambda)*
     $           smcutm**((-lambda-1.d0)/lambda)*ds
            hmax=hcut-dhp*smcutm
c           write(45,*) 'linear fit',ds,dhp,hmax
            
         endif
         prop1=hmax
         dprop11=0.d0
      endif	

      end subroutine brooks_corey_cap
