      subroutine vgcap_fit3(iflg,slr,slm,slcut,smcut,fac,alpha
     &                     ,alamda,c1,c2,c3,c4,hmin)

      implicit none

      integer iflg,ncut 
      real*8 slr,slm,slcut,alpha,alamda,c1,c2,c3,c4 
      real*8 hmax,hmin,hcut,smcut,alpi,ds,dhp,slope
      real*8 fac,fac_min,fac_use,hcut2,smcut2,ancut
      parameter(fac_min=2.0,ncut = 1)

      if(iflg.eq.0) then
       if(fac.eq.0.0) then
c zero slope, zero curvature at sl=0
c h = a*x**3 + d
c this is most common fehm way
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
            ds = 1.0/(slm-slr)
            dhp = 1.0/alpha*(1.0-alamda)/
     *              (1.0/smcut**alpi-1.0)**alamda
     *              *(-alpi/smcut**(alpi+1.0))*ds
            fac=1.0/3.0
       		hmax=hcut-fac*dhp*slcut
       		c1=dhp/(3.0*slcut**2)
       		c2=0.0
       		c3=0.0
       		c4=hmax
       else if(fac.gt.0.0) then
c linear fit from slcut to sl=0.0
c with hmax=fac*hcut
            fac_use=max(fac,fac_min)
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
       		hmax=fac_use*hcut
       		slope=-(hmax-hcut)/slcut
       		c3=slope
       		c4=hmax
       else 
c linear fit from slcut to sl=0.0
c with hmax: linear extension of slope at slcut
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
            ds = 1.0/(slm-slr)          
            dhp = 1.0/alpha*(1.0-alamda)/
     *              (1.0/smcut**alpi-1.0)**alamda
     *              *(-alpi/smcut**(alpi+1.0))*ds
       		hmax=hcut-dhp*slcut
       		slope=-(hmax-hcut)/slcut
       		c3=slope
       		c4=hmax
      endif
      elseif(iflg.eq.4) then          
c linear fit from smcut to star= 1.00
c with hmax: linear extension of slope at slcut
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
       slope=(hcut-hmin)/(1.0-slcut)
       c3=-slope
       c4=slope+hmin
       ancut = (1./slcut**ncut-1.0)
       c3 = hcut/ancut
       c4 = -c3 + hmin
      continue
      endif
      return
      end
