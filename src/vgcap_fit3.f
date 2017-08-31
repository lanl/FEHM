      subroutine vgcap_fit3(iflg,slr,slm,slcut,smcut,fac,alpha
     &                     ,alamda,c1,c2,c3,c4,hmin)
!*************************************************************************
! Copyright  2015.   Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
! National  Laboratory  (LANL),  which is operated by  Los Alamos National
! Security, LLC  for the U. S. Department of Energy.  The U. S. Government
! has rights to use, reproduce, and distribute this software.  Neither the
! U. S. Government nor Los Alamos National Security, LLC or persons acting
! on their behalf,  make any warranty,  express or implied, or assumes any
! liability for the accuracy, completeness, or usefulness of the software,
! any information pertaining to the software,  or  represents that its use
! would not infringe privately owned rights.

! The  software  being licensed  may  be Export Controlled.  It may not be
! distributed  or  used by individuals  or entities prohibited from having
! access to the software package, pursuant to United States export control
! laws and regulations. An export control review and determination must be
! completed before LANS will provide access to the identified Software.
!*************************************************************************

      implicit none

      integer iflg,ncut 
      real*8 slr,slm,slcut,alpha,alamda,c1,c2,c3,c4 
      real*8 hmax,hmin,hcut,smcut,alpi,ds,dhp,slope
      real*8 fac,fac_min,fac_use,hcut2,smcut2,ancut
      parameter(fac_min=2.0,ncut = 1)
		c1=0.;c2=0.;c3=0.;c4=0.
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
