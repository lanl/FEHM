      subroutine vgrlpa( sl, star, alpha, beta, hmin,
     2     hp, dhp, rpa1, rpa2, rpa3, rpa4, rpa5,
     3     slr, slm, rl, drls, rv, drvs )
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To compute the liquid relative permeability and derivative for the
!D1 van Genuchten model.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/vgrlpa.f_a  $
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.4.4 Relative-permeability and capillary-pressure functions
!D3
!**********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      implicit none

      real*8 sl  
      real*8 star
      real*8 alpha
      real*8 beta
      real*8 hp
      real*8 dhp
      real*8 rl
      real*8 drls
      real*8 rv
      real*8 drvs
      real*8 hmin


      real*8 alamda
      real*8 dahp
      real*8 term1
      real*8 ahp
      real*8 al2
      real*8 dterm1
      real*8 ratio
      real*8 term2
      real*8 terma
      real*8 dratio
      real*8 dterma
      real*8 dterm2
      real*8 termhp1
      real*8 termhp2
      real*8 ahp1
      real*8 termahp1
      real*8 termahp2
      real*8 rpa1, rpa2, rpa3, rpa4, rpa5
      real*8 term3,term4,term5,term6
      real*8 slr, slm, hpc
c     relative perm functions from Touma and Vauclin (1986)
c     as reported by Celia and Binney (1992)
c
c check for minimum capillary pressure
c
      if (hp .lt. hmin ) hp=hmin
      if(sl.gt.slr.and.sl.lt.slm) then
c     calculate the liquid relative permeability
         rl = rpa5*sl**rpa1
         drls=rpa5*rpa1*sl**(rpa1-1.0d00)           
      else if(sl.le.slm) then
c     lower residual cutoff
         rl = 0.0
         drls= 0.0
c     upper residual cutoff
      else
         rl = 1.0
         drls= 0.0
      endif
c     calculate the gas relative permeability
c     convert head from meters to centimeters
      hpc = hp*100.d00
c     1.0d02 factor incorporates 100.d00 term above 
      rv = rpa2*rpa3/(rpa3 + hpc**rpa4)
      drvs = -(rpa2*rpa3/(rpa3 + 1.d00/hpc**(-rpa4))**2)*
     &        rpa4*1.0d02/hpc**(1.d00-rpa4)*dhp
      return
      end
