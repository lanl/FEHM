      subroutine thomeercap(iflg,i,it,sl,star,pc,dpcs,rl,drls,rv,drvs)
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
!D1 Calculate Thomeer relative permeabilities and capillary pressures
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/thomeercap.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.4.4 Relative-permeability and capillary-pressure functions
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comhi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comki
      implicit none

      integer iflg,i,it
      real*8 sl,star,hp,dhp,rl,drls,rv,drvs
      real*8 rp1,rp2,rp3,rp4,denom,ds
      real*8 pd1,fg1,term1,dterm1,pc,pcd,dpcs,cp1,cp2,cp3,apc

      if(iflg.eq.0) then
c initialization and input
c   rp20f- storage for gridblock max saturations
c   rp21f- storage for gridblock local minima
       if(allocated(rp20f)) then
        deallocate (rp20f,rp21f)
        allocate(rp20f(n0))
        allocate(rp21f(n0))
        rp20f = 0.0
        rp21f = 0.0
       else
        allocate(rp20f(n0))
        allocate(rp21f(n0))
        rp20f = 0.0
        rp21f = 0.0
       endif
      else if(iflg.eq.-2) then
c     two-phase mixture (test)
c     linear capillary function
c     
c     linear forsythe(1989) model
c     max cap p = 1 Mpa
c     cap = 0 at sl = 1.  
                  cp1=1.0
                  cp3=1.0
                  pc=cp1*(cp3-sl)
                  dpcs=-cp1

      else if(iflg.eq.-3) then
c     two-phase mixture (test)
c     linear rel  function
               rp1 = 0.
               rp2 = 0.
               rp3 = 1.
               rp4 = 1.
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  sv = 1.0-sl
                  if(sl.ge.rp1) then
                     rl = 1.0
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom
                        drls = 1.0/denom
                     endif
                  endif
                  if(sv.ge.rp2) then
                     rv = 1.0
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom
                        drvs=-1.0/denom
                     endif
                  endif
      else if(iflg.eq.1) then
c     Thomeer capillary relationship
c calculate max and min capillary pressures
               cp1 = cp1f(it)
               cp2 = cp2f(it)
               pd1 = rp3f(it)
               fg1 = rp4f(it)
               term1 = -fg1/(log(1.0-cp1))
               pcd = pd1*10.**term1
               rp5f(it) = pcd
               term1 = -fg1/(log(1.0-cp2))
               pcd = pd1*10.**term1
               rp6f(it) = pcd
               rp7f(it) = -pcd/(1.0-cp2)

      else if(iflg.eq.2) then
c     Thomeer capillary relationship
c     cp1-cutoff lower saturation
c     cp2-cutoff upper saturation
c     rp3-Pd1 parameter in Boitnott's formulation
c     rp4-Fg1 parameter in Boitnott's formulation
c     rp5-Pc at lower saturation
c     rp6-Pc at upper saturation
               cp1 = cp1f(it)
               cp2 = cp2f(it)
               pd1 = rp3f(it)
               fg1 = rp4f(it)

                  if(sl.lt.cp1) then 
                   pc = rp5f(it) 
                   dpcs = 0.0
                  elseif(sl.gt.cp2) then
                   apc = rp7f(it)
                   pc = apc*(1.-sl)+rp6f(it) 
                   dpcs = -apc 
                  else
                   term1 = -fg1/(log(1.0-sl))
                   dterm1 = -fg1/(log(1.0-sl))**2*(1.0/(1.0-sl))
                   pc = pd1*10.**term1
                   dpcs = 2.3*pc*dterm1
c                  dpcs = pd1*term1*10.**(term1-1.)*dterm1
                   continue
                  endif

      else if(iflg.eq.3) then
c
c calculated relative permeabilities
c corey relationships(geothermal)
c rp1 is residual liquid saturation
c rp2 is residual gas saturation
c
                  rp1 = rp1f(it)
                  rp2 = rp2f(it)
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  denom = 1.0-rp1-rp2
                  star=(sl-rp1)/denom
                  ds = 1.0/denom
                  if(star.le.0.0) then
                     rv = 1.0
                  else if(star.ge.1.0) then
                     rl = 1.0
                  else
                     rl = star**4
                     drls = 4.0*star**3*ds
                     rv=(1.0-star)**2*(1.0-star**2)
                     drvs=(-2.0*(1.0-star)*(1.0-star**2)+
     2                    (1.0-star)**2*(-2.0*star))*ds
                  endif


      endif
      return
      end
