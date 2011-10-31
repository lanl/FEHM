      subroutine stressperm_5(i)
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

c     thermal contribution based formulation
c     perm model 5 displacement based model-explicitly coupled
c     
c     perm changes for individual node
c     assumes call has been made to allocate memory
      use comai
      use combi
      use comdi
      use comsi

      implicit none
      integer i
      integer kbx1,kbx2,kby1,kby2,kbz1,kbz2,idir,iispmd

      real*8  ksx, ksy, ksz
      real*8  epsxx12, epsyy12, epszz12 
      real*8  epsxy12, epsxz12, epsyz12
      real*8  sisjx, sisjy, sisjz           
      real*8 dispx12, dispy12, dispz12
      real*8 disptx12, dispty12, disptz12
      real*8 dispxy12, dispxz12, dispyx12, dispyz12, dispzx12, dispzy12
      real*8 frac_bx, frac_by, frac_bz, frac_tol
      real*8 amultx, amulty, amultz
      real*8 amultxy, amultxz, amultyz
      real*8 disx,disy,disz,alpi,dt,efac,phid,pi

      parameter (frac_tol=0.00001, pi=3.1415926535)
c     
c     x direction (orthogonal pieces)
c     
c     identify parameters (model 5 uses only the rock strengths)
c     
c     
      iispmd = ispm(i)

      strx_min = spm1f(iispmd)
      stry_min = spm2f(iispmd) 
      strz_min = spm3f(iispmd)
      phid = spm6f(iispmd)
c     
c     the parameters in the fully coupled model are replaced 
c     by initial permeabilities
c     
      
      frac_bx = max(spm4f(iispmd),frac_tol) 
      frac_by = max(spm4f(iispmd),frac_tol) 
      frac_bz = max(spm4f(iispmd),frac_tol) 
c     perx_m = spm7f(iispmd)
c     pery_m = spm8f(iispmd) 
c     perz_m = spm9f(iispmd) 
      
      kbx1 = ipermx(i,1)
      kbx2 = ipermx(i,2)
      kby1 = ipermy(i,1)
      kby2 = ipermy(i,2)
      kbz1 = ipermz(i,1)
      kbz2 = ipermz(i,2)
c     
c     identify displacements
c     
      disx = (cord(kbx2,1)-cord(kbx1,1))/2.
      disy = (cord(kby2,2)-cord(kby1,2))/2.
      disz = (cord(kbz2,3)-cord(kbz1,3))/2.
c     
c     
c     different definitions from model 3
c     
      amultx = 1./(disx*frac_bx)
      amulty = 1./(disy*frac_by)
      amultz = 1./(disz*frac_bz)
      amultxy = 1.
      amultxz = 1.
      amultyz = 1.              
c     
c     displacement terms (orthogonal)
c     
      dispx12 = 
     &     (du(kbx2)-du(kbx1)-(du_ini(kbx2)-du_ini(kbx1)))/2.
      dispy12 = 
     &     (dv(kby2)-dv(kby1)-(dv_ini(kby2)-dv_ini(kby1)))/2.  
      dispz12 = 
     &     (dw(kbz2)-dw(kbz1)-(dw_ini(kbz2)-dw_ini(kbz1)))/2. 
      epsxx12 = dispx12/disx
      epsyy12 = dispy12/disy
      epszz12 = dispz12/disz 
c     
c     displacement terms (thermal)
c     
      alpi = alp(i)     
      dt = t(i)-tini(i)
      efac = 3.d0*e2(i)+2.d0*e3(i)        
c     tx12 = (t(i)-tini(i))*alpi*disx
c     ty12 = (t(i)-tini(i))*alpi*disy
c     tz12 = (t(i)-tini(i))*alpi*disz 
      
c     
c     determine the net contribution  
c     
      
      dispx12= dispx12
     &     +(e1(i)*epsxx12+e2(i)*(epsyy12+epszz12)
     &     -alpi*efac*dt)/knx_stressperm      
      
      dispy12 = dispy12
     &     +(e1(i)*epsyy12+e2(i)*(epsxx12+epszz12)
     &     -alpi*efac*dt)/kny_stressperm
      
      dispz12 = dispz12
     &     +(e1(i)*epszz12+e2(i)*(epsxx12+epsyy12)
     &     -alpi*efac*dt)/knz_stressperm
      
      disptx12 = max(dispx12,0.0d0)
      dispty12 = max(dispy12,0.0d0) 
      disptz12 = max(dispz12,0.0d0)              
c     
c     displacement terms (shear)
c     
      
      epsxy12 = abs((du(kby2)-du(kby1))/2)/disy
     &     +abs((dv(kbx2)-dv(kbx1))/2)/disx          
      epsxz12 = abs((du(kbz2)-du(kbz1))/2)/disz
     &     +abs((dw(kbx2)-dw(kbx1))/2)/disx
      epsyz12 = abs((dv(kbz2)-dv(kbz1))/2)/disz
     &     +abs((dw(kby2)-dw(kby1))/2)/disy
      
      if(str_x(i)-phi(i).gt.0) then
         dispyx12 = (disy+e3(i)/ksy)*epsxy12
         dispzx12 = (disz+e3(i)/ksz)*epsxz12
      else
         dispyx12 = 0.
         dispzx12 = 0.
      endif
      
      if(str_y(i)-phi(i).gt.0) then
         dispxy12 = (disx+e3(i)/ksx)*epsxy12
         dispzy12 = (disz+e3(i)/ksz)*epsyz12  
      else
         dispxy12 = 0.
         dispzy12 = 0.
      endif
      
      if(str_z(i)-phi(i).gt.0) then
         dispxz12 = (disx+e3(i)/ksx)*epsxz12
         dispyz12 = (disy+e3(i)/ksy)*epsyz12
      else
         dispxz12 = 0.
         dispyz12 = 0.
      endif
      
c     dispxy12 = (disx+e3(i)/ksx)*epsxy12
c     dispxz12 = (disx+e3(i)/ksx)*epsxz12
c     dispyx12 = (disy+e3(i)/ksy)*epsxy12
c     dispyz12 = (disy+e3(i)/ksy)*epsyz12
c     dispzx12 = (disz+e3(i)/ksz)*epsxz12
c     dispzy12 = (disz+e3(i)/ksz)*epsyz12  
      
c     
c     perm enhancement                               
c     
      
      if(frac_by.eq.0.and.frac_bz.eq.0) then
         rlxs(i) = 1.
      else
         rlxs(i) = 1/(frac_by**3+frac_bz**3)*((frac_by+dispty12
     &        +(dispxy12+dispzy12)*tan(phid/180*pi))**3+
     &      (frac_bz+disptz12+(dispyz12+dispxz12)*tan(phid/180*pi))**3)
      endif
      
      if(frac_bx.eq.0.and.frac_bz.eq.0) then
         rlys(i) = 1.
      else
         rlys(i) = 1/(frac_bx**3+frac_bz**3)*((frac_bx+disptx12
     &        +(dispyx12+dispzx12)*tan(phid/180*pi))**3+
     &      (frac_bz+disptz12+(dispyz12+dispxz12)*tan(phid/180*pi))**3)
      endif
      
      if(frac_bx.eq.0.and.frac_by.eq.0) then
         rlzs(i) = 1.
      else
         rlzs(i) = 1/(frac_bx**3+frac_by**3)*((frac_bx+disptx12
     &        +(dispyx12+dispzx12)*tan(phid/180*pi))**3+
     &      (frac_by+dispty12+(dispxy12+dispzy12)*tan(phid/180*pi))**3)
      endif
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccCHECK          
      check(i)= (dispxy12+dispzy12)*tan(phid/180*pi)
      
c     rlxs(i) = amultyz*(1. + amulty*dispty12)**3*
c     &              (1. + amultz*disptz12)**3
      drlxs(i,2) = 3.*amultyz*(1. + amulty*dispty12)**2*
     &     (1. + amultz*disptz12)**3*amulty
      drlxs(i,1) = -3.*amultyz*(1. + amulty*dispty12)**2*
     &     (1. + amultz*disptz12)**3*amulty
      drlxs(i,4) = 3.*amultyz*(1. + amulty*dispty12)**3*
     &     (1. + amultz*disptz12)**2*amultz
      drlxs(i,3) = -3.*amultyz*(1. + amulty*dispty12)**3*
     &     (1. + amultz*disptz12)**2*amultz
      
c     rlys(i) = amultxz*(1. + amultx*disptz12)**3*
c     &              (1. + amultz*disptz12)**3
      drlys(i,2) = 3.*amultxz*(1. + amultx*disptz12)**2*
     &     (1. + amultz*disptz12)**3*amultx
      drlys(i,1) = -3.*amultxz*(1. + amultx*disptz12)**2*
     &     (1. + amultz*disptz12)**3*amultx
      drlys(i,4) = 3.*amultxz*(1. + amultx*disptz12)**3*
     &     (1. + amultz*disptz12)**2*amultz
      drlys(i,3) = -3.*amultxz*(1. + amultx*disptz12)**3*
     &     (1. + amultz*disptz12)**2*amultz  
      
c     rlzs(i) = amultxy*(1. + amultx*disptz12)**3*
c     &              (1. + amulty*dispty12)**3
      drlzs(i,2) = 3.*amultxy*(1. + amultx*disptz12)**2*
     &     (1. + amulty*dispty12)**3*amultx
      drlzs(i,1) = -3.*amultxy*(1. + amultx*disptz12)**2*
     &     (1. + amulty*dispty12)**3*amultx
      drlzs(i,4) = 3.*amultxy*(1. + amultx*disptz12)**3*
     &     (1. + amulty*dispty12)**2*amulty
      drlzs(i,3) = -3.*amultxy*(1. + amultx*disptz12)**3*
     &     (1. + amulty*dispty12)**2*amulty     
      
c     
c     thermal contribution based formulation
c     
c     
c     now change absolute permeabilities         
      
      pnx(i) = pnx0(i)*rlxs(i)
      pny(i) = pny0(i)*rlys(i)
      pnz(i) = pnz0(i)*rlzs(i)
      
      return
      
      end
c....................................................................
