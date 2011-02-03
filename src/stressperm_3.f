      subroutine stressperm_3(i)
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

c     perm model 3 fully coupled model
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

      real*8 frac_bx,frac_by,frac_bz,disx,disy,disz
      real*8 amultx,amulty,amultz,amultxy,amultxz,amultyz
      real*8 dispx12,dispy12,dispz12, alpi,tx12,ty12,tz12
      real*8 disptx12,dispty12,disptz12,dispxy12,dispxz12,dispyx12
      real*8 dispyz12,dispzx12,dispzy12
      
c     x direction (orthogonal pieces)
c     
c     identify parameters
c     
c     
      strx_min = spm1f(iispmd)
      stry_min = spm2f(iispmd) 
      strz_min = spm3f(iispmd) 
      frac_bx = spm4f(iispmd) 
      frac_by = spm5f(iispmd) 
      frac_bz = spm6f(iispmd)
      perx_m = spm7f(iispmd)
      pery_m = spm8f(iispmd) 
      perz_m = spm9f(iispmd) 
      
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
      amultx = frac_bx**3
      amulty = frac_by**3
      amultz = frac_bz**3
      amultxy = 1./(amultx + amulty)
      amultxz = 1./(amultx + amultz)
      amultyz = 1./(amulty + amultz)
c     
c     displacement terms (orthogonal)
c     
      dispx12 = du(kbx2)-du(kbx1)
      dispy12 = dv(kby2)-dv(kby1)  
      dispz12 = dw(kbz2)-dw(kbz1)  
c     
c     displacement terms (thermal)
c     
      alpi = alp(i)     
      tx12 = (t(i)-tini(i))*alpi*disx/2. 
      ty12 = (t(i)-tini(i))*alpi*disy/2. 
      tz12 = (t(i)-tini(i))*alpi*disz/2. 
c     
c     determine the net contribution  
c     
      disptx12 = max(dispx12-tx12,0.0d0)
      dispty12 = max(dispy12-ty12,0.0d0) 
      disptz12 = max(dispz12-tz12,0.0d0)              
c     
c     displacement terms (shear)
c     
      dispxy12 = du(kby2)-du(kby1)
      dispxz12 = du(kbz2)-du(kbz1)  
      dispyx12 = dv(kbx2)-dv(kbx1)
      dispyz12 = dv(kbz2)-dv(kbz1)
      dispzx12 = dw(kbx2)-dw(kbx1)
      dispzy12 = dw(kby2)-dw(kby1)                         
      
      
      rlxs(i) = amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3
      drlxs(i,2) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**2*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
      drlxs(i,1) = -3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**2*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amulty
      drlxs(i,4) = 3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
      drlxs(i,3) = -3.*amultyz*(1. + amulty*(dv(kby2)-dv(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
      
      rlys(i) = amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3
      drlys(i,2) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
      drlys(i,1) = -3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**3*amultx
      drlys(i,4) = 3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz
      drlys(i,3) = -3.*amultxz*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amultz*(dw(kbz2)-dw(kbz1)))**2*amultz  
      
      rlzs(i) = amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amulty*(dv(kbz2)-dv(kbz1)))**3
      drlzs(i,2) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &     (1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
      drlzs(i,1) = -3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**2*
     &     (1. + amulty*(dv(kbz2)-dv(kbz1)))**3*amultx
      drlzs(i,4) = 3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty
      drlzs(i,3) = -3.*amultxy*(1. + amultx*(du(kby2)-du(kby1)))**3*
     &     (1. + amulty*(dv(kbz2)-dv(kbz1)))**2*amulty     
      
      return
      
      end    
cc....................................................................
