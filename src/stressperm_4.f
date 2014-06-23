      subroutine stressperm_4(i)
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

c     perm model 4 - keita model
c     differs from perm model 2 by allowing cubic variation 
c     of perm with stress  
c     lagged permeability (only after each time step)
c     only good for 2-D models
      use comai
      use combi
      use comdi
      use comsi
      
      implicit none
      integer i,iispmd

      real*8  ipmd4_p1,ipmd4_p2,ipmd4_p3,ipmd4_p4
      real*8  ipmd4_k
c     
c     calculate components of volume strain
c     
c     dt = t(i) - tini(i)
c     alpv = alp(i)
c     vol_temp(i) = alpv*dt*sx1(i) - vol_strain(i)*sx1(i) 
c     ipchk = 0
      
c     kx_fac = 1. + stry_min*abs(str_y(i)**3)
      ipmd4_fx = spm1f(iispmd)
      ipmd4_br = spm2f(iispmd) 
      ipmd4_bmx = spm3f(iispmd) 
      ipmd4_alx = spm4f(iispmd) 
      ipmd4_aly = spm5f(iispmd) 
      ipmd4_fdx = spm6f(iispmd)
      ipmd4_dmx = spm7f(iispmd)
      ipmd4_gmx = spm8f(iispmd) 
      ipmd4_kc = spm9f(iispmd)
      
      ipmd4_fy = spm10f(iispmd)
      ipmd4_btx = spm11f(iispmd)
      ipmd4_bty = spm12f(iispmd)
      ipmd4_fdy = spm13f(iispmd)
      ipmd4_gmy = spm14f(iispmd)
      
c********************have to decide how to do with perx_m etc. input?   
      perx_m = 100
      pery_m = 100     
      
c***********************************NORMAL
      ipmd4_p1 = ipmd4_alx*(str_x(i)-pho(i))+ipmd4_aly*(str_y(i)-pho(i))
c     maximum aperture dilation	   
      if(ipmd4_p1.lt.-5) ipmd4_p1 = -5
      
c     normal dilation component	   
      ipmd4_p2 = ipmd4_fx/12*(ipmd4_br+ipmd4_bmx*exp(-ipmd4_p1))**3
      
      ipmd4_p3 = ipmd4_btx*(str_x(i)-pho(i))+ipmd4_bty*(str_y(i)-pho(i))
      if(ipmd4_p3.lt.-5) ipmd4_p4 = -5
      
      ipmd4_p4 = ipmd4_fy/12*(ipmd4_br+ipmd4_bmx*exp(-ipmd4_p3))**3
      
c***********************************SHEAR
c     ipmd4_k is the ratio of min stress to max stress
c     if(str_x(i).gt.str_y(i)) then
c     ipmd4_k = str_x(i)/str_y(i)
c     else 
c     ipmd4_k = str_y(i)/str_x(i)   
c     endif
c     if(ipmd4_k.gt.0) then
c     ipmd4_p3 = ipmd4_dmx*(1-exp(-ipmd4_gmx*(ipmd4_k-ipmd4_kc)))
c     else
c     ipmd4_p3 = 0
c     endif
c***********************************SHEAR
      
      pnx(i) = ipmd4_p2*1e6
      pny(i) = ipmd4_p4*1e6
      
      pnx(i) = min(perx_m*pnx0(i),pnx(i))
      pny(i) = min(pery_m*pny0(i),pny(i))
            
      return
      
      end    
c....................................................................
