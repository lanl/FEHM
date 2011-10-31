      subroutine stressperm_31(i)
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

c     change mech and flow properties with ice formation
c     based on model 2
c     changes the Youngs modulus with temperature 
c     modified from kayla Lewis' YoungsMod.f90
c     only good for 2D or 3D models
      use comai
      use combi
      use comdi
      use comsi
      implicit none

      integer i,iispmd,ipchk

      real*8 t_ramp,ratio_frac,ratio
      real*8 temp_node,temp_low,temp_high
      real*8 DeltaPoros, DeltaHeight
      
c     
c     calculate components of volume strain
c     
      iispmd = ispm(i)     
     
c     spm1f is the youngs modulus at ambient conditions
c     spm2f is is the youngs modulus at frozen conditions
c     spm3f is temperture differnce for the ramp 
c    
c     
      if(icnl.ne.0) then   
c     2D x-y version  (y can be vertical)  
         ratio =  spm1f(iispmd)
         t_ramp = spm2f(iispmd)
         temp_low = -0.5*t_ramp
         temp_high = 0.5*t_ramp 
         temp_node = t(i) 
         ratio_frac =  (temp_node-temp_low)/t_ramp
c     
         ipchk = 0
           if(temp_node.lt.temp_low) then
            e1(i) = ratio*e10(i)
            e2(i) = ratio*e20(i)
            e3(i) = ratio*e30(i)
           else if(temp_node.gt.temp_high) then
            e1(i) = e10(i)
            e2(i) = e20(i)
            e3(i) = e30(i)
           else
            e1(i) = (1.-ratio_frac)*ratio*e10(i)
            e2(i) = (1.-ratio_frac)*ratio*e20(i)
            e3(i) = (1.-ratio_frac)*ratio*e30(i)              
           endif
  
      else if(icnl.eq.0) then   
c     
c     3D  version  (same as the 2D version)
c     
         ratio =  spm1f(iispmd)
         t_ramp = spm2f(iispmd)
         temp_low = -0.5*t_ramp
         temp_high = 0.5*t_ramp 
         temp_node = t(i) 
         ratio_frac =  (temp_node-temp_low)/t_ramp
c     
         ipchk = 0
           if(temp_node.lt.temp_low) then
            e1(i) = ratio*e10(i)
            e2(i) = ratio*e20(i)
            e3(i) = ratio*e30(i)
           else if(temp_node.gt.temp_high) then
            e1(i) = e10(i)
            e2(i) = e20(i)
            e3(i) = e30(i)
           else
            e1(i) = ratio_frac*ratio*e10(i)
            e2(i) = ratio_frac*ratio*e20(i)
            e3(i) = ratio_frac*ratio*e30(i)              
           endif
      endif
c
c get porosity
c  
c          ps(i) = DeltaPoros(t(i),psini(i)) 
c
c get height of cell
c   
          dzrg_new(i) = DeltaHeight(t(i),dzrg(i),psini(i))
      return

      end
c.....................................................................
