      subroutine stressperm_23(jpt)
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

c using the hydraulic aperture along the direction of shear movement
c a provided by Stuart Walsh from LANL, 6 Nov 2012.
c direction of shear slip assumed to be vertical. mohr-coulomb failure 
c criteria on the vertical plane 

      use comai
      use comdi
      use comsi
      implicit none

      integer jpt,fail_flag,iispmd
      real*8 rm(3,3)
      real*8 frac_perm

      fail_flag = 0
      frac_perm = 1.e-18

      iispmd = ispm(jpt) 
      call stressperm_23_failure(jpt,fail_flag,rm)
      if(fail_flag.eq.1) then
         if(itemp_perm23(jpt).eq.0) then
            write(91,*)jpt
            itemp_perm23(jpt) = 1
         endif
         call stressperm_23_perm(jpt,frac_perm)
c fehm normalizes permeabilities by 10^6
         frac_perm = frac_perm*1.e+6
c         if(frac_perm.lt.pnz(jpt)) then
c            write(*,*) l, jpt,frac_perm,pnz(jpt) 
c         endif
         pnz(jpt)=frac_perm

      endif
      

      return
      end

c.......................................................................

      subroutine stressperm_23_failure(jpt, fail_flag,rm21)
c mohr-coulomb failure criteria on a verticalplane
      use comai , only: iptty,ierr
      use comdi
      use comsi

      implicit none
      integer jpt,fail_flag
      integer i,j,k,ijk, iispmd
      real*8 sig_norm,sig_xp,sig_yp,sig_xyp,taumax
      real*8 xnorm,s_yf,s_xf,tau_xyf
      real*8 pi, shear_max_inplane, angle_min,angle_max
      real*8  perm1,perm2,perm3
      real*8 rm21(3,3), sig(3,3), sig21(3,3), sum21
      real*8 rm31s,rm33s,rm32s
      real*8 pp_fac

      pi=3.1415926535

      iispmd = ispm(jpt)
      pp_fac = spm3f(iispmd)

      rm21=0.
      sig=0.
      sig21=0.

c assemble the symmetric streess matrix
      sig(1,1)=str_x(jpt)
      sig(2,2)=str_y(jpt)
      sig(3,3)=str_z(jpt)
      sig(1,2)=str_xy(jpt)
      sig(1,3)=str_xz(jpt)
      sig(2,3)=str_yz(jpt)
      sig(2,1)=sig(1,2)
      sig(3,1)=sig(1,3)
      sig(3,2)=sig(2,3)

c Normal stress accross the specified plane. rm21(i,j) is the rotation matrix
c take z-prime normal to the plane. rminv=inverse(rm21)=transpose(rm21)
c Jaeger&Cook-pp25, eq 6.  direction cosins of the normal to the vertical 
c are (lambda,mue,0) with lambda^2+mue^2 = 1. assume mue = 0.
      rm31s=1.
      rm31s=sqrt(rm31s)
      rm21(3,1)= 1.
      rm21(3,2)= 0.
      rm21(3,3)= 0.
c within the plane x-prime and y-prime can be rotated, for simplicity
c we choose , notating z-prime=(l,m,n), x-prime=(m, -l,0)/normx and
c y-prime=(l*n,m*n,-normx)/normx, where normx=sqrt(l*l+m*m)
      xnorm=sqrt(rm21(3,1)*rm21(3,1)+rm21(3,2)*rm21(3,2))
      if(xnorm.gt.0.) then
         rm21(1,1)=+rm21(3,2)/xnorm
         rm21(1,2)=-rm21(3,1)/xnorm
         rm21(2,1)=+rm21(3,1)*rm21(3,3)/xnorm
         rm21(2,2)=+rm21(3,2)*rm21(3,3)/xnorm
         rm21(2,3)=-xnorm
      else
         rm21(1,1)=1.0
         rm21(2,2)=1.
      endif
c stresses in the coordinate system of the fault plane
      do i=1,3
         do j=1,3
            sum21=0.
            do k=1,3
               do ijk=1,3
                 sum21=sum21+rm21(i,ijk)*rm21(j,k)*sig(ijk,k)
               enddo
            enddo
            sig21(i,j)=sum21            
         enddo
      enddo
c effective stress normal to the plane
      sig_norm=sig21(3,3)-pp_fac*phi(jpt)
c stresses within the plane of the fault, in the new coordinates
      s_xf=sig21(1,1)
      s_yf=sig21(2,2)
      tau_xyf=sig21(1,2)
c calculate in-plane principal values and angle of maximum shear
      shear_max_inplane=sqrt(tau_xyf*tau_xyf+0.25*(s_xf-s_yf)**2.)
      angle_min=datan(2*tau_xyf/(s_xf-s_yf))
      angle_max=angle_min+pi/4.
c check for failure using mohr-coulomb, effective stress
      excess_shear(jpt)=abs(shear_max_inplane)-spm2f(iispmd)
      excess_shear(jpt)=excess_shear(jpt)-spm1f(iispmd)*(sig_norm)
      if(excess_shear(jpt).gt.0.)then
         fail_flag=1
      endif

      return
      end
c..........................................................................

      subroutine stressperm_23_perm(jpt,frac_perm)
c using the hydraulic aperture along the direction of shear movement
c a provided by Stuart Walsh from LANL, 6 Nov 2012. for now, assuming
c shear movement to be in the vertical
      use comdi
      use comsi

      implicit none
      integer jpt
      integer i,j,k,l,m,n, iispmd
      real*8  perm1,perm2,perm3
      real*8 shear_mod, shear_disp, sfac
      real*8 haperture,frac_perm   

      iispmd = ispm(jpt)    

      shear_mod = spm4f(iispmd)
      shear_disp = excess_shear(jpt)/shear_mod
      sfac=min (1.d0, (shear_disp/spm5f(iispmd)))
      haperture = max(1.d-6, abs(spm6f(iispmd)*sfac)) 
      frac_perm = max(1.d-18,(haperture**3)/(12.0*spm7f(iispmd)))
      
      return
      end

c.....................................................................
