      subroutine stressperm_21(jpt)
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

c mohr-coulomb failure criteria on a user-specified plane
c spm1f,spm2f,spm3f=direction cosins of the normal to fault plane
c spm4f:friction coefficient of shear in the fault plane
c spm5f:shear strength of the fault plane
c spm6f:factor in effective stress=sigma-(pp_fac*pore pressure)
c spm7f:range of excess shear stress over which damage is ramped
c spm8f:maximum multiplier for young's modulus  in x-prime direction
c spm9f:maximum multiplier for young's modulus  in y-prime direction
c spm10f:maximum multiplier for young's modulus  in z-prime direction
c spm11f: maximum multiplier  permeability x-prime direction
c spm12f: maximum multiplier  permeability y-prime direction
c spm13f: maximum multiplier  permeability z-prime direction
      use comsi
      implicit none
      integer jpt, fail_flag
      real*8 rm(3,3)
      fail_flag=0

      call stressperm_21_failure(jpt, fail_flag,rm)
      if(fail_flag.eq.1) then
         call stressperm_21_perm(jpt,rm)
         call stressperm_21_emod(jpt,rm)
      endif

      return
      end
c.......................................................................

      subroutine stressperm_21_failure(jpt, fail_flag,rm21)
c mohr-coulomb failure criteria on a user-specified plane
      use comai , only: iptty,ierr
      use comdi
      use comsi

      implicit none
      integer jpt,fail_flag
      integer i,j,k,l,m,n, iispmd
      real*8 sig_norm,sig_xp,sig_yp,sig_xyp,taumax
      real*8 xnorm,s_yf,s_xf,tau_xyf
      real*8 pi, shear_max_inplane, angle_min,angle_max
      real*8  perm1,perm2,perm3
      real*8 rm21(3,3), sig(3,3), sig21(3,3), sum21
      real*8 rm31s,rm33s,rm32s
      real*8 pp_fac

      pi=3.1415926535

      iispmd = ispm(jpt)
      pp_fac = spm6f(iispmd)

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
c Jaeger&Cook-pp25, eq 6. Make sure the direction cosins of the normal
c (specified in the input file) are normalized
      rm31s=spm1f(iispmd)*spm1f(iispmd)+spm2f(iispmd)*spm2f(iispmd)
     &     +spm3f(iispmd)*spm3f(iispmd)
      rm31s=sqrt(rm31s)
      if(rm31s.le.0.) then
         write(iptty,*)'error in input. check permmodel 21 dir cosins.'
         write(iptty,*)'STOP'
         write(ierr,*)'error in input. check permmodel 21 dir cosins.'
         write(ierr,*)'STOP'
         stop
      endif
      rm21(3,1)=spm1f(iispmd)/rm31s
      rm21(3,2)=spm2f(iispmd)/rm31s
      rm21(3,3)=spm3f(iispmd)/rm31s
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
               do l=1,3
                 sum21=sum21+rm21(i,l)*rm21(j,k)*sig(l,k)
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
      excess_shear(jpt)=abs(shear_max_inplane)-spm5f(iispmd)
      excess_shear(jpt)=excess_shear(jpt)-spm4f(iispmd)*(sig_norm)
      if(excess_shear(jpt).gt.0.)then
         fail_flag=1
      endif

      return
      end
c..........................................................................

      subroutine stressperm_21_perm(jpt,rm)
c rotate back the tensor of multiplication factors from the system 
c of the fault plane to the original axis, multiply permeabilities
      use comdi
      use comsi

      implicit none
      integer jpt
      integer i,j,k,l,m,n, iispmd
      real*8 perm1,perm2,perm3
      real*8 rm(3,3), rminv(3,3), sum21, sfac
      real*8 fac(3,3),fac_p(3,3)

c material failed, modify perm and young's modulus. fac1,fac2,fac3 are
c user specified multiplers, face1 and fac2 are in the plane of the
c fault, fac1 along the maximum shear, fac2 normal to fac1, and fac3
c normal to the fracture plane. These are treated as the principal
c values of a tensor and projected on x,y,z,axis to obtain multipliers
c for Kxx, Kyy, and Kzz. Crossterms in permeability are egnored for now
c we are rotating back, using the inverse of rm

      iispmd = ispm(jpt)    

      do i=1,3
         do j=1,3
            rminv(i,j)=rm(j,i)
         enddo
      enddo
      fac_p=0.
      fac=0.
      sfac=excess_shear(jpt)/spm7f(iispmd)
      sfac=min(1.d0,sfac)
      fac_p(1,1)=1.0 + (spm11f(iispmd)-1)*sfac
      fac_p(2,2)=1.0 + (spm12f(iispmd)-1)*sfac
      fac_p(3,3)=1.0 + (spm13f(iispmd)-1)*sfac
      do i=1,3
         do j=1,3
            sum21=0.
            do k=1,3
               do l=1,3
                  sum21=sum21+rminv(i,l)*rminv(j,k)*fac_p(l,k)
               enddo
            enddo
            fac(i,j)=sum21
         enddo
      enddo
      
      pnx(jpt)=pnx0(jpt)*fac(1,1)
      pny(jpt)=pny0(jpt)*fac(2,2)
      pnz(jpt)=pnz0(jpt)*fac(3,3)
      
      return
      end

c.....................................................................
      subroutine stressperm_21_emod(jpt,rm)

      use comsi
      implicit none

      integer jpt,iispmd
      real*8 rm(3,3)

      iispmd = ispm(jpt)    

      return
      end

c.....................................................................
