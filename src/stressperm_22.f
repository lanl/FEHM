      subroutine stressperm_22(jpt)
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

c mohr-coulomb failure criteria on the plane of maximum shear

      use comdi
      use comsi
      implicit none

      integer jpt,fail_flag,iispmd
      real*8  shear_max,stress_norm
c gaz 052017      
c      real*8 eigenvec(3,3),alambda(3),rm(3,3)
      real*8 rm(3,3)
      real*8 friction,strength
      real*8 fac(3,3),fac_E(3), fac_por
      real*8 porosity_damage_factor

      fail_flag = 0
      fac = 1.
      fac_E = 1.
      fac_por=1.
      iispmd = ispm(jpt) 
      call stressperm_22_failure(jpt,fail_flag,rm)

      if(fail_flag.eq.1) then
         if(itemp_perm22(jpt).eq.0) then
            write(91,*)jpt
            itemp_perm22(jpt) = 1
         endif
         call stressperm_22_perm(jpt,rm, fac, fac_por)
         call stressperm_22_emod(jpt,rm,fac_E)
      endif
      
         pnx(jpt)=pnx0(jpt)*fac(1,1)
         pny(jpt)=pny0(jpt)*fac(2,2)
         pnz(jpt)=pnz0(jpt)*fac(3,3)
      
         porosity_damage_factor = 1.
         if(spm11f(iispmd).gt.0.) 
     &     porosity_damage_factor=1.0 + (spm11f(iispmd)-1.)*fac_por
         ps(jpt) = min(1.d0, porosity_damage_factor * psini(jpt))
      
         e1(jpt)=e10(jpt)*fac_E(1)
         e2(jpt)=e20(jpt)*fac_E(1)
         e3(jpt)=e30(jpt)*fac_E(1)
      

      return
      end

c.....................................................................
      subroutine stressperm_22_failure(jpt,fail_flag,rm)
c mohr-coulomb failure criteria on the plane that miximizes
c abs(shear)-friction*normal stress
c spm1f:friction coefficient of shear in the fault plane
c spm2f:shear strength of the fault plane
c spm3f:factor in effective stress=sigma-(pp_fac*pore pressure)
c spm4f:range of excess shear stress over which damage is ramped
c spm5f:maximum multiplier for young's modulus  in x-prime direction
c spm6f:maximum multiplier for young's modulus  in y-prime direction
c spm7f:maximum multiplier for young's modulus  in z-prime direction
c spm8f: maximum multiplier  permeability x-prime direction
c spm9f: maximum multiplier  permeability y-prime direction
c spm10f: maximum multiplier  permeability z-prime direction
c  here z-prime is along tehnormal to the plane of failure, and
c  y-prime is along the median principal stress

      use comsi
      implicit none

      integer jpt,iispmd,fail_flag
      real*8  shear_max,stress_norm
c gaz 052017      
c      real*8 eigenvec(3,3),alambda(3),rm(3,3)
      real*8 rm(3,3)
      real*8 friction,strength, pi, pp_fac, cossh, sinsh
      real*8 stress_factor_initial, shear_threshold

      pi=3.1415926535

      iispmd = ispm(jpt) 
      friction = spm1f(iispmd)
      strength = spm2f(iispmd)
      pp_fac = spm3f(iispmd)
      stress_factor_initial = 0.
      shear_threshold = 0.
      if(incremental_shear_permmodel.eq.1) then
         stress_factor_initial = spm12f(iispmd)
         shear_threshold = spm13f(iispmd)
      endif

c find the principal stresses. eigenvec(j=1 to 3,k) contains the x,y,z
c components of the eigenvector for the k-th eigenvalue.  
c gaz 052017
      if(incremental_shear_permmodel.eq.1) then
c         call principal_incremental_stress_3D(jpt,alambda,eigenvec
c     &        ,stress_factor_initial)          
         call principal_incremental_stress_3D(jpt,
     &        stress_factor_initial)
      else
c gaz 052017                   
c         call principal_stress_3D(jpt,alambda,eigenvec)
                  call principal_stress_3D(jpt)          
      endif
c max shear stress,and normal stress on this plane (45deg to principal)
c gaz 052017      
c      call max_excess_shear(jpt,friction,strength,alambda,
c     &     pp_fac,stress_factor_initial, eigenvec)
      call max_excess_shear(jpt,friction,strength,
     &     pp_fac,stress_factor_initial)      
      cossh=dcos(shear_angle(jpt))
      sinsh=dsin(shear_angle(jpt))
c choose z-prime=normal to the failure plane, and 
c y-prime=median principal stress. 
c rm(column3)= direction cos of (z-prime) etc.
      rm(1,3)=cossh*eigenvec(1,3) + sinsh*eigenvec(1,1)
      rm(2,3)=cossh*eigenvec(2,3) + sinsh*eigenvec(2,1)
      rm(3,3)=cossh*eigenvec(3,3) + sinsh*eigenvec(3,1)
      rm(1,2)=eigenvec(1,2)
      rm(2,2)=eigenvec(2,2)
      rm(3,2)=eigenvec(3,2)
      rm(1,1)= rm(2,3)*rm(3,2)-rm(3,3)*rm(2,2)
      rm(2,1)=-rm(1,3)*rm(3,2)+rm(3,3)*rm(1,2)
      rm(3,1)= rm(1,3)*rm(2,2)-rm(2,3)*rm(1,2)

      if(excess_shear(jpt).gt.shear_threshold) then
         fail_flag = 1
      endif

      return
      end

c.....................................................................

      subroutine stressperm_22_perm(jpt,rm, fac,fac_por)
c rotate back the tensor of multiplication factors from the system 
c of the fault plane to the original axis, multiply permeabilities
      use comdi
      use comsi

      implicit none
      integer jpt
      integer i,j,k,l,m,n, iispmd
      real*8  perm1,perm2,perm3
      real*8 rm(3,3), rminv(3,3), sum21, sfac
      real*8 fac(3,3),fac_p(3,3), fac_por

c material failed, modify perm and young's modulus. fac1,fac2,fac3 are
c user specified multiplers, face1 and fac2 are in the plane of the
c fault, fac1 along the maximum shear, fac2 normal to fac1, and fac3
c normal to the fracture plane. These are treated as the principal
c values of a tensor and projected on x,y,z,axis to obtain multipliers
c for Kxx, Kyy, and Kzz. Crossterms in permeability are egnored for now
c we are rotating back, using the inverse of rm

      fac_p=0.
      iispmd = ispm(jpt)    

      do i=1,3
         do j=1,3
            rminv(i,j)=rm(j,i)
         enddo
      enddo
      sfac=abs(excess_shear(jpt)/spm4f(iispmd))
      sfac=min(1.d0,sfac)
      fac_por=sfac
      fac_p(1,1)=1.0 + (spm8f(iispmd)-1.)*sfac
      fac_p(2,2)=1.0 + (spm9f(iispmd)-1.)*sfac
      fac_p(3,3)=1.0 + (spm10f(iispmd)-1.)*sfac
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
      
      return
      end

c.....................................................................
      subroutine stressperm_22_emod(jpt,rm,fac_E)

      use comai, only: iptty,ierr
      use comsi
      implicit none

      integer jpt,iispmd
      real*8 fac_E(3), sfac, rm(3,3)

      iispmd = ispm(jpt)    

c s kelkar Aug 11 2011. For now, simply multiplying youngs mod in every
c direction, keeping poissons ratio fixed and
c not taking care of the tensor nature of D
      if(spm5f(iispmd).le.0..or.spm6f(iispmd).le.0.
     &     .or.spm7f(iispmd).le.0.) then
         write(iptty,*)'Error in stressperm_22_emod. STOP'
         write(iptty,*)'spm(4,5, and 6)f must be',
     &        'gt 0. in the strs-permmodel=22'
         write(ierr,*)'Error in stressperm_22_emod. STOP'
         write(ierr,*)'spm(4,5,and 6)f must be',
     &        'gt 0. in the strs-permmodel=22'
         stop
      endif

      sfac=abs(excess_shear(jpt)/spm4f(iispmd))
      sfac=min(1.d0,sfac)
      fac_E(1)=1.0 - sfac + sfac/spm5f(iispmd)
      fac_E(2)=1.0 - sfac + sfac/spm6f(iispmd)
      fac_E(3)=1.0 - sfac + sfac/spm7f(iispmd)

      return
      end

c.....................................................................
c      subroutine max_excess_shear(jpt,friction,strength,alambda,
c     &     pp_fac,stress_factor_initial, eigenvec)
      subroutine max_excess_shear(jpt,friction,strength,
     &     pp_fac,stress_factor_initial)

      use comsi
      use comdi
      implicit none

      integer jpt,iispmd
      real*8 shear_max ,stress_norm  
c gaz 052017      
c      real*8 eigenvec(3,3),alambda(3)
      real*8 friction,strength, pi, pp_fac
      real*8 stress_factor_initial

      pi=3.1415926535

c      iispmd = ispm(jpt) 
c      friction = spm1f(iispmd)
c      strength = spm2f(iispmd)

c Jaeger & Cook, pp 95-96
      shear_max = 0.5*(alambda(3)-alambda(1))
      stress_norm = 0.5*(alambda(3)+alambda(1))
     &     -pp_fac*(phi(jpt)-stress_factor_initial*pho(jpt))
c minimum value of (shear - friction*normal stress)
      excess_shear(jpt) = shear_max*sqrt(1.0d0+friction*friction)
     &     - friction*stress_norm
      excess_shear(jpt) = excess_shear(jpt) - strength
      shear_angle(jpt) = pi*0.5-0.5*atan(1./friction)

      return
      end

c.....................................................................
c gaz 052017
c      subroutine principal_incremental_stress_3D(i,alambda,eigenvec
c     &        ,stress_factor_initial)
      subroutine principal_incremental_stress_3D(i,
     &        stress_factor_initial)      
      use comai
      use comsi

      implicit none
      integer i
c gaz 052017      
c      real*8 AMAT(3,3), eigenvec(3,3),alambda(3)
      real*8 AMAT(3,3)
      real*8 AI1,AI2,AI3
      real*8 stress_factor_initial
      real*8 str_init_temp(3,3)

      AMAT=0.0
      str_init_temp = 0.0

      if (allocated(str_x0_perm)) then
        str_init_temp(1,1) = stress_factor_initial*str_x0_perm(i) 
        str_init_temp(2,2) = stress_factor_initial*str_y0_perm(i) 
        str_init_temp(3,3) = stress_factor_initial*str_z0_perm(i) 
        str_init_temp(1,2) = stress_factor_initial*str_xy0_perm(i) 
        str_init_temp(1,3) = stress_factor_initial*str_xz0_perm(i) 
        str_init_temp(2,3) = stress_factor_initial*str_yz0_perm(i)
        str_init_temp(2,1) = str_init_temp(1,2)
        str_init_temp(3,1) = str_init_temp(1,3)
        str_init_temp(3,2) = str_init_temp(2,3)
      endif

      AMAT(1,1) = str_x(i) 
      AMAT(2,2) = str_y(i) 
      AMAT(3,3) = str_z(i) 
      AMAT(1,2) = str_xy(i) 
      AMAT(1,3) = str_xz(i) 
      AMAT(2,3) = str_yz(i) 
      AMAT(2,1) = AMAT(1,2)
      AMAT(3,2) = AMAT(2,3)
      AMAT(1,3) = AMAT(3,1)

      AMAT = AMAT - str_init_temp
      
! Compute the invariants for the matrix
      AI1 = AMAT(1,1) + AMAT(2,2) + AMAT(3,3)
      AI2 = AMAT(1,1)*AMAT(2,2) + AMAT(2,2)*AMAT(3,3) 
      AI2 = AI2 + AMAT(3,3)*AMAT(1,1) - AMAT(1,2)*AMAT(1,2) 
     &     - AMAT(2,3)*AMAT(2,3)  - AMAT(3,1)*AMAT(3,1)
      AI3 = AMAT(1,1)*(AMAT(2,2)*AMAT(3,3) - AMAT(2,3)*AMAT(3,2))
      AI3 = AI3 + AMAT(1,2)*(AMAT(2,3)*AMAT(3,1) - AMAT(2,1)*AMAT(3,3))
      AI3 = AI3 + AMAT(1,3)*(AMAT(2,1)*AMAT(3,2) - AMAT(3,1)*AMAT(2,2))
      
! Solve the characteristic equation for the eigenvalues
      call solve_cubic(-AI1, AI2,-AI3,alambda)
      
! Make the matrix tridiagonal using Householder transformation
      call reduce_to_tridiagonal(AMAT)
      
! Solve for the eigenvectors
      call eigenvectors(AMAT, alambda, eigenvec)
      
      return

      end

c.......................................................................
