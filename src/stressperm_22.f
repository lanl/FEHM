      subroutine stressperm_22(jpt)
c mohr-coulomb failure criteria on the plane of maximum shear

      use comsi
      implicit none

      integer jpt,fail_flag
      real*8 excess_shear, shear_max,shear_angle,stress_norm
      real*8 eigenvec(3,3),alambda(3),rm(3,3)
      real*8 friction,strength, pp_fac

      fail_flag = 0
      call stressperm_22_failure(jpt,excess_shear,shear_angle
     &     ,fail_flag,rm)
      if(fail_flag.eq.1) then
         call stressperm_22_perm(jpt,excess_shear,rm)
         call stressperm_22_emod(jpt,excess_shear,rm)
      endif


      return
      end

c.....................................................................
      subroutine stressperm_22_failure(jpt,excess_shear,shear_angle
     &     ,fail_flag,rm)
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
      real*8 excess_shear, shear_max,shear_angle,stress_norm
      real*8 eigenvec(3,3),alambda(3),rm(3,3)
      real*8 friction,strength, pi, pp_fac,cossh,sinsh

      pi=3.1415926535

      iispmd = ispm(jpt) 
      friction = spm1f(iispmd)
      strength = spm2f(iispmd)
      pp_fac = spm3f(iispmd)

c find the principal stresses. eigenvec(j=1 to 3,k) contains the x,y,z
c components of the eigenvector for the k-th eigenvalue.  
      call principal_stress_3D(jpt,alambda,eigenvec)
c max shear stress,and normal stress on this plane (45deg to principal)
      call max_excess_shear(jpt,friction,strength,alambda,
     &     pp_fac,eigenvec,excess_shear,shear_angle )
      cossh=dcos(shear_angle)
      sinsh=dsin(shear_angle)
c choose z-prime=normal to the failure plane, and 
c y-prime=median principal stress. 
c rm(column3)= direction cos of (z-prime) etc.
      rm(1,3)=cossh*eigenvec(1,3) + sinsh*eigenvec(1,1)
      rm(2,3)=cossh*eigenvec(2,3) + sinsh*eigenvec(2,1)
      rm(3,3)=cossh*eigenvec(3,3) + sinsh*eigenvec(3,1)
      rm(1,2)=eigenvec(1,2)
      rm(2,2)=eigenvec(2,2)
      rm(3,2)=eigenvec(3,2)
      rm(1,1)=rm(2,3)*rm(3,2)-rm(3,3)*rm(2,2)
      rm(2,1)=rm(1,3)*rm(3,1)-rm(3,3)*rm(1,2)
      rm(3,1)=rm(1,3)*rm(2,2)-rm(2,3)*rm(1,2)

      if(excess_shear.gt.0.0d0) then
         fail_flag = 1
      endif

      return
      end

c.....................................................................

      subroutine stressperm_22_perm(jpt,excess_shear,rm)
c rotate back the tensor of multiplication factors from the system 
c of the fault plane to the original axis, multiply permeabilities
      use comdi
      use comsi

      implicit none
      integer jpt
      integer i,j,k,l,m,n, iispmd
      real*8 excess_shear, perm1,perm2,perm3
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
      sfac=abs(excess_shear/spm4f(iispmd))
      sfac=min(1.,sfac)
      fac_p(1,1)=1.0 + (spm8f(iispmd)-1)*sfac
      fac_p(2,2)=1.0 + (spm9f(iispmd)-1)*sfac
      fac_p(3,3)=1.0 + (spm10f(iispmd)-1)*sfac
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
      subroutine stressperm_22_emod(jpt,excess_shear,rm)

      use comsi
      implicit none

      integer jpt,iispmd
      real*8 excess_shear, rm(3,3)

      iispmd = ispm(jpt)    

      return
      end

c.....................................................................
      subroutine max_excess_shear(jpt,friction,strength,alambda,
     &     pp_fac,eigenvec,excess_shear,shear_angle )

      use comsi
      use comdi, only : phi
      implicit none

      integer jpt,iispmd
      real*8 excess_shear,shear_angle,shear_max ,stress_norm  
      real*8 eigenvec(3,3),alambda(3)
      real*8 friction,strength, pi, pp_fac

      pi=3.1415926535

c      iispmd = ispm(jpt) 
c      friction = spm1f(iispmd)
c      strength = spm2f(iispmd)

c Jaeger & Cook, pp 95-96
      shear_max = 0.5*(alambda(3)-alambda(1))
      stress_norm = 0.5*(alambda(3)+alambda(1))-pp_fac*phi(jpt)
c minimum value of (shear - friction*normal stress)
      excess_shear = shear_max*sqrt(1.0d0+friction*friction)
     &     - friction*stress_norm
      excess_shear = excess_shear - strength
      shear_angle = pi*0.5-0.5d0*datan(1.0d0/friction)

      return
      end

c.....................................................................
