      subroutine stressperm_24(jpt)
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

c d dempsey, Feb 2013

c mohr-coulomb failure criteria applied to an ensemble of fractures with
c normal orientations evenly distributed about a stereonet. post-failure
c increase in permeability occurs according to an empirical model based
c on Lee and Cho (2002).

c INPUTS:
c Ks = shear fracture toughness (spm1f)
c mu_s = static friction coefficient (spm2f)
c mu_d = dynamic friction coefficient (spm3f)
c frac_num = number of fractures per control volume (spm4f)
c d1 = shear displacement at which perm enhancement begins (spm5f)
c Dd = shear displacement interval to complete perm enhancement (spm6f)
c dk = total perm enhancement (in log(perm)) (spm7f)
c cohesion = fracture cohesion (spm8f)
c frac_den = fracture density in control volume (spm9f)
      use comai
      use comdi
      use comsi
      implicit none
      
      integer i,jpt,fail_flag,iispmd,frac_num
      real*8 dc_l1,dc_m1,dc_n1
      real*8 dc_l2,dc_m2,dc_n2
      real*8 dc_l3,dc_m3,dc_n3
      real*8 theta(1000), sprl_phi(1000)
      real*8 str_g(3,3), str_l(3,3)
      real*8 perm_frac,frac_den,tau_ex,perm0
      real*8 perm_mult(3,3), perm_denom(3,3)

c get current and initial global stress components
      str_g(1,1) = str_x(jpt)
      str_g(2,2) = str_y(jpt)
      str_g(3,3) = str_z(jpt)
      str_g(1,2) = str_xy(jpt)
      str_g(1,3) = str_xz(jpt)
      str_g(2,3) = str_yz(jpt)
      str_g(2,1) = str_g(1,2)
      str_g(3,1) = str_g(1,3)
      str_g(3,2) = str_g(2,3)

c zero permeability modifying matrix
      perm_mult(1,1)=0.
      perm_mult(2,2)=0.
      perm_mult(3,3)=0.
      perm_mult(1,2)=0.
      perm_mult(1,3)=0.
      perm_mult(2,3)=0.
      perm_mult(2,1) = perm_mult(1,2)
      perm_mult(3,1) = perm_mult(1,3)
      perm_mult(3,2) = perm_mult(2,3)
      
	  perm_denom(1,1)=0.
      perm_denom(2,2)=0.
      perm_denom(3,3)=0.
      perm_denom(1,2)=0.
      perm_denom(1,3)=0.
      perm_denom(2,3)=0.
      perm_denom(2,1) = perm_denom(1,2)
      perm_denom(3,1) = perm_denom(1,3)
      perm_denom(3,2) = perm_denom(2,3)
	  
      iispmd = ispm(jpt)
      frac_num = spm4f(iispmd)
      frac_den = spm9f(iispmd)
      call stressperm_24_spiral(jpt,frac_num,theta,sprl_phi)
      do i=2,frac_num+1
        call stressperm_24_local(i,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2
     &       ,dc_l3,dc_m3,dc_n3,theta,sprl_phi,str_g)
        call stressperm_24_rotate(dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,
     &       dc_l3,dc_m3,dc_n3,str_g,str_l)
        perm_frac=0.
        call stressperm_24_failure(jpt,fail_flag,str_l)
        if(fail_flag.eq.1) then
	      call stressperm_24_postfailure(jpt,str_l,perm_frac)
        endif
        call stressperm_24_perm(jpt,perm_mult,perm_denom,perm_frac
     &       ,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2)
      enddo
c calculate 'per fracture' permeability and multiply by fracture density
      perm_mult(1,1) = (perm_mult(1,1)/perm_denom(1,1)-1)*frac_den+1
      perm_mult(2,2) = (perm_mult(2,2)/perm_denom(2,2)-1)*frac_den+1
      perm_mult(3,3) = (perm_mult(3,3)/perm_denom(3,3)-1)*frac_den+1
c      perm_mult(1,2) = (perm_mult(1,2)/perm_denom(1,2)-1)*frac_den+1
c      perm_mult(1,3) = (perm_mult(1,3)/perm_denom(1,3)-1)*frac_den+1
c      perm_mult(2,3) = (perm_mult(2,3)/perm_denom(2,3)-1)*frac_den+1
c      pnx(jpt)=pnx0(jpt)*perm_mult(1,1)
c      pny(jpt)=pny0(jpt)*perm_mult(2,2)
c      pnz(jpt)=pnz0(jpt)*perm_mult(3,3)
      perm0 = (pnx0(jpt)+pny0(jpt)+pnz(jpt))/2
      pnx(jpt)=perm_denom(1,1)*perm_mult(1,1)*perm0/frac_num
      pny(jpt)=perm_denom(2,2)*perm_mult(2,2)*perm0/frac_num
      pnz(jpt)=perm_denom(3,3)*perm_mult(3,3)*perm0/frac_num
c      write(*,'(a,F8.4,f8.4)') 'dkx,kx = ', perm_mult(1,1),pnx(jpt)
c      write(*,'(a,F8.4,f8.4)') 'dky,ky = ', perm_mult(2,2),pny(jpt)
c      write(*,'(a,F8.4,f8.4)') 'dkz,kz = ', perm_mult(3,3),pnz(jpt)
      return
      end
c........................................................................
      subroutine stressperm_24_spiral(jpt,frac_num,theta,sprl_phi)
c calculate a set of evenly distributed fracture orientations using the 
c Saff and Kuijlaars (1997) spiral walk algorithm.
      use comai
      use comdi
      use comsi
      implicit none
      
      integer i,jpt,iispmd,frac_num
      real*8 theta(1000),sprl_phi(1000),sprl_par(1000)
      real*8 pi,sprl_par1,sprl_par2,i2,hi,frac_num_real
      
      pi = 3.1415926536
	  
      iispmd = ispm(jpt)
c because we are excluding polar orientations, the fracture set should
c be supplemented by an additional two members in place.
      frac_num_real = frac_num+2.
      sprl_par1 = 1.-1./(frac_num_real-3.)
      sprl_par2 = 0.5*(frac_num_real+1.)/(frac_num_real-3.)
      sprl_par(1) = 0.
      theta(1) = pi
      sprl_phi(1) = 0.
      do i=2,frac_num+1
        i2 = sprl_par1*(i-1)+sprl_par2
        hi = 2.*(i2-1.)/(frac_num_real-1.)-1.
        sprl_par(i) = sqrt(1-hi*hi)
c calculate polar angles of next point on the walk
        theta(i) = dacos(hi)
        sprl_phi(i) = sprl_phi(i-1)+7.2/sqrt(frac_num_real)
     &/(sprl_par(i-1)+sprl_par(i))
      enddo
      return
      end
c........................................................................
      subroutine stressperm_24_local(i,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2
     &       ,dc_n2,dc_l3,dc_m3,dc_n3,theta,sprl_phi,str_g)
c for a given fracture orientation with normal corresponding to z-prime,
c calculate x- and y- prime directions. note, x-prime is aligned along 
c the fracture shear stress. shear stress along y-prime is zero
      use comai
      use comdi
      use comsi
      
      implicit none
      integer i
      real*8 dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,dc_l3,dc_m3,dc_n3,A,B
      real*8 theta(1000),sprl_phi(1000)
      real*8 str_g(3,3)

c calculate direction cosines of fracture normal from polar angles
      dc_l3 = dsin(sprl_phi(i))*dsin(theta(i))
      dc_m3 = dcos(sprl_phi(i))*dsin(theta(i))
      dc_n3 = dcos(theta(i))
c calculate direction cosines of fracture plane vector perpendicular
c to the resolved shear stress
      A = -dc_m3*(str_g(1,1)+dc_n3/dc_l3*str_g(1,3)+dc_m3/dc_l3
     &*str_g(1,2))+dc_m3*str_g(2,2)+dc_n3*str_g(2,3)+dc_l3*str_g(1,2)
      B = -dc_n3*(str_g(1,1)+dc_n3/dc_l3*str_g(1,3)+dc_m3/dc_l3
     &*str_g(1,2))+dc_n3*str_g(3,3)+dc_m3*str_g(2,3)+dc_l3*str_g(1,3)
      dc_m2 = 1./sqrt(((A/B*dc_n3-dc_m3)/dc_l3)**2+(A/B)**2+1.)
      dc_n2 = -A/B*dc_m2
      dc_l2 = (-dc_m2*dc_m3-dc_n2*dc_n3)/dc_l3
c calculate direction cosines of fracture plane vector parallel to 
c the resolved shear stress - cross product of other two basis vectors
      dc_l1 = dc_m2*dc_n3 - dc_m3*dc_n2
      dc_m1 = dc_n2*dc_l3 - dc_l2*dc_n3
      dc_n1 = dc_l2*dc_m3 - dc_m2*dc_l3
      return
      end
c........................................................................
      subroutine stressperm_24_rotate(dc_l1,dc_m1,dc_n1,dc_l2,dc_m2
     &       ,dc_n2,dc_l3,dc_m3,dc_n3,str_g,str_l)
c calculate stress components in the local coordinate system
      use comai
      use comdi
      use comsi
      implicit none
      
      real*8 dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,dc_l3,dc_m3,dc_n3
      real*8 str_g(3,3), str_l(3,3)
	  
c sigma_zz	  
      str_l(3,3)=dc_l3**2*str_g(1,1)+dc_m3**2*str_g(2,2)
     &+dc_n3**2*str_g(3,3)+2*dc_m3*dc_n3*str_g(2,3)
     &+2*dc_n3*dc_l3*str_g(1,3)+2*dc_l3*dc_m3*str_g(1,2)
c sigma_yy
c      str_l(2,2)=dc_l2**2*str_g(1,1)+dc_m2**2*str_g(2,2)
c     &+dc_n2**2*str_g(3,3)+2*dc_m2*dc_n2*str_g(2,3)
c     &+2*dc_n2*dc_l2*str_g(1,3)+2*dc_l2*dc_m2*str_g(1,2)
c sigma_xx
c      str_l(1,1)=dc_l1**2*str_g(1,1)+dc_m1**2*str_g(2,2)
c     &+dc_n1**2*str_g(3,3)+2*dc_m1*dc_n1*str_g(2,3)
c     &+2*dc_n1*dc_l1*str_g(1,3)+2*dc_l1*dc_m1*str_g(1,2)
c sigma_xy	  
c      str_l(1,2)=dc_l1*dc_l2*str_g(1,1)+dc_m1*dc_m2*str_g(2,2)
c     &+dc_n1*dc_n2*str_g(3,3)+(dc_m1*dc_n2+dc_m2*dc_n1)*str_g(2,3) 
c     &+(dc_n1*dc_l2+dc_n2*dc_l1)*str_g(1,3)
c     &+(dc_l1*dc_m2+dc_l2*dc_m1)*str_g(1,2)
c sigma_xz
      str_l(1,3)=dc_l1*dc_l3*str_g(1,1)+dc_m1*dc_m3*str_g(2,2)
     &+dc_n1*dc_n3*str_g(3,3)+(dc_m1*dc_n3+dc_m3*dc_n1)*str_g(2,3) 
     &+(dc_n1*dc_l3+dc_n3*dc_l1)*str_g(1,3)
     &+(dc_l1*dc_m3+dc_l3*dc_m1)*str_g(1,2)
c sigma_yz
c      str_l(2,3)=dc_l2*dc_l3*str_g(1,1)+dc_m2*dc_m3*str_g(2,2)
c     &+dc_n2*dc_n3*str_g(3,3)+(dc_m2*dc_n3+dc_m3*dc_n2)*str_g(2,3)
c     &+(dc_n2*dc_l3+dc_n3*dc_l2)*str_g(1,3)
c     &+(dc_l2*dc_m3+dc_l3*dc_m2)*str_g(1,2)
      return
      end
c........................................................................
      subroutine stressperm_24_failure(jpt,fail_flag,str_l)
c determine if Mohr-Coulomb failure has occurred on the fracture plane
      use comai
      use comdi
      use comsi
      implicit none
      
      integer jpt,fail_flag,iispmd
      real*8 str_l(3,3)
      real*8 tau_ex,mu_s,cohesion
      
      iispmd = ispm(jpt)
      mu_s = spm2f(iispmd)
      cohesion = spm8f(iispmd)
      tau_ex = abs(str_l(1,3))-mu_s*(str_l(3,3)-phi(jpt))-cohesion
      fail_flag = 0
      if (tau_ex.gt.0) then
        fail_flag = 1
      endif
      return
      end
c........................................................................
      subroutine stressperm_24_postfailure(jpt,str_l,perm_frac)
c calculate post-failure permeability enhancement according to the 
c Lee and Cho (2002) empirical model
      use comai
      use comdi
      use comsi
      implicit none
      
      integer jpt,iispmd
      real*8 str_l(3,3)
      real*8 perm_frac,mu_d,d1,Dd,dk,Ks,d_shear,l19,tau_ex,cohesion,mu_s
      iispmd = ispm(jpt)
      
      Ks = spm1f(iispmd)
      mu_s = spm2f(iispmd)
      mu_d = spm3f(iispmd)
      d1 = spm5f(iispmd)
      Dd = spm6f(iispmd)
      dk = spm7f(iispmd)
      cohesion = spm8f(iispmd)
	  
      tau_ex=abs(str_l(1,3))-mu_d*(str_l(3,3)-phi(jpt))
      d_shear = tau_ex/Ks
      l19 = log(19.)
      perm_frac = perm_frac+dk/(1.+exp(1.-2.*((d_shear-d1)/Dd))**l19)
      return
      end
c........................................................................
      subroutine stressperm_24_perm(jpt,perm_mult,perm_denom,perm_frac
     &           ,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2)
c accumulate fracture modified permeability components in global frame
      use comai
      use comdi
      use comsi
      implicit none
      
      integer jpt
      real*8 perm_frac,frac_perm,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2
      real*8 perm_mult(3,3), perm_denom(3,3), perm_term
	  
c diagonal terms	  
      perm_mult(1,1)=perm_mult(1,1)+(dc_l1**2+dc_l2**2)*10**perm_frac
      perm_mult(2,2)=perm_mult(2,2)+(dc_m1**2+dc_m2**2)*10**perm_frac
      perm_mult(3,3)=perm_mult(3,3)+(dc_n1**2+dc_n2**2)*10**perm_frac
c off-diagonal terms	  
c      perm_mult(1,2)=perm_mult(1,2)+
c     &(dc_l1*dc_m1+dc_l2*dc_m2)*10**perm_frac
c      perm_mult(1,3)=perm_mult(1,3)+
c     &(dc_l1*dc_n1+dc_l2*dc_n2)*10**perm_frac
c      perm_mult(2,3)=perm_mult(2,3)+
c     &(dc_m1*dc_n1+dc_m2*dc_n2)*10**perm_frac
c denominator	 
      perm_denom(1,1)=perm_denom(1,1)+(dc_l1**2+dc_l2**2)
      perm_denom(2,2)=perm_denom(2,2)+(dc_m1**2+dc_m2**2)
      perm_denom(3,3)=perm_denom(3,3)+(dc_n1**2+dc_n2**2)
c      perm_denom(1,2)=perm_denom(1,2)+(dc_l1*dc_m1+dc_l2*dc_m2)
c      perm_denom(1,3)=perm_denom(1,3)+(dc_l1*dc_n1+dc_l2*dc_n2)
c      perm_denom(2,3)=perm_denom(2,3)+(dc_m1*dc_n1+dc_m2*dc_n2)
	  
      return
      end
c........................................................................
