      subroutine stressperm_25(jpt)
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
c orientations distributed according to fracture orienation data. 
c pre-failure increase in permeability occurs according to the 
c Bai et al. (1999) model. post-failure increase in permeability occurs
c according to an empirical model based on Lee and Cho (2002).

c INPUTS:
c Ks = shear fracture toughness (spm1f)
c mu_s = static friction coefficient (spm2f)
c mu_d = dynamic friction coefficient (spm3f)
c frac_num = number of fractures per control volume (spm4f)
c d1 = shear displacement at which perm enhancement begins (spm5f)
c Dd = shear displacement interval to complete perm enhancement (spm6f)
c dk = total perm enhancement (in log(perm)) (spm7f)
c cohesion = fracture cohesion (spm8f)
      use comai
      use comdi
      use comsi
      implicit none
      
      integer i,jpt,fail_flag,iispmd,frac_num
      real*8 dc_l1,dc_m1,dc_n1
      real*8 dc_l2,dc_m2,dc_n2
      real*8 dc_l3,dc_m3,dc_n3
      real*8 sprl_theta(1000), sprl_phi(1000)
      real*8 str_g(3,3), str_l(3,3)
      real*8 perm_frac,frac_den,tau_ex
      real*8 perm_mult(3,3), perm_denom(3,3)
      real*8 pf0,perm0
	  	  
      pf0 = str_pf0_perm(jpt)

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
      call stressperm_25_random(jpt)
      do i=1,frac_num
        call stressperm_25_local(jpt,i,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2
     &       ,dc_n2,dc_l3,dc_m3,dc_n3,str_g)
        call stressperm_25_rotate(dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,
     &       dc_l3,dc_m3,dc_n3,str_g,str_l)
        perm_frac=0.
        call stressperm_25_failure(jpt,fail_flag,str_l)
        if(fail_flag.eq.1) then
	      call stressperm_25_postfailure(jpt,str_l,perm_frac)
        endif
        call stressperm_25_perm(jpt,perm_mult,perm_denom,perm_frac
     &       ,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2)
      enddo
c calculate 'per fracture' permeability and multiply by fracture density
      perm_mult(1,1) = (perm_mult(1,1)/perm_denom(1,1)-1)+1
      perm_mult(2,2) = (perm_mult(2,2)/perm_denom(2,2)-1)+1
      perm_mult(3,3) = (perm_mult(3,3)/perm_denom(3,3)-1)+1
      perm_mult(1,2) = (perm_mult(1,2)/perm_denom(1,2)-1)+1
      perm_mult(1,3) = (perm_mult(1,3)/perm_denom(1,3)-1)+1
      perm_mult(2,3) = (perm_mult(2,3)/perm_denom(2,3)-1)+1
	  
      perm0 = (pnx0(jpt)+pny0(jpt)+pnz0(jpt))/2
      if(perm_mult1(jpt).lt.perm_mult(1,1)) then
         perm_mult1(jpt)=perm_mult(1,1)
      endif
      if(perm_mult2(jpt).lt.perm_mult(2,2)) then
         perm_mult2(jpt)=perm_mult(2,2)
      endif
      if(perm_mult3(jpt).lt.perm_mult(3,3)) then
         perm_mult3(jpt)=perm_mult(3,3)
      endif

      pnx(jpt)=perm_denom(1,1)*perm_mult1(jpt)*perm0/frac_num
      pny(jpt)=perm_denom(2,2)*perm_mult2(jpt)*perm0/frac_num
      pnz(jpt)=perm_denom(3,3)*perm_mult3(jpt)*perm0/frac_num
	  
      return
      end
c........................................................................
      subroutine stressperm_25_random(jpt)
c calculate a random set of fracture orienations based on supplied 
c fracture data. synthetic data has the same orienations density profile
c as field data.
      use comai
      use comdi
      use comsi
c      use compart
      implicit none
      
      integer i,j,jpt,iispmd,N_sectors,N_allocated,str25_N_syn,rseed
      integer N_allocated_orig, i_correct
      real ran_sp
      real*8 sprl_theta(1000),sprl_phi(1000),sprl_par(1000)
      real*8 pi,sprl_par1,sprl_par2,i2,hi,frac_num_real
      real*8 dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,mid_x,mid_y,mid_z
      real*8 dir_x,dir_y,dir_z,mag_dir,tan_x,tan_y,tan_z,tan_r
      real*8 dc_l3(1000),dc_m3(1000),dc_n3(1000)
      real*8 R_sector, N_density, rand_num
      integer I_allocated(1500)

      if (str25_xfrac(jpt,1).lt.1.) return
      pi = 3.1415926536
	  
      rseed=466201
      iispmd = ispm(jpt)
c number of synthetic fracture orientations to generate
      str25_N_syn = spm4f(iispmd)
      if (str25_N_syn.gt.1000) then
        write(*,'(A,I5,A)') 'Too many synthetic fracture orienations (',
     &  str25_N_syn, '), maximum 1000.'
        stop
      endif
c ************ decompose stereonet into sectors *****************
c decompose stereonet into sectors, defined by mid-points
      N_sectors = 200
      R_sector = sqrt(10./float(N_sectors))
c because we are excluding polar orientations, the fracture set should
c be supplemented by an additional two members in place.
      frac_num_real = N_sectors+2.
      sprl_par1 = 1.-1./(frac_num_real-3.)
      sprl_par2 = 0.5*(frac_num_real+1.)/(frac_num_real-3.)
      sprl_par(1) = 0.
      sprl_theta(1) = pi
      sprl_phi(1) = 0.
      do i=2,N_sectors+1
        i2 = sprl_par1*(i-1)+sprl_par2
        hi = 2.*(i2-1.)/(frac_num_real-1.)-1.
        sprl_par(i) = sqrt(1-hi*hi)
c calculate polar angles of next point on the walk
        sprl_theta(i) = dacos(hi)
        sprl_phi(i) = sprl_phi(i-1)+7.2/sqrt(frac_num_real)
     &/(sprl_par(i-1)+sprl_par(i))
c        write(*,'(a,f8.4,f8.4)') 'theta,phi=',sprl_theta(i),sprl_phi(i)
      enddo
c     calculate and save sector mid-points
      do i=2,N_sectors+1
c       sector mid-point
        dc_l3(i) = dsin(sprl_phi(i))*dsin(sprl_theta(i))
        dc_m3(i) = dcos(sprl_phi(i))*dsin(sprl_theta(i))
        dc_n3(i) = dcos(sprl_theta(i))
c        write(*,'(a,f8.4,f8.4,f8.4)') 'mid_x,y,z=',dc_l3(i),
c     &dc_m3(i),dc_n3(i)
      enddo
c ************ calculate fracture density data *****************
      if (str25_density(1).eq.0.) then
c       assign non-zero background fracture density
        do j=2,N_sectors+1
          str25_density(j)=0.05
        enddo
c       for each observation
        do j=1,str25_N_obs  
c         direction cosines of observed fracture pole
          dc_l2 = dsin(str25_azi(j)/180.*pi)*dsin(str25_dip(j)/180.*pi)
          dc_m2 = dcos(str25_azi(j)/180.*pi)*dsin(str25_dip(j)/180.*pi)
          dc_n2 = dcos(str25_dip(j)/180.*pi)
c         for each sector
          do i=2,N_sectors+1
c           if angular distance between fracture pole and sector mid-point
c           less than tolerance, add to sector
            if (acos(dc_l2*dc_l3(i)+dc_m2*dc_m3(i)+dc_n2*dc_n3(i))
     &          <R_sector) then 
              str25_density(i) = str25_density(i) + 1.
            endif
          enddo
        enddo
c       convert density vector to pdf
        N_density = 0.
        do j=2,N_sectors+1
          N_density = N_density + str25_density(j)
        enddo
        do j=2,N_sectors+1
          str25_density(j) = str25_density(j)/N_density
        enddo
      endif
c ************ calculate synthetic fracture data *****************
c distribute synthetic fractures among sectors
      N_allocated = 0
c     while not all fractures have been allocated
101   if (N_allocated.lt.str25_N_syn) then
c       for each sector
        do i=2,N_sectors+1
c         test pdf to see if fracture in sector
          rand_num = ran_sp(rseed)
c          write(*,'(f8.4)') rand_num
          if (rand_num<str25_density(i)) then
            N_allocated = N_allocated+1
            I_allocated(N_allocated)=i
c            write(*,'(i4)') i
          endif
        enddo
        goto 101
      endif
      N_allocated_orig = N_allocated
c     check if too many fractures have been allocated
102   if (N_allocated.gt.str25_N_syn) then
        I_allocated(int(ran_sp(rseed)*float(N_allocated-2)+1.))=0
        N_allocated = N_allocated - 1
        goto 102
      endif
c     for each synthetic fracture, calculate its position within the sector
      i_correct = 0
      do i=1,N_allocated_orig
c       sector midpoint
         if(I_allocated(i).gt.0) then

            i_correct = i_correct + 1

            mid_x=dc_l3(I_allocated(i))
            mid_y=dc_m3(I_allocated(i))
            mid_z=dc_n3(I_allocated(i))
c     random direction
            dir_x=2.*ran_sp(rseed)-1.
            dir_y=2.*ran_sp(rseed)-1.
            dir_z=2.*ran_sp(rseed)-1.
c     normalise direction
            mag_dir = sqrt(dir_x**2+dir_y**2+dir_z**2)
            dir_x=dir_x/mag_dir
            dir_y=dir_y/mag_dir
            dir_z=dir_z/mag_dir
c     calculate cross product, tangent to stereonet hemisphere
            tan_x = dir_y*mid_z - mid_y*dir_z
            tan_y = dir_z*mid_x - dir_x*mid_z
            tan_z = dir_x*mid_y - dir_y*mid_x
c     write(*,'(F8.4,F8.4,F8.4)') tan_x,tan_y,tan_z
c     calculate fracture pole
            tan_r = ran_sp(rseed)*R_sector
            dir_x = mid_x + tan_r*tan_x
            dir_y = mid_y + tan_r*tan_y
            dir_z = mid_z + tan_r*tan_z		
c     normalise pole, save
            mag_dir = sqrt(dir_x**2+dir_y**2+dir_z**2)
            str25_xfrac(jpt,i_correct)=dir_x/mag_dir
            str25_yfrac(jpt,i_correct)=dir_y/mag_dir
            str25_zfrac(jpt,i_correct)=dir_z/mag_dir
c     write(*,'(F8.4,F8.4,F8.4)') str25_xfrac(jpt,i),
c     &  str25_yfrac(jpt,i),str25_zfrac(jpt,i)
         endif
      enddo
      return
      end
c........................................................................
      subroutine stressperm_25_local(jpt,i,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2
     &       ,dc_n2,dc_l3,dc_m3,dc_n3,str_g)
c for a given fracture orientation with normal corresponding to z-prime,
c calculate x- and y- prime directions. note, x-prime is aligned along 
c the fracture shear stress. shear stress along y-prime is zero
      use comai
      use comdi
      use comsi
      
      implicit none
      integer jpt,i
      real*8 dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2,dc_l3,dc_m3,dc_n3,A,B
      real*8 str_g(3,3)

c calculate direction cosines of fracture normal from polar angles
      dc_l3 = str25_xfrac(jpt,i)
      dc_m3 = str25_yfrac(jpt,i)
      dc_n3 = str25_zfrac(jpt,i)
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
      subroutine stressperm_25_rotate(dc_l1,dc_m1,dc_n1,dc_l2,dc_m2
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
      str_l(2,2)=dc_l2**2*str_g(1,1)+dc_m2**2*str_g(2,2)
     &+dc_n2**2*str_g(3,3)+2*dc_m2*dc_n2*str_g(2,3)
     &+2*dc_n2*dc_l2*str_g(1,3)+2*dc_l2*dc_m2*str_g(1,2)
c sigma_xx
      str_l(1,1)=dc_l1**2*str_g(1,1)+dc_m1**2*str_g(2,2)
     &+dc_n1**2*str_g(3,3)+2*dc_m1*dc_n1*str_g(2,3)
     &+2*dc_n1*dc_l1*str_g(1,3)+2*dc_l1*dc_m1*str_g(1,2)
c sigma_xy	  
      str_l(1,2)=dc_l1*dc_l2*str_g(1,1)+dc_m1*dc_m2*str_g(2,2)
     &+dc_n1*dc_n2*str_g(3,3)+(dc_m1*dc_n2+dc_m2*dc_n1)*str_g(2,3) 
     &+(dc_n1*dc_l2+dc_n2*dc_l1)*str_g(1,3)
     &+(dc_l1*dc_m2+dc_l2*dc_m1)*str_g(1,2)
c sigma_xz
      str_l(1,3)=dc_l1*dc_l3*str_g(1,1)+dc_m1*dc_m3*str_g(2,2)
     &+dc_n1*dc_n3*str_g(3,3)+(dc_m1*dc_n3+dc_m3*dc_n1)*str_g(2,3) 
     &+(dc_n1*dc_l3+dc_n3*dc_l1)*str_g(1,3)
     &+(dc_l1*dc_m3+dc_l3*dc_m1)*str_g(1,2)
c sigma_yz
      str_l(2,3)=dc_l2*dc_l3*str_g(1,1)+dc_m2*dc_m3*str_g(2,2)
     &+dc_n2*dc_n3*str_g(3,3)+(dc_m2*dc_n3+dc_m3*dc_n2)*str_g(2,3)
     &+(dc_n2*dc_l3+dc_n3*dc_l2)*str_g(1,3)
     &+(dc_l2*dc_m3+dc_l3*dc_m2)*str_g(1,2)

      return
      end
c........................................................................
      subroutine stressperm_25_failure(jpt,fail_flag,str_l)
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
       subroutine stressperm_25_postfailure(jpt,str_l,perm_frac)
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
      if (cohesion.ge.0) then
        perm_frac = perm_frac+dk/(1.+exp(1.-2.*((d_shear-d1)/Dd))**l19)
      endif
	
      return
      end
c........................................................................
      subroutine stressperm_25_perm(jpt,perm_mult,perm_denom,perm_frac
     &           ,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2)
c accumulate fracture modified permeability components in global frame
      use comai
      use comdi
      use comsi
      implicit none
      
      integer jpt,iispmd
      real*8 perm_frac,frac_perm,dc_l1,dc_m1,dc_n1,dc_l2,dc_m2,dc_n2
      real*8 perm_mult(3,3), perm_denom(3,3), perm_term
	  
      iispmd = ispm(jpt)
	  
      perm_mult(1,1)=perm_mult(1,1)+(dc_l1**2+dc_l2**2)*10**perm_frac
      perm_mult(2,2)=perm_mult(2,2)+(dc_m1**2+dc_m2**2)*10**perm_frac
      perm_mult(3,3)=perm_mult(3,3)+(dc_n1**2+dc_n2**2)*10**perm_frac
      perm_mult(1,2)=perm_mult(1,2)+
     &(dc_l1*dc_m1+dc_l2*dc_m2)*10**perm_frac
      perm_mult(1,3)=perm_mult(1,3)+
     &(dc_l1*dc_n1+dc_l2*dc_n2)*10**perm_frac
      perm_mult(2,3)=perm_mult(2,3)+
     &(dc_m1*dc_n1+dc_m2*dc_n2)*10**perm_frac
      perm_denom(1,1)=perm_denom(1,1)+(dc_l1**2+dc_l2**2)
      perm_denom(2,2)=perm_denom(2,2)+(dc_m1**2+dc_m2**2)
      perm_denom(3,3)=perm_denom(3,3)+(dc_n1**2+dc_n2**2)
      perm_denom(1,2)=perm_denom(1,2)+(dc_l1*dc_m1+dc_l2*dc_m2)
      perm_denom(1,3)=perm_denom(1,3)+(dc_l1*dc_n1+dc_l2*dc_n2)
      perm_denom(2,3)=perm_denom(2,3)+(dc_m1*dc_n1+dc_m2*dc_n2)
	  
      return
      end
