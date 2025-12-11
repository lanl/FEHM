      subroutine random_walk(iflag,np1, current_node, current_model,
     1	edt,x_ts, y_ts, z_ts,itime,istep,time_step_flag)

!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1 	
!D1 PURPOSE
!D1 	
!D1 Perform random walk calculation in the streamline particle tracking 
!D1 algorithm.	
!D1	
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: ?, Programmer: ?
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/random_walk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Replaced random_seed / random_number with ran_sp
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:20   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************

      use comai , only : ierr
      use comdi
      use comsptr
      use compart
      use comsk, only : omr_flag
      
      implicit none
      
      integer np1, current_node, current_model
      integer i,j,k
      integer iflag,itime,istep
      integer time_step_flag
      
C     real rand_return
      real ran_sp
      real*8 x_ts, y_ts, z_ts
      
      real*8 time_random, edt
      real*8 rnum(3),b(3),r(3,3),dx(3),bdisp(3,3),dx_rand(3)
      real*8 coef1, coef2,at,rnum1,rnum2,rnum3
      real*8 termx1,termx2,termx3,termy1,termy2,termy3
      real*8 termz1,termz3
      real*8  ux,uy,uz,uxux,uyuy,uzuz,uxuz,beta,betau
      real*8 u,uu,al,alv,alh,ath,atv,dtsqrt,uxuy,uyuz,sbeta
      real*8 dterm(3),dxterm(3),dispterm(3),dphiterm(3)
      
      real*8 dm,dmm,a1,a2,a3,a4,asx,asy,asz,sig,sig2,ap,bp,cp,dp
      real*8 term,root2,zeta1,zeta2,root1i,fac1,fac2
      real*8 eig1,eig2,eig3,vca1,vca2,vca3,avca
      
      if(tprpflag(current_model).eq.2.or.
     $     tprpflag(current_model).eq.4.or.
     $     tprpflag(current_model).eq.13.or.
     $     tprpflag(current_model).eq.14) then
         
         
c     initialize random walk displacement to zero.
c     unless isotropic_flag = 0,1,2,4,5,6;  on return dx=0
c     and x2-in = x2-out
         
         dx = 0.
         dx_rand=0.
         
         dm=dispersivity1(current_model)
         
         
         ddxv(np1) = ddx(ijkv(np1))
         ddyv(np1) = ddy(ijkv(np1))
         ddzv(np1) = ddz(ijkv(np1))
         
         ux=(ggg(current_node,1)-ggg(current_node,-1))*.5
         uy=(ggg(current_node,2)-ggg(current_node,-2))*.5
         uz=(ggg(current_node,3)-ggg(current_node,-3))*.5
         
         uxux=ux*ux
         uyuy=uy*uy
         uzuz=uz*uz
         uxuz=ux*uz
         uxuy=ux*uy
         uyuz=uy*uz
         
         vx(np1)=ux
         vy(np1)=uy
         vz(np1)=uz
         
         uu=uxux+uyuy+uzuz
         u = sqrt(uu)
         if(u.gt.1.d-30) then
            
            
c.......................................
            
c     pre-dec-2000 way of doing the random tensor
            
            if(abs(itensor).eq.5) then
               
c     
               beta=dsqrt((uu+2.*uxuz+uyuy))
               betau=beta*u
               r=0.
               r(1,1)=ux/u
               r(2,1)=uy/u
               r(3,1)=uz/u
               r(1,2)=-uy/beta
               r(2,2)=(ux+uz)/beta
               r(3,2)=r(1,2)
               r(1,3)=-(uyuy+uzuz+uxuz)/betau
               r(2,3)=uy*(ux-uz)/betau
               r(3,3)=(uxux+uyuy+uxuz)/betau
               
               al =dispersivity1(current_model)
               ath=dispersivity2(current_model) 
               atv=dispersivity3(current_model)
               dm =dispersivity4(current_model)
               
               b=0.
               b(1)=sqrt(2.*(al*u+dm))
               b(2)=sqrt(2.*((ath*(uxux+uyuy)+atv*uzuz)/u+dm))
               b(3)=sqrt(2.*(atv*u+dm))
               
c     bdisp used below for either the time step determination
c     or the random walk determination
               
               do i=1,3
                  do j=1,3
                     bdisp(j,i)=r(j,i)*b(i)
                  enddo
               enddo
               
               
c     do i=1,3
c     !		    call random_number(rand_return)
c     !		    rnum(i) = (rand_return - .5)*dtsqrt
C     rnum(i) = (ran_sp(rseed) - .5)*dtsqrt
c     enddo
               
c     do i=1,3
c     dx(i)=0.
c     do k=1,3
c     dx(i)=dx(i)+r(i,k)*b(k)*rnum(k)
c     enddo
c     enddo
               
               
            elseif(abs(itensor).eq.1) then
c........................
               
c     generalized form of the dispersion tensor, axisymmetric 
c     formulation, 
c     assumed that axis of symmetry oblique to verticle and to the flow
c     Lichtner et al 99, eq 41-47 and my notes dated 10/25/99
c     here a1,a2,a3,a4 are the 4 generatized dispersivities from eq 41
c     and asx, asy,asz are the direction cosins of the axis of symmetry 
c     denoted by
c     Lambda in eq41
c     vxlambda is an eigen vector of D, the other two are constructed 
c     from linear
c     combination of v and lambda
c     the columns of U (denoted here by r(i,j)) are made up of the
c     eigenvectors. 
c     vas = unit vector in velocity direction. dot. axis os symmetry
               
               a1=dispersivity2(current_model)
               a2=dispersivity3(current_model)
               a3=dispersivity4(current_model)
               a4=dispersivity6(current_model)
               asx=dispersivity7(current_model)
               asy=dispersivity8(current_model)
               asz=sqrt(1.-asx*asx-asy*asy)
               sig=(asx*ux+asy*uy+asz*uz)/u
               sig2=sig*sig
               if((1.-sig2).lt.1.e-20) then
c     v parallel to axis of symmetry, requires a different formulation
                  write(ierr,*)'v parallel to axis of symmetry, '
                  write(ierr,*)' requires a different formulation. stop'
                  stop
               else
                  ap=a1+a2+sig2*a3+sig*a4
                  bp=sig*a3+0.5*a4
                  cp=(1.-sig2)*bp
                  dp=a1+(1.-sig2)*a3
                  
                  term=(ap-dp)**2.+4.*bp*cp
                  if(term.lt.0.) then
                     write(ierr,*)'error, term.le.0 in random_walk.'
                     write(ierr,*)'check dispersiv coeff in sptr macro.'
                     write(ierr,*)'STOP'
                     stop
                  endif
                  term=sqrt(term)
                  
                  eig1=0.5*((ap+dp)+term)
                  eig2=0.5*((ap+dp)-term)
                  
                  root2=-cp/(ap-eig2)
c     11/18/99 s kelkar NOTE: refere to eq 37-40 in lochtner et al
c     writing xsi-2=(root2*v-hat+omega-hat)/norm but
c     xsi-1=(v-hat+inv.root1*omega-hat)/norm
c     in order to facilitate the limit of the axissymetrc case
c     with alpha-3 = alpha-4 =0.
c     also note that the equations bellow are cast in terms of lamda-hat 
c     instead of omega-hat
                  
                  root1i=-bp/(dp-eig1)
                  
                  zeta2=sqrt(root2*root2+1.-sig2)
                  zeta1=sqrt(1.+(1.-sig2)*root1i*root1i)
                  
                  eig3=a1
                  
                  b(1)=sqrt(2.*eig1*u)
                  b(2)=sqrt(2.*eig2*u)
                  b(3)=sqrt(2.*eig3*u)
                  
                  fac1=(1.-root1i*sig)/u
                  fac2=(root2-sig)/u
                  
                  r(1,1)=(fac1*ux+root1i*asx)/zeta1
                  r(2,1)=(fac1*uy+root1i*asy)/zeta1
                  r(3,1)=(fac1*uz+root1i*asz)/zeta1
                  r(1,2)=(fac2*ux+asx)/zeta2
                  r(2,2)=(fac2*uy+asy)/zeta2
                  r(3,2)=(fac2*uz+asz)/zeta2
                  vca1=+uy*asz-uz*asy
                  vca2=-ux*asz+uz*asx
                  vca3=+ux*asy-uy*asx
                  avca=sqrt(vca1*vca1+vca2*vca2+vca3*vca3)
                  r(1,3)=vca1/avca
                  r(2,3)=vca2/avca
                  r(3,3)=vca3/avca
                  
                  do i=1,3
                     do j=1,3
                        bdisp(j,i)=r(j,i)*b(i)
                     enddo
                  enddo
                  
                  
               endif
               
c     do i=1,3
c     !		    call random_number(rand_return)
c     !		    rnum(i) = (rand_return-.5)*dtsqrt
C     rnum(i) = (ran_sp(rseed) - .5)*dtsqrt
c     enddo
               
c     do i=1,3
c     dx(i)=0.
c     do k=1,3
c     dx(i)=dx(i)+bdisp(i,k)*rnum(k)
c     enddo
c     enddo
               
               
            elseif(abs(itensor).eq.2) then
               
c..............................
               
c     Burnett and Frind tensor, axisymmetric 
c     formulation, 
c     assumed that axis of symmetry is along z axis.
c     Lichtner et al 99, eq 59-62 and 76-78
               r=0.
               
               beta=(uxux+uyuy)
               
               if(beta.lt.1.e-34) then
                  r(1,3)=-1.
                  r(2,2)=+1.
                  r(3,1)=+1.
                  sig2 = 1.
               else
                  sbeta=sqrt(beta)
                  sig=uz/u
                  sig2=sig*sig
                  
                  r=0.
                  
                  r(1,1)=+ux/u
                  r(2,1)=+uy/u
                  r(3,1)=+uz/u
                  
                  r(1,3)=-uxuz/(u*sbeta)
                  r(2,3)=-uyuz/(u*sbeta)
                  r(3,3)=+sbeta/u
                  
                  r(1,2)=-uy/sbeta
                  r(2,2)=+ux/sbeta
                  r(3,2)=0.
                  
               endif
               
               al  =dispersivity2(current_model)
               ath =dispersivity3(current_model) 
               atv =dispersivity4(current_model)
               
               a1=ath+sig2*(atv-ath)
               
               b=0.
               b(1)=sqrt(2.*(al*u+dm))
               b(3)=sqrt(2.*(atv*u+dm))
               b(2)=sqrt(2.*(a1*u+dm))
               
               do i=1,3
                  do j=1,3
                     bdisp(j,i)=r(j,i)*b(i)
                  enddo
               enddo
               
               
c     do i=1,3
c     !		    call random_number(rand_return)
c     !		    rnum(i) = (rand_return-.5)*dtsqrt
C     rnum(i) = (ran_sp(rseed) - .5)*dtsqrt
c     enddo
               
c     do i=1,3
c     dx(i)=0.
c     do k=1,3
c     dx(i)=dx(i)+r(i,k)*b(k)*rnum(k)
c     enddo
c     enddo
               
            elseif(abs(itensor).eq.3) then
               
c................................................
               
c     modified Burnett and friend, Lichtner et al eq 57-61,
c     formulation, Lamda assumed to be along z axis 
c     assumed that axis of symmetry oblique to verticle and to the flow
c     here a1,a2,a3,a4 are the 4 generatized dispersivities from eq 41
c     vxlambda is an eigen vector of D, the other two are constructed 
c     from linear
c     combination of v and lambda
c     the columns of U (denoted here by r(i,j)) are made up of the
c     eigenvectors. 
c     vas = unit vector in velocity direction. dot. axis os symmetry
               
               alh=dispersivity2(current_model)
               alv=dispersivity3(current_model)
               ath=dispersivity4(current_model)
               atv=dispersivity6(current_model)
c     asx=dispersivity7(current_model)
c     asy=dispersivity8(current_model)
c     asz=sqrt(1.-asx*asx-asy*asy)
               
               asx=0.
               asy=0.
               asz=1.
               
               beta=(uxux+uyuy)
               sig=(asx*ux+asy*uy+asz*uz)/u
               sig2=sig*sig
               
               if((1.-sig2).lt.1.e-20) then
c     v parallel to axis of symmetry, requires a different formulation
                  r=0.
                  r(1,3)=-1.
                  r(2,2)=+1.
                  r(3,1)=+1.
                  
                  b=0.
                  b(1)=sqrt(2.*(alv*u+dm))
                  b(2)=sqrt(2.*(ath*u+dm))
                  b(3)=sqrt(2.*(ath*u+dm))
                  
               else
                  
                  al=alh+sig2*(alv-alh)
                  at=atv+sig2*(ath-atv)
                  
                  sbeta=sqrt(beta)
                  
                  r=0.
                  
                  r(1,1)=+ux/u
                  r(2,1)=+uy/u
                  r(3,1)=+uz/u
                  
                  r(1,2)=-uxuz/(u*sbeta)
                  r(2,2)=-uyuz/(u*sbeta)
                  r(3,2)=+sbeta/u
                  
                  r(1,3)=-uy/sbeta
                  r(2,3)=+ux/sbeta
                  r(3,3)=0.
                  
                  b=0.
                  b(1)=sqrt(2.*(al*u+dm))
                  b(2)=sqrt(2.*(at*u+dm))
                  b(3)=sqrt(2.*(ath*u+dm))
                  
                  do i=1,3
                     do j=1,3
                        bdisp(j,i)=r(j,i)*b(i)
                     enddo
                  enddo
                  
               endif
               
c     do i=1,3
c     !		    call random_number(rand_return)
c     !		    rnum(i) = (rand_return-.5)*dtsqrt
C     rnum(i) = (ran_sp(rseed) - .5)*dtsqrt
c     enddo
               
c     do i=1,3
c     dx(i)=0.
c     do k=1,3
c     dx(i)=dx(i)+bdisp(i,k)*rnum(k)
c     enddo
c     enddo
               
               
            elseif(abs(itensor).eq.4) then
               
c..............................
               
c     isotropic formulation of the tensor as in thompson's report
               
               al=dispersivity2(current_model)
               at=dispersivity3(current_model) 
               
c     coef1=dtsqrt*dsqrt(2.*(al*u+dm))
               beta =sqrt(uu+2.*uxuz+uyuy)
               coef2=sqrt(at/al)/beta
               
c     !		 call random_number(rand_return)
c     !		 rnum1 = (rand_return-.5)
C     rnum1 = (ran_sp(rseed) - .5)
c     !		 call random_number(rand_return)
c     !		 rnum2 = (rand_return-.5)
C     rnum2 = (ran_sp(rseed) - .5)
c     !		 call random_number(rand_return)
c     !		 rnum3 = (rand_return-.5)
C     rnum3 = (ran_sp(rseed) - .5)
c     x-displacement
c     termx1 = +rnum1*ux/u
c     termx2 = -rnum2*coef2*uy
c     termx3 = -rnum3*coef2*(uyuy+uzuz+uxuz)/u
c     dx(1) = coef1*(termx1+termx2+termx3)
               
               b(1) = dsqrt(2.*(al*u+dm))
               b(2) = dsqrt(2.*(al*u+dm))
               b(3) = dsqrt(2.*(al*u+dm))
               r(1,1) = ux/u
               r(1,2) = -coef2*uy
               r(1,3) = -coef2*(uyuy+uzuz+uxuz)/u
               r(2,1) = uy/u
               r(2,2) = coef2*(ux+uz)
               r(2,3) = coef2*uy*(ux-uz)/u
               r(3,1) = uz/u
               r(3,2) = -coef2*uy
               r(3,3) = coef2*(uxux+uyuy+uxuz)/u
               
               do i=1,3
                  do j=1,3
                     bdisp(j,i)=r(j,i)*b(i)
                  enddo
               enddo
               
c     y-displacement
c     termy1 = +rnum1*uy/u
c     termy2 = +rnum2*coef2*(ux+uz)
c     termy3 = +rnum3*coef2*uy*(ux-uz)/u
c     dx(2) = coef1*(termy1+termy2+termy3)
c     z-displacement
c     term2 for dx and dz are the same
c     termz1 = +rnum1*uz/u
c     c       termz2 = -rnum2*coef2*uy    !  Line commented out before
c     c                                      switching to bdisp formulation
c     termz3 = +rnum3*coef2*(uxux+uyuy+uxuz)/u
c     dx(3) = coef1*(termz1+termx2+termz3)
            end if
c     End of itensor options, this else is for
c     no velocity (i.e. diffusion only)
         else
c     if velocity is zero, only molecular dispersion
            
c     !	      call random_number(rand_return)
c     !	      rnum1 = (rand_return-.5)
C     rnum1 = (ran_sp(rseed) - .5)
c     !	      call random_number(rand_return)
c     !	      rnum2 = (rand_return-.5)
C     rnum2 = (ran_sp(rseed) - .5)
c     !	      call random_number(rand_return)
c     !	      rnum3 = (rand_return-.5)
C     rnum3 = (ran_sp(rseed) - .5)
c     dmm=dtsqrt*sqrt(2.*dm)
            
c     bdisp is diagonal to reproduce molecular diffusion
c     random walk below
            bdisp = 0.
            bdisp(1,1) = sqrt(2.*dm)
            bdisp(2,2) = sqrt(2.*dm)
            bdisp(3,3) = sqrt(2.*dm)
c     dx(1)=dmm*rnum1
c     dx(2)=dmm*rnum2
c     dx(3)=dmm*rnum3
            
         endif
         
c     iflag is used to set what is being computed
c     0 - determine characteristic time step for establishing
c     minimum time step for calculation
c     not 0 - determine dx(1:3) = random walk terms
         if(iflag .eq. 0) then
            
c     Time step determination (no random walk term)
            do i=1,3
               dx(i)=0.
               do k=1,3
                  dx(i)=dx(i)+bdisp(i,k)
               enddo
            enddo
            
            x_ts=(ddxv(np1)/max(1.d-30,abs(dx(1))))**2
            y_ts=(ddyv(np1)/max(1.d-30,abs(dx(2))))**2
            z_ts=(ddzv(np1)/max(1.d-30,abs(dx(3))))**2
            
         else
            
c     Determine random walk terms and gradd if called for
c     Make time step not any larger than the selected time step
c     to prevent it from getting large when edt is large enough
c     for particles to traverse several cells
            
            time_random =dmin1(edtmax(np1),edt,dt(np1))
            dtsqrt=dsqrt(time_random)*sqrt(3.)*2.
            
c     the sqrt(3.)*2. is to make the random numbers have <z^2>=1         
            
c     Random numbers (3 independent values)
c     dtsqrt included for each
            
            do i=1,3
               
                                ! ZVD implement with ran_sp
c     !		 call random_number(rand_return)
c     !		 rnum(i) = (rand_return-.5)*dtsqrt
               
               rnum(i) = (ran_sp(rseed) - .5)*dtsqrt
               
            enddo
            
c     random walk displacement in each direction computed
            
            do i=1,3
               dx_rand(i)=0.
               do k=1,3
                  dx_rand(i)=dx_rand(i)+bdisp(i,k)*rnum(k)
               enddo
            enddo
            
            
c..............................
            
c     add the grad-dispersion term here. gradd was calculated in 
c     subroutine gradd_BF which was called from ptrac3 before the 
c     loop on all particles.
            
            
c     The sign of itensor is the flag that says if we are to
c     use the gradd term or not
            
            if(itensor.gt.0) then
               
c...  8/29/01, 3/19/02 11/1/02  s kelkar..............
c     add the div-dispersion and grad(ln(porosity)) terms here.
c we are using the bilinear interpolation of Darcy velocities in
c     calculating the Divergence of Dispersion tensor term
c     the gradient of diffusivity is being ignored for the present.
               
               do i=1,3
               enddo
               
               do i=1,3
                  if(omr_flag) then
                     dispterm(i)=divd_omr(current_node,i)
     $                    *time_random
                     dphiterm(i)=dm*dpordx_omr(current_node,i)
     $                    *time_random
                  else
                     dispterm(i)=divd(i)*time_random
                     dphiterm(i)=dm*dpordx(i)*time_random
                  endif
               enddo
               
c...............................................................
               
               do i=1,3
                  dx(i)=dx(i)+dx_rand(i)
               enddo
               
               do i=1,3
                  dx(i)=dx(i)+dispterm(i)
               enddo
               
               do i=1,3
                  dx(i)=dx(i)+dphiterm(i)
               enddo
               
               
            end if
            
c     New x2, y2, z2 location
            
            x2(np1)=x1(np1)+dx(1)
            y2(np1)=y1(np1)+dx(2)
            z2(np1)=z1(np1)+dx(3)
            
c     Assign default values for these (not used for iflag=1)
            x_ts = 1.d30
            y_ts = 1.d30
            z_ts = 1.d30
            
         end if
         
         
c.........................
         
      else
         
c     No dispersion, set time step limits large (iflag=0)
c     or x2=x1 etc., for iflag .ne. 0
         if(iflag .eq. 0) then
            x_ts = 1.d30
            y_ts = 1.d30
            z_ts = 1.d30
         else
            x2(np1)=x1(np1)
            y2(np1)=y1(np1)
            z2(np1)=z1(np1)
         end if
         
      end if
      
      return
      
      end
      
c........................................................................
