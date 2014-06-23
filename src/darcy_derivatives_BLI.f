      subroutine darcy_derivatives_BLI(i,i1,j1,k1,xp,yp,zp,
     $                    dvd,dvxd,dvyd,dvzd,v,vv,pormax)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Calculate velocity derivatives needed for evaluating the
!D1 dispersion gradient, using a trilinear interpolation of the 
!D1 velocities within the subvolume that contains the particle
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 
!D2 Initial implementation:Sep 02, Programmer: S. Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/darcy_derivatives_BLI.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:01:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!D4 grad.D term from Tompson and Gelhar, WRR, Oct 90, eq. 8
!D4 the flux derivatives are calculated using a trilinear 
!D4 interpolation which is an extention of the method outlined by 
!D4 LeBolle et al, WRR, March 1996, Table 1, pp 588 and then 
!D4 converted to velocity derivatives. This avoids having to calculate
!D4 the derivatives of porosity
!**********************************************************************

      use combi
      use comci
      use comdi
      use comsptr

      implicit none
      integer i
      real*8 dvd(3),dvxd(3),dvyd(3),dvzd(3),v(3),vv

      integer i1,j1,k1,ix,iy,iz
      integer iisk,iyz,izx,iyx,ixy

      real*8 xp,yp,zp
      real*8 ux1,ux2,ux3,ux4,ux5,ux6,ux7,ux8
      real*8 uy1,uy2,uy3,uy4,uy5,uy6,uy7,uy8
      real*8 uz1,uz2,uz3,uz4,uz5,uz6,uz7,uz8
      real*8 dvdx, dvdy, dvdz
      real*8 dx,dvxdx,dvxdy,dvxdz
      real*8 dvxdyr,dvxdyl,dvxdzr,dvxdzl
      real*8 dvydxr,dvydxl,dvydzr,dvydzl
      real*8 dvzdyr,dvzdyl,dvzdxr,dvzdxl
      real*8 dy,dvydx,dvydy,dvydz
      real*8 dz,dvzdx,dvzdy,dvzdz
      real*8 fx,fy,fz
      real*8 length_factor,pormax

      real*8 gi(-3:3),gix(-3:3),giy(-3:3),giz(-3:3),giyz(-3:3)
      real*8 gizx(-3:3),giyx(-3:3),gixy(-3:3)

      parameter(length_factor=2.)

c      open(unit=98,file='dispersion_debug.out')
c      write(98,*)' using darcy_derivatives_BLI.f'
c      close(98)

      ix =irray(i ,  i1)
c account for boundary nodes, s kelkar, 12/19 02
      if(ix.le.0) ix =i 
      iy =irray(i ,2*j1)
      if(iy.le.0) iy = i
      iz =irray(i ,3*k1)
      if(iz.le.0) iz = i
      iyz=irray(iy,3*k1)
      if(iyz.le.0) iyz = iy
      izx=irray(iz,  i1)
      if(izx.le.0) izx = iz
      iyx=irray(iy,  i1)
      if(iyx.le.0) iyx = iy
      ixy=irray(ix,2*j1)
      if(ixy.le.0) ixy = ix


c .. ....skelkar  may 22 02
c find the max porosity amongst the connected nodes. This will be used in 
c dispersion_divergence for calculating (DIv(Phi.D))/pormax. this is being 
c done to avoide large jumps in (DIv(Phi.D))/porosity resulting from
c those in ps(I)
c zvd 13-May-08 change ps(*) to ps_trac(*)
      pormax=ps_trac(i)
      if(ps_trac(ix ).gt.pormax) pormax=ps_trac(ix )
      if(ps_trac(iy ).gt.pormax) pormax=ps_trac(iy )
      if(ps_trac(iz ).gt.pormax) pormax=ps_trac(iz )
      if(ps_trac(iyz).gt.pormax) pormax=ps_trac(iyz)
      if(ps_trac(izx).gt.pormax) pormax=ps_trac(izx)
      if(ps_trac(iyx).gt.pormax) pormax=ps_trac(iyx)
      if(ps_trac(ixy).gt.pormax) pormax=ps_trac(ixy)
c...............................


c calculate the darcy velocities at nodes of interest
      do iisk=-3,3,+1
         gi(iisk)=ggg(i,iisk)*ps_trac(i)
         gix(iisk)=ggg(ix,iisk)*ps_trac(ix)
         giy(iisk)=ggg(iy,iisk)*ps_trac(iy)
         giz(iisk)=ggg(iz,iisk)*ps_trac(iz)
      enddo
      giyz(+1)=ggg(iyz,+1)*ps_trac(iyz)
      giyz(-1)=ggg(iyz,-1)*ps_trac(iyz)
      gizx(+2)=ggg(izx,+2)*ps_trac(izx)
      gizx(-2)=ggg(izx,-2)*ps_trac(izx)
      giyx(+3)=ggg(iyx,+3)*ps_trac(iyx)
      giyx(-3)=ggg(iyx,-3)*ps_trac(iyx)
      gixy(+3)=ggg(ixy,+3)*ps_trac(ixy)
      gixy(-3)=ggg(ixy,-3)*ps_trac(ixy)
c
      dx=(cord(ix,1)-cord(i,1))/length_factor
      dy=(cord(iy,2)-cord(i,2))/length_factor
      dz=(cord(iz,3)-cord(i,3))/length_factor
c account for boundary nodes, s kelkar, 12/19 02
      if(i.eq.ix) dx=(cord(irray(i,-i1),1)-cord(i,1))/length_factor
      if(i.eq.iy) dy=(cord(irray(i,-2*j1),2)-cord(i,2))/length_factor
      if(i.eq.iz) dz=(cord(irray(i,-3*k1),3)-cord(i,3))/length_factor
      
      fx=(xp-cord(i,1))/dx
      fy=(yp-cord(i,2))/dy
      fz=(zp-cord(i,3))/dz

c     coeffecients of interpolation for the x component of Darcy velocity
      ux1=0.5*(gi(+1)-gi(-1))
      ux2=i1*gi(i1)
      ux3=0.5*(ux2+i1*giy(i1))
      ux4=0.5*(ux1+0.5*(giy(+1)
     1     -giy(-1)))
      ux5=0.5*(ux1+0.5*(giz(+1)
     1     -giz(-1)))
      ux6=0.5*(ux2+i1*giz(i1))
      ux7=0.5*(ux3+0.5*(i1*giyz(i1)
     $     +i1*giz(i1)))
      ux8=0.5*(ux4+0.25*(giyz(+1)
     1     -giyz(-1)
     2     +giz(+1)
     3     -giz(-1)))
      
c     coeffecients of interpolation for the y component of Darcy velocity
      uy1=0.5*(gi(+2)-gi(-2))
      uy2=0.5*(uy1+0.5*(gix(+2)
     1     -gix(-2)))
      uy3=0.5*(j1*gix(2*j1)
     1     +j1*gi(2*j1))
      uy4=j1*gi(2*j1)
      uy5=0.5*(uy1+0.5*(giz(+2)
     1     -giz(-2)))
      uy6=0.5*(uy2+0.25*(giz(+2)
     1     -giz(-2)
     2     +gizx(+2)
     3     -gizx(-2)))
      uy7=0.5*(uy3+0.5*(j1*giz(2*j1)
     1     +j1*gizx(2*j1)))
      uy8=0.5*(uy4+j1*giz(2*j1))
      
c     coeffecients of interpolation for the z component of Darcy flux
      uz5=k1*gi(3*k1)
      uz6=0.5*(uz5+k1*gix(3*k1))
      uz7=0.5*(uz6+k1*giy(3*k1)
     1     +k1*giyx(3*k1))
      uz8=0.5*(uz5+k1*giy(3*k1))
      uz1=0.5*(uz5-k1*gi(-3*k1))
      uz2=0.5*(uz6+0.5*(-k1*gi(-3*k1)
     1     -k1*gix(-3*k1)))
      uz3=0.5*(uz7+0.25*(-k1*gi(-3*k1)
     1     -k1*gix(-3*k1)
     2     -k1*giy(-3*k1)
     3     -k1*gixy(-3*k1)))
      uz4=0.5*(uz8+0.5*(-k1*gi(-3*k1)
     1     -k1*giy(-3*k1)))

c interpolate Darcy velocity at the particle location
      v(1)=(1.-fx)*(1.-fy)*(1.-fz)*ux1 +  fx    *(1.-fy)*(1.-fz)*ux2
     1    + fx    *fy     *(1.-fz)*ux3 + (1.-fx)*fy     *(1.-fz)*ux4
     2    +(1.-fx)*(1.-fy)*fz     *ux5 + fx     *(1.-fy)*fz     *ux6
     3    +fx     *fy     *fz     *ux7 + (1.-fx)*fy     *fz     *ux8

      v(2)=(1.-fx)*(1.-fy)*(1.-fz)*uy1 +  fx    *(1.-fy)*(1.-fz)*uy2
     1    + fx    *fy     *(1.-fz)*uy3 + (1.-fx)*fy     *(1.-fz)*uy4
     2    +(1.-fx)*(1.-fy)*fz     *uy5 + fx     *(1.-fy)*fz     *uy6
     3    +fx     *fy     *fz     *uy7 + (1.-fx)*fy     *fz     *uy8

      v(3)=(1.-fx)*(1.-fy)*(1.-fz)*uz1 +  fx    *(1.-fy)*(1.-fz)*uz2
     1    + fx    *fy     *(1.-fz)*uz3 + (1.-fx)*fy     *(1.-fz)*uz4
     2    +(1.-fx)*(1.-fy)*fz     *uz5 + fx     *(1.-fy)*fz     *uz6
     3    +fx     *fy     *fz     *uz7 + (1.-fx)*fy     *fz     *uz8

c derivatives of the x-darcy velocity
      dvxd(1)=(     - (1.-fy)*(1.-fz)*ux1      +(1.-fy)*(1.-fz)*ux2
     1              +      fy*(1.-fz)*ux3      -     fy*(1.-fz)*ux4
     2              -      (1.-fy)*fz*ux5      +(1.-fy)*fz     *ux6
     3              +      fy*fz     *ux7      -     fy*fz     *ux8)/dx
      dvxd(2)=(-(1.-fx)      *(1.-fz)*ux1-fx           *(1.-fz)*ux2
     1         +fx           *(1.-fz)*ux3+(1.-fx)      *(1.-fz)*ux4
     2         -(1.-fx)           *fz*ux5-fx                *fz*ux6
     3         +fx                *fz*ux7+(1.-fx)           *fz*ux8)/dy
      dvxd(3)=(-(1.-fx)*(1.-fy)      *ux1-fx     *(1.-fy)      *ux2
     1         -fx     *fy           *ux3-(1.-fx)*fy           *ux4
     2          +(1.-fx)*(1.-fy)     *ux5+fx     *(1.-fy)      *ux6
     3          +fx     *fy          *ux7+(1.-fx)*fy           *ux8)/dz

c derivatives of the y-darcy velocity
      dvyd(1)=(      -(1.-fy)*(1.-fz)*uy1      +(1.-fy)*(1.-fz)*uy2
     1               +fy     *(1.-fz)*uy3      -fy     *(1.-fz)*uy4
     2               -(1.-fy)*fz     *uy5      +(1.-fy)     *fz*uy6
     3               +fy     *fz     *uy7      -fy          *fz*uy8)/dx
      dvyd(2)=(-(1.-fx)      *(1.-fz)*uy1-fx           *(1.-fz)*uy2
     1       +fx             *(1.-fz)*uy3+(1.-fx)      *(1.-fz)*uy4
     2       -(1.-fx)        *fz     *uy5-fx           *fz     *uy6
     3       +fx             *fz     *uy7+(1.-fx)      *fz    *uy8)/dy
      dvyd(3)=(-(1.-fx)*(1.-fy)      *uy1-fx    *(1.-fy)       *uy2
     1         -fx     *fy           *uy3-(1.-fx)*fy           *uy4
     2         +(1.-fx)*(1.-fy)      *uy5+fx     *(1.-fy)      *uy6
     3         +fx     *fy           *uy7+(1.-fx)*fy           *uy8)/dz

c derivatives of the z-darcy velocity
      dvzd(1)=(      -(1.-fy)*(1.-fz)*uz1+      (1.-fy)*(1.-fz)*uz2
     1               +fy     *(1.-fz)*uz3-      fy     *(1.-fz)*uz4
     2               -(1.-fy)*fz     *uz5+      (1.-fy)*fz     *uz6
     3               +fy     *fz     *uz7-      fy     *fz     *uz8)/dx
      dvzd(2)=(-(1.-fx)      *(1.-fz)*uz1-fx           *(1.-fz)*uz2
     1         +fx           *(1.-fz)*uz3+(1.-fx)      *(1.-fz)*uz4
     2         -(1.-fx)           *fz*uz5-fx                 *fz*uz6
     3         +fx           *fz     *uz7+(1.-fx)           *fz*uz8)/dy
      dvzd(3)=(-(1.-fx)*(1.-fy)      *uz1-fx     *(1.-fy)      *uz2
     1         -fx     *fy           *uz3-(1.-fx)*fy           *uz4
     2         +(1.-fx)*(1.-fy)      *uz5+fx     *(1.-fy)      *uz6
     3         +fx     *fy           *uz7+(1.-fx)     *fy      *uz8)/dz

c........................................................

      vv = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      v(1)=v(1)/max(1.d-30,vv)
      v(2)=v(2)/max(1.d-30,vv)
      v(3)=v(3)/max(1.d-30,vv)


      
c     derivatives of the abs. value of darcy velocity , v
      
      dvd(1)=v(1)*dvxd(1)+v(2)*dvyd(1)+v(3)*dvzd(1)
      dvd(2)=v(1)*dvxd(2)+v(2)*dvyd(2)+v(3)*dvzd(2)
      dvd(3)=v(1)*dvxd(3)+v(2)*dvyd(3)+v(3)*dvzd(3)
c     
      
      return

      end
      
      
