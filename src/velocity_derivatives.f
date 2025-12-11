      subroutine velocity_derivatives(i,dvd,dvxd,dvyd,dvzd,v,vv)
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
!D1 dispersion gradient.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 25-AUG-99, Programmer: S. Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/velocity_derivatives.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:32:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
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
!D4 grad.D term from Tompson and Gelhar, WRR, Oct 90, eq. 8
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
**********************************************************************

      use comsptr
      implicit none

      integer i
      real*8 dvd(3),dvxd(3),dvyd(3),dvzd(3),v(3),vv

      integer ix1, ix2, iy1, iy2, iz1, iz2
      real*8   vxir,vxil,vyir,vyil,vzir,vzil
      real*8   xy1r,xy1l,xy2r,xy2l,xz1r,xz1l,xz2r,xz2l
      real*8   yx1r,yx1l,yx2r,yx2l,yz1r,yz1l,yz2r,yz2l
      real*8   zx1r,zx1l,zx2r,zx2l,zy1r,zy1l,zy2r,zy2l
      real*8 dvdx, dvdy, dvdz
      real*8 dx,dx1,dx2,dvxdx,dvxdy,dvxdz
      real*8 dvxdyr,dvxdyl,dvxdzr,dvxdzl
      real*8 dvydxr,dvydxl,dvydzr,dvydzl
      real*8 dvzdyr,dvzdyl,dvzdxr,dvzdxl
      real*8 dy,dy1,dy2,dvydx,dvydy,dvydz
      real*8 dz,dz1,dz2,dvzdx,dvzdy,dvzdz
      real*8 length_factor
      parameter(length_factor=1.)

c velocities are defined at cell-faces along the direction normal
c to each face.
c dvx/dx is formed within the cell by refering to the two y-z
c faces normal to the x axis. However, to get the cross derivatives 
c have to refere to the neighbouring cells.
c eg: consider dvx/dy. Take the x-y plane going thru the node. In this
c plane, there are 6 x-velocities defined in the neighbourhood of this
c node, ie one on the face to the right, one on the face to the left, 
c then taking the nodes above and bellow, to the right and left of each
c of those. Then dvxdy is formed as the average centered difference.

c i is the current node, iy1 is the next one along +y axis,
c iy2 is the next one along -y axis etc.
         ix1 = irray(i,+1)
         ix2 = irray(i,-1)
         iy1 = irray(i,+2)
         iy2 = irray(i,-2)
         iz1 = irray(i,+3)
         iz2 = irray(i,-3)
c take care of the boundary and corner nodes. If i is on a face
c and there no node next to it, then the next node is defined as
c i itself, so that the terms in the derivative expressions involving
c that node drop out.
         if(ix1.eq.0) ix1=i
         if(ix2.eq.0) ix2=i
         if(iy1.eq.0) iy1=i
         if(iy2.eq.0) iy2=i
         if(iz1.eq.0) iz1=i
         if(iz2.eq.0) iz2=i
c form distances for the node centered cells
         dx=ddx(i)/length_factor
         dx1=0.5*(ddx(i)+ddx(ix1))/length_factor
         dx2=0.5*(ddx(i)+ddx(ix2))/length_factor
         dy=ddy(i)/length_factor
         dy1=0.5*(ddy(i)+ddy(iy1))/length_factor
         dy2=0.5*(ddy(i)+ddy(iy2))/length_factor
         dz=ddz(i)/length_factor
         dz1=0.5*(ddz(i)+ddz(iz1))/length_factor
         dz2=0.5*(ddz(i)+ddz(iz2))/length_factor
c recall velocities at the current cell faces
c note the -ve sign in front of the gg terms refering to - faces
         vxir=+ggg(i,+1)
         vxil=-ggg(i,-1)
         vyir=+ggg(i,+2)
         vyil=-ggg(i,-2)
         vzir=+ggg(i,+3)
         vzil=-ggg(i,-3)
c dvx/dx.
         dvxdx=(vxir-vxil)/dx
c dvx/dy
c velocities are named as follows: xy1r is the x velocity for the
c cell corrosponding to node y1 on the +y face etc.
         xy1r=+ggg(iy1,+1)
         xy1l=-ggg(iy1,-1)
         xy2r=+ggg(iy2,+1)
         xy2l=-ggg(iy2,-1)
         dvxdyl=0.5*(xy1l-vxil)/dy1
     $         +0.5*(vxil-xy2l)/dy2
         dvxdyr=0.5*(xy1r-vxir)/dy1
     $         +0.5*(vxir-xy2r)/dy2
         dvxdy=0.5*(dvxdyl+dvxdyr)
c dvx/dz
         xz1r=+ggg(iz1,+1)
         xz1l=-ggg(iz1,-1)
         xz2r=+ggg(iz2,+1)
         xz2l=-ggg(iz2,-1)
         dvxdzl=0.5*(xz1l-vxil)/dz1
     $         +0.5*(vxil-xz2l)/dz2
         dvxdzr=0.5*(xz1r-vxir)/dz1
     $         +0.5*(vxir-xz2r)/dz2
         dvxdz=0.5*(dvxdzl+dvxdzr)

         dvxd(1)=dvxdx
         dvxd(2)=dvxdy
         dvxd(3)=dvxdz

c dvy/dy.
         dvydy=(vyir-vyil)/dy
c dvy/dx
         yx1r=+ggg(ix1,+2)
         yx1l=-ggg(ix1,-2)
         yx2r=+ggg(ix2,+2)
         yx2l=-ggg(ix2,-2)
         dvydxl=0.5*(yx1l-vyil)/dx1
     $         +0.5*(vyil-yx2l)/dx2
         dvydxr=0.5*(yx1r-vyir)/dx1
     $         +0.5*(vyir-yx2r)/dx2
         dvydx=0.5*(dvydxl+dvydxr)
c dvy/dz
         yz1r=+ggg(iz1,+2)
         yz1l=-ggg(iz1,-2)
         yz2r=+ggg(iz2,+2)
         yz2l=-ggg(iz2,-2)
         dvydzl=0.5*(yz1l-vyil)/dz1
     $         +0.5*(vyil-yz2l)/dz2
         dvydzr=0.5*(yz1r-vyir)/dz1
     $         +0.5*(vyir-yz2r)/dz2
         dvydz=0.5*(dvydzl+dvydzr)

         dvyd(1)=dvydx
         dvyd(2)=dvydy
         dvyd(3)=dvydz

c dvz/dz.
         dvzdz=(vzir-vzil)/dz
c dvz/dx
         zx1r=+ggg(ix1,+3)
         zx1l=-ggg(ix1,-3)
         zx2r=+ggg(ix2,+3)
         zx2l=-ggg(ix2,-3)
         dvzdxl=0.5*(zx1l-vzil)/dx1
     $         +0.5*(vzil-zx2l)/dx2
         dvzdxr=0.5*(zx1r-vzir)/dx1
     $         +0.5*(vzir-zx2r)/dx2
         dvzdx=0.5*(dvzdxl+dvzdxr)
c dvz/dy
         zy1r=+ggg(iy1,+3)
         zy1l=-ggg(iy1,-3)
         zy2r=+ggg(iy2,+3)
         zy2l=-ggg(iy2,-3)
         dvzdyl=0.5*(zy1l-vzil)/dy1
     $         +0.5*(vzil-zy2l)/dy2
         dvzdyr=0.5*(zy1r-vzir)/dy1
     $         +0.5*(vzir-zy2r)/dy2
         dvzdy=0.5*(dvzdyl+dvzdyr)

         dvzd(1)=dvzdx
         dvzd(2)=dvzdy
         dvzd(3)=dvzdz

c derivatives of the abs. value of velocity , v

         v(1)=0.5*(vxir+vxil)
         v(2)=0.5*(vyir+vyil)
         v(3)=0.5*(vzir+vzil)
         vv = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         v(1)=v(1)/max(1.d-30,vv)
         v(2)=v(2)/max(1.d-30,vv)
         v(3)=v(3)/max(1.d-30,vv)

         dvdx=v(1)*dvxdx+v(2)*dvydx+v(3)*dvzdx
         dvdy=v(1)*dvxdy+v(2)*dvydy+v(3)*dvzdy
         dvdz=v(1)*dvxdz+v(2)*dvydz+v(3)*dvzdz
         dvd(1)=dvdx
         dvd(2)=dvdy
         dvd(3)=dvdz
c

         end

