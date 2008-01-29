      subroutine area2d_tri(ntri,
     *                      x1,y1,z1,x2,y2,z2,x3,y3,z3,
     *                      xarea_tri,xvol1,xvol2,xvol3,
     *                      xarea12,xarea23,xarea31)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1      This routine calculates the total area of a triangle and
CD1         also the three Voronoi areas along each of the three
CD1         sides. The vertices of the triangle are assumed to be
CD1         ordered counter-clockwise.
CD1
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2  $Log:   /pvcs.config/fehm90/src/area2d_tri.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Mar 31 08:29:28 1997   gaz
CD2 harolds volumes and areas
CD2 
CD2    Rev 1.1   Thu Jun 27 11:29:02 1996   gaz
CD2 added prolog
CD2 
CD2    Rev 1.0   Fri Apr 26 14:57:56 1996   gaz(H. Trease)
CD2 Initial revision.
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
C***********************************************************************
C
C     INPUT ARGUMENTS -
C
C        ntri - Number of triangles
C        x1(1:ntri) - The X-coordinates of all the 1st vertices.
C        y1(1:ntri) - The Y-coordinates of all the 1st vertices.
C        z2(1:ntri) - The Z-coordinates of all the 1st vertices.
C        x2(1:ntri) - The X-coordinates of all the 2nd vertices.
C        y2(1:ntri) - The Y-coordinates of all the 2nd vertices.
C        z2(1:ntri) - The Z-coordinates of all the 2nd vertices.
C        x3(1:ntri) - The X-coordinates of all the 3rd vertices.
C        y3(1:ntri) - The Y-coordinates of all the 3rd vertices.
C        z3(1:ntri) - The Z-coordinates of all the 3rd vertices.
C
C     OUTPUT ARGUMENTS -
C
C        xarea_tri(1:ntri) - The total area of the triangles.
C        xvol1(1:ntri)     - The area associated with vertex 1.
C        xvol2(1:ntri)     - The area associated with vertex 2.
C        xvol3(1:ntri)     - The area associated with vertex 3.
C        xarea12(1:ntri) - The area of the Voronoi edge associated
C                              with the 1-2 connection divided by
C                              the length of the 1-2 connection.
C        xarea23(1:ntri) - The area of the Voronoi edge associated
C                              with the 2-3 connection divided by
C                              the length of the 2-3 connection.
C        xarea31(1:ntri) - The area of the Voronoi edge associated
C                              with the 3-1 connection divided by
C                              the length of the 3-1 connection.
C
C
C#######################################################################
C
C
C     integer npts
C
      implicit none
      integer ntri
      real*8 crosx, crosy, crosz, dotpr, vecmag, ds13, ds23
      real*8 xa11, ya11, za11, xa12, ya12, za12, xa21, ya21, za21
      real*8 xa22, ya22, za22, xa31, ya31, za31, xa32, ya32, za32
      real*8 a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, a123x
      real*8 a123y, a123z, area_vor, area_tri, area_error, a1_sign
      real*8 a2_sign, a3_sign, area1, area2, area3
      real*8 x1(ntri), y1(ntri), z1(ntri)
      real*8 x2(ntri), y2(ntri), z2(ntri)
      real*8 x3(ntri), y3(ntri), z3(ntri)
C
      real*8 xarea_tri(ntri),xarea12(ntri),xarea23(ntri),xarea31(ntri)
      real*8 xvol1(ntri), xvol2(ntri), xvol3(ntri)
C
      integer i1
C
      real*8 a, b, c, d, e, f
      real*8 xa, ya, za
      real*8 xb, yb, zb
      real*8 xd, yd, zd
      real*8 xm, ym, zm
      real*8 xn, yn, zn, rn
      real*8 xn1, yn1, zn1
      real*8 xv, yv, zv
      real*8 xl, yl, zl
      real*8 ql
      real*8 x12, y12, z12
      real*8 x23, y23, z23
      real*8 x13, y13, z13
      real*8 dot3, dotb3, xdot, xdot1, xdot2, xdot3
      real*8 xarea1, xarea2, xarea3
      real*8 rb3
      real*8 ds1, ds2, ds3, ds12
      real*8 cmx, cmy, cmz
      real*8 cvx, cvy, cvz
C
C#######################################################################
c
c     Statement functions for cross product, dot product and vector magnitude
c
      crosx(a,b,c,d,e,f)=b*f-c*e
      crosy(a,b,c,d,e,f)=c*d-a*f
      crosz(a,b,c,d,e,f)=a*e-b*d
      dotpr(a,b,c,d,e,f)=a*d + b*e + c*f
      vecmag(a,b,c) = 0.5d0*sqrt(a**2 + b**2 + c**2)
c
C
C#######################################################################
C
C
      do i1=1,ntri
         xa=x1(i1)
         ya=y1(i1)
         za=z1(i1)
         xb=(x2(i1)-xa)
         yb=(y2(i1)-ya)
         zb=(z2(i1)-za)
         xd=(x3(i1)-xa)
         yd=(y3(i1)-ya)
         zd=(z3(i1)-za)
         xm=(x1(i1)+x2(i1)+x3(i1))/3.0d+00
         ym=(y1(i1)+y2(i1)+y3(i1))/3.0d+00
         zm=(z1(i1)+z2(i1)+z3(i1))/3.0d+00
C
C        ***************************************************************
C        Calculate the Voronoi point
C
         xn1=crosx(xb,yb,zb,xd,yd,zd)
         yn1=crosy(xb,yb,zb,xd,yd,zd)
         zn1=crosz(xb,yb,zb,xd,yd,zd)
         xn=crosx(xb,yb,zb,xn1,yn1,zn1)
         yn=crosy(xb,yb,zb,xn1,yn1,zn1)
         zn=crosz(xb,yb,zb,xn1,yn1,zn1)
         rn=1.0/sqrt(xn*xn+yn*yn+zn*zn)
         xn=xn*rn
         yn=yn*rn
         zn=zn*rn
         dotb3=xb*xd+yb*yd+zb*zd
         dot3=dotb3/(xd*xd+yd*yd+zd*zd)
         rb3=1.0/(xb*xb+yb*yb+zb*zb)
         ql=(1.0-dot3)/(1.0-dot3*dotb3*rb3+1.0d-30)
         xl=0.5*(ql*(xd-dotb3*rb3*xb)+xb)
         yl=0.5*(ql*(yd-dotb3*rb3*yb)+yb)
         zl=0.5*(ql*(zd-dotb3*rb3*zb)+zb)
         ds1=sqrt((xl)**2+(yl)**2+(zl)**2)
         ds2=sqrt((xl-xb)**2+(yl-yb)**2+(zl-zb)**2)
         ds3=sqrt((xl-xd)**2+(yl-yd)**2+(zl-zd)**2)
         xv=xl+xa
         yv=yl+ya
         zv=zl+za
C
C        ***************************************************************
C
         x12=0.5*(x1(i1)+x2(i1))
         y12=0.5*(y1(i1)+y2(i1))
         z12=0.5*(z1(i1)+z2(i1))
         x13=0.5*(x1(i1)+x3(i1))
         y13=0.5*(y1(i1)+y3(i1))
         z13=0.5*(z1(i1)+z3(i1))
         x23=0.5*(x2(i1)+x3(i1))
         y23=0.5*(y2(i1)+y3(i1))
         z23=0.5*(z2(i1)+z3(i1))
         xdot1=(x2(i1)-x1(i1))*(x3(i1)-x1(i1))+
     *         (y2(i1)-y1(i1))*(y3(i1)-y1(i1))+
     *         (z2(i1)-z1(i1))*(y3(i1)-y1(i1))
         xdot2=(x1(i1)-x3(i1))*(x2(i1)-x3(i1))+
     *         (y1(i1)-y3(i1))*(y2(i1)-y3(i1))+
     *         (z1(i1)-z3(i1))*(y2(i1)-y3(i1))
         xdot3=(x3(i1)-x2(i1))*(x1(i1)-x2(i1))+
     *         (y3(i1)-y2(i1))*(y1(i1)-y2(i1))+
     *         (z3(i1)-z2(i1))*(y1(i1)-y2(i1))
         xarea1=sqrt((xv-x12)*(xv-x12)+
     *               (yv-y12)*(yv-y12)+
     *               (zv-z12)*(zv-z12))
         xarea2=sqrt((xv-x13)*(xv-x13)+
     *               (yv-y13)*(yv-y13)+
     *               (zv-z13)*(zv-z13))
         xarea3=sqrt((xv-x23)*(xv-x23)+
     *               (yv-y23)*(yv-y23)+
     *               (zv-z23)*(zv-z23))
         ds12=sqrt((x2(i1)-x1(i1))**2+
     *             (y2(i1)-y1(i1))**2+
     *             (z2(i1)-z1(i1))**2)
         cmx=  (y2(i1)-y1(i1))*(zm-z1(i1))-(ym-y1(i1))*(z2(i1)-z1(i1))
         cmy=-((x2(i1)-x1(i1))*(zm-z1(i1))-(xm-x1(i1))*(z2(i1)-z1(i1)))
         cmz=  (x2(i1)-x1(i1))*(ym-y1(i1))-(xm-x1(i1))*(y2(i1)-y1(i1))
         cvx=  (y2(i1)-y1(i1))*(zv-z1(i1))-(yv-y1(i1))*(z2(i1)-z1(i1))
         cvy=-((x2(i1)-x1(i1))*(zv-z1(i1))-(xv-x1(i1))*(z2(i1)-z1(i1)))
         cvz=  (x2(i1)-x1(i1))*(yv-y1(i1))-(xv-x1(i1))*(y2(i1)-y1(i1))
         xdot=cmx*cvx+cmy*cvy+cmz*cvz
         if(xdot.lt.0.0) then
            xarea12(i1)=-xarea1/ds12
         else
            xarea12(i1)= xarea1/ds12
         endif
         ds13=sqrt((x3(i1)-x1(i1))**2+
     *             (y3(i1)-y1(i1))**2+
     *             (z3(i1)-z1(i1))**2)
         cmx=  (y1(i1)-y3(i1))*(zm-z3(i1))-(ym-y3(i1))*(z1(i1)-z3(i1))
         cmy=-((x1(i1)-x3(i1))*(zm-z3(i1))-(xm-x3(i1))*(z1(i1)-z3(i1)))
         cmz=  (x1(i1)-x3(i1))*(ym-y3(i1))-(xm-x3(i1))*(y1(i1)-y3(i1))
         cvx=  (y1(i1)-y3(i1))*(zv-z3(i1))-(yv-y3(i1))*(z1(i1)-z3(i1))
         cvy=-((x1(i1)-x3(i1))*(zv-z3(i1))-(xv-x3(i1))*(z1(i1)-z3(i1)))
         cvz=  (x1(i1)-x3(i1))*(yv-y3(i1))-(xv-x3(i1))*(y1(i1)-y3(i1))
         xdot=cmx*cvx+cmy*cvy+cmz*cvz
         if(xdot.lt.0.0) then
            xarea31(i1)=-xarea2/ds13
         else
            xarea31(i1)= xarea2/ds13
         endif
         ds23=sqrt((x3(i1)-x2(i1))**2+
     *             (y3(i1)-y2(i1))**2+
     *             (z3(i1)-z2(i1))**2)
         cmx=  (y3(i1)-y2(i1))*(zm-z2(i1))-(ym-y2(i1))*(z3(i1)-z2(i1))
         cmy=-((x3(i1)-x2(i1))*(zm-z2(i1))-(xm-x2(i1))*(z3(i1)-z2(i1)))
         cmz=  (x3(i1)-x2(i1))*(ym-y2(i1))-(xm-x2(i1))*(y3(i1)-y2(i1))
         cvx=  (y3(i1)-y2(i1))*(zv-z2(i1))-(yv-y2(i1))*(z3(i1)-z2(i1))
         cvy=-((x3(i1)-x2(i1))*(zv-z2(i1))-(xv-x2(i1))*(z3(i1)-z2(i1)))
         cvz=  (x3(i1)-x2(i1))*(yv-y2(i1))-(xv-x2(i1))*(y3(i1)-y2(i1))
         xdot=cmx*cvx+cmy*cvy+cmz*cvz
         if(xdot.lt.0.0) then
            xarea23(i1)=-xarea3/ds23
         else
            xarea23(i1)= xarea3/ds23
         endif
C
         xarea_tri(i1)=0.5d+00*sqrt(xn1**2+yn1**2+zn1**2)
C
         xa11= crosx((x12-x1(i1)),(y12-y1(i1)),(z12-z1(i1)),
     *               (xv -x1(i1)),(yv -y1(i1)),(zv -z1(i1)))
         ya11= crosy((x12-x1(i1)),(y12-y1(i1)),(z12-z1(i1)),
     *               ( xv-x1(i1)),( yv-y1(i1)),( zv-z1(i1)))
         za11= crosz((x12-x1(i1)),(y12-y1(i1)),(z12-z1(i1)),
     *               ( xv-x1(i1)),( yv-y1(i1)),( zv-z1(i1)))
         xa12=-crosx((x13-x1(i1)),(y13-y1(i1)),(z13-z1(i1)),
     *               ( xv-x1(i1)),( yv-y1(i1)),( zv-z1(i1)))
         ya12=-crosy((x13-x1(i1)),(y13-y1(i1)),(z13-z1(i1)),
     *               ( xv-x1(i1)),( yv-y1(i1)),( zv-z1(i1)))
         za12=-crosz((x13-x1(i1)),(y13-y1(i1)),(z13-z1(i1)),
     *               ( xv-x1(i1)),( yv-y1(i1)),( zv-z1(i1)))
         xa21= crosx((x23-x2(i1)),(y23-y2(i1)),(z23-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         ya21= crosy((x23-x2(i1)),(y23-y2(i1)),(z23-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         za21= crosz((x23-x2(i1)),(y23-y2(i1)),(z23-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         xa22=-crosx((x12-x2(i1)),(y12-y2(i1)),(z12-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         ya22=-crosy((x12-x2(i1)),(y12-y2(i1)),(z12-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         za22=-crosz((x12-x2(i1)),(y12-y2(i1)),(z12-z2(i1)),
     *               ( xv-x2(i1)),( yv-y2(i1)),( zv-z2(i1)))
         xa31= crosx((x13-x3(i1)),(y13-y3(i1)),(z13-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
         ya31= crosy((x13-x3(i1)),(y13-y3(i1)),(z13-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
         za31= crosz((x13-x3(i1)),(y13-y3(i1)),(z13-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
         xa32=-crosx((x23-x3(i1)),(y23-y3(i1)),(z23-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
         ya32=-crosy((x23-x3(i1)),(y23-y3(i1)),(z23-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
         za32=-crosz((x23-x3(i1)),(y23-y3(i1)),(z23-z3(i1)),
     *               ( xv-x3(i1)),( yv-y3(i1)),( zv-z3(i1)))
c
c     Area vector of triangle
c
      xa= crosx((x2(i1)-x1(i1)),(y2(i1)-y1(i1)),
     *    (z2(i1)-z1(i1)),(x3(i1)-x1(i1)),(y3(i1)-y1(i1)),
     *    (z3(i1)-z1(i1)))
      ya= crosy((x2(i1)-x1(i1)),(y2(i1)-y1(i1)),
     *    (z2(i1)-z1(i1)),(x3(i1)-x1(i1)),(y3(i1)-y1(i1)),
     *    (z3(i1)-z1(i1)))
      za= crosz((x2(i1)-x1(i1)),(y2(i1)-y1(i1)),
     *    (z2(i1)-z1(i1)),(x3(i1)-x1(i1)),(y3(i1)-y1(i1)),
     *    (z3(i1)-z1(i1)))
c
c     Vector sum of Voronoi area vector for each of three nodes in the triangle
c
      a1x = xa11 + xa12
      a1y = ya11 + ya12
      a1z = za11 + za12
      a2x = xa21 + xa22
      a2y = ya21 + ya22
      a2z = za21 + za22
      a3x = xa31 + xa32
      a3y = ya31 + ya32
      a3z = za31 + za32
c
c     Total area vector of triangle formed by summing 6 voronoi contributions
c
      a123x = a1x + a2x + a3x
      a123y = a1y + a2y + a3y
      a123z = a1z + a2z + a3z
c
c     Comparison of area calculated two different ways
c
      area_vor  =0.5*sqrt(a123x**2+a123y**2+a123z**2)
      area_tri  =0.5*sqrt(xa**2+ya**2+za**2)
      area_error=0.5*sqrt((a123x-xa)**2+(a123y-ya)**2+(a123z-za)**2)
      
      if(area_error .gt. 1.e-9*area_tri)then
        print *,'Error Calculating Voronoi Area'
        print *,area_vor,area_tri,area_error
      endif
c
c     Dot product of voronoi area vector with triangle area vector. This
c     quantity will be positive for positive area contributions and negative
c     for negative area contributions
c
      a1_sign = dotpr(a1x,a1y,a1z,xa,ya,za)
      a2_sign = dotpr(a2x,a2y,a2z,xa,ya,za)
      a3_sign = dotpr(a3x,a3y,a3z,xa,ya,za)
c
c     Find the absolute magnitude of the voronoi contributions
c 
      area1   = vecmag(a1x,a1y,a1z)
      area2   = vecmag(a2x,a2y,a2z)
      area3   = vecmag(a3x,a3y,a3z)
c
c     Give the area contributions their correct sign
c      
      area1 = sign(area1,a1_sign)
      area2 = sign(area2,a2_sign)
      area3 = sign(area3,a3_sign)
c
c load areas into arrays
c
      xvol1(i1) = area1
      xvol2(i1) = area2
      xvol3(i1) = area3
      enddo
C
      goto 9999
9999  continue
C
      return
      end
