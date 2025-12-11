      subroutine stress_2D_post_fem()
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
! 
! Integrates the stresses and strains from the gausspoints onto the nodes
! when using 'fem' computations
! 
! Author : Sai Rapaka
!

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsi
      use comfem
      implicit none

      integer i
      integer icd
      integer ii1
      integer ii2
      integer idg 
      integer iq
      integer jmi
      integer jml
      integer jmia
      integer jm
      integer neqp1
      integer ij
      integer ij1
      integer ij2
      integer iz
      integer kb
      integer neighc
      integer iau
      integer ial
      integer kz
      
      integer nmatavw
      real*8 sx1d
      real*8 axi
      real*8 ayi
      real*8 azi
      real*8 alxi
      real*8 alyi
      real*8 alzi
      real*8 avxi
      real*8 avyi
      real*8 avzi
      real*8 pvii
      real*8 phii
      real*8 swi
      real*8 dlpi
      real*8 dlpkb
      real*8 dvpi
      real*8 dvpkb
      real*8 dlei
      real*8 dvei
      real*8 dlekb
      real*8 dvekb
      real*8 dlapi
      real*8 dvapi
      real*8 dlapkb
      real*8 dvapkb
      real*8 dili
      real*8 dilkb
      real*8 divi
      real*8 divkb
      real*8 axkb
      real*8 aykb
      real*8 azkb
      real*8 alxkb
      real*8 alykb
      real*8 alzkb
      real*8 sx2c
      real*8 sx4d

      real*8 pvikb
      real*8 phikb
      real*8 radi
      real*8 radkb
      real*8 fid
      real*8 fid1
      real*8 axyd
      real*8 axy
      real*8 axyf
      real*8 heatc
      real*8 pvxy
      real*8 pxy
      real*8 pxyh
      real*8 pxyi
      real*8 sx4h
      real*8 vxyd
      real*8 vxy
      real*8 vxyf
      real*8 dpvti
      real*8 dilpi
      real*8 dilei
      real*8 divpi
      real*8 divei
      real*8 dilpkb
      real*8 divekb
      real*8 dilekb
      real*8 sx2t
      real*8 sx3t
      real*8 sxzt
      real*8 dlaei
      real*8 dlaekb
      real*8 dvaei
      real*8 dvaekb
      real*8 ti
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 redu
      real*8 divpkbction_factor
      real*8 grav_air
      real*8 grav_term

      real*8 e1i,e2i,e3i,e1kb,e2kb,e3kb
      real*8 e1bar,e2bar,e3bar
	real*8 efac,bulki,bulkkb,bulkb,alpi,alpkb,alphab
      real*8 dui,dvi,dwi,dukb,dvkb,dwkb
      real*8 xxx,xyy,xzz,xyx,xxy,xzx,yxy
      real*8 yyx,yyy,yxx,yzz,yzy,yyz,zxz
      real*8 zzx,zyz,zzy,zzz,zyy,zxx,xxz
      real*8 xx,xy,xz,yy,yx,yz,zx,zy,zz
      real*8 ddx,ddy,ddz,pxx,pxz,pyx
	real*8 pyy,pyz,pzx,pzy,pzz
	real*8 fad,fbd,fcd
      real*8 pdumi,tdumi,pdumkb,tdumkb,pdumt,tdumt,roli
	real*8 tdumx,tdumy,tdumz
	real*8 bforcex,bforcey,bforcez
	real*8 dtdumxp,dtdumxt
	real*8 dtdumyp,dtdumyt
	real*8 dtdumzp,dtdumzt
      real*8 sixsjx,siysjy,sizsjz,sixsjy
      real*8 siysjx,sixsjz,sizsjx,siysjz
      real*8 sizsjy,sjsix,sjsiy,sjsiz
	real*8 sisjx,sisjy,sisjz,termx,termy,termz
      real*8 vol,sxd,syd,tauxy,tauyz,tauzx,ctherm,biot,erat
	real*8 epi,eti,dpd,dt,shti,shpi
       real*8 fac, onedV

c      integer, allocatable ::   itstress(:)
      integer iws, j

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd,k
      integer    :: ContrElem, el
      integer    :: node(8)
      real*8 val
      real*8 :: Dmat(3,3)
      real*8 :: intsig(3), intstrain(3), cvol(neq)
      real*8, allocatable :: kmat(:)
      real*8 deltaT, deltaP, alphadT, betadP
      real*8 kmat_gp
      real*8 xdV, ydV, zdV

      if(ifem.eq.1) then

        if(icnl.ne.1) then
          write(iout,*) 'ERROR: 2D FEM Stress computations only 
     &                   enabled in XY plane'
          write(iptty,*) 'ERROR: 2D FEM Stress computations only 
     &                   enabled in XY plane'
          stop
        endif

        str_x = 0.0d0
        str_y = 0.0d0
        str_z = 0.0d0
        str_xy = 0.0d0
        str_yz = 0.0d0
        str_xz = 0.0d0

        strain_xx = 0.0d0
        strain_yy = 0.0d0
        strain_zz = 0.0d0
        strain_xy = 0.0d0
        strain_yz = 0.0d0
        strain_zx = 0.0d0

        cvol = 0.0d0

        do el=1,nei

          do k=1,ns
            node(k) = elnode(el, k)
          enddo

          intsig = 0.0d0
          intstrain = 0.d0
          onedV = 0.0d0

          do j=1,numgausspoints
            fac = detJ(el, j)*gpweight(j)
            intsig = intsig + fem_stress(el, j, :)*fac
            intstrain = intstrain + fem_strain(el, j, :)*fac
            onedV = onedV + fac
          enddo

          do k=1,ns
            cvol(node(k)) = cvol(node(k)) + onedV
            str_x(node(k)) = str_x(node(k)) + intsig(1)
            str_y(node(k)) = str_y(node(k)) + intsig(2)
            !str_z(node(k)) = str_z(node(k)) + intsig(3)
            str_xy(node(k)) = str_xy(node(k)) + intsig(3)
            !str_yz(node(k)) = str_yz(node(k)) + intsig(5)
            !str_xz(node(k)) = str_xz(node(k)) + intsig(6)

            strain_xx(node(k)) = strain_xx(node(k)) + intstrain(1)
            strain_yy(node(k)) = strain_yy(node(k)) + intstrain(2)
            !strain_zz(node(k)) = strain_zz(node(k)) + intstrain(3)
            strain_xy(node(k)) = strain_xy(node(k)) + intstrain(3)
            !strain_yz(node(k)) = strain_yz(node(k)) + intstrain(5)
            !strain_zx(node(k)) = strain_zx(node(k)) + intstrain(6)

          enddo
        enddo

        ! Make xx, yy, zz stresses negative to change to compression
        ! positive convention! 
        str_x = -str_x/cvol
        str_y = -str_y/cvol
c        str_z = -str_z/cvol
        str_xy = str_xy/cvol
c        str_yz = str_yz/cvol
c        str_xz = str_xz/cvol

        strain_xx = strain_xx/cvol
        strain_yy = strain_yy/cvol
c        strain_zz = strain_zz/cvol
        strain_xy = strain_xy/cvol
c        strain_yz = strain_yz/cvol
c        strain_zx = strain_zx/cvol
        
         if(residual_stress) then
          str_x = str_x + str_x0
          str_y = str_y + str_y0
c          str_z = str_z + str_z0
          str_xy =str_xy + str_xy0
c          str_yz =str_yz + str_yz0
c          str_xz =str_xz + str_xz0
         endif

        return
      else
        write(iout, *) 'ERROR! Stress_3D_post_fem called without 
     &     using fem macro! '
        stop
      endif

      return
      end


