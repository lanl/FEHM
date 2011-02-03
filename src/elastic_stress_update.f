      subroutine elastic_stress_update(i)
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
! Computes the stresses from strains when using 'plastic' submacro using
! the FEHM algorithm. Currently, just calls stress_3D_post
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
      real*8, dimension(1:6) :: tmp_eps, tmp_str
      
c      integer, allocatable ::   itstress(:)
      integer iws
      real*8, dimension(1:6, 1:6)  :: d_i

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    

      call stress_3D_post(i)
  
c changed by avw -- entered here by seh
!      neqp1=neq+1
!      ldna=nelm(neqp1)-neqp1
!
!      ! call compute_strains()
!      
!      ! do i=1,neq
!
!        call elastic_stiffness(i, d_i)
!        
!        alpi = alp(i)
!        bulki = bulk(i)
!        
!        tdumi = t(i) - tini(i)
!        pdumi = phi(i) - phini(i)
!        
!        tmp_eps(1) = strain_xx(i) - alpi*tdumi - bulki*pdumi
!        tmp_eps(2) = strain_yy(i) - alpi*tdumi - bulki*pdumi
!        tmp_eps(3) = strain_zz(i) - alpi*tdumi - bulki*pdumi
!        tmp_eps(4) = strain_xy(i)
!        tmp_eps(5) = strain_yz(i)
!        tmp_eps(6) = strain_zx(i)
!        
!        tmp_str = matmul(d_i, tmp_eps)
!        
!        write(iout, *) i,tmp_str(1:3)
!        
!        str_x(i) = -tmp_str(1)
!        str_y(i) = -tmp_str(2)
!        str_z(i) = -tmp_str(3)
!        str_xy(i) = -tmp_str(4)
!        str_yz(i) = -tmp_str(5)
!        str_xz(i) = -tmp_str(6)
!        
!      ! enddo
      
      return
      end


