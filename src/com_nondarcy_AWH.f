	module com_nondarcy_AWH
!     module 
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
!D1  PURPOSE
!D1
!D1  Include file for nondarcy variables added to files to eventually simplify
!D1  calling parameters
!D1
!***********************************************************************

      logical nd_flow
      real*8, allocatable :: nd_beta(:)
      real*8 d_vnd_pi
      real*8 vel_tol, vel_tol_min, v_tol
      real*8 visc_corr
      integer nr_iter_max, nd_test, nd_test_write

      real*8  muij,dmuijpi,dmuijpj,dmuijei,dmuijej
      real*8  muvij,dmuvijpi,dmuvijpj,dmuvijei,dmuvijej
      real*8  rolij,drolijpi,drolijpj,drolijei,drolijej
      real*8  rovij,drovijpi,drovijpj,drovijei,drovijej
      real*8  enlij,denlijpi,denlijpj,denlijei,denlijej
      real*8  envij,denvijpi,denvijpj,denvijei,denvijej
c gaz 131025 
      real*8  dlcpkb_nd,dlcpi_nd,dlcekb_nd,dlcei_nd  
      real*8  axyf,aexyf,acxyf,acxy_nd
      real*8  enli,enlkb,cnli,cnlkb
      real*8  dilpi,dilpkb,dilei,dilekb,dilci,dilckb
      real*8  deli,delkb,delei,delekb,delci,delckb
      real*8  dcli,dclkb,dclei,dclekb,dclci,dclckb
      real*8  dleckb_nd,dleci_nd,dlcckb_nd,dlcci_nd
      real*8  dlackb_nd,dlaci_nd
      real*8  dili,dilkb
c gaz 141025
      real*8  cnvkb,cnvi
c gaz 180925      
      real*8 daxydci, daxydckb
      real*8 fid,fid1,xrl
      real*8 pxyi,pxy
c gaz 020425
      real*8 xrl_nd, xrv_nd
c gaz 110822 debug
c gaz 120225 nd_flow liq and gas phase
      real*8 vxy_nd, dvapi_nd, dvapkb_nd, vel_nd
      real*8 axy_nd, dlapi_nd, dlapkb_nd
      real*8 aexy_nd,dlaei_nd, dlaekb_nd, dvaei_nd, dvaekb_nd
c gaz  
      real*8 dlepi_nd, dlepkb_nd, dvepi_nd, dvepkb_nd
      real*8 dleei_nd, dleekb_nd, dveei_nd, dveekb_nd
      real*8 dlei_nd, dlekb_nd, dlpi_nd, dlpkb_nd
      real*8  dvpi_nd, dvpkb_nd, dvei_nd, dvekb_nd
c gaz 021025
      real*8  dcndpi,dcndpj,dandpi,dandpj,dbndpi,dbndpj  
      real*8  dcndei,dcndej,dandei,dandej,dbndei,dbndej 
      real*8  dcndci,dcndcj,dbndci,dbndcj,dandci,dandcj
      real*8  dvelpi,dvelpj,dvelei,dvelej,dvelci,dvelcj,s_i,s_j
      real*8 velij_nd,aij,den_term
      real*8 axyd_nd,vxyd_nd,kij,kij_tol
      real*8 sx4d
c gaz 301025
      real*8 dvevpi,dvevpj,dvevei,dvevej,dvevci,dvevcj
      real*8 vxyf,vexyf,vcxyf,vexy_nd,vcxy_nd
      real*8 dvcpkb_nd,dvcpi_nd,dvcekb_nd,dvcei_nd  
      real*8 envi,envkb
      real*8 divpi,divpkb,divei,divekb,divci,divckb 
      real*8 devi,devkb,devei,devekb,devci,devckb 
      real*8 dcvi,dcvkb,dcvei,dcvekb,dcvci,dcvckb 
      real*8 dveckb_nd,dveci_nd,dvcckb_nd,dvcci_nd 
      real*8 dvackb_nd,dvaci_nd 
      real*8 divi,divkb

      parameter(vel_tol = 1.d-6, nr_iter_max = 15, visc_corr = 1.d-6)
      parameter(vel_tol_min = 1.d-22, nd_test_write = 0)
c gaz 100125
      parameter(kij_tol=1.d-24)
      end module com_nondarcy_AWH
                      