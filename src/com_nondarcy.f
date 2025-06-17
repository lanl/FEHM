	module com_nondarcy
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
!D1  Include file for nondarcy variables
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
c gaz 110624
      real*8 phi_grad,g_term
      real*8 dg_termpi,dg_termpkb,dg_termei,dg_termekb
      real*8 daxydpi, daxydpkb, daxydei, daxydekb
c gaz 020425
      real*8 xrl_nd, xrv_nd
      parameter(vel_tol = 1.d-6, nr_iter_max = 15, visc_corr = 1.d-6)
      parameter(vel_tol_min = 1.d-22, nd_test_write = 0)
c
	end module com_nondarcy

