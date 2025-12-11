      module comwt
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
***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global file for wtsi
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/comwt.f_a  $
!D2
!***********************************************************************

      integer n_wt_cols, ic, move_wt
      integer wtsi_isot, column_read, isconwt
      integer, allocatable :: dry_zone(:), move_type(:)
      integer, allocatable :: wcol(:), col(:,:), n_col(:)
      integer, allocatable :: col_out(:)
      integer, allocatable :: node_above(:)
      integer, allocatable :: pointer_above(:)

      integer wt_flag, iad_up_wtsi, irich
      real*8  wt_elev

      real*8 sattol, zfac_ani, head_ck_first, head_id
c      parameter (sattol = 1.d-1)

      real*8 dry_tol, rlptol
      parameter (dry_tol=1.d-10)

c move_type =0 do not move source, 1 move source with water table, < 0 
c  move source within well to uppermost saturated node
c
	end module comwt
