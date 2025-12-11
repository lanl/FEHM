      module comriv 
!     comriv
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
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
!D1
!D1 PURPOSE
!D1
!D1 Global file for array variables and pointers river and well modules
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 13-Dec-06, Programmer: R. Pawar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comriv.f_a  $
!D2
!***********************************************************************

      integer npoint_riv,iriver,maxlay,nbc_riv,npoint_riv_old,nic_old
      integer maxriver_nodes,n_well_prod,n_well_seg, nodes_well2_added
      integer iwell_phy_all, max_seg_div,nnelm_riv
      integer, allocatable :: isriverf(:) 
      integer, allocatable :: ifriverf(:)
      integer, allocatable :: iriverf(:)
      integer, allocatable :: inriverf(:)
      integer, allocatable :: iwsp(:)
      integer, allocatable :: isriver(:) 
      integer, allocatable :: ibc_river(:) 
      integer, allocatable :: iriver_list(:)
      integer, allocatable :: iriver_con_node(:) 	 
      integer, allocatable :: izone_river(:) 
      integer, allocatable :: izone_river_nodes(:) 
      integer, allocatable :: river_nodes_local(:)
      integer, allocatable :: river_nodes_global(:)
      integer, allocatable :: iriver_out_node(:)
      integer, allocatable :: iriver_first_node(:)
      integer, allocatable :: mdnodes_riv(:) 
      integer, allocatable :: mdnode_riv(:,:)
      integer, allocatable :: rivbegin(:)
      integer, allocatable :: rivend(:) 

      real*8, allocatable ::  wgt_river(:)
      real*8, allocatable ::  bc_river(:,:)
      real*8, allocatable ::  coor_riv(:,:)
      real*8, allocatable  :: coor_dum(:,:)
      real*8, allocatable ::  area_riv(:,:)
      real*8, allocatable ::  csg_th(:,:)
      real*8, allocatable ::  perm_riv(:,:)  
      real*8, allocatable ::  river01(:)
      real*8, allocatable ::  river02(:)
      real*8, allocatable ::  river03(:)
      real*8, allocatable  :: vol1(:)
      real*8, allocatable ::  mod_dis(:,:)
c arrays for welltype = 2
	integer, allocatable :: iwell_geom(:)
	integer, allocatable :: iwell_phys(:) 
	integer, allocatable :: iwell_prod(:) 
	integer, allocatable :: iwell_seg(:,:) 
	integer, allocatable :: nwell2_prim(:) 
	integer, allocatable :: new_node_well2_segid(:) 
	integer, allocatable :: new_node_well2(:) 
	integer, allocatable :: neigh_well2(:,:)
	integer, allocatable :: neigh_well2_count(:) 
	integer, allocatable :: neigh_well2_new(:,:) 
	integer, allocatable :: nwell2_int(:) 
	integer, allocatable :: iwell2_int(:,:) 
	integer, allocatable :: izone_well2(:)		
	integer, allocatable :: izlabelp(:)
	integer, allocatable :: izlabels(:) 
	integer, allocatable :: nelm_riv(:,:)	
	integer, allocatable :: idir_seg(:) 	
	integer, allocatable :: iwdum(:) 	
	integer, allocatable :: iwell_end(:) 			 						
	
	real*8, allocatable ::  coor_well2(:,:)	
	real*8, allocatable ::  coor_new_well2(:,:)	
	real*8, allocatable ::  well_rad(:)	
	real*8, allocatable ::  well_dz(:)
	real*8, allocatable ::  well_connect_fac(:)
	real*8, allocatable ::  delxw2(:)
	real*8, allocatable ::  delyw2(:)
	real*8, allocatable ::  delzw2(:)
	
	real*8, allocatable ::  sx_w(:,:)				 
      real*8, allocatable ::  seg_cos(:,:)
      
      parameter (maxriver_nodes = 100000)

      end module comriv
