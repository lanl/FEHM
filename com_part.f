      module com_part 
!     com_part
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
!D1 Global include file for array variables for partitioning
!D1 (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 21-MAY-03    P. Fasel       22      Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/com_part.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 N/A
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 None
!D4 
!***********************************************************************
!D5
!D5 GLOBAL OBJECTS
!D5
!D5 Global Constants
!D5
!D5   None
!D5
!D5 Global Types
!D5
!D5   None
!D5
!D5 Global Variables
!D5
!D5                            COMMON
!D5   Identifier      Type     Block  Description
!D5
!D5   ***** COMMON Block fbb pointers and associated variables *****
!D5   nelm_part       INT      fbb    Node connectivity for a partition
!D5   bound_part      INT      fbb    Boundary nodes of partition
!D5   bound_conn_part INT      fbb    Connection nodes to boundary nodes
!D5   index_part      INT      fbb    New node index to old node number
!D5   rindex_part     INT      fbb    Old node index to new node number
!D5
!D5 Global Subprograms
!D5
!D5   None
!D5
!***********************************************************************

      integer nzone_para,nzone_para0,npart_type,npall,igrow
      integer, allocatable ::  nelm_part(:,:)
      integer, allocatable ::  nelmdg_part(:,:)

      integer, allocatable ::  bound_part(:,:)
      integer, allocatable ::  bound_istrw(:,:)

      integer, allocatable ::  index_part(:,:)
      integer, allocatable ::  rindex_part(:,:)

      integer, allocatable ::  neq_part(:)
      integer, allocatable ::  iad_part(:)
      integer, allocatable ::  itert_part(:)
      integer, allocatable ::  itotal_part(:)
      integer, allocatable ::  itotals_part(:)
      integer, allocatable ::  bound_pos(:)
      integer, allocatable ::  bound_zone(:)	


      integer, allocatable ::  izone_para(:)
      integer, allocatable ::  igrow_part(:)
	integer, allocatable ::  ipart_ini(:)
      integer, allocatable ::  izone_para_temp(:)
      integer, allocatable ::  igrow_part_temp(:)
      integer, allocatable ::  ipart_ini_temp(:)
      integer, allocatable ::  index_part_temp(:,:)
      integer, allocatable ::  rindex_part_temp(:,:)
      integer, allocatable ::  zone_level(:)

      integer, allocatable ::  decomp_part(:)
      integer, allocatable ::  decomp_part_temp(:)
	integer, allocatable ::  decomp_part_temp1(:)

      integer, allocatable ::  nmat_part(:,:)
      integer, allocatable ::  nb_part(:,:)
      integer, allocatable ::  nrhs_part(:,:)

      integer, allocatable ::  npvt_part(:,:)
      integer, allocatable ::  istrw_part(:,:)

      real*8, allocatable ::  bp_part(:,:)
      real*8, allocatable ::  fdum_part(:)
      real*8, allocatable ::  f0_part(:)
      real*8, allocatable ::  timing_part(:)
      real*8, allocatable ::  repeat_part(:)
	real*8, allocatable ::  flux_l_part(:,:)
	real*8, allocatable ::  flux_v_part(:,:)
      

      end module com_part

