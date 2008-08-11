      module comzone
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
!D1 Global file for making zone 
!D1
!***********************************************************************

      integer node_azones
      integer hflag, pflag, tflag, eflag, ozflag, carbflag, carbflag2
      real*8, allocatable :: avg_values(:,:), zone_volume(:)

      type :: integer_value
         integer node_number
         type (integer_value), pointer :: nnp
      end type

      type (integer_value), pointer :: node_head, node_head_num
      type (integer_value), pointer :: node_tail, node_tail_num
      type (integer_value), pointer :: node_ptr, node_ptr_num

      type (integer_value), pointer :: zone_head
      type (integer_value), pointer :: zone_tail
      type (integer_value), pointer :: zone_ptr

      end module comzone
