      subroutine  inhstz
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To define zones for average parameter value output.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inhstz.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3 2.7 Provide Restart Capability
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai, only : inpt
      use combi, only : izonef
      use comdti, only : n0
      use comzone
      implicit none

      integer i, j, k, l, istat
      character*80 input_string
      logical null1
      
      do
         read (inpt, '(a80)') input_string
         if (null1(input_string)) exit
         read (input_string, *) i, j, k
         node_azones = node_azones + 1
         if (.not. associated(node_head_num)) then
            allocate(node_head_num, stat=istat)
            node_tail_num => node_head_num
            nullify (node_tail_num%nnp)
            node_tail_num%node_number = 0
         else
            allocate(node_tail_num%nnp, stat=istat)
            node_tail_num => node_tail_num%nnp
            nullify (node_tail_num%nnp)
            node_tail_num%node_number = 0
         end if
         if (.not. associated(zone_head)) then
            allocate(zone_head, stat=istat)
            zone_tail => zone_head
            nullify (zone_tail%nnp)
            zone_tail%node_number = 0
         else
            allocate(zone_tail%nnp, stat=istat)
            zone_tail => zone_tail%nnp
            nullify (zone_tail%nnp)
            zone_tail%node_number = 0
         end if

         if (i .gt. 0) then
            zone_tail%node_number = -node_azones
            if (j .eq. 0 .and. k .eq. 0) then
               j = n0
               k = 1
            end if
            do l = i, j, k
               if (.not. associated(node_head)) then
                  allocate(node_head, stat=istat)
                  node_tail => node_head
                  nullify (node_tail%nnp)
                  node_tail%node_number = l
               else
                  allocate(node_tail%nnp, stat=istat)
                  node_tail => node_tail%nnp
                  nullify (node_tail%nnp)
                  node_tail%node_number = l
               end if
               node_tail_num%node_number = 
     &              node_tail_num%node_number + 1
            end do
         else
            zone_tail%node_number = -i
            do l = 1, n0
               if (izonef(l) .eq. abs(i)) then
                  if (.not. associated(node_head)) then
                     allocate(node_head, stat=istat)
                     node_tail => node_head
                     nullify (node_tail%nnp)
                     node_tail%node_number = l
                  else
                     allocate(node_tail%nnp, stat=istat)
                     node_tail => node_tail%nnp
                     nullify (node_tail%nnp)
                     node_tail%node_number = l
                  end if
                  node_tail_num%node_number = 
     &                 node_tail_num%node_number + 1
               end if
            end do
         end if
      end do  
    
      end
