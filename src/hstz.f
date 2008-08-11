      subroutine  hstz
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
!D1 To determine average parameter values for a zone
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/hstz.f_a  $
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

      use combi, only : sx1
      use comco2
      use comdi, only : head, t, phi, qh, pnx, ifree, rlxyf
      use comwt, only : rlptol, sattol, wt_flag, wt_elev, head_id
      use comzone
      implicit none

      integer i, ij, j, k, count, num_zone
      integer, allocatable :: zn(:)
      real*8 headdum
      real*8, allocatable :: zv(:)

      allocate (zv(node_azones), zn(node_azones))
!      zv = zone_volume
      zv = 0.
      if (.not. allocated(avg_values)) 
     &     allocate(avg_values(node_azones,ozflag))
      node_ptr_num => node_head_num
      node_ptr => node_head
      avg_values = 0
      num_zone = 0
      outer: do
         if (.not. associated(node_ptr_num)) exit
         num_zone = num_zone + 1
! Number of nodes in current zone
         k = node_ptr_num%node_number
         count = 0
         if (k .ne. 0 ) then
            sum_values: do 
               if (.not. associated(node_ptr)) exit
! Number of current node in the zone
               i = node_ptr%node_number
               zn(num_zone) = i
! Sum values (weighted by node volume)
! head
               if (hflag .ne. 0) then
                  call headctr(4,i,phi(i),headdum)
c gaz 7-22-05
                  if (ifree .ne. 0) then
                     if(rlxyf(i).lt.rlptol+sattol) headdum = head_id
                  end if
                  if(headdum.gt.0.) then
                     avg_values(num_zone,hflag) = sx1(i)*pnx(i) * 
     &                    headdum + avg_values(num_zone,hflag)
                     zv(num_zone)=zv(num_zone)+sx1(i)*pnx(i)
!                  else
!                     zv(num_zone)=zv(num_zone)-sx1(i)*pnx(i)
                  endif
               end if
! Added permeability weighting for p,t,e for consistency with
! zone volume calculation (zvd 05/02/2007)
! pressure
               if (pflag .ne. 0) then
                  avg_values(num_zone,pflag) = sx1(i)*pnx(i) * phi(i) + 
     &                 avg_values(num_zone,pflag)
               end if
! temperature
               if (tflag .ne. 0) then
                  avg_values(num_zone,tflag) = sx1(i)*pnx(i) * t(i) + 
     &                 avg_values(num_zone,tflag)
               end if
! enthalpy
               if (eflag .ne. 0) then
                  avg_values(num_zone,eflag) = sx1(i)*pnx(i) * qh(i) + 
     &                 avg_values(num_zone,eflag)
               end if
c     RJP added below for outputing total carbon mass in different zones
               if (carbflag.ne.0) then
                  avg_values(num_zone,carbflag) = 
     &                 avg_values(num_zone,carbflag)+denco2h(i)*sx1(i)
               endif
c     RJP 08/09/07 added below
               if (carbflag2.ne.0) then
                  avg_values(num_zone,carbflag2) = 
     &                 avg_values(num_zone,carbflag2)+fl(i)
               endif
               node_ptr =>node_ptr%nnp
               count = count + 1
               if (count .eq. k) exit
            end do sum_values
         end if
         node_ptr_num =>node_ptr_num%nnp
      end do outer
      do i = 1, num_zone
         do j = 1, ozflag
            if(j.ne.carbflag) then
               if (zone_volume(i) .eq. 0.) then
! No nodes in zone
                  avg_values(i,j) = -9.99e6
               else
                  if (j .ne. hflag ) then
! parameter being averaged is not head
                     avg_values(i,j) = avg_values(i,j)/zone_volume(i)
                  else
! head is being averaged
                     if (zv(i) .ne. 0.) then
                        avg_values(i,j) = avg_values(i,j)/zv(i)
                     else
                        ij = -1*zn(i)
                        if (wt_flag .ne. 0 .and. ij .ne. 0) then
! Find water_table in column
                           call wtsictr(ij)
                           avg_values(i,j) = -1.*wt_elev
                        else
! All nodes in zone dry
                           avg_values(i,j) = -1.11e6
                        end if
                     end if
                  end if
               end if
            end if
         end do
      end do
      deallocate (zv, zn)

      end
