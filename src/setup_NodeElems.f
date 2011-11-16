      subroutine Setup_NodeElems()
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
! store elements connected to each node
!
! Author : Sai Rapaka

      use comai, only: nei, neq, ns
      use comfem

      implicit none

      integer, dimension(neq)        :: elcount
      integer node, nsize
      integer i_start, i_end 
      logical i_done
      integer i, j, k

      ! Currently written only for hex elements in 3-D
      ! Easy to add conditions for other elements
      elcount = 0
      do i=1,nei
        do j=1,ns
          node = elnode(i,j)
          elcount(node) = elcount(node) + 1
        enddo
      enddo

      ! Compute size of NodeElems
      nsize = neq + 1
      do i=1,neq
        nsize = nsize + elcount(i)
      enddo

      allocate(NodeElems(nsize))
      NodeElems = 0

      do i=1,neq+1
        if(i.eq.1) then
          NodeElems(i) = (neq+1)
        else
          NodeElems(i) = NodeElems(i-1) + elcount(i-1)
        endif
      enddo

      do i=1,nei
        do j=1,ns
          node = elnode(i,j)
          i_start = NodeElems(node) + 1
          i_end = NodeElems(node+1)
          i_done = .false.
          do k=i_start,i_end
            if((NodeElems(k).eq.0).and.(.not.i_done)) then
              NodeElems(k) = i
              i_done = .true.
            endif
          enddo
        enddo
      enddo

      end subroutine Setup_NodeElems

