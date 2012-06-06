      subroutine setup_edgePointers_3D()
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
! Author : Sai Rapaka 
!
      
      use comai, only: nei, neq
      use combi, only: nelm
      use comfem
      
      integer                            :: el, j, node_I, node_J
      integer, parameter                 :: numEdges = 28
      
      if(.not. allocated(edges)) then
        allocate(edges(numEdges, 2))
      endif

!      edges( 1,1) = 1; edges( 1,2) = 2
!      edges( 2,1) = 2; edges( 2,2) = 3
!      edges( 3,1) = 3; edges( 3,2) = 4
!      edges( 4,1) = 4; edges( 4,2) = 1
!      edges( 5,1) = 5; edges( 5,2) = 6
!      edges( 6,1) = 6; edges( 6,2) = 7
!      edges( 7,1) = 7; edges( 7,2) = 8
!      edges( 8,1) = 8; edges( 8,2) = 5
!      edges( 9,1) = 3; edges( 9,2) = 7
!      edges(10,1) = 6; edges(10,2) = 2
!      edges(11,1) = 4; edges(11,2) = 8
!      edges(12,1) = 5; edges(12,2) = 1
      
      edges( 1,1) = 1; edges( 1,2) = 2
      edges( 2,1) = 1; edges( 2,2) = 3
      edges( 3,1) = 1; edges( 3,2) = 4
      edges( 4,1) = 1; edges( 4,2) = 5
      edges( 5,1) = 1; edges( 5,2) = 6
      edges( 6,1) = 1; edges( 6,2) = 7
      edges( 7,1) = 1; edges( 7,2) = 8
      edges( 8,1) = 2; edges( 8,2) = 3
      edges( 9,1) = 2; edges( 9,2) = 4
      edges(10,1) = 2; edges(10,2) = 5
      edges(11,1) = 2; edges(11,2) = 6
      edges(12,1) = 2; edges(12,2) = 7
      edges(13,1) = 2; edges(13,2) = 8
      edges(14,1) = 3; edges(14,2) = 4
      edges(15,1) = 3; edges(15,2) = 5
      edges(16,1) = 3; edges(16,2) = 6
      edges(17,1) = 3; edges(17,2) = 7
      edges(18,1) = 3; edges(18,2) = 8
      edges(19,1) = 4; edges(19,2) = 5
      edges(20,1) = 4; edges(20,2) = 6
      edges(21,1) = 4; edges(21,2) = 7
      edges(22,1) = 4; edges(22,2) = 8
      edges(23,1) = 5; edges(23,2) = 6
      edges(24,1) = 5; edges(24,2) = 7
      edges(25,1) = 5; edges(25,2) = 8
      edges(26,1) = 6; edges(26,2) = 7
      edges(27,1) = 6; edges(27,2) = 8
      edges(28,1) = 7; edges(28,2) = 8

      if(.not. allocated(edgeNum1)) then
        allocate(edgeNum1(nei, numEdges))
      endif
      
      if(.not. allocated(edgeNum2)) then
        allocate(edgeNum2(nei, numEdges))
      endif
      
      if(.not. allocated(numelems)) then
        allocate(numelems(nelm(neq+1)))
      endif

      if(.not. allocated(permfactor)) then
        allocate(permfactor(nelm(neq+1), 3))
        permfactor = 1.00d0
        allocate(permtmp(nelm(neq+1), 3))
        permtmp = 1.00d0
      endif

      do el=1,nei        ! Loop over elements
        do j=1,numEdges  ! Loop over edges
          node_I = elnode(el, edges(j,1))
          node_J = elnode(el, edges(j,2))
      
          I_begin = nelm(node_I)+1
          I_end   = nelm(node_I + 1)
          do n = I_begin, I_end
            if(nelm(n).eq.node_J) then
              edgeNum1(el, j) = n
            endif
          enddo
      
          J_begin = nelm(node_J) + 1
          J_end   = nelm(node_J + 1)
          do n=J_begin, J_end
            if(nelm(n).eq.node_I) then
              edgeNum2(el, j) = n
            endif
          enddo
        enddo
      enddo
      
      numelems = 0

      do el = 1, nei
        do n=1, numEdges
          n1 = edgeNum1(el, n)
          n2 = edgeNum2(el, n) 
          numelems(n1) = numelems(n1) + 1
          numelems(n2) = numelems(n2) + 1
        enddo
      enddo

      do n=1, nelm(neq+1)
        if(numelems(n).eq.0) then
          numelems(n) = -1
        endif
      enddo

      end subroutine setup_edgePointers_3D
      
