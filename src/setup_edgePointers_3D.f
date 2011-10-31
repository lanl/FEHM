      subroutine setup_edgePointers_3D()
      
      use comai, only: nei, neq
      use combi, only: nelm
      use comfem
      
      integer                            :: el, j, node_I, node_J
      
      if(.not. allocated(edges)) then
        allocate(edges(12, 2))
      endif

      edges( 1,1) = 1; edges( 1,2) = 2
      edges( 2,1) = 2; edges( 2,2) = 3
      edges( 3,1) = 3; edges( 3,2) = 4
      edges( 4,1) = 4; edges( 4,2) = 1
      edges( 5,1) = 5; edges( 5,2) = 6
      edges( 6,1) = 6; edges( 6,2) = 7
      edges( 7,1) = 7; edges( 7,2) = 8
      edges( 8,1) = 8; edges( 8,2) = 5
      edges( 9,1) = 3; edges( 9,2) = 7
      edges(10,1) = 6; edges(10,2) = 2
      edges(11,1) = 4; edges(11,2) = 8
      edges(12,1) = 5; edges(12,2) = 1
      
      if(.not. allocated(edgeNum1)) then
        allocate(edgeNum1(nei, 12))
      endif
      
      if(.not. allocated(edgeNum2)) then
        allocate(edgeNum2(nei, 12))
      endif
      
      if(.not. allocated(numElems)) then
        allocate(numElems(nelm(neq+1)))
      endif

      if(.not. allocated(permFactor)) then
        allocate(permFactor(nelm(neq+1)))
      endif

      if(.not. allocated(permFactor_nodal)) then
        allocate(permFactor_nodal(neq))
      endif

      do el=1,nei        ! Loop over elements
        do j=1,12        ! Loop over edges
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
      
      numElems = 0

      do el = 1, nei
        do n=1, 12
          n1 = edgeNum1(el, n)
          n2 = edgeNum2(el, n) 
          numElems(n1) = numElems(n1) + 1
          numElems(n2) = numElems(n2) + 1
        enddo
      enddo

      do n=1, nelm(neq+1)
        if(numElems(n).eq.0) then
          numElems(n) = -1
        endif
      enddo

      open(37, file='edgeNumbers.dat', status='unknown')
      write(37, *) 'Edge Number information :'
      
      do el=1, nei
        write(37, '(13I8)') el, edgeNum1(el,:)
      enddo
 
      write(37,*) 
      write(37,*) 

      do n=1, nelm(neq+1)
        write(37, *) n, numElems(n)
      enddo

      close(37)

      end subroutine setup_edgePointers_3D
      
