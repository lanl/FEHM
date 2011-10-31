      subroutine compute_average_stress()
      
      use comai, only: nei, neq
      use combi, only: nelm
      use comdi, only: phi
      use comsi, only: ihms, bulk, e2, e3
      use comfem
      
      implicit none
      
      integer                            :: el, j, n1, n2, n
      integer                            :: node_I, node_J

      real*8, dimension(6)               :: intsig
      real*8                             :: onedV, fac, biot
      real*8                             :: mean_stress, temp, P_ij
      real*8                             :: bulki, bulkj, bulkbar
      real*8                             :: harmonic_mean

      integer  counter_nodal(neq)
      
      if(.not.allocated(permFactor_nodal)) then
         allocate(permFactor_nodal(neq))
      endif

      if(.not. allocated(avg_stress)) then
        allocate(avg_stress(nei, 6))
      endif
      
      if(.not. allocated(edgeNum1)) then
        call setup_edgePointers_3D()
      endif

      permFactor = 0.0d0
      permFactor_nodal = 0.0d0
      counter_nodal = 0

      do el=1,nei
      
        avg_stress(el,:) = 0.0d0
        intsig = 0.0d0
        onedV = 0.0d0
      
        do j=1,numgausspoints
          fac = detJ(el,j)*gpweight(j)
          intsig = intsig + fem_stress(el,j,:)*fac
          onedV = onedV + fac
        enddo
      
        avg_stress(el,:) = intsig/onedV
        mean_stress = sum(avg_stress(el, 1:3))/3.0d0

        do n=1, 12
          node_I = elnode(el, edges(n, 1))
          node_J = elnode(el, edges(n, 2))

          P_ij = 0.5d0*(phi(node_I) + phi(node_J))
          biot = 0.5d0*(bulk(node_I) + bulk(node_J))

          bulki = 3.0d0*e2(node_I) + 2.0d0*e3(node_I)
          bulkj = 3.0d0*e2(node_J) + 2.0d0*e3(node_J)
          bulkbar = harmonic_mean(bulki, bulkj)

          temp = mean_stress + bulkbar*biot*P_ij

          if(temp.ge.0.0d0) then
            n1 = edgeNum1(el, n)
            n2 = edgeNum2(el, n)
            fac = min(1000.0d0, 1.0d0 + temp*99.0d0)
            permFactor(n1) = permFactor(n1) + fac
            permFactor(n2) = permFactor(n2) + fac
          else
            n1 = edgeNum1(el, n)
            n2 = edgeNum2(el, n)
            permFactor(n1) = permFactor(n1) + 1.0d0
            permFactor(n2) = permFactor(n2) + 1.0d0
          endif

          permfactor_nodal(node_I) = permfactor_nodal(node_I) + fac
          counter_nodal(node_I) = counter_nodal(node_I) + 1
          permfactor_nodal(node_J) = permfactor_nodal(node_J) + fac
          counter_nodal(node_J) = counter_nodal(node_J) + 1

        enddo
      enddo
      
      permFactor = permFactor/numElems

      do j=1,neq
         if(counter_nodal(j).le.0) then
            write(*,*)'in comupte_average_stress'
            write(*,*)'counter_nodal=0 for j=', j
            stop
         endif
      enddo

      permfactor_nodal = permfactor_nodal/counter_nodal

      open(38, file='avg_stress.dat', status='unknown')
      do el = 1,nei
        write(38, *) el, avg_stress(el, 1:3)
      enddo

      do n=1, nelm(neq+1)
        if(numElems(n).gt.0) then
          write(38,*) n, permFactor(n)
        endif
      enddo  
    
      close(38)
      
      end subroutine compute_average_stress
