      subroutine fem_elastic_stiffness(i, j, D)
      
      use comai, only: iout, iptty, ns
      use comsi, only: modelNumber, plasticModel, elastic_mod, poisson
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6,6)      :: D

      integer                      :: k, itmp, iModel, node
      real*8                       :: e1bar, e2bar, e3bar
      real*8                       :: E, nu, lambda, G

      itmp = modelNumber(elnode(i, 1))
      iModel = plasticModel(itmp)

      do k=2,ns
        itmp = modelNumber(elnode(i, k))
        if(iModel.ne.plasticModel(itmp)) then
          write(iout, *) 'Multiple plastic models being used ! 
     &        Not supported at this time! '
          write(iptty, *) 'Multiple plastic models being used ! 
     &        Not supported at this time! '
          stop
        endif
      enddo

      D = 0.0d0
      if(iModel.eq.1) then

        e1bar = 0.0d0
        e2bar = 0.0d0
        e3bar = 0.0d0

        do k=1,ns
          node = elnode(i, k)
          E = elastic_mod(node)
          nu = poisson(node)
          lambda = E*nu/((1 + nu)*(1 - 2.0d0*nu))
          G = E/(2.0d0*(1 + nu))
          
          e1bar = e1bar + Psi(i, j, k)*(lambda + 2.0d0*G)
          e2bar = e2bar + Psi(i, j, k)*lambda
          e3bar = e3bar + Psi(i, j, k)*G
        enddo

        D(1,1) = e1bar
        D(2,2) = e1bar
        D(3,3) = e1bar
        
        D(1,2) = e2bar
        D(1,3) = e2bar
        D(2,1) = e2bar
        D(2,3) = e2bar
        D(3,1) = e2bar
        D(3,2) = e2bar

        D(4,4) = e3bar
        D(5,5) = e3bar
        D(6,6) = e3bar
      else
        write(iout,*) '***ERROR: elastic stiffness routine called for
     &   a node not set to be elastic'
      endif
      
      end subroutine fem_elastic_stiffness
