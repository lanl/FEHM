      subroutine interpolate_q2n_3r(val, f, k)

      ! TODO : add a description of function - only interpolates from quadrature
      ! points to nodes

      use comfem

      implicit none
      real*8,  dimension(8)        :: f
      real*8                       :: val
      integer                      :: k
      
      real*8                       :: sq3, zeta, eta, mu
      real*8,  dimension(8)        :: N
      integer                      :: i

      sq3 = sqrt(3.0d0)
      
      zeta = gpcord(k,1)*sq3*sq3
      eta  = gpcord(k,2)*sq3*sq3
      mu   = gpcord(k,3)*sq3*sq3

      N(1) = 0.125d0*(1 - zeta)*(1 - eta)*(1 - mu)
      N(2) = 0.125d0*(1 + zeta)*(1 - eta)*(1 - mu)
      N(3) = 0.125d0*(1 + zeta)*(1 + eta)*(1 - mu)
      N(4) = 0.125d0*(1 - zeta)*(1 + eta)*(1 - mu)
      N(5) = 0.125d0*(1 - zeta)*(1 - eta)*(1 + mu)
      N(6) = 0.125d0*(1 + zeta)*(1 - eta)*(1 + mu)
      N(7) = 0.125d0*(1 + zeta)*(1 + eta)*(1 + mu)
      N(8) = 0.125d0*(1 - zeta)*(1 + eta)*(1 + mu)
      
      val = 0.0d0
      do i=1,8
        val = val + N(i)*f(i)
      enddo

      end subroutine interpolate_q2n_3r
 
