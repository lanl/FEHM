      subroutine compute_ipsi()

      use comfem, only: iPsi

      implicit none

      real*8   :: xc(8), yc(8), zc(8)
      real*8   :: zeta, eta, mu
      real*8   :: xmult, solved
      integer i,j,k
      real*8   :: Nmat(8, 8), rhs(8, 8)
      real*8   :: sol(8)

      xc = 1.0d0
      yc = 1.0d0
      zc = 1.0d0
      xc=(/-1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)
      yc=(/-1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0/)
      zc=(/-1.0d0, -1.0d0, -1.0d0, -1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
      do i=1,8
        do j=1,8
          zeta = xc(i)/dsqrt(3.0d0)
          eta = yc(i)/dsqrt(3.0d0)
          mu = zc(i)/dsqrt(3.0d0)
          Nmat(i, j) = (1.0d0/8.0d0)*(1 + xc(j)*zeta)*
     &                 (1 + yc(j)*eta)*(1 + zc(j)*mu)
        enddo
      enddo

      rhs = 0.d00
      do i=1,8
        rhs(i,i) = 1.0d0
      enddo

      do k=1,7
        do i=k+1,8
          xmult = Nmat(i,k)/Nmat(k,k)
          do j=1,8
            Nmat(i,j) = Nmat(i,j) - xmult*Nmat(k,j)
            rhs(i,j) = rhs(i,j) - xmult*rhs(k,j)
          enddo
        enddo
      enddo

      do k=1,8
        ! Solve for the k-th column of rhs
        sol(8) = rhs(8, k)/Nmat(8,8)
        do j=7,1,-1
          solved = 0.0d0
          do i=j+1,8
            solved = solved + sol(i)*Nmat(j,i)
          enddo
          sol(j) = (rhs(j, k) - solved)/Nmat(j, j)
        enddo
        do j=1,8
          iPsi(j, k) = sol(j)
        enddo
      enddo

      return
      end subroutine compute_ipsi

