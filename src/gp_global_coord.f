      subroutine gp_global_coord(i, j, x, y, z)

      use comfem
      use comai, only: ns
      use combi, only: cord

      implicit none
      integer           :: i, j
      real*8            :: x, y, z

      integer           :: node, k

      x = 0.0d0
      y = 0.0d0
      z = 0.0d0

      do k=1,ns
        node = elnode(i, k)
        x = x + cord(node, 1)*Psi(i, j, k)
        y = y + cord(node, 2)*Psi(i, j, k)
        z = z + cord(node, 3)*Psi(i, j, k)
      enddo
        
      end subroutine gp_global_coord
