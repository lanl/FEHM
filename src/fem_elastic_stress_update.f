      subroutine fem_elastic_stress_update(i, j, gp_stress, gp_strain)

      use comsi, only: iPlastic, plasticModel, modelNumber
      use comai, only: iout, iptty, ns
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6)        :: gp_stress, gp_strain

      real*8,  dimension(6,6)      :: D

      call fem_elastic_stiffness(i, j, D)
      gp_stress = matmul(D, gp_strain)

      end subroutine fem_elastic_stress_update

