      subroutine fem_material_stress_update(i, j, gp_stress, gp_strain,
     &     gp_strain_mech,DSai)
      
      use comsi, only: iPlastic, plasticModel, modelNumber
      use comai, only: iout, iptty, ns
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6)        :: gp_stress, gp_strain
      real*8,  dimension(6)        :: gp_strain_mech
      real*8 xgp,ygp,zgp
      real*8,  dimension(6, 6)     :: DSai

      integer                      :: itmp, iModel, k,i1

      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stress update routine called 
     &   without plastic flag being set! '
        write(iptty,*) '***ERROR : material stress update routine called
     &   without plastic flag being set! '
      endif
      
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

      if(iModel.eq.1) then
        ! Linear, isotropic, elastic rock
        call fem_elastic_stress_update(i, j, gp_stress, gp_strain)
      else if(iModel.eq.2) then
        ! von Mises material
        ! material properties such as yield_stress etc are available
        ! as global variables
        call fem_vonMises_stress_update(i, j, gp_stress, gp_strain)
      endif
      
      end subroutine fem_material_stress_update
