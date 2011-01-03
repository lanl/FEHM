      subroutine fem_material_stiffness(i, j, D)
      
      use comai, only: iout, iptty, ns
      use comsi, only: iPlastic, modelNumber, plasticModel, e1, e2, e3
      use comfem

      implicit none
      real*8,  dimension(6,6)      :: D
      integer                      :: i, j, itmp
      integer                      :: iModel, k
      real *8                      :: e1bar, e2bar, e3bar

      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stiffness routine called 
     &   without plastic flag being set! '
        write(iptty,*) '***ERROR : material stiffness routine called 
     &   without plastic flag being set! '
        stop
      endif
      
      if(ifem.eq.0) then
        write(iout,*) '***ERROR : FEM material stiffness called without
     &   setting fem flag! '
        write(iptty,*) '***ERROR : FEM material stiffness called without
     &   setting fem flag! '
        stop
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
        call fem_elastic_stiffness(i, j, D)
      else if(iModel.eq.2) then
        ! von Mises material
        ! material properties such as yield_stress etc are available
        ! as global variables
        call fem_vonMises_stiffness(i, j, D)
      endif
      
      end subroutine fem_material_stiffness
