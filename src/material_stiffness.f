      subroutine material_stiffness(i, d_i)
      
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsi
      
      real*8, dimension(1:6, 1:6)  :: d_i
      integer i, itmp
      
      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stiffness routine called 
     &   without plastic flag being set! '
      endif
      
      itmp = modelNumber(i)
      iModel = plasticModel(itmp)
      
      if(iModel.eq.1) then
        ! Linear, isotropic, elastic rock
        call elastic_stiffness(i, d_i)
      else if(iModel.eq.2) then
        ! von Mises material
        ! material properties such as yield_stress etc are available
        ! as global variables
        call vonMises_stiffness(i, d_i)
      endif
      
      end subroutine material_stiffness