      subroutine material_stress_update()
      
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
      
      integer i
            
      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stress update routine called 
     &   without plastic flag being set! '
      endif
      
      call compute_strains()
      
      do i=1,n0
      
        itmp = modelNumber(i)
        iModel = plasticModel(itmp)
      
        if(iModel.eq.1) then
          ! Linear, isotropic, elastic rock
          call elastic_stress_update(i)
        else if(iModel.eq.2) then
          ! von Mises material
          ! material properties such as yield_stress etc are available
          ! as global variables
          call vonMises_stress_update(i)
        endif
      
      enddo  
      
      end subroutine material_stress_update