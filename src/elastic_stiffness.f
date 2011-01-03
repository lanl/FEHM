      subroutine elastic_stiffness(i, d_i)
      
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
      integer i, itmp, iModel
      real*8 E0, nu0, lambda0, G0
      
      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : elastic stiffness routine called 
     &   without plastic flag being set '
      endif
      
      itmp = modelNumber(i)
      iModel = plasticModel(itmp)
      
      d_i = 0.0d0
      if(iModel.eq.1) then
        E0 = elastic_mod(i)
        nu0 = poisson(i)
        lambda0 = E0*nu0/((1 + nu0)*(1 - 2.0d0*nu0))
        G0 = E0/(2.0d0*(1 + nu0))
        
        d_i(1:3, 1:3) = lambda0
        d_i(1,1) = d_i(1,1) + 2.0d0*G0
        d_i(2,2) = d_i(2,2) + 2.0d0*G0
        d_i(3,3) = d_i(3,3) + 2.0d0*G0
        d_i(4,4) = G0
        d_i(5,5) = G0
        d_i(6,6) = G0
      else
        write(iout,*) '***ERROR: elastic stiffness routine called for
     &   a node not set to be elastic'
      endif
      
      end subroutine elastic_stiffness