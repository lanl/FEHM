      subroutine vonMises_stiffness(i, d_i)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************
! 
! Returns the 6x6 material stiffness matrix for von Mises material when
! using the FEHM algorithm
! 
! Author : Sai Rapaka
! 
      
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
      
      end subroutine vonMises_stiffness
