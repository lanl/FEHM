      subroutine material_stress_update()
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
! Top level code to call the correct stress update routine when using the
! FEHM algorithm
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
      
      integer i
            
      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stress update routine called 
     &   without plastic flag being set! '
        stop
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
