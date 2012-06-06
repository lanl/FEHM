      subroutine fem_material_stress_update(i, j, gp_stress, gp_strain,
     & iUnload)

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
! Top level code to call the correct stress update routine when the 
! 'plastic' submacro is used with 'fem' computations
! 
! Author : Sai Rapaka
!
      
      use comsi, only: iPlastic, plasticModel, modelNumber
      use comai, only: iout, iptty, ns
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6)        :: gp_stress, gp_strain
      real*8 xgp,ygp,zgp
      logical                      :: iUnload
      integer                      :: itmp, iModel, k,i1

      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : material stress update routine called 
     &   without plastic flag being set! '
        write(iptty,*) '***ERROR : material stress update routine called
     &   without plastic flag being set! '
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
        iUnload = .false.
        call fem_elastic_stress_update(i, j, gp_stress, gp_strain)
      else if(iModel.eq.2) then
        ! von Mises material
        ! material properties such as yield_stress etc are available
        ! as global variables
        call fem_vonMises_stress_update(i, j, gp_stress, gp_strain,
     & iUnload )
      else if(iModel.eq.3) then
         ! s karra, 17 May 2012, Drucker Prager model (without cap)
        call fem_DruckerPrager_stress_update(i, j, gp_stress,
     &        gp_strain, iUnload)
      endif
      
      end subroutine fem_material_stress_update
