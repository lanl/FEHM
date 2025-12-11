      subroutine fem_material_stiffness_2D(i, j, D)
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
! Top level code to call the correct material stiffness routine when the 
! 'plastic' submacro is used with 'fem' computations
!
! Author : Sai Rapaka
!
      
      use comai, only: iout, iptty, ns
      use comsi, only: iPlastic, modelNumber, plasticModel, e1, e2, e3
      use comfem

      implicit none
      real*8,  dimension(3,3)      :: D
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
      
      end subroutine fem_material_stiffness_2D
