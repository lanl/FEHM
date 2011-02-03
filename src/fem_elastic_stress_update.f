      subroutine fem_elastic_stress_update(i, j, gp_stress, gp_strain)
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
! Compute the stresses from strains for a linear, isotropic, elastic solid
! when the 'plastic' submacro is used with 'fem' computations
! 
! Author : Sai Rapaka
!

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

