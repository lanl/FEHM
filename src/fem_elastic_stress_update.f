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
      use comai, only: iout
      use comsi, only: alp, bulk, iPlastic, plasticModel, modelNumber,
     &     stress_anisotropy_in
      use comai, only: iout, iptty, iad, ns
      use comdi, only: t, phi, tini, phini
      use comfem

      implicit none
      integer                      :: i, j, k, node
      real*8,  dimension(6)        :: gp_stress, gp_strain

      real*8,  dimension(6,6)      :: D
      real*8                       :: alpha, deltaT, beta, deltaP
      real*8                       :: alphadeltaT, betadeltaP

      fem_strain(i,j,:) = conv_strain(i,j,:) + gp_strain

      alphadeltaT = 0.0d0
      betadeltaP = 0.0d0

      do k=1,ns
        node = elnode(i,k)
        alpha = alp(node)
        deltaT = t(node) - tini(node)
        beta = bulk(node)
        deltaP = phi(node) - phini(node)
        alphadeltaT = alphadeltaT + 
     &                 Psi(i,j,k)*alpha*deltaT
        betadeltaP = betadeltaP + 
     &                 Psi(i,j,k)*beta*deltaP
      enddo

      gp_strain = fem_strain(i,j,:)
      gp_strain(1) = gp_strain(1) - alphadeltaT - betadeltaP
      gp_strain(2) = gp_strain(2) - alphadeltaT - betadeltaP
      gp_strain(3) = gp_strain(3) - alphadeltaT - betadeltaP

      if (stress_anisotropy_in) then
c s karra 17May2012
         call fem_transverse_isotropy_elastic_stiffness(i, j, D)
      else
         call fem_elastic_stiffness(i, j, D)
      endif

      gp_stress = matmul(D, gp_strain)

      fem_stress(i,j,:) = gp_stress

      end subroutine fem_elastic_stress_update

