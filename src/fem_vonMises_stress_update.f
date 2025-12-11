      subroutine fem_vonMises_stress_update(i, j, gp_stress, gp_strain,
     & iUnload )
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
! Computes the stresses from strains for the von Mises plasticity model
! 
! Author : Sai Rapaka
!

      use comsi, only: iPlastic, plasticModel, modelNumber, e1, e2, e3
      use comsi, only: isPlastic, plasticParam1, plastic_strain
      use comsi, only: du, dv, dw, alp, bulk
      use comdi, only: t, phi, tini, phini
      use comai, only: iout, iptty, iad, ns, nei
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6)        :: gp_stress, gp_strain, dstrain
      real*8,  dimension(6)        :: tot_strain, dev_strain
      integer, dimension(8)        :: node
      real*8,  dimension(6, 24)    :: B
      real*8,  dimension(24)       :: disp
      logical                      :: iUnload
      real*8,  dimension(8)        :: alpha, beta, deltaT, deltaP
      real*8,  dimension(8)        :: lambda, G
      real*8,  dimension(6,6)      :: D
      real*8                       :: alphadeltaT, betadeltaP
      real*8 vol_strain, pressure
      real*8 bulk_mod, lambda_bar, G_bar, J2, Q, deltaGamma
      real*8 phi_trial, phi_tol, fac, sigma_y
      integer itmp, iModel
      integer k

      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : von Mises stress update called
     &   without plastic flag being set '
      endif

      if(.not.allocated(conv_pstrain)) then
        allocate(conv_pstrain(nei, numgausspoints))
        conv_pstrain = 0.0d0
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

      phi_tol = 1.0d-6
      isPlastic(i,j) = 0

      bulk_mod = 0.0d0
      lambda_bar = 0.0d0
      G_bar = 0.0d0

      sigma_y = plasticParam1(itmp)

      alphadeltaT = 0.0d0
      betadeltaP  = 0.0d0

      do k=1, ns
        node(k) = elnode(i, k)
        lambda(k) = e2(node(k))
        G(k) = 0.5d0*(e1(k) - e2(k))
        lambda_bar = lambda_bar + lambda(k)*Psi(i,j,k)
        G_bar = G_bar + G(k)*Psi(i,j,k)
        alpha(k) = alp(node(k))
        deltaT(k) = t(node(k)) - tini(node(k))
        beta(k) = bulk(node(k))
        deltaP(k) = phi(node(k)) - phini(node(k))
        alphadeltaT = alphadeltaT + 
     &                 Psi(i,j,k)*alpha(k)*deltaT(k)
        betadeltaP = betadeltaP + 
     &                 Psi(i,j,k)*beta(k)*deltaP(k)

      enddo

!      if(iad.eq.0) then
!        fem_strain(i,j,:) = 0.0d0
!      endif

      ! Change gp_strain from incremental strain to total new strain
      tot_strain = conv_strain(i, j, :) + gp_strain

      ! tot_strain(1) = tot_strain(1) - alphadeltaT - betadeltaP
      ! tot_strain(2) = tot_strain(2) - alphadeltaT - betadeltaP
      ! tot_strain(3) = tot_strain(3) - alphadeltaT - betadeltaP

      bulk_mod = lambda_bar + (2.0d0/3.0d0)*G_bar
      vol_strain = tot_strain(1) + tot_strain(2) + tot_strain(3)
      ! vol_strain = vol_strain - 3.0d0*(alphadeltaT + betadeltaP)
      pressure = bulk_mod*(vol_strain -3.d0*(alphadeltaT + betadeltaP))

      dev_strain = tot_strain
      dev_strain(1) = dev_strain(1) - vol_strain/3.0d0
      dev_strain(2) = dev_strain(2) - vol_strain/3.0d0 
      dev_strain(3) = dev_strain(3) - vol_strain/3.0d0
      dev_strain(4) = dev_strain(4)/2.0d0
      dev_strain(5) = dev_strain(5)/2.0d0
      dev_strain(6) = dev_strain(6)/2.0d0

      call J2Invariant(J2, 2.0d0*G_bar*dev_strain)
      Q = sqrt(3.0d0*J2)
      phi_trial = Q - sigma_y

      if((phi_trial.ge.phi_tol).and.(iUnload.eqv..false.)) then
        isPlastic(i, j) = 1
        deltaGamma = phi_trial/(3.0d0*G_bar)
        plastic_strain(i,j) = conv_pstrain(i,j) + deltaGamma

        !  phi_trial = Q - 3.0d0*G_bar*deltaGamma - sigma_y/sqrt(3.0d0)

        fac = 2.0d0*G_bar*(1.0d0 - 3.0d0*G_bar*deltaGamma/Q)
        gp_stress(1) = fac*dev_strain(1) + pressure
        gp_stress(2) = fac*dev_strain(2) + pressure
        gp_stress(3) = fac*dev_strain(3) + pressure
        gp_stress(4) = fac*dev_strain(4)
        gp_stress(5) = fac*dev_strain(5)
        gp_stress(6) = fac*dev_strain(6)

        fac = fac/(2.0d0*G_bar)
        tot_strain(1) = fac*dev_strain(1) + vol_strain/3.0d0
        tot_strain(2) = fac*dev_strain(2) + vol_strain/3.0d0
        tot_strain(3) = fac*dev_strain(3) + vol_strain/3.0d0
        tot_strain(4) = 2.0d0*fac*dev_strain(4)
        tot_strain(5) = 2.0d0*fac*dev_strain(5)
        tot_strain(6) = 2.0d0*fac*dev_strain(6)

      else
        fac = 2.0d0*G_bar
        gp_stress(1) = fac*dev_strain(1) + pressure
        gp_stress(2) = fac*dev_strain(2) + pressure
        gp_stress(3) = fac*dev_strain(3) + pressure
        gp_stress(4) = fac*dev_strain(4)
        gp_stress(5) = fac*dev_strain(5)
        gp_stress(6) = fac*dev_strain(6)
      endif

      fem_strain(i,j,:) = tot_strain
      fem_stress(i,j,:) = gp_stress

      end subroutine fem_vonMises_stress_update

