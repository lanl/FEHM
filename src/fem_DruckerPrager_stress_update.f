      subroutine fem_DruckerPrager_stress_update(i, j, gp_stress, 
     & gp_strain, iUnload)
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
! Computes the stresses from strains for the drucker-prager (without cap)
! plasticity model
! 
! Author : Satish Karra
! Date : March 26 2012
!

      use comsi, only: iPlastic, plasticModel, modelNumber, e1, e2, e3
      use comsi, only: isPlastic, plasticParam1, plastic_strain
      use comsi, only: plasticParam2, plasticParam3
      use comsi, only: du, dv, dw, alp, bulk
      use comdi, only: t, phi, tini, phini
      use comai, only: iout, iptty, iad, ns, nei
      use comfem

      implicit none
      integer                                :: i, j, k   
      double precision,  dimension(6)        :: gp_stress, gp_strain
      double precision,  dimension(6)        :: dstrain
      double precision,  dimension(6)        :: tot_strain, dev_strain
      integer, dimension(8)                  :: node
      logical                                :: iUnload
      double precision,  dimension(8)        :: alpha, beta, deltaT
      double precision,  dimension(8)        :: lambda, G, deltaP
      double precision,  dimension(6,6)      :: D
      double precision alphadeltaT, betadeltaP, vol_strain, pressure
      double precision bulk_mod, lambda_bar, G_bar, J2, deltaGamma
      double precision phi_trial, phi_tol, fac, eta, xi, rad
      integer itmp, iModel


      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : Drucker-Prager stress update called
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

      alphadeltaT = 0.0d0
      betadeltaP  = 0.0d0

!   Assigning Drucker Prager parameters for now
      eta = plasticParam1(itmp) 
      xi = plasticParam2(itmp)
      rad = plasticParam3(itmp)

      do k = 1, ns
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

      ! Change gp_strain from incremental strain to total new strain
      tot_strain = conv_strain(i, j, :) + gp_strain

      bulk_mod = lambda_bar + (2.0d0/3.0d0)*G_bar
      vol_strain = tot_strain(1) + tot_strain(2) + tot_strain(3)
      vol_strain = vol_strain - 3.0d0*(alphadeltaT + betadeltaP)
      pressure = bulk_mod*(vol_strain)

      dev_strain = tot_strain
      dev_strain(1) = dev_strain(1) - vol_strain/3.0d0
      dev_strain(2) = dev_strain(2) - vol_strain/3.0d0 
      dev_strain(3) = dev_strain(3) - vol_strain/3.0d0
      dev_strain(4) = dev_strain(4)/2.0d0
      dev_strain(5) = dev_strain(5)/2.0d0
      dev_strain(6) = dev_strain(6)/2.0d0

      call J2Invariant(J2, 2.0d0*G_bar*dev_strain)
      phi_trial = sqrt(J2) + eta*pressure - xi*rad 

      if((phi_trial.ge.phi_tol).and.(iUnload.eqv..false.)) then
        isPlastic(i, j) = 1

        deltaGamma = phi_trial/(G_bar + bulk_mod*eta**2)

        if (sqrt(J2) - G_bar*deltaGamma.ge.phi_tol) then
! going to the smooth part of the cone
          fac = (1.d0 - G_bar*deltaGamma/sqrt(J2))*2.d0*G_bar
          pressure = pressure - bulk_mod*deltaGamma*eta
          gp_stress(1) = fac*dev_strain(1) + pressure
          gp_stress(2) = fac*dev_strain(2) + pressure
          gp_stress(3) = fac*dev_strain(3) + pressure
          gp_stress(4) = fac*dev_strain(4)
          gp_stress(5) = fac*dev_strain(5)
          gp_stress(6) = fac*dev_strain(6)
          plastic_strain(i,j) = conv_pstrain(i,j) + eta*deltaGamma
          fac = fac/(2.0d0*G_bar)
          tot_strain(1) = fac*dev_strain(1) + pressure/(3.0d0*bulk_mod)
          tot_strain(2) = fac*dev_strain(2) + pressure/(3.0d0*bulk_mod)
          tot_strain(3) = fac*dev_strain(3) + pressure/(3.0d0*bulk_mod)
          tot_strain(4) = 2.d0*fac*dev_strain(4)
          tot_strain(5) = 2.d0*fac*dev_strain(5) 
          tot_strain(6) = 2.d0*fac*dev_strain(6) 
        else
! going to return to the apex
          deltaGamma = (pressure - xi/eta*rad)/(eta*bulk_mod)
          gp_stress(1) = pressure - bulk_mod*eta*deltaGamma
          gp_stress(2) = pressure - bulk_mod*eta*deltaGamma
          gp_stress(3) = pressure - bulk_mod*eta*deltaGamma
          gp_stress(4) = 0.d0
          gp_stress(5) = 0.d0
          gp_stress(6) = 0.d0
          plastic_strain(i,j) = conv_pstrain(i,j) +
     &                          deltaGamma*sqrt(0.5d0 + eta**2/3.d0)
          tot_strain(1) = gp_stress(1)/(3.0d0*bulk_mod)
          tot_strain(2) = gp_stress(2)/(3.0d0*bulk_mod)
          tot_strain(3) = gp_stress(3)/(3.0d0*bulk_mod)
          tot_strain(4) = 0.d0
          tot_strain(5) = 0.d0
          tot_strain(6) = 0.d0
        endif

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

      end subroutine fem_DruckerPrager_stress_update

