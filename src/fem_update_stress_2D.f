      subroutine fem_update_stress_2D(i, j)
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
! Computes stresses from strains in two-dimensional 'fem' computations
! 
! Author : Sai Rapaka
!

      use comai, only: ns,iout
      use comdi, only: t, phi, tini, phini
      use comsi, only: e1, e2, e3, du, dv, dw, alp, bulk, iPlastic
      use comfem

      implicit none

      integer                      :: i,j
      integer                      :: k
      integer, dimension(4)        :: node
      real*8,  dimension(4)        :: alpha, beta, deltaT, deltaP
      real*8,  dimension(3, 8)     :: B
      real*8,  dimension(3, 3)     :: D
      real*8,  dimension(3)        :: gp_stress, gp_strain
      real*8,  dimension(6)        :: gp_strain_mech
      real*8,  dimension(8)        :: disp

      real*8                       :: e1bar, e2bar, e3bar
      real*8                       :: alphadeltaT, betadeltaP
      logical                      :: iUnload


      ! first compute the strain
      do k=1,ns
        node(k) = elnode(i,k)
        disp(2*k-1) = du(node(k))
        disp(2*k) = dv(node(k))
        alpha(k) = alp(node(k))
        deltaT(k) = t(node(k)) - tini(node(k))
        beta(k) = bulk(node(k))
        deltaP(k) = phi(node(k)) - phini(node(k))
      enddo

      B = 0.0d0
      do k=1,ns
        B(1,2*(k-1) + 1) = dPsidX(i, j, k)
        B(2,2*(k-1) + 2) = dPsidY(i, j, k)
        B(3,2*(k-1) + 1) = dPsidY(i, j, k)
        B(3,2*(k-1) + 2) = dPsidX(i, j, k)
      enddo
          
      gp_strain = matmul(B, disp)
      fem_strain(i,j,:) = gp_strain

      e1bar = 0.0d0
      e2bar = 0.0d0
      e3bar = 0.0d0

      alphadeltaT = 0.0d0
      betadeltaP = 0.0d0

      do k=1,ns
        e1bar = e1bar + Psi(i, j, k)*e1(node(k))
        e2bar = e2bar + Psi(i, j, k)*e2(node(k))
        e3bar = e3bar + Psi(i, j, k)*e3(node(k))
        
        alphadeltaT = alphadeltaT + 
     &                 Psi(i,j,k)*alpha(k)*deltaT(k)
        betadeltaP = betadeltaP + 
     &                 Psi(i,j,k)*beta(k)*deltaP(k)

      enddo

      gp_strain(1) = gp_strain(1) - alphadeltaT - betadeltaP
      gp_strain(2) = gp_strain(2) - alphadeltaT - betadeltaP

      if(iPlastic.eq.1) then
        call fem_material_stress_update(i, j, gp_stress, gp_strain, 
     &        iUnload)
      else
        D = 0.0d0
        D(1,1) = e1bar
        D(1,2) = e2bar
        D(2,1) = e2bar
        D(2,2) = e1bar
        D(3,3) = e3bar
 
        gp_stress = matmul(D, gp_strain)
      endif

      fem_stress(i,j,:) = gp_stress

      end subroutine fem_update_stress_2D
