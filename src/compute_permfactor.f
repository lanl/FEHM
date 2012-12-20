      subroutine compute_permfactor(el, node_k, duu, dvv, dww
     &     , recompute_stress)
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
! compute elemental permbeality factors
!
! Author : Sai Rapaka

      use comai, only: nei, neq, ns, iout
      use combi, only: nelm
      use comdi, only: t, tini, phi, phini
      use comsi, only: ispm, ipermstr2, perx_m, pery_m, perz_m
      use comsi, only: strx_min, stry_min, strz_min, e1, e2, e3
      use comsi, only: spm1f, spm2f, spm3f
      use comsi, only: spm7f, spm8f, spm9f, du, dv, dw, alp, bulk
      use comfem

      implicit none

      double precision                   :: duu, dvv, dww
      double precision                   :: mean_str(6)
      double precision                   :: D(6,6)
      double precision                   :: B(6, 24)
      integer                            :: node(8)
      double precision                   :: disp(24)
      double precision                   :: gp_strain(6), gp_stress(6)
      double precision                   :: d_strain(6)
      double precision                   :: fac, perm_fac
      double precision                   :: str_max, str_tol
      double precision                   :: alphadeltaT, betadeltaP
      double precision                   :: e1bar, e2bar, e3bar
      integer, parameter                 :: numEdges = 28
      integer el, node_k, i, j, k
      integer node_I, node_J, edge_1, edge_2
      logical recompute_stress
      integer iispmd

      real*8 epsilon_perm

      epsilon_perm = 1.e-18


      str_max = 1.0d0
      mean_str = 0.0
c      permfactor = 1.0

      if(ipermstr2.ne.0) then
        strx_min = spm1f(1)
        stry_min = spm2f(1)
        strz_min = spm3f(1)
        perx_m  = spm7f(1)
        pery_m  = spm8f(1)
        perz_m  = spm9f(1)
      endif

c      recompute = .false.
c      do j=1,ns
c        if(elnode(el,j).eq.node_k) then
c          recompute=.true.
c        endif
c      enddo
      
      if(recompute_stress) then
        ! Recompute the stress in the element
        do j=1,numgausspoints
          ! first compute the strain
          alphadeltaT = 0.0d0
          betadeltaP = 0.0d0
          e1bar = 0.0d0
          e2bar = 0.0d0
          e3bar = 0.0d0

          !TODO: Should be setting up disp outside gausspoint loop
          disp = 0.0d0
          do k=1,ns
            node(k) = elnode(el, k)
            if(node(k).eq.node_k) then
              disp(3*k-2) = duu
              disp(3*k-1) = dvv
              disp(3*k  ) = dww
            endif
!            if(node(k).eq.node_k) then
!              disp(3*k-2) = disp(3*k-2) + duu
!              disp(3*k-1) = disp(3*k-1) + dvv
!              disp(3*k  ) = disp(3*k  ) + dww
!            endif
            alphadeltaT = alphadeltaT + alp(node(k))*
     &            Psi(el, j, k)*(t(node(k)) - tini(node(k)))
            betadeltaP = betadeltaP + bulk(node(k))*
     &            Psi(el, j, k)*(phi(node(k)) - phini(node(k)))
            e1bar = e1bar + e1(node(k))*Psi(el, j, k)
            e2bar = e2bar + e2(node(k))*Psi(el, j, k)
            e3bar = e3bar + e3(node(k))*Psi(el, j, k)
          enddo

          B = 0.0d0
          do k=1,ns
            B(1,3*(k-1) + 1) = dPsidX(el, j, k)
            B(2,3*(k-1) + 2) = dPsidY(el, j, k)
            B(3,3*(k-1) + 3) = dPsidZ(el, j, k)
            B(4,3*(k-1) + 1) = dPsidY(el, j, k)
            B(4,3*(k-1) + 2) = dPsidX(el, j, k)
            B(5,3*(k-1) + 2) = dPsidZ(el, j, k)
            B(5,3*(k-1) + 3) = dPsidY(el, j, k)
            B(6,3*(k-1) + 1) = dPsidZ(el, j, k)
            B(6,3*(k-1) + 3) = dPsidX(el, j, k)
          enddo
          
          !! Computing incremental change in strain
          d_strain = matmul(B, disp)

          !! Add incremental strain to elastic strain at gausspoint
          gp_strain = fem_strain(el, j, :) + d_strain 
          gp_strain(1) = gp_strain(1) - alphadeltaT - betadeltaP
          gp_strain(2) = gp_strain(2) - alphadeltaT - betadeltaP
          gp_strain(3) = gp_strain(3) - alphadeltaT - betadeltaP

          D = 0.0d0
          !! fem_material_stiffness returns elasto-plastic matrix for
          !! plastic, otherwise elastic matrix
          call fem_material_stiffness(el, j, D)
!          D(1,1) = e1bar
!          D(1,2) = e2bar
!          D(1,3) = e2bar
!          D(2,1) = e2bar
!          D(2,2) = e1bar
!          D(2,3) = e2bar
!          D(3,1) = e2bar
!          D(3,2) = e2bar
!          D(3,3) = e1bar
!          D(4,4) = e3bar
!          D(5,5) = e3bar
!          D(6,6) = e3bar

          gp_stress = matmul(D, gp_strain)
          mean_str = mean_str + gp_stress
        enddo
        mean_str = mean_str/numgausspoints
      else
        mean_str = 0.0d0
        do j=1,numgausspoints
          mean_str = mean_str + fem_stress(el, j, :)
        enddo
        mean_str = mean_str/numgausspoints
      endif

      !! Effective stress
!      write(iout, *) 'Eff'
!      write(iout, *) mean_str(1:3)
!      write(iout, *) betadeltaP
!      write(iout, *) 3.0*e2bar + 2.0*e3bar
!      write(iout, *) (3.0*e2bar + 2.0*e3bar)*betadeltaP

      do j=1, numEdges
        node_I = elnode(el, edges(j,1))
        node_J = elnode(el, edges(j,2))
        edge_1 = edgeNum1(el, j)
        edge_2 = edgeNum2(el, j)
        
        !! Explicitly considering model 2 here
        !! Need to add a loop to look at model numbers and 
        !! break it into different routines
        if(mean_str(1).gt.strx_min) then
          fac = ((mean_str(1)-strx_min)/str_max)
          perm_fac = fac*(perx_m - 1.0d0) + 1.0d0
          permfactor(edge_1, 2) = permfactor(edge_1, 2) + perm_fac
          permfactor(edge_1, 3) = permfactor(edge_1, 3) + perm_fac
          permfactor(edge_2, 2) = permfactor(edge_2, 2) + perm_fac
          permfactor(edge_2, 3) = permfactor(edge_2, 3) + perm_fac
        endif
        if(mean_str(2).gt.stry_min) then
          fac = ((mean_str(2)-stry_min)/str_max)
          perm_fac = fac*(pery_m - 1.0d0) + 1.0d0
          permfactor(edge_1, 3) = permfactor(edge_1, 3) + perm_fac
          permfactor(edge_1, 1) = permfactor(edge_1, 1) + perm_fac
          permfactor(edge_2, 3) = permfactor(edge_2, 3) + perm_fac
          permfactor(edge_2, 1) = permfactor(edge_2, 1) + perm_fac
        endif
        if(mean_str(3).gt.strz_min) then
          fac = ((mean_str(3)-strz_min)/str_max)
          perm_fac = fac*(perz_m - 1.0d0) + 1.0d0
          permfactor(edge_1, 1) = permfactor(edge_1, 1) + perm_fac
          permfactor(edge_1, 2) = permfactor(edge_1, 2) + perm_fac
          permfactor(edge_2, 1) = permfactor(edge_2, 1) + perm_fac
          permfactor(edge_2, 2) = permfactor(edge_2, 2) + perm_fac
        endif
      enddo

      end subroutine compute_permfactor
