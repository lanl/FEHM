      subroutine compute_permfactor_vonMises(el, node_k, duu, dvv, dww
     &, dpp, flag_u_pp)
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
! similar to compute_permfactor but use effective stress
!
! Author : Satish Karra
! Date updated: 1/26/12

      use comai, only: nei, neq, ns, iout
      use combi, only: nelm
      use comdi, only: t, tini, phi, phini
      use comsi, only: ipermstr2, perx_m, pery_m, perz_m
      use comsi, only: e1, e2, e3
      use comsi, only: ipermstr22
      use comsi, only: spm1f, spm2f, spm3f, spm4f, du, dv, dw, alp, bulk
      use comfem
      use comsi, only: modelNumber, plasticParam1, plasticModel
      use comai, only: iptty
      use comsi, only: plastic_strain, ispmd

      implicit none

      double precision                   :: duu, dvv, dww
      double precision                   :: dpp
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
      double precision                   :: p_eff,mean_str_eff
      double precision                   :: dev_stress(6)
      double precision                   :: G(8)
      double precision                   :: trace, sigma_y, phi_trial 
      double precision                   :: Q, J2, deltaGamma, G_bar
      double precision                   :: accum_pstrain
      double precision                   :: phi_tol, mean_accum_pstrain
      double precision                   :: accum_pstrain_max
      double precision                   :: norm_per, change_pstrain
      
      integer, parameter                 :: numEdges = 28
      integer el, node_k, i, j, k
      integer node_I, node_J, edge_1, edge_2
      logical recompute
      integer flag_u_pp
      integer itmp, iModel
      logical iUnload

      str_max = 1.0d0
      mean_str = 0.0
      mean_str_eff = 0.0
      p_eff = 0.0
      mean_accum_pstrain = 0.0
      
      if(ispmd.eq.100) then
        perx_m  = spm1f(1)
c        pery_m  = spm2f(1)
c        perz_m  = spm3f(1)
        accum_pstrain_max = spm4f(1)
      else
          write(iout, *) 'Only permmodel 100 is supported ! 
     &        with von Mises! '
          write(iptty, *) 'Only permmodel 100 is supported ! 
     &        with von Mises! '
          stop
      endif

      recompute = .false.
      do j=1,ns
         if(flag_u_pp.eq.1) then
            if(elnode(el,j).eq.node_k) then
               recompute=.true.
            endif
         endif
      enddo
      
      if(recompute) then
        ! Recompute the stress in the element
        do j=1,numgausspoints
          ! first compute the strain
          alphadeltaT = 0.0d0
          betadeltaP = 0.0d0
          e1bar = 0.0d0
          e2bar = 0.0d0
          e3bar = 0.0d0

          disp = 0.0d0
          do k=1,ns
            node(k) = elnode(el, k)
            if(node(k).eq.node_k) then
              disp(3*k-2) = duu
              disp(3*k-1) = dvv
              disp(3*k  ) = dww
            endif

            alphadeltaT = alphadeltaT + alp(node(k))*
     &            Psi(el, j, k)*(t(node(k)) - tini(node(k)))
            betadeltaP = betadeltaP + bulk(node(k))*
     &            Psi(el, j, k)*(phi(node(k)) - phini(node(k)))
            e1bar = e1bar + e1(node(k))*Psi(el, j, k)
            e2bar = e2bar + e2(node(k))*Psi(el, j, k)
            e3bar = e3bar + e3(node(k))*Psi(el, j, k)
            G(k) = 0.5d0*(e1(k) - e2(k))
            G_bar = G_bar + G(k)*Psi(el ,j, k)


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

          gp_stress = matmul(D, gp_strain)

          
          ! Assuming that when in plastic region, whatever is the change in 
          ! strain there is, it is due to plastic strain 
          change_pstrain = sqrt(sum(gp_strain(1:6)**2))
          accum_pstrain = plastic_strain(el, j) + change_pstrain
!             trace = (gp_stress(1) + gp_stress(2) + gp_stress(3))/3
!             dev_stress(1) = gp_stress(1) - trace
!             dev_stress(2) = gp_stress(2) - trace
!             dev_stress(3) = gp_stress(3) - trace
!             dev_stress(4) = gp_stress(4) 
!             dev_stress(5) = gp_stress(5)          
!             dev_stress(6) = gp_stress(6)             
! 
!           !! Going through von Mises calculations           
!             call J2Invariant(J2, dev_stress)
!             itmp = modelNumber(elnode(el, 1))
!             sigma_y = plasticParam1(itmp)
!             iModel = plasticModel(itmp)
! 
!             Q = sqrt(3.0d0*J2)
!             phi_trial = Q - sigma_y
!             phi_tol = 1.0d-6
!             iUnload = .false.
!             if((phi_trial.ge.phi_tol).and.(iUnload.eqv..false.)) then
!               deltaGamma = phi_trial/(3.0d0*G_bar)
!               accum_pstrain = plastic_strain(el,j) + deltaGamma
!             else
!               write(iout, *) 'In the elastic region ! 
!      &          Should not be checking for von Mises criterion! '
!               write(iptty, *) 'In the elastic region ! 
!      &          Should not be checking for von Mises criterion! '
!            endif
          
          mean_str = mean_str + gp_stress
          mean_accum_pstrain = mean_accum_pstrain + accum_pstrain
        enddo
        mean_str = mean_str/numgausspoints
        mean_accum_pstrain = mean_accum_pstrain/numgausspoints
      else
        mean_str = 0.0d0
        mean_accum_pstrain = 0.0d0
        do j=1,numgausspoints
          mean_str = mean_str + fem_stress(el, j, :)
          mean_accum_pstrain = mean_accum_pstrain + 
     &                            plastic_strain(el,j)
        enddo
        mean_str = mean_str/numgausspoints
        mean_accum_pstrain = mean_accum_pstrain/numgausspoints
      endif
     
      do j=1, numEdges
        node_I = elnode(el, edges(j,1))
        node_J = elnode(el, edges(j,2))
        edge_1 = edgeNum1(el, j)
        edge_2 = edgeNum2(el, j)
      
        !! Perm dependence on accumulated plastic strain
        !! ramp function upto a maximum value in acc. plastic strain
        !! maximum acc. plastic strain -- Karra
        if(mean_accum_pstrain.gt.0.0) then
          if(mean_accum_pstrain.lt.accum_pstrain_max) then
            fac = mean_accum_pstrain/accum_pstrain_max
            else 
            fac = 1.d0 
          endif
        else 
          fac = 0.d0 
        endif

c        norm_per = sqrt(perx_m*perx_m + pery_m*pery_m +
c     &                  perz_m*perz_m)
        norm_per = perx_m

        perm_fac = fac*(norm_per - 1.0d0) + 1.0d0

        write(iout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        write(iout,*) 'mean_accum_pstrain, fac,norm_per,perm_fac'
        write(iout,*) mean_accum_pstrain, fac,norm_per,perm_fac
        write(iout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
 
        permfactor(edge_1, :) = permfactor(edge_1, :) + perm_fac
        permfactor(edge_2, :) = permfactor(edge_2, :) + perm_fac
    
       enddo

      end subroutine compute_permfactor_vonMises
