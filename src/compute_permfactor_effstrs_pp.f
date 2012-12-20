      subroutine compute_permfactor_effstrs_pp(el,node_J, duu, dvv, dww
     &, dpp, flag_u_pp, recompute_stress)
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
! Author : S Kelkar

      use comai, only: nei, neq, ns, iout
      use combi, only: nelm
      use comdi, only: t, tini, phi, phini
      use comsi, only: ipermstr2, perx_m, pery_m, perz_m
      use comsi, only: strx_min, stry_min, strz_min, e1, e2, e3
      use comsi, only: spm1f, spm2f, spm3f
      use comsi, only: spm7f, spm8f, spm9f, du, dv, dw, alp, bulk
      use comfem

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
      double precision                   :: efac_1,efac_2

      integer, parameter                 :: numEdges = 28
      integer el, node_J, i, j, k, j_edge
      integer node_1, node_2, edge_1, edge_2
      logical recompute_stress
      integer flag_u_pp


      str_max = 1.0d0
      mean_str = 0.0
      mean_str_eff = 0.0
      p_eff=0.0


      if(ipermstr2.ne.0) then
        strx_min = spm1f(1)
        stry_min = spm2f(1)
        strz_min = spm3f(1)
        perx_m  = spm7f(1)
        pery_m  = spm8f(1)
        perz_m  = spm9f(1)

      endif
c
c      recompute = .false.
c      do j=1,ns
c         if(flag_u_pp.eq.1) then
c            if(elnode(el,j).eq.node_J) then
c               recompute=.true.
c            endif
c         endif
c      enddo
c      
      if(recompute_stress) then
        ! Recompute the stress in the element
         disp = 0.0d0
         do k=1,ns
            node(k) = elnode(el, k)
            if(node(k).eq.node_J) then
               disp(3*k-2) = duu
               disp(3*k-1) = dvv
               disp(3*k  ) = dww
            endif
         enddo

        do j=1,numgausspoints
          ! first compute the strain
          alphadeltaT = 0.0d0
          betadeltaP = 0.0d0
          e1bar = 0.0d0
          e2bar = 0.0d0
          e3bar = 0.0d0

          do k=1,ns
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

      do j_edge=1, numEdges
        node_1 = elnode(el, edges(j_edge,1))
        node_2 = elnode(el, edges(j_edge,2))
        edge_1 = edgeNum1(el, j_edge)
        edge_2 = edgeNum2(el, j_edge)
        efac_1 = 3.d0*e2(node_1)+2.d0*e3(node_1)
        efac_2 = 3.d0*e2(node_2)+2.d0*e3(node_2)
        if(node_1.eq.node_J) then
       p_eff=0.5*(efac_1*bulk(node_1)*(phi(node_1)-phini(node_1)+dpp)
     &          + efac_2*bulk(node_2)*(phi(node_2)-phini(node_2)))     
        elseif(node_2.eq.node_J) then
       p_eff=0.5*(efac_1*bulk(node_1)*(phi(node_1)-phini(node_1))
     &          + efac_2*bulk(node_2)*(phi(node_2)-phini(node_2)+dpp))    
        else
       p_eff=0.5*(efac_1*bulk(node_1)*(phi(node_1)-phini(node_1))
     &          + efac_2*bulk(node_2)*(phi(node_2)-phini(node_2)))     
        endif
        !! Explicitly considering model 2 here, but with effective stress
        !! Need to add a loop to look at model numbers and 
        !! break it into different routines
c NOTE: sign convention here is -ve in compression, opposite to general
c fehm convention
c s kelkar 5 june 2012. If strx_min <or= 0 model is interpreted to be that
c for local failure driven by pore pressure only. 
c  If strx_min>0 , it is taken to be the usual model 2
        if(strx_min.gt.0.) then
           mean_str_eff = mean_str(1) + p_eff
           if(mean_str_eff.gt.strx_min) then
              fac = ((mean_str_eff-strx_min)/str_max)
              perm_fac = fac*(perx_m - 1.0d0) + 1.0d0
              permfactor(edge_1, 2) = permfactor(edge_1, 2) + perm_fac
              permfactor(edge_1, 3) = permfactor(edge_1, 3) + perm_fac
              permfactor(edge_2, 2) = permfactor(edge_2, 2) + perm_fac
              permfactor(edge_2, 3) = permfactor(edge_2, 3) + perm_fac
           endif
           mean_str_eff = mean_str(2) + p_eff
           if(mean_str_eff.gt.stry_min) then
              fac = ((mean_str_eff-stry_min)/str_max)
              perm_fac = fac*(pery_m - 1.0d0) + 1.0d0
              permfactor(edge_1, 3) = permfactor(edge_1, 3) + perm_fac
              permfactor(edge_1, 1) = permfactor(edge_1, 1) + perm_fac
              permfactor(edge_2, 3) = permfactor(edge_2, 3) + perm_fac
              permfactor(edge_2, 1) = permfactor(edge_2, 1) + perm_fac
           endif
           mean_str_eff = mean_str(3) + p_eff
           if(mean_str_eff.gt.strz_min) then
              fac = ((mean_str_eff-strz_min)/str_max)
              perm_fac = fac*(perz_m - 1.0d0) + 1.0d0
              permfactor(edge_1, 1) = permfactor(edge_1, 1) + perm_fac
              permfactor(edge_1, 2) = permfactor(edge_1, 2) + perm_fac
              permfactor(edge_2, 1) = permfactor(edge_2, 1) + perm_fac
              permfactor(edge_2, 2) = permfactor(edge_2, 2) + perm_fac
           endif
        else
           fac = p_eff + strx_min
           if(fac.gt.0.0) then
              fac=fac/str_max
              perm_fac = fac*(perx_m - 1.0d0) + 1.0d0
              permfactor(edge_1, 1) = permfactor(edge_1, 1) + perm_fac
              permfactor(edge_1, 2) = permfactor(edge_1, 2) + perm_fac
              permfactor(edge_1, 3) = permfactor(edge_1, 3) + perm_fac
              permfactor(edge_2, 1) = permfactor(edge_2, 1) + perm_fac
              permfactor(edge_2, 2) = permfactor(edge_2, 2) + perm_fac
              permfactor(edge_2, 3) = permfactor(edge_2, 3) + perm_fac
           endif
        endif
      enddo

      end subroutine compute_permfactor_effstrs_pp
