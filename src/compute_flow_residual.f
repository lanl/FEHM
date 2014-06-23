      subroutine compute_flow_residual(Ri,node_I,node_k,duu,dvv,dww
     &     , dpp, flag_u_pp)
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
!compute residuals to allow numerical differentiation
! 
!
! Author : Sai Rapaka

      use comai, only: nei, neq, upwgt, dnwgt, sx_mult
      use comai, only: igrav, grav, iout
      use combi, only: nelm, sx, sx1, cord, istrw, istrw_itfc
      use combi, only: red_factor, isox, isoy, isoz
      use comci, only: deni, div, dil, rovf, rolf, enlf, envf, denei
      use comdi, only: sk, phi, pcp, pnx, pny, pnz, thx, thy, thz, t, qh
      use comsi, only: perx_m, pery_m, perz_m, flag_permmodel
      use comfem

      implicit none
      double precision                   :: duu, dvv, dww
      double precision                   :: dpp
      double precision                   :: Ri(2)
      real*8                             :: harmonic_mean
      real*8                             :: sx2c
      real*8                             :: trans_liq, trans_vap
      real*8                             :: rho_liq_bar, rho_vap_bar
      real*8                             :: kx_bar, ky_bar, kz_bar
      real*8                             :: tx_bar, ty_bar, tz_bar
      real*8                             :: k_bar, t_bar
      real*8                             :: reduction_factor
      real*8                             :: P_Liq_I, P_Liq_J
      real*8                             :: P_Vap_I, P_Vap_J
      real*8                             :: h_l, h_v
      real*8                             :: delx2, dely2, delz2, dis2
      real*8                             :: delZ
      real*8                             :: fid, gid
      real*8                             :: grav_air
      real*8                             :: delP_liq, delP_vap
      real*8                             :: flow_liq, flow_vap
      real*8                             :: dis_tol
      real*8                             :: pmf(3)
      integer i, j, iw
      integer node_I, node_J, node_K
      integer j_begin, j_end, el
      integer flag_u_pp      
      logical recompute_stress
      parameter(dis_tol = 1.0d-12)

      Ri = 0.0d0
      j_begin = nelm(node_I)+1
      j_end   = nelm(node_I+1)

      if(node_k.gt.0. and. flag_permmodel.eq.1) then
        ! save permfactor
        permtmp = permfactor
        ! overwrite with updated values
        permfactor = 0.0d0
c loop over elements connected to the node
        do j= NodeElems(node_I)+1,NodeElems(node_I+1)
           recompute_stress = .True.
           el = NodeElems(j)
c           call compute_permfactor(el, node_k, duu, dvv, dww
c     &          , recompute_stress)
c           call compute_permfactor_effstrs(el, node_k, duu, dvv, dww
c     &          ,dpp, flag_u_pp, recompute_stress)
           call compute_permfactor_effstrs_pp(el,node_k,duu, dvv, dww
     &          ,dpp, flag_u_pp, recompute_stress)
c     if(ipermstr2.ne.0) then 
c     &         call fem_permfactor_2(el, node_k, duu, dvv, dww
c     &         ,dpp, flag_u_pp, recompute_stress)
c     if(ipermstr22.ne.0) then
c     &            call fem_permfactor_2(el, node_k, duu, dvv, dww
c     &         ,dpp, flag_u_pp, recompute_stress)
c     if(ipermstr41.ne.0) then
c     &            call fem_permfactor_2(el, node_k, duu, dvv, dww
c     &         ,dpp, flag_u_pp, recompute_stress)
        enddo 

        do j=j_begin, j_end
          if(numelems(j).gt.0) then
            permfactor(j,:) = permfactor(j,:)/numelems(j)
          endif
          permfactor(j,1) = min(permfactor(j,1), perx_m)
          permfactor(j,2) = min(permfactor(j,2), pery_m)
          permfactor(j,3) = min(permfactor(j,3), perz_m)
          permfactor(j,1) = max(permfactor(j,1), 1.0d0)
          permfactor(j,2) = max(permfactor(j,2), 1.0d0)
          permfactor(j,3) = max(permfactor(j,3), 1.0d0)
        enddo
      endif

      ! start the loop over the list of neighbors
      do j=j_begin, j_end

        node_J = nelm(j)

        if(node_J.eq.node_I) goto 10

        ! effective permeabilities        
        kx_bar = harmonic_mean(pnx(node_I), pnx(node_J))
        ky_bar = harmonic_mean(pny(node_I), pny(node_J))
        kz_bar = harmonic_mean(pnz(node_I), pnz(node_J))

        ! effective thermal conductivities
        tx_bar = harmonic_mean(thx(node_I), thx(node_J))
        ty_bar = harmonic_mean(thy(node_I), thy(node_J))
        tz_bar = harmonic_mean(thz(node_I), thz(node_J))

        ! multiply by permeability factors
        pmf = 1.0
        if(allocated(permfactor)) then
           pmf(1) = permfactor(j, 1)
           pmf(2) = permfactor(j, 2)
           pmf(3) = permfactor(j, 3)
        endif
        kx_bar = kx_bar*pmf( 1)
        ky_bar = ky_bar*pmf( 2)
        kz_bar = kz_bar*pmf( 3)

        reduction_factor = red_factor(istrw_itfc(j - neq - 1))

        P_Vap_I = phi(node_I)
        P_Vap_J = phi(node_J)
        P_Liq_I = P_Vap_I - pcp(node_I)
        P_Liq_J = P_Vap_J - pcp(node_J)

        iw = istrw(j - neq - 1)
        sx2c = sx(iw, isox) + sx(iw, isoy) + sx(iw, isoz)

        delx2=(cord(node_J,1)-cord(node_I,1))**2
        dely2=(cord(node_J,2)-cord(node_I,2))**2
        delz2=(cord(node_J,3)-cord(node_I,3))**2
        dis2=delx2+dely2+delz2

        if(dis2.gt.dis_tol) then
          k_bar=sx2c*dis2/(delx2/kx_bar+dely2/ky_bar+delz2/kz_bar)
          t_bar=sx2c*dis2/(delx2/tx_bar+dely2/ty_bar+delz2/tz_bar)
        else
          k_bar=sx2c*sx_mult*max(kx_bar,ky_bar,kz_bar)
          t_bar=sx2c*sx_mult*max(tx_bar,ty_bar,tz_bar)
        endif

        k_bar = k_bar*reduction_factor

        delP_vap = P_Vap_J - P_Vap_I
        delP_liq = P_Liq_J - P_Liq_I

        rho_liq_bar = 0.5d0*(rolf(node_I) + rolf(node_J))
        rho_vap_bar = 0.5d0*(rovf(node_I) + rovf(node_J))

        delZ = cord(node_J, igrav) - cord(node_I, igrav)
        flow_liq = k_bar*(delP_liq - rho_liq_bar*grav*delZ)
        flow_vap = k_bar*(delP_vap - rho_vap_bar*grav_air*delZ)

        ! Set upweighting parameter
        if(flow_liq.lt.0.0d0) then
          fid = dnwgt
        elseif(flow_liq.gt.0.0d0) then
          fid = upwgt
        else
          fid = 0.5d0
        endif

        if(flow_vap.lt.0.0d0) then
          gid = dnwgt
        elseif(flow_vap.gt.0.0d0) then
          gid = upwgt
        else
          gid = 0.5d0
        endif

        ! rel_perm*rho/mu
        trans_liq = fid*dil(node_J) + (1.0d0 - fid)*dil(node_I)
        trans_vap = gid*div(node_J) + (1.0d0 - gid)*div(node_I)
        ! rho*u_l*h_l
        h_l = fid*dil(node_J)*enlf(node_J) + 
     &        (1.0d0 - fid)*dil(node_I)*enlf(node_I)
        h_v = gid*div(node_J)*envf(node_J) +
     &        (1.0d0 - gid)*div(node_I)*envf(node_I)

        Ri(1) = Ri(1) + flow_liq*trans_liq
        Ri(1) = Ri(1) + flow_vap*trans_vap

        Ri(2) = Ri(2) + flow_liq*h_l
        Ri(2) = Ri(2) + flow_vap*h_v

        ! Add heat conduction terms
        Ri(2) = Ri(2) + t_bar*(t(node_J) - t(node_I))

10      continue
      enddo
        
      Ri(1) = Ri(1) + sx1(node_I)*deni(node_I) + sk(node_I)
      Ri(2) = Ri(2) + sx1(node_I)*denei(node_I) + qh(node_I)

      if(node_k.gt.0. and. flag_permmodel.eq.1) then
        ! recover permfactors
        permfactor = permtmp
      endif

      end subroutine compute_flow_residual
