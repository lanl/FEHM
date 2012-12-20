      subroutine update_permfactors()
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
! Author : Sai Rapaka
!!    Update the permfactors for all elements (entire domain)
!!    Calls the routine compute_permfactor for each element.

      use comai, only: nei, neq, iout
      use combi, only: nelm
      use comsi, only: ispm, ispmd, perx_m, pery_m, perz_m, ipermstr2
      use comsi, only: spm7f, spm8f, spm9f, spm10f
      use comsi, only: spm1f, spm2f, spm3f
      use comsi, only: iPlastic, modelNumber, plasticModel
      use comsi, only: flag_pstrain_perm_coupling, idof_stress
      use comfem

      implicit none
      double precision                   :: duu, dvv, dww
      double precision                   :: dpp
      double precision                   :: x_max, y_max, z_max
      integer id, i, j, node_j
      integer i_begin, i_end, el
      integer flag_u_pp
      integer itmp, iModel
      logical recompute_stress
      integer iispmd

      permfactor = 0.0d0
c forcing the entire domain to have the same plastic model
      if(iPlastic.eq.1) then
         iModel = modelNumber (1)
         itmp = plasticModel(iModel)
      endif

      ! Read the maximum permeability multipliers from the permmodel input
      if(ipermstr2.ne.0) then
         perx_m = spm7f(1)
         pery_m = spm8f(1)
         perz_m = spm9f(1)
      endif

      if(ispmd.eq.22) then
        perx_m = spm8f(1)
        pery_m = spm9f(1)
        perz_m = spm10f(1)
      endif

      if(ispmd.eq.100) then
        perx_m = spm1f(1)
        pery_m = spm2f(1)
        perz_m = spm3f(1)
      endif

      !! Accumulate permfactor for all the connections, element
      !! by element
      do el=1,nei
         recompute_stress = .False.
        node_j = 0
        duu = 0.0d0; dvv = 0.0d0; dww = 0.0d0; dpp = 0.0d0
        flag_u_pp = 0
        if(iPlastic.eq.0) then
c        call compute_permfactor(el, node_j, duu, dvv, dww
c     &          , recompute_stress)
c        call compute_permfactor_effstrs(el, node_j, duu, dvv, dww
c     &       ,dpp,flag_u_pp, recompute_stress)
           call compute_permfactor_effstrs_pp(el,node_j,duu, dvv, dww
     &          ,dpp,flag_u_pp, recompute_stress)
        elseif(flag_pstrain_perm_coupling.eq.1) then
c           elseif(itmp.eq.3) then
c              call compute_permfactor_DruckerPrager(el, node_j, duu, dvv, 
c     &             dww, dpp, flag_u_pp)
c the same routine handles V-M and D-P models
           if(itmp.eq.2.or.itmp.eq.3) then
              call compute_permfactor_vonMises(el, node_j, duu, dvv, dww
     &           ,dpp, flag_u_pp)
           endif
        endif
      enddo 

      !! Normalize by number of elements sharing a connection
      do i=1,nelm(neq+1)
        if(numelems(i).gt.0) then
          permfactor(i,1) = permfactor(i,1)/numelems(i)
          permfactor(i,2) = permfactor(i,2)/numelems(i)
          permfactor(i,3) = permfactor(i,3)/numelems(i)

          permfactor(i,1) = min(permfactor(i,1), perx_m)
c          permfactor(i,2) = min(permfactor(i,2), pery_m)
c          permfactor(i,3) = min(permfactor(i,3), perz_m)
          permfactor(i,2) = min(permfactor(i,2), perx_m)
          permfactor(i,3) = min(permfactor(i,3), perx_m)

          permfactor(i,1) = max(permfactor(i,1), 1.0d0)
          permfactor(i,2) = max(permfactor(i,2), 1.0d0)
          permfactor(i,3) = max(permfactor(i,3), 1.0d0)
        else
          !! Edge from node i to itself
          permfactor(i,:) = 1.0d0
        endif
        if(permfactor(i,1).le.0.0) permfactor(i,1) = 1.0d0
        if(permfactor(i,2).le.0.0) permfactor(i,1) = 1.0d0
        if(permfactor(i,3).le.0.0) permfactor(i,1) = 1.0d0
      enddo

      end subroutine update_permfactors

