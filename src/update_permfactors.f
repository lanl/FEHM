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
      use comsi, only: perx_m, pery_m, perz_m, ipermstr2
      use comsi, only: spm7f, spm8f, spm9f
      use comfem

      implicit none
      double precision                   :: duu, dvv, dww
      double precision                   :: dpp
      double precision                   :: x_max, y_max, z_max
      integer id, i, j, node_j
      integer i_begin, i_end, el
      integer flag_u_pp

      permfactor = 0.0d0

      ! Read the maximum permeability multipliers from the permmodel input
      if(ipermstr2.ne.0) then
        x_max = spm7f(1)
        y_max = spm8f(1)
        z_max = spm9f(1)
      endif

      !! Accumulate permfactor for all the connections, element
      !! by element
      do el=1,nei
        node_j = 0
        duu = 0.0d0; dvv = 0.0d0; dww = 0.0d0; dpp = 0.0d0
        flag_u_pp = 0
c        call compute_permfactor(el, node_j, duu, dvv, dww)
        call compute_permfactor_effstrs(el, node_j, duu, dvv, dww
     &       ,dpp,flag_u_pp)
      enddo 

      !! Normalize by number of elements sharing a connection
      do i=1,nelm(neq+1)
        if(numelems(i).gt.0) then
          permfactor(i,:) = permfactor(i,:)/numelems(i)
          permfactor(i,1) = min(permfactor(i,1), x_max)
          permfactor(i,2) = min(permfactor(i,2), y_max)
          permfactor(i,3) = min(permfactor(i,3), z_max)
          permfactor(i,1) = max(permfactor(i,1), 1.0d0)
          permfactor(i,2) = max(permfactor(i,2), 1.0d0)
          permfactor(i,3) = max(permfactor(i,3), 1.0d0)
        else
          !! Edge from node i to itself
          permfactor(i,:) = 1.0d0
        endif
      enddo

      end subroutine update_permfactors

