      subroutine compute_average_stress()
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
!
      
      use comai, only: nei, neq
      use combi, only: nelm
      use comdi, only: phi
      use comsi, only: bulk, e2, e3, ispm, ispmt
      use comfem
      
      implicit none
      
      integer                            :: el, j, n1, n2, n
      integer                            :: node_I, node_J
      integer                            :: iispmd, jispmd
      integer                            :: model_I, model_J
      integer                            :: permmodel

      real*8, dimension(6)               :: intsig
      real*8, dimension(6)               :: avg_stress
      real*8                             :: onedV, fac, biot
      real*8                             :: mean_stress, P_ij
      real*8                             :: eff_stress, xperm
      real*8                             :: bulki, bulkj, bulkbar
      real*8                             :: harmonic_mean
      
      if(.not. allocated(edgeNum1)) then
        call setup_edgePointers_3D()
      endif

      permFactor = 0.0d0

      do el=1,nei
      
        avg_stress = 0.0d0
        intsig = 0.0d0
        onedV = 0.0d0
      
        do j=1,numgausspoints
          fac = detJ(el,j)*gpweight(j)
          intsig = intsig + fem_stress(el,j,:)*fac
          onedV = onedV + fac
        enddo
      
        avg_stress = intsig/onedV
        mean_stress = sum(avg_stress(1:3))/3.0d0

        do n=1, 12
          node_I = elnode(el, edges(n, 1))
          node_J = elnode(el, edges(n, 2))

          P_ij = 0.5d0*(phi(node_I) + phi(node_J))
          biot = 0.5d0*(bulk(node_I) + bulk(node_J))

          bulki = 3.0d0*e2(node_I) + 2.0d0*e3(node_I)
          bulkj = 3.0d0*e2(node_J) + 2.0d0*e3(node_J)
          bulkbar = harmonic_mean(bulki, bulkj)

          eff_stress = mean_stress + bulkbar*biot*P_ij

        enddo
      enddo
      
      end subroutine compute_average_stress
