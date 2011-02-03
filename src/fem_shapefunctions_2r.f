       subroutine fem_shapefunctions_2r
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
! Computes the shapefunctions and their derivatives at each gausspoint
! for each quadrilateral element in two dimensions. Called when 'fem'
! computations are used in 2-D.
! 
! Author : Sai Rapaka
!

       use comai, only: icnl, nei, ns
       use combi, only: cord
       use comfem

       implicit none

       real*8, dimension(4)   :: N, dNdZeta, dNdEta, dNdMu
       real*8, dimension(4)   :: xx, yy
       real*8, dimension(2,2) :: invJ

       real*8 zeta, eta, mu
       real*8 dXdZeta, dXdEta
       real*8 dYdZeta, dYdEta
       real*8 m11, m12, m13, m21, m22, m23, m31, m32, m33
       real*8 o4, dum
       integer i,j,k

       if(icnl.ne.1) return
       if(ns.ne.4)   return
       !if(ifem.ne.1) return

       do i=1,nei

         do k=1,4
           xx(k) = cord(elnode(i,k),1)
           yy(k) = cord(elnode(i,k),2)
         enddo

         do j=1,numgausspoints
           zeta = gpcord(j,1)
           eta  = gpcord(j,2)

           o4 = 1.0d0/4.0d0
           N(1) = o4*(1 - zeta)*(1 - eta)
           N(2) = o4*(1 + zeta)*(1 - eta)
           N(3) = o4*(1 + zeta)*(1 + eta)
           N(4) = o4*(1 - zeta)*(1 + eta)

           dNdZeta(1) = -o4*(1 - eta)
           dNdZeta(2) = +o4*(1 - eta)
           dNdZeta(3) = +o4*(1 + eta)
           dNdZeta(4) = -o4*(1 + eta)

           dNdEta(1) = -o4*(1 - zeta)
           dNdEta(2) = -o4*(1 + zeta)
           dNdEta(3) = +o4*(1 + zeta)
           dNdEta(4) = +o4*(1 - zeta)

           dXdZeta = 0.0d0
           dXdEta  = 0.0d0
           dYdZeta = 0.0d0
           dYdEta  = 0.0d0
           
           do k=1,4
             dXdZeta = dXdZeta + xx(k)*dNdZeta(k)
             dXdEta  = dXdEta  + xx(k)*dNdEta(k)

             dYdZeta = dYdZeta + yy(k)*dNdZeta(k)
             dYdEta  = dYdEta  + yy(k)*dNdEta(k)
           enddo

           dum = dXdZeta*dYdEta - dXdEta*dYdZeta
           detJ(i,j) = dum

           ! minors
           m11 = dYdEta
           m12 = dXdEta
           m21 = dYdZeta
           m22 = dXdZeta

           ! inverse of jacobian matrix
           invJ(1,1) = +m11
           invJ(1,2) = -m21
           invJ(2,1) = -m12
           invJ(2,2) = +m22

           invJ = invJ/detJ(i,j)

           do k = 1,ns
             Psi(i,j,k) = N(k)
             dPsidX(i,j,k) = invJ(1,1)*dNdZeta(k) + invJ(1,2)*dNdEta(k)
             dPsidY(i,j,k) = invJ(2,1)*dNdZeta(k) + invJ(2,2)*dNdEta(k)
           enddo

         enddo

       enddo

       end subroutine fem_shapefunctions_2r
