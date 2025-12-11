       subroutine fem_shapefunctions_3r
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
! for each hex element in three dimensions. Called when 'fem'
! computations are used in 3-D.
! 
! Author : Sai Rapaka
!

       use comai, only: icnl, nei, ns
       use combi, only: cord
       use comfem

       implicit none

       real*8, dimension(8)   :: N, dNdZeta, dNdEta, dNdMu
       real*8, dimension(8)   :: xx, yy, zz
       real*8, dimension(3,3) :: invJ

       real*8 zeta, eta, mu
       real*8 dXdZeta, dXdEta, dXdMu
       real*8 dYdZeta, dYdEta, dYdMu
       real*8 dZdZeta, dZdEta, dZdMu
       real*8 m11, m12, m13, m21, m22, m23, m31, m32, m33
       real*8 o8, dum
       integer i,j,k

       if(icnl.ne.0) return
       if(ns.ne.8)   return
       !if(ifem.ne.1) return

       do i=1,nei

         do k=1,8
           xx(k) = cord(elnode(i,k),1)
           yy(k) = cord(elnode(i,k),2)
           zz(k) = cord(elnode(i,k),3)
         enddo

         do j=1,numgausspoints
           zeta = gpcord(j,1)
           eta  = gpcord(j,2)
           mu   = gpcord(j,3)

           o8 = 1.0d0/8.0d0
           N(1) = o8*(1 - zeta)*(1 - eta)*(1 - mu)
           N(2) = o8*(1 + zeta)*(1 - eta)*(1 - mu)
           N(3) = o8*(1 + zeta)*(1 + eta)*(1 - mu)
           N(4) = o8*(1 - zeta)*(1 + eta)*(1 - mu)
           N(5) = o8*(1 - zeta)*(1 - eta)*(1 + mu)
           N(6) = o8*(1 + zeta)*(1 - eta)*(1 + mu)
           N(7) = o8*(1 + zeta)*(1 + eta)*(1 + mu)
           N(8) = o8*(1 - zeta)*(1 + eta)*(1 + mu)

           dNdZeta(1) = -o8*(1 - eta)*(1 - mu)
           dNdZeta(2) = +o8*(1 - eta)*(1 - mu)
           dNdZeta(3) = +o8*(1 + eta)*(1 - mu)
           dNdZeta(4) = -o8*(1 + eta)*(1 - mu)
           dNdZeta(5) = -o8*(1 - eta)*(1 + mu)
           dNdZeta(6) = +o8*(1 - eta)*(1 + mu)
           dNdZeta(7) = +o8*(1 + eta)*(1 + mu)
           dNdZeta(8) = -o8*(1 + eta)*(1 + mu)

           dNdEta(1) = -o8*(1 - zeta)*(1 - mu)
           dNdEta(2) = -o8*(1 + zeta)*(1 - mu)
           dNdEta(3) = +o8*(1 + zeta)*(1 - mu)
           dNdEta(4) = +o8*(1 - zeta)*(1 - mu)
           dNdEta(5) = -o8*(1 - zeta)*(1 + mu)
           dNdEta(6) = -o8*(1 + zeta)*(1 + mu)
           dNdEta(7) = +o8*(1 + zeta)*(1 + mu)
           dNdEta(8) = +o8*(1 - zeta)*(1 + mu)

           dNdMu(1) = -o8*(1 - zeta)*(1 - eta)
           dNdMu(2) = -o8*(1 + zeta)*(1 - eta)
           dNdMu(3) = -o8*(1 + zeta)*(1 + eta)
           dNdMu(4) = -o8*(1 - zeta)*(1 + eta)
           dNdMu(5) = +o8*(1 - zeta)*(1 - eta)
           dNdMu(6) = +o8*(1 + zeta)*(1 - eta)
           dNdMu(7) = +o8*(1 + zeta)*(1 + eta)
           dNdMu(8) = +o8*(1 - zeta)*(1 + eta)

           dXdZeta = 0.0d0
           dXdEta  = 0.0d0
           dXdMu   = 0.0d0
           dYdZeta = 0.0d0
           dYdEta  = 0.0d0
           dYdMu   = 0.0d0
           dZdZeta = 0.0d0
           dZdEta  = 0.0d0
           dZdMu   = 0.0d0
           
           do k=1,8
             dXdZeta = dXdZeta + xx(k)*dNdZeta(k)
             dXdEta  = dXdEta  + xx(k)*dNdEta(k)
             dXdMu   = dXdMu   + xx(k)*dNdMu(k)

             dYdZeta = dYdZeta + yy(k)*dNdZeta(k)
             dYdEta  = dYdEta  + yy(k)*dNdEta(k)
             dYdMu   = dYdMu   + yy(k)*dNdMu(k)

             dZdZeta = dZdZeta + zz(k)*dNdZeta(k)
             dZdEta  = dZdEta  + zz(k)*dNdEta(k)
             dZdMu   = dZdMu   + zz(k)*dNdMu(k)
           enddo

           dum = dXdZeta*(dYdEta*dZdMu - dZdEta*dYdMu)
           dum = dum + dYdZeta*(dZdEta*dXdMu - dXdEta*dZdMu)
           dum = dum + dZdZeta*(dXdEta*dYdMu - dYdEta*dXdMu)
           detJ(i,j) = dum

           ! minors
           m11 = dYdEta*dZdMu - dZdEta*dYdMu
           m12 = dXdEta*dZdMu - dZdEta*dXdMu
           m13 = dXdEta*dYdMu - dYdEta*dXdMu
           m21 = dYdZeta*dZdMu - dZdZeta*dYdMu
           m22 = dXdZeta*dZdMu - dZdZeta*dXdMu
           m23 = dXdZeta*dYdMu - dYdZeta*dXdMu
           m31 = dYdZeta*dZdEta - dZdZeta*dYdEta
           m32 = dXdZeta*dZdEta - dZdZeta*dXdEta
           m33 = dXdZeta*dYdEta - dYdZeta*dXdEta

           ! inverse of jacobian matrix
           invJ(1,1) = +m11
           invJ(1,2) = -m21
           invJ(1,3) = +m31
           invJ(2,1) = -m12
           invJ(2,2) = +m22
           invJ(2,3) = -m32
           invJ(3,1) = +m13
           invJ(3,2) = -m23
           invJ(3,3) = +m33

           invJ = invJ/detJ(i,j)

           do k = 1,ns
             Psi(i,j,k) = N(k)
             dPsidX(i,j,k) = invJ(1,1)*dNdZeta(k) + invJ(1,2)*dNdEta(k)
     &                       + invJ(1,3)*dNdMu(k)
             dPsidY(i,j,k) = invJ(2,1)*dNdZeta(k) + invJ(2,2)*dNdEta(k)
     &                       + invJ(2,3)*dNdMu(k)
             dPsidZ(i,j,k) = invJ(3,1)*dNdZeta(k) + invJ(3,2)*dNdEta(k)
     &                       + invJ(3,3)*dNdMu(k)
           enddo

         enddo

       enddo

       end subroutine fem_shapefunctions_3r
