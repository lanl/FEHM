      subroutine geneq_stress_fem_2D_XY()
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
! Assembles the mechanical equations for 'fem' computations in two dimensions
! 
! Author : Sai Rapaka
!

      use comai, only: nei, neq, ns, igrav, grav
      use combi, only: nelm, nelmdg
      use comdi, only: denr
      use comei, only: a
      use comgi, only: bp
      use comsi, only: e1, e2, e3, iPlastic, ibodyforce
      use comfem
      use davidi,only: nmat, nrhs

      implicit none

      integer                      :: i,j,k,li,m,ii1,ii2
      integer                      :: nodei, nodej, jmia, iau
      integer, dimension(4)        :: node
      real*8,  dimension(3, 8)     :: B
      real*8,  dimension(3, 3)     :: D
      real*8,  dimension(3, 8)     :: DB
      real*8,  dimension(8, 8)     :: BtDB
      real*8,  dimension(2, 2)     :: nodalK
      real*8,  dimension(3)        :: gp_stress
      real*8,  dimension(8)        :: Btsigma
      real*8,  dimension(4)        :: bforcex, bforcey, bforcez

      real*8                       :: e1bar, e2bar, e3bar
      real*8                       :: fac

      do i=1,nei
        do k=1,ns
          node(k) = elnode(i,k)
        enddo

        do j=1,numgausspoints

          call fem_update_stress_2D(i, j)

          B = 0.0d0
          do k=1,ns
            B(1,2*(k-1) + 1) = dPsidX(i, j, k)
            B(2,2*(k-1) + 2) = dPsidY(i, j, k)
            B(3,2*(k-1) + 1) = dPsidY(i, j, k)
            B(3,2*(k-1) + 2) = dPsidX(i, j, k)
          enddo
          
          !!! Get D matrix from material module
          if(iPlastic.eq.1) then
            call fem_material_stiffness_2D(i, j, D)
          else
            e1bar = 0.0d0
            e2bar = 0.0d0
            e3bar = 0.0d0
 
            do k=1,ns
              e1bar = e1bar + Psi(i, j, k)*e1(node(k))
              e2bar = e2bar + Psi(i, j, k)*e2(node(k))
              e3bar = e3bar + Psi(i, j, k)*e3(node(k))
            enddo

            D = 0.0d0
            D(1,1) = e1bar
            D(1,2) = e2bar
            D(2,1) = e2bar
            D(2,2) = e1bar
            D(3,3) = e3bar
          endif

          DB = matmul(D, B)
          BtDB = matmul(transpose(B), DB)

          fac = detJ(i,j)*gpweight(j)

          do k=1,ns
            nodei = node(k)
            do li=1,ns
              if(li.ne.k) then
                nodej = node(li)
                ! diagonal entry
                jmia = nelmdg(nodei) - (neq + 1)
                ! (i,j) entry
                ii1 = nelm(nodei) + 1
                ii2 = nelm(nodei + 1)
                iau = 0
                do m=ii1,ii2
                  if(nelm(m) .eq. nodej) then
                    iau  = m - (neq + 1)
                  endif
                enddo
  
                nodalK = BtDB(2*k-1:2*k, 2*li-1:2*li)
                a(iau + nmat(1)) = a(iau + nmat(1)) + nodalK(1,1)*fac
                a(iau + nmat(2)) = a(iau + nmat(2)) + nodalK(1,2)*fac
                a(iau + nmat(3)) = a(iau + nmat(3)) + nodalK(2,1)*fac
                a(iau + nmat(4)) = a(iau + nmat(4)) + nodalK(2,2)*fac
  
                a(jmia + nmat(1)) = a(jmia + nmat(1)) - nodalK(1,1)*fac
                a(jmia + nmat(2)) = a(jmia + nmat(2)) - nodalK(1,2)*fac
                a(jmia + nmat(3)) = a(jmia + nmat(3)) - nodalK(2,1)*fac
                a(jmia + nmat(4)) = a(jmia + nmat(4)) - nodalK(2,2)*fac
              endif
            enddo ! inner loop over neighboring node
          enddo   ! loop over node i

          !body forces
          bforcex = 0.0d0
          bforcey = 0.0d0

          if(ibodyforce.ne.0) then
            do k=1,ns
              if(igrav.eq.1) then
                bforcex(k) = -denr(k)*grav*Psi(i,j,k)*fac
              elseif(igrav.eq.2) then
                bforcey(k) = -denr(k)*grav*Psi(i,j,k)*fac
              endif
              bp(node(k) + nrhs(1)) = bp(node(k) + nrhs(1)) + bforcex(k)
              bp(node(k) + nrhs(2)) = bp(node(k) + nrhs(2)) + bforcey(k)
            enddo
          endif

          gp_stress = fem_stress(i,j,:)
          Btsigma = matmul(transpose(B), gp_stress)

          do k=1,ns
            bp(node(k) + nrhs(1)) = bp(node(k) + nrhs(1)) + 
     &         Btsigma(2*k-1)*fac
            bp(node(k) + nrhs(2)) = bp(node(k) + nrhs(2)) + 
     &         Btsigma(2*k  )*fac
          enddo

        enddo     ! loop over gausspoints
      enddo       ! loop over elements

      end subroutine geneq_stress_fem_2D_XY
