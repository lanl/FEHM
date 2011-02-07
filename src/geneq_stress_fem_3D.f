      subroutine geneq_stress_fem_3D()

      use comai, only: nei, neq, ns, igrav, grav
      use combi, only: nelm, nelmdg
      use comdi, only: denr
      use comei, only: a
      use comgi, only: bp
      use comsi, only: e1, e2, e3, iPlastic,ibodyforce
      use comfem
      use davidi,only: nmat, nrhs

      implicit none

      integer                      :: i,j,k,li,m,ii1,ii2
      integer                      :: nodei, nodej, jmia, iau
      integer, dimension(8)        :: node
      real*8,  dimension(6, 24)    :: B
      real*8,  dimension(6, 6)     :: D
      real*8,  dimension(6, 24)    :: DB
      real*8,  dimension(24,24)    :: BtDB
      real*8,  dimension( 3, 3)    :: nodalK
      real*8,  dimension(6)        :: gp_stress
      real*8,  dimension(24)       :: Btsigma
      real*8,  dimension(8)        :: bforcex, bforcey, bforcez

      real*8                       :: e1bar, e2bar, e3bar
      real*8                       :: fac

      real*8,  dimension(neq)      :: sx_sai

      sx_sai = 0.0d0

      do i=1,nei
        do k=1,ns
          node(k) = elnode(i,k)
        enddo

        do j=1,numgausspoints

          call fem_update_stress(i, j)

          B = 0.0d0
          do k=1,ns
            B(1,3*(k-1) + 1) = dPsidX(i, j, k)
            B(2,3*(k-1) + 2) = dPsidY(i, j, k)
            B(3,3*(k-1) + 3) = dPsidZ(i, j, k)
            B(4,3*(k-1) + 1) = dPsidY(i, j, k)
            B(4,3*(k-1) + 2) = dPsidX(i, j, k)
            B(5,3*(k-1) + 2) = dPsidZ(i, j, k)
            B(5,3*(k-1) + 3) = dPsidY(i, j, k)
            B(6,3*(k-1) + 1) = dPsidZ(i, j, k)
            B(6,3*(k-1) + 3) = dPsidX(i, j, k)
          enddo
          
          !!! Get D matrix from material module
          if(iPlastic.eq.1) then
            call fem_material_stiffness(i, j, D)
          else
            e1bar = 0.0d0
            e2bar = 0.0d0
            e3bar = 0.0d0
 
            do k=1,ns
              e1bar = e1bar + Psi(i, j, k)*e1(node(k))
              e2bar = e2bar + Psi(i, j, k)*e2(node(k))
              e3bar = e3bar + Psi(i, j, k)*e3(node(k))
            enddo

!            e1bar = e1(node(j))
!            e2bar = e2(node(j))
!            e3bar = e3(node(j))

            D = 0.0d0
            D(1,1) = e1bar
            D(1,2) = e2bar
            D(1,3) = e2bar
            D(2,1) = e2bar
            D(2,2) = e1bar
            D(2,3) = e2bar
            D(3,1) = e2bar
            D(3,2) = e2bar
            D(3,3) = e1bar
            D(4,4) = e3bar
            D(5,5) = e3bar
            D(6,6) = e3bar
          endif

          DB = matmul(D, B)
          BtDB = matmul(transpose(B), DB)

          fac = detJ(i,j)*gpweight(j)

          do k=1,8
            nodei = node(k)
            do li=1,8
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
  
                nodalK = BtDB(3*k-2:3*k, 3*li-2:3*li)
                a(iau + nmat(1)) = a(iau + nmat(1)) + nodalK(1,1)*fac
                a(iau + nmat(2)) = a(iau + nmat(2)) + nodalK(1,2)*fac
                a(iau + nmat(3)) = a(iau + nmat(3)) + nodalK(1,3)*fac
                a(iau + nmat(4)) = a(iau + nmat(4)) + nodalK(2,1)*fac
                a(iau + nmat(5)) = a(iau + nmat(5)) + nodalK(2,2)*fac
                a(iau + nmat(6)) = a(iau + nmat(6)) + nodalK(2,3)*fac
                a(iau + nmat(7)) = a(iau + nmat(7)) + nodalK(3,1)*fac
                a(iau + nmat(8)) = a(iau + nmat(8)) + nodalK(3,2)*fac
                a(iau + nmat(9)) = a(iau + nmat(9)) + nodalK(3,3)*fac
  
                a(jmia + nmat(1)) = a(jmia + nmat(1)) - nodalK(1,1)*fac
                a(jmia + nmat(2)) = a(jmia + nmat(2)) - nodalK(1,2)*fac
                a(jmia + nmat(3)) = a(jmia + nmat(3)) - nodalK(1,3)*fac
                a(jmia + nmat(4)) = a(jmia + nmat(4)) - nodalK(2,1)*fac
                a(jmia + nmat(5)) = a(jmia + nmat(5)) - nodalK(2,2)*fac
                a(jmia + nmat(6)) = a(jmia + nmat(6)) - nodalK(2,3)*fac
                a(jmia + nmat(7)) = a(jmia + nmat(7)) - nodalK(3,1)*fac
                a(jmia + nmat(8)) = a(jmia + nmat(8)) - nodalK(3,2)*fac
                a(jmia + nmat(9)) = a(jmia + nmat(9)) - nodalK(3,3)*fac
              endif
            enddo ! inner loop over neighboring node
          enddo   ! loop over node i

          !body forces
          bforcex = 0.0d0
          bforcey = 0.0d0
          bforcez = 0.0d0

          if(ibodyforce.ne.0) then
            do k=1,8
              if(igrav.eq.1) then
                bforcex(k) = -denr(k)*grav*Psi(i,j,k)*fac
              elseif(igrav.eq.2) then
                bforcey(k) = -denr(k)*grav*Psi(i,j,k)*fac
              elseif(igrav.eq.3) then
                bforcez(k) = -denr(k)*grav*Psi(i,j,k)*fac
              endif
              bp(node(k) + nrhs(1)) = bp(node(k) + nrhs(1)) + bforcex(k)
              bp(node(k) + nrhs(2)) = bp(node(k) + nrhs(2)) + bforcey(k)
              bp(node(k) + nrhs(3)) = bp(node(k) + nrhs(3)) + bforcez(k)
            enddo
          endif

          gp_stress = fem_stress(i,j,:)
          Btsigma = matmul(transpose(B), gp_stress)

          do k=1,8
            bp(node(k) + nrhs(1)) = bp(node(k) + nrhs(1)) + 
     &         Btsigma(3*k-2)*fac
            bp(node(k) + nrhs(2)) = bp(node(k) + nrhs(2)) + 
     &         Btsigma(3*k-1)*fac
            bp(node(k) + nrhs(3)) = bp(node(k) + nrhs(3)) + 
     &         Btsigma(3*k  )*fac
          enddo

        enddo     ! loop over gausspoints
      enddo       ! loop over elements

      end subroutine geneq_stress_fem_3D
