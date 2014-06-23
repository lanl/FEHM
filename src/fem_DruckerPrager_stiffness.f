      subroutine fem_DruckerPrager_stiffness(i, j, D)
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
! Returns the 6x6 material stiffness (D) matrix for Drucker Prager without cap
! model 
!
! Author : Satish Karra
! Date : March 28, 2012
!
      use comsi, only: iPlastic, plasticModel, modelNumber, e1, e2, e3
      use comsi, only: isPlastic, plasticParam1, plastic_strain, du
      use comai, only: iout, iptty, ns
      use comfem
      
      implicit none
      real*8, dimension(1:6, 1:6)  :: D
      real*8, dimension(6)         :: dev_stress, N, DeN
      real*8, dimension(8)         :: lambda, G
      integer, dimension(8)        :: node
      real*8 lambda_bar, G_bar,  H, NDeN
      real*8 pressure, J2, eta
      integer i, j, k, l, m, itmp, iModel

      lambda_bar = 0.0
      G_bar = 0.0
      NDeN = 0.0
      pressure = 0.0
      


      if(iPlastic.eq.0) then
        write(iout,*) '***ERROR : elastic stiffness routine called 
     &   without plastic flag being set '
      endif
      
      itmp = modelNumber(elnode(i, 1))
      iModel = plasticModel(itmp)

      eta = plasticParam1(itmp)

      do k = 2,ns
        itmp = modelNumber(elnode(i, k))
        if(iModel.ne.plasticModel(itmp)) then
          write(iout, *) 'Multiple plastic models being used ! 
     &        Not supported at this time! '
          write(iptty, *) 'Multiple plastic models being used ! 
     &        Not supported at this time! '
          stop
        endif
      enddo

      do k=1,ns
        node(k) = elnode(i, k)
        lambda(k) = e2(node(k))
        G(k) = 0.5d0*(e1(k) - e2(k))
        lambda_bar = lambda_bar + lambda(k)*Psi(i,j,k)
        G_bar = G_bar + G(k)*Psi(i,j,k)
      enddo

      D = 0.0d0
      if(iModel.eq.3) then
        ! First initialize D to the elastic matrix
        D(1:3, 1:3) = lambda_bar
        D(1,1) = D(1,1) + 2.0d0*G_bar
        D(2,2) = D(2,2) + 2.0d0*G_bar
        D(3,3) = D(3,3) + 2.0d0*G_bar
        D(4,4) = G_bar
        D(5,5) = G_bar
        D(6,6) = G_bar

        if(isPlastic(i,j).eq.1) then
          ! Compute deviatoric stress
          dev_stress = fem_stress(i, j, :)
          pressure = sum(dev_stress(1:3))/3.0d0
          dev_stress(1) = dev_stress(1) - pressure
          dev_stress(2) = dev_stress(2) - pressure
          dev_stress(3) = dev_stress(3) - pressure

          call J2Invariant(J2, dev_stress)

          N(1) = 1/2.d0*sqrt(J2)*dev_stress(1) + eta/3
          N(2) = 1/2.d0*sqrt(J2)*dev_stress(2) + eta/3
          N(3) = 1/2.d0*sqrt(J2)*dev_stress(3) + eta/3
          N(4) = 1/2.d0*sqrt(J2)*dev_stress(4) 
          N(5) = 1/2.d0*sqrt(J2)*dev_stress(5) 
          N(6) = 1/2.d0*sqrt(J2)*dev_stress(6) 

!          N(4) = N(4)*2.0d0
!          N(5) = N(5)*2.0d0
!          N(6) = N(6)*2.0d0

          DeN = matmul(D, N)
          NDeN = sum(N*DeN)

          do l=1,6
            do m=1,6
                D(l,m) = D(l,m) - DeN(l)*DeN(m)/NDeN
            enddo
          enddo

        endif
      else
        write(iout,*) '***ERROR: Drucker Prager stiffness routine called for
     &   a node not set to be Drucker Prager'
      endif
      
      end subroutine fem_DruckerPrager_stiffness
