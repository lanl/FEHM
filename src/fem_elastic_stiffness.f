      subroutine fem_elastic_stiffness(i, j, D)
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
! Returns the 6x6 material stiffness (D) matrix for a linear, isotropic,
! elastic solid when the 'plastic' submacro is used with 'fem' computations
!
! Author : Sai Rapaka
!
      
      use comai, only: iout, iptty, ns
      use comsi, only: modelNumber, plasticModel, elastic_mod, poisson
      use comsi, only: iPlastic
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6,6)      :: D

      integer                      :: k, itmp, iModel, node
      real*8                       :: e1bar, e2bar, e3bar
      real*8                       :: E, nu, lambda, G

      if(iPlastic.ne.0) then
        itmp = modelNumber(elnode(i, 1))
        iModel = plasticModel(itmp)

        do k=2,ns
          itmp = modelNumber(elnode(i, k))
          if(iModel.ne.plasticModel(itmp)) then
            write(iout, *) 'Multiple plastic models being used ! 
     &          Not supported at this time! '
            write(iptty, *) 'Multiple plastic models being used ! 
     &          Not supported at this time! '
            stop
          endif
        enddo
      endif

      D = 0.0d0

      e1bar = 0.0d0
      e2bar = 0.0d0
      e3bar = 0.0d0

      do k=1,ns
        node = elnode(i, k)
        E = elastic_mod(node)
        nu = poisson(node)
        lambda = E*nu/((1 + nu)*(1 - 2.0d0*nu))
        G = E/(2.0d0*(1 + nu))
        
        e1bar = e1bar + Psi(i, j, k)*(lambda + 2.0d0*G)
        e2bar = e2bar + Psi(i, j, k)*lambda
        e3bar = e3bar + Psi(i, j, k)*G
      enddo

      D(1,1) = e1bar
      D(2,2) = e1bar
      D(3,3) = e1bar
      
      D(1,2) = e2bar
      D(1,3) = e2bar
      D(2,1) = e2bar
      D(2,3) = e2bar
      D(3,1) = e2bar
      D(3,2) = e2bar

      D(4,4) = e3bar
      D(5,5) = e3bar
      D(6,6) = e3bar
      
      end subroutine fem_elastic_stiffness
