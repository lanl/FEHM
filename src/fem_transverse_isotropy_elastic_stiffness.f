      subroutine fem_transverse_isotropy_elastic_stiffness(i, j, D)
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
! Returns the 6x6 material stiffness (D) matrix for a transverse isotropic
! elastic solid when the 'plastic' submacro is used with 'fem' computations
!
! Author : Satish Karra
!
      
      use comai, only: iout, iptty, ns
      use comsi, only: modelNumber, plasticModel, elastic_mod, poisson
      use comsi, only: iPlastic
      use comsi, only: elastic_mod_t, poisson_t, shearmod_t
      use comfem

      implicit none
      integer                      :: i, j
      real*8,  dimension(6,6)      :: D

      integer                      :: k, itmp, iModel, node
      real*8                       :: e1bar, e2bar, e3bar, e4bar, ezzbar
      real*8                       :: young_p, young_t, pois_p, pois_t
      real*8                       :: shear_t, shear_t_bar, pois_sq
      real*8                       :: fac1, fac2, fac3, e1, e2, e3 
      real*8                       :: e4, ezz
    
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
      e4bar = 0.0d0
      shear_t_bar = 0.0d0
      
      do k = 1,ns
!     Notation
!     e1=c11, e2=c12=c21, e3=c66=Gp, e4=c13, and ezz=c33
!     these goto isotropic limit when Ep=Et and Nue-p=Nue-t
        young_p = elastic_mod(node)
        young_t = elastic_mod_t(node)
        pois_p = poisson(node)
        pois_t = poisson_t(node)
        shear_t = shearmod_t(node)
        pois_sq = pois_t*pois_t
        fac1 = young_p/young_t
        fac2 = 1.0d0 - pois_p - 2.0d0*fac1*pois_sq
        fac3 = fac2*(1.0d0 + pois_p)
        e1 = young_p*(1.0d0 - fac1*pois_sq)/fac3
        e2 = young_p*(pois_p + fac1*pois_sq)/fac3
        e3 = 0.5d0*young_p/(1.d0 + pois_p)
        e4 = young_p*pois_t/fac2
        ezz = young_t*(1.0d0 - pois_p)/fac2
        
        e1bar = e1bar + Psi(i, j, k)*e1
        e2bar = e2bar + Psi(i, j, k)*e2
        e3bar = e3bar + Psi(i, j, k)*e3
        e4bar = e4bar + Psi(i, j, k)*e4
        ezzbar = ezzbar + Psi(i, j, k)*ezz
        shear_t_bar = shear_t_bar + Psi(i, j, k)*shear_t

      enddo

      D(1,1) = e1bar
      D(2,2) = e1bar

      D(3,3) = ezzbar
      
      D(1,2) = e2bar
      D(2,1) = e2bar

      D(1,3) = e4bar
      D(3,1) = e4bar
      D(2,3) = e4bar
      D(3,2) = e4bar

      D(4,4) = e3bar

      D(5,5) = shear_t_bar
      D(6,6) = shear_t_bar

      end subroutine fem_transverse_isotropy_elastic_stiffness
