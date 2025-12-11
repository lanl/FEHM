      module comrlp_param
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
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
!D1
!D1 PURPOSE
!D1
!D1 Global parameters for VG regularization
!D1 GAZ 030324 Initial Coding    
!D1
!***********************************************************************

c gaz 041524      
c      real(8) :: scutm,hmin,darcyf,tol_l,tol_u,su_cut,su_cut5
      real(8) :: scutm,hmin,darcyf,tol_l,tol_u,su_cut
c gaz 040824 need scutm = 1.d-03 to pass VV  
      parameter(scutm = 1.d-03)
c      parameter(scutm = 1.d-05)
c gaz debug test hmin
c      parameter(hmin = 1.d-8)
      parameter(hmin = 1.d-10)
      parameter(darcyf = 1.d12)
c      parameter(tol_l  = 1.d-5)
c      parameter(tol_u  = 1.d-5)
c gaz debug test su_cut    
      parameter(su_cut = 0.99d00)
c      parameter(su_cut = 0.999d00)
c gaz 041524 su_cut5    not used 
c      parameter (su_cut5 = 0.99d0)
      integer :: ireg0
      integer :: ireg = 1
      
      real*8 alpha,beta,alamda,alpi,smcut,slcut,fac,ds,dhp
      real*8 rp1,rp2,rp3,rp4,denom,star,hp,rl1,rv1,hp1,dhp1
      real(8) :: rl = 1., rv = 1., drls = 0., drvs = 0.
      real*8 drls1,drvs1,akf,akm,porf,permb,sl
      real*8 smcutm,smcutf,alpham,alamdam,facf
      real*8 rpa1, rpa2, rpa3, rpa4, rpa5

      end module comrlp_param
