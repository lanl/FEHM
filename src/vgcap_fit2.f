      subroutine vgcap_fit2(iflg,mm,iphase)
!***********************************************************************
! Copyright 2009 Los Alamos National Security, LLC  All rights reserved
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
!***********************************************************************

C***********************************************************************
c iflg =0 normal
c iflg =1 fracture
c mm= model
c iphase = phase
c
      use comai, only : nrlp
      use comrlp

      implicit none

      integer iflg,mm,iphase,k,ncut
      real*8 slr,slm,slcut,alpha,alamda,c1,c2,c3,c4 
      real*8 hmax,hmin,hcut,smcut,alpi,ds,dhp,slope
      real*8 fac,fac_min,fac_use
      real*8 su_cut,ss,smcutm,smcutf
      real*8 alpham, alamdam, facm, ancut
      parameter(hmin = 1.d-18)
      parameter(su_cut = 0.99d00)
c      parameter(su_cut = 0.70d00)
      parameter(fac_min=2.0)
      parameter(ncut = 1)

      k = (iphase - 1) * max_cp
      c1 = 0.
      c2 = 0.
      c3 = 0.
      c4 = 0.
      if(iflg.eq.0) then
! regular
         if (.not. allocated(vg1)) then
            allocate (vg1(nrlp,2),vg2(nrlp,2),vg3(nrlp,2),vg4(nrlp,2))
            allocate (cp1s(nrlp,2),cp2s(nrlp,2))
         end if
         alamda = 1. - 1.0 / cap_param(mm, k + 4)
         alpha = cap_param(mm, k + 3)
         slr = cap_param(mm, k + 1)
         slm = cap_param(mm, k + 2)
         smcut = cap_param(mm, k + 6)
         fac = cap_param(mm, k + 5)
         alpham = alpha
         alamdam = alamda
         facm = fac
         smcutm = smcut
      else
! fracture
         if (.not. allocated(vg5)) then
            allocate (vg5(nrlp,2),vg6(nrlp,2),vg7(nrlp,2),vg8(nrlp,2))
            allocate (cp3s(nrlp,2),cp4s(nrlp,2))
         end if
         alamdam = 1. - 1.0 / cap_param(mm, k + 4)
         alamda = 1. - 1.0 / cap_fparam(mm, k + 4)
         alpham = cap_param(mm, k + 3)
         alpha = cap_fparam(mm, k + 3)
         slr = cap_fparam(mm, k + 1)
         slm = cap_fparam(mm, k + 2)
         smcutm = cap_param(mm, k + 6)
         smcut = cap_fparam(mm, k + 6)
         facm = cap_param(mm, k + 5)
         fac = cap_fparam(mm, k + 5)
         if (facm .ge. 0.d0) then
            call vgcap_match(2, smcutm, alpham, alamdam, facm,
     &           smcut, alpha, alamda, fac)
            cap_fparam(mm, k + 5) = fac
         end if
      endif

      slcut=smcut*(slm-slr)+slr
		
      if (fac .eq. 0.) then
c zero slope, zero curvature at sl=0
c h = a*x**3 + d
c this is most common fehm way
         alpi = 1.0/alamda
         hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
         ds = 1.0/(slm-slr)
         dhp = 1.0/alpha*(1.0-alamda)/
     *        (1.0/smcut**alpi-1.0)**alamda
     *        *(-alpi/smcut**(alpi+1.0))*ds
         fac=1.0/3.0
         hmax=hcut-fac*dhp*slcut
         c1=dhp/(3.0*slcut**2)
         c2=0.0
         c3=0.0
         c4 = hmax
      else if (fac .gt. 0.) then
c linear fit from slcut to sl=0.0
c with hmax=fac*hcut
         fac_use=max(fac,fac_min)
         alpi = 1.0/alamda
         hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
         hmax=fac_use*hcut
         slope=-(hmax-hcut)/slcut
         c3=slope
         c4=hmax
      else if (fac .lt. 0.) then
c linear fit from slcut to sl=0.0
c with hmax: linear extension of slope at slcut
         alpi = 1.0/alamda
         hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
         ds = 1.0/(slm-slr)          
         dhp = 1.0/alpha*(1.0-alamda)/
     *        (1.0/smcut**alpi-1.0)**alamda
     *        *(-alpi/smcut**(alpi+1.0))*ds
         hmax=hcut-dhp*slcut
         slope=-(hmax-hcut)/slcut
         c3=slope
         c4=hmax
      endif

      if(iflg.eq.0) then
         vg1(mm,iphase)=c1
         vg2(mm,iphase)=c2
         vg3(mm,iphase)=c3
         vg4(mm,iphase)=c4
      else
         vg5(mm,iphase)=c1
         vg6(mm,iphase)=c2
         vg7(mm,iphase)=c3
         vg8(mm,iphase)=c4
      endif
c linear fit from smcut to star= 1.00
c with hmax: linear extension of slope at slcut
      slcut = su_cut * (slm  - slr) + slr
      alpi = 1.0/alamda
      hcut = 1.0/alpha*(1.0/su_cut**alpi-1.0)**(1.0-alamda)
      slope=(hcut-hmin)/(1.0-slcut)
c      c3=-slope
c      c4=slope+hmin
      ancut = (1./slcut**ncut-1.0)
      c3 = hcut/ancut
      c4 = -c3 + hmin
      if(iflg.eq.0) then
         cp1s(mm,iphase)=c3
         cp2s(mm,iphase)=c4
      else
         cp3s(mm,iphase)=c3
         cp4s(mm,iphase)=c4
      endif

      end subroutine vgcap_fit2
