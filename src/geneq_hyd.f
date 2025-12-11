      subroutine  geneq_hyd ( i )
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1  PURPOSE
!D1
!D1  This subroutine generates the equations for 3-dimensional heat
!D1  and mass transfer for methane componet of a mixture 
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 25-FEB-03    G. Zyvoloski           Initial implementation.
!D2
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/geneq_hyd.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.3.2 Heat- and mass-transfer equations
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comflow
      use davidi
      use comji
      use comhi
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use commeth
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz
      integer neighc, neqp1, nmatavw
      integer imd,iwd
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8  aexy, aexyf, alxi, alxkb, alyi, alykb, alzi, alzkb
      real*8  avxi, avyi, avzi, axi, axkb, axy, axyd, axyf
      real*8  ayi, aykb, azi, azkb
      real*8  delei, delekb, deli, delkb, devei, devekb, devi, devkb 
      real*8  dilei, dilekb, dili, dilkb, dilpi, dilpkb 
      real*8  divei, divekb, divi, divkb, divpi, divpkb
      real*8  dlaei, dlaekb, dlapi, dlapkb, dleei, dleekb 
      real*8  dlei, dlekb, dlepi, dlepkb, dlpi, dlpkb, dpvti  
      real*8  dvaei, dvaekb, dvapi, dvapkb, dveei, dveekb
      real*8  dvei, dvekb, dvepi, dvepkb, dvpi, dvpkb
      real*8  fid, fid1, grav_air, heatc, enli, enlkb, envi, envkb
      real*8  phii, phikb, pvii, pvikb, pvxy, pxy, pxyi, pxyh
      real*8  radi, radkb
      real*8  swi, sx1d, sx2c, sx2t, sx3c, sx3t, sx4d, sx4h, sxzt
      real*8  thxi, thxkb, thyi, thykb, thzi, thzkb, ti
      real*8  vexy, vexyf, vxy, vxyd, vxyf

      real*8  divmi, divmkb, dilmi, dilmkb                      
      real*8  dlami, dlamkb, dlemi, dlemkb                      
      real*8  dvami, dvamkb, dvemi, dvemkb                      

c
c
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif

      sx1d=sx1(i)
      phii=phihyd(i)
      enli=enlfhyd(i)
      ti=tmeth(i)
      swi=frachyd(i)
c
c form constants for i>neq
c
      if(i.gt.neq.and.idualp.eq.0) then
      icd=neq
      else
      icd=0
      endif
      iz=i-icd
c
      iz4m1 = 4*(iz-1) +1
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1

c add accumulation terms
c note no dependence on water frac (variable #3) 
c pressure is pmeth, but = pwater
c no source/sink terms but dissociation
c t=thyd

      bp(iz+nrhs(5))=bp(iz+nrhs(5))+sx1d*denhydi(i)+skhyd(i)
      bp(iz+nrhs(6))=bp(iz+nrhs(6))+sx1d*denehydi(i)+qhhyd(i)
      a(jmia+nmat(25))=a(jmia+nmat(25))+sx1d*dmpf(i)+dskhyd1(i)
      a(jmia+nmat(26))=a(jmia+nmat(26))+sx1d*dmef(i)+dskhyd2(i)
      a(jmia+nmat(31))=a(jmia+nmat(31))+sx1d*depf(i)+dqhhyd1(i)
      a(jmia+nmat(32))=a(jmia+nmat(32))+sx1d*deef(i)+dqhhyd2(i)

c  derivatives wrt hydrate fraction

      a(jmia+nmat(30))=a(jmia+nmat(30))+sx1d*dmmf(i) + dskhyd4(i)
      a(jmia+nmat(36))=a(jmia+nmat(36))+sx1d*demf(i) + dqhhyd4(i)

c     derivatives wrt water fraction
c      dmwf(i)=0.0
c      dewf(i)=0.0
      a(jmia+nmat(29))=a(jmia+nmat(29))+sx1d*dmwf(i) + dskhyd3(i)
      a(jmia+nmat(35))=a(jmia+nmat(35))+sx1d*dewf(i)

      r e t u r n
      e    n    d
