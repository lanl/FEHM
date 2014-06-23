      subroutine  geneq_hyd_equil ( i )
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine generates the equations for 3-dimensional heat
CD1  and mass transfer for methane componet of a mixture 
CD1  Equilibrium model
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 25-FEB-03    G. Zyvoloski           Initial implementation.
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/geneq_hyd_equil.f_a  $
!D2
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

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
      real*8 dis2,dis_tol,sx_min,sx1dl,sx1dv
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

c modify water mass and energy balance
      sx1dl = sx1d*fracl

      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1dl*denhydi(i)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1dl*denehydi(i)
c  derivatives wrt pressure and temperature (water mass eq)

      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1dl*dmpf(i)
      a(jmia+nmat(2))=a(jmia+nmat(2))+sx1dl*dmef(i)
c  derivatives wrt pressure and temperature (water energy eq)

      a(jmia+nmat(7))=a(jmia+nmat(7))+sx1dl*depf(i)
      a(jmia+nmat(8))=a(jmia+nmat(8))+sx1dl*deef(i)
c  derivatives wrt water fraction

      a(jmia+nmat(5))=a(jmia+nmat(5))+sx1dl*dmwf(i) 
      a(jmia+nmat(11))=a(jmia+nmat(11))+sx1dl*dewf(i)
c  derivatives wrt hydrate fraction

      a(jmia+nmat(6))=a(jmia+nmat(6))+sx1dl*dmmf(i) 
      a(jmia+nmat(12))=a(jmia+nmat(12))+sx1dl*demf(i) 

c modify methane mass and energy balance
      sx1dv = sx1d*fracv

      bp(iz+nrhs(3))=bp(iz+nrhs(3))+sx1dv*denhydi(i)
      bp(iz+nrhs(4))=bp(iz+nrhs(4))+sx1dv*denehydi(i)
c  derivatives wrt pressure and temperature (methane mass eq)

      a(jmia+nmat(15))=a(jmia+nmat(15))+sx1dv*dmpf(i)
      a(jmia+nmat(16))=a(jmia+nmat(16))+sx1dv*dmef(i)
c  derivatives wrt pressure and temperature (methane energy eq)

      a(jmia+nmat(21))=a(jmia+nmat(21))+sx1dv*depf(i)
      a(jmia+nmat(22))=a(jmia+nmat(22))+sx1dv*deef(i)
c  derivatives wrt water fraction

      a(jmia+nmat(17))=a(jmia+nmat(17))+sx1dv*dmwf(i) 
      a(jmia+nmat(23))=a(jmia+nmat(23))+sx1dv*dewf(i)

c  derivatives wrt hydrate fraction

      a(jmia+nmat(18))=a(jmia+nmat(18))+sx1dv*dmmf(i) 
      a(jmia+nmat(24))=a(jmia+nmat(24))+sx1dv*demf(i) 
      r e t u r n
      e    n    d
