      subroutine  gdpm_geneqh ( i )
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
CD1  This subroutine generates the equations for 2-3-dimensional heat
CD1  transfer FOR RATE-LIMITED PROCESSES.  
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.

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
      implicit none

      logical bit
      integer i, iz4m1
      integer ial, iau, icd, ii1, ii2, idg, ij, ij1, ij2, isl, iq, iz
      integer jm, jmi, jmia, jml, kb, kz, ig, kbg
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
      parameter(dis_tol=1.d-12)
c
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif

       thxi=val_conh
       thyi=val_conh
       thzi=val_conh

      ti=t(i)
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

      do jm=ii1,ii2
        kb=nelm(jm)+icd
        if(kb.ne.i) then
	   iq=iq+1
         it8(iq)=kb
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
         it11(iq)=jm-neqp1
         ij1=nelm(nelm(jm))+1
         ij2=nelm(nelm(jm)+1)
         do ij=ij1,ij2
          if(nelm(ij).eq.iz) then
           it12(iq)=ij-neqp1
          endif
         enddo
	  endif
      enddo
      if(icnl.eq.0) then
c
c 3-d geometry
c
         do jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw=it10(jm)
            iwd=it10(jm)
            iw =abs(iwd)
         
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
	      if(igdpm_rate_nodes(kb).ne.0) then
	       thxkb=val_conh
             thykb=val_conh
             thzkb=val_conh
	      else
             thxkb=thx(kb)
             thykb=thy(kb)
             thzkb=thz(kb)
	       thxkb=val_conh
             thykb=val_conh
             thzkb=val_conh
            endif

            sx2t=2.*thxi*thxkb/(thxi+thxkb)
            sx3t=2.*thyi*thykb/(thyi+thykb)
            sxzt=2.*thzi*thzkb/(thzi+thzkb)
           
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
             
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t+
     &           delz2/sxzt)
            else
              sx3c=sx2c*sx_mult*max(sx2t,sx3t,sxzt)
            endif
            
            
            t5(neighc)=sx3c                    
            
         enddo
      else if ( icnl.ne.0 )  t h e n

c
c 2-d geometry
c
         radi=cord(iz,3)
         do jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw=it10(jm)
            iwd=it10(jm)
            iw =abs(iwd)
         
             radkb=0.5*(radi+cord(kz,3))
             sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
	      if(igdpm_rate_nodes(kb).ne.0) then
	       thxkb=val_conh
             thykb=val_conh
	      else
             thxkb=dis_tol
             thykb=dis_tol
            endif
             sx2t=2.*thxi*thxkb/(thxi+thxkb)
             sx3t=2.*thyi*thykb/(thyi+thykb)
           
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
             
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               sx3c=sx2c*dis2/
     &          (delx2/sx2t+dely2/sx3t) 
            else
               sx3c=sx2c*sx_mult*max(sx2t,sx3t)
            endif
            
            
            t5(neighc)=sx3c                    
            
         enddo

      endif

c
c add heat conduction
c
      do 66 jm=1,iq
      kb=it8(jm)
      kz=kb-icd
      neighc=it9(jm)
      iau=it11(jm)
      ial=it12(jm)
      jml=nelmdg(kz)-neqp1
      heatc=t5(neighc)
      bp(iz+nrhs(2))=bp(iz+nrhs(2))+heatc*(t(kb)-ti)      
      a(jmia+nmat(4))=a(jmia+nmat(4))-heatc*dtpae(i)     
      a(iau+nmat(4))=a(iau+nmat(4))+heatc*dtpae(kb)
	if(igdpm_rate_nodes(kb).eq.0) then
c add terms to other gridblocks
       bp(kz+nrhs(2))=bp(kz+nrhs(2))-heatc*(t(kb)-ti) 
       a(jml+nmat(4))=a(jml+nmat(4))-heatc*dtpae(kb)
	 a(ial+nmat(4))=a(ial+nmat(4))+heatc*dtpae(i)     
	endif      
66    continue

      r e t u r n
      e    n    d
