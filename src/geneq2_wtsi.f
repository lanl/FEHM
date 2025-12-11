      subroutine geneq2_wtsi(i) 
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To generate equations isothermal air-water solution at each node.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 22-May-02, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/geneq2_wtsi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.3.3 Noncondensible gas flow equations        
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************
c
c generate equations for 3-d schemes'isothermal air-water ,
c full derivatives'
c
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comii
      use comdti
      use comai
      use comwt
      use comriv
      implicit none

      integer i,imm
      integer icd
      integer ii1
      integer ii2
      integer idg
      integer iq
      integer jmi
      integer jml
      integer jmia
      integer jm
      integer neqp1
      integer ij
      integer ij1
      integer ij2
      integer iz
      integer kb,iunsat,kb_unsat
      integer neighc
      integer iau
      integer ial
      integer kz
      integer nmatavw,dry_nodes
      real*8 sx1d
      real*8 axi
      real*8 ayi
      real*8 azi
      real*8 alxi
      real*8 alyi
      real*8 alzi
      real*8 avxi
      real*8 avyi
      real*8 avzi
      real*8 pvii
      real*8 phii
      real*8 swi
      real*8 dlpi
      real*8 dlpkb
      real*8 dvpi
      real*8 dvpkb
      real*8 dlei,df
      real*8 dvei
      real*8 dlekb
      real*8 dvekb
      real*8 dlapi
      real*8 dvapi
      real*8 dlapkb
      real*8 dvapkb
      real*8 dili
      real*8 dilkb
      real*8 divi,dz
      real*8 divkb
      real*8 axkb
      real*8 aykb
      real*8 azkb
      real*8 alxkb
      real*8 alykb
      real*8 alzkb
      real*8 sx2c
      real*8 sx4d
c     real*8 sxzc
      real*8 pvikb
      real*8 phikb
      real*8 radi
      real*8 radkb
      real*8 fid
      real*8 fid1
      real*8 axyd,axyds
      real*8 axy
      real*8 axyf
      real*8 heatc
      real*8 pvxy
      real*8 pxy
      real*8 pxyh
      real*8 pxyi
      real*8 sx4h
      real*8 vxyd
      real*8 vxy
      real*8 vxyf
      real*8 dpvti
      real*8 dilpi
      real*8 dilei
      real*8 divpi
      real*8 divei
      real*8 dilpkb
      real*8 divpkb
      real*8 divekb
      real*8 dilekb
      real*8 sx2t
      real*8 sx3t
      real*8 sxzt
      real*8 dlaei
      real*8 dlaekb
      real*8 dvaei
      real*8 dvaekb
      real*8 ti
      real*8 dis2,dis_tol,sx_min,dxy2
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 reduction_factor
      real*8 grav_air,s_wt_fac,ds_fac_wti,ds_fac_wtkb
      real*8 rlpfree,dfid,dfid1,phi_0,cap_fac
      real*8 rlxyd, drlxyi, drlxykb, dfi, dfkb, tiny
      parameter(dis_tol=1.d-12, cap_fac = 0.0d00, tiny = 1.d-20)

      logical bit, test_bit
c     integer isl,isw_term
      integer isl
      integer iz4m1
      integer imd,iwd

c
c additional terms 
c
      real*8 tcap1(nn),tcap2(nn),dtcap1(nn),dtcap2(nn)
      real*8 pcai, axycap, pcaij,pcakb
c
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref 
c gaz 112021 pref, now in comai and possible change with air table      
      phi_0 = pref + phi_inc
c     water column height correction
      pcai = cap_fac*dzrg(i)*(rlxyf(i)-0.5)

c changed by avw -- entered here by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif
      if(icons.le.abs(maxit)) then
       grav_air=0.0
      else
       grav_air=grav
      endif

c
c storage for upwind
c

      sx1d=sx1(i)
      axi=pnx(i)
      ayi=pny(i)
      azi=pnz(i)
      alxi=axi
      avxi=axi
      alyi=ayi
      avyi=ayi
      alzi=azi
      avzi=azi
      pvii=phi(i)
      dili=dil(i)
      dilpi=dilp(i)
      if(irdof.ne.13) then
        phii=pvii-pcp(i)
        dpvti=dpcef(i)
        divi=div(i)
        divpi=divp(i)
        divei=dive(i)
        dilei=dile(i)
      else
        phii=pvii
      endif
      ti=t(i)
      swi=s(i)
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
      iz4m1 = 4*(iz-1)+1
c
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1

c Take care of source/sink term by putting it into empty slot (from 
c same node to same node). 
c If this is an isothermal air-water simulation
      a_axy(jmia+nmatavw)=sk(i)
c Take care of source/sink term 
c If this is an isothermal air-water simulation
      if (irdof .ne. 13) then
         a_vxy(jmia+nmatavw)=qh(i)
      end if

      do 58 jm=jmi+1,ii2
         kb=nelm(jm)+icd
         if(iriver.eq.2.and.kb.gt.neq_primary) go to 58
         iq=iq+1
         it8(iq)=kb
         it9(iq)=jm-ii1+1
         it10(iq)=istrw(jm-neqp1)
c gaz 11-18-2001
c         if(imdnode.ne.0) then
c           imd = mdnodes(i) + mdnodes(kb)
c           if(imd.lt.2) it10(iq) = -abs(it10(iq))
c         endif
         it11(iq)=jm-neqp1
         ij1=nelm(nelm(jm))+1
         ij2=nelmdg(nelm(jm))-1
         do 68 ij=ij1,ij2
            if(nelm(ij).eq.iz) then
               it12(iq)=ij-neqp1
            endif
 68      continue
 58   continue
c
c 3-d geometry
c
      if(icnl.eq.0) then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iwd=it10(jm)
            iw =abs(iwd)
            axkb=pnx(kb)
            aykb=pny(kb)
            azkb=pnz(kb)
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
            reduction_factor = red_factor(istrw_itfc(it11(jm)))
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)

            pvikb=phi(kb)
            if (irdof .ne. 13) then
               phikb=pvikb-pcp(kb)
            else
               phikb=pvikb
            end if
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
            endif
	    pxy = pxy*reduction_factor
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue
c     
c     
c 2-d geometry 
c
      elseif(icnl.ne.0) then
         radi=cord(iz,3)
         do 69 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iwd=it10(jm)
            iw = abs(iwd)
            neighc=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
            reduction_factor = red_factor(istrw_itfc(it11(jm)))
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
            pvikb=phi(kb)
            if (irdof .ne. 13) then
               phikb=pvikb-pcp(kb)
            else
               phikb=pvikb
            end if
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,igrav)-cord(iz,igrav))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.and.iwd.gt.0) then
              pxy=sx2c*dis2/(delx2/perml(1)+
     &           dely2/perml(2))
            else
              pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
	    pxy = pxy*reduction_factor            
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 69      continue
      endif
c     
c     liquid phase calculations
       

         do 70 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            pxyi=t1(neighc)
		  pxy=t3(neighc)
            sx4d=t6(neighc)
            pcakb = cap_fac*dzrg(kb)*(rlxyf(kb)-0.5)
            axyd=pxyi+0.5*sx4d
     &           *(rolf(i)+rolf(kb))
     &           *((cord(kz,igrav)-cord(iz,igrav))   
     &           +(pcakb-pcai))        
            t8(neighc)=axyd
      
 70      continue
c     
c     determine upwind nodes and if liquid phase exists
c     wtsi_isot = 1 isotropic; wtsi_isot = 0 anisotropic
         isl=1
         do 71 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            if(wtsi_isot.ne.0) then
c calculate information needed for direction cosines
               delx2=(cord(kz,1)-cord(iz,1))**2
               if(icnl.eq.0) then
c 3-d problem
                  dely2=(cord(kz,2)-cord(iz,2))**2
               else
c 2-d problem
                  dely2=0.
               endif
               delz2=(cord(kz,igrav)-cord(iz,igrav))**2
               dxy2=delx2+dely2
			   dis2=delx2+dely2+delz2+dis_tol		
            endif
            if(iad.le.iad_up) then
               axyd=t8(neighc)
               if(axyd.gt.0.0) then
                  imm=iz
               else
                  imm=kz
               endif

               if(wtsi_isot.eq.2) then
c
c calculate harmonic weighted relative permeability
c   

                  rlxyd = 2.*rlxyf(iz)*rlxyf(kz)/(rlxyf(iz)+rlxyf(kz)
     &                 +tiny)
                  drlxyi = 2.*drlxyf(iz)*rlxyf(kz)/(rlxyf(iz)+rlxyf(kz)
     &                 +tiny) - 2.*rlxyf(iz)*rlxyf(kz)/(rlxyf(iz)
     &                 +rlxyf(kz)+tiny)**2*drlxyf(iz)
                  drlxykb= 2.*drlxyf(kz)*rlxyf(iz)/(rlxyf(iz)+rlxyf(kz)
     &                 +tiny) - 2.*rlxyf(iz)*rlxyf(kz)/(rlxyf(iz)
     &                 +rlxyf(kz)+tiny)**2*drlxyf(kz)
c not correct
                  rlpfree=rlxyd*(dxy2/dis2)+rlzf(imm)*(delz2/dis2)
                  dfi=drlxyi*(dxy2/dis2)+drlzf(imm)*(delz2/dis2)
                  dfkb=drlxykb*(dxy2/dis2)+drlzf(imm)*(delz2/dis2)
               else if(wtsi_isot.eq.1) then
	
                  rlpfree=rlxyf(imm)*(dxy2/dis2)+rlzf(imm)*(delz2/dis2)
                  df=drlxyf(imm)*(dxy2/dis2)+drlzf(imm)*(delz2/dis2)

               else
                  rlpfree=rlxyf(imm)	
                  df = drlxyf(imm)	
               endif

               if(axyd.gt.0.0) then
c water is flowing out of node i
                  if(wtsi_isot.eq.2) then
                     t91(neighc)=rlpfree
                     t9(neighc)=0.0	         
                     dfid1f(neighc)=dfi
                     dfidf(neighc)=dfkb
                  else
                     t91(neighc)=rlpfree
                     t9(neighc)=0.0
                     dfid1f(neighc)=df
                     dfidf(neighc)=0.0
                  end if
               else
c water is flowing into node i
                  if(wtsi_isot.eq.2) then
                     t9(neighc)=rlpfree
                     t91(neighc)=0.0  
                     dfid1f(neighc)=dfi
                     dfidf(neighc)=dfkb 
                  else 
                     t9(neighc)=rlpfree
                     t91(neighc)=0.0  
                     dfid1f(neighc)=0.0
                     dfidf(neighc)=df
                  end if
               end if

               call setbit(nbits,neighc,upwind_v(iz4m1),fid)
            else
c placeholder if_block:  we should never go here
               test_bit = bit(nbits,neighc,upwind_v(iz4m1))
               if(.not.test_bit) then 
c placeholder if_block:  we should never go here
               endif
            endif
 71      continue
c     
c     form equations
c     
         if(isl.ne.0) then
            do 72 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               iau=it11(jm)
               ial=it12(jm)
               jml=nelmdg(kz)-neqp1
               heatc=t5(neighc)
               axyd=t8(neighc)

               fid=t9(neighc)
               fid1=t91(neighc)
               dfid=dfidf(neighc)
               dfid1=dfid1f(neighc)
               pxyi=t1(neighc)
               pxy=t3(neighc)
               sx4d=t6(neighc)
               dilkb=dil(kb)
               dilpkb=dilp(kb)
               pcakb = cap_fac*dzrg(kb)*(rlxyf(kb)-0.5)
	         axyd = axyd 
	         dlpi=-pxy+0.5*sx4d*(dglp(i))
     &              *((cord(kz,igrav)-cord(iz,igrav))+(pcakb-pcai)) 
     &              -0.5*sx4d*(rolf(i)+rolf(kb))
     &              *cap_fac*dzrg(i)*drlxyf(i)

               dlpkb=pxy+0.5*sx4d*(dglp(kb))
     &              *((cord(kz,igrav)-cord(iz,igrav))+(pcakb-pcai))
     &              +0.5*sx4d*(rolf(i)+rolf(kb))
     &              *cap_fac*dzrg(kb)*drlxyf(kb)

c
c put derivative wrt free surface in
c

               axyf=(fid*dilkb+fid1*dili)
               axy=axyd*axyf 
c
               dlapi=(dlpi*axyf+axyd*fid1*dilpi
     &         + axyd*dili*dfid1)

               dlapkb=(dlpkb*axyf+axyd*fid*dilpkb
     &         + axyd*dilkb*dfid) 
               a_axy(iau+nmatavw)=axy
               a_axy(ial+nmatavw)=-axy

               bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
               bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy
               a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
               a(ial+nmat(1))=a(ial+nmat(1))-dlapi
               a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
               a(jml+nmat(1))=a(jml+nmat(1))-dlapkb  
               if(kb.gt.neq_primary) then
c                write(*,*) i, kb, jmia, ial, iau,jml, axy, dlapi,dlapkb
c                pause
               endif
 72         continue
	endif
      return
      end
