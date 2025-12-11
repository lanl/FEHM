      subroutine geneq2_wtsi_well(i) 
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

      neqp1=neq+1
c define diagonal term for connectivity
      jmi=nelmdg(i-icd)
c define diagonal term for jacobian matrix 
      jmia=jmi-neqp1
      
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
c this loop contains upper and lower diagonal connections      
      iq=0
      do jm=ii1,ii2
	 kb=nelm(jm)+icd
	 if(kb.ne.i) then
        iq=iq+1
        it8(iq)=kb
        it9(iq)=jm-ii1+1
        it10(iq)=istrw(jm-neqp1)
        it11(iq)=jm-neqp1
c find connection to/from primary grid  
c it12(iq) = 0 will mean non primary grid correction 
c need to find geometric coefficient from primary grid to wellbore
        it12(iq) = 0     
        if(kb.lt.i) then
         ij1 = nelmdg(kb)+1
         ij2 = nelm(kb+1)
         do ij = ij1,ij2
          if(nelm(ij).eq.i) then
           it12(iq) = ij-neqp1
           it10(iq) = istrw(ij-neqp1)
          endif
         enddo   
        endif
	 endif
      enddo
c
c 3-d geometry,added well node(iriver=2)
c
      if(icnl.eq.0) then
         do  jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iwd=it10(jm)
            iw =abs(iwd)
c pny is along the well,pnx is radial (often = 0)            
            axkb=pnx(kb)
            alykb=pny(kb)          
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)            
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            if(kb.gt.neq_primary) then
c along wellbore            
             pxy = sx2c*perml(2)
            else
c connection to primary grid (radial)          
             pxy = sx2c*alxi
            endif
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            t3(neighc)=pxy
            t4(neighc)=pxy                   
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
        enddo
c     
c     liquid phase calculations       
c
         do 70 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            pxyi=t1(neighc)
            sx4d=t6(neighc)
            axyd=pxyi+0.5*sx4d
     &           *(rolf(i)+rolf(kb))
     &           *(cord(kz,igrav)-cord(iz,igrav))        
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

               bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
               a(jmia+nmat(1))=a(jmia+nmat(1))+dlapi
               a(iau+nmat(1))=a(iau+nmat(1))+dlapkb
c
c connection to primary grid
c 
      if(kb.lt.neq_primary) then
       a_axy(ial+nmatavw)=-axy
       bp(kz+nrhs(1))=bp(kz+nrhs(1))-axy  
       a(ial+nmat(1))=a(ial+nmat(1))-dlapi
       a(jml+nmat(1))=a(jml+nmat(1))-dlapkb
                        
      endif              

 72         continue
	endif
	endif
      return
      end
