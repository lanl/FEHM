      subroutine geneq_stress_uncoupled_2D(i)
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

C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To generate equations for stress at each node.
CD1 This simple version does not use symmetry
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 9-08-05     G. Zyvoloski   00022   Initial implementation.
CD2
CD3
CD3
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
C developer notes
C 3-21-06 gaz
C array stencil symmetry not used yet- will put in later
C need to think on how upper and lower diagonal are used 
C*******************************

***************************************
c
c generate equations for 3-d stress with finite elements ,
c full derivatives
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
      use comdti
      use comai
      use comsi
      implicit none

      integer i
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
      integer kb
      integer neighc
      integer iau
      integer ial
      integer kz
      
      integer nmatavw
      real*8 sx1d
      real*8 pvii
      real*8 phii
      real*8 swi
      real*8 dlpi
      real*8 dlpkb
      real*8 dvpi
      real*8 dvpkb
      real*8 dlei
      real*8 dvei
      real*8 dlekb
      real*8 dvekb
      real*8 dlapi
      real*8 dvapi
      real*8 dlapkb
      real*8 dvapkb
      real*8 dili
      real*8 dilkb
      real*8 divi
      real*8 divkb
      real*8 axkb
      real*8 aykb
      real*8 azkb
      real*8 alxkb
      real*8 alykb
      real*8 alzkb
      real*8 sx2c
      real*8 sx4d

      real*8 pvikb
      real*8 phikb
      real*8 radi
      real*8 radkb
      real*8 fid
      real*8 fid1
      real*8 axyd
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
      real*8 dis2,dis_tol,sx_min
      real*8 delx2
      real*8 dely2
      real*8 delz2
      real*8 redu
      real*8 divpkbction_factor
      real*8 grav_air
      real*8 grav_term

      real*8 e1i,e2i,e3i,e1kb,e2kb,e3kb
      real*8 e1bar,e2bar,e3bar
	real*8 efac,efaci, efackb, bulki,bulkkb,bulkb,alpi,alpkb,alphab
	real*8 ealphai, ealphakb, ebulki, ebulkkb, ealphabar, ebulkbar
      real*8 dui,dvi,dwi,dukb,dvkb,dwkb
      real*8 xxx,xyy,xzz,xyx,xxy,xzx,yxy
      real*8 yyx,yyy,yxx,yzz,yzy,yyz,zxz
      real*8 zzx,zyz,zzy,zzz,zyy,zxx,xxz
      real*8 xx,xy,xz,yy,yx,yz,zx,zy,zz
      real*8 ddx,ddy,ddz,pxx,pxz,pyx
	real*8 pyy,pyz,pzx,pzy,pzz
	real*8 fad,fbd,fcd
      real*8 pdumi,tdumi,pdumkb,tdumkb,pdumt,tdumt,roli
	real*8 tdumx,tdumy,tdumz
	real*8 bforcex,bforcey,bforcez
	real*8 dtdumxp,dtdumxt
	real*8 dtdumyp,dtdumyt
	real*8 dtdumzp,dtdumzt
      real*8 sixsjx,siysjy,sizsjz,sixsjy
      real*8 siysjx,sixsjz,sizsjx,siysjz
      real*8 sizsjy,sjsix,sjsiy,sjsiz
c      integer, allocatable ::   itstress(:)
      integer iws

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    

  


c changed by avw -- entered here by seh
      neqp1=neq+1
      ldna=nelm(neqp1)-neqp1
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
c zero out temporary storage
c
c      if(.not.allocated(itstress)) then
c       allocate(itstress(200))
c      endif

c
c fluid and grid properties at node i
c
      sx1d=sx1(i)
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
      e1i = e1(i)
      e2i = e2(i)
      e3i = e3(i)
c recall displacements at the node i
      dui = du(i)
      dvi = dv(i)

c pore pressure and thermal expansion terms
      bulki=bulk(i)
      alpi=alp(i)
      tdumi=t(i)-tini(i)
      pdumi=phi(i) - phini(i)
c reference rock density
      roli=denr(i)
c      
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
c
      neqp1=neq+1
c define diagonal term for connectivity
      jmi=nelmdg(i-icd)
c define diagonal term for jacobian matrix 
      jmia=jmi-neqp1
      
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      
      iq=0
      do jm=ii1,ii2
	 kb=nelm(jm)+icd
	 if(kb.ne.i) then
        iq=iq+1
        it8(iq)=kb
        it9(iq)=jm-ii1+1
        it10(iq)=istrw(jm-neqp1)
        itstress(iq) = istrws(jm-neqp1) 
        it11(iq)=jm-neqp1
	 endif
      enddo
c
c body forces
c
c  gravity term (grav) is times 1.e-6 to account for stress in Mpa
c  also grav is negative so this should be considered in force balance
c 
   	   bforcex = 0.0d0
	   bforcey = 0.0d0
         if(ibodyforce.ne.0) then
          if(igrav.eq.1) then
           bforcex = -roli*grav*sx1d
	    else if(igrav.eq.2) then
           bforcey = -roli*grav*sx1d
	    else if(igrav.eq.3) then
           bforcey = -roli*grav*sx1d
	    endif
	   endif
	   bp(iz+nrhs(1)) = bforcex
	   bp(iz+nrhs(2)) = bforcey
c
c 2-d geometry 
c
         do jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iw = it10(jm)
            iws=itstress(jm)    
            iau=it11(jm)	     
          
c     
c recall shape function integrals, calculated in gencof
c
c   xx term
         sixsjx=sxs(iws,5)
c   yy term
         siysjy=sxs(iws,6)
c   xy term
         sixsjy=sxs(iws,1) 
c   yx term
         siysjx=sxs(iws,2)
                            
            e1kb = e1(kb)
            e2kb = e2(kb)
            e3kb = e3(kb)
            e1bar=2.*e1i*e1kb/(e1i+e1kb + dis_tol)
            e2bar=2.*e2i*e2kb/(e2i+e2kb + dis_tol)
            e3bar=2.*e3i*e3kb/(e3i+e3kb + dis_tol)
            
            dukb=du(kb)
            dvkb=dv(kb)
            
c
c compute the terms in the stiffness matrix  k(i,j)
c
         xxx = e1bar*sixsjx
         xyy = e3bar*siysjy
         
         xyx=e3bar*siysjx
         xxy=e2bar*sixsjy
         
         
         yxy=e3bar*sixsjy
         yyx=e2bar*siysjx
         yyy=e1bar*siysjy
         yxx=e3bar*sixsjx
         
         
         
         
c
c form stiffness matrix
c
         xx=xxx+xyy
         xy=xyx+xxy
         
         yx=yxy+yyx
         yy=yyy+yxx
         
         
         ddx=(dukb-dui)
         ddy=(dvkb-dvi)
        

         pxx=xx*ddx
         pxy=xy*ddy
         
         pyx=yx*ddx
         pyy=yy*ddy
       
       

c         
c form residuals of the stress balance equations
c
         fad=pxx+pxy
         fbd=pyx+pyy
         
c         
c form residuals of the stress balance equations
c
	      
         bp(iz+nrhs(1))=bp(iz+nrhs(1))+fad
         bp(iz+nrhs(2))=bp(iz+nrhs(2))+fbd 



c x equation derivatives

      a(jmia+nmat(1))=a(jmia+nmat(1))-xx
      a(jmia+nmat(2))=a(jmia+nmat(2))-xy
  
     
      a(iau+nmat(1))=a(iau+nmat(1))+xx   
      a(iau+nmat(2))=a(iau+nmat(2))+xy
     
                
      
c y equation derivatives

      a(jmia+nmat(3))=a(jmia+nmat(3))-yx
      a(jmia+nmat(4))=a(jmia+nmat(4))-yy
        
      a(iau+nmat(3))=a(iau+nmat(3))+yx
      a(iau+nmat(4))=a(iau+nmat(4))+yy
                  
      enddo
c
c add thermal expansion and biot terms
c

	do jm = ii1,ii2
      kb = nelm(jm)
	iws = istrws(jm-neqp1)

c x term for pore pressure and thermal expansion term
         sjsix=sxs(iws,7)
c y term for pore pressure and thermal expansion term
         sjsiy=sxs(iws,8)
    
            e1kb = e1(kb)
            e2kb = e2(kb)
            e3kb = e3(kb)
            e1bar=2.*e1i*e1kb/(e1i+e1kb + dis_tol)
            e2bar=2.*e2i*e2kb/(e2i+e2kb + dis_tol)
            e3bar=2.*e3i*e3kb/(e3i+e3kb + dis_tol)
            alpkb=alp(kb)
            alphab=2.*alpi*alpkb/(alpi+alpkb + dis_tol)
c boit term
            bulkkb=bulk(kb)
            bulkb=2.*bulkkb*bulki/(bulkkb+bulki + dis_tol)

c            efac = 3.d0*e2bar + 2.d0*e3bar
! Sai 10/08/2010 : Modified averaging which seems to produce analytically correct
! solutions for simple geometries
            efaci = 3.0d0*e2i + 2.0d0*e3i
            ealphai = efaci*alpi
            ebulki  = efaci*bulki

            efackb = 3.0d0*e2kb + 2.0d0*e3kb
            ealphakb = efackb*alpkb
            ebulkkb  = efackb*bulkkb

            ealphabar = 2.0d0*ealphai*ealphakb
            ealphabar = ealphabar/(ealphai + ealphakb + dis_tol)
            ebulkbar = 2.0d0*ebulki*ebulkkb
            ebulkbar = ebulkbar/(ebulki + ebulkkb + dis_tol)

            ealphabar = 2.0d0*ealphabar - ealphai
            ebulkbar  = 2.0d0*ebulkbar  - ebulki

            if(istrs_coupl.eq.-99) then
             tdumt=t(kb)-tini(kb)
             pdumt=phi(kb)-phini(kb)
	      else
             tdumt=t(kb)-tini(kb)
             pdumt=phi(kb)-phini(kb)
	      endif
c

! Sai 10/08/2010 
c           tdumx=sjsix*(tdumt*alphab+pdumt*bulkb)*efac
c           tdumy=sjsiy*(tdumt*alphab+pdumt*bulkb)*efac
           tdumx=sjsix*(tdumt*ealphabar+pdumt*ebulkbar)
           tdumy=sjsiy*(tdumt*ealphabar+pdumt*ebulkbar)
          
c         
c form residuals of the stress balance equations
c
	      
         bp(iz+nrhs(1))=bp(iz+nrhs(1))+tdumx
         bp(iz+nrhs(2))=bp(iz+nrhs(2))+tdumy 

      enddo
      
      return
      end
      subroutine stress_boun2(iflg,iz)
c
c apply bounday conditions for stress equations
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
      use comdti
      use comai
      use comsi

      implicit none
      integer i1,i2,jj,kb,jjkb
      integer i,iflg,iz,neqp1,jji,krdu,krdv,krdw
	real*8 fac, delu, delv
	real*8 du_boun, dv_boun
      real*8 forcx_boun, forcy_boun, xboun, yboun
	
      if(iflg.eq.0) then

	else if(iflg.eq.1) then
	 neqp1 = neq+1
	 do i = 1,neq
c check for x fixed displacement	     
           i1 = nelm(i)+1
           i2 = nelm(i+1)

           do jj= i1,i2
            kb = nelm(jj)
            if(kr(kb,1).eq.1) then
             jjkb = jj-neqp1
             delu = disp(kb,1) - du(kb)
             bp(i+nrhs(1)) = bp(i+nrhs(1)) - a(jjkb+nmat(1))*delu
             bp(i+nrhs(2)) = bp(i+nrhs(2)) - a(jjkb+nmat(3))*delu
             a(jjkb+nmat(1)) = 0.0d0
             a(jjkb+nmat(3)) = 0.0d0
            endif
           enddo
          jji = nelmdg(i)-neqp1

          krdu = kr(i,1)
	    if(krdu.eq.1) then
            do jj= i1,i2
             jjkb = jj-neqp1
             a(jjkb+nmat(1)) = 0.0d0 
             a(jjkb+nmat(2)) = 0.0d0               
            enddo
            jji = nelmdg(i)-neqp1
            a(jji+nmat(1)) = 1.
!		  bp(i+nrhs(1)) = 0.0d0		     
                  bp(i+nrhs(1)) = (du(i) - disp(i,1))
		else if(krdu.eq.-1) then
c    added force
	     forcx_boun = forc(i,1)
	     bp(i+nrhs(1)) = bp(i+nrhs(1)) - forcx_boun
		endif 
    
c check for y fixed displacement	 
           i1 = nelm(i)+1
           i2 = nelm(i+1)
           do jj= i1,i2
            kb = nelm(jj)
            if(kr(kb,2).eq.2) then
             jjkb = jj-neqp1
             delv = disp(kb,2) - dv(kb)
             bp(i+nrhs(1)) = bp(i+nrhs(1)) - a(jjkb+nmat(2))*delv
             bp(i+nrhs(2)) = bp(i+nrhs(2)) - a(jjkb+nmat(4))*delv
             a(jjkb+nmat(2)) = 0.0d0
             a(jjkb+nmat(4)) = 0.0d0
         
            endif
           enddo
          jji = nelmdg(i)-neqp1
	    
          krdv = kr(i,2)
	    if(krdv.eq.2) then
            do jj= i1,i2
             jjkb = jj-neqp1
             a(jjkb+nmat(3)) = 0.0d0 
             a(jjkb+nmat(4)) = 0.0d0 
                   
            enddo
           jji = nelmdg(i)-neqp1
           a(jji+nmat(4)) = 1.	
!		 bp(i+nrhs(2)) = 0.0d0		       		     
                 bp(i+nrhs(2)) = (dv(i) - disp(i,2))
		else if(krdv.eq.-2) then
c    added force
	     forcy_boun = forc(i,2)
	     bp(i+nrhs(2)) = bp(i+nrhs(2)) - forcy_boun
		endif  
       enddo
      else if(iflg.eq.2) then

	endif
  
	return
	end

