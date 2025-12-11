      subroutine stress_3D_post(i)
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
      use comfem
      implicit none

      integer i, li
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
	real*8 efac,bulki,bulkkb,bulkb,alpi,alpkb,alphab
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
	real*8 sisjx,sisjy,sisjz,termx,termy,termz
      real*8 vol,sxd,syd,tauxy,tauyz,tauzx,ctherm,biot,erat
	real*8 epi,eti,dpd,dt,shti,shpi

c s kelkar 12/11/09 axisymmetric anisotropy
	real*8  ezzi,ezzkb,ezzbar,efacxy,efacz 
	real*8 shtixy,shtiz,shpixy,shpiz, e4i,e4kb,e4bar
        real*8 sheari, shearkb,shearbar,efac_ks,betat,efac_betat  
c..................................................................
c gaz 052017        
c      real*8 eigenvec(3,3),alambda(3)
      real*8 onedV, fac
      integer j

c      integer, allocatable ::   itstress(:)
      integer iws

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    
c.......................................
      integer k
      integer    :: ContrElem, el
      integer :: node
      real*8 val
      real*8 :: Dmat(6,6), sol(8)
      real*8 :: rhs1(8), rhs2(8), rhs3(8), rhs4(8), rhs5(8), rhs6(8)
      real*8, allocatable :: kmat(:)
      real*8 deltaT, deltaP, alphadT, betadP
      real*8 kmat_gp
      real*8 xdV, ydV, zdV

      if(ifem.eq.1) then
        allocate(kmat(neq))
        kmat = 0.d00

        ContrElem = 0
        str_x(i) = 0.0d0
        str_y(i) = 0.0d0
        str_z(i) = 0.0d0
        str_xy(i) = 0.0d0
        str_yz(i) = 0.0d0
        str_xz(i) = 0.0d0

        strain_xx(i) = 0.0d0
        strain_yy(i) = 0.0d0
        strain_zz(i) = 0.0d0
        strain_xy(i) = 0.0d0
        strain_yz(i) = 0.0d0
        strain_zx(i) = 0.0d0

        onedV = 0.0d0

        do el=1,nei
          do k=1,ns 
            if(elnode(el,k).eq.i) then
              ContrElem = ContrElem + 1

              do j=1,numgausspoints
                fac = detJ(el,j)*gpweight(j)

                e1bar = 0.0d0
                e2bar = 0.0d0
                e3bar = 0.0d0

                ! Get alpha*deltaT and beta*deltaP at gausspoint
                alphadT = 0.0d0
                betadP  = 0.0d0
                do li=1,ns
                  node = elnode(el, li)
                  deltaT = t(node) - tini(node)
                  deltaP = phi(node) - phini(node)

                  alphadT = alphadT + 
     &                      Psi(el, j, li)*alp(node)*deltaT
                  betadP  = betadP  +
     &                      Psi(el, j, li)*bulk(node)*deltaP
                  e1bar = e1bar + Psi(el, j, li)*e1(node)
                  e2bar = e2bar + Psi(el, j, li)*e2(node)
                  e3bar = e3bar + Psi(el, j, li)*e3(node)

                enddo

                if(iPlastic.eq.1) then
                  call fem_material_stiffness(el, j, Dmat)
                else
                  Dmat = 0.0d0
                  Dmat(1,1) = e1bar
                  Dmat(1,2) = e2bar
                  Dmat(1,3) = e2bar
                  Dmat(2,1) = e2bar
                  Dmat(2,2) = e1bar
                  Dmat(2,3) = e2bar
                  Dmat(3,1) = e2bar
                  Dmat(3,2) = e2bar
                  Dmat(3,3) = e1bar
                  Dmat(4,4) = e3bar
                  Dmat(5,5) = e3bar
                  Dmat(6,6) = e3bar
                endif
                  
                kmat_gp = (Dmat(1,1) + Dmat(1,2) + Dmat(1,3))
                kmat(i) = kmat(i) + kmat_gp*fac
                onedV = onedV + 1.0d0*fac

                rhs1(j) = fem_stress(el, j, 1) + 
     &                    kmat_gp*(alphadT + betadP)
                rhs2(j) = fem_stress(el, j, 2) +
     &                    kmat_gp*(alphadT + betadP)
                rhs3(j) = fem_stress(el, j, 3) +
     &                    kmat_gp*(alphadT + betadP)
                rhs4(j) = fem_stress(el, j, 4)
                rhs5(j) = fem_stress(el, j, 5)
                rhs6(j) = fem_stress(el, j, 6)
              enddo 

              sol = matmul(iPsi, rhs1)
              str_x(i) = str_x(i) + sol(k)
              sol = matmul(iPsi, rhs2)
              str_y(i) = str_y(i) + sol(k)
              sol = matmul(iPsi, rhs3)
              str_z(i) = str_z(i) + sol(k)
              sol = matmul(iPsi, rhs4)
              str_xy(i) = str_xy(i) + sol(k)
              sol = matmul(iPsi, rhs5)
              str_yz(i) = str_yz(i) + sol(k)
              sol = matmul(iPsi, rhs6)
              str_xz(i) = str_xz(i) + sol(k)

            endif
          enddo
        enddo
        
        kmat(i) = kmat(i)/onedV
        
        tdumt = alp(i)*(t(i) - tini(i))
        pdumt = bulk(i)*(phi(i) - phini(i))

        str_x(i) = str_x(i)/ContrElem - kmat(i)*(tdumt + pdumt)
        str_y(i) = str_y(i)/ContrElem - kmat(i)*(tdumt + pdumt)
        str_z(i) = str_z(i)/ContrElem - kmat(i)*(tdumt + pdumt)
        str_xy(i) = str_xy(i)/ContrElem
        str_yz(i) = str_yz(i)/ContrElem
        str_xz(i) = str_xz(i)/ContrElem

!        str_x(i) = str_x(i)/onedV
!        str_y(i) = str_y(i)/onedV
!        str_z(i) = str_z(i)/onedV
!        str_xy(i) = str_xy(i)/onedV
!        str_yz(i) = str_yz(i)/onedV
!        str_xz(i) = str_xz(i)/onedV

!        strain_xx(i) = strain_xx(i)/onedV
!        strain_yy(i) = strain_yy(i)/onedV
!        strain_zz(i) = strain_zz(i)/onedV
!        strain_xy(i) = strain_xy(i)/onedV
!        strain_yz(i) = strain_yz(i)/onedV
!        strain_zx(i) = strain_zx(i)/onedV

        deallocate(kmat)

c     s kelkar change signs of diagonal stress components to be 
c     compatible with standard FEHM
        str_x(i) = -str_x(i)
        str_y(i) = -str_y(i)
        str_z(i) = -str_z(i)

      else
c......................................
c     changed by avw -- entered here by seh
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
c     zero out temporary storage
c     
c     if(.not.allocated(itstress)) then
c     allocate(itstress(200))
c     endif
         
c     
c     fluid and grid properties at node i
c     
         sx1d=sx1(i)
         vol = sx1d
         pvii=phi(i)     
         phii=pvii
         ti=t(i)
         swi=s(i)
c     recall displacements at the node i
         dui = du(i)
         dvi = dv(i)
         dwi = dw(i)
         
c     reference rock density
         roli=denr(i)
         bforcex = 0.0d0
         bforcey = 0.0d0
         bforcez = 0.0d0
c     
c     
c     form constants for i>neq
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
         
         neqp1=neq+1
c     define diagonal term for connectivity
         jmi=nelmdg(i-icd)
c     define diagonal term for jacobian matrix 
         jmia=jmi-neqp1
         
         ii1=nelm(i-icd)+1
         ii2=nelm(i-icd+1)
         idg=nelmdg(i-icd)-ii1+1
         
         iq=0
         do jm=ii1,ii2
            kb=nelm(jm)+icd
            
            iq=iq+1
            it8(iq)=kb
            itstress(iq) = istrws(jm-neqp1) 
            
         enddo
         
         termx=0.
         termy=0.
         termz=0.
         tauxy=0.
         tauyz=0.
         tauzx=0.
         
         
         ctherm=alp(i)
c     ctherm=coef of temp expansion,  biot=1./3h from biot's paper
         biot=bulk(i)
         e1i = e1(i)
         e2i = e2(i)
         e3i = e3(i)
         erat=e2i/e1i
         dpd=phi(i)-phini(i)
         dt=t(i)-tini(i)
c     s kelkar 12/11/09 axisymmetric anisotropy.........................
         if(stress_anisotropy_use) then
            ezzi= ezz(i)
            e4i=e4(i)
            sheari=shearmod_t(i)
            efacxy= (e1i+e2i+e4i)
            efacz = 2.0d0*e4i + ezzi
c     pore pressure terms. See Keita's notes dated 2/25/2010, Here 
c     efacxy is defined as efacxy=3*Hp=(C11+C12+C13) and  
c     bulk()=beta_p/3Hp or  beta_p=bulk*efacxy. Then
c     Ks=Hp/(1-beta_p)=(efacxy/3)/(1-efacxy*bulkb). then beta_t 
c     calculated from
c     betat=1-Ht/Ks=1-(efacz/3)/Ks where efacz=3*Ht=(2C13+C33)
c     also define efac_betat=betat/efacz for consistancy of notation
            efac_ks=(efacxy/3.0)/(1.0-efacxy*biot)
            betat=1.0-(efacz/3.0)/efac_ks
            efac_betat=1.0/efacz-1.0/(3.0*efac_ks)
            shtixy=efacxy*ctherm*tdumt
            shpixy=efacxy*biot*dpd
            shtiz=efacz*ctherm*tdumt
            shpiz=efacz*efac_betat*dpd
            str_x(i)=(shpixy+shtixy)
            str_y(i)=(shpixy+shtixy)
            str_z(i)=(shpiz+shtiz)
c.......................................................
         else
            efac=3.d0*e2i+2.d0*e3i
c     stress due to temp and pore pressure changes
            epi=efac*biot
            eti=efac*ctherm
c     changed 5-25-08 gaz
            shti=(eti*dt)
            shpi=(epi*dpd)
            str_x(i)=(shpi+shti)
            str_y(i)=(shpi+shti)
            str_z(i)=(shpi+shti)
         endif
         
c     
c     3-d geometry 
c     
         do jm=1,iq
            kb=it8(jm)
            kz=kb-icd          
            
            iws=itstress(jm)    
            
            
c     recall geometric integrals, calculated in gencof.
c     might be 7 8 9 gaz 11-2-2006
            sisjx=sxs(iws,7)
            sisjy=sxs(iws,8)
            sisjz=-sxs(iws,9)
c     average material properties                     
            e1kb = e1(kb)
            e2kb = e2(kb)
            e3kb = e3(kb)
            e1bar=2.*e1i*e1kb/(e1i+e1kb + dis_tol)
            e2bar=2.*e2i*e2kb/(e2i+e2kb + dis_tol)
            e3bar=2.*e3i*e3kb/(e3i+e3kb + dis_tol)
c     s kelkar 12/5/09 axisymmetric anisotropy.........................
            if(stress_anisotropy_use) then
               ezzkb = ezz(kb)
               e4kb = e4(kb)
               shearkb=shearmod_t(kb)
               ezzbar= 2.*ezzi*ezzkb/(ezzi+ezzkb + dis_tol)
               e4bar= 2.*e4i*e4kb/(e4i+e4kb + dis_tol)
               shearbar= 2.*sheari*shearkb
               shearbar=shearbar/(sheari+shearkb+dis_tol)
            endif
c     thermal conductivity
c     alpkb=alp(kb)
c     alphab=2.*alpi*alpkb/(alpi+alpkb + dis_tol)
c     biot term
c     bulkkb=bulk(kb)
c     bulkb=2.*bulkkb*bulki/(bulkkb+bulki + dis_tol)
c     efac = 3.d0*e2bar + 2.d0*e3bar
            
c     tdumkb=t(kb)-tini(kb)
c     pdumkb=phi(kb)-phini(kb)
c     tdumt=0.5*(tdumkb+tdumi)
c     pdumt=0.5*(pdumkb+pdumi)
c     
c     tdumx=sjsix*(tdumt*alphab+pdumt*bulkb)*efac
c     tdumy=sjsiy*(tdumt*alphab+pdumt*bulkb)*efac
c     tdumz=sjsiz*(tdumt*alphab+pdumt*bulkb)*efac
            
c     calculate stresses 
            
            ddx=du(kb)-dui
            ddy=dv(kb)-dvi
            ddz=dw(kb)-dwi
c     s kelkar 12/5/09 axisymmetric anisotropy.........................
c     e1=c11, e2=c12=c21, e3=c66, e4=c13, and ezz=c33
c     these goto isotropic limit when Ep=Et and Nue-p=Nue-t
            if(stress_anisotropy_use) then
               xx=e1bar*sisjx*ddx
               xy=e2bar*sisjy*ddy
               xz=e4bar*sisjz*ddz
               yx=e2bar*sisjx*ddx
               yy=e1bar*sisjy*ddy
               yz=e4bar*sisjz*ddz
               zx=e4bar*sisjx*ddx
               zy=e4bar*sisjy*ddy
               zz=ezzbar*sisjz*ddz
               xyx=e3bar*sisjy*ddx
               xyy=e3bar*sisjx*ddy
               yzy=shearbar*sisjz*ddy
               yzz=shearbar*sisjy*ddz
               zxx=shearbar*sisjz*ddx
               zxz=shearbar*sisjx*ddz
               
            else
               xx=e1bar*sisjx*ddx
               xy=e2bar*sisjy*ddy
               xz=e2bar*sisjz*ddz
               yx=e2bar*sisjx*ddx
               yy=e1bar*sisjy*ddy
               yz=e2bar*sisjz*ddz
               zx=e2bar*sisjx*ddx
               zy=e2bar*sisjy*ddy
               zz=e1bar*sisjz*ddz
               xyx=e3bar*sisjy*ddx
               xyy=e3bar*sisjx*ddy
               yzy=e3bar*sisjz*ddy
               yzz=e3bar*sisjy*ddz
               zxx=e3bar*sisjz*ddx
               zxz=e3bar*sisjx*ddz
            endif
c     
            termx=termx+(xx+xy+xz)
            termy=termy+(yx+yy+yz)
            termz=termz+(zx+zy+zz)
            tauxy=tauxy+xyx+xyy
            tauyz=tauyz+yzy+yzz
            tauzx=tauzx+zxx+zxz 
            
         enddo
         str_x(i)=(str_x(i)+termx/vol)
         str_y(i)=(str_y(i)+termy/vol)
         str_z(i)=(str_z(i)+termz/vol)
         str_xy(i)=tauxy/vol
         str_yz(i)=tauyz/vol
         str_xz(i)=tauzx/vol
         
       if(residual_stress) then
        str_x(i)=str_x(i)+str_x0(i)
        str_y(i)=str_y(i)+str_y0(i)
        str_z(i)=str_z(i)+str_z0(i)
        str_xy(i)=str_xy(i)+str_xy0(i)
        str_yz(i)=str_yz(i)+str_yz0(i)
        str_xz(i)=str_xz(i)+str_xz0(i)
      endif        
         
         if(flag_principal.eq.1) then
c gaz 052017             
c           call principal_stress_3D(i,alambda,eigenvec)
            call principal_stress_3D(i)
c     save the eigenvlaues in str_z,str)y,str_x in decreasing order
c     str_z is the max principal stress
            str_x(i)= alambda(1)
            str_y(i)= alambda(2)
            str_z(i)= alambda(3)
c     save two direction cosins of max principal direction 
c     in str_xz,str_yz
            str_xz(i) = eigenvec(1,3)
            str_yz(i) = eigenvec(2,3)
c     save the first dir cos for the min principal direction in str_x
c     the rest of the dir cos values can be obtained from 
c     ortho-normality
            str_xy(i) = eigenvec(1,1)
         endif
         
      endif
      
      return
      end
      
c................................................................

c gaz 052017      
c      subroutine principal_stress_3D(i,alambda,eigenvec)
      subroutine principal_stress_3D(i)

      use comai
      use comsi

      implicit none
      integer i
c gaz 052017      
c      real*8 AMAT(3,3), eigenvec(3,3),alambda(3)
      real*8 AMAT(3,3)
      real*8 AI1,AI2,AI3

      AMAT=0.0
      AMAT(1,1) = str_x(i)
      AMAT(2,2) = str_y(i)
      AMAT(3,3) = str_z(i)
      AMAT(1,2) = str_xy(i)
      AMAT(1,3) = str_xz(i)
      AMAT(2,3) = str_yz(i)
      AMAT(2,1) = AMAT(1,2)
      AMAT(3,2) = AMAT(2,3)
      AMAT(1,3) = AMAT(3,1)
      
! Compute the invariants for the matrix
      AI1 = AMAT(1,1) + AMAT(2,2) + AMAT(3,3)
      AI2 = AMAT(1,1)*AMAT(2,2) + AMAT(2,2)*AMAT(3,3) 
      AI2 = AI2 + AMAT(3,3)*AMAT(1,1) - AMAT(1,2)*AMAT(1,2) 
     &     - AMAT(2,3)*AMAT(2,3)  - AMAT(3,1)*AMAT(3,1)
      AI3 = AMAT(1,1)*(AMAT(2,2)*AMAT(3,3) - AMAT(2,3)*AMAT(3,2))
      AI3 = AI3 + AMAT(1,2)*(AMAT(2,3)*AMAT(3,1) - AMAT(2,1)*AMAT(3,3))
      AI3 = AI3 + AMAT(1,3)*(AMAT(2,1)*AMAT(3,2) - AMAT(3,1)*AMAT(2,2))
      
! Solve the characteristic equation for the eigenvalues
      call solve_cubic(-AI1, AI2,-AI3,alambda)
      
! Make the matrix tridiagonal using Householder transformation
      call reduce_to_tridiagonal(AMAT)
      
! Solve for the eigenvectors
      call eigenvectors(AMAT, alambda, eigenvec)
      
      return

      end

c.......................................................................
