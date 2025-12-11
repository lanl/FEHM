      subroutine geneq_plastic_3D()
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
! Assembles the equations for plastic deformations using FEHM algorithm
! 
! Author : Sai Rapaka
!

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
      integer jmj, jmja, ijdg, jjm, kbi, jau
      
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
      real*8 efac,bulki,bulkkb,bulkb,alpi,alpkb,alphab
      real*8 dui,dvi,dwi,dukb,dvkb,dwkb
      real*8 xxx,xyy,xzz,xyx,xxy,xzx,yxy
      real*8 yyx,yyy,yxx,yzz,yzy,yyz,zxz
      real*8 zzx,zyz,zzy,zzz,zyy,zxx,xxz
      real*8 xx,xy,xz,yy,yx,yz,zx,zy,zz
      real*8 ddx,ddy,ddz,pxx,pxz,pyx
      real*8 pyy,pyz,pzx,pzy,pzz
      real*8 fad,fbd,fcd
      real*8 fad2,fbd2,fcd2
      real*8 pdumi,tdumi,pdumkb,tdumkb,pdumt,tdumt,roli
      real*8 tdumx,tdumy,tdumz
      real*8 bforcex,bforcey,bforcez
      real*8 dtdumxp,dtdumxt
      real*8 dtdumyp,dtdumyt
      real*8 dtdumzp,dtdumzt
      real*8 sixsjx,siysjy,sizsjz,sixsjy
      real*8 siysjx,sixsjz,sizsjx,siysjz
      real*8 sizsjy,sjsix,sjsiy,sjsiz
      real*8 sixsj, siysj, sizsj
      real*8 sxx_bar, syy_bar, szz_bar, sxy_bar, syz_bar, szx_bar
      real*8 harmonic_mean
      
      real*8, dimension(1:3,1:3) :: nodalk
      real*8, dimension(1:6,1:6) :: d_i, d_j
      real*8 avg_d

      integer iws

      parameter(dis_tol=1.d-12)

      logical bit
      integer isl
      integer iz4m1
      integer imd,iwd    
      integer mm,qq,p
  

c code below needed when dpdp is implemented with stress
c changed by avw -- entered here by seh
      ! neqp1=neq+1
      ! ldna=nelm(neqp1)-neqp1
      ! if(i.gt.neq) then
      !    nmatavw=ldna
      ! else
      !    nmatavw=0
      ! endif
      if(icons.le.abs(maxit)) then
       grav_air=0.0
      else
       grav_air=grav
      endif

!      call elastic_stress_update()

      call material_stress_update()
            
      do i=1,neq

        !fluid and grid properties at node i
        sx1d=sx1(i)
       
        !pore pressure and thermal expansion terms
        bulki=bulk(i)
        alpi=alp(i)
        
        !reference rock density
        roli=denr(i)

        !form constants for i>neq
        !needed for dpdp - however, loop will also have to be changed
        if(i.gt.neq.and.idualp.eq.0) then
          icd=neq
        else
          icd=0
        endif
        iz=i-icd

        iz4m1 = 4*(iz-1)+1

        neqp1=neq+1
        !define diagonal term for connectivity
        jmi=nelmdg(i-icd)
        !define diagonal term for jacobian matrix 
        jmia=jmi-neqp1
      
        ii1=nelm(i-icd)+1
        ii2=nelm(i-icd+1)
        idg=nelmdg(i-icd)-ii1+1
      
        !body forces
        !
        !gravity term (grav) is times 1.e-6 to account for stress in Mpa
        !also grav is negative so this should be considered in force balance 
        bforcex = 0.0d0
        bforcey = 0.0d0
        bforcez = 0.0d0
        if(ibodyforce.ne.0) then
          if(igrav.eq.1) then
            bforcex = -roli*grav*sx1d
          else if(igrav.eq.2) then
            bforcey = -roli*grav*sx1d
          else if(igrav.eq.3) then
            bforcez = -roli*grav*sx1d
          endif
        endif
        bp(iz+nrhs(1)) = bforcex
        bp(iz+nrhs(2)) = bforcey
        bp(iz+nrhs(3)) = bforcez

        call material_stiffness(i, d_i)
        
       !Start loop over the connected neighbors
        do jm=ii1,ii2
          kb=nelm(jm)+icd
          kz=kb-icd
          iws = istrws(jm-neqp1)
          iau = jm - neqp1
          
          jmj = nelmdg(kb - icd)
          jmja = jmj - neqp1

          ij1=nelm(kb-icd)+1
          ij2=nelm(kb-icd+1)
          ijdg=nelmdg(kb-icd)-ij1+1
      
          do jjm=ij1,ij2
            kbi=nelm(jjm)+icd
            if(kbi.eq.i) then
              jau=jjm-neqp1
            endif
          enddo

         ! Computation of F_int = B^{T}\sigma
         ! F_int_x = dNidxNj*sigma_xx(j) + dNidyNj*sigma_xy(j) + dNidzNj*sigma_xz(j)
          sixsj =  sxs(iws,13)
          siysj =  sxs(iws,14)
          sizsj = -sxs(iws,15)

          ! sxx_bar = harmonic_mean(str_x(i), str_x(kb))
          ! syy_bar = harmonic_mean(str_y(i), str_y(kb))
          ! szz_bar = harmonic_mean(str_z(i), str_z(kb))
          ! sxy_bar = harmonic_mean(str_xy(i), str_xy(kb))
          ! syz_bar = harmonic_mean(str_yz(i), str_yz(kb))
          ! szx_bar = harmonic_mean(str_xz(i), str_xz(kb))

          sxx_bar = str_x(kb)
          syy_bar = str_y(kb)
          szz_bar = str_z(kb)
          sxy_bar = str_xy(kb)
          syz_bar = str_yz(kb)
          szx_bar = str_xz(kb)

          ! sxx_bar = 0.5d0*(str_x(i) + str_x(kb))
          ! syy_bar = 0.5d0*(str_y(i) + str_y(kb))
          ! szz_bar = 0.5d0*(str_z(i) + str_z(kb))
          ! sxy_bar = 0.5d0*(str_xy(i) + str_xy(kb))
          ! syz_bar = 0.5d0*(str_yz(i) + str_yz(kb))
          ! szx_bar = 0.5d0*(str_xz(i) + str_xz(kb))

          fad2 = sixsj*sxx_bar + siysj*sxy_bar + sizsj*szx_bar
          fbd2 = sixsj*sxy_bar + siysj*syy_bar + sizsj*syz_bar
          fcd2 = sixsj*szx_bar + siysj*syz_bar + sizsj*szz_bar

          bp(iz + nrhs(1)) = bp(iz + nrhs(1)) + fad2
          bp(iz + nrhs(2)) = bp(iz + nrhs(2)) + fbd2
          bp(iz + nrhs(3)) = bp(iz + nrhs(3)) + fcd2

          e1kb = e1(kb)
          e2kb = e2(kb)
          e3kb = e3(kb)
          alpkb = alp(kb)
          bulkkb = bulk(kb)

          e1bar = harmonic_mean(e1i, e1kb)
          e2bar = harmonic_mean(e2i, e2kb)
          e3bar = harmonic_mean(e3i, e3kb)
          alphab = harmonic_mean(alpi, alpkb)
          bulkb = harmonic_mean(bulki, bulkkb)

          efac = 3.d0*e2bar + 2.d0*e3bar

          tdumt=t(kb)-tini(kb)
          pdumt=phi(kb)-phini(kb)
          
          tdumx=sixsj*(tdumt*alphab+pdumt*bulkb)*efac
          tdumy=siysj*(tdumt*alphab+pdumt*bulkb)*efac
          tdumz=sizsj*(tdumt*alphab+pdumt*bulkb)*efac

          !bp(iz+nrhs(1))=bp(iz+nrhs(1))+tdumx
          !bp(iz+nrhs(2))=bp(iz+nrhs(2))+tdumy 
          !bp(iz+nrhs(3))=bp(iz+nrhs(3))+tdumz 

          if(kb.ne.i) then
          
            call material_stiffness(kb, d_j)
           !temporarily specify the material properties matrix manually            
            ! e1kb = e1(kb)
            ! e2kb = e2(kb)
            ! e3kb = e3(kb)

            ! d_j = 0.0d0
            ! d_j(1,1) = e1kb
            ! d_j(1,2) = e2kb
            ! d_j(1,3) = e2kb
            ! d_j(2,1) = e2kb
            ! d_j(2,2) = e1kb
            ! d_j(2,3) = e2kb
            ! d_j(3,1) = e2kb
            ! d_j(3,2) = e2kb
            ! d_j(3,3) = e1kb
            ! d_j(4,4) = e3kb
            ! d_j(5,5) = e3kb
            ! d_j(6,6) = e3kb

           !Assemble the 3x3 stiffness matrix corresponding to nodes i and j
            do mm=1,3
              do qq=1,3
                nodalk(mm,qq) = 0.0d0
                do p=1,9
                  avg_d = harmonic_mean(d_i(row(mm,qq,p),col(mm,qq,p)),
     &                    d_j(row(mm,qq,p),col(mm,qq,p)))
                  nodalk(mm,qq) = nodalk(mm,qq) + dnidnj(iws,p)*avg_d
                enddo
              enddo
            enddo

           !Assemble the nodal stiffness matrix into the global stiffness matrix
            a(iau+nmat(1))=a(iau+nmat(1))+nodalk(1,1)
            a(iau+nmat(2))=a(iau+nmat(2))+nodalk(1,2)
            a(iau+nmat(3))=a(iau+nmat(3))+nodalk(1,3)
                  
            a(iau+nmat(4))=a(iau+nmat(4))+nodalk(2,1)
            a(iau+nmat(5))=a(iau+nmat(5))+nodalk(2,2)
            a(iau+nmat(6))=a(iau+nmat(6))+nodalk(2,3)
                    
            a(iau+nmat(7))=a(iau+nmat(7))+nodalk(3,1)
            a(iau+nmat(8))=a(iau+nmat(8))+nodalk(3,2)
            a(iau+nmat(9))=a(iau+nmat(9))+nodalk(3,3)
         
           !Now add corresponding contributions to the diagonal block to ensure
           !zero stresses for uniform displacement i.e., [1 1 1 ...] is in null space
            a(jmia+nmat(1))=a(jmia+nmat(1))-nodalk(1,1)
            a(jmia+nmat(2))=a(jmia+nmat(2))-nodalk(1,2)
            a(jmia+nmat(3))=a(jmia+nmat(3))-nodalk(1,3)
       
            a(jmia+nmat(4))=a(jmia+nmat(4))-nodalk(2,1)
            a(jmia+nmat(5))=a(jmia+nmat(5))-nodalk(2,2)
            a(jmia+nmat(6))=a(jmia+nmat(6))-nodalk(2,3)
       
            a(jmia+nmat(7))=a(jmia+nmat(7))-nodalk(3,1)
            a(jmia+nmat(8))=a(jmia+nmat(8))-nodalk(3,2)
            a(jmia+nmat(9))=a(jmia+nmat(9))-nodalk(3,3)
          endif

        enddo
      enddo

      return
      end
