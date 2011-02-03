      subroutine compute_strains() 
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
! Computes the strain tensor using the FEHM approach. Only called when the
! 'fem' submacro is not used.
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
      
      if(icons.le.abs(maxit)) then
       grav_air=0.0
      else
       grav_air=grav
      endif

      do i=1,neq 

        vol=sx1(i)

!        e1i = e1(i)
!        e2i = e2(i)
!        e3i = e3(i)

        dui = du(i)
        dvi = dv(i)
        dwi = dw(i)

        if(i.gt.neq.and.idualp.eq.0) then
           icd=neq
        else
           icd=0
        endif
        iz=i-icd

        iz4m1 = 4*(iz-1)+1

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

        do jm=1,iq
          kb=it8(jm)
          kz=kb-icd          
        
          iws=itstress(jm)    
          
          !! Signs in shape function integrals was changed to make compression positive
          !! So, using negative sign below
          sisjx= -sxs(iws,7)
          sisjy= -sxs(iws,8)
          sisjz=  sxs(iws,9)

          ddx=du(kb)-dui
          ddy=dv(kb)-dvi
          ddz=dw(kb)-dwi
c      
          termx = termx + sisjx*ddx
          termy = termy + sisjy*ddy
          termz = termz + sisjz*ddz

          tauxy = tauxy + sisjy*ddx + sisjx*ddy
          tauyz = tauyz + sisjz*ddy + sisjy*ddz
          tauzx = tauzx + sisjx*ddz + sisjz*ddx

        enddo
      
        strain_xx(i)=termx/vol
        strain_yy(i)=termy/vol
        strain_zz(i)=termz/vol
        strain_xy(i)=tauxy/vol  
        strain_yz(i)=tauyz/vol
        strain_zx(i)=tauzx/vol
     
      enddo
 
      return
      end


