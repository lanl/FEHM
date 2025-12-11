      subroutine dualta(spec_num,matnum)
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
CD1 To compute the Jacobian and residual terms of the concentration
CD1 equation at each node for a dual porosity solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dualta.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Mon Jan 29 15:42:54 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.6   08/18/95 10:20:20   llt
CD2 iw was already defined, removed for cray
CD2 
CD2    Rev 1.5   08/07/95 11:55:30   awolf
CD2 Commented out stored_residual terms and stored_derivative terms at very end on advice of BAR
CD2 
CD2    Rev 1.4   04/25/95 09:36:46   llt
CD2 retrieved log history information
CD2 
CD2    Rev 1.3   01/28/95 13:54:34   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.2   12/09/94 16:02:28   llt
CD2 Changes to fix  radi (made by gaz).
CD2
CD2    Rev 1.1   03/18/94 15:50:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:32   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3                    
CD3 spec_num              int      The current species number
CD3 matnum                int      The current submatrix number
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 lenreal, nn, neq, nelmdg, sx1, apuv1, tmch, npn, phi, pcp, 
CD4 pnx, pny, pnz, dil, div, s, anl, anv, danl, danv, rolf, rovf,
CD4 displx, disply, displz, dispvx, dispvy, dispvz, nsp, icns, igrav,
CD4 cord, dnwgta, upwgta, rc, aw, ay, denci, dencj, akc, drc, a, bp,
CD4 wb11, tb11, a21mpf, a32mpf, rb2mf, r3mf, stored_residual,
CD4 stored_derivative, npt, perml, permv, zero_t
CD4 
CD4 Global Subprograms
CD4
CD4 Name       Type     Description
CD4 
CD4 mmgetblk   N/A      Allocates memory for an array
CD4 mmrelblk   N/A      Frees array space
CD4 
C**********************************************************************
CD5 
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 itr          int         Solute unknown number
CD5 sx1d         real*8      Volume associated with the current node
CD5 pvii         real*8      Pressure of vapor
CD5 phii         real*8      Pressure of liquid
CD5 swi          real*8      Saturation
CD5 dili         real*8      Liquid transmissibility
CD5 dilkb        real*8      Liquid transmissibility, connected node
CD5 divi         real*8      Vapor transmissibility
CD5 divkb        real*8      Vapor transmissibility, connected node
CD5 anli         real*8      Liquid concentration
CD5 anlkb        real*8      Liquid concentration, connected node
CD5 anvi         real*8      Vapor concentration
CD5 anvkb        real*8      Vapor concentration
CD5 danli        real*8      Derivative of liquid concentration
CD5 danlkb       real*8      Derivative of liquid concentration
CD5 danvi        real*8      Derivative of vapor concentration
CD5 danvkb       real*8      Derivative of vapor concentration
CD5 anlri        real*8      Density times liquid concentration
CD5 anvri        real*8      Density times vapor concentration
CD5 danlri       real*8      Density times derivative of liquid
CD5                          concentration
CD5 danvri       real*8      Density times derivative of vapor
CD5                          concentration
CD5 id           int         Integer index parameter
CD5 idg          int         Integer index parameter
CD5 jm           int        Integer index parameter
CD5 neqp1        int         Number of equations plus 1
CD5 neighc       int         Integer index parameter
CD5 axkb         real*8      Permeability in x-direction, connected node
CD5 sx2c         real*8      Real parameter used in calculation
CD5 sx4d         real*8      Real parameter used in calculation
CD5 sx4h         real*8      Real parameter used in calculation
CD5 pvikb        real*8      Vapor pressure at connected node
CD5 phikb        real*8      Liquid pressure at connected node
CD5 radi         real*8      Parameter used in calculation
CD5 radkb        real*8      Parameter used in calculation
CD5 fid          real*8      Parameter used in calculation
CD5 fid1         real*8      Parameter used in calculation
CD5 axyd         real*8      Parameter used in calculation
CD5 axy          real*8      Parameter used in calculation
CD5 axyf         real*8      Parameter used in calculation
CD5 vxy          real*8      Parameter used in calculation
CD5 vxyf         real*8      Parameter used in calculation
CD5 dlaei        real*8      Derivative term used in calculation
CD5 dlaekb       real*8      Derivative term used in calculation
CD5 dvaei        real*8      Derivative term used in calculation
CD5 dvaekb       real*8      Derivative term used in calculation
CD5 pxyh         real*8      Derivative term used in calculation
CD5 heatc        real*8      Derivative term used in calculation
CD5 vxyd         real*8      Derivative term used in calculation
CD5 icode        int         Returned flag from call to mmgetblk
CD5 neq2         int         Twice the number of equations
CD5 idum         int         Index for arrays
CD5 r1           real*8      Parameter in mass balance
CD5 r2           real*8      Parameter in mass balance
CD5 r3           real*8      Parameter in mass balance
CD5 a11          real*8      Derivative parameter in mass balance
CD5 a12          real*8      Derivative parameter in mass balance
CD5 a21          real*8      Derivative parameter in mass balance
CD5 a22          real*8      Derivative parameter in mass balance
CD5 a23          real*8      Derivative parameter in mass balance
CD5 a32          real*8      Derivative parameter in mass balance
CD5 a33          real*8      Derivative parameter in mass balance
CD5 tot          real*8      Total volume at a node
CD5 idl          int         Unknown number of first matrix node
CD5 idll         int         Unknown number of second matrix node
CD5 frac0        real*8      Volume fraction of primary node
CD5 frac1        real*8      Volume fraction of first matrix node
CD5 frac2        real*8      Volume fraction of second matrix node
CD5 alen         real*8      Volume to area ratio
CD5 area         real*8      Area term used in calculation
CD5 al0          real*8      Length scale for primary node
CD5 al1          real*8      Length scale for first matrix node
CD5 al2          real*8      Length scale for second matrix node
CD5 dist01       real*8      Average length scale used in calculation
CD5 dist12       real*8      Average length scale of matrix nodes
CD5 sx1d0        real*8      Volume of primary node
CD5 sx1d1        real*8      Volume of first matrix node
CD5 sx1d2        real*8      Volume of second matrix node
CD5 coef1        real*8      Length scale used in calculation
CD5 coef2        real*8      Length scale used in calculation
CD5 i            int         Index for first matrix node number
CD5 kb           int         Index for second matrix node number
CD5 iz           int         Index for primary node number
CD5 kz           int         Index for primary node number
CD5 iw           int         Index for accessing array
CD5 axi          real*8      Permeability in x-direction
CD5 ayi          real*8      Permeability in y-direction
CD5 azi          real*8      Permeability in z-direction
CD5 alxi         real*8      Permeability in x-direction
CD5 alyi         real*8      Permeability in y-direction
CD5 alzi         real*8      Permeability in z-direction
CD5 avxi         real*8      Permeability in x-direction
CD5 avyi         real*8      Permeability in y-direction
CD5 avzi         real*8      Permeability in z-direction
CD5 diskbl       real*8      Dispersion coefficient for liquid in matrix
CD5 diskbv       real*8      Dispersion coefficient for vapor in matrix
CD5 sx2tl        real*8      Liquid dispersion parameter used in
CD5                             calculation
CD5 sx2tv        real*8      Vapor dispersion parameter used in
CD5                             calculation
CD5 pxy          real*8      Parameter used in calculation
CD5 pvxy         real*8      Parameter used in calculation
CD5 pxyi         real*8      Parameter used in calculation
CD5 ktr          int         Integer index parameter
CD5 r1c          real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 r2c          real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 r3c          real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 ab11         real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 ab22         real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 rb2          real*8      Term used in correction of residual of
CD5                             matrix nodes
CD5 spec_num     int         The current species number
CD5 id_sol       int         The array index for the residual
CD5 idum_sol     int         The array index for the a matrix
CD5 matnum       int         The current submatrix number
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6     
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.3.4 Solute-transport equations
CD9 2.4.7 Dual-porosity formulation
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN dualta
CPS 
CPS mmgetblk - allocate space in temporary arrays
CPS 
CPS FOR each node
CPS 
CPS   Initialize parameters for matrix mass balance calculations
CPS   
CPS   IF there is no volume in the second set of matrix nodes
CPS     Set parameters accordingly
CPS   ENDIF
CPS   Compute terms used later in the calculation
CPS   
CPS   IF this is a liquid phase or Henry's Law solute
CPS 
CPS     Compute liquid advection term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     Compute liquid dispersion term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     
CPS   ENDIF
CPS   
CPS   IF this is a vapor phase or Henry's Law solute
CPS 
CPS     Compute vapor advection term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     Compute vapor dispersion term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     
CPS   ENDIF
CPS   
CPS   Add accumulation term for primary node
CPS   
CPS   Initialize parameters for matrix mass balance calculations
CPS   Compute terms used later in the calculation
CPS   
CPS   IF this is a liquid phase or Henry's Law solute
CPS 
CPS     Compute liquid advection term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     Compute liquid dispersion term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     
CPS   ENDIF
CPS   
CPS   IF this is a vapor phase or Henry's Law solute
CPS 
CPS     Compute vapor advection term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     Compute vapor dispersion term for matrix nodes
CPS     Add to residual for matrix nodes
CPS     
CPS   ENDIF
CPS   
CPS   Form matrix for solution
CPS   Compute matrices to back substitute to get solution
CPS   Save arrays needed to get solution
CPS   Save residuals and derivatives in arrays
CPS   
CPS ENDFOR
CPS 
CPS mmrelblk - free temporary array space
CPS 
CPS END dualta
CPS 
C**********************************************************************

      use comcouple
      use comrxnb
      use comrxni
      use comji
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      real*8 sx1d,pvii,phii,swi,dili,dilkb,divi,divkb,anli,anlkb,anvi
      real*8 anvkb,danli,danlkb,danvi,danvkb,anlri,anvri,danlri,danvri
      real*8 axkb,sx2c,sx4d,sx4h,pvikb,phikb,radi,radkb,fid,fid1,axyd
      real*8 axyf,axy,vxy,vxyf,dlaei,dlaekb,dvaei,dvaekb,pxyh,heatc,vxyd
      real*8 r1,r2,r3,a11,a12,a21,a22,a23,a32,a33,tot,frac0,frac1,frac2
      real*8 alen,area,al0,al1,al2,dist01,dist12,sx1d0,sx1d1,sx1d2,coef1
      real*8 axi,coef2,ayi,azi,alxi,alyi,alzi,avxi,avyi,avzi,diskbl
      real*8 diskbv,sx2tl,sx2tv,pxy,pvxy,pxyi,r1c,r2c,r3c,ab11,ab22,rb2
      integer id,itr,idg,jm,neqp1,neighc,icode,neq2,idum,idl,idll
      integer i,kb,iz,kz,ktr,spec_num,id_sol,idum_sol,matnum

      neq2=neq+neq
      neqp1=neq+1
c     
c modify geometric terms
c save nelm terms for later use
c
c loop on nodes
      do id=1,neq
c identify diagonal member
         idg=nelmdg(id)
         idum=idg-neqp1
         id_sol = id + nrhs_sol(spec_num)
         idum_sol = idum + nmat_sol(matnum)

         r1=0.0
         r2=0.0
         r3=0.0
         a11=0.0
         a12=0.0
         a21=0.0
         a22=0.0
         a23=0.0
         a32=0.0
         a33=0.0
c     
c form equations at id+neq2
c     
         idl=id+neq
         idll=id+neq2
c
c generate transport terms
c
         tot =  sx1(id)+sx1(idl)+sx1(idll)
         frac0 =  sx1(id)/tot
         frac1 =  sx1(idl)/tot
         frac2 =  sx1(idll)/tot
         alen =  apuv1(id)
         area=tot/alen
         al0 =  frac0*alen
         al1 =  frac1*alen
         al2 =  frac2*alen
         dist12=0.5*(al1+al2)
         coef2=-area/dist12
         sx1d0=sx1(id)
         sx1d1=sx1(idl)
         sx1d2=sx1(idll)
         if(frac2.le.tmch) then
            coef2=0.0
            sx1d2=1.0
         endif

         i=idl
         kb=idll
         iz=id
         kz=id
         iw=1

         itr=i+npn
         ktr=kb+npn
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
         phii=pvii-pcp(i)
         swi=s(i)
         dili=dil(i)
         divi=div(i)
         anli=anl(itr)
         anvi=anv(itr)
         danli=danl(itr)
         danvi=danv(itr)
         anlri=anli*rolf(i)
         anvri=anvi*rovf(i)
         danlri=danli*rolf(i)
         danvri=danvi*rovf(i)

c form constants for i>neq
c     
         axkb =  max( pnx(idll  ),pny(idll  ),
     *        pnz(idll  ),zero_t )
         diskbl =  max( displx(idll  ),disply(idll  ),
     *        displz(idll  ),zero_t )
         diskbv =  max( dispvx(idll  ),dispvy(idll  ),
     *        dispvz(idll  ),zero_t )

c     2-d geometry
c     
         radi = 1.
         jm=1
         neighc=1
         perml(1)=axkb
         permv(1)=perml(1)
         radkb=radi
         sx2c=radkb*coef2
         sx2tl=sx2c*diskbl
         sx2tv=sx2c*diskbv
         pvikb=phi(kb)
         phikb=pvikb-pcp(kb)
         pxy=sx2c*perml(1)
         pvxy=sx2c*permv(1)
         pxyi=pxy*(phikb-phii)
         pxyh=pvxy*(pvikb-pvii)
         t1(neighc)=pxyi
         t2(neighc)=pxyh
         t3(neighc)=pxy
         t4(neighc)=pvxy
         t5(neighc)=sx2tl
         t5v(neighc)=sx2tv
         t6(neighc)=0.0
         t7(neighc)=0.0

c     choose vapour or liquid tracer
c     
******
         if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2) then
c     
c     liquid phase calculations
c     
            jm=1
            neighc=1
            pxyi=t1(neighc)
            sx4d=t6(neighc)
            axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=axyd
c     
c     determine upwind nodes and if liquid phase exists
c     
            jm=1
            neighc=1
            fid=.5
            axyd=t8(neighc)
            if(axyd.lt.0.0) fid=dnwgta
            if(axyd.gt.0.0) fid=upwgta
            t9(neighc)=fid

c     form equations
c     
            jm=1
            neighc=1
            axyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            dilkb=dil(kb)
            anlkb=anl(ktr)
            danlkb=danl(ktr)
            axyf=(fid*dilkb*anlkb+fid1*dili*anli)
            axy=axyd*axyf
            dlaei=axyd*fid1*dili*danli
            dlaekb=axyd*fid*dilkb*danlkb

            r2=r2+axy
            r3=r3-axy
            a22=a22+dlaei
            a32=a32-dlaei
            a23=a23+dlaekb
            a33=a33-dlaekb
c     
c add dispersion term for liquid
c
            jm=1
            neighc=1
            heatc=t5(neighc)
            axy=heatc*(anl(ktr)*rolf(kb)-anlri)
            dlaei=-heatc*danlri
            dlaekb=heatc*rolf(kb)*danl(ktr)
            r2=r2+axy
            r3=r3-axy
            a22=a22+dlaei
            a32=a32-dlaei
            a23=a23+dlaekb
            a33=a33-dlaekb
         endif
         if (icns(nsp).eq.-1.or.abs(icns(nsp)).eq.2) then

c     
c     vapour phase calculations
c     
            jm=1
            neighc=1
            pxyh=t2(neighc)
            sx4h=t7(neighc)
            vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=vxyd
c     
c     determine upwind nodes and if vapour phase exists
c     
            jm=1
            neighc=1
            fid=.5
            vxyd=t8(neighc)
            if(vxyd.lt.0.0) fid=dnwgta
            if(vxyd.gt.0.0) fid=upwgta
            t9(neighc)=fid

c     form equations
c     
            jm=1
            neighc=1
            vxyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            divkb=div(kb)
            anvkb=anv(ktr)
            danvkb=danv(ktr)
            vxyf=(fid*divkb*anvkb+fid1*divi*anvi)
            vxy=vxyd*vxyf
            dvaei=vxyd*fid1*divi*danvi
            dvaekb=vxyd*fid*divkb*danvkb
            r2=r2+vxy
            r3=r3-vxy
            a22=a22+dvaei
            a32=a32-dvaei
            a23=a23+dvaekb
            a33=a33-dvaekb

c     add dispersion term
c     
            jm=1
            neighc=1
            heatc=t5v(neighc)
            vxy=heatc*(anv(ktr)*rovf(kb)-anvri)
            dvaei=-heatc*danvri
            dvaekb=heatc*rovf(kb)*danv(ktr)

            r2=r2+vxy
            r3=r3-vxy
            a22=a22+dvaei
            a32=a32-dvaei
            a23=a23+dvaekb
            a33=a33-dvaekb

         endif
****  rc is used ****
c
c     add accumulation terms for second and first matrix levels
c     
         sx1d=sx1d1
         r2=r2+sx1d*(aw*denci(itr)+ay*dencj(itr))+rc(itr)
         a22=a22+sx1d*(aw*akc(itr))+drc(itr)

         sx1d=sx1d2
         r3=r3+sx1d*(aw*denci(ktr)+ay*dencj(ktr))+rc(ktr)
         a33=a33+sx1d*(aw*akc(ktr))+drc(ktr)
c     
c     form equations at node id
c     
c     geometry
c
         dist01=0.5*(al0+al1)
         coef1=-area/dist01
         if(frac1.le.0.0) then
            coef1=0.0
         endif

         i=id
         kb=idl
         iz=id
         kz=id

         itr=i+npn
         ktr=kb+npn
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
         phii=pvii-pcp(i)
         dili=dil(i)
         divi=div(i)
         anli=anl(itr)
         anvi=anv(itr)
         danli=danl(itr)
         danvi=danv(itr)
         anlri=anli*rolf(i)
         anvri=anvi*rovf(i)
         danlri=danli*rolf(i)
         danvri=danvi*rovf(i)

c form constants for i>neq
c
         axkb =  max( pnx(idl  ),pny(idl  ),
     *        pnz(idl  ),zero_t )
         diskbl =  max( displx(idl  ),disply(idl  ),
     *                                      displz(idl  ),zero_t )
         diskbv =  max( dispvx(idl  ),dispvy(idl  ),
     *                                      dispvz(idl  ),zero_t )

c     2-d geometry
c     
         radi = 1.
         jm=1
         neighc=1
         perml(1)=axkb
         permv(1)=perml(1)
         radkb=radi
         sx2c=radkb*coef1
         sx2tl=sx2c*diskbl
         sx2tv=sx2c*diskbv
         pvikb=phi(kb)
         phikb=pvikb-pcp(kb)
         pxy=sx2c*perml(1)
         pvxy=sx2c*permv(1)
         pxyi=pxy*(phikb-phii)
         pxyh=pvxy*(pvikb-pvii)
         t1(neighc)=pxyi
         t2(neighc)=pxyh
         t3(neighc)=pxy
         t4(neighc)=pvxy
         t5(neighc)=sx2tl
         t5v(neighc)=sx2tv
         t6(neighc)=0.0
         t7(neighc)=0.0

c     choose vapour or liquid tracer
c     
         if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2) then
c     
c     liquid phase calculations
c
            jm=1
            neighc=1
            pxyi=t1(neighc)
            sx4d=t6(neighc)
            axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=axyd
c
c     determine upwind nodes and if liquid phase exists
c
            jm=1
            neighc=1
            fid=.5
            axyd=t8(neighc)
            if(axyd.lt.0.0) fid=dnwgta
            if(axyd.gt.0.0) fid=upwgta
            t9(neighc)=fid

c form equations
c
            jm=1
            neighc=1
            axyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            dilkb=dil(kb)
            anlkb=anl(ktr)
            danlkb=danl(ktr)
            axyf=(fid*dilkb*anlkb+fid1*dili*anli)
            axy=axyd*axyf
            dlaei=axyd*fid1*dili*danli
            dlaekb=axyd*fid*dilkb*danlkb

            r1=r1+axy
            r2=r2-axy
            a11=a11+dlaei
            a21=a21-dlaei
            a12=a12+dlaekb
            a22=a22-dlaekb
c
c add dispersion term for liquid
c
            jm=1
            neighc=1
            heatc=t5(neighc)
            axy=heatc*(anl(ktr)*rolf(kb)-anlri)
            dlaei=-heatc*danlri
            dlaekb=heatc*rolf(kb)*danl(ktr)
            r1=r1+axy
            r2=r2-axy
            a11=a11+dlaei
            a21=a21-dlaei
            a12=a12+dlaekb
            a22=a22-dlaekb
            endif
            if (icns(nsp).eq.-1.or.abs(icns(nsp)).eq.2) then

c     vapour phase calculations
c     
            jm=1
            neighc=1
            pxyh=t2(neighc)
            sx4h=t7(neighc)
            vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *           *(cord(kz,igrav)-cord(iz,igrav))
            t8(neighc)=vxyd
c     
c determine upwind nodes and if vapour phase exists
c
            jm=1
            neighc=1
            fid=.5
            vxyd=t8(neighc)
            if(vxyd.lt.0.0) fid=dnwgta
            if(vxyd.gt.0.0) fid=upwgta
            t9(neighc)=fid

c     form equations
c     
            jm=1
            neighc=1
            vxyd=t8(neighc)
            fid=t9(neighc)
            fid1=1.0-fid
            divkb=div(kb)
            anvkb=anv(ktr)
            danvkb=danv(ktr)
            vxyf=(fid*divkb*anvkb+fid1*divi*anvi)
            vxy=vxyd*vxyf
            dvaei=vxyd*fid1*divi*danvi
            dvaekb=vxyd*fid*divkb*danvkb
            r1=r1+vxy
            r2=r2-vxy
            a11=a11+dvaei
            a21=a21-dvaei
            a12=a12+dvaekb
            a22=a22-dvaekb

c add dispersion term
c
            jm=1
            neighc=1
            heatc=t5v(neighc)
            vxy=heatc*(anv(ktr)*rovf(kb)-anvri)
            dvaei=-heatc*danvri
            dvaekb=heatc*rovf(kb)*danv(ktr)

            r1=r1+vxy
            r2=r2-vxy
            a11=a11+dvaei
            a21=a21-dvaei
            a12=a12+dvaekb
            a22=a22-dvaekb
         endif

c     form matrix ab22=a22-a23*a33i*a32
c     
*         rc is used here
         r1c=r1
         r2c=r2
         r3c=r3
         ab22=a22-a23*a32/a33
         rb2=r2c-a23*r3c/a33
c     
c     form matrix ab11=a11-a12*ab22i*a21
c     
         ab11=a11-a12*a21/ab22
         a(idum)=a(idum)+ab11
         bp(id)=bp(id)+r1c-a12*rb2/ab22
c     
c     save matrices to back out solution
c     ab22i=wb,a33i=tb,a21,a32
c     
         wb11(id)=1.0/ab22
         tb11(id)=1.0/a33
         a21mpf(id)=a21
         a32mpf(id)=a32
c     
c save vectors neccessary to get solution
c
         rb2mf(id)=rb2
         r3mf(id)=r3c

      enddo

      return
      end
