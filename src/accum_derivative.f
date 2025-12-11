      subroutine accum_derivative(mid,iflg) 
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
!D1 PURPOSE
!D1
!D1 To calculate accumulation terms and derivatives for isothermal
!D1 air-water system.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/accum_derivative.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:27:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of thrair 08-Feb-02
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      use comrxni
      use comsplitts
      use davidi
      use comii
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer mid,mi,ieosd,kq,iflg
c gaz 110819 removed tref, pref (now global)         
      real*8 dtin,tempc,drocp,drocp0,rolref,xvisl,comw,xvisv
      real*8 xrl,drl,drlp,xrv,drv,drvp,pl,sl,svd,pwl,dpwlp,dpwls
      real*8 qdis,qwdis,qadis,dqws,dqas,dqwp,dqap,por,sflux,permsd
      real*8 pflowd,area,uperm,roc,drocs,rol,drolp,dporpl
      real*8 rcomd,drols,dena,ddenp,ddens,ddenap,ddenas,dql,dqv
      real*8 pcl0, rolref_b    
      real*8 tl_last,viln,vild,vil
      real*8 viln1,viln2,viln3,vild1,vild2,vild3
      real*8 vla0,vlpa1,vlpa2,vlpa3,vlta1,vlta2,vlta3,vlpta,vlp2ta
      real*8 vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3,vlptb,vlp2tb
      real*8 vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tlx,tl2x,tlx2
c cden_correction is a real function
      real*8 cden_correction
      integer i_mem_rlp
      save i_mem_rlp
      parameter(pcl0 = 0.101325)
c gaz debug 011114 roc0 now on comai
c      parameter(roc0 = 1.292864)
     
c     calculate variable porosity if enabled
c      if(iporos.ne.0) call porosi(1)
c     
c     dependent variables vap p and sl
c     
c     misc. constants
c gaz 110819 removed tref, pref (now global)      
c      tref = crl(6,1)
      tempc=(273.0)/(tref+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
      xvisl=crl(2,1)
      comw=crl(3,1)
c gaz 110819 pref, tref (global) read in scanin      
c      pref=crl(4,1)
      xvisv=crl(5,1)

c
c     liquid viscosity
c
c     numerator coefficients
c     
      dtin=1.0/dtot
c     
c     generate coef and derivatives
c     
c called  for  every grid block
C         mi=mid+ndummy
         mi=mid
         ieosd=2    
       if(cden)then
        rolref_b= rolref+cden_correction(mi)
       else if(iden_vis.gt.0) then
c spatially variable density and viscosity
        rolref_b = den_spatial(mi)
        xvisl = vis_spatial(mi)
        if(comp_spatial(mi).gt.0.0) then
          comw = comp_spatial(mi)
        endif
       else
        rolref_b= rolref
       endif
       rcomd=comw*rolref_b
      
c     
         pl=phi(mi)
         sl=s(mi)
         por = ps(mi)
c     
          if(irdof.ne.13) then
            pwl=pl-pcp(mi)
            dpwlp=1.0
            dpwls=-dpcef(mi)
          else
            pwl=pl
            dpwlp=1.0d00
            dpwls=0.0d00      
          endif
c     
c
c add new coding for sv
c
c next line used to carry sv info for irdof=14
c
        if(abs(irdof).eq.14) then
         svd=denj(mi)   
        else
         svd=1.0-s(mi)
        endif
****  Average molecular weight is equal to the molecular
****  weight of air, since no water is present in the vapor phase
****  in this routine ****
         avgmolwt(mi) = mw_air
c     density of air
c GAZ had to move it here for constant density approx
c added the derivative of porosity wrt pressure
         dporpl=dporp(mi)
         drocp=drocp0                     
         roc=drocp0*pl
         drocs=0.0
c     water density

          rol=rolref_b*(1.0+comw*(pl-pref))
          drolp=rcomd
          drols=0.0
c     accumulation terms
          den=por*rol*sl
          dena=por*roc*svd
          ddenp=por*sl*drolp+dporpl*rol*sl 
          ddens=por*rol
          ddenap=por*drocp*svd + dporpl*roc*svd
          ddenas=por*(drocs*svd-roc)         
          if(ifree.ne.0) then
           ddenp = ddenp + por*rol*drlxyf(mi)
          endif
c     
         deni(mi)=(den-denh(mi))*dtin
         dmpf(mi)=ddenp*dtin
         if(icons.lt.maxit) then
c 
c make densities constant for transport terms
c
c     density of air
         roc=drocp0*pref  
         drocs=0.0
         drocp=0.0
c     water density
         rol=rolref
         drolp=0.0
         drols=0.0		 
         endif
         if(irdof.ne.13) then
           rovf(mi)=roc
           denei(mi)=(dena-deneh(mi))*dtin
           dmef(mi)=ddens*dtin
           depf(mi)=ddenap*dtin
           deef(mi)=ddenas*dtin
         endif
      return
      end
