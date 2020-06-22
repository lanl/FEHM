      subroutine accum_term(mid, ndummy)       
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
!D1 To calculate water and air accumulation terms.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/accum_term.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:28:10   pvcs
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

      use comdti
      use comai
c      use combi
      use comci
      use comdi
c      use comei
c      use comfi
c      use comgi
      use comii
      use davidi
c      use comrxni
      implicit none
c
C      integer ndummy,mid,mi,ieosd,kq,iflg
      integer ndummy,mid,mi,ieosd
C      real*8 dtin,tempc,drocp,drocp0,rolref,xvisl,comw,pref,xvisv
c gaz 110819 pref now global      
      real*8 tempc,drocp0,rolref,xvisl,comw,xvisv
C      real*8 xrl,drl,drlp,xrv,drv,drvp,pl,sl,svd,pwl,dpwlp,dpwls
      real*8 pl,sl,svd,pwl
C      real*8 qdis,qwdis,qadis,dqws,dqas,dqwp,dqap,por,sflux,permsd
      real*8 por
C      real*8 pflowd,perml,dpls,dpas,perma,roc,drocs,rol,drolp,dporpl
      real*8 roc,rol
C      real*8 rcomd,drols,dena,ddenp,ddens,ddenap,ddenas,dql,dqv
      real*8 rcomd,dena
C      real*8 pcl0,roc0,wimped,airmobile,conwa,conaw
c      real*8 pcl0,roc0
      real*8 pcl0
      real*8 cden_correction, cden_cor
      parameter(pcl0 = 0.101325)
C      parameter(airmobile = 10.0)
C      integer i_mem_rlp
C      save i_mem_rlp
     

c     
c     dependent variables vap p and sl
c     
c     misc. constants
c gaz 110819 (t() can be spatially variable-use t(mi)
c for now use roc0 and pcl0 at atmospheric definitions
      tempc=(273.0)/(t(mi)+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
      xvisl=crl(2,1)
      comw=crl(3,1)
c gaz 110819 pref, tref (global) read in scanin        
c      pref=crl(4,1)
      xvisv=crl(5,1)
      rcomd=comw*rolref
c     
c     generate coef and derivatives
c     
         mi=mid
         ieosd=2          
         pl=phi(mi)
         sl=s(mi)
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
c
         pwl=pl-pcp(mi)
         por=ps(mi)
c     density of air
         roc=drocp0*phi(mi)
c     water density
         rol=rolref*(1.0+comw*(pl-pref))
         if(cden) then
            cden_cor = cden_correction(mi)
            rol = rol + cden_cor
         end if
c     accumulation terms
         den=por*rol*sl
         dena=por*roc*svd
c     
         denh(mi)=den                  
         deneh(mi)=dena                  
      return
      end
