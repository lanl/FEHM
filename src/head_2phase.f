      subroutine head_2phase(iflg)
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
CD1  This subroutine sets air and water pressure with 2 phase conditions
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 FEHM Version 2.20
CD2 
CD2 Initial implementation: ?, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/head_2phase.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.3.2 Heat- and mass-transfer equations
CD3 2.3.3 Noncondensible gas flow equations 
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  Finite Element Ouput Format
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************
c sets air and water pressure with 2 phase conditions

      use comai
      use combi
      use comdi
      use comdti
      use comfi
      use comgi
      use comii
      use comrxni
      implicit none

      integer iflg, neq2
      integer mi,mid
      real*8 rho2grav,hmax,drocp0
c gaz 110819 removed tref, pref (now global)        
      real*8 pcl0, rolref, tempc,roc
      parameter(pcl0 = 0.101325)
c gaz debug 011014
c      parameter(roc0 = 1.292864)
c gaz 110819 tref  read in scanin        
c      tref = crl(6,1)
      tempc=(273.0)/(tref+273.0)
      drocp0=roc0*tempc/pcl0
      rolref=crl(1,1)
c gaz 110819 pref  read in scanin     
c      pref=crl(4,1)
      roc = drocp0*pref
      rho2grav = roc*(-grav)
      neq2 = 2*neq

      if(iflg.eq.0.and.ifree.eq.0) then
c     correct for negative pressures (head solution)
c     first find max height
         hmax= -1.e30
         mi = 0
         do mid = 1,neq
            if(cord(mid,igrav).ge.hmax) then
               hmax = cord(mid,igrav)
               mi = mid
            endif
         enddo
         do mi = 1,n
            if(mi.gt.neq2) then
               mid=mi-neq2
            elseif(mi.gt.neq) then
               mid=mi-neq
            else
               mid=mi
            endif
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref             
            if(pho(mi).lt.pref.and.l.eq.0) then
c     correct initial pressure
               pho(mi)= (hmax-cord(mid,igrav))*rho2grav + pref
c     change to seepage face condition if specified head
            endif
            if(pflow(mi).lt.pref) then
               pflow(mi)= (hmax-cord(mid,igrav))*rho2grav + pref
c     if(ka(mi).lt.0) ka(mi) = -3
            else 
c     s(mi) = 1.0
            endif
            
         enddo
      endif
      return
      end
