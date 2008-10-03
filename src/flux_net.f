      subroutine flux_net(i)
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
!D1 To generate net flux for isothermal air-water solution at a node.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 09-FEB-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/flux_net.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:05:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 
!D4 
!**********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comgi
      use comji
      use comflow
      use davidi
      implicit none

      integer i
      integer icd
      integer i1
      integer i2
      integer neqp1
      integer iz
      integer jj
      integer kb
      integer kz
      integer nmatavw
      real*8 axy,vxy
      real*8 sx1d    

c changed by avw -- entered here by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif
c
c form constants for i>neq
c
      if(i.gt.neq.and.idualp.eq.0) then
         icd=neq
      else
         icd=0
      endif
      iz=i-icd
      sx1d = sx1(iz)
c
      i1=nelm(iz)+1
      i2=nelm(iz+1)
c
      if(irdof.ne.11) then
c     
c liquid phase calculations
c
            do jj = i1,i2
              kb = nelm(jj)
              if(kb.ne.iz)then
               axy=a_axy(jj-neqp1+nmatavw)
               bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
              endif
            enddo
      endif
      if(irdof.ne.13.and.jswitch.eq.0) then   
c     
c vapour phase calculations
c
c
c form equations
c
            do jj = i1,i2
              kb = nelm(jj)
              if(kb.ne.iz)then
               vxy=a_vxy(jj-neqp1+nmatavw)
               bp(iz+nrhs(2))=bp(iz+nrhs(2))+vxy
              endif
            enddo
      endif
c     
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)
      if(irdof.ne.13) then
        bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)
      endif
      
      return
      end
