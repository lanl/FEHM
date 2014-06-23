      subroutine cappr_hyd(iz,ndummy)
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Calculate capillary pressure for methane hydrate.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/cappr_hyd.f_a  $
!D2 
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 2.4.4 Relative-permeability and capillary-pressure functions
!D9
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8  Requirements from SDN: 10086-RD-2.20-00
!D8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
!D8    FEHM Application Version 2.20
!D8
!***********************************************************************
c
c subroutine to calculate the capillary pressure
c
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comdti
      use comai
      use comki
      use commeth

      implicit none

      integer i,iz,ndummy,mid,it,mi,itp,itperm,itpperm,im
      real*8 cp1,cp3

c read in data
      if(icapp.ne.0) then
         if(iz.ne.0) then
c     
c     load nodal capillary pressures
c     
            do mid=1,neq
               mi=mid+ndummy
               it=icap(mi)
               itperm = irlp(mi)
               itpperm = irlpt(itperm)
               if(itpperm .ge. 3 ) then
                  itp = 0
               else
                  itp=icapt(it)
               end if
               if(itp.eq.1) then
c     
c     linear forsythe(1989) model
c     
                  cp1=cp1f(it)
                  cp3=cp3f(it)
                  if(ieos(mi).eq.2) then
                     pcp(mi)=cp1*(cp3-fracw(mi)-frachyd(mi))
                     dpcp3(mi)=-cp1
                     dpcp4(mi)=-cp1
                     if(pcp(mi).lt.0.0) then
                        pcp(mi)=0.0
                        dpcp3(mi)=0.0
                        dpcp4(mi)=0.0
                     endif
                   else if(ieos(mi).eq.3) then
                     pcp(mi)=cp1*cp3
                     dpcp3(mi)=0.0
                     dpcp4(mi)=0.0
                   else if(ieos(mi).eq.1) then
                     pcp(mi)=0.0        
                     dpcp3(mi)=0.0
                     dpcp4(mi)=0.0
                  endif
               endif
            enddo
         end if
      endif

      return
      end
