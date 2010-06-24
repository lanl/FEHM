      subroutine area_flow_bc(iflg)
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
!D1 Calculate areas for flow through BC
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.0
!D2 
!D2 Initial implementation: 11-Nov-2009, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/area_flow_bc.f_a  $
!D2
!***********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai
      use combi
      use comci      
      use comdi
      use comdti
      use comki
      use davidi
      implicit none

      integer i,icode,iflg,neqp1,i1,i2,jj,kb
      integer i3,i4,kk,idir
      real*8 area_mult,areat,sx2c,dis,disz2,cosz,rlp_min
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 dil_dum,perm
      parameter (rlp_min=1.d-2)
      dil_dum = (997./1.e-3)
      neqp1 = neq +1
      if(iflg.eq.1) then 
         do i = 1, n0
            if (ka(i).eq.-23.or.ka(i).eq.-24) then  
c     free drainage condition         
               areat = 0.0 
               area_mult = wellim(i) 
               cord1 = cord(i,1)
               cord2 = cord(i,2)
               cord3 = 0.0
               if(icnl.eq.0) then
                  cord3 = cord(i,3)
               endif
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  if(kb.lt.i) then
                     i3 = nelmdg(kb)+1
                     i4 = nelm(kb+1)
                     do kk = i3,i4
                        if(nelm(kk).eq.i) then
                           iw = istrw(kk-neqp1)
                           go to 100
                        endif
                     enddo
                  else if(kb.gt.i) then
                     iw = istrw(jj-neqp1)
                  else
                     go to 200
                  endif
 100              cord1j = cord(kb,1)
                  cord2j = cord(kb,2)
                  cord3j = 0.0
                  if(icnl.eq.0) then
                     cord3j = cord(kb,3)
                  endif 
                  dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &                 (cord3-cord3j)**2)                   
                  disz2 = cord(kb,igrav)-cord(i,igrav)
                  cosz = abs(disz2)/dis
                  if(disz2.ge.0.0) then
                     sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
                     areat = max(sx2c*dis*cosz,areat)
                  endif
 200              continue
               enddo
               wellim(i)=areat*area_mult
            else if (ka(i).le.-11.and.ka(i).ge.-13) then  
c     constant head or pressure, physically based impedance         
               areat = 0.0 
               area_mult = wellim(i)
               idir = abs(ka(i))-10 
               if(idir.eq.1) then
                  perm = pnx(i)
               else if(idir.eq.2) then
                  perm = pny(i)
               else
                  perm = pnz(i)
               endif
               cord1 = cord(i,1)
               cord2 = cord(i,2)
               cord3 = 0.0
               if(icnl.eq.0) then
                  cord3 = cord(i,3)
               endif
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  if(kb.lt.i) then
                     i3 = nelmdg(kb)+1
                     i4 = nelm(kb+1)
                     do kk = i3,i4
                        if(nelm(kk).eq.i) then
                           iw = istrw(kk-neqp1)
                           go to 110
                        endif
                     enddo
                  else if(kb.gt.i) then
                     iw = istrw(jj-neqp1)
                  else
                     go to 210
                  endif
 110              cord1j = cord(kb,1)
                  cord2j = cord(kb,2)
                  cord3j = 0.0
                  if(icnl.eq.0) then
                     cord3j = cord(kb,3)
                  endif 
                  dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &                 (cord3-cord3j)**2)                   
                  disz2 = cord(kb,idir)-cord(i,idir)
                  cosz = abs(disz2)/dis
                  if(disz2.ge.0.0) then
                     sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
                     areat = max(sx2c*cosz,areat)
                  endif
 210              continue
               enddo
               dil_dum = dil(i)
               if(irdof.ne.13) then
                  if(rlf(i).lt.1.0) dil_dum = 
     &                 dil_dum/max(rlf(i),rlp_min)
               endif
               wellim(i)=areat*area_mult*perm*dil_dum
            endif            
         end do
      endif
      

      end
