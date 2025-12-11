      subroutine fd_calc_heat(iflg)
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
!D1 Calculate fd coefficients for use with heat conduction terms of 
!D1 anisotropic perms 
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.0
!D2 
!D2 Initial implementation: Programmer: George Zyvoloski
!D2 1-May-2010
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
      use comji
      use comdti
      use comki
      use davidi
      implicit none

      integer i,icode,iflg,neqp1,i1,i2,jj,kb,kc,ic,ik,kd
      integer i3,i4,kk,idir,imodel,j,kcmin,nsize,nsize_sx
      integer kb_adv_x, kb_adv_y, kb_adv_z
      real*8 area_i,area_kb,sx2c,dism
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 cord1k,cord2k,cord3k
      real*8 dis2,dis2i,disx2, disy2, disz2
      real*8 cosxi0,cosyi0,coszi0,cosx,cosy,cosz
      real*8 cosmin,cosikc,areat,dis2min
      real*8 disikb,dis_totalkb,dis_totali
    
      neqp1 = neq +1
      if(iflg.eq.1) then 
c find areas for diffrence terms to augment  anisotropy
       ic = 0 
       nsize = nelm(neqp1)-neqp1
       if(.not.allocated(istrw))then
        allocate(istrw(nsize))
        isoy = 1
        isoz = 1
calculate storage requirements
        nsize_sx = 3*neq
        allocate(sx(nsize_sx,3))
        istrw = 0
       else
c error condition
        if(iout.ne.0)  
     &   write(iout,*) 'error condition:istrw previously allocated'    
        if(iptty.ne.0)  
     &   write(iptty,*) 'error condition:istrw previously allocated' 
       endif
       do i = 1, n0  
            areat = 0.0 
            cord1 = cord(i,1)
            cord2 = cord(i,2)
            cord3 = 0.0
            if(icnl.eq.0) then
             cord3 = cord(i,3)
            endif
            do ik = 1,3
             kd = ncon_adv(i,ik)
             if(kd.gt.0) then
             i1 = nelmdg(i)+1
             i2 = nelm(i+1)
             do jj = i1,i2
              kb = nelm(jj)
              if(kb.eq.kd) go to 1400
             enddo
1400         continue             
             if(kb.ne.i) then
              cord1j = cord(kb,1)
              cord2j = cord(kb,2)
              cord3j = 0.0
              if(icnl.eq.0) then
               cord3j = cord(kb,3)
              endif 
c calculate distance between i-kb  
              disx2 = (cord1-cord1j)**2
              disy2 = (cord2-cord2j)**2
              disz2 = (cord3-cord3j)**2
              dis2i = disx2 + disy2 + disz2
              disikb =sqrt(dis2i)
c find cosine squared              
              cosxi0 = disx2/dis2i
              cosyi0 = disy2/dis2i
              coszi0 = disz2/dis2i
c find maximum length in i-kb  (different from i-kb) direction for each gridblock
c gridblock i  (use dot product)  
            cosmin = 0.1
            kcmin = 0
            i3 = i1
            i4 = i2
            do kk = i3,i4 
             kc = nelm(kk)
             if(kc.ne.i.and.kc.ne.kb) then
              cord1k = cord(kc,1)
              cord2k = cord(kc,2)
              cord3k = 0.0
              if(icnl.eq.0) then
               cord3k = cord(kc,3)
              endif            
              disx2 = (cord1-cord1k)**2
	        disy2 = (cord2-cord2k)**2
	        disz2 = (cord3-cord3k)**2
	        dis2 = disx2 + disy2 + disz2 
c find cosine squared              
	       cosx = disx2/dis2
	       cosy = disy2/dis2
             cosz = disz2/dis2 
c compare cosines
             cosikc = abs(cosx-cosxi0)+abs(cosy-cosyi0)+abs(cosz-coszi0)
               if(cosikc.lt.cosmin) then
                cosmin = cosikc
                kcmin = kc
                dis2min = dis2
               endif
             endif
             enddo
c found distance calculate total distance
             dis_totali = 0.5*(sqrt(dis2min) +sqrt(dis2i))
c gridblock kb  (use dot product)  
            cosmin = 0.1
            kcmin = 0
            i3 = nelm(kb)+1
            i4 = nelm(kb+1)
            do kk = i3,i4 
             kc = nelm(kk)
             if(kc.ne.i.and.kc.ne.kb) then
              cord1k = cord(kc,1)
              cord2k = cord(kc,2)
              cord3k = 0.0
              if(icnl.eq.0) then
               cord3k = cord(kc,3)
              endif            
            disx2 = (cord1j-cord1k)**2
	      disy2 = (cord2j-cord2k)**2
	      disz2 = (cord3j-cord3k)**2
	      dis2 = disx2 + disy2 + disz2 
c find cosine squared              
	       cosx = disx2/dis2
	       cosy = disy2/dis2
             cosz = disz2/dis2 
c compare cosines
             cosikc = abs(cosx-cosxi0)+abs(cosy-cosyi0)+abs(cosz-coszi0)
               if(cosikc.lt.cosmin) then
                cosmin = cosikc
                kcmin = kc
                dis2min = dis2
               endif
             endif
c found distance calculate total distance
             dis_totalkb = 0.5*(sqrt(dis2min) +sqrt(dis2i))
            enddo    
c with distances and volumes can calculate areas for control volumes 
            area_i = sx1(i)/dis_totali/disikb
            area_kb = sx1(kb)/dis_totalkb/disikb
c    
            if(kb.gt.i) then  
c the "divide by 3" is for isotropic coefficients            
             ic = ic +1  
             istrw(jj-neqp1) = ic   
             sx(ic,1) = -max(area_i,area_kb)/3.
            endif
           endif  
           endif
           enddo       
       enddo
      endif
      return
      end
