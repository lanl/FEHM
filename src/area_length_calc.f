      subroutine area_length_calc(iflg)
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
!D1 Calculate areas and lengths for GDPM calcs 
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.0
!D2 
!D2 Initial implementation: ?, Programmer: George Zyvoloski
!D2 11-Nov-2009
!D2 $Log:   /pvcs.config/fehm90/src/inflo3.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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
      use comriv
      use davidi
      implicit none

      integer i,id,icode,iflg,neqp1,i1,i2,jj,kb,kc
      integer i3,i4,kk,idir,imodel,j,n_loop
      real*8 area_mult,areat,sx2c,dis,cosz,rlp_min,rad
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 dil_dum,perm, disx1, disy1, disz1, disx2, disy2, disz2
      parameter (rlp_min=1.d-2)
c isothermal default density/viscosity 
c might want to make a function of fluid and conditions     
      dil_dum = 997.d0/1.d-3
      neqp1 = neq +1
      if(iflg.eq.1) then 
c find areas for BCs, physically based impedances      
       do i = 1, n0
        if((ka(i).eq.-1.or.ka(i).eq.-2).and.ianpe.eq.0) then
c general fixed pressure or head BC       
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
100            cord1j = cord(kb,1)
               cord2j = cord(kb,2)
               cord3j = 0.0
               if(icnl.eq.0) then
                cord3j = cord(kb,3)
               endif 
               dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &              (cord3-cord3j)**2)                   
                sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
                areat = max(sx2c*dis,areat)
 200          continue
            enddo
            wellim(i)=areat*area_mult
          else if((ka(i).eq.-1.or.ka(i).eq.-2).and.ianpe.ne.0) then
c general fixed pressure or head BC with anisotropy
c bc impedance is multiple of the volume              
           area_mult = wellim(i) 
           wellim(i)=area_mult*sx1(i)            
          else if((ka(i).eq.-23.or.ka(i).eq.-24).and.ianpe.eq.0) then
c free drainage condition         
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
                   go to 101
                  endif
                enddo
              else if(kb.gt.i) then
                iw = istrw(jj-neqp1)
              else
               go to 201
              endif
101            cord1j = cord(kb,1)
               cord2j = cord(kb,2)
               cord3j = 0.0
               if(icnl.eq.0) then
                cord3j = cord(kb,3)
               endif 
               dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &              (cord3-cord3j)**2)                   
               disz2 = cord(kb,igrav)-cord(i,igrav)
               cosz = abs(disz2)/dis
               if(cosz.gt.0.0) then
                sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
                areat = max(sx2c*dis*cosz,areat)
               endif
 201          continue
            enddo
            wellim(i)=areat*area_mult  
         else if((ka(i).eq.-23.or.ka(i).eq.-24).and.ianpe.ne.0) then
c free drainage with anisotropy
c bc impedance is multiple of the volume              
           area_mult = wellim(i) 
           wellim(i)=area_mult*sx1(i)  
         else if ((ka(i).le.-11.and.ka(i).ge.-13).and.ianpe.eq.0) then  
c constant head or pressure, physically based impedance         
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
110            cord1j = cord(kb,1)
               cord2j = cord(kb,2)
               cord3j = 0.0
               if(icnl.eq.0) then
                cord3j = cord(kb,3)
               endif 
               dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &              (cord3-cord3j)**2)                   
               disz2 = cord(kb,idir)-cord(i,idir)
               cosz = abs(disz2)/dis
               if(cosz.ge.0.0) then
                sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
                areat = max(sx2c*cosz,areat)
               endif
 210          continue
            enddo
            dil_dum = dil(i)
            if(irdof.ne.13) then
             if(rlf(i).lt.1.0) dil_dum = dil_dum/max(rlf(i),rlp_min)
            endif
            wellim(i)=areat*area_mult*perm*dil_dum
         else if ((ka(i).le.-11.and.ka(i).ge.-13).and.ianpe.ne.0) then    
c bc impedance is multiple of the volume              
           area_mult = wellim(i) 
           wellim(i)=area_mult*sx1(i)  
         endif               
       end do
      else if((iflg.eq.2.and.gdkm_flag.ne.0).and.ianpe.eq.0) then
c calculate some areas and lengths for GDPM and GDKM calculations   
       do i = 1,neq_primary
        imodel= igdpm(i)
        if(imodel.ne.0) then
          cord1 = cord(i,1)
	    cord2 = cord(i,2)
	    cord3 = 0.0
	    if(icnl.eq.0) then
           cord3 = cord(i,3)
          endif
c find areas of all connecting (complete) gridblocks
          i1 = nelm(i)+1
          i2 = nelm(i+1)
          j = 0
          do jj = i1,i2
           kb = nelm(jj)
           iw = istrw(jj-neqp1)
           j = j+1
           cord1j = cord(kb,1)
           cord2j = cord(kb,2)
           cord3j = 0.0
           if(icnl.eq.0) then
            cord3j = cord(kb,3)
           endif 
           sx2c=abs(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
           dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &          (cord3-cord3j)**2)      
           it10(j) = kb
           t1(j) = sx2c*dis
          enddo  
c examine connections to neighbors               
          do j = 1, ngdpm_layers(imodel)
c note special numbering for gdpm nodes           
           kc =  kc + 1
           cord1j = cord(kb,1)
	     cord2j = cord(kb,2)
           cord3j = cord(kb,3)
         enddo
        endif
       enddo
      else if((iflg.eq.2.and.gdkm_flag.ne.0).and.ianpe.ne.0) then 
       if(iout.ne.0) then
        write(iout,*)'error:area_length_calc.f(no capability for anpe)'
        write(iout,*)'stopping'
       endif
       if(iptty.ne.0) then
        write(iptty,*)'error:area_length_calc.f(no capability for anpe)'
        write(iptty,*)'stopping'
       endif    
       stop   
      else if(iflg.eq.3.and.ianpe.eq.0) then
c calculate cell lengths for GDPM and GDKM calculations   
       if(.not.allocated(dzrg)) allocate(dzrg(n))
       if(.not.allocated(dyrg)) allocate(dyrg(n))
       if(.not.allocated(dxrg)) allocate(dxrg(n))   
       if(gdkm_flag.ne.0) then   
        n_loop = neq_primary
       else if(idpdp.ne.0) then
        n_loop = neq_primary
       else
        n_loop = n
       endif
       do i = 1, n_loop
          cord1 = cord(i,1)
	    cord2 = cord(i,2)
	    cord3 = 0.0
	    if(icnl.eq.0) then
           cord3 = cord(i,3)
          endif
c find lengths of all connecting primary gridblocks
          disx1 = 0.0
          disy1 = 0.0
          disz1 = 0.0      
          disx2 = 1.e20
          disy2 = 1.e20
          disz2 = 1.e20
          i1 = nelm(i)+1
          i2 = nelm(i+1)
          do jj = i1,i2
           kb = nelm(jj)
           cord1j = cord(kb,1)
           cord2j = cord(kb,2)
           cord3j = 0.0
           if(icnl.eq.0) then
            cord3j = cord(kb,3)
           endif 
              disx1=max(cord1j-cord1,disx1)
              disx2=min(cord1j-cord1,disx2)
              disy1=max(cord2j-cord2,disy1)
              disy2=min(cord2j-cord2,disy2) 
              disz1=max(cord3j-cord3,disz1)
              disz2=min(cord3j-cord3,disz2)                           
          enddo  
            if(ivf.eq.-1) then
               if(disx1.eq.0.0.and.disx2.eq.0.0)
     &             disx1 = abs(cord1-x_orig)*2.
               if(disy1.eq.0.0.and.disy2.eq.0.0) 
     &             disy1 = abs(cord2-y_orig)*2.
               if(disz1.eq.0.0.and.disz1.eq.0.0)
     &             disz1 = abs(cord3-z_orig)*2.
               dzrg(i) = max(disz1,abs(disz2))
               dyrg(i) = max(disy1,abs(disy2))
               dxrg(i) = max(disx1,abs(disx2))
               continue               
            else
               dzrg(i) = abs(disz1-disz2)/2.	
               dyrg(i) = abs(disy1-disy2)/2.    
               dxrg(i) = abs(disx1-disx2)/2.                        
            endif      
       enddo 
       if(iriver.eq.2) then
        do i = 1,nodes_well2_added
c identify segment and radius        
         id = new_node_well2_segid(i)
         rad = well_rad(id)
         dzrg(i+neq_primary) = 2*rad
        enddo
       endif
       continue
c       
c  lengths for gdkm set (for now) equal to that of primary node
c
       if(gdkm_flag.ne.0) then
          do i = 1,neq_gdkm
           imodel= igdpm(i)
c     first gdkm node is given last connection of primary node
            if(imodel.ne.0) then
             kb = nelm(nelm(i+1))
             corz(kb,1) = cord(i,1)
             corz(kb,2) = cord(i,2)
             corz(kb,3) = cord(i,3)
            endif
          enddo
       endif  
      
      else if(iflg.eq.3.and.ianpe.ne.0) then 
       if(iout.ne.0) then
        write(iout,*)'warning:area_length_calc.f(no capability w anpe)'
        write(iout,*)'gdkm and anpe not compatible'
       endif
       if(iptty.ne.0) then
        write(iptty,*)'warning:area_length_calc.f(no capability w anpe)'
        write(iptty,*)'gdkm and anpe not compatible'
       endif 
c       stop          
      else if(iflg.eq.-3.and.gdkm_flag.ne.0) then
c calculate cell lengths for GDPM and GDKM calculations  
c deallocations 
       if(allocated(dzrg)) deallocate(dzrg)
       if(allocated(dyrg)) deallocate(dyrg)
       if(allocated(dxrg)) deallocate(dxrg)   
      endif
      end
