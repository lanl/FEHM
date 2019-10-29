      subroutine gdkm_connect(iflg)
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
!D1 Identify connections and geometry for gdkm 
!D1 Called from add_gdpm to connect 'last' gdpm node to a surrounding 
!D1    primary grid matrix block
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.0
!D2 
!D2 Initial implementation: Programmer: George Zyvoloski
!D2 8-May-2010
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

      integer i,icode,iflg,neqp1,i1,i2,jj,kb,kc,ii
      integer i3,i4,kk,idir,imodel,j,kcmin,iparchek
      real*8 area_i,area_kb,sx2c,dism,areat,tol_cos
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 dis2,dis2i,disx2, disy2, disz2
      real*8 cosxi0,cosyi0,coszi0,cosx,cosy,cosz
      real*8 cosmin,cosikc
      real*8 cosneqx, cosposx
      integer kb_neg,kb_pos,i_pos,kb_pos_neg,i_neg,jk
      integer ic,ngdkm,nlayerm,kblast,nsize,kbmax,kbmin
      integer i_gdkm,i1_gdkm,i2_gdkm, iww
      real*8 cord_kb_pos_neg, cord_i_pos, area_kb_pos
      real*8  cord_i_neg, area_kb_neg
      real*8 disx_max,disx_min
c gaz 060117 added volume fraction cals      
      real*8 a11,vfrac,vfrac2,sx1_primary,sx1_total,vf_gdpm
      integer, allocatable :: idum1(:)
      integer, allocatable :: idump(:)
      integer, allocatable :: idumn(:)
      integer, allocatable :: nelm_temp(:)
      integer, allocatable :: istrw_temp(:)
      real*8, allocatable :: sx_12(:,:)
      real*8, allocatable :: sx_temp(:,:)
      parameter(iparchek = 1)
      parameter(tol_cos = 1.d-8)
c
      neqp1 = neq +1
      if(iflg.eq.1) then 
c find connections in directions gdkm
c extra storage used later
       allocate(iconn_gdkm(neq_primary,4))
        if(gdpm_flag.eq.11) then       
        do i = 1,neq_primary
            areat = 0.0 
            cord1 = cord(i,1)
            cord2 = cord(i,2)
            cord3 = 0.0
            if(icnl.eq.0) then
             cord3 = cord(i,3)
            endif
c     loop on neighbors 
            cosneqx =0.0
            cosposx =0.0  
            kb_neg = 0
            kb_pos = 0       
            i1 = nelm(i)+1
            i2 = nelm(i+1)
            do jj = i1,i2
             kb = nelm(jj)
             if(kb.ne.i.and.kb.le.neq_primary) then
              cord1j = cord(kb,1)
              cord2j = cord(kb,2)
              cord3j = 0.0
              if(icnl.eq.0) then
               cord3j = cord(kb,3)
              endif 
c calculate distance i-kb  
              disx2 = (cord1-cord1j)**2
              disy2 = (cord2-cord2j)**2
              disz2 = (cord3-cord3j)**2
              dis2i = disx2 + disy2 + disz2
c find cosine squared              
              cosxi0 = disx2/dis2i
c find min cosine connection in + and - directions 
              if(cosxi0.gt.cosneqx+tol_cos
     &         .and.cord1j-cord1.gt.tol_cos) then
               cosneqx = cosxi0
               kb_pos = kb
              endif
              if(cosxi0.gt.cosposx+tol_cos
     &         .and.cord1-cord1j.gt.tol_cos) then
	         cosposx = cosxi0
	         kb_neg = kb
              endif
             endif
            enddo
            iconn_gdkm(i,1) = kb_pos
            iconn_gdkm(i,2) = kb_neg 
        if(igdpm(i).gt.0) then    
c gdkm nodes
c iconn_gdkm(3) has a different meaning   
       i_gdkm = nelm(nelm(i+1))-1      
       do kk = 1,ngdpm_layers(igdpm(i))
         i_gdkm = i_gdkm +1
         i1_gdkm = nelm(i_gdkm)+1
         i2_gdkm = nelm(i_gdkm+1)
         kb_neg = 0
         kb_pos = 0 
         disx_max = tol_cos
         disx_min = -tol_cos
         do ii =  i1_gdkm, i2_gdkm   
            kb = nelm(ii)
            cord1j = cord(kb,1)
            cord2j = cord(kb,2)
            cord3j = 0.0
              if(icnl.eq.0) then
               cord3j = cord(kb,3)
              endif    
c find min and max distances in x directions
              disx2 = cord1j-cord1
              if(disx2.gt.disx_max) then
               disx_max = disx2
               kb_pos = kb
              else if(disx2.lt.disx_min) then
               disx_min = disx2
               kb_neg = kb  
              endif      
         enddo    
       enddo
        iconn_gdkm(i,3) = kb_pos
        iconn_gdkm(i,4) = kb_neg 
       endif
       enddo
       endif
      else if(iflg.eq.2) then
c make and break connections 
       if(gdpm_flag.eq.11) then 
c model 11 has x direction nodes    
c start counter at number at total nodes   
c estimate size of new nelm array  
c assume all gdkm nodes are double sided    
       allocate(idump(neq))
       allocate(idumn(neq))
       idump = 0
       idumn = 0
       neqp1 = neq+1 
       ngdkm = neq-neq_primary
       nsize = nelm(neqp1)+5*ngdkm +ngdkm
       allocate(nelm_temp(nsize))
       nelm_temp(1) = neqp1
       allocate(istrw_temp(nsize-neqp1))
       allocate(idum1(neq))  
       idum1 = 0          
       allocate(sx_12(neq,2))
       if(isoy.eq.1) then
        allocate(sx_temp(nr+2*neq,1))
       else
        allocate(sx_temp(nr+2*neq,3))
       endif      
       ic = neq  
c save original connectivity (and stor file pointer)      
        do i = 1,neq_primary
         kbmin = neq
         kbmax = 0
         do jj = nelm(i)+1,nelm(i+1)
          kb = nelm(jj)
          kbmin = min(kbmin,kb)
          kbmax = max(kbmax,kb)
          idum1(kb) = istrw(jj-neqp1)
         enddo
         if(igdpm(i).ne.0) then
c kb_pos is global connection in  + x direction    
c kb_neg is global connection in  - x direction        
          kb_pos = iconn_gdkm(i,1)
          kb_neg = iconn_gdkm(i,2)         
          if(kb_pos.gt.0) then   
c since i is a gdpm node  find the +x and -x last gdpm nodes
c to these nodes, a new connection will be added
c i_pos is the last gdpm node attached to node i in the +x direction
c i_pos will connect to kb_pos (or a gdpm node attached to kb_pos)
           i_pos = iconn_gdkm(i,3)  
           cord_i_pos = cord(i_pos,1)
           area_i = areat_gdpm(i)
c connection to gdpm nodes should be last one: nelm(kb_pos+1)
c connetion in + x direction
c note different meaning of iconn_gdkm for gdpm nodes
          if(igdpm(kb_pos).ne.0) then
           kb_pos_neg = iconn_gdkm(nelm(kb_pos+1),2)
c find first gdkm node in -x direction
           if(kb_pos_neg.ne.0) then
c identify last gdkm node
            cord_kb_pos_neg = cord(kb_pos_neg,1)
c identify area of secondary material             
            area_kb_pos = areat_gdpm(kb_pos)
c primary connection is kb_pos
            idump(i_pos) = kb_pos_neg
            sx_12(i,1) = -((area_i+area_kb_pos)/2.)
     &         /(cord_kb_pos_neg-cord_i_pos)
           endif
          else
c kb_pos has no gdkm nodes- connect to primary nodel 
c remember the reciprocal relationship
           idump(i_pos) = kb_pos
           idumn(kb_pos) = i_pos
           sx_12(kb_pos,1) = -(area_i)
     &         /(cord(kb_pos,1)-cord_i_pos)                    
         endif    
         endif
         if(kb_neg.gt.0) then   
c since i is a gdpm node  find the +x and -x last gdpm nodes
           i_neg = iconn_gdkm(i,4)  
           cord_i_neg = cord(i_neg,1)
           area_i = areat_gdpm(i)    
         endif
        endif
       enddo
c now add new connection
c will leave connection from i to kb_pos (may need it later- for stress)
       ii = neqp1
       iww = nr
       if(isoy.eq.1) then
        sx_temp(1:nr,1) = sx(1:nr,1)
       else
        sx_temp(1:nr,1) = sx(1:nr,1)
        sx_temp(1:nr,2) = sx(1:nr,2)
        sx_temp(1:nr,3) = sx(1:nr,3)
       endif
c loop on all nodes  
       idum1 = 0
       do i = 1, neq
        kbmin = neq
        kbmax = 0
        i1 = nelm(i)+1
        i2 = nelm(i+1)
        do jj = i1,i2
         kb = nelm(jj)
         kbmin = min(kbmin,kb)
         kbmax = max(kbmax,kb)
         if(kb.ne.i) then
          idum1(kb) = istrw(jj-neqp1)
         else
          idum1(kb) = -1
         endif
        enddo
        if(i.gt.neq_primary.and.idump(i).ne.0) then
         iww = iww + 1
         sx_temp(iww,1) = sx_12(i,1)
         kb = idump(i)
         kbmin = min(kbmin,kb)
         kbmax = max(kbmax,kb)
         idum1(kb) = iww
        endif   
        if(i.lt.neq_primary.and.idumn(i).ne.0) then
         iww = iww + 1
         sx_temp(iww,1) = sx_12(i,1)         
         kb =  idumn(i)
         kbmin = min(kbmin,kb)
         kbmax = max(kbmax,kb)
         idum1(kb) = iww     
        endif             
        do kb = kbmin, kbmax
         iw = idum1(kb)
         idum1(kb) = 0
         if(iw.ne.0) then
          ii = ii +1
          nelm_temp(ii)= kb
          istrw_temp(ii-neqp1) = max(iw,0)
          if(kb.eq.i) nelmdg(i) = ii
         endif
        enddo
        nelm_temp(i+1) = ii
       enddo 
c copy nelm_temp to nelm  
       deallocate(sx,nelm,istrw) 
       if(isoy.eq.1)then
        allocate(sx(iww,1))
       else
        allocate(sx(iww,3))
       endif
       allocate(istrw(ii-neqp1))
       allocate(nelm(ii))
       if(isoy.eq.1) then
        sx(1:iww,1) = sx_temp(1:iww,1)
       else
        sx(1:nr,1) = sx_temp(1:nr,1)
        sx(1:nr,2) = sx_temp(1:nr,2)
        sx(1:nr,3) = sx_temp(1:nr,3)
        sx(nr+1:iww,1) = 3.*sx_temp(nr+1:iww,1)
        sx(nr+1:iww,2) = 0.0
        sx(nr+1:iww,2) = 0.0
       endif
       nelm(1:ii) = nelm_temp(1:ii)
       istrw(1:ii-neqp1) = istrw_temp(1:ii-neqp1)
       deallocate(idump,idumn,nelm_temp,istrw_temp,sx_temp)       
      neqp1=neq+1
c test of coefficients
       if(iparchek.eq.1) then
       write(ierr,*) 'fe coef. neighbors'
       write(ierr,*)
       do i = 1,neq
        write(ierr,*)' vol node ', i,' = ', sx1(i)
        i1 = nelm(i)+1
        i2 = nelmdg(i)-1
        do j = i1,i2
         kb = nelm(j)
         i3 = nelmdg(kb)+1
         i4 = nelm(kb+1)
         do jj = i3,i4
         if(nelm(jj).eq.i) then
          iw = istrw(jj-neqp1)
          if(iw.gt.0) then
           a11 = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
          else
           a11 = 0.0
          endif
          write(ierr,355)i,kb,iw,a11
          go to 899
         endif
        enddo
899     continue         
       enddo        
        i1 = nelmdg(i)+1
        i2 = nelm(i+1)        
        do j = i1,i2
         kb = nelm(j)
         iw = istrw(j-neqp1)
          if(iw.gt.0) then
           a11 = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
          else
           a11 = 0.0
          endif
         write(ierr,355)i,kb,iw,a11
        enddo        
       enddo   
355   format(i6,1x,i6,1x,i6,1x,f12.3)        
       stop
       endif   
       endif 
      else if(iflg.eq.3) then   
c combine intrinsic permeability with fracture volume
       do i = 1,neq_primary
        imodel= igdpm(i)
c first gdpm node is given last connection of primary node
c check this logic
        kb = nelm(nelm(i+1))
        if(imodel.ne.0) then
          sx1_primary = sx1(i)
	    vfrac = vfrac_primary(imodel)
c multiply pemeabilities by volume fractions	    
	     pnx(i) = pnx(i)*vfrac
	     pny(i) = pny(i)*vfrac
	     pnz(i) = pnz(i)*vfrac
c save gdkm volume fraction for node 
           gdkm_volume_fraction(i) = vfrac           
           sx1_total = sx1_primary/vfrac
           do jk = 1, ngdpm_layers(imodel)
            kc = kb+jk-1
            vfrac2 = sx1(kc)/sx1_total
	      pnx(kc) = pnx(kc)*vfrac2
	      pny(kc) = pny(kc)*vfrac2
	      pnz(kc) = pnz(kc)*vfrac2  
c save gdkm volume fraction for node 
            gdkm_volume_fraction(kc) = vfrac2
           enddo
        endif
       enddo      
      else if(iflg.eq.-1) then    
       deallocate(iconn_gdkm)
      endif
      return
      end
