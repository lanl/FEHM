      subroutine gdkm_calc(iflg) 
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
!D1 Calculate areas and lengths for GDKM calcs 
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
      use davidi
      implicit none

      integer i,icode,iflg,neqp1,i1,i2,jj,kb,kc,ii,kb1,i3s,i4s,kks
      integer i3,i4,kk,idir,imodel,j,nsize,ngdkm_con,ic,iws,kbs
      integer jk,kbmodel,idp,kbm,kcm,jm,minkb,maxkb,i1s,i2s,jjs
      integer icw,i1_mem
      integer, allocatable :: istrw_gdkm(:)  
      integer, allocatable :: idum(:) 
      integer, allocatable :: nelm_temp(:)
      integer, allocatable :: nconi(:) 
      integer, allocatable :: nneigh(:,:) 
      real*8,  allocatable :: vf_gdpm(:,:)
      real*8,  allocatable :: sxdum(:)
      real*8,  allocatable :: sx_gdkm(:,:)
      real*8,  allocatable :: sx_dum(:,:)
      real*8 area_mult,areat,sx2c,dis,dis2,cosz,rlp_min
      real*8 dis2min,aream
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 cord1jg,cord2jg,cord3jg
      real*8 fac_nop
      real*8 sx1_primary,vfrac,sx1_total,vol_frac,perm
      parameter (rlp_min = 1.d-2)
      parameter (ngdkm_con = 6)
     
      neqp1 = neq +1
      if(gdkm_flag.eq.0) return
      if(iflg.eq.1) then
c     calculate some areas and lengths for GDKM calculations 
c     establish connectivity for gdkm calcs
c     need estimate for  size of nelm_gdkm
c     at this call nelm should be populated gdpm connectivity
         allocate (idum(neq))
         allocate (sxdum(neq))
         i1 = (neq-neq_primary)*ngdkm_con
         i2 = neq-neq_primary
         allocate (sx_gdkm(i2,27))
         allocate (nconi(i2))
         allocate (nneigh(i2,27))
         sx_gdkm = 0.0d0
         sxdum = 0.0d0
         idum = 0
         nconi = 0
         nneigh = 0
c     ic is a counter for forming the connectivity matrix       
         ic = i2
c     calculate volume_fractions of all gdpm nodes 
         allocate(vf_gdpm(neq_primary,maxgdpmlayers))
         do i = 1,neq_primary
            imodel= igdpm(i)
c     first gdpm node is given last connection of primary node
c     check this logic
            kb = nelm(nelm(i+1))
            if(imodel.ne.0) then
               sx1_primary = sx1(i)
               vfrac = vfrac_primary(imodel)
               sx1_total = sx1_primary/vfrac
               do jk = 1, ngdpm_layers(imodel)
                  vf_gdpm(i,jk) = sx1(kb+jk-1)/sx1_total
               enddo
            endif
         enddo
c     now finish calculation   
c     this logic is only strictly valid for parallel fracture model
         do i = 1,neq_primary
            imodel= igdpm(i)
            if(imodel.ne.0) then
               cord1 = cord(i,1)
               cord2 = cord(i,2)
               cord3 = 0.0
               if(icnl.eq.0) then
                  cord3 = cord(i,3)
               endif
c     find areas of all connecting (primary) gridblocks
c     parts lt i and gt i
c     first lt i (if i2<i1 then no nodes lt)
               i1 = nelm(i)+1
               i2 = nelmdg(i)-1
               j = 0
               do jj = i1,i2
                  kb = nelm(jj)
                  i3 = nelmdg(kb)+1
                  i4 = nelm(kb+1)
                  do kk = i3,i4
                     kb1 = nelm(kk)
                     if(kb1.eq.i) then
                        iw = istrw(kk-neqp1)
                        go to 100
                     endif
                  enddo
 100              continue           
                  j = j+1     
                  cord1j = cord(kb,1)
                  cord2j = cord(kb,2)
                  cord3j = 0.0
                  if(icnl.eq.0) then
                     cord3j = cord(kb,3)
                  endif 
                  sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                  dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &                 (cord3-cord3j)**2)      
                  it10(j) = kb
                  t1(j) = sx2c
               enddo  
c     next gt i (but only primary nodes; hence the "-1"
               i1 = nelmdg(i)+1
               i2 = nelm(i+1)-1
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
                  sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                  dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2 +
     &                 (cord3-cord3j)**2)      
                  it10(j) = kb
                  t1(j) = sx2c
               enddo  
c     first find total volume of primary and secondary porosity of node i
               sx1_primary = sx1(i)
               vfrac = vfrac_primary(imodel)
               sx1_total = sx1_primary/vfrac
c     look at secondary porosity nodes of node i  
c     examine connections to neighbors (idp 1st first gdpm or secondary node)
               idp = nelm(nelm(i+1))   
               do jk = 1, ngdpm_layers(imodel) 
c     ii is the secondary porosity node number          
                  ii = idp +(jk-1) 
                  maxkb = 0
                  minkb = neq 
                  cord1 = cord(ii,1)
                  cord2 = cord(ii,2)
                  cord3 = 0.0
                  if(icnl.eq.0) then
                     cord3 = cord(ii,3)
                  endif
c     remember volume frac of secondary node           
                  vol_frac = sx1(ii)/sx1_total
c     save all areas and nodes in secondary porosity node ii 
                  idum(ii) = 1  
c     neighbor nodes greater than ii
                  i1s = nelmdg(ii)+1
                  i2s = nelm(ii+1)
                  do jjs = i1s,i2s 
                     kbs = nelm(jjs)
                     iws = istrw(jjs-neqp1)
                     sx2c=sx(iws,isox)+sx(iws,isoy)+sx(iws,isoz)
                     idum(kbs) = 1
                     sxdum(kbs) = sx2c
                  enddo 
c     neighbor nodes less than ii            
                  i1s = nelm(ii)+1
                  i2s = nelmdg(ii)-1
                  do jjs = i1s,i2s 
                     kbs = nelm(jjs)
                     i3s = nelmdg(kbs)+1
                     i4s = nelm(kbs+1) 
                     do kks = i3s,i4s
                        if(nelm(kks).eq.ii) then
                           iws = istrw(kks-neqp1)
                           sx2c=sx(iws,isox)+sx(iws,isoy)+sx(iws,isoz)
                           idum(kbs) = 1
                           sxdum(kbs) = sx2c
                           maxkb = max(kbs,maxkb)
                           minkb = min(kbs,minkb)
                           go to 50
                        endif
                     enddo
 50                  continue            
                  enddo
c     
c     at this point we have, for every primary gridblock, the primary connections
c     and the areas, these are in it10(j) and t1(j) respectively
c     we also have the 
c     
                  do jj = 1,j
c     loop to find other primary node connections (with gdpm nodes)   
                     kb = it10(jj)
c     areat = t1(jj)  this is the A/d formulation   
                     areat = t1(jj)
                     cord1j = cord(kb,1)
                     cord2j = cord(kb,2)
                     cord3j = 0.0
                     if(icnl.eq.0) then
                        cord3j = cord(kb,3)
                     endif
c     find node kb matrix nodes   
                     kbmodel= igdpm(kb)
                     if(kbmodel.ne.0) then
c     identify first matrix node            
                        kbm = nelm(nelm(kb+1))
                        dis2min = 1.d20
                        do jm = 1, ngdpm_layers(kbmodel)
c     note special numbering for gdpm nodes  
                           kc =  kbm+jm-1
                           cord1jg = cord(kc,1)
                           cord2jg = cord(kc,2)
                           cord3jg = 0.0
                           if(icnl.eq.0) then
                              cord3jg = cord(kc,3)
                           endif
c     now find closest node kb gdpm node to i gdpm node  
c     can work with square of distance
                           dis2 = (cord1-cord1jg)**2 + 
     &                          (cord2-cord2jg)**2 + (cord3-cord3jg)**2
                           if(dis2.lt.dis2min) then
                              kcm = kc
                              dis2min = dis2
                              aream = areat
                           endif                  
                        enddo
                        maxkb = max(kcm,maxkb)
                        minkb = min(kcm,minkb)
                        idum(kcm) = idum(kcm) + 1
c     calculated area/dis term             
                        sxdum(kcm) = aream
                     else
c     connect to primary node  
                        idum(kb) = 1
                        dis = sqrt((cord1-cord1j)**2 + (cord2-cord2j)**2
     &                       + (cord3-cord3j)**2)    
                        sxdum(kb) = areat
                     endif
                  enddo
c     add connection to secondary node (ii) associated with primary node i
                  idum(ii) = 1
                  maxkb = max(ii,maxkb)
                  minkb = min(ii,minkb)
c     adjust connectivity and area coefficients           
                  do jj = minkb,maxkb
                     if(idum(jj).ne.0) then
                        nconi(ii-neq_primary) = nconi(ii-neq_primary) +1
                        nneigh(ii-neq_primary,nconi(ii-neq_primary)) = 
     &                       jj
                        sx_gdkm(ii-neq_primary,nconi(ii-neq_primary)) =
     &                       sxdum(jj)
                        idum(jj) = 0
                        sxdum(jj) = 0.0
                     endif
                  enddo
c     
               enddo
            endif
         enddo
      endif
c     insure connection symmetry    
      do i = neq_primary+1, neq
         i1 = nconi(i-neq_primary)
         do jj = 1,i1
c     identify i-kb connection       
            kb = nneigh(i-neq_primary,jj)      
c     look for kb-i connection 
            if(kb.gt.neq_primary.and.kb.ne.i) then 
               i2 = nconi(kb-neq_primary)     
               do kk = 1,i2
                  kc = nneigh(kb-neq_primary,kk)
                  if(kc.eq.i) then
c     connection found 
                     go to 150        
                  endif
               enddo
c     no connection found (add it)
               nconi(kb-neq_primary) = nconi(kb-neq_primary) + 1   
               nneigh(kb-neq_primary,nconi(kb-neq_primary)) = i  
               sx_gdkm(kb-neq_primary,nconi(kb-neq_primary)) = 1.d20
 150           continue 
            endif       
         enddo
      enddo 
c     now all the connections are included
c     complete area terms (always choose smallest A/d)
      do i = neq_primary+1, neq
         i1 = nconi(i-neq_primary)
         do jj = 1,i1
c     identify i-kb connection       
            kb = nneigh(i-neq_primary,jj)      
c     look for kb-i connection 
            areat = abs(sx_gdkm(i-neq_primary,jj)) 
            if(kb.gt.neq_primary.and.kb.ne.i) then
               i2 = nconi(kb-neq_primary)     
               do kk = 1,i2
                  kc = nneigh(kb-neq_primary,kk)
                  if(kc.eq.i) then
c     connection found 
                     aream = abs(sx_gdkm(kb-neq_primary,kk))
                     go to 250        
                  endif
               enddo
c     error condition
               write (ierr,*) 'missing connection i = ', i, 'kb = ',kb
               write (ierr,*) 'stopping in gdkm_calc'
               stop
 250           continue
               aream = -min(aream,areat)
               sx_gdkm(kb-neq_primary,kk) = aream
               sx_gdkm(i-neq_primary,jj) = aream    
            endif             
         enddo
      enddo
c     now add up extra space required for nelm
      ic = 0
      do i = neq_primary+1, neq
         ic = ic + nconi(i-neq_primary)
      enddo  
c     correct size for nelm is nsize + ic
      nsize = nelm(neq+1)
c     this is the gdpm count:right total nodes and one connection from secondary nodes  
c     ic contains secodary node connection back to primary plus connections 2nd to 2nd
      i2 = nsize + ic 
      allocate(nelm_temp(nsize))
      nelm_temp(1:nsize) = nelm(1:nsize) 
      deallocate (nelm)
      allocate (nelm(i2)) 
      ic = neq+1
      nelm(1) = ic
      do i = 1,neq_primary
         i1 = nelm_temp(i)+1
         i2 = nelm_temp(i+1)
         do j = i1,i2
            nelm(j) = nelm_temp(j)
            ic = ic +1
         enddo
         nelm(i+1) = ic
      enddo
c     
      do i =  neq_primary+1, neq 
         do j = 1, nconi(i-neq_primary)
            ic = ic + 1
            kb = nneigh(i-neq_primary,j)
            nelm(ic) = kb
            if(kb.eq.i) nelmdg(i) = ic
         enddo
         nelm(i+1) = ic
      enddo       
c     construct new coefficient file
c     first calculate the new size (icw) 
c     adjust for symmetry and no diag
c     
      icw = 0
      do i = 1, neq
         i1 = nelmdg(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
            icw = icw+1
         enddo
      enddo
      nsize = nelm(neq+1)-(neq+1)
c     
      allocate(istrw_gdkm(nsize))
      istrw_gdkm = 0
      if(isoy.eq.1) then
         allocate(sx_dum(icw,1)) 
         icw = 0
         do i = 1, neq_primary
            i1 = nelmdg(i)+1
            i2 = nelm(i+1)
            do jj = i1,i2
               iw = istrw(jj-neqp1)
               icw = icw+1
               istrw_gdkm(jj-(neq+1)) = icw
               if(iw.ne.0) then
                  sx2c = sx(iw,1)
                  sx_dum(icw,1) = sx2c
               else
                  sx_dum(icw,1) = 0.
               endif
            enddo
         enddo 
      else
         allocate(sx_dum(icw,3))
         icw = 0
         do i = 1, neq_primary
            i1 = nelmdg(i)+1
            i2 = nelm(i+1)
            do jj = i1,i2
               iw = istrw(jj-neqp1)
               sx2c = sx(iw,1)
               icw = icw+1
               istrw_gdkm(jj-(neq+1)) = icw
               sx_dum(icw,1) = sx(iw,1)
               sx_dum(icw,2) = sx(iw,2)
               sx_dum(icw,3) = sx(iw,3)
            enddo
         enddo 
      endif 
      deallocate (sx)

c     fill in symmetric coefficient matrix
      deallocate(istrw)
      ic = icw
c     for connection with secondary nodes
c     insure that the newly calculated coefficients are used 
c     so must search below the diagonal
      sxdum = 0.0d0
      do i = neq_primary+1, neq  
         do j = 1,nconi(i-neq_primary)
            sxdum(nneigh(i-neq_primary,j)) = sx_gdkm(i-neq_primary,j)
         enddo
         i1 = nelmdg(i)+1
         i2 = nelm(i+1)
         do jj = i1,i2
            kb = nelm(jj)
            ic = ic +1
            istrw_gdkm(jj-neqp1) = ic
            if(isoy.eq.1) then
               sx_dum(ic,isox) = sxdum(kb)/3.
            else
               sx_dum(ic,isox) = sxdum(kb)
               sx_dum(ic,isoy) = 0.0
               sx_dum(ic,isoz) = 0.0
            endif
         enddo
c     adjust primary to secondary connection 
c     primary node must be first in connectivity list for node i      
         i1 = nelm(i)+1
         kb = nelm(i1)
c     secondary node is last in neighbor list for primary node        
         kk = nelm(kb+1)  
         iw = istrw_gdkm(kk-neqp1) 
         if(isoy.eq.1) then
            sx_dum(iw,1) = sxdum(kb)/3.
         else
            sx_dum(iw,isox) = sxdum(kb)
            sx_dum(iw,isoy) = 0.0
            sx_dum(iw,isoz) = 0.0          
         endif          
      enddo   
      allocate(istrw(nsize))
      istrw = istrw_gdkm
      deallocate(istrw_gdkm)
c     
c     change definition of nr
      nr = ic + 1
c     
      if(isoy.eq.1) then
         allocate(sx(nr,1)) 
         do j = 1, ic
            sx(j,1) = sx_dum(j,1)
         enddo 
      else
         allocate(sx(nr,3))
         do j = 1, ic
            sx(j,1) = sx_dum(j,1) 
            sx(j,2) = sx_dum(j,2)
            sx(j,3) = sx_dum(j,3)
         enddo
      endif    
      if(isoy.eq.1) then
         sx(nr,1) = 0.0
      else
         sx(nr,1) = 0.0
         sx(nr,2) = 0.0
         sx(nr,3) = 0.0
      endif

      do i = 1,neq_primary
         imodel= igdpm(i)
c     first gdkm node is given last connection of primary node
         if(imodel.ne.0) then
            kb = nelm(nelm(i+1))
            corz(kb,1) = cord(i,1)
            corz(kb,2) = cord(i,2)
            corz(kb,3) = cord(i,3)
         endif
      enddo
      gdpm_flag = 0 
      neq_gdkm = neq_primary
      neq_primary = neq              
c     deallocate memory      
      deallocate(idum,nconi,nneigh)
      deallocate(vf_gdpm,sxdum,sx_dum,sx_gdkm)
      if(allocated(vfrac_primary)) deallocate(vfrac_primary)
c     endif 
c     
      end
