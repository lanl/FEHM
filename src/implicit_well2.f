      subroutine implicit_well(iflg)
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!D1 
!D1  PURPOSE
!D1
!D1      Create a simple well 
!D1      Using mdnodes technology
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 14-Aug-08, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/implicit_well.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.2 Finite-Element Coefficient Generation
!D3
C***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 GAZ initial implementation 1-17-05
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
C***********************************************************************
C  WE want to implement a river or stream that is connected to the aquifer
C  this is done by defining a stream by adding nodes to the grid and connecting 
C  them to the nearest gridblock
C
C  gdpm and river macros are incompatible
C
C***********************************************************************
C
C     INPUT ARGUMENTS -
C coordinates and label of river sections
C cross sectional area of sections
C river bottom are of sections
C permeability od sections 
C head-impedance data 
C
C#######################################################################
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comrxni
      use comii
      use comxi
      use davidi
      use comriv
c     
      implicit none
      integer iflg,i,j,i1,i2,ik,izone,inode,ii,ij,kb,iafile,neqp1_old
      integer neqp1,iii,node,neq_new, id, nbc, kk, ic, icsx, nsize_sx
      integer imr,ipr,imrl,iprl,kb2,j1,j2,i3, ncont, ipiv, iskp1, iskp2
      integer neqpi, iexpfl, icoef,nr_old,kbmin,kbmax, npiv, iowell
      integer open_file, inwel,max_con,max_well_con,n_ncon_old
      integer maxiarea,iareap, jj, n_ncon, irnode, isnode, nseg,k,kc

      real*8  vol0, vol2, pi, pi2, disa, raddif,disr,rad2,amult,area_tol
      real*8  sx_max,area_dum,area_new,sx_dum,vol_well_dum,sx_well_seg
      real*8  seg,seg2,dis,xc,yc,zc,segx,segy,segz,segx2,segy2,segz2
	real*8  disx,disy,disz,sx_well_to_grid,sx_well,rad,pd
	real*8  delx,dely,delz,xi,yi,zi,term1,term2,radkb
	real*8  rthick,rw,r0,pnxa,pnya,aperm,peaceman_term,onethird

      real*8, allocatable  :: sx_new(:,:)
      real*8, allocatable  :: dum(:)
      real*8, allocatable  :: dum1(:,:)
      real*8, allocatable  :: seg_min(:)

      integer, allocatable :: idum(:)
      integer, allocatable :: idum1(:,:)
      integer, allocatable :: idum2(:)
      integer, allocatable :: ncon_new(:)
      integer, allocatable :: istrw_new(:)
      integer, allocatable :: n_seg_min(:)
      integer, allocatable :: in_seg_yes(:)
c     
      parameter (area_tol=1.e-10, pi= 3.1415927, pi2=6.2831850)
	parameter (max_con = 100,max_well_con= 10)
	parameter (onethird = 1.0/3.0)
c     

      logical null1
      logical  used
	character*80 dum_string

      integer i4,ic1,node_seg,iseg,ntot2,nodew
      allocate(seg_min(max_con))
      allocate(n_seg_min(max_con))
	allocate(in_seg_yes(max_con))
C     
C#######################################################################
C     
c     
c     return if no wells 
c
      if(nriver.eq.0) return
	if(iflg.eq.0) then
      else if(iflg.eq.1) then
c 
c find primary grid nodes associated with well intervals
c
       allocate(nwell2_prim(n_well_prod))
c       
c nwell2_prim is the primary grid node closest to the segment definition nodes 
c      
       do i = 1, n_well_prod
	   xc = coor_well2(i,1)
	   yc = coor_well2(i,2)
	   zc = coor_well2(i,3)
	   call near3(xc, yc, zc, nodew, 0)
	   nwell2_prim(i) = nodew
	 enddo
c 
c find number of intersection of wells
c

        allocate(nwell2_int(n_well_prod))
	  allocate(iwell2_int(n_well_prod,max_well_con))  
	  allocate(seg_cos(n_well_seg,3)) 
	  allocate(idir_seg(n_well_seg))
	  nwell2_int = 0
	  iwell2_int = 0  
        do i = 1, n_well_seg
	   j1 = iwell_seg(i,1)
         j2 = iwell_seg(i,2)
         nwell2_int(j1) = nwell2_int(j1) + 1
         nwell2_int(j2) = nwell2_int(j2) + 1
         iwell2_int(j1,nwell2_int(j1)) = i
         iwell2_int(j2,nwell2_int(j2)) = i
c find orientation of each segment         
         segx = (coor_well2(j1,1)-coor_well2(j2,1))
         segy = (coor_well2(j1,2)-coor_well2(j2,2))  
         segz = (coor_well2(j1,3)-coor_well2(j2,3))
	   segx2 = segx**2
	   segy2 = segy**2
	   segz2 = segz**2
	   dis = segx2 + segy2 + segz2
	   seg_cos(i,1) = segx2/dis
	   seg_cos(i,2) = segy2/dis
	   seg_cos(i,3) = segz2/dis
	   disa = max(seg_cos(i,1),seg_cos(i,2),seg_cos(i,3))
	   if(disa-seg_cos(i,1).lt.area_tol) then
	    idir_seg(i) = 1
	   else if(disa-seg_cos(i,2).lt.area_tol) then
	    idir_seg(i) = 2
	   else
	    idir_seg(i) = 3
	   endif
	  enddo
c
c label new nodes 
c
c seg is segment length 
c dis is desired node spacing
c
       allocate(new_node_well2(n_well_seg*max_seg_div))
       allocate(new_node_well2_segid(n_well_seg*max_seg_div))
       allocate(coor_new_well2(n_well_seg*max_seg_div,3))
	 allocate(neigh_well2_new(npoint_riv,max_con))
	 allocate(idum2(npoint_riv))
	 allocate(idum(npoint_riv))
       allocate(nelm_riv(nnelm_riv,2))
c idum2 is a counter used in adding connections	 
	 idum2 = 3
	 idum = 3
	 neigh_well2_new = 0
       iwell_prod = 0 
	 ntot2 = neq_primary
       ic1 = 0 
       nnelm_riv = 0
	 do i = 1, n_well_seg
	  j1 = iwell_seg(i,1)
c
c search for existing node with same(j1) position
c
        xc = coor_well2(j1,1)
        yc = coor_well2(j1,2)
        zc = coor_well2(j1,3)
        iskp1 = 0
        do ii = 1, ic1  
         segx = (coor_new_well2(ii,1)-xc)
         segy = (coor_new_well2(ii,2)-yc)  
         segz = (coor_new_well2(ii,3)-zc)
	   segx2 = segx**2
	   segy2 = segy**2
	   segz2 = segz**2
	   dis = segx2 + segy2 + segz2
	   if(dis.lt.area_tol) then
c found match - don't add node at beginning of segment
          iskp1 = 1
          idum2(ii) = idum2(ii)+1
          neigh_well2_new(ii,idum2(ii)) = ic1 + 1 + neq_primary  
          neigh_well2_new(ic1+1,2) = ii + neq_primary   
          go to 121 
	   endif
        enddo	  
121     j2 = iwell_seg(i,2)
c
c search for existing node with same(j1) position
c
        xc = coor_well2(j2,1)
        yc = coor_well2(j2,2)
        zc = coor_well2(j2,3)
        iskp2 = 0
        do jj = 1, ic1
         segx = (coor_new_well2(jj,1)-xc)
         segy = (coor_new_well2(jj,2)-yc)  
         segz = (coor_new_well2(jj,3)-zc)
	   segx2 = segx**2
	   segy2 = segy**2
	   segz2 = segz**2
	   dis = segx2 + segy2 + segz2
	   if(dis.lt.area_tol) then
c found match - don't add node at end of segment
          iskp2 = 1
          idum(jj) = idum(jj)+1  
          go to 122       
	   endif
        enddo	        
122     segx = (coor_well2(j2,1)-coor_well2(j1,1))
        segy = (coor_well2(j2,2)-coor_well2(j1,2))  
        segz = (coor_well2(j2,3)-coor_well2(j1,3))
	  segx2 = segx**2
	  segy2 = segy**2
	  segz2 = segz**2
	  seg = sqrt(segx2 + segy2 + segz2)
	  dis = well_dz(i)
	  node_seg = seg/dis + 1
c "1" insures spacing is le than desired spacing
	  disx = segx/(node_seg-1)
	  disy = segy/(node_seg-1)
	  disz = segz/(node_seg-1)
c
c  fill in additional nodes along well
c   
        if(iwell_prod(j1).eq.0.and.iskp1.eq.0) then   
         ic1 = ic1 +1   
	   iwell_prod(j1) = ic1 
	   new_node_well2(ic1) = ic1 + ntot2
c set porosity of wellbore node = 1.	   
	   ps(ic1 + ntot2) = 1.0
	   pnx(ic1 + ntot2) = 1.e-30
	   rad = well_rad(i)
	   pny(ic1 + ntot2) = rad*rad/8.0
	   coor_new_well2(ic1,1) = coor_well2(j1,1)
	   coor_new_well2(ic1,2) = coor_well2(j1,2)
	   coor_new_well2(ic1,3) = coor_well2(j1,3)
	   new_node_well2_segid(ic1) = i
         neigh_well2_new(ic1,2) = 0
	   neigh_well2_new(ic1,3) = ic1 + neq_primary
	   neigh_well2_new(ic1,4) = ic1+1 + neq_primary
c identify zone of node with seqment i	   	   
	   izonef(ic1+neq_primary) = izlabels(i)
	   nnelm_riv = nnelm_riv + 1
	   nelm_riv(nnelm_riv,1) = ic1 + neq_primary
	   nelm_riv(nnelm_riv,2) = ic1 + 1 + neq_primary
c	  else if(iwell_prod(j1).eq.0.and.iskp1.eq.1) then 
	  else if(iskp1.eq.1) then 
	    neigh_well2_new(ic1+1,2) = ii + neq_primary   
	    nnelm_riv = nnelm_riv + 1	    
	    nelm_riv(nnelm_riv,1) =  ii + neq_primary
	    nelm_riv(nnelm_riv,2) = ic1 + 1 + neq_primary
	  endif
        do j = 2,node_seg-1
	   ic1 = ic1 + 1
         new_node_well2(ic1) = ic1 + ntot2
c set porosity of wellbore node = 1.	   
	   ps(ic1 + ntot2) = 1.0
	   pnx(ic1 + ntot2) = 1.e-30
	   rad = well_rad(i)
	   pny(ic1 + ntot2) = rad*rad/8.0      
	   coor_new_well2(ic1,1) = coor_new_well2(ic1-1,1) + disx
	   coor_new_well2(ic1,2) = coor_new_well2(ic1-1,2) + disy
	   coor_new_well2(ic1,3) = coor_new_well2(ic1-1,3) + disz
	   new_node_well2_segid(ic1) = i
	   nnelm_riv = nnelm_riv + 1
	   nelm_riv(nnelm_riv,1) =  ic1 + neq_primary
	   nelm_riv(nnelm_riv,2) =  ic1 + 1 + neq_primary
	   if(iskp1.eq.0) then
          neigh_well2_new(ic1,2) = ic1-1 + neq_primary
         else
          iskp1 = 0
         endif 
	   neigh_well2_new(ic1,3) = ic1 + neq_primary

	 	 if(iskp2.eq.0.or.j.ne.node_seg-1) then
	    neigh_well2_new(ic1,4) = ic1+1 + neq_primary
         else 
          neigh_well2_new(ic1,4) = jj + neq_primary   
          neigh_well2_new(jj,idum(jj)) = ic1 + neq_primary
         endif 	   
c identify zone of node with seqment i	   
	   izonef(ic1+neq_primary) = izlabels(i)
	  enddo  
        if(iwell_prod(j2).eq.0.and.iskp2.eq.0) then   
	   ic1 = ic1 +1  
c	   nnelm_riv = nnelm_riv + 1
c	   nelm_riv(nnelm_riv,1) =  ic1 + neq_primary
c	   nelm_riv(nnelm_riv,2) =  ic1 + 1 + neq_primary	  
	   iwell_prod(j2) = ic1   
	   new_node_well2(ic1) = ic1 + ntot2
c set porosity of wellbore node = 1.	   
	   ps(ic1 + ntot2) = 1.0
	   pnx(ic1 + ntot2) = 1.e-30
	   rad = well_rad(i)
	   pny(ic1 + ntot2) = rad*rad/8.0	   
	   coor_new_well2(ic1,1) = coor_well2(j2,1)
	   coor_new_well2(ic1,2) = coor_well2(j2,2)
	   coor_new_well2(ic1,3) = coor_well2(j2,3)
	   new_node_well2_segid(ic1) = i	
         neigh_well2_new(ic1,2) = ic1-1 + neq_primary
	   neigh_well2_new(ic1,3) = ic1 + neq_primary
	   neigh_well2_new(ic1,4) = 0	  
c identify zone of node with seqment i	   	   
	   izonef(ic1+neq_primary) = izlabels(i)	  
	  else if(iwell_prod(j2).eq.0.and.iskp2.eq.1) then 
c	    nnelm_riv = nnelm_riv + 1
	    nelm_riv(nnelm_riv,1) =  ic1 + neq_primary
	    nelm_riv(nnelm_riv,2) =  jj + neq_primary	  
	    neigh_well2_new(ic1,4) = jj + neq_primary 
	    neigh_well2_new(jj,4) = jj + neq_primary         	         
	  endif  	   
       enddo

       deallocate(idum2,idum)
c
c resize arrays
c
       npoint_riv = ic1
       n = neq_primary + npoint_riv
       n0 = neq_primary + npoint_riv
       if(iout.ne.0) then
        write(iout,*)'n0 adjusted for river nodes to ', n0 
       endif
       if(iptty.ne.0) then
        write(iptty,*)'n0 adjusted for river nodes to ', n0 
       endif
       allocate (idum(ic1))
       idum = 0
	 idum(1:ic1) = new_node_well2(1:ic1)
	 deallocate(new_node_well2)
	 allocate(new_node_well2(ic1))
       new_node_well2(1:ic1) = idum(1:ic1)

	 idum(1:ic1) = new_node_well2_segid(1:ic1)
	 deallocate(new_node_well2_segid)
	 allocate(new_node_well2_segid(ic1))
       new_node_well2_segid(1:ic1) = idum(1:ic1)

	 allocate (dum1(ic1,3))
	 dum1(1:ic1,1:3) = coor_new_well2(1:ic1,1:3)
	 deallocate(coor_new_well2)
	 allocate(coor_new_well2(ic1,3))
       coor_new_well2(1:ic1,1:3) = dum1(1:ic1,1:3)
c
c find well neighbors in well grid
c

c
c find neighbors in primary grid
c need to divede areas and and volumes bsed on the number
c of interscting well nodes
c
        allocate(neigh_well2(neq_primary,max_con))
        allocate(neigh_well2_count(neq_primary))

        neigh_well2_count = 0
        neigh_well2 = 0
        do i = 1, ic1
	   xc = coor_new_well2(i,1)
	   yc = coor_new_well2(i,2)
	   zc = coor_new_well2(i,3)
	   call near3(xc, yc, zc, nodew, 0)
         neigh_well2_count(nodew) = neigh_well2_count(nodew)+1
         neigh_well2(nodew,neigh_well2_count(nodew)) = i
	   neigh_well2_new(i,1) = nodew
	  enddo
        nodes_well2_added = ic1
	  neq = neq_primary + nodes_well2_added
c	  
c fill in coordinates of new well nodes
c
        do i = 1,ic1
	   xc = coor_new_well2(i,1)
	   yc = coor_new_well2(i,2)
	   zc = coor_new_well2(i,3)  
	   cord(i+neq_primary,1)= xc   
	   cord(i+neq_primary,2)= yc 
	   cord(i+neq_primary,3)= zc   
         corz(i+neq_primary,1)= xc   
	   corz(i+neq_primary,2)= yc 
	   corz(i+neq_primary,3)= zc        	        
        enddo
c allocate space for well coefficients        
        allocate(sx_w(ic1,max_con))
        sx_w = 0.0
calculate  coeficients               
	 deallocate(idum,dum1)
c 
	else if(iflg.eq.2) then
c	
	else if(iflg.eq.3) then
c
c  modify volumes
c  this model assumes wellbore volumes are small compared to gridblock 
c  volumes
c
	else if(iflg.eq.4) then
c     modify connectivity array 
c     called from riverctr(iflag= 4)
c     riverctr called from anonp
c     estimate ncon from previous connectivity by added "nodes_well2_added"
c     old size + nodes_well2_added (for row starts)
c     + nodes_well2_added (connections added to primary nodes)+ indexing for rows
c     + 7*nodes_well2_added (conection in wellbore)
         
         neqp1_old = neq_primary + 1
	   neq_new = neq_primary + nodes_well2_added
c
c  well_rad(n_well_seg)-well_radius for each segment
c  well_dz(n_well_seg)- distance increment for gridblock along well segment
c
        vol_well_dum = 0.0d0
        do ii = 1,nodes_well2_added
c identify segment and radius        
         id = new_node_well2_segid(ii)
         rad = well_rad(id)
         area_dum= pi*rad**2
         i = ii + neq_primary
         sx1(i) = 0.0
         xc = cord(i,1)  
	   yc = cord(i,2)
	   zc = cord(i,3)  
c only consider well segments here	   
         sx_w(ii,1) = 0.0      
         do j = 2,4
          kb = neigh_well2_new(ii,j)          
          if(kb.ne.0.and.kb.ne.i) then
c identify segment and radius        
           id = new_node_well2_segid(kb-neq_primary)
           radkb = well_rad(id)
           area_new= pi*radkb**2          
           segx = (xc-cord(kb,1))**2  
	     segy = (yc-cord(kb,2))**2
	     segz = (zc-cord(kb,3))**2  
	     dis = sqrt(segx+segy+segz)    
           sx_w(ii,j) = -0.5*(area_dum+area_new)/dis
           sx1(i) = sx1(i) + area_dum*dis/2.
c sx_w(ii,1) will contain the surface area of the segment ii          
           sx_w(ii,1) = sx_w(ii,1)+pi2*(rad+radkb)/2.*dis
          endif
         enddo  
         vol_well_dum = vol_well_dum + sx1(i)
        enddo   
        if(iout.ne.0) then
          write(iout,*)'>>>>> wellbore volume ',vol_well_dum,' <<<<<<'
        endif    
        if(iptty.ne.0) then
          write(iptty,*)'>>>>> wellbore volume ',vol_well_dum,' <<<<<<'
        endif  
c	   
c
c this estimate for n_ncon is too large - can improve
c  
         n_ncon_old= nelm(neq_primary)  
         n_ncon = n_ncon_old + nodes_well2_added*7 
         neqp1 =  neq_new +1  
         allocate(istrw_new(n_ncon-neqp1))
c estimate for size of sx_new is n_ncon/2-neqp1
         nsize_sx = max(n_ncon-neqp1, nr-1)
         allocate(sx_new(nsize_sx,3)) 
         sx_new = 0.0        
         allocate(ncon_new(n_ncon))
         allocate(idum(2*neq_new))
         allocate(idum2(neq_primary))
         allocate(dum(n_ncon-neqp1))
         allocate(dum1(neq,3))
         dum = 0.0
         j=neqp1
         icoef = nr
         nr_old = nr-1 
         icsx = nr_old     
         ncon_new(1)=neqp1
         idum2(1:neq_primary)= nelmdg(1:neq_primary)
         deallocate (nelmdg)
         allocate (nelmdg(neq_new))
         ic1 = neqp1
	   idum = 0
	   dum1 = 0.0
c load sx_new with primary grid sx 	
         do i=1,nr_old
            sx_new(i,isox) = sx(i,isox)
            sx_new(i,isoy) = sx(i,isoy)
            sx_new(i,isoz) = sx(i,isoz)
         enddo  
c do loop i start         
         do i=1,neq_primary
c            idum(i)=1
            kbmin=neq_new
            kbmax=0
            kbmin=min(kbmin,i) 
            kbmax=max(kbmax,i)
            i1=nelm(i)+1     
            i2=nelm(i+1) 
c do loop jj start            
            do jj=i1,i2
              kb=nelm(jj)
              kbmin=min(kbmin,kb)
              kbmax=max(kbmax,kb)
              idum(kb)=1   
c upper diagonal only for sx coefficients                         
              if(kb.gt.i) then             
               iw = istrw(jj-neqp1_old)
               idum(neq_new+kb)=iw
               if(iw.eq.nr) idum(neq_new+kb) = -99999999
               dum1(kb,isox) = sx(iw,isox)
               dum1(kb,isoy) = sx(iw,isoy)
               dum1(kb,isoz) = sx(iw,isoz) 
              endif              
            enddo
c do loop jj end              
c
c  added additional connections to well from primary grid
c  zero out first (1.e-9)
    
           sx_well = sx1(i)**0.333
c need to find orientation and use cosines to figure delx delz            
c do loop j start            
           do j = 1, neigh_well2_count(i) 
	      kb = neigh_well2(i,j) + neq_primary
            id = new_node_well2_segid(kb-neq_primary)    	      
            kbmin=min(kbmin,kb)
            kbmax=max(kbmax,kb)
            idum(kb) = 1
            dum1(kb,1) = -sx_w(id,1)/sx_well
	     enddo
c do loop j end  	     
c
c do loop j start 
	     do j = kbmin,kbmax  
	      if(idum(j).ne.0) then 
	       ic1 = ic1 +1
	       ncon_new(ic1) = j
	       if(j.eq.i) nelmdg(i) = ic1
	       idum(j) = 0
	       if (j.le.neq_primary.and.j.ne.i) then
	         iw = idum(neq_new+j)
	         istrw_new(ic1-neqp1)	= iw
	         idum(j+neq_new) = 0
	       else if(j.gt.neq_primary) then
c always add here new coefficients  (j > i here)    
	        icsx = icsx + 1
	        istrw_new(ic1-neqp1) = icsx	
	        sx_new(icsx,1) = dum1(kb,1)
	        dum1(j,1) = 0.0             
	      endif
	     endif
	    enddo
c do loop j end  	    
           ncon_new(i+1) = ic1
	   enddo
c do loop i end	   
c now complete up to primary nodes
c    (with wellbore connections)    
c
c now add  new wellbore nodes and connections 
c need to add neq_primary to the wellbore nodes
c assume for start volumes are calculated from input data
c radius aand segment length calculated previously
c
c
c do loop ii start
c
c contribution from just 1 wellbore node to istrw_new
c
         do ii=1,nodes_well2_added
            i = ii + neq_primary
            kbmin = neq_new
            kbmax = 0
c do loop j start            
	     do j = 1,max_well_con
            kb = neigh_well2_new(ii,j) 
	      if(kb.ne.0) then
	       idum(kb) = 1
	       idum(neq_new+kb) = -kb
	       kbmin=min(kbmin,kb)
             kbmax=max(kbmax,kb)
 	      endif
	     enddo
c do loop j end
c do loop j start	     
	     do j = kbmin,kbmax
	      if(idum(j).ne.0) then
	       idum(j) = 0
	       ic1 = ic1 +1
	       ncon_new(ic1) = j	
	       if(j.gt.i) then   
	        icsx = icsx +1    
	        istrw_new(ic1-neqp1) = icsx
c 4th position is the wellbore connection with larger node number	        
	        sx_new(icsx,1) = sx_w(i-neq_primary,4)
	       endif
	       if(j.eq.i) nelmdg(i) = ic1
	      endif
	     enddo
c do loop j end	     
            ncon_new(i+1) = ic1	     
	    enddo
c do loop ii end	    
c
c     area coefficient calculation
c     now calculate the coefficients
c     deallocate and allocate some arrays
c     note:nelmdg already loaded
c     
         icoef = icsx +1
         j = ic1
         deallocate(nelm,istrw,istrw_itfc)
         allocate (nelm(j))
         allocate (istrw(j-neqp1))
         allocate (istrw_itfc(j-neqp1))
         do i=1,j
            nelm(i)=ncon_new(i) 
         enddo
c         
c     (RJP 01/08/07 added following)- GAZ obtained from river_ctr routine
c  -99999999  below indicates last position
         ncont=ncon_new(neq_new+1)
         nelmd = ncont
         do i=1,j-neqp1 
            if(istrw_new(i).eq.-99999999) istrw_new(i) = icoef+1
            istrw(i)=istrw_new(i)
         enddo
         deallocate(sx,sx_w)
         if(isoy.ne.1) then
            allocate (sx(icoef+1,3))
         else
            allocate (sx(icoef+1,1))
         endif
         
         
         do i = 1,icoef
            if(isoy.ne.1) then
               sx(i,1) = sx_new(i,1)            
               sx(i,2) = sx_new(i,2)           
               sx(i,3) = sx_new(i,3)         
            else
               sx(i,1) = sx_new(i,1)             
            endif
         enddo
 
         istrw_itfc  = 0
         
         nr = icoef + 1
c gaz 060209 be careful         
c         neq_primary = neq
c
c adjust area coeffcients for wellbore nodes
c
c
c Peaceman model 
c
         allocate(delxw2(nodes_well2_added))
         allocate(delyw2(nodes_well2_added))
         allocate(delzw2(nodes_well2_added))
         allocate(iwdum(nodes_well2_added))
         do ii = 1,nodes_well2_added
c  i is primary node           
            i = neigh_well2_new(ii,1)
c  k is the well node            
            kc = ii + neq_primary
             delx = 0.0
             dely = 0.0
             delz = 0.0
c find gridblock dimensions of primary grid that holds well node
             xi = cord(i,1)
             yi = cord(i,2)
             zi = cord(i,3)
             i1 = nelm(i)+1
             i2 = nelm(i+1)
             do k = i1,i2
              kb = nelm(k)
              if(kb.eq.kc) iwdum(ii) = k
              delx = max(abs(xi-cord(kb,1)),delx)
              dely = max(abs(yi-cord(kb,2)),dely)
              delz = max(abs(zi-cord(kb,3)),delz)
             enddo
            delxw2(ii) = delx
            delyw2(ii) = dely
            delzw2(ii) = delz  
          enddo
c          
         do ii = 1,nodes_well2_added
c  i is primary node           
            i = neigh_well2_new(ii,1)
c  k is the well node            
            kc = ii + neq_primary
             delx = delxw2(ii)
             dely = delyw2(ii)
             delz = delzw2(ii)
c find gridblock dimensions  
c jj is the connection primary to wellbore
            jj = iwdum(ii)           
            id = new_node_well2_segid(ii)
	      iw = istrw(jj-neqp1)
	      delx = delxw2(ii)
            dely = delyw2(ii)
            delz = delzw2(ii)
	      if(idir_seg(id).eq.1) then
c   aligned with x axis
            else if(idir_seg(id).eq.2) then
c   aligned with y axis                    
            else 
c   aligned with z axis               
               rthick = delz
               rw = well_rad(id)
   	         r0 = 0.14*sqrt(delx*delx+dely*dely)
   	         term2 = pi2*rthick/log(r0/rw)   
              if(isoy.eq.1) then 	         
   	          sx(iw,1) = -onethird*term2
   	         else
   	          sx(iw,1) = -term2
   	          sx(iw,2) = 0.0
   	          sx(iw,3) = 0.0
   	         endif   	         	            
            endif
         enddo
c
c set peaceman type connection to primary grid block
c only for nodes that represent segment ends
c     
       do ii = 1, n_well_prod
       	  i = iwell_end(ii)
       	  j = i
       	  kk = i - neq_primary     	  
201     nodew = nwell2_prim(ii)
             delx = delxw2(ii)
             dely = delyw2(ii)
             delz = delzw2(ii)
c find gridblock dimensions  
          jj = iwdum(kk)           
          id = new_node_well2_segid(i-neq_primary)
	    iw = istrw(jj-neqp1)
	    if(idir_seg(id).eq.1) then
c   aligned with x axis
          else if(idir_seg(id).eq.2) then
c   aligned with y axis                    
          else
c   aligned with z axis   
               rthick = delz
               rw = well_rad(id)
   	         pnxa = 1.e-6*pnx(nodew)
   	         pnya = 1.e-6*pny(nodew)
   	         aperm = sqrt(pnxa*pnya)
   	         term1 = (pnya/pnxa)**0.25 + (pnxa/pnya)**0.25
   	   term2 = sqrt(sqrt(pnya/pnxa)*delx*delx + sqrt(pnxa/pnya)*
     &	               dely*dely) 
    	         r0 = 0.28*term2/term1	         
   	         aperm = sqrt(pnxa*pnya)
   	         peaceman_term = 1.e06*pi2*aperm*rthick/log(r0/rw)   	
               if(isoy.eq.1) then 	         
   	          sx(iw,1) = -onethird*peaceman_term
   	         else
   	          sx(iw,1) = -peaceman_term
   	          sx(iw,2) = 0.0
   	          sx(iw,3) = 0.0
   	         endif
   	      endif
c for a producer set perm x (connection to primary grid) = 1.   	         
   	      pnx(i) = well_connect_fac(ii)   	  
	 enddo
         deallocate(idum,dum,istrw_new,ncon_new,sx_new)  
         deallocate(seg_cos,idir_seg,iwdum,iwell_end)
         deallocate(delxw2,delyw2,delzw2)   
	else if(iflg.eq.5) then
c
c assign zone number to potential source/sink nodes   
c     
       allocate(iwell_end(n_well_prod))
       do ii = 1, n_well_prod
	   xc = coor_well2(ii,1)
	   yc = coor_well2(ii,2)
	   zc = coor_well2(ii,3)
	   do j = neq_primary+1,neq
	    dis = abs(cord(j,1)-xc)+abs(cord(j,2)-yc)+abs(cord(j,3)-zc)
	    if(dis.lt.area_tol) then
	     izonef(j) = izlabelp(ii)
	     iwell_end(ii) = j
	    endif
	   enddo
	 enddo          
	else if(iflg.eq.-33) then
c
c    
c
       continue
	else if(iflg.eq.33) then
c
c    
c
       continue
	else if(iflg.eq.6) then
c
c organize output 
c  
c for now, called when contour plot is called 
c
      allocate(dum(n_well_seg))
      dum = 0.0
      inquire(unit = nufilb(30), opened = used)
      if(.not.used) then
       open(unit = nufilb(30), file = nmfil(30),form = cform(30))
      endif
      iowell = nufilb(30)
      write(iowell,*)
      write(iowell,*) 'Output for well module; time(days) = ', days
      write(iowell,500) 
500   format(/,1x,'well node ',' pg node ',' well id  ',5x,'p or h ',7x,
     &'t ',5x,'src ',7x,'src e ',10x,'x coor ',6x,'y coor ',3x,'z coor')
      jj = neq-nodes_well2_added        
      do ii = 1,nodes_well2_added
       i = ii + jj
       id = new_node_well2_segid(ii)
       nodew = neigh_well2_new(ii,1)
       xc = cord(i,1)
       yc = cord(i,2)
       zc = cord(i,3)
       dum(id) = dum(id) + sk(i)
       if(ihead.eq.0) then
        pd=phi(i)
       else
        call headctr(4,i,phi(i),pd)
       endif
       write(iowell,501) i,nodew,id,pd,t(i),sk(i),qh(i),xc,yc,zc
      enddo
501   format(1x,3(1x,i8),1x,f12.4,1x,f8.3,1p,2(1x,g12.4),0p,2(1x,f12.2),
     & f10.3)  
c 
c print out net mass flow
c   
      write(iowell,*) 
      write(iowell,*)'Net source/sink for segment(seg id net massflow)'
      write(iowell,502) (i,dum(i), i = 1,n_well_seg)  
502   format(1x,1p,5('(',i8,1x,g12.5,')'))              
      endif
      return
      end
