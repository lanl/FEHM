      subroutine river_ctr(iflg)
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
!D1      Set up surface flow (rivers and streams)
!D1      Using mdnodes technology
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 17-Jan-05, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/river_ctr.f_a  $
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
C  10-26-05 GAZ
C  n_ncon needs a 2*ii which I don't understand (OK but watch) 
C  10-27-05 GAZ
C  cludged well overlying nodes 
C  management of multiple segments needs work
C  node gets out of bounds for nwell 
C  set istrw_itfc = 0 (probably negates this capability)
C question of what coorinates to use for radial elements 
C 11-05-2005
C anisotropic perm can get off center and weight law perm inexpectedly
C 11-11-05  
C divided sx(ik,1) by 3 for isotropic calculation in geneq1 etc.
C 11-16-05
C made connection from primary node to well mult_well*largest connection
C
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
      integer neqp1,iii,node,neq_new, id, nbc, kk
      integer imr,ipr,imrl,iprl,kb2,j1,j2,i3, ncont
      integer isection,nmodels,ic,mb,mc,k,jk,kbmax,kbmin,icoef,nr_old
      real*8 x1,x2,y1,y2,area_total,perm,xriv,yriv,zriv,area_tol, area0
      real*8  rad0, rad1, rad01, rad12, dis0, disx, disy, disz, radmax
      real*8  vol0, vol2, pi, pi2, disa, raddif,disr,rad2,amult,rad02
      real*8  mult_well,sx_max,area_dum,sx_dum,outrad,outarea,ht
      real*8 area_new, rads
      character*10 macro1
      character*4 macro
      integer i4
      integer neqpi, iexpfl
      real*8, allocatable  :: sx_new(:,:)
      real*8, allocatable  :: dum(:)
      real*8, allocatable  :: dum1(:)
      real*8, allocatable  :: cv_rad(:)

      integer, allocatable :: idum(:)
      integer, allocatable :: idum1(:)
      integer, allocatable :: ncon_new(:)
      integer, allocatable :: istrw_new(:)
      integer maxiarea, iareap, jj, n_ncon, irnode, isnode
c     
      integer neq_riv, maxiriver, maxc, istart
c     RJP 1/9/07 changed following	     
      parameter (maxiriver=100000, maxc=5)
      parameter (area_tol=1.e-10, pi= 3.1415927, pi2=6.2831850)
c     
      parameter (mult_well=1000.0d0)
      integer open_file, inwel
      logical null1
c     RJP 1/9/07 added river_zone
      integer iroot, river_zone(maxiriver), icdum
      character*120 well_file_name, well_root
	character*80 dum_string
c RJP 05/21/08 added following
      integer ic1
c     RJP 1/9/07 changed following
      save n_ncon, river_zone, iexpfl, neq_new, well_root, iroot, inwel
           
      real*8 well_connect_tol
      parameter(well_connect_tol = 1.d-30)
      
      maxlay = 300

c     
C     
C#######################################################################
C     
C     
c     return if no rivers or wells 
      if(nriver.eq.0) return
      if(iflg.eq.0) then
c     
c     read input 
c     izone_area: zone on which to calculate areas etc for this macro
c     iriver = type  of surface flow
c     iriver=0 : no routing, ponding, groundwater connection
c     iriver=1 : simple fluid routing, no groundwater connection
c     iriver=2 : simple stream definition, no groundwater connection
c     iriver=3 : simple stream definition, with groundwater connection
c     
c     iriver=4 : simple stream definition, with groundwater connection (with ponding)
c     iriver=5 : 
c     iriver=6 : 
c     iriver=7 : 
c     iriver=8 : 
c     
c     iriver01 = type of connection with other river segments
c     iriver02 : simple fluid routing, no groundwater connection
c     iriver=2 : simple stream definition, no groundwater connection
c     iriver=3 : simple stream definition, with groundwater connection
c     
c     wgt_river= river or river weight
c     
c     river_nodes_local(j) local river number of global gridblock j (size = all nodes)
c     river_nodes_global(j) global gridblock number of local river number j (size = river nodes)
c     
c     RJP 11/21/05 Added output file
c     zvd 06/14/07 Added root_name selection
         if (null1(root_name)) then
            if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ')
     &           then
               call file_prefix(nmfil(5), iroot)
               if (iroot .gt. 100) iroot = 100
               well_root(1:iroot) = nmfil(5)(1:iroot)
            else 
               if (nmfil(2)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil2 file: ', nmfil(2),
     &                 ' unable to determine well file prefix'
                  stop
               else
                  call file_prefix(nmfil(2), iroot)
                  if (iroot .gt. 100) iroot = 100
                  well_root(1:iroot) = nmfil(2)(1:iroot)
               end if
            endif
         else
            iroot = len_trim (root_name)
            if (iroot .gt. 100) iroot = 100
            well_root(1:iroot) = root_name(1:iroot)
         end if
c     open(1001, file = well_file_name, status='unknown',
c     *           form = 'formatted')

c         if(iriver.eq.0) then
c     
c     iflg = 0 and iflg = 1 called from incoord
c     at this point neq_primary has been set in scanin
c     n0 = neq_primary + npoint_riv
c     npoint_riv has been estimated in scanin
c     

c            isriver = 0
            nbc_riv = 0
            iriver = 1
            macro = "well"
c
c read in "sub macros" for well
      
100     continue
        read (inpt, '(a80)') wdd1
        if (wdd1(1:9) .eq. 'wellend') go to 200 
        if (wdd1(1:9) .eq. '       ') go to 200 
        if (wdd1(1:1) .eq. '#') go to 40 
        read (wdd1, '(a10)') macro1
     	  write(iout, 50) macro1
	  if (iptty .gt. 0) write(iptty, 50) macro1 
50         format(3x, '**** well sub macro : ', a10,' **** ' ) 
        if(macro1.eq.'text      ') then
c     
c     read text here     
c
        else if(macro1.eq.'wellmodel ') then           

c     nriver is the number of models for this call (all have to be of the same type)
            read(inpt,'(a80)') dum_string
	      read(dum_string,*,end = 140) nriver, iriver
	      go to 141
140         iriver = 1
141         continue
            allocate(iwell_geom(nriver))
            do i = 1,nriver
	       iwell_geom(i) = iriver
	      enddo
            nall = n0

            if(.not.allocated(izone_river)) then
               allocate(izone_river(max(1,nriver)))
               allocate(iriverf(max(1,nriver)))
            end if
c         endif
         nmodels = 0
c     
c     "i" is the river or well number
c 
       if(iriver.eq.1) then 
c	
          npoint_riv = 0      
          if(.not.allocated(isriverf)) then   
            allocate(isriverf(maxiriver))
            allocate(ifriverf(maxiriver))
            allocate(inriverf(maxiriver))
            allocate(iwsp(maxiriver))
            allocate(bc_river(maxiriver,2))
            allocate(ibc_river(maxiriver))
            allocate(rivbegin(maxiriver))
            allocate(rivend(maxiriver))
c     allocate(wgt_river(n0))
            allocate(isriver(maxriver_nodes))
            isriver = 0
            allocate(coor_riv(maxriver_nodes,3))
            allocate(area_riv(maxriver_nodes,2))
c     RJP 07/20/07 added a new array to specify casing & cement thickness
            allocate(csg_th(maxriver_nodes,2))
            allocate(perm_riv(maxriver_nodes,2))
          endif
c     
c     zone number rules
c     
         mc = 0
         
         do i = 1, nriver
c     RJP 07/21/07 added a new flag stored in 'iwsp' array, to specify 
c     whether details of casing thickness and cement thickness behind 
c     casing are specified. If the flag is zero then no casing and cement 
c     is present. The thickness are specified in subsequent row.
            read(inpt,*)   inriverf(i),isriverf(i),ifriverf(i),iwsp(i)
            npoint_riv =   npoint_riv + ifriverf(i)
            nmodels = nmodels+1
c     eliminate duplicate models           
            do ii = 1, nriver-1
               if(isriverf(i).eq.isriverf(ii)) nmodels = nmodels - 1
            enddo
            
c     inriverf(i) section number(id) of section i 
c     isriverf(i) number of layers for the  section 
c     isriver(j) is the section type of river node j 
c     ifriverf(i) number of coordinate points in the section
c     
c     for all models
c     
c     can skip nodes (with a negative number as in other macros)
c     
c     read a list of node numbers and weights

c     well type simple
c     read for all section points: x,y,z, radius of well, radius of largest layer,
               ic = npoint_riv - ifriverf(i) +1 
               do ii = 1,ifriverf(i)+1
                  read (inpt, '(a80)') wdd1
                  if (null1(wdd1)) go to 10 
                  backspace inpt
c     RJP 07/16/07 changed area_riv below to go from 2 to 4.
c     The default is that the first three numbers
c     are specified explicitly by the user in the input file. 
c     If the user does not want to specify casing and annular cement 
c     explicitly, then the second and third number should be zeros
c     and fourth number should be the outside of wellbore region
                  read (inpt, *) mb, (coor_riv(iabs(mb)+ic-1,j),j=1,3),
     &                 area_riv(iabs(mb)+ic-1,1), 
     &                 (csg_th(iabs(mb)+ic-1,j),j=1,2),
     &                 area_riv(iabs(mb)+ic-1,2)       
                  if (mb .lt. 0) then
                     call riv_interpolate(1,ic,mb,mc,i,isriver(i))
                     call area_interpolate(1,ic,mb,mc,i,isriver(i))
c     go to 10
                  endif
                  mc = iabs(mb)
               enddo 
 10            continue
c     mc = mc + iabs(mb)

               isriver(ic:npoint_riv) = isriverf(i)
               continue
         enddo
            else if(iriver.eq.2) then
c
c simple list of nodes defined - then filled in
c
c input production coordinates.
c
c n_well_prod 
c iwell_prod
c coor_well2 coordinates for simple wells
c n_well_seg pipe lengths (length connected to iwell_prod nodes)
c
c  well_rad(n_well_seg)-well_radius for each segment
c  well_dz(n_well_seg)- distance increment for gridblock along well segment
c  well_connect_fac - multiplier of well to primary grid 
c      
             n_well_prod = nriver
	       allocate(iwell_prod(n_well_prod))
	       allocate(coor_well2(n_well_prod,3))
	       allocate(well_connect_fac(n_well_prod))
	       allocate(izlabelp(n_well_prod))	       
             do i = 1,n_well_prod
	        read(inpt,*) ii,
     &			coor_well2(ii,1),coor_well2(ii,2),coor_well2(ii,3),
     &          well_connect_fac(ii),izlabelp(ii)
                iwell_prod(i) = ii   
                well_connect_fac(ii) = 
     &            max(abs(well_connect_fac(ii)),well_connect_tol)
                
             enddo
             read(inpt,*) n_well_seg
	       allocate(iwell_seg(n_well_seg,2))
	       allocate(well_rad(n_well_seg))
	       allocate(well_dz(n_well_seg))
	       allocate(izlabels(n_well_seg))		       
             do i = 1,n_well_seg
	        read(inpt,*) j,iwell_seg(j,1),iwell_seg(j,2),
     &			well_rad(j),well_dz(j),izlabels(j)	   
             enddo
            continue
            else if(iriver.eq.3) then
c high resolution 3D embedded patch of cells 
            else if(iriver.eq.4) then
c surface grid nodes     
            else
            endif
c     
c     finished reading  river section and attributes      
c     end loop on sections 
c     
        else if(macro1.eq.'wellphysi ') then 
c
c  read in well identification
c
         read(inpt,'(a80)') dum_string(1:80)
	   if(dum_string(1:8).eq.'allwells') then
	    
c apply to all well nodes
	     iwell_phy_all = 0
	     read (inpt,*) iwell_phy_all
	    else if(dum_string(1:8).eq.'wellzone') then
c identify wells
         endif
	     
        else
         write(iout,*) 'ERROR IN WELL INPUT(STOPPING)'
         write(*,*) 'ERROR IN WELL INPUT(STOPPING)'
	   stop
        end if
40     continue
      go to 100
200   continue         
      else if(iflg.eq.1) then  
c
c  need a loop on the number of types of wells
c
	  if(iriver.eq.1) then
c     
c     add nodes to fehm list (connectivity will be established later)
c     called after have coordinates (in startup)
c     
c     complete area (radius information)    
c     
c     load mdnodes array
c     
c     following code assumes all river or wells have the same number of layers
c     
         allocate (mdnodes(npoint_riv*maxlay))
         allocate (mdnode(npoint_riv*maxlay,maxc))
         allocate (coor_dum(npoint_riv,3))
c     iriver_list will contain the global number of the well or river
c     not including the intermediate nodes
         allocate (iriver_list(npoint_riv))
c     iriver_con_node will contain the global number of the primary node
c     that is closest to the river/well node
         allocate (iriver_con_node(npoint_riv*maxlay))
         allocate (dum(npoint_riv*maxlay))
         allocate (dum1(maxlay))
         allocate (iriver_out_node(maxiriver+neq))
         allocate (iriver_first_node(maxiriver+neq))
c         coor_dum = coor_riv
         coor_dum(1:npoint_riv,1:3) = coor_riv(1:npoint_riv,1:3)
         npoint_riv_old = npoint_riv
         ic = 0
c     RJP 11/30/07 added following
         i4 = 0
         do ii = 1, nriver
            rivbegin(ii) = ic+1
            do i2 = 1, ifriverf(ii)
               ic = ic+1
               i4 = i4 + 1
               isection = isriver(i4)
               do i = 1, isection
                  ic = ic+1
               enddo
            enddo
            rivend(ii) = ic
         enddo
         ic = 0
         i3 = 0
c     RJP 11/30/07 added following
         i4 = 0
         maxlay = 0
         iexpfl = 0
         do ii = 1, nriver
            do i2 = 1, ifriverf(ii)
c     RJP 11/30/07 added following
               i4 = i4 + 1
c     RJP 11/30/07 changed isriver(i2) below to isriver(i4)
               isection = isriver(i4)
               maxlay = max(maxlay,isection+1)
               ic = ic +1
               i3 = i3 +1
               iriver_list(i3) = ic
               call primary_connect(1,i3,ic,isection,node)
               mdnodes(ic) = 4
               mdnode(ic,1) = ic-isection-1
               if((ic-isection-1).lt.rivbegin(ii)) mdnode(ic,1) = 0
               mdnode(ic,2) = ic
               mdnode(ic,3) = ic+1
               mdnode(ic,4) = ic+isection+1
               if((ic+isection+1).gt.rivend(ii)) mdnode(ic,4) = 0
               dum(ic) = area_riv(i3,1)
               radmax = area_riv(i3,2)
c     call is with isection +1 because wellbore is not in the section count
c     call mult(isection+1,dum(ic),radmax,amult)
c     RJP 11/28/05 changed the first argument in above statement back to 
c     isection to get the multiplier correct
c     Changed the multiplier specifier below. The default is that the first 
c     three numbers are specified explicitly by the user in the input file. 
c     If the user does not want to specify casing and annular cement 
c     explicitly, then the second and third number should c be zeros and 
c     fourth number should be the outside of the wellbore region. The 
c     number of onion skins or detailed regions are changed based on user 
c     input on the casing and annular cement. Overall radial sections are 
c     kept the same as the second number in the first line of wellbore 
c     specification.
               if (iwsp(ii).eq.1) then
                  dum1(2) = csg_th(i3,1)/2.d0
                  dum1(3) = csg_th(i3,2)/2.d0
                  rads = dum(ic)+2.d0*(dum1(2)+dum1(3))
                  call mult(isection-1,rads,radmax,amult)

                  do i = 4, isection+1
                     dum1(i) = rads*amult**(i-3)
                  enddo

                  dum(ic) = 0.d0
                  vol0 = dum(ic)
                  dum(ic+1) = area_riv(i3,1)+dum1(2)
                  vol0 = vol0+dum(ic+1)
                  dum(ic+2) = dum(ic+1)+dum1(2)+dum1(3)
                  vol0 = dum(ic+2)+dum1(3)
                  jj = ic+2
                  do i = 1, isection-2
                     jj = jj +1
                     dum(jj) = vol0 +dum1(i+3)/2.d0
                     vol0 = vol0 +dum1(i+3)
                  enddo
                  dum(jj) = vol0	
               else
                  rads = 2.d0*dum(ic)
                  call mult(isection,rads,radmax,amult)
                  do i = 2, isection
                     dum1(i) = 2.d0*dum(ic)*amult**(i-1)
                  enddo
                  jj = ic+1
                  dum(ic+1) = 2.d0*dum(ic)
                  dum(ic) = 0.d0
                  vol0 = dum(ic+1)+dum(ic)
                  do i = 1, isection-1
                     jj = jj +1
                     dum(jj) = vol0 +dum1(i+1)
                     vol0 = dum(jj)
                  enddo	

               endif

               do j = 1,isection-1
                  ic = ic + 1
                  mdnodes(ic) = 5
                  mdnode(ic,1) = ic-1
                  mdnode(ic,2) = ic 
                  mdnode(ic,3) = ic+1
                  mdnode(ic,4) = ic-isection-1
                  mdnode(ic,5) = ic+isection+1
                  if ((ic-isection-1).lt.rivbegin(ii)) mdnode(ic,4) = 0
                  if ((ic+isection+1).gt.rivend(ii)) mdnode(ic,5) = 0
               enddo
               ic = ic +1
               mdnodes(ic) = 5
               mdnode(ic,1) = ic-1
               mdnode(ic,2) = ic
               mdnode(ic,3) = -node
               mdnode(ic,4) = ic-isection-1
               mdnode(ic,5) = ic+isection+1
               if ((ic-isection-1).lt.rivbegin(ii)) mdnode(ic,4) = 0
               if ((ic+isection+1).gt.rivend(ii)) mdnode(ic,5) = 0

               iriver_con_node(ic-isection:ic) = node
            enddo
         enddo
         mdnode(1,1) = 0
         mdnode(ic-isection,4) = 0

         npoint_riv = ic

c     RJP 11/28/05 changed below
c     do ii = 1,ic
c     area_riv(ii,1) = dum(ii)
c     enddo
         do ii = 1,ic
            area_riv(ii,1) = dum(ii)
         enddo
c     do ii = 1, npoint_riv
c     coor_riv(ii,1) = coor_riv(ii,1)+area_riv(ii,1)
c     coor_riv(ii,2) = coor_riv(ii,2)+area_riv(ii,1)
c     enddo
         deallocate(dum)
         allocate (idum(npoint_riv))
         idum = 0
         do i = 1, npoint_riv_old 
            idum(iriver_list(i)) = isriver(i)
         enddo
         deallocate(isriver)
         allocate(isriver(npoint_riv))
         isriver = idum
         deallocate (idum)        
         deallocate (coor_dum)
         allocate (coor_dum(neq_primary+npoint_riv,3))
         coor_dum(1:npoint_riv,1:3) = coor_riv(1:npoint_riv,1:3)
         deallocate (coor_riv)
         allocate (coor_riv(npoint_riv,3))
c	 coor_riv = coor_dum
         coor_riv(1:npoint_riv,1:3) = coor_dum(1:npoint_riv,1:3)
         do j= 1,3
            do i = 1,neq_primary
               coor_dum(i,j) = cord(i,j)
            enddo
         enddo  
         deallocate(cord)
c     add additioanal coordinates and identify local and global nodes
c     count number of models
         neq_new = neq_primary+npoint_riv
         allocate(river_nodes_local(neq_new))
         allocate(river_nodes_global(npoint_riv))
         allocate(izone_river_nodes(neq_new))   
         izone_river_nodes = 0
         ic = 0
         nmodels = 0
         do inode = neq_primary+1, neq_new 
            ic = ic +1
            river_nodes_local(inode) = ic
            river_nodes_global(ic)   = inode
            izone_river_nodes(inode) = isriver(ic)
            if (mod(ic,isection+1).eq.0) then
               node = iriver_con_node(ic)
               iriver_out_node(node) = river_nodes_global(ic)
            endif
            do j = 1,3
               coor_dum(inode,j) = coor_riv(ic,j)
            enddo
         enddo

         ic = 0
         do inode = neq_primary+1, neq_new
            ic = ic + 1
            if (isriver(ic).ne.0) then
               node = iriver_con_node(ic)
               iriver_first_node(node) = river_nodes_global(ic)
               iriver_out_node(node) = iriver_first_node(node)
     &              +isriver(ic)
            endif
         enddo
c     resize arrays 
c     write(*,*) "Here 2"
         
         allocate (cord(neq_new,3))
         cord = coor_dum
c     release  memory         
         deallocate (coor_dum,coor_riv)
c
       else if(iriver.eq.2) then
c      generate connectivity for implicit well system
        call implicit_well (1)
c
	 endif

      else if(iflg.eq.2) then
c     
c     redefine mdnodes to global
c     
c     
c     mdnodes is local to river or well nodes
c    
        if(iriver.eq.1) then 
         neq_new =  neq
         neqp1 = neq_primary + 1
         allocate(mdnodes_riv(neq_new))
         allocate(mdnode_riv(neq_new,maxlay))
         mdnode_riv = 0 
         mdnodes_riv = 0
         ii = 0
         do iii = 1, npoint_riv
            j = mdnodes(iii) 
            inode = river_nodes_global(iii)  
            mdnodes_riv(inode) = j
            ii = ii+j
            do jj = 1,j
               ij = mdnode(iii,jj)
               if(ij.gt.0) then
                  ik = river_nodes_global(ij)
                  mdnode_riv(inode,jj) = ik
               else if(ij.lt.0) then
c     special case where abs(ij) is the primary node                    
                  ik = -ij
                  mdnode_riv(inode,jj) = ik
                  mdnodes_riv(ik)= mdnodes_riv(ik)+1
                  mdnode_riv(ik,mdnodes_riv(ik)) = inode
               endif
            enddo
         enddo
         n_ncon = 2*ii +nelm(neqp1) + npoint_riv
       else if(iriver.eq.2) then
       endif
      else if(iflg.eq.3) then
c     
c     resize volumes
c     
c     called after connectivity 
c     but is needed for sx (coefficient)  calculations        
c     
c     fill in volumes and areas (really areas divided by distance)
c     
c     
       if(iriver.eq.1) then
c
         allocate(river01(neq_new))
         allocate(river02(neq_new)) 
         allocate(river03(neq_new)) 
         allocate(vol1(npoint_riv))
         allocate(cv_rad(npoint_riv))
         allocate (coor_dum(neq_new,1))
         river02 = 0.0d0 	 
c     kk is a count of the wellbore gridblocks  
         kk = 0  
         do ii = 1, nriver
            do jj = rivbegin(ii), rivend(ii) 
               ik = river_nodes_global(jj)
               id = isriver(jj)
               isection = id
               
               if(id.ne.0) then
c     this is a river node
                  kk = kk +1
                  imr = mdnode_riv(ik,1)
                  ipr = mdnode_riv(ik,4)
c     dis0 is length of river or well gridblock
                  if (jj.eq.rivbegin(ii)) then
c     first node is 1/2 cell length, dis to next cell center is ful
                     disx = (cord(ik,1) - cord(ipr,1))
                     disy = (cord(ik,2) - cord(ipr,2))
                     disz = (cord(ik,3) - cord(ipr,3))
                     disa = sqrt(disx**2+disy**2+disz**2)
                     dis0 = disa/2.0
c                    dis0 = disa
c     last node is 1/2 cell
                  else if (jj.eq.rivend(ii)-isection) then
                     disx = (cord(imr,1) - cord(ik,1))
                     disy = (cord(imr,2) - cord(ik,2))
                     disz = (cord(imr,3) - cord(ik,3))
                     disa = sqrt(disx**2+disy**2+disz**2)
                     dis0 = disa/2.0
c                    dis0 = disa
                  else
c     middle cell are full cell length, full connection
                     disx = (cord(imr,1) - cord(ipr,1))
                     disy = (cord(imr,2) - cord(ipr,2))
                     disz = (cord(imr,3) - cord(ipr,3))
                     dis0 = sqrt(disx**2+disy**2+disz**2)/2.0
                     disx = (cord(ik,1) - cord(ipr,1))
                     disy = (cord(ik,2) - cord(ipr,2))
                     disz = (cord(ik,3) - cord(ipr,3))
                     disa = sqrt(disx**2+disy**2+disz**2) 
                  endif
c     rad0 is outer radius in well gridblock

                  ic = ik-1
c     RJP 7/20/07 changed below
                  if(iwsp(ii).eq.1) then
                     do j = 1, isection+1
                        ic = ic+1
                        area_riv(river_nodes_local(ic),2)=dis0
                        if (j .eq. 1) then 
                           cv_rad(j) = area_riv(river_nodes_local(ic+1),

     &				1)-csg_th(kk,1)/2.d0
                        elseif (j .eq. 2) then
                           cv_rad(j) = area_riv(river_nodes_local(ic+1),

     &				1)-csg_th(kk,2)/2.d0
                        elseif (j .eq. 3) then
                           cv_rad(j) = area_riv(river_nodes_local(ic),1)
     &				+csg_th(kk,2)/2.d0
                        elseif (j.lt.isection+1) then
                           rad0 = area_riv(river_nodes_local(ic),1)
                           rad1 = area_riv(river_nodes_local(ic+1),1)
                           cv_rad(j) = (rad0+rad1)/2.0
                        endif
                     enddo
                  else
                     do j = 1, isection+1
                        ic = ic+1
                        area_riv(river_nodes_local(ic),2)=dis0
                        if (j.lt.isection+1) then
                           rad0 = area_riv(river_nodes_local(ic),1)
                           rad1 = area_riv(river_nodes_local(ic+1),1)
                           cv_rad(j) = (rad0+rad1)/2.0
                        endif
                     enddo
                  end if
                  ic = ik-1
                  do j = 1, isection+1
                     ic = ic+1
                     coor_dum(ic,1)= area_riv(river_nodes_local(ic),1)
c     if (j.eq.1) then
c     coor_dum(ic,1) = 0.0
c     else
c     coor_dum(ic,1) = cv_rad(j-1)
c     endif
                  enddo
                  ic = ik-1
                  do j = 1, isection+1
                     ic = ic+1
                     rad1 =  area_riv(river_nodes_local(ic),1)
                     if (j.ne.isection+1) then
                        rad2 =  area_riv(river_nodes_local(ic+1),1)
                     endif
                     if (j.eq.1) then
                        river03(ic) = pi*cv_rad(j)**2/disa
                        vol0 = river03(ic)*disa*dis0
                     elseif (j.eq.isection+1) then
                        river03(ic) = pi*(rad1**2-cv_rad(j-1)**2)/disa
                        vol0 = river03(ic)*disa*dis0
                     else
                        river03(ic) = pi*(cv_rad(j)**2-cv_rad(j-1)**2)/
     &                       disa
                        vol0 = river03(ic)*disa*dis0
                     endif
                     river01(ic) = vol0
                     vol1(kk) = vol1(kk)+vol0
                     disr = rad2 - rad1
                     if (j.eq.isection+1) disr = rad2-
     &                    area_riv(river_nodes_local(ic-1),1)
                     if (j.eq.1) then
                        area0 = 0.5*dis0*pi2*((cv_rad(j)/2.)+
     &                       (cv_rad(j)+cv_rad(j+1))/2.)
                     elseif (j.eq.isection+1) then
                        area0 = 0.5*dis0*pi2*(((cv_rad(j)+cv_rad(j-1))/
     &                       2.)+(cv_rad(j)/2.))
                     elseif (j.eq.isection) then
                        radmax =  area_riv(river_nodes_local(ic+1),1)
                        area0 = 0.5*dis0*pi2*(((cv_rad(j)+cv_rad(j-1))/
     &                       2.)+(cv_rad(j)+radmax)/2.)
                     else
                        area0 = 0.5*dis0*pi2*(((cv_rad(j)+cv_rad(j-1))/
     &                       2.)+(cv_rad(j)+cv_rad(j+1))/2.)
                     endif
                     river02(ic) = area0/disr
c     arbitrarily set x coordiate of well to be maximum radial distance away
                  enddo

               endif
            enddo
         enddo   
         do i = neq_primary+1,neq_new
            cord(i,1) = cord(i,1)+coor_dum(i,1)
c     cord(i,2) = 0.0
         enddo 
         continue
         deallocate(coor_dum)
c     
c     modify volumes
c     
         allocate(dum(neq_new))  
         dum(1:neq_primary) = sx1(1:neq_primary) 
         deallocate (sx1)
c         open(unit=1002,file="well_connectivity.dat")
         well_file_name = ''
         well_file_name(1:iroot) = well_root(1:iroot)
         well_file_name(iroot+1:iroot+13)='.connectivity'
         inwel = open_file(well_file_name, 'unknown')

         do i = neq_primary + 1, neq_new
            ii = river_nodes_local(i)
            node = iriver_con_node(ii)  
c     node is the primary gridblock asociated with river(well node)
c     subtract river volume from primary node volume
            dum(node) = dum(node) - river01(i)
            dum(i) = river01(i)
            write(inwel,1003) i, ii, node, dum(node), river01(i)
         enddo
 1003    format(i6,i6,i6,1x,f14.4,1x,f14.4)   
         allocate (sx1(neq_new))
         sx1 = dum
         close(inwel)
         deallocate(dum,vol1)
       else if(iriver.eq.2) then
c
c  
c	   
	 endif
      else if(iflg.eq.-3) then
c     
c     printout some well or river properties
c     RJP 11/21/05
c     
       if(iriver.eq.1) then
         well_file_name = ''
         well_file_name(1:iroot) = well_root(1:iroot)
         well_file_name(iroot+1:iroot+4)='.wel'
         inwel = open_file(well_file_name, 'unknown')

         write(inwel,*) "WellNode Zone GlobaNode GlobalX+radius GlobalZ"
         do i = 1, npoint_riv
            corz(river_nodes_global(i),1)=cord(river_nodes_global(i),1)
            corz(river_nodes_global(i),2)=cord(river_nodes_global(i),2)
            corz(river_nodes_global(i),3)=cord(river_nodes_global(i),3)
            write(inwel,666) i, river_zone(i),river_nodes_global(i), 
     &           cord(river_nodes_global(i),1),
     &           cord(river_nodes_global(i),3)
         enddo
         close (inwel)
 666     FORMAT(I5,2x,I5,2x,I8,2x,g12.6,2x,g12.6)
       else if(iriver.eq.2) then
        call implicit_well(-3)
       endif
      else if(iflg.eq.4) then
C     
C     modify connectivity
C     
C     estimate array sizes
C     
C     are the additional connections
c     
c     add total number of connections (ii +ii) due to new connections in primary and river nodes
c     
        if(iriver.eq.1) then
c
         neqp1_old = neq_primary + 1
         neqp1 =  neq_new +1  
         allocate(istrw_new(n_ncon-neqp1))
         allocate(ncon_new(n_ncon))
         allocate(idum(2*neq_new))
         allocate(idum1(neq_primary))
         allocate(dum(n_ncon-neqp1))
         j=neqp1
         icoef = nr
         nr_old = nr-1      
         ncon_new(1)=neqp1
         idum1(1:neq_primary)= nelmdg(1:neq_primary)
         deallocate (nelmdg)
         allocate (nelmdg(neq_new))
         do i=1,neq_primary
            idum = 0
            idum(i)=1
            kbmin=neq
            kbmax=0
            kbmin=min(kbmin,i) 
            kbmax=max(kbmax,i)
            i1=nelm(i)+1     
            i2=nelm(i+1) 
            do jj=i1,i2
               kb=nelm(jj)
               kbmin=min(kbmin,kb)
               kbmax=max(kbmax,kb)
               idum(kb)=1
               idum(neq_new+kb)=istrw(jj-neqp1_old)
               if(idum(neq_new+kb).eq.nr) 
     &              idum(neq_new+kb) = -99999999
            enddo
            sx_max= 0.0d0
            i1 = idum1(i)+1
            do jj=i1,i2
               ii = istrw(jj-neqp1_old)
               vol0 = abs(sx(ii,isox)+sx(ii,isoy)+sx(ii,isoz))
               sx_max = max(sx_max,vol0)
            enddo
C     RJP Modify later for non-centered grid
            sx_max = mult_well*sx_max
C     identify additional connections                  
            if(mdnodes_riv(i).ne.0) then
               do k=1,mdnodes_riv(i)
                  kb=mdnode_riv(i,k)
                  kbmin=min(kbmin,kb)
                  kbmax=max(kbmax,kb)
                  idum(kb)=1
                  if(kb.gt.i) then
                     icoef = icoef +1
                     idum(neq_new+kb) = icoef
                     dum(icoef) = -sx_max
                  else if(kb.lt.i) then
                     do jk = nelmdg(kb)+1, ncon_new(kb+1)
                        if(ncon_new(jk).eq.i) then
                           idum(neq_new+kb) = istrw_new(jk-neqp1)
                        endif
                     enddo
                  endif
               enddo
            endif
            do k=kbmin,kbmax
               if(idum(k).ne.0) then
                  j=j+1
                  ncon_new(j)=k
                  istrw_new(j-neqp1)=idum(k+neq_new)
                  if(k.eq.i) nelmdg(i)=j
                  idum(k)=0
                  idum(k+neq_new)=0
               endif
            enddo
            ncon_new(i+1)=j
         enddo
C     
C     now add connections for river nodes
C     
         ic = j 
         do i = neq_primary+1, neq_new
            idum = 0
            idum(i)=1
            kbmin=neq
            kbmax=0
            kbmin=min(kbmin,i)
            kbmax=max(kbmax,i)
            j = mdnodes_riv(i)
            do jj = 1,j
               kb = mdnode_riv(i,jj)
               if (kb.ne.0) then
                  kbmin=min(kbmin,kb)
                  kbmax=max(kbmax,kb)
                  idum(kb) = 1
                  if(kb.gt.i+1) then
                     icoef = icoef +1
c     note reference to node i for this coefficient
                     if(abs(river03(i)).gt.area_tol) then
                        idum(neq_new+kb) = icoef
                     else
                        idum(neq_new+kb) = -icoef
                     endif
                     dum(iabs(icoef)) = -river03(i)                 
                  else if(kb.eq.i+1) then
c     note reference to node i for this coefficient
                     icoef = icoef +1
                     if(abs(river02(i)).gt.area_tol) then
                        idum(neq_new+kb) = icoef
                     else
                        idum(neq_new+kb) = -icoef
                     endif
                     dum(iabs(icoef)) = -river02(i)
                  else if(kb.lt.i) then
                     do jk = nelmdg(kb)+1, ncon_new(kb+1)
                        if(ncon_new(jk).eq.i) then
                           idum(neq_new+kb) = istrw_new(jk-neqp1)
                        endif
                     enddo
                  endif
               endif
            enddo 
            do k = kbmin,kbmax
               if(idum(k).ne.0) then
                  ic = ic +1
                  ncon_new(ic) = k
                  istrw_new(ic-neqp1)=idum(k+neq_new)
                  if(k.eq.i) nelmdg(i) = ic
               endif
            enddo
            ncon_new(i+1) = ic
         enddo 
c     
c     deallocate and allocate some arrays
c     note:nelmdg already loaded
c     
         j = ic
         deallocate(nelm,istrw,istrw_itfc)
         allocate (nelm(j))
         allocate (istrw(j-neqp1))
         allocate (istrw_itfc(j-neqp1))
         do i=1,j
            nelm(i)=ncon_new(i)
         enddo
c     RJP 01/08/07 added following
         ncont=ncon_new(neq+1)
         nelmd = ncont
         do i=1,j-neqp1 
            if(istrw_new(i).eq.-99999999) istrw_new(i) = icoef+1
            istrw(i)=istrw_new(i)
         enddo
         nr = icoef+1
         
         if(isoy.ne.1) then
            allocate (sx_new(icoef+1,3))
         else
            allocate (sx_new(icoef+1,1))
         endif
         
         do i=1,nr_old
            sx_new(i,isox) = sx(i,isox)
            sx_new(i,isoy) = sx(i,isoy)
            sx_new(i,isoz) = sx(i,isoz)
         enddo
         
         deallocate(sx)
         if(isoy.ne.1) then
            allocate (sx(icoef+1,3))
         else
            allocate (sx(icoef+1,1))
         endif
         
         do i=1,nr_old
            sx(i,isox) = sx_new(i,isox)
            sx(i,isoy) = sx_new(i,isoy)
            sx(i,isoz) = sx_new(i,isoz)
         enddo
         
         do i = nr_old+1,icoef
            if(isoy.ne.1) then
               sx(i,1) = dum(i)             
               sx(i,2) = 0.0d00          
               sx(i,3) = 0.0d00          
            else
               sx(i,1) = dum(i)/3.             
            endif
         enddo

         istrw_itfc  = 0
         deallocate(sx_new,istrw_new)

         neq = neq_new
         neqpi = neq
         allocate(istrw_new(j-neqp1))
         istrw_new = istrw
c     
c     need to modify interface pointer
c     (for now)
c     
         if (isoy.ne.1) then
            allocate(sx_new(2*icoef,3))
         else
            allocate(sx_new(2*icoef,1))
         endif
         do i = 1, icoef
            sx_new(i,isox)=sx(i,isox)
            sx_new(i,isoy)=sx(i,isoy)
            sx_new(i,isoz)=sx(i,isoz)
         enddo
         ic = icoef
         nic_old = icoef
c RJP 05/21/08 added following to save modified distances for calculations in Geneqs
         ic1 = 0
         allocate(mod_dis(npoint_riv_old*10,2))
         mod_dis = 0.d0
         do i = 1, neq_primary
            if (mdnodes_riv(i).ne.0) then
               irnode = iriver_out_node(i)
               outrad = area_riv(river_nodes_local(irnode),1)
               ht = area_riv(river_nodes_local(irnode),2)
               outarea = pi2*outrad*ht/4.
               isnode=  iriver_first_node(i)
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do j = i1, i2
                  kb = nelm(j)
                  iw = istrw(j-neqpi-1)
                  if ((iw.ne.0).and.(kb.ne.irnode))then
                     sx_dum = sx(iw,isox)
c     calculate distance between the two nodes
                     disx = (cord(isnode,1) - cord(kb,1))
                     disy = (cord(isnode,2) - cord(kb,2))
                     disz = (cord(isnode,3) - cord(kb,3))
                     disa = sqrt(disx**2+disy**2+disz**2)
c     calculate area from distance and sx here it is assumed that 
c     sx(isox), sx(isoy), sx(isoz) are all equal.
c     calculate the new area as average between
c     two areas
                     if(disz.eq.0) then 
                        area_dum = -3.*sx_dum*disa
                        area_new = (area_dum+outarea)/2.
                        area_new = outarea*0.8
                        sx_dum = -area_new/(disa-outrad)
                        ic = ic+1
                        istrw_new(j-neqpi-1)=ic
                        sx_new(ic,isox)=sx_dum
c RJP 05/21/08 calculating and storing new distances for calculations in GENEQs
c Added a new array to store the distances
                        ic1 = ic1+1
                        if(disx.eq.0.d0) then
                           mod_dis(ic1,1) = disx
                        else
                           mod_dis(ic1,1)=dabs(disx)-outrad
                        endif
                        if(disy.eq.0.d0) then
                           mod_dis(ic1,2) = disy
                        else
                           mod_dis(ic1,2)=dabs(disy)-outrad
                        endif
 
                        if (kb.lt.i) then
                           j1 = nelm(kb)+1
                           j2 = nelm(kb+1)
                           do k = j1, j2
                              kb2 = nelm(k)
                              if (kb2.eq.i) then
                                 istrw_new(k-neqpi-1)=ic
                              endif
                           enddo
                        endif 
                     endif
                  end if
               enddo
            endif
         enddo
         deallocate(sx)
         if (isoy.ne.1) then
            allocate(sx(ic,3))
         else
            allocate(sx(ic,1))
         endif

         do i = 1, icoef
            sx(i,isox) = sx_new(i,isox)
            sx(i,isoy) = sx_new(i,isoy)
            sx(i,isoz) = sx_new(i,isoz)
         enddo

         do i = icoef+1, ic
            if (isoy.ne.1) then
               sx(i,1) = sx_new(i,isox)
               sx(i,2) = 0.0d00
               sx(i,3) = 0.0d00
            else
               sx(i,1) = sx_new(i,1)/3.
            endif
         enddo
         nr = ic + 1
         istrw = istrw_new
c         neq_primary = neq
         deallocate(idum,dum,istrw_new,ncon_new,sx_new)
         deallocate (river01,river02,river03)
       elseif(iriver.eq.2) then
c simple peaceman type wellbore
        call implicit_well(4)
       endif
      else if(iflg.eq.-4) then
c     
c     printout connectivities
c     

      else if(iflg.eq.5) then
         
c     
c     assign zones to well or river plus layers
c       
       if(iriver.eq.1) then
c
         ic = neq-npoint_riv
         do i = 1,nriver
            i1 = inriverf(i)
            jj = isriverf(i)
            i2 = ifriverf(i)   	      
            
            do kk = 1,i2
               ic = ic + 1
               if(kk.eq.1) then

                  izonef(ic) = 1000*(i1) 
                  do j = 1,jj
                     ic = ic + 1
                     izonef(ic) = izonef(iriver_con_node(ic - neq
     &                    +npoint_riv))
                  enddo
               else if(kk.eq.i2) then
                  izonef(ic) = 1000*(i1) + 900
                  do j = 1,jj
                     ic = ic + 1
                     izonef(ic) = izonef(iriver_con_node(ic - neq
     &                    +npoint_riv))
                  enddo		  	
               else 
                  izonef(ic) = 1000*(i1) + kk - 1    
                  do j = 1,jj
                     ic = ic + 1
                     izonef(ic) = izonef(iriver_con_node(ic - neq
     &                    +npoint_riv))
                  enddo
               endif

            enddo

         enddo
         continue

         do i = 1, npoint_riv
            river_zone(i) = izonef(neq-npoint_riv+i)
         enddo
	 else if(iriver.eq.2) then
c
	  call implicit_well(5)
c
	 endif

      else if(iflg.eq.-5) then  
c     
c     printout zonation
c     
      else if(iflg.eq.33) then  
c     
c     printout zonation
c     
c     n = neq-npoint_riv
c
	 if(iriver.eq.1) then
         do i = neq-npoint_riv+1, neq
            if(pnx(i).eq.zero_t) then
               pnx(i) = pnx(iriver_con_node(i-neq+npoint_riv))
               ps(i) = ps(iriver_con_node(i-neq+npoint_riv))
               denr(i) = denr(iriver_con_node(i-neq+npoint_riv))
               cpr(i) = cpr(iriver_con_node(i-neq+npoint_riv))
            endif
            if(pny(i).eq.zero_t) then
               pny(i) = pny(iriver_con_node(i-neq+npoint_riv))
            endif
            if(pnz(i).eq.zero_t) then
               pnz(i) = pnz(iriver_con_node(i-neq+npoint_riv))
            endif
         enddo

	 else if(iriver.eq.2) then
c
	  call implicit_well(33)
c
	 endif

      elseif (iflg.eq.-33) then
c
	 if(iriver.eq.1) then
c
         n = neq-npoint_riv
         do i = n+1, n+npoint_riv
            pnx(i) = -9999
            pny(i) = -9999
            pnz(i) = -9999
            ps(i) = -9999
            denr(i) = -9999
            cpr(i) = -9999
         enddo
	 else if(iriver.eq.2) then
c
	  call implicit_well(-33)
c
	 endif
      else if(iflg.eq.6) then
c     
c     output well or river information
c     timesteps
c     
c     implement later
       if(iriver.eq.1) then
c     
	 else if(iriver.eq.2) then
c
	  call implicit_well(6)
c
	 endif        

C     *******************************************************
      endif

      return
      end

      subroutine riv_interpolate(iz,ic,mb,mc,id,ie)
c     
c     interpolate quantities for river macro
c     
      
      use combi      
      use comci      
      use comdti      
      use comai 
      use comriv    
      implicit none            
      integer iz,ie,ic
      integer mb,mc,id,io,i,ir,mdif,mdifm,mn
      real*8 dx,dy,dz

      mdif   =  iabs(mb)-mc
      mdifm  =  mdif-1
      mn     =  mc+ic-1
c     
c**** interpolate on coordinates ****
c     
      dx  =  ( coor_riv(iabs(mb)+ic-1,1)-coor_riv(mn,1) )/mdif
      dy  =  ( coor_riv(iabs(mb)+ic-1,2)-coor_riv(mn,2) )/mdif
      dz  =  ( coor_riv(iabs(mb)+ic-1,3)-coor_riv(mn,3) )/mdif
      do io = 1, mdifm
         coor_riv(mn+io,1) =  coor_riv(mn,1)+io*dx
         coor_riv(mn+io,3) =  coor_riv(mn,3)+io*dz
         coor_riv(mn+io,2) =  coor_riv(mn,2)+io*dy
      end do
      return
      end

      subroutine area_interpolate(iz,ic,mb,mc,id,ie)
c     
c     interpolate quantities for river macro
c     
      
      use combi      
      use comci      
      use comdti      
      use comai 
      use comriv    
      implicit none            
      integer iz,ie,ic
      integer mb,mc,id,io,i,ir,mdif,mdifm,mn
      real*8 area1,area2,area3,area4

      mdif   =  iabs(mb)-mc
      mdifm  =  mdif-1
      mn     =  mc+ic-1
c     
c**** interpolate on coordinates ****
c     
      area1 = ( area_riv(iabs(mb)+ic-1,1)-area_riv(mn,1) )/mdif
      area2 = ( area_riv(iabs(mb)+ic-1,2)-area_riv(mn,2) )/mdif
      area3 = ( csg_th(iabs(mb)+ic-1,1)-csg_th(mn,1) )/mdif
      area4 = ( csg_th(iabs(mb)+ic-1,2)-csg_th(mn,2) )/mdif

      do io = 1, mdifm
         area_riv(mn+io,1) =  area_riv(mn,1)+io*area1
         area_riv(mn+io,2) =  area_riv(mn,2)+io*area2
c         area_riv(ic+io,1) =  area_riv(mn,1)+io*area1
c         area_riv(ic+io,2) =  area_riv(mn,2)+io*area2
         csg_th(mn+io,1) =  csg_th(mn,1)+io*area3
         csg_th(mn+io,2) =  csg_th(mn,2)+io*area4
      end do
      
      return
      end

      subroutine perm_interpolate(iz,mb,mc,id,ie)
c     
c     interpolate quantities for river macro
c     
      
      use combi      
      use comci      
      use comdti      
      use comai 
      use comriv    
      implicit none            
      integer iz,ie
      integer mb,mc,id,io,i,ir,mdif,mdifm
      real*8 perm1,perm2

      mdif   =  iabs(mb)-mc
      mdifm  =  mdif-1

c     
c**** interpolate on coordinates ****
c     
      perm1     =  ( perm_riv(iabs(mb),1)-perm_riv(mc,1) )/mdif
      perm2     =  ( perm_riv(iabs(mb),1)-perm_riv(mc,2) )/mdif

      do io = 1, mdifm
         perm_riv(mc+io,1) =  perm_riv(mc,1)+io*perm1
      end do
      
      do io = 1, mdifm
         perm_riv(mc+io,2) =  perm_riv(mc,2)+io*perm2
      end do
      
      return
      end

      subroutine primary_connect(iz,ii,ie,isect,node)
c     
c     connect river nodes to primary nodes
c     
      
      use combi      
      use comci      
      use comdti      
      use comai 
      use comriv
      
      implicit none            
      integer iz,i,ii,iii,isect,node,ie
      real*8 xriv,yriv,zriv,xdis,ydis,zdis,frac
c     
c     find closest primary node
c     
      if(iz.eq.1) then

         xriv = coor_dum(ii,1)
         yriv = coor_dum(ii,2)
         zriv = coor_dum(ii,3)
         call near3(xriv, yriv, zriv, node, 0)  
c     
c     divide distance based on section number
c     
         xdis = cord(node,1)-xriv
         ydis = cord(node,2)-yriv
         zdis = cord(node,3)-zriv
	 if(inriverf(1).ne.0) then
c     radial well model in xy direction
            xdis = 0.0d0
            ydis = 0.0d0
            zdis = 0.0d0
            if(isect.gt.1) then
            endif
c     
c     add layers if necessary
c     
            coor_riv(ie,1) = xriv 
            coor_riv(ie,2) = yriv 
            coor_riv(ie,3) = zriv     
            do i = 1, isect 
               iii = ie +i  
               frac = float(i)/float(isect)       
               coor_riv(iii,1) = xriv + xdis*frac
               coor_riv(iii,2) = yriv + ydis*frac
               coor_riv(iii,3) = zriv + zdis*frac
            enddo
         endif
      else if(iz.eq.2) then
         
      endif
      
      return
      end

      subroutine mult(n,delx,xmax,amult)
c     
c     subroutine to calculate multiplication factor for geometric spacing
c     
c     n is the number of sections (input)
c     delx is initial radius (input)
c     xmax is max radius (input)
c     amult is the grid scacing multipler(output) 

      use comai, only : ierr, iout, iptty
      implicit none
      real*8 delx,xmax,amult,tol,r,drda
      integer maxit,i,n
      parameter(tol=1.d-6,maxit=100)
      amult=2.
      do 10 i=1,maxit
         r=delx*(1.0-amult**n)/(1.0-amult)-xmax
         drda=delx*(-n*amult**(n-1)/(1.0-amult)+
     &        (1.0-amult**n)/(1.0-amult)**2)
         amult=amult-r/drda
         if(abs(r).le.tol) return
 10   continue
      write(ierr, *) 'no convergence in mult,stopping'
      if (iout .ne. 0) write(iout,*) 'no convergence in mult,stopping'
      if (iptty .ne. 0) write(iptty,*) 'no convergence in mult,stopping'
      stop
      end
