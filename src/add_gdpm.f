      subroutine add_gdpm
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
!D1 To add nodes to the connectivity array and determine the
!D1 area/distance terms and store them in the proper location.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: 03-NOV-99, Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/add_gdpm.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 08:54:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2
!D2    Rev 2.3   14 Nov 2001 13:04:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2
!D2    Rev 2.2   06 Jun 2001 13:28:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2
!D2    Rev 2.1   30 Nov 2000 11:55:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.8 Generalized dual-porosity formulation
!D3
!**********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comdti
      use comai
      use combi
      use comdi
      implicit none

      integer, allocatable :: istrw_new(:)
      integer ncont_primary
      integer ncont
      integer ngdpm_primary, i1, i2, jj
      integer i, j, n_n_n, new_place, jcounter, j_old
      integer nposition, connected_node, primary_node
      integer gdpm_counter,iw_c,kb
      integer imodel
      integer nmodels
      integer inew_gdpm
      integer inewgdpm
      integer nr_primary

      real*8, allocatable :: vol_gdpm(:)
      real*8,  allocatable :: sx_dum(:,:)
      real*8 area_factor
      real*8 delx_gdpm, sx_1,sx_2,sx_3
      real*8 vfac, vfackb, vfacmin, vfrac_pri
      real*8 sx1save, sx1gdpm, wgtgdpm
      real*8 length_total, coord_primary
      real*8 length_cell, length_pri, length_sec, gdpm_len
      real*8 gdpm_left, gdpm_right
      real*8 dis_p, dis_s

c     Set actual size of neq
c     gdpm_flag: if nonzero, determines the geometry of the matrix. 
c     1: parallel plate fractures,
c     2: spherical geometry
c     3: 1-d column geometry (edge)
c     4: 1-d column geometry (edge,based on total area)
c     5: 1-d column geometry (block-centered)
c     6: 1-d column geometry (block-centered,based on total area)

c This has already been done in setparams
c      neq = neq_primary + ngdpm_actual

c	Determine the sizes of arrays

      ncont_primary = nelm(neq_primary+1)

c	Find the number of primary nodes that have gdpm nodes

      ngdpm_primary = 0
      nmodels = 0
      do i = 1, neq_primary
         imodel = igdpm(i)
         nmodels = max(nmodels,imodel)
         if(imodel.ne.0) ngdpm_primary = ngdpm_primary + 1
      end do

c	Determine the size of the new nelm array
      ncont = ncont_primary+ngdpm_actual+ngdpm_primary+
     2     3*(ngdpm_actual-ngdpm_primary)+2*ngdpm_primary
c     New value of nr is required because it's used in sx_combine
c     to break connections of two nodes that are both flow 
c     boundary nodes. 
c     Recall that nr is the number of connections +1
c     Therefore,the new size is ncont-(neq+1)+1 = ncont-neq
      nr_primary = nr
c     nr = ncont-neq
c GAZ 4-18-2001
c the number of new connections should be equal to the number
c of new nodes since they are 1-D and symmetric
      nr = nr_primary + ngdpm_actual + 1

c	Store old nelm array, allocate space for new array
      allocate(nelm_primary(ncont_primary))
      nelm_primary = nelm
      deallocate(nelm)
      allocate(nelm(ncont))
      nelm = 0
      nelmd = ncont
c	Store old istrw arrays, allocate space for new arrays
c RJP 12/13/06 modified below for wellbore
      allocate(istrw_primary(ncont_primary-(neq_primary+1)))
      istrw_primary(1:ncont_primary-(neq_primary+1)) = 
     &              istrw(1:ncont_primary-(neq_primary+1))
      deallocate(istrw)
      allocate(istrw(ncont))
      istrw = 0
c      deallocate(istrw_itfc,istrw_cold)

c     No possibility of dpdp run if gdpm is invoked, so
c     the sizes of these arrays are ncont, not 2*ncont

c      allocate(istrw_itfc(ncont),istrw_cold(ncont))
c      istrw_itfc = 0
c      istrw_cold = 0

c	Store old sx values, allocate space for new array
      if(isoy.eq.1) then
c        allocate(sx_primary(ncont_primary-neq_primary,1))
c gaz 4-18-2001
c
         allocate(sx_primary(nr_primary,1))                    
         sx_primary = 0.
         sx_primary(1:nr_primary,1) = sx(1:nr_primary,1)
      else
c        allocate(sx_primary(ncont_primary-neq_primary,3))
c gaz 4-18-2001
c
         allocate(sx_primary(nr_primary,3))                   
         sx_primary = 0.
         sx_primary(1:nr_primary,1) = sx(1:nr_primary,1)
         sx_primary(1:nr_primary,2) = sx(1:nr_primary,2)
         sx_primary(1:nr_primary,3) = sx(1:nr_primary,3)
      end if
      deallocate(sx)
c gaz 4-18-2001
c     if(isoy.eq.1) then
c        allocate(sx(ncont-neq,1))
c     else
c        allocate(sx(ncont-neq,3))
c     end if
      if(isoy.eq.1) then
         allocate(sx(nr,1))
      else
         allocate(sx(nr,3))
      end if
      sx = 0.
c	reallocate nelmdg
      allocate(nelmdg_primary(neq_primary))
c nelmdg was allocated n0 which was already modified for gdpm nodes in
c setparams where neq was also set equal to n0
      nelmdg_primary=nelmdg(1:neq_primary)
c      deallocate(nelmdg)
c      allocate(nelmdg(neq))
      nelmdg=0

c	Reconstruct the pointer part of nelm, with new gdpm nodes
c	added

      inewgdpm = 0
      maxgdpmlayers = 0
      nelm(1) = nelm_primary(1) + ngdpm_actual
      do i = 1, neq_primary
         imodel = igdpm(i)
         if(ngdpm_layers(imodel).ne.0) then
c	This node has gdpm nodes attached
            maxgdpmlayers = 
     &      max(maxgdpmlayers,ngdpm_layers(imodel))
            inewgdpm = 1
         else
            inewgdpm = 0
         end if

c	Determine pointer position for this primary node
         nelm(i+1) = nelm(i) + 
     2   (nelm_primary(i+1)-nelm_primary(i))+inewgdpm
            
      end do
c	Set pointer array for the portion containing the new
c	gdpm nodes. Also set coordinate locations

c	n_n_n is the node number of the current new node
      n_n_n = neq_primary
      do i = 1, neq_primary
         coord_primary = cord(i,1)
         imodel = igdpm(i)
         if(ngdpm_layers(imodel).ne.0) then
            do j = 1, ngdpm_layers(imodel)-1
               n_n_n = n_n_n + 1
               cord(n_n_n,1) = coord_primary + gdpm_x(imodel,j)
               cord(n_n_n,2) = cord(i,2)
               cord(n_n_n,3) = cord(i,3)
c	3 nodes in this portion of nelm array
               nelm(n_n_n+1) = nelm(n_n_n)+3
            end do
c	2 nodes in this portion because it is the last gdpm node
               n_n_n = n_n_n + 1
               cord(n_n_n,1) = coord_primary + 
     2              gdpm_x(imodel,ngdpm_layers(imodel))
               cord(n_n_n,2) = cord(i,2)
               cord(n_n_n,3) = cord(i,3)
               nelm(n_n_n+1) = nelm(n_n_n)+2

         end if
      end do
c	This section fills the new nelm array with node numbers
c	start with the primary nodes

      do i = 1, neq_primary
         new_place = nelm(i)
         do j = nelm_primary(i)+1,nelm_primary(i+1)
            new_place = new_place + 1
            nelm(new_place) = nelm_primary(j)
         end do
      end do

c	Now fill in portion of the new array containing new nodes

      n_n_n = neq_primary

      do i = 1, neq_primary
         imodel = igdpm(i)
         if(ngdpm_layers(imodel).ne.0) then
            do j = 1, ngdpm_layers(imodel)
               n_n_n = n_n_n + 1
               if(j.eq.1) then
c	Fill in first gdpm node in primary node location
                  nelm(nelm(i+1)) = n_n_n
c	Now fill in connections for gdpm
c       first one is the primary node
                  nelm(nelm(n_n_n)+1) = i
                  nelm(nelm(n_n_n)+2) = n_n_n
                  if(ngdpm_layers(imodel).gt.1)
     2            nelm(nelm(n_n_n)+3) = n_n_n + 1
               else
                  nelm(nelm(n_n_n)+1) = n_n_n - 1
                  nelm(nelm(n_n_n)+2) = n_n_n
                  if(j.lt.ngdpm_layers(imodel))
     2            nelm(nelm(n_n_n)+3) = n_n_n + 1
                  
               end if
            end do
         end if
      end do

c	Determine nelmdg

       do i = 1, neq
          nelmloop: do j = nelm(i)+1,nelm(i+1)
             if(nelm(j).eq.i) then
                nelmdg(i)=j
                exit nelmloop
             end if
          end do nelmloop
       end do

c	Section determines the volumes (sx1)

      n_n_n = neq_primary
      
c Check for model type 4, the 1-d column based on total area
c The volumes associated with an individual model are summed
c The area associated with an individual column is the volume 
c fraction times the total area
      if(gdpm_flag.eq.4.or.gdpm_flag.eq.6) then
     
       allocate(vol_gdpm(nmodels))
       vol_gdpm = 0.0
       do i = 1,neq_primary
          imodel = igdpm(i)
          if(imodel.ne.0) then
           vol_gdpm(imodel) = vol_gdpm(imodel) + sx1(i)
          endif
       enddo
      else if(gdpm_flag.eq.11) then
       allocate(areat_gdpm(neq_primary))
      endif

      do i = 1, neq_primary
         imodel = igdpm(i)
         if(ngdpm_layers(imodel).gt.0) then
c	The node has gdpm nodes, divide up the volume
            sx1save = sx1(i)
            length_total = gdpm_x(imodel,ngdpm_layers(imodel))

c	Volume of primary node
c   Don't change primary gridblock volume if model 4 or 6 (linear column model)
            if(gdpm_flag.eq.4.or.gdpm_flag.eq.6) then
              area_factor = sx1(i)/vol_gdpm(imodel)
              sx1(i) = sx1save
            else if(gdpm_flag.le.2) then
              sx1(i) = sx1save*vfrac_primary(imodel)
            endif	
c	gdpm volume
            if(gdpm_flag.ne.11) then
             vfrac_pri = vfrac_primary(imodel)
             sx1gdpm = sx1save*(1.-vfrac_pri)
            else if(gdpm_flag.eq.11) then
c physically-based fracture model 
c coordinates in x direction may be wrong
c will have 22(y) and 33(z) models as well
c             areat_gdpm(i) = sx1save/length_total
             areat_gdpm(i) = sx1save/dxrg(i)
             vfrac_pri = vfrac_primary(imodel)
             length_sec = dxrg(i)-vfrac_pri
             gdpm_len = gdpm_vol(imodel,1)*length_sec/wgt_length(imodel)
             gdpm_left = gdpm_len
             gdpm_x(imodel,1) = gdpm_left/2.0
             do j = 2, ngdpm_layers(imodel)
              gdpm_len = gdpm_vol(imodel,j)
     &         *length_sec/wgt_length(imodel)
              gdpm_right = gdpm_len
              gdpm_x(imodel,j) = gdpm_right/2. + gdpm_left
              gdpm_left = gdpm_right + gdpm_left    
             enddo                    
             sx1(i) = sx1save*vfrac_pri/dxrg(i)
             sx1gdpm = sx1save*(dxrg(i)-vfrac_pri)/dxrg(i)
            endif
c	Loop to determine the gdpm node volumes

            do j = 1, ngdpm_layers(imodel)
               n_n_n = n_n_n + 1

               call gdpm_volume

            end do

         end if

      end do

c	Section adjusts the sx array and adds new values for the
c	gdpm nodes. Does the istrw array too

c	First fill in the old sx values and istrw

      do i = 1, neq_primary

         jcounter = 0
c	This loop is over the new nelm, j_old is the position
c	in the old nelm, and the first part of sx is the same
c	in the new as in the old, since it's the area/distance
c	of the primary grid, and all of the node numbers are
c	the same
c	Only go up to nelm(i+1)-1, because nelm(i+1) may be a
c	gdpm node. We handle that case after this loop
c     we store only uppper diagonal area/dis
         do j = nelm(i)+1, nelm(i+1)-1
            jcounter = jcounter + 1
            j_old = nelm_primary(i)+jcounter
            istrw(j-neq-1) = istrw_primary(j_old-neq_primary-1)
            iw = abs(istrw(j-neq-1))
c gaz 112409   if(nelm(j).ne.i) then            
            if(nelm(j).gt.i) then
               sx(iw,isox) = 
     2            sx_primary(iw,isox)
               sx(iw,isoy) =
     2            sx_primary(iw,isoy)
               sx(iw,isoz) =
     2            sx_primary(iw,isoz)
            end if
         end do
c	If the last node in this section is a primary node (i.e.
c	there are no gdpm nodes), then we do the same thing as
c	in the loop above. Otherwise we skip over it for now.

         j = nelm(i+1)
         if(nelm(j).le.neq_primary) then
            jcounter = jcounter + 1
            j_old = nelm_primary(i)+jcounter
            istrw(j-neq-1) = istrw_primary(j_old-neq_primary-1)
            iw = abs(istrw(j-neq-1))
            if(nelm(j).ne.i) then
               sx(iw,isox) = sx_primary(iw,isox)
               sx(iw,isoy) = sx_primary(iw,isoy)
               sx(iw,isoz) = sx_primary(iw,isoz)
            end if
         end if

      end do

c	Now fill in sx and istrw for the new connections
c	These include gdpm nodes connected to each other
c	and to a primary node

c	Every gdpm node has one node of number lower than
c	itself. We identify this connection, and then fill
c	in the position in two places in istrw, then compute sx

c    nposition = ncont_primary-neq_primary-1
c gaz 4-18-2001 new positions in sx start from nr_primary
      nposition = nr_primary

      do i = neq_primary+1, neq

         nposition = nposition + 1

c	The node number of the node that is lower in number
c	than node i is nelm(nelm(i)+1)

c	The position in the nelm array for the connection
c	from connected_node to i is nelm(connected_node+1)

c	The corresponding position of the connection from i to
c	connected_node is nelm(i)+1
         connected_node = nelm(nelm(i)+1)
         istrw(nelm(connected_node+1)-neq-1) = nposition
         istrw(nelm(i)-neq) = nposition
c  identify original gridblock volume         
         sx1save =sx1(connected_node) +sx1(i)
c     Compute area/distance term, stor in sx array


         call gdpm_area

      end do
c
c correct fracture connections by volume fraction
c this is only valid for model 11 because
c vfrac_primary is the width of the fracture
c probably have to redefine the whole sx and istrw arrays 
c (this will enlarge the array in most cases)
c scale the sx (area/d) by the primary 
c
      if(gdpm_flag.eq.11) then 
c find neighbors in x direction       
       call gdkm_connect(1) 
c adjust spacing and connections   
       call gdkm_connect(2)        
       i1 = nelm(neq+1) + 1
       if(isoy.ne.1) then
        allocate(sx_dum(i1,3))
        allocate(istrw_new(i1-neq-1))
       else
        allocate(sx_dum(i1-neq-1,1))
        allocate(istrw_new(i1-neq-1))
       endif
        sx_dum = 0.0
        istrw_new = 0   
        iw_c = 0 
        do i = 1,neq
        if(i.le.neq_primary) then
         imodel = igdpm(i)
         if(imodel.gt.0) then
          vfac = vfrac_primary(imodel)
         else
          vfac = 1.0
         endif
        else
         imodel = 0
         vfac = 1.0
        endif
         i1 = nelmdg(i)+1     
         i2 = nelm(i+1)
          do jj = i1,i2
           kb = nelm(jj)
           iw = istrw(jj-neq-1)
           if(kb.le.neq_primary) then
c primariy  to primary connection         
            imodel = igdpm(kb)
            if(imodel.gt.0) then
             vfackb = vfrac_primary(imodel)
            else
             vfackb = 1.0
            endif  
            vfacmin = min(vfac,vfackb)          
           else
            imodel = 0
            vfacmin = 1.
           endif
           iw_c = iw_c + 1
           istrw_new(jj-neq-1) = iw_c
           sx_1 = sx(iw,isox)
           sx_2 = sx(iw,isoy)
           sx_3 = sx(iw,isoz)
           sx_dum(iw_c,isox) = vfacmin*sx_1 
           sx_dum(iw_c,isoy) = vfacmin*sx_2
           sx_dum(iw_c,isoz) = vfacmin*sx_3           
          enddo
         enddo
      nr = iw_c +1
      deallocate(sx)
      if(isoy.ne.1) then
       allocate(sx(nr,3)) 
       do i = 1,nr
        sx(i,1) = sx_dum(i,1)
        sx(i,2) = sx_dum(i,2)
        sx(i,3) = sx_dum(i,3)
       enddo
      else
       allocate(sx(nr,3))
       do i = 1,nr
        sx(i,1) = sx_dum(i,1)
       enddo       
      endif
      i1 = nelm(neq+1) -(neq+1)
      do i = 1,i1
       istrw(i) = istrw_new(i)
      enddo 
      deallocate(sx_primary)
      deallocate(istrw_new)
      deallocate(sx_dum)
      call gdkm_connect(-1)
      endif
c    deallocate memory      
      if(allocated(vol_gdpm)) deallocate(vol_gdpm)
      if(gdkm_flag.eq.0) then
       if(allocated(gdpm_x)) deallocate(gdpm_x)
      endif
      if(allocated(gdpm_vol)) deallocate(gdpm_vol)
      if(allocated(vfrac_primary).and.gdkm_flag.eq.0) 
     &  deallocate(vfrac_primary)
      if(allocated(wgt_length)) deallocate(wgt_length)
      
      return

      contains

      subroutine gdpm_volume
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To compute control volumes for gdpm nodes..
!D1
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.2 Finite-Element Coefficient Generation
!D3 2.4.8 Generalized dual-porosity formulation
!D3
!**********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      implicit none
      real*8 delx_low, delx_high, x_centered
      real*8 r_low, r_high
      real*8 frac_num
      if(gdpm_flag.eq.6.or.gdpm_flag.eq.4) then

c     Linear 1-d column model (based on total area) 
c     Input distances are assumed block centered
c     GAZ 4-20-2001

c     Volume calculated by difference between adjacent nodes
c     times the Area which is read in palace of the primary fraction
c     vfrac_primary(imodel) contains the area
c     assumes primary node starts at x=0

         if(j.eq.1.and.ngdpm_layers(imodel).eq.1) then
          sx1(n_n_n) = gdpm_x(imodel,1)* 
     2                 vfrac_primary(imodel)*area_factor         
         else if(j.eq.1) then
          sx1(n_n_n) = (gdpm_x(imodel,j+1)+gdpm_x(imodel,j))/2.*
     2                 vfrac_primary(imodel)*area_factor
         else if(j.eq.ngdpm_layers(imodel)) then
c     assumes  last node in model is at end of the column
          sx1(n_n_n) = (gdpm_x(imodel,j)-gdpm_x(imodel,j-1))*
     2                 vfrac_primary(imodel)*area_factor
         else
          sx1(n_n_n) = (gdpm_x(imodel,j+1)-gdpm_x(imodel,j-1))/2.*
     2                 vfrac_primary(imodel)*area_factor
         endif
      

      else if(gdpm_flag.eq.5.or.gdpm_flag.eq.3) then

c     Linear 1-d column model
c     Input distances are assumed block centered
c     GAZ 4-20-2001

c     Volume calculated by difference between adjacent nodes
c     times the Area which is read in palace of the primary fraction
c     vfrac_primary(imodel) contains the area

         if(j.eq.1) then
c     assumes primary node starts at x=0
          sx1(n_n_n) = (gdpm_x(imodel,j+1)+gdpm_x(imodel,j))/2.*
     2                 vfrac_primary(imodel)
         else if(j.eq.ngdpm_layers(imodel)) then
c     assumes  last node in model is at end of the column
          sx1(n_n_n) = (gdpm_x(imodel,j)-gdpm_x(imodel,j-1))*
     2                 vfrac_primary(imodel)
         else
          sx1(n_n_n) = (gdpm_x(imodel,j+1)-gdpm_x(imodel,j-1))/2.*
     2                 vfrac_primary(imodel)
         endif
      
      else if(gdpm_flag.eq.2) then


c     Spherical model, diffusion from outer surface
c     into the sphere

c     Determine the radial positions of the point half-way between
c     the current node and the ones on either side of it


         if(j.eq.1) then
c     First node, coordinate position is outermost radius of sphere
            r_high = gdpm_x(imodel,ngdpm_layers(imodel))
         else
c     First determine average distance from outer surface,
c     then get radius by subtracting from outermost radius
            r_high = 0.5*
     2           (gdpm_x(imodel,j) + gdpm_x(imodel,j-1))
            r_high = gdpm_x(imodel,ngdpm_layers(imodel)) - r_high
         end if

         if(j.eq.ngdpm_layers(imodel)) then
c     Last node, no delx_high
            r_low = 0.
         else
c     Same type of calculation as r_high
            r_low = 0.5*
     2           (gdpm_x(imodel,j+1)+gdpm_x(imodel,j))
            r_low = gdpm_x(imodel,ngdpm_layers(imodel)) - r_low
         end if

c     Volume of cell is total volume of GDPM nodes times
c     ratio of shell volume to total sphere volume

         sx1(n_n_n) = sx1gdpm*(r_high**3-r_low**3)/
     2        gdpm_x(imodel,ngdpm_layers(imodel))**3
     
      else  if(gdpm_flag.eq.11) then


c     Physical Fracture Model (x direction)
    
c      length-weighted fractions of total volume

         sx1(n_n_n) = gdpm_vol(imodel,j)*sx1save
        continue
  
      else if(gdpm_flag.eq.1) then


c     Parallel Fracture Model


c     Determine half-lengths between adjacent nodes
         if(j.eq.1) then
c     First node, coordinate position is 0, use entire length
            delx_low = gdpm_x(imodel,1)
         else
            delx_low = 0.5*
     2           (gdpm_x(imodel,j) - gdpm_x(imodel,j-1))
         end if
         if(j.eq.ngdpm_layers(imodel)) then
c     Last node, no delx_high
            delx_high = 0.
         else
            delx_high = 0.5*
     2           (gdpm_x(imodel,j+1)-gdpm_x(imodel,j))
         end if
         
c     Compute gdpm node volume
         x_centered = delx_high+delx_low
         sx1(n_n_n) = sx1gdpm*x_centered/length_total
         
      end if

      end subroutine gdpm_volume


c     Subroutine for determining area/distance term
c     for GDPM connections

      subroutine gdpm_area
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To compute area/distance terms for gdpm connections.
!D1
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.2 Finite-Element Coefficient Generation
!D3 2.4.8 Generalized dual-porosity formulation
!D3
!**********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      implicit none
      real*8 sumr_gdpm, r_high, r_low, delr_gdpm
      real*8 frac_num, tolf, frac_num_dpdp
      parameter (tolf= 1.d-12)
      integer int, itest
      parameter( itest =7 )
      

c     Linear 1-D column model (based on total area)

      if(gdpm_flag.eq.6.or.gdpm_flag.eq.4) then


c     Determine if the connected node is a primary node, if
c     it is, change the primary node to this one

         if(connected_node.le.neq_primary) then
            primary_node = connected_node
            imodel = igdpm(primary_node)
            area_factor = sx1(primary_node)/vol_gdpm(imodel)
            gdpm_counter = 1
            delx_gdpm = gdpm_x(imodel,1)
         else
            gdpm_counter = gdpm_counter + 1
            delx_gdpm = gdpm_x(imodel,gdpm_counter)
     2           -gdpm_x(imodel,gdpm_counter-1)
         end if
         
c     Compute sx, stor in sx(iposition,isox)
         sx(nposition,isox) = 
     2    -vfrac_primary(imodel)*area_factor/delx_gdpm
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if

c     Linear 1-D column model

      else if(gdpm_flag.eq.5.or.gdpm_flag.eq.3) then


c     Determine if the connected node is a primary node, if
c     it is, change the primary node to this one

         if(connected_node.le.neq_primary) then
            primary_node = connected_node
            imodel = igdpm(primary_node)
            gdpm_counter = 1
            delx_gdpm = gdpm_x(imodel,1)
         else
            gdpm_counter = gdpm_counter + 1
            delx_gdpm = gdpm_x(imodel,gdpm_counter)
     2           -gdpm_x(imodel,gdpm_counter-1)
         end if
         
c     Compute sx, stor in sx(iposition,isox)
         sx(nposition,isox) = -vfrac_primary(imodel)/delx_gdpm
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if

c     Spherical model

      else if(gdpm_flag.eq.2) then

c     Determine if the connected node is a primary node, if
c     it is, change the primary node to this one

         if(connected_node.le.neq_primary) then
            primary_node = connected_node
            imodel = igdpm(primary_node)
c     Need to compute gdpm volume
c     sx1(primary_node) is the volume of only the 
c     primary portion
            sx1gdpm = sx1(primary_node)*(1.-vfrac_primary(imodel))
     2           /vfrac_primary(imodel)
c            sx1gdpm = sx1save*(1.-vfrac_primary(imodel))
            gdpm_counter = 1
            r_high = gdpm_x(imodel,ngdpm_layers(imodel))
            r_low = gdpm_x(imodel,ngdpm_layers(imodel)) -
     2           gdpm_x(imodel,1)
            sumr_gdpm = r_low + r_high
            delr_gdpm = r_high - r_low
         else
            gdpm_counter = gdpm_counter + 1
            r_low = gdpm_x(imodel,ngdpm_layers(imodel)) -
     2           gdpm_x(imodel,gdpm_counter)
            r_high = gdpm_x(imodel,ngdpm_layers(imodel)) -
     2           gdpm_x(imodel,gdpm_counter-1)
            sumr_gdpm = r_low + r_high
            delr_gdpm = r_high - r_low
         end if
         
c     Compute sx, stor in sx(iposition,isox)
         sx(nposition,isox) = -0.75*sumr_gdpm**2*sx1gdpm/
     2        (delr_gdpm*gdpm_x(imodel,ngdpm_layers(imodel))**3)
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if


c     Parallel fracture model for gdkm
      else if(gdkm_flag.eq.1) then
c     retired :Physical fracture model (model 11)
c gaz 090816 need to retire this model (replaced by gdkm models
c 1,2,3, 0
          if(connected_node.le.neq_primary) then
            primary_node = connected_node
c gaz 120822 frac number now input          
            imodel = igdpm(primary_node) 
             if(gdkm_dir(imodel).eq.1) then
              length_total = dxrg(primary_node)
              frac_num = gdpm_x(imodel,1)
             elseif(gdkm_dir(imodel).eq.2) then
              length_total = dyrg(primary_node) 
              frac_num = gdpm_x(imodel,1)
             elseif(gdkm_dir(imodel).eq.3) then 
              length_total = dzrg(primary_node)
              frac_num = gdpm_x(imodel,1)
             elseif(gdkm_dir(imodel).eq.4) then
c generic fracture model  with num of fractures/gridblock              
                 if(icnl.eq.0) then
                  length_total = sx1save**(0.333333)
                 else
                  length_total = sx1save**(0.5)
                 endif
                delx_gdpm =  gdpm_x(imodel,1)
                frac_num = max(int(length_total/delx_gdpm+tolf),1)
                continue

              elseif(gdkm_dir(imodel).eq.0) then                                                 
c generic fracture model original  
                 dis_p = 0.0
                 dis_s = 0.0
                 if(icnl.eq.0) then
                  length_total = sx1save**(0.333333)
                 else
                  length_total = sx1save**(0.5)
                 endif
                 delx_gdpm = gdpm_x(imodel,1)  
                 frac_num = 0
             endif
            if(gdkm_dir(imodel).ne.0) then
             length_pri = vfrac_primary(imodel) 

c             frac_num = delx_gdpm  gaz 031317
c gaz 120822              
c            frac_num_dpdp = max(int(length_total/delx_gdpm+tolf),1)
c A/d for primary-secondary connection is
c Ax = Volume_total/gridblock_length in x dir (Ax = sx1save/length_total)          
c Compute sx, stor in sx(iposition,isox) (- sign for FEHM sign convention)
c 1/2 distance primary volume(dis_p) =   length_pri*length_total/2
c 1/2 distance secondary volume(dis_s)  =   (1.-length_pri)*length_total/2    
c assuming sec volume is centered, divide dis_p again by 2 for the center(2*2 = 4)   
c delx_gdpm = gdpm_x(imodel,1) not used as frac spacing (inputs fracture number)

             dis_p =   length_pri*length_total/((frac_num+1)*2)

             dis_s =  (1.-length_pri)*length_total/(frac_num+1)            
c gaz 121322 distance does not need dis_s
             dis_s = 0.0
c multiply by 2 because centered secondary volume has 2 area faces 
c multiply by frac_num because of added surface area         
c this is affected by number of fractures present (length_total/gdpm_x(imodel,1))    
c gaz 042218 
             sx(nposition,isox) = -2.*(sx1save/length_total)/
     &                      (dis_p+dis_s)  
          continue
          if(primary_node.eq.itest) write(ierr,901)  
          write(ierr,902) primary_node,sx1( primary_node), dis_p,dis_s,
     &    nposition, sx(nposition,isox),gdkm_dir(imodel),
     &    frac_num, gdpm_x(imodel,1)    
901       format (t1,'primary_node',t15, 'sx1(pri)',
     &    t30,'dis_p',t45,'dis_s',t59,'nposition',t72,'sx(npos,isox)')
902     format (t1,i10,t15,g12.4,t29,g12.4,t43,g12.4,t59,i10,t72,g12.4,
     &    ' gdkm_dir', i4,' frac_num ', f8.3,' delx_gdpm ', f9.3)
          else                                                                                
c older generic fracture model - like dpdp
           sx1gdpm = sx1( primary_node)*(1.-vfrac_primary(imodel)) 
     &           /vfrac_primary(imodel)
c           sx1gdpm = sx1(primary_node)/vfrac_primary(imodel)
c           sx1gdpm = sx1(primary_node)*(vfrac_primary(imodel))
c     &           /(1.-vfrac_primary(imodel))          
c            note that length_total and delx_gdpm are the same (gives dpdp result)
            length_total = gdpm_x(imodel,ngdpm_layers(imodel)) 
            gdpm_counter = 1
c  use older model          
            delx_gdpm = gdpm_x(imodel,1)
c            delx_gdpm = gdpm_x(imodel,1)/2.0
c gaz debug 031718  
c need to use full volume "sx1(primary_node)/vfrac_primary(imodel)" for symmetry            
c            sx(nposition,isox) = -sx1gdpm/(delx_gdpm*length_total)          
c        sx(nposition,isox) = -(sx1(primary_node)/vfrac_primary(imodel))
c    &       /(delx_gdpm*length_total)
         sx(nposition,isox) = -sx1save/(delx_gdpm*length_total)

          write(ierr,902) primary_node,sx1( primary_node), dis_p,dis_s,
     &    nposition, sx(nposition,isox), gdkm_dir(imodel),frac_num,
     &    gdpm_x(imodel,1)
          endif
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if
          endif  
      else if(gdpm_flag.eq.11) then
         if(connected_node.le.neq_primary) then
            primary_node = connected_node
            imodel = igdpm(primary_node)
c     Need to compute gdpm volume
c     sx1(primary_node) is the volume of only the 
c     primary portion
c     length_total is length of secondary material
c     length_pri is length of primary material
c            sx1gdpm = sx1(primary_node)*(1.-vfrac_primary(imodel))
c     2           /vfrac_primary(imodel)
 
             length_total = dxrg(primary_node)
             length_pri = vfrac_primary(imodel)
            
            gdpm_counter = 1
            delx_gdpm = gdpm_x(imodel,1)*length_total+0.5*length_pri
         else
            gdpm_counter = gdpm_counter + 1
            delx_gdpm = (gdpm_x(imodel,gdpm_counter)
     2           -gdpm_x(imodel,gdpm_counter-1))*length_total
         end if
         
c     Compute sx, stor in sx(iposition,isox)
         sx(nposition,isox) = -sx1save/(delx_gdpm*length_total)
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if


c     Parallel fracture model

      else if(gdpm_flag.eq.1) then

c     Determine if the connected node is a primary node, if
c     it is, change the primary node to this one

         if(connected_node.le.neq_primary) then
            primary_node = connected_node
            imodel = igdpm(primary_node)
c     Need to compute gdpm volume
c     sx1(primary_node) is the volume of only the 
c     primary portion
            sx1gdpm = sx1(primary_node)*(1.-vfrac_primary(imodel))
     2           /vfrac_primary(imodel)
c            sx1gdpm = sx1save*(1.-vfrac_primary(imodel))
            length_total = gdpm_x(imodel,ngdpm_layers(imodel))
            gdpm_counter = 1
            delx_gdpm = gdpm_x(imodel,1)
         else
            gdpm_counter = gdpm_counter + 1
            delx_gdpm = gdpm_x(imodel,gdpm_counter)
     2           -gdpm_x(imodel,gdpm_counter-1)
         end if        
         
c     Compute sx, stor in sx(iposition,isox)         
         sx(nposition,isox) = -sx1gdpm/(delx_gdpm*length_total)
         if(isoy.eq.1) then
            if(icnl.eq.0) then
               sx(nposition,isox) = sx(nposition,isox)/3.
            else
               sx(nposition,isox) = sx(nposition,isox)/2.
            end if
         else
            sx(nposition,isoy) = 0.
            sx(nposition,isoz) = 0.
         end if

      end if

      end subroutine gdpm_area

      end subroutine add_gdpm

