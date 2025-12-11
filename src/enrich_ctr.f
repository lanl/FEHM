       subroutine enrich_ctr(iflg)
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
!!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 enrich the grid in defined zones and in coordinate directions
!D1
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 29-Nov-09, Programmer: George Zyvoloski
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 3.0
!D4 
!**********************************************************************
c NOTES gaz 08010
c need to call scanin to count models etc.
c need to skip these elements and nodes in anonp etc
c direction follow local coordinates , not global coordinates
c At present the formulation is only for hexahedral elements
c Could use this directly with anisotropic K 
c First generate all nodes 
c check for repeated nodes then generate elements

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai

      implicit none
c    
      integer, allocatable :: lrp(:) 
      integer, allocatable :: nrp(:)  
      integer i,j,ii,ij,jj,kb,ie,iflg,iieosd
      integer node_add,nei_add,ietype
      integer nlayers,nrepeat
      integer max_iter,imodel
      character*4 vard
      character*80 dummy_string
      logical null1
      real*8 dis_tol
      parameter (max_iter = 10, dis_tol = 1.e-8)
c      
      if(ienrich.eq.0) return
c      
      if(iflg.eq.0) then
 1000 continue
      read(inpt,'(a80)') dummy_string
      if(.not. null1(dummy_string)) then
         imodel = imodel + 1
         backspace inpt
         read(inpt,*) ienrich_dir(imodel),ienrich_layers(imodel)
         goto 1000
      end if      
c
c input similiar to gdpm
c
c     Set flag to identify which nodes with enrichment model
      narrays = 1
      itype(1) = 4
      default(1) = 0
      igroup = 2
      macro = 'enri'

      call initdata2 (inpt, ischk, neq_primary, narrays, itype, 
     *     default, macroread(10), macro, igroup, ireturn,
     *     i4_1=ienrich_model(1:neq_primary)) 
c       
      else if(iflg.eq.1) then       
c
c  identify elements connected to enrichment nodes
c
      if(ns.ne.8) then
c
c error condition
c
      if(iptty.ne.0) write(iptty,*) 
     & 'stopping in enrich_ctr because ns ne 8'
      endif
      allocate (ipl_enrich(n0),nop_enrich(n0,8))
      allocate (nei_enrich_id(nei))
      nei_enrich_id = 0
      ipl_enrich=0
      do ie=1,nei
        do j=1,ns
           i=nelm((ie-1)*ns+j)
           if (ienrich_model(i).ne.0) then
              nei_enrich_id(ie) = ienrich_model(i)
              ipl_enrich(i)=ipl_enrich(i)+1
              nop_enrich(i,ipl_enrich(i))=ie
              noo_enrich(i,ipl_enrich(i))=j 
           endif
        enddo
      enddo
c      
      else if(iflg.eq.2) then  
c      
c generate additional resolution in identified elements     
c      
      node_add = 0
      nei_add = 0
      do ie = 1,nei
       if(nei_enrich_id(ie).ne.0) then
        imodel = nei_enrich_id(ie)
        ietype = ienrich_dir(imodel)
        nei_enrich = nei_enrich + 1
        nei_enrich_list(nei_enrich) = ie
        if(ietype.eq.1) then
c enrichment in the local x direction
         nlayers = ienrich_layers(imodel)
         call add_resolution(1,ie,itype,nlayers,nei_add,node_add)          
        endif
       endif
      enddo      
c eliminate duplicate nodes (this could be very slow)
      nrepeat = 0
      do ii = neq_primary + 1,neq_primary + node_add
       do jj =  ii + 1, neq_primary + node_add
        if(abs(cord(jj,1)-cord(ii,1)).lt.dis_tol) then
         if(abs(cord(jj,2)-cord(ii,2)).lt.dis_tol) then
          if(abs(cord(jj,3)-cord(ii,3)).lt.dis_tol) then
            nrepeat = nrepeat + 1
            lrp(nrepeat) = jj
            nrp(nrepeat) = ii
          endif
         endif
        endif
      enddo
      enddo
      
      else if(iflg.eq.-1) then
c      
c print out error condition      
c
               if(ii .gt. max_iter) then
                  if(iout.ne.0) 
     *                 write(iout,*) 'failed in enrich_ctr'
                  if(iptty.ne.0)
     *                 write(iptty,*) 'failed in enrich_ctr'
               end if
       
      endif   
      return 
      end
      subroutine add_resolution(iflg,ie,imtype,nlayers,nei_add,node_add)
c
c add directional resolution to an element
c
      use combi
      use comdti     
      use comki 
      use comai
      implicit none
      integer iflg,ie,ii,nlayers,nei_add,node_add
c      integer n1,n2,n3,n4,n5,n6,n7,n8
      integer jj,node_add0,new_tot,imtype
      integer n1,nei_t(8)
      real*8 delx,dely,delz
      real*8 x1,y1,z1,x2,y2,z2
c      
      if(iflg.eq.1) then
       if(imtype.eq.1) then
         node_add0 = node_add
         new_tot = (nlayers-1)*4
         ii=(ie-1)*ns
         n1 = nelm(ii + 1)
         n2 = nelm(ii + 2)
         n3 = nelm(ii + 3)
         n4 = nelm(ii + 4)  
         n5 = nelm(ii + 5)
         n6 = nelm(ii + 6)
         n7 = nelm(ii + 7)
         n8 = nelm(ii + 8) 
c         
c  identify x pairs
c  generate new nodes and elements 
c n1-n2
         x1 = cord(n1,1)
         y1 = cord(n1,2)
         z1 = cord(n1,3)
         x2 = cord(n2,1)
         y2 = cord(n2,2)
         z2 = cord(n2,3)   
         delx = (x2-x1)/nlayers
         dely = (y2-y1)/nlayers
         delz = (z2-z1)/nlayers
         do jj = 1, nlayers-1
          node_add = node_add+1
          cord(node_add,1) = x1 + jj*delx
          cord(node_add,2) = y1 + jj*dely
          cord(node_add,3) = z1 + jj*delz
         enddo    
c n4-n3
         x1 = cord(n4,1)
         y1 = cord(n4,2)
         z1 = cord(n4,3)
         x2 = cord(n3,1)
         y2 = cord(n3,2)
         z2 = cord(n3,3)   
         delx = (x2-x1)/nlayers
         dely = (y2-y1)/nlayers
         delz = (z2-z1)/nlayers
         do jj = 1, nlayers-1
          node_add = node_add+1
          cord(node_add,1) = x1 + jj*delx
          cord(node_add,2) = y1 + jj*dely
          cord(node_add,3) = z1 + jj*delz
         enddo  
c n5-n6
         x1 = cord(n5,1)
         y1 = cord(n5,2)
         z1 = cord(n5,3)
         x2 = cord(n6,1)
         y2 = cord(n6,2)
         z2 = cord(n6,3)   
         delx = (x2-x1)/nlayers
         dely = (y2-y1)/nlayers
         delz = (z2-z1)/nlayers
         do jj = 1, nlayers-1
          node_add = node_add+1
          cord(node_add,1) = x1 + jj*delx
          cord(node_add,2) = y1 + jj*dely
          cord(node_add,3) = z1 + jj*delz
         enddo
c n8-n7
         x1 = cord(n8,1)
         y1 = cord(n8,2)
         z1 = cord(n8,3)
         x2 = cord(n7,1)
         y2 = cord(n7,2)
         z2 = cord(n7,3)   
         delx = (x2-x1)/nlayers
         dely = (y2-y1)/nlayers
         delz = (z2-z1)/nlayers
         do jj = 1, nlayers-1
          node_add = node_add+1
          cord(node_add,1) = x1 + jj*delx
          cord(node_add,2) = y1 + jj*dely
          cord(node_add,3) = z1 + jj*delz
         enddo         
c
c new coordinates generated
c now generate new elements
c remember to replace old element with first new element
c last generated element will have 'old' nodes
c
         nei_t(1) = n1 
         nei_t(2) = node_add0 + 1
         nei_t(3) = n4
         nei_t(4) = node_add0 + nlayers-1 + 1
         nei_t(5) = n5
         nei_t(6) = node_add0 + 2*(nlayers-1) + 1 
         nei_t(7) = n8 
         nei_t(8) = node_add0 + 3*(nlayers-1) + 1 
c replace original element with first new element         
         nelm(ii+1:ii+8) = nei_t(1:8)
c         
         do jj = 1,nlayers-2
          nei_add = nei_add + 1
	  ii=(nei_add-1)*ns
	  nelm(ii+1) = nei_t(2)
	  nelm(ii+4) = nei_t(3)
	  nelm(ii+5) = nei_t(6)
	  nelm(ii+8) = nei_t(7)
	  nelm(ii+2) = nelm(ii+1) + 1
	  nelm(ii+3) = nelm(ii+4) + 1 + (nlayers-1)
	  nelm(ii+6) = nelm(ii+5) + 1 + 2*(nlayers-1)
	  nelm(ii+7) = nelm(ii+8) + 1 + 3*(nlayers-1) 
	  nei_t(2) = nelm(ii+2)
	  nei_t(3) = nelm(ii+3)
	  nei_t(6) = nelm(ii+6)
	  nei_t(7) = nelm(ii+7)	  
         enddo 
c now do last new element separetly
          nei_add = nei_add + 1
	  ii=(nei_add-1)*ns
	  nelm(ii+1) = nei_t(2)
	  nelm(ii+4) = nei_t(3)
	  nelm(ii+5) = nei_t(6)
	  nelm(ii+8) = nei_t(7)
	  nelm(ii+2) = n2
	  nelm(ii+3) = n3
	  nelm(ii+6) = n6
	  nelm(ii+7) = n7
       endif  
      endif
      return
      end
