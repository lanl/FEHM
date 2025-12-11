      subroutine setup_cliffnodes(icliff1,icliff2,icliff3)
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
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
!D1 Setup cliff nodes.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/cliffnodes.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comsk
      implicit none

      integer cliff_count, icliff1, icliff2, icliff3, istep_count

c     s kelkar march 30, 05 generate list of cliff nodes
      if (icliff1 .ne. 0) then
         call generate_cliff_list(cliff_count, icliff1, icliff2)
         if(cliff_count.gt.0) then
            call cliff_setflag(cliff_count, istep_count)
         endif
         call write_cliff (icliff3, cliff_count, istep_count)
      else 
         call read_cliff (icliff3)
      end if
      if (allocated(list_case)) deallocate (list_case)
      if (allocated(node_save)) deallocate (node_save)
      if (allocated(node_flag1)) deallocate (node_flag1)
      if (allocated(ncase)) deallocate (ncase)
      if (allocated(node_global)) deallocate (node_global)
      if (allocated(cords)) deallocate (cords)
      if (allocated(ne)) deallocate (ne)
      if (allocated(rnodecol)) deallocate (rnodecol)
      if (allocated(relemcol)) deallocate (relemcol)
      if (allocated(tri_med)) deallocate (tri_med)
      if (allocated(tri_norm)) deallocate (tri_norm)
      if (allocated(jtet)) deallocate (jtet)
      if (allocated(jtetoff)) deallocate (jtetoff)
      if (allocated(edge_node_tri)) deallocate (edge_node_tri)

      return

      end subroutine setup_cliffnodes

c......................................................................

      subroutine read_surface_mesh(icliff1, icliff2, nedges, num_node, 
     &     num_element)

      use comsk
      use comai, only : ierr, iptty 

      implicit none

      character*4 dumm
      character*5 dumm5
      character*7 dumm7
      character*14 dumm14
      character*17 dumm17

      character*160 line

      integer maxdim,maxcorner,num_node,num_element
      integer nelem_att,node_att2, nelem_att2,node_att,ii,node
      integer nelm,it,jt,ntets,nelmnef,mbndry,nef_cmo,nedges,i
      integer jf,itemp,n1,j1,jj,junk,i1,iii,inodecol,idum(100)
      integer icolxmed,icolymed,icolzmed
      integer icolxnorm,icolynorm,icolznorm
      integer icliff1,icliff2

      real*8 elem,rnode,xmed,ymed,zmed,xnorm,ynorm,znorm
      real*8 xn1,yn1,zn1,xn2,yn2,zn2,pi,cj,di,pj,ptemp

      maxdim=1000000
      maxcorner=8

c     read node coordinates
      read(icliff1,*) num_node, num_element, node_att, nelem_att,idum(1)
      allocate (cords(num_node,3))
      node_att2=1
      nelem_att2=0
      do ii=1,num_node
         read(icliff1,*)node,(cords(node,jj),jj=1,3)
      enddo
c     read triangles defined by 3 nodes
      allocate (ne(num_element,3))
      do ii=1,num_element
         read(icliff1,*)nelm,junk,dumm,(ne(nelm,jj),jj=1,3)
      enddo
c     skip line for formatting
      read(icliff1,*)
      idum=1
c identify column with idnode0 attribute
      do ii=1,node_att
         read(icliff1,4147) dumm7
 4147    format(a7)
         if(dumm7.eq.'idnode0') inodecol=ii
      enddo
c     read the node numbering in the original 3-D mesh corrosponding
c     to the surface numbering.
      allocate (rnodecol(node_att), node_global(num_node))
      do ii=1,num_node
         read(icliff1,*)node,(rnodecol(iii),iii=1,node_att)
         node_global(node)=int(rnodecol(inodecol))
      enddo
c     skip line for formatting
      read(icliff1,*)
      allocate (relemcol(nelem_att))
      allocate (tri_med(num_element,3),tri_norm(num_element,3))
c identify columns with needed element attributes attribute
      do ii=1,nelem_att
         read(icliff1,4142) dumm5
 4142    format(a5)
         if(dumm5.eq.'xmed,') icolxmed=ii
         if(dumm5.eq.'ymed,') icolymed=ii
         if(dumm5.eq.'zmed,') icolzmed=ii
         if(dumm5.eq.'xnorm') icolxnorm=ii
         if(dumm5.eq.'ynorm') icolynorm=ii
         if(dumm5.eq.'znorm') icolznorm=ii
      enddo
c     for each element, read cords of the median point and components
c     of the outward normal
      do ii=1,num_element
         read(icliff1,*)nelm,(relemcol(iii),iii=1,nelem_att)
         tri_norm(nelm,1)=relemcol(icolxnorm)
         tri_norm(nelm,2)=relemcol(icolynorm)
         tri_norm(nelm,3)=relemcol(icolznorm )
         tri_med(nelm,1)=relemcol(icolxmed)
         tri_med(nelm,2)=relemcol(icolymed)
         tri_med(nelm,3)=relemcol(icolzmed)
      enddo

c     read jtetoff, pointer to the element edge information sotred in
c     the array jtet
c     read(icliff2,*)
 4141 format(a160)
      allocate (jtetoff(num_element),jtet(num_element*3))
      dumm17='Attribute: jtetoff'
      do
         read(icliff2,4141,end=9111)line
         if(line(1:17).eq.dumm17(1:17)) exit
      enddo
      do ii=1,num_element
         read(icliff2,*)it,jt
         jtetoff(it)=jt
      enddo

      dumm14='Attribute: jtet'
      do 
         read(icliff2,4141,end=9112)line
         if(line(1:14).eq.dumm14(1:14)) exit
      enddo
      do ii=1,3*num_element
         read(icliff2,*)it,jt
         jtet(ii)=jt
      enddo

      ntets=num_element
      nelmnef=3
      mbndry=0
      nef_cmo=3
      nedges=0
      allocate(edge_node_tri(int(num_element*1.5),4))
      do it=1,ntets
         do i=1,nelmnef
c     check if element face is on an external boundry
            if(jtet(jtetoff(it)+i).eq.mbndry) then
               jt=0
               jf=0
c     check if element face is on an internal boundry
            elseif(jtet(jtetoff(it)+i).gt.mbndry)then
               jt=1+(jtet(jtetoff(it)+i)-mbndry-1)/nef_cmo
               jf=jtet(jtetoff(it)+i)-mbndry-nef_cmo*(jt-1)
C     Volume element
            else
               jt=1+(jtet(jtetoff(it)+i)-1)/nef_cmo
               jf=jtet(jtetoff(it)+i)-nef_cmo*(jt-1)
            endif
c     it is assumed that jtet is organized in order of non-decreasing
c     values of it. Note that every triangle has 3 edges, and every edge 
c     has 2 triangls associated with it, because the the surface mesh
c     was created as a wrap-around of a 3-D object, it has no external
c     edges. Also, the way jtet is constructed, the various values of jt
c     corresponding to an it are all distinct.
            if(it.lt.jt) then
               nedges=nedges+1
               edge_node_tri(nedges,3)=it
               edge_node_tri(nedges,4)=jt
c     find the two nodes common to the triangles it and jt
               itemp=0
               do i1=1,3
                  n1=ne(it,i1)
                  do j1=1,3
                     if(n1.eq.ne(jt,j1)) then
                        itemp=itemp+1
                        edge_node_tri(nedges,itemp)=n1
                     endif
                  enddo
               enddo
            endif 
            
         enddo
      enddo

      close (icliff1)
      close (icliff2)
      return
 9111 write(ierr,*)'wrong format in triangle2.jtet'
      write(ierr,*)' cant find string "Attribute: jtetoff"'
      write(ierr,*)'STOP.'
      if (iptty .ne. 0) then
         write(iptty,*)'wrong format in triangle2.jtet'
         write(iptty,*)' cant find string "Attribute: jtetoff"'
         write(iptty,*)'STOP.'
      end if
      stop
 9112 write(ierr,*)'wrong format in triangle2.jtet'
      write(ierr,*)' cant find string "Attribute: jtet"'
      write(ierr,*)'STOP.'
      if (iptty .ne. 0) then
         write(iptty,*)'wrong format in triangle2.jtet'
         write(iptty,*)' cant find string "Attribute: jtet"'
         write(iptty,*)'STOP.'
      end if
      stop

      end subroutine read_surface_mesh
c......................................................................

      subroutine generate_cliff_list(cliff_count, icliff1, icliff2)
c     s kelkar Mar 22, 2005
c     using a triangular mesh representing the surface nodes of a 3-D 
c     mesh and the jtet array for this triangular mesh, generate a list 
c     of cliff nodes that have an internal angle at external faces.

      use comsk
      implicit none

      integer num_node, num_element, node_att, nelem_att,maxdim
      integer node,junk,maxcorner, node_att2, nelem_att2
      integer mbndry,nef_cmo,nedges,icliff1, icliff2
      integer itemp,n1,i1,j1,cliff_count,ie,cliffcode
      integer node1,node2,tri1,tri2,tritemp,i,j,k,ii,jj,nelm,iii
      integer it,jt,ntets,nelmnef,nodeg,ncs,jf,list_code
      integer idum(100)

      real*8 elem,rjunk,rnode,xmed,ymed,zmed,xnorm,ynorm,znorm
      real*8 xn1,yn1,zn1,xn2,yn2,zn2,pi,cj,di,pj,ptemp
      
      call read_surface_mesh(icliff1, icliff2, nedges, num_node, 
     &     num_element)
      
      allocate (node_flag1(num_node))
      allocate (ncase(num_node))
      allocate (list_case(num_node,8))
      allocate (node_save(num_node))

      ncase=0
      list_case=0
      cliff_count=0

      do ie=1,nedges

c     recall the nodes and triangle #s associated with the edge
         node1=edge_node_tri(ie,1)
         node2=edge_node_tri(ie,2)
         tri1=edge_node_tri(ie,3)
         tri2=edge_node_tri(ie,4)

c     recall the outward normals to the two triangles
         xn1=tri_norm(tri1,1)
         yn1=tri_norm(tri1,2)
         zn1=tri_norm(tri1,3)
         xn2=tri_norm(tri2,1)
         yn2=tri_norm(tri2,2)
         zn2=tri_norm(tri2,3)

c     find the axis index corresponding to the normals
         if(abs(xn1).eq.1.) then
            i=1
            pi=xn1
         elseif(abs(yn1).eq.1.) then
            i=2
            pi=yn1
         elseif(abs(zn1).eq.1.) then
            i=3
            pi=zn1
         endif

         if(abs(xn2).eq.1.) then
            j=1
            pj=xn2
         elseif(abs(yn2).eq.1.) then
            j=2
            pj=yn2
         elseif(abs(zn2).eq.1.) then
            j=3
            pj=zn2
         endif

         if(i.ne.j) then
c     if i and j in the same plane then its not a corner
c     find the 3rd axis, k so that (i,j,k)=(1,2,3) (not in any order)
            if(i*j.eq.2) then
               k=3
            elseif(i*j.eq.3) then
               k=2
            elseif(i*j.eq.6) then
               k=1
            endif

c     using the projections coordinates of the midpoints of the triangles
c     onto the i-j plane, if the point of intersection of the lines
c     thru the projected points and parallel to the respective projections
c     of the normals is in the same direction as the normals, then the 
c     point is outside and the grid has an internal angle a the edge
            do ii=1,2
               node=edge_node_tri(ie,ii)
               cj=tri_med(tri1,j)-cords(node,j)
               di=tri_med(tri2,i)-cords(node,i)
               if(pi*di.gt.0..and.pj*cj.gt.0.) then
                  if(node_flag1(node1).eq.0) then
                     cliff_count=cliff_count+1
                     node_flag1(node1)=1
                     node_save(cliff_count)=node1
                  endif
                  if(node_flag1(node2).eq.0) then
                     cliff_count=cliff_count+1
                     node_flag1(node2)=1
                     node_save(cliff_count)=node2
                  endif
                  call get_list_case(i,j,pi,pj,list_code)
                  do iii=1,ncase(node1)
                     if(list_case(node1,iii).eq.list_code) goto 9191
                  enddo
                  ncase(node1)=ncase(node1)+1
                  list_case(node1,ncase(node1))=list_code
 9191             continue
                  do iii=1,ncase(node2)
                     if(list_case(node2,iii).eq.list_code) goto 9192
                  enddo
                  ncase(node2)=ncase(node2)+1
                  list_case(node2,ncase(node2))=list_code
 9192             continue
               endif
            enddo
         endif
         
      enddo

      do ii=1,cliff_count
         node=node_save(ii)
         nodeg=node_global(node)
         ncs=ncase(node)
      enddo
      
      do ii=1,num_node
         do jj=1,cliff_count
            if(ii.eq.node_save(jj)) then
               cliffcode=ncase(jj)
               exit
            endif
         enddo
         cliffcode=0
      enddo


      return

      end subroutine generate_cliff_list

c.......................................................................

      subroutine  get_list_case(iin,jin,piin,pjin,list_code)

c     if i is a surface node on a stair-step on an internal corner 
c     of a 3D mesh
c     then it does not have a brickshape cc. create reflecting
c     boundaries for this. Such nodes are tagged in Tag_stairstep
c     as a part of the mesh generating. 
c     n_stairstep= # of such nodes. For these nodes,
c     for non-OMR node, irray(i,0)=-1e+8-icounter, 
c     (later, in subroutine , for OMR nodes, this is changed to 
c     irray(i,0)=-2e+8-icounter)
c     where icounter points to the
c     slot in the array istep_cases. istep_cases(icounter) is in turn
c     another pointer that points to the slot in the array istep_cases
c     where information for the node i is saved. simillar to nelm()
c     ip=-irray(i,0)-1e+8
c     i1=istep_cases(ip)+1
c     i2=istep_cases(ip+1)
c     ncases=i2-i1+1 = # of cliff cases intersecting at node i
c     istep_cases(i1:i1)= case # given below.
c     Each node can have upto 3 corners cut out from its cc. 
c     so 0<ncases<=3. These are classified according to
c     the direction of the outward normals defing the two faces of the 
c     cc that are cut out. There are 6*4/2=12 cases as follows:
c     Case  1 : normal1 = +x, normal2 = +y
c     Case  2 : normal1 = +x, normal2 = -y
c     Case  3 : normal1 = +x, normal2 = +z
c     Case  4 : normal1 = +x, normal2 = -z
c     Case  5 : normal1 = -x, normal2 = +y
c     Case  6 : normal1 = -x, normal2 = -y
c     Case  7 : normal1 = -x, normal2 = +z
c     Case  8 : normal1 = -x, normal2 = -z
c     Case  9 : normal1 = +y, normal2 = +z
c     Case 10 : normal1 = +y, normal2 = -z
c     Case 11 : normal1 = -y, normal2 = +z
c     Case 12 : normal1 = -y, normal2 = -z

      use comai, only : ierr
      use comsk
      implicit none

      integer i,j,list_code,itemp,jtemp,iin,jin

      real*8 pi,pj,ptemp,piin,pjin

c     arrange i<j
      if(iin.gt.jin) then
         i=jin
         pi=pjin
         j=iin
         pj=piin
      elseif(iin.lt.jin) then
         i=iin
         pi=piin
         j=jin
         pj=pjin
      else
         write(ierr,*)'iin=jin in get_list_case.STOP'
         stop
      endif

      if(i.eq.1) then
         if(pi.eq.+1.) then
c     now j =2 or 3 since 1<=i<j<=3
            if(j.eq.2) then
               if(pj.eq.+1.) then
                  list_code=1
               elseif(pj.eq.-1.) then
                  list_code=2
               endif
            elseif(j.eq.3) then
               if(pj.eq.+1.) then
                  list_code=3
               elseif(pj.eq.-1.) then
                  list_code=4
               endif
            endif
         elseif(pi.eq.-1.) then
            if(j.eq.2) then
               if(pj.eq.+1.) then
                  list_code=5
               elseif(pj.eq.-1.) then
                  list_code=6
               endif
            elseif(j.eq.3) then
               if(pj.eq.+1.) then
                  list_code=7
               elseif(pj.eq.-1.) then
                  list_code=8
               endif
            endif
         endif
      elseif(i.eq.2) then
c     now j must =3 since 1<=i<j<=3
         if(pi.eq.+1.) then
            if(pj.eq.+1.) then
               list_code=9
            elseif(pj.eq.-1.) then
               list_code=10
            endif
         elseif(pi.eq.-1.) then
            if(pj.eq.+1.) then
               list_code=11
            elseif(pj.eq.-1.) then
               list_code=12
            endif

         endif
      endif

      return

      end subroutine get_list_case

c................................................................

      subroutine cliff_setflag(cliff_count, istep_count)
c     s kelkar june 16, 05
c     set flag in irray(i,0) for cliff nodes

      use comai, only : ierr, iptty 
      use comsk
      use comsptr

      implicit none

      integer j,cliff_count,cliff_pointer_zero
      integer i1,i2,ij0,ij,ncjj,jj,jj_global,istep_count

      allocate(istep_cases(cliff_count*4))

      cliff_pointer_zero = 100000000

c     -200000000 < irray(i,0) < -100000000 : non-OMR cliff node 
c     irray(i,0) < -200000000 : OMR cliff node
c     if the node is a cliff node, but has a specified boundary
c     outflow at it, remove cliff tag and mark as a regular
c     boundry node in load_omr_flux_array

      i1=cliff_count+1

      do j=1,cliff_count
         jj=node_save(j)
         jj_global=node_global(jj)
c if the node is already defined to be a spring node, dont
c reset the flag, otherwise reset irray(,0)
         if(irray(jj_global,0).ne.-(jj_global+2000000)) then
            irray(jj_global,0)=-(j+cliff_pointer_zero)
         endif
c     (later, in subroutine struct_conn_array (ptrac1.f),
c     for OMR nodes, this is changed to 
c     irray(i,0)=-2e+8-icounter)

         istep_cases(j)=i1
         ncjj=ncase(jj)
         do ij=1,ncjj
            ij0=i1+ij-1
            if(ij0.gt.4*cliff_count) then
               write(iptty,*)'increase dim of istep_cases'
               write(iptty,*)'in cliff_setflag. STOP'
               write(ierr,*)'increase dim of istep_cases'
               write(ierr,*)'in cliff_setflag. STOP'
               stop
            endif
            istep_cases(ij0)=list_case(jj,ij)
         enddo
         i1=ij0+1
      enddo
      istep_count = i1 - 1

      return

      end subroutine cliff_setflag

c................................................................

      subroutine write_cliff(icliff3, cliff_count, istep_count)

      use comsptr, only : irray
      use comsk

      implicit none
      integer cliff_count, icliff3, istep_count, j, jj

      write (icliff3, *) cliff_count, istep_count
      do j = 1, cliff_count
         jj = node_global(node_save(j))
         write (icliff3, *) jj
      end do
      
      do j = 1, istep_count
         write (icliff3, *) istep_cases(j)
      end do

      close (icliff3)

      end subroutine write_cliff
         
c................................................................

      subroutine read_cliff(icliff3)

      use comsptr, only : irray
      use comsk, only : istep_cases

      implicit none
      integer cliff_count, icliff3, istep_count, j, jj
      integer cliff_pointer_zero

      cliff_pointer_zero = 100000000

      read (icliff3, *) cliff_count, istep_count
      allocate (istep_cases(istep_count))

      do j = 1, cliff_count
         read (icliff3, *) jj
         irray(jj,0)=-(j+cliff_pointer_zero)
      end do
      
      do j = 1, istep_count
         read (icliff3, *) istep_cases(j)
      end do

      close (icliff3)

      end subroutine read_cliff
 
