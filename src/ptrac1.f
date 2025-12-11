      subroutine ptrac1
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
!D1 Perform initial setup functions for the streamline particle 
!D1 tracking calculations. 
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/ptrac1.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:14   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
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

      use comai
      use combi
      use comci
      use comdi
      use comsptr
      use compart
      use comsk
      use comxi, only : nmfil
      use davidi
      use comwt

      implicit none
      integer neqp1,n50
      integer i3,ii1,ii2,i1,i33,ix,iy,iz
      integer ip,is,i5,kbp,kbm,kb,i,np1
      integer current_node
      integer  n_porosi0,itemp_col,itemp_node, kb2,flag_box,inp1
      real*8 dx,dy,dz,ep,ep5,x60,x33
      real*8 rprime, spacing
      real*8 xcoordw, ycoordw, zcoordw, del_plus, del_minus
      real*8 ps_print
      real*8 s_print
      integer connect_flag, upper_limit
      integer iprcount, iprint
      real*8 tol_c
      parameter(tol_c=1.d-20)
      integer position_in_string, final_position
      integer n_written, iprops, jprops, ijkv_find
      
c......dec 4 01 s kelkar insert omr changes ...................

      integer inode,iwsk

      real*8 aread,aread_max
      integer flag_sk

      integer iomr_flag
      integer idbg,kdbg
      real*8 gotcord

      integer npart_ist2

      real*8 epsilonwt
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Called from insptr now
!      if(cliff_flag) call setup_cliffnodes

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      epsilonwt=1.e-12
      count_steps = 0
      if (omr_flag) iboulist = 0
      if (omr_flag) ipc_save = 0

! Added to read static arrays for calculations once they have been computed zvd 20-Apr-04
      if (save_omr) then
         inquire ( file = nmfil(23), exist = sptr_exists )
         if (sptr_exists) then
! Read arrays, no need to recompute
            call sptr_save (0)
         end if
      else
         sptr_exists = .false.
      end if

c!!!for debuggin 5/13/04
c         do i=1,neq
c           if((i).eq.1..or.ddy(i).eq.1..or.ddz(i).eq.1.) then
c              write(*,*)'i,ddx(i),ddy(i),ddz(i);',i,ddx(i),ddy(i),ddz(i)
c           endif
c         enddo
c!!!!!!

c     Initialize parameters

      neqp1=neq+1
      ep=1.e-7
      ep5=.5*ep

c     Initialize particle time arrays

      ttt=0
      if (irsttime .ne. 0) then
         tt1 = rsttime
c     ZVD - 10-Dec-09, time can't be less than the simulation start time
         do i = 1, num_part
            if (ttp1(i) .lt. tt1) ttp1(i) = tt1
         end do
      else
         tt1= 0.d0
      end if
      tt = 0
! ZVD - 14-Oct-05, initialized in insptr (for particle restarts)
!      ttp1=0

c     Initialize omr arrays if data not read in
      if (.not. sptr_exists .and. omr_flag) then
         node_omr = 0
         isave_omr = 0
      end if

c 5-Nov-02 zvd Add default nsp and icns values (liquid only)
      nsp = 1
      icns = 1


      if(.not.unstruct_sptr) then

c     Determine setup arrays for structured grids

c     Call routine to construct the connectivity array for
c     structured grids
c !!!!!!!!!!!!!!!!!5/13/04 next two lines taken out for debugging
         if (.not. sptr_exists) then

            call struct_conn_array(iomr_flag)

c.. s kelkar sep 28 04, moved call to geom_array inside conn_array
c     Call routine to determine geometric sizes of structured grid cells
c            call struct_geom_array

         endif

         if(node_count.eq.0) omr_flag=.false.

      else

c     Setup for unstructured grids

         call unstruct_arrays

      end if

c     Determine initial state of particles (positions, nodes)

      if(ist.eq.2) then
         call particle_patch(npart_ist2)
      end if

c  *******initial state *********
c     Find cell location where particle starts
      if(abs(ist).ge.1) then

         call find_particle_cell(npart_ist2,ierr,iptty)

c..............................................................
c...Oct 15 2008 s kelkar if porosity<=0 move down the column
c     routine wtsi_column sorts nodes in vertical columns
c     wcol(node)=column# corresponding to the node
c     n_col(column#)=# of nodes in the column
c     col(kb,wcol(node))=node #s in that column
         n_porosi0= 0 
         do i=1,num_part
            if (ijkv(i) .ne. 0) then
               if(ps(ijkv(i)).le.0.) then
                  if (ifree .ne. 0) then
c     Only try to move particle if wtsi problem
                     n_porosi0 =  n_porosi0 +1
                     if( n_porosi0.eq.1) then
                        if(.not.allocated(wcol)) then
                           call wtsi_column
                        endif
                     endif
                     inp1=ijkv(i)
                     itemp_col = wcol(inp1)
                     do kb2 = 1, n_col(itemp_col)
                        itemp_node=col(itemp_col,kb2)
                        if(inp1.eq.itemp_node) then
                           do kb = kb2+1, n_col(itemp_col)
                              itemp_node=col(itemp_col,kb)
                              if(ps(itemp_node).gt.0.) then
                                 ijkv(i) = itemp_node
                                 z1(i)=cord(itemp_node,3)-
     &                                corn(itemp_node,3)
                                 goto 96969
                              endif
                           enddo
c     did not find a porosity>0 node in the column. Do a neighbor search
                           call tree_search_porosity(ijkv(i),5,flag_box)
                           if(flag_box.gt.0) then
                              ijkv(i)=flag_box
                              z1(i) = cord(flag_box,3)
                           else
                              write(ierr, 222) i, ijkv(i) 
                              ijkv(i) = 0
                              istop(i) = 1                           
c                           write(ierr,*)"error in ptrac1. can't find"
c                           write(ierr,*)'neighbor with porosity>0. ',
c     &                          'STOP.'
c                           stop
                           endif
                        endif
                     enddo
96969                continue
                  else
c     Just remove the particle 
                     write(ierr, 223) ijkv(i), i
                     istop(i) = 1
                     ijkv(i) = 0
                  end if              
               endif
            end if
         enddo
 222     format ("Error in ptrac1: can't find neighbor with porosity>0",
     &        ' for particle ', i8, 'at node ', i8)
 223     format ('Error in ptrac1: Invalid particle start ',
     &        'at 0 porosity node ', i8, ' for particle number ', i8)
c..............................................................

         if (ist .eq. 2) then
            do i = 1, num_part
               part_id(i,1) = i
               part_id(i,2) = ijkv(i)
            end do
         end if
      end if
c zvd 06-21-07 Set x3,y3,z3 to initial particle location
      do i = 1, num_part
         if (abs(ijkv(i)) .ne. 0) then
            x3(i) = x1(i) + corn(ijkv(i), 1)
            y3(i) = y1(i) + corn(ijkv(i), 2)
            z3(i) = z1(i) + corn(ijkv(i), 3)
c zvd 03-16-2010 Set x3, y3, z3 to initial particle location in insptr for ist = 1
         else if (ist .ne. 1) then
            x3(i) = x1(i)
            y3(i) = y1(i)
            z3(i) = z1(i)
         end if
      end do

c...for a quick fix of bouncing particles s kelkar 1/16/02
      do np1=1,num_part
         oldnode(np1)=ijkv(np1)
         oldnode2(np1)=ijkv(np1)
         oldnode_rand(np1)=ijkv(np1)
      enddo
c......................................

c     Output of particle information changed: BAR 6-15-99
c     
      if(iprto.ne.0) then
         call output_info
      end if

c     Subroutine to initialize particle tracking transport parameters
c s kelkar may 28 09 moved call to fehmn.f where ptrac1 used to be called

c      call init_sptr_params
      
c s kelkar may 20 09 moved call to fehmn.f where ptrac1 used to be called
c      if (.not. compute_flow) then
c         if (.not. sptr_exists) call load_omr_flux_array
c         if(.not.random_flag) then
c            if(allocated(sx)) deallocate(sx)
c            if(allocated(istrw)) deallocate(istrw)
c         end if
c      end if

! Added to save static arrays for calculations once they have been computed zvd 20-Apr-04
c zvd 19-Nov-2010
! Moved to fehmn, needs to be called after call to load_omr_flux_array
!      if (.not. sptr_exists .and. save_omr) call sptr_save (1)

      return

      contains

******************************************************************
******************************************************************

      subroutine struct_conn_array(iomr_flag)

      use comsk
      use comai, only : ierr, iout, iptty
      implicit none
      integer ibou, iflag_boundry, iomr_flag,idum,irray0,idumm


c **** get irray *****

      if(icnl.eq.0) then
         upper_limit = 3
      else
         upper_limit = 2
      end if

      do i=1,neq
        do i3=-3,3
           ggg(i,i3)=0.
        enddo
      enddo

c ...dec 4 01, s kelkar  nov 1 01, OMR stuff ...............................
      node_count=0
c.................................................................

      do i=1,neq
        do i3=-3,3
c     do not zero out irray(i,0) because it has particle capture 
c     information read in insptr
           if(i3.ne.0) irray(i,i3)=0
        enddo
        
c     ...dec 4 01, s kelkar  nov 1 01, OMR stuff ......................
        iomr_count=0
        iomr_flag=0
        
        isave_omr=0

        do i1=1,iomrmax
           iomr_neighbour(i1)=0
        enddo
c.................................................................

        ii1=nelm(i)+1
        ii2=nelm(i+1)

c.....dec 4, 01 s kelkar  9/21/01 find the max area connection at this node
c     for flagging connection-areas that are zero (less than epsilon)
        aread_max=0.
        aread=0.
        do i1=ii1,ii2
           kb=nelm(i1)
           if(kb.ne.i) then
c  istrw is a pointer array for sx, corresponding to the connection
c  i-i1. Generally (but not always) it is a scalar, the negative of the 
c  magnitude of (cross sectional area0/3. divided by the inter-nodal 
c  distance
              iwsk=istrw(i1-neq-1)
              aread=-sx(iwsk,1)
              if(aread.gt.aread_max) aread_max=aread
           endif
        enddo
        aread_max=aread_max*1.e-8
c...................................................

        do i1=ii1,ii2
          kb=nelm(i1)
          if(kb.ne.i) then
c     Do loop filters out any connections that aren't 
c     oriented only in the x, y, or z directions

             connect_flag = 0
             do i3=1,upper_limit
                x33=cord(kb,i3)-cord(i,i3)
                if(x33.gt. tol_c ) connect_flag = connect_flag + 1
                if(x33.lt.-tol_c ) connect_flag = connect_flag + 1
             end do
             
             if(connect_flag.lt.2) then
c     if connect_flag =1 then its a good connection, store the node
c     number in the correct slot of irray
                i33=0
                do i3=1,upper_limit
                   x33=cord(kb,i3)-cord(i,i3)
                   if(x33.gt. tol_c) then
                      i33= i3
                   endif
                   if(x33.lt.-tol_c) then
                      i33=-i3
                   endif
                enddo
c s kelkar, may 25,04, check for -ve porosity at kb, these are treated 
c as no flow connections. Flag these with -ve sign
! Use rock matirx pososity here to account for nodes that have been
! eliminated using negative porosities
                if(ps(kb).gt.0.) then
                   irray(i,i33)=+kb
                else
                   irray(i,i33)=-kb
                endif
             endif
             
             if(connect_flag.gt.1) then

c     ........... dec 4 01 s kelkar   9/21/01..........................
c     if connect_flag > 1 then the connection is not lined up
c     with any of the axis. If this connection has a nonzero area
c     then it signals change from a structured to unstructured 
c     part of the grid. Howver, zero area connections that are
c     not lined up with axis can also occure
c     in a structured part of the grid, ie a diagonal in a square.
c     these have to be filtered out.
             
                iwsk=istrw(i1-neq-1)
                aread=-sx(iwsk,1)
                if(aread.gt.aread_max) then
                
c     ...s kelkar  nov 1 01, OMR stuff ...............................
c     here is how the information is stored:
c     node_count counts the number of nodes that have at least one bad 
c     connection with any one of its neighbours, and the node number of 
c     is stored in node_omr(node_count).
c     for each such node i, the # of its neighbours that have a bad 
c     connection with it are counted in iomr_count, and the pointers
c     (i1) to the  neighbour 
c     node numbers are stored temporarily in iomr_neighbour(iomr_count)
c     When the loop ovver the neighbours 'kb' is finished, then for each
c     node i with a bad connection, the subroutine 'subomr' is called,
c     which sorts out the faces on which bad connections occur.
c     komr_count(k) is a temporary counter for the number of neighbours
c     of 'i' on a particular face (k). and the pointers (i1) to the node
c     numbers of such neighbours are saved in 
c     isave_omr(face#,komr_count(k)) 
                
                   if(iomr_flag.eq.0) then
                      iomr_flag=1
                      node_count=node_count+1
                      if(node_count.gt.omr_nodes) then
                         write (ierr, 1001) node_count, omr_nodes
                         if (iptty .ne. 0) 
     .                        write (iptty, 1001) node_count, omr_nodes
                         call exit_ptrac1
                      endif
 1001                 format ( 'ERROR in PTRAC1: count ', i8, 
     &                     ' greater than number of omr nodes ',
     &                     i8, /, 'STOPPING')
                      node_omr(node_count)=i
                   endif
                   iomr_count=iomr_count+1
                   if(iomr_count.gt.iomrmax) then
                      write(ierr,*)'iomr_count.gt.iomrmax in ptrac1'
                      write(ierr,*)'increase dimension of ',
     .                     'iomr_neighbour'
                      write(ierr,*)'STOPPING'
                      if (iptty .ne. 0) then
                         write(iptty,*)'iomr_count.gt.iomrmax in ptrac1'
                         write(iptty,*)'increase dimension of ',
     .                        'iomr_neighbour'
                         write(iptty,*)'STOPPING'
                      end if                         
                      call exit_ptrac1
                   endif
c     NOTE: saving the pointer i1 rather than the node # kb because
c     i1 can be directly used as a pointer to sx and a_axy arrays
c     and kb can be retrieved from nelm(i1)
                   iomr_neighbour(iomr_count)=i1
                   
c     close(87)
                   
                endif
             endif
c.................................................................
          
          endif

        enddo

        if(iomr_flag.eq.1) then
c OMR node can not be specified as a well-capture node. At this stage
c this could have only come from insptr- sptr input file must be
c modified
           irray0=irray(i,0)
           if(irray0.eq.-i.or.irray0.eq.-(i+1000000)
     1          .or.(irray0.lt.-(10000000).and.irray0.gt.-(100000000))
     2          )then
              write(ierr,*)'OMR node can not specified as a '
              write(ierr,*)'well-capture node.sptr input file'
              write(ierr,*)'must be modified. check keyword cpatur'
              write(ierr,*)'SUBROUTINE struct_conn_array (ptrac1)'
              write(ierr,*)'Node Number=',i
              write(ierr,*)'STOP'
              call exit_ptrac1
           endif

           if(irray0.lt.-100000000) then
              idumm=-(irray0+100000000)
              irray0=-200000000-idumm
           elseif(irray(i,0).ne.-(i+2000000)) then
              irray0 = 0
           endif
           irray(i,0) = irray0

c omr node allowed to be a spring node

c....s kelkar  Jan 27march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular, interior node, not a source/sink
c     irray(i,0) = -i-2000 : regular node on a external boundary
c     irray(i,0) = -i : regular interior node that is a sink/source
c                        but not explicitly specified in sptr macro
c     -100000000 < irray(i,0) < -10000000 : regular interior node that is 
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c     -200000000 < irray(i,0) < -100000000 : non-OMR cliff node 
c     irray(i,0) < -200000000 : OMR cliff node
c if the node is a cliff node, but has a specified boundary
c outflow at it, remove cliff tag and mark as a regular
c boundry node in load_omr_flux_array
c                = -(i+2000000) : spring node
c                = -(i+1000000) : well-capture node on external bound
c                   similar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c
c flag OMR nodes that are on the exterior boundary, and save the 
c exterior faces in iboulist. Also set irray(i,0)=-i-1000 for
c OMR nodes on a external boundary

           call boundary_nodes(i,+1)
c...................................................

        else
c flag non-OMR nodes that are on the exterior boundary, and save the 
c exterior faces in iboulist.
           call boundary_nodes(i,0)
           irray0=irray(i,0)
           if(irray0.ne.-i.and.irray0.ne.-(i+1000000).and.
     1          irray0.ne.-(i+2000000).and.irray0.ne.-(i+2000).
     2          and.irray0.gt.-10000000) then
              irray(i,0) = i
           endif
        endif

c ...dec 4 01 s kelkar  nov 1 01, OMR stuff ..........................
        if(iomr_flag.eq.1) call subomr2(i)
c.................................................................
c... s kelkar sep 28, 04  reduce storage by calling geom_array here
c and changin the dim of iomr_save to (la,k)

c calculate del_plus, del_minus- the distances to the 
c control volume boundaries, for use in ddx_corn_array
        call struct_geom_array(i)

      enddo

c set up the ddx,ddy,ddz and corn arrays
      call ddx_corn_array
c      
      do i=1,neq
        do i3=-3,3
           ggg(i,i3)=0.
        enddo
      enddo

      write(iout, 1003) omr_nodes, node_count
      if (iptty .ne. 0) write(iptty, 1003) omr_nodes, node_count
 1003 format ('Number of OMR nodes set to ', i8, /,
     .     'Actual number of OMR nodes ', i8)

c **** got irray *****

      return

      end subroutine struct_conn_array

******************************************************************
******************************************************************

      subroutine output_info

      implicit none
      real*8 sptr_time
      character*200 sptr_heading, sptr_prop_values

c     Sets up output of particle tracking info

c     Only used for regular output with transient particle start times, 
c     trans_flag is used for minimal output options
      pstart_out = .true.

c     Robinson added minimal write option 3-11-02
      if(iprto.lt.0) then
         if (xyz_flag) then
c ZVD modified minimal write option to include coordinates 12-08-2005
            if (iprto .eq. -1) then
               write(isptr2,100)
            else
               write (isptr2) 'XYZ'
            end if
            do np1 = 1, num_part
               current_node = ijkv(np1)
               xcoordw = x1(np1) + corn(current_node,1)
               ycoordw = y1(np1) + corn(current_node,2)
               zcoordw = z1(np1) + corn(current_node,3)
               sptr_time = ttp1(np1)
               if (iprto .eq. -1) then
                  write(isptr2,105) part_id(np1,1),sptr_time,
     &                 current_node,xcoordw, ycoordw, zcoordw
               else
                  write (isptr2) part_id(np1,1),sptr_time,current_node,
     &                 xcoordw, ycoordw, zcoordw 
               end if
            end do
         else if (ip_flag .or. trans_flag) then
c ZVD added option to write initial position to abbreviated output file
            if (iprto .eq. -1) then
               if (ip_flag) then
                  write(isptr2, 110) ''
               else
                  write(isptr2, 110) 'TRA : '
               end if
            else
               if (trans_flag) write (isptr2) 'TRA'
            end if
            do np1 = 1, num_part
c ZVD added option for transient where particle start time is saved but
c     initial node is set to 0
c ZVD 07-Feb-2011 negative of starting node is now output 
c (this way particles that have been excluded can be distinguished from
c particles that have a delayed start time) 
               if (trans_flag) then
c                  current_node = 0
                  current_node = -ijkv(np1)
               else
                  current_node = ijkv(np1)
               end if
               sptr_time = ttp1(np1)
               if (iprto .eq. -1) then
                  write(isptr2,105) part_id(np1,1),sptr_time,
     &                 current_node
               else
                  write (isptr2) part_id(np1,1),sptr_time,current_node
               end if
            end do
         end if            
 100     format ('XYZ : Part_no    time_days     cell_leaving',
     &        '       X             Y             Z')
 105     format(1x,i8,1x,g21.14,1x,i8,3(1x,g16.9))
 110     format (a, 'Part_no    time_days     cell_leaving')
      elseif(iprto.eq.1) then
         pstart_out = .false.
         do iprint = 1, 200
            sptr_heading(iprint:iprint) = ' '
         end do
         sptr_prop_values = ''

         sptr_heading(1:58) =
     2      ' particle_number  x(m)      y(m)      z(m)      time(days)'
c     2      ' particle_number  x         y         z         time      '
         position_in_string = 60
c     Determine how many property columns are being written

         n_written = 0
         do iprops = 1, nsptrprops
            if(write_prop(iprops).ne.0) then
               n_written = n_written + 1
            end if
         end do

c     Find which one is written next, write it to string
         do iprops = 1, n_written

            inner: do jprops = 1, nsptrprops
            if(write_prop(jprops).eq.iprops) then

               if(jprops.eq.1) then
                  sptr_heading(position_in_string:position_in_string+9)
     2                 = 'porosity  '
                  position_in_string = position_in_string + 10
                  exit inner
c     Fluid saturation
               elseif(jprops.eq.2) then
                  sptr_heading(position_in_string:position_in_string+11)
     2                 = 'saturation  '
                  position_in_string = position_in_string + 12
                  exit inner
c     Permeability
               elseif(jprops.eq.3) then
                  sptr_heading(position_in_string:position_in_string+13)
     2                 = 'permeability  '
                  position_in_string = position_in_string + 14
                  exit inner
c     Rock density
               elseif(jprops.eq.4) then
                  sptr_heading(position_in_string:position_in_string+13)
     2                 = 'rock_density  '
                  position_in_string = position_in_string + 14
                  exit inner
c     Pressure
               elseif(jprops.eq.5) then
                  sptr_heading(position_in_string:position_in_string+9)
     2                 = 'pressure  '
                  position_in_string = position_in_string + 10
                  exit inner
c     Temperature
               elseif(jprops.eq.6) then
                  sptr_heading(position_in_string:position_in_string+12)
     2                 = 'temperature  '
                  position_in_string = position_in_string + 13
                  exit inner
c     Zone number
               elseif(jprops.eq.7) then
                  sptr_heading(position_in_string:position_in_string+5)
     2                 = 'zone  '
                  position_in_string = position_in_string + 6
                  exit inner
c     Particle ID
               elseif(jprops.eq.8) then
                  sptr_heading(position_in_string:position_in_string+3)
     2                 = 'ID  '
                  position_in_string = position_in_string + 4
                  exit inner
               end if
               
            end if
            end do inner

         end do

         sptr_heading(position_in_string:position_in_string+17)
     2        = 'old_node  new_node'
         final_position = position_in_string+17
         
         write(isptr2,'(a)') trim(sptr_heading)
         do np1 = 1, num_part
            current_node = ijkv(np1)
            xcoordw = x1(np1) + corn(current_node,1)
            ycoordw = y1(np1) + corn(current_node,2)
            zcoordw = z1(np1) + corn(current_node,3)
            sptr_prop = 0.
            if(current_node.eq.0) then
               iprcount = 0
               do iprint = 1, nsptrprops
                  if(write_prop(iprint).ne.0) then
                     iprcount = iprcount + 1
                     sptr_prop(write_prop(iprint)) = 1.d-30
                  end if
               end do
c     ps_print = 1.e-30
c     s_print = 1.e-30
            else
               iprcount = 0
               
c     Porosity
               if(write_prop(1).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(1)) = ps_trac(current_node)
               end if
c     Fluid saturation
               if(write_prop(2).ne.0) then
                  iprcount = iprcount + 1
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     sptr_prop(write_prop(2)) = s(current_node)
                  else
                     sptr_prop(write_prop(2)) = 1.0d0
                  end if
               end if
c     Permeability
               if(write_prop(3).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(3)) = 1.d-6*pnx(current_node)
               end if
c     Rock density
               if(write_prop(4).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(4)) = denr(current_node)
               end if
c     Pressure
               if(write_prop(5).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(5)) = phi(current_node)
               end if
c     Temperature
               if(write_prop(6).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(6)) = t(current_node)
               end if
c     Zone number
               if(write_prop(7).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(7)) = izonef(current_node)
               end if
c     Particle ID
               if(write_prop(8).ne.0) then
                  iprcount = iprcount + 1
                  sptr_prop(write_prop(8)) = part_id(np1, 2)
               end if

               position_in_string = 1

               do i = 1, iprcount
                  if (i .eq. write_prop(7) .or. i .eq. write_prop(8)) 
     &                 then
                     write (sptr_prop_values(position_in_string:
     &                    position_in_string+9), '(i8,2x)') 
     &                    int(sptr_prop(i))
                     position_in_string = position_in_string+10
                  else
                     write (sptr_prop_values(position_in_string:
     &                    position_in_string+17), '(g16.9,2x)') 
     &                    sptr_prop(i)
                     position_in_string = position_in_string+18
                  end if
               end do
               write (sptr_prop_values(position_in_string:
     &              position_in_string+17), '(2(i8,2x))')
     &              current_node, current_node

               sptr_time = ttp1(np1)
c Don't output if particle time is greater than starting time
c               if (sptr_time .gt. days) sptr_time = days
               if (sptr_time .le. days) then
                  pstart_out(np1) = .true.
                  position_in_string = len_trim(sptr_prop_values)

                  write(isptr2,8001) part_id(np1,1), xcoordw, ycoordw, 
     2                 zcoordw, sptr_time, 
     3                 sptr_prop_values(1: position_in_string)
               end if
            end if

         end do
      end if
 8001 format(1x, i8, 3(1x,g16.9), 1x, g21.14, 200a)
      return
      
      end subroutine output_info

******************************************************************
******************************************************************


      subroutine struct_geom_array(i)
      
      use comsptr
      use comsk
      
      implicit none

      integer i,j,kb,kb_omr,jab

      real*8 gotcord, delkb,del

c NOTE: using ggg as scratch storage for saving del+ and del-
c temporarily

      do j=-3,3
         delkb=0.
         if(j.ne.0) then
            jab=abs(j)
            kb= abs (irray(i, j))
            if((irray(i,0).eq.0).or.(irray(i,0).eq.-(i+1000)).or.
     1          (irray(i,0).lt.-200000000)) then
               kb_omr=0
               call getcord(i,kb_omr,j,gotcord)
               del = 0.5*abs((gotcord-cord(i,jab)))
c     on return from gotcord, kb_omr =0 only for boundary nodes
c     use expected symmetry of ggg to set ggg(kb_omr,-j)
c     Also if del(i,j) and del(kb,-j)
c     are not equal, then use the greater of the two- this way we may 
c     increase overlap, but we reduce chances of holes.
c     for omr nodes, kb can be 0 either bcs its a boundry node or bcs
c     grid refinement has a missing node there. In the case of a
c     missing node, return value of kb_omr can be a legitimate node.
c      in that case, use symmetry of ggg and irray. ggg(kb_omr could be
c     nonzero thru this process even if kb_omr lt i. Hence when i 
c     becomes value of kb_omr, that del is recalculated and the  
c     greater value is used.

                  ggg(i,j)=del
                  if(irray(i,j).eq.0 .and. kb_omr .gt. 0) then
                     if (ps(kb_omr).gt.0.) then
                        irray(i,j)=+kb_omr
                     else
                        irray(i,j)=-kb_omr
                     endif
                  endif

c               if(kb_omr.gt.0) then 
c                  delkb=ggg(kb_omr,-j)
c                  if(delkb.gt.del) del=delkb
c                  ggg(kb_omr,-j)=del
c                  ggg(i,j)=del
c                  irray(i,j)=kb_omr
c                  irray(kb_omr,-j)=i
c               else
c                  ggg(i,j)=del
c               endif

               else

c            elseif(ggg(i,j).eq.0.) then
c     non-OMR node.
c     there is no need to recalculate ggg if it is already nonzero
               kb_omr=abs(kb)
               call getcord(i,kb_omr,j,gotcord)
               del = 0.5*abs((gotcord-cord(i,jab)))
               ggg(i,j)=del
c               if (kb_omr.gt.0) ggg(kb_omr,-j)=del
            endif
         endif

      enddo

      return

      end subroutine struct_geom_array

c******************************************************************
c     Subroutine is empty - not yet implemented

      subroutine unstruct_arrays

      return
      end subroutine unstruct_arrays

******************************************************************
******************************************************************

      subroutine particle_patch(npart_ist2)
      implicit none

      integer npart_ist2

c  *******initial state section *********

      if(nx.eq.1) then
         dx=0
      else
         dx=xdim/(1.e-20+nx-1)
      endif
      if(ny.eq.1) then
         dy=0
      else
         dy=ydim/(1.e-20+ny-1)
      endif
      if(icnl.ne.0) nz = 1
      if(nz.eq.1) then
         dz=0
      else
         dz=zdim/(1.e-20+nz-1)
      endif
      do ix=1,nx
         do iy=1,ny
            do iz=1,nz
               ip=ix+(iy-1)*nx+(iz-1)*nx*ny
               x1(ip)=x10+(ix-1)*dx
               y1(ip)=y10+(iy-1)*dy
               z1(ip)=z10+(iz-1)*dz
            enddo
         enddo
      enddo

      npart_ist2=ip

c....dec 4 01, s kelkar  sep 20 2001
      if(ist.eq.3) then
         do ip=1,num_part
            inode=ijkv(ip)
            x1(ip)=cord(inode,1)
            y1(ip)=cord(inode,2)
            z1(ip)=cord(inode,3)
         enddo
      endif
c.....................................

      return
      end subroutine particle_patch

******************************************************************

      subroutine find_particle_cell(npart_ist2,ierr,iptty)
      implicit none

      integer flag_box,npart_ist2,ijkv_last,ierr,iptty
      integer nout,save_out_face(6)

c *** fudged initial
      xo=x1
      yo=y1
      zo=z1
c     If the search algorithm starts at node 1 for all nodes,
c     it can't get inside a locally structured part of the grid
c     from outside of it. The following search gets each particle
c     close or even at its starting node, so as long as that
c     coordinate is inside a locally structured part of the
c     grid and inside the model domain, the subsequent search
c     should work.  BAR 3-9-00
c     
      
      if(abs(ist).eq.1) then
! If we have read in the corresponding nodes use them otherwise search
         if (.not. sptr_snode) then
c     Particles could be anywhere, do search on each one
            do i1 = 1, num_part
               call near3(xo(i1),yo(i1),zo(i1),ijkv(i1),0)
            end do
         end if
      elseif(ist.eq.2) then
c     Particles are in a cluster, do search on first
c     one only, use that as starting location
         call near3(xo(1),yo(1),zo(1),ijkv_find,0)
         ijkv(1)=ijkv_find
      
c..s kelkar Jan 21,05 replacing with tree-search      
c     *** n50 must exceed max flow path length/smallest nodal separation
c      n50=100000
c      
c      
c      do is=1,n50 
c         ijkvss=ijkv
c         x1=xo-corn(ijkv,1)
c         y1=yo-corn(ijkv,2)
c         z1=zo-corn(ijkv,3)
c         ddxv=ddx(ijkv)
c         ddyv=ddy(ijkv)
c         ddzv=ddz(ijkv)
c         
c         ijkvs=ijkv
c         where(x1/ddxv.gt.1.)
c            ijkv=irray(ijkv, 1)
c         endwhere
c         where(x1/ddxv.lt.0.)
c            ijkv=irray(ijkv,-1)
c         endwhere
c         where(ijkv.eq.0) ijkv=ijkvs
c         
c         ijkvs=ijkv
c         where(y1/ddyv.gt.1.)
c            ijkv=irray(ijkv, 2)
c         endwhere
c         where(y1/ddyv.lt.0.)
c            ijkv=irray(ijkv,-2)
c         endwhere
c         where(ijkv.eq.0) ijkv=ijkvs
c         
c         ijkvs=ijkv
c         where(z1/ddzv.gt.1.)
c            ijkv=irray(ijkv, 3)
c         endwhere
c         where(z1/ddzv.lt.0.)
c            ijkv=irray(ijkv,-3)
c         endwhere
c         where(ijkv.eq.0) ijkv=ijkvs
c         
c         i5=0
c         do i1=1,num_part
c            if(ijkv(i1).ne.ijkvss(i1)) i5=1
c         enddo
c         if(i5.eq.0) go to 201
c         
c      enddo
c 201  continue
c........................................................
         ijkv_last=ijkv_find
         do is=2,npart_ist2     
            call tree_search(ijkv_last,20,ierr,iptty,
     $           xo(is),yo(is),zo(is),flag_box,nout,save_out_face)
            if(flag_box.gt.0) then
               ijkv(is)=flag_box
               ijkv_last=flag_box
            else
               call tree_search(ijkv(1),20,ierr,iptty,
     $              xo(is),yo(is),zo(is),flag_box,nout,save_out_face)
               if(flag_box.gt.0) then
                  ijkv(is)=flag_box
                  ijkv_last=flag_box
               else
c     tree-search failed, do a global search
                  call near3(xo(is),yo(is),zo(is),flag_box,0)
                  if(flag_box.gt.0) then
                     ijkv(is)=flag_box
                     ijkv_last=flag_box
                  else
                     write(ierr,*)'Error in find_particle_cell. Ist=2 '
                     write(ierr,*)'Cant find the CC for particle #', is
                     write(ierr,*)'xo,yo,zo=',xo(is),yo(is),zo(is)
                     write(ierr,*)'STOP'
                     if (iptty. ne. 0) then
                    write(iptty,*)'Error in find_particle_cell. Ist=2 '
                    write(iptty,*)'Cant find the CC for particle #', is
                        write(iptty,*)'xo,yo,zo=',xo(is),yo(is),zo(is)
                        write(iptty,*)'STOP'
                     end if
                     call exit_ptrac1
                  end if
               endif
            end if
         enddo
      end if

         ijkvs=ijkv

c***  
c     s kelkar aug 29, 06
c     for water table nodes, if S(ijkv(np1))>Smin then set the initial
c     position below ddz*S, if S<Smin then move the particle vertically
c     downward until a node with S>Smin is encountered.

c     zvd added to time loop (ptrac3), to move particle down to wt 
c     when it starts to move aug 27, 2007
c      if (ifree.ne.0) then
c         if(deltawt.gt.epsilonwt) then
c            call wtsi_ptrac1_init
c         endif
c      endif
c****


      x1=xo-corn(ijkv,1)
      y1=yo-corn(ijkv,2)
      z1=zo-corn(ijkv,3)
      
c     ********at final relative initial state *********
      
      ddxv=ddx(ijkv)
      ddyv=ddy(ijkv)
      ddzv=ddz(ijkv)
      
c     **** is the initial state valid? *****
cc      where((ijkv.lt.1).or.(ijkv.gt.neq)) ijkv=0
c      where((x1.gt.ddxv).or.(x1.lt.0.))     ijkv=0
c      where((y1.gt.ddyv).or.(y1.lt.0.))     ijkv=0
c      if(icnl.eq.0) then
c         where((z1.gt.ddzv).or.(z1.lt.0.))     ijkv=0
c      end if

      do is=1,num_part
c     Check to see if particle should be excluded if out side the model domain
         if (exclude_particle) then
            if(x1(is).gt.ddxv(is) .or. x1(is).lt.0. .or. 
     &           y1(is).gt.ddyv(is) .or. y1(is).lt.0. .or.
     &           ps(ijkv(is)) .le. 0.) then
               istop(is) = 1
               ijkv(is) = 0
            end if
            if (icnl.eq.0) then
               if(z1(is).gt.ddzv(is) .or. z1(is).lt.0.) then
                  istop(is) = 1
                  ijkv(is) = 0
               end if
            end if           
         else
            if(x1(is).gt.ddxv(is)) x1(is)=ddxv(is)
            if(x1(is).lt.0.) x1(is)=0.
            if(y1(is).gt.ddyv(is)) y1(is)=ddyv(is)
            if(y1(is).lt.0.) y1(is)=0.
            if(icnl.eq.0) then
               if(z1(is).gt.ddzv(is)) z1(is)=ddzv(is)
               if(z1(is).lt.0.) z1(is)=0.
            end if
         end if
      enddo      

c....dec4 01  s kelkar  nov  11 01 commented out next 3 lines to
c    allow omr nodes as  initial locations
c      do np1 = 1, num_part
c         if(irray(ijkv(np1),0).lt.0) ijkv(np1) = 0
c      end do
c...................................
      
c     *** print initial state ******
c     ***** stop if initial state of any particle is invalid: ijkv=0 ****
      istop=0
      do i1=1,num_part
         if(ijkv(i1).eq.0) then
            write(ierr, 224) i1
!            call exit_ptrac1
            istop (i1)=1
         end if
      enddo
 224  format ('Error in ptrac1: Initial state of particle is invalid',
     &     ' for particle number ', i8)

c     **** set istop=1 if point out of domain****
c      where(ijkv.eq.0) istop=1
      
c     *** move initial points off element boundaries***
      ddxv=ddx(ijkv)
      ddyv=ddy(ijkv)
      ddzv=ddz(ijkv)
      where(x1.eq.0. ) x1=ddxv*ep
      where(x1.eq.ddxv) x1=ddxv*(1.-ep)
      where(y1.eq.0. ) y1=ddyv*ep
      where(y1.eq.ddyv) y1=ddyv*(1.-ep)
      if(icnl.eq.0) then
         where(z1.eq.0. ) z1=ddzv*ep
         where(z1.eq.ddzv) z1=ddzv*(1.-ep)
      end if
               
      return
      end subroutine find_particle_cell
               
c******************************************************************

      end subroutine ptrac1

c***********************************************************************

      subroutine subomr2(i)

c s kelkar 11 jul 05
c this is  a modified version (and simplified) of 'subomr'
c  each connected node with a 'bad' connection, from array
c iomr_neighbour(:) is counted and stored for every direction
c that the connection if off-axis in the array isave_omr(:,:)
c connections are not classified as type I or II (that is done
c in subomr, but no longer needed. In the current version of ptrac1,
c that classification leads to wrong ddx,ddy,ddz values)

      use comai
      use combi
      use comci
      use comdi
      use comflow
      use comsptr
      use comsk
      use compart

      implicit none

      integer i,i1,k,ia,kb,ikb,ipos

      real*8 epsilon,xkb,xia,d

      epsilon=1.e-8

      do i1=-3,3
         komr_count(i1)=0
      enddo
      
      do i1=1,iomr_count
         
c     NOTE: in iomr_neighbour,saved ikb rather than the node # kb 
c     because ikb can be directly used as a pointer to sx and a_axy 
c     arrays and kb can be retrieved from nelm(ikb)
         ikb=iomr_neighbour(i1)
         kb=nelm(ikb)
         do k=1,3            
c            ia=abs(irray(i,k))
            d=cord(kb,k)-cord(i,k)
            if(d.gt.epsilon) then
               komr_count(k)=komr_count(k)+1
               if(komr_count(k).gt.komrmax) then
                  write(ierr,*)'komr_count(k).gt.komrmax in ',
     &                 'subomr2'
                  write(ierr,*)'change dimension in comomr.'
                  write(ierr,*)'i,iomr_count',i,iomr_count
                  write(ierr,*)'STOP'
                  if (iptty .ne. 0) then
                     write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                    'subomr2'
                     write(iptty,*)'change dimension in comomr.'
                     write(iptty,*)'i,iomr_count',i,iomr_count
                     write(iptty,*)'STOP'
                  end if
                  call exit_ptrac1
               endif
               isave_omr(k,komr_count(k))=ikb
            endif
            
c            ia=abs(irray(i,-k))
            d=cord(kb,k)-cord(i,k)
            if(d.lt.-(epsilon)) then
               komr_count(-k)=komr_count(-k)+1
               if(komr_count(-k).gt.komrmax) then
                  write(ierr,*)'komr_count(k).gt.komrmax in ',
     &                 'subomr2'
                  write(ierr,*)'change dimension in comomr.'
                  write(ierr,*)'STOP'
                  if (iptty .ne. 0) then
                     write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                    'subomr2'
                     write(iptty,*)'change dimension in comomr.'
                     write(iptty,*)'STOP'
                  end if
                  call exit_ptrac1
               endif
               isave_omr(-k,komr_count(-k))=ikb
            endif            
         enddo
            
      enddo
      
      return

      end
      
c.......................................................................
      
      subroutine getcord(i,kb,la,gotcord)
c     s kelkar, modified Feb 10,05

      use comai, only : ierr, iptty
      use combi
      use comsptr
      use comsk

      implicit none

      integer i,kb,k,l,la,j,jk,jkb,jkbmax,ibou,j1,j2,jkc
      integer itempf, jtemp(200), jtemp_count, lasign,jjjj
      integer i_augment,i1,i2,iii,jkd, jkbb,jkb1,jkb2,i_dir
      integer l_perp,irray0

      real*8 gotcord,dmax,d, djtemp, djkb,dtemp1,dtemp2

      l=abs(la)
      lasign=isign(1,la)

      if(kb.gt.0) then
         gotcord=cord(kb,l)

      elseif(kb.eq.0) then

         irray0=irray(i,0)
         if((irray0.eq.i).or.(irray0.eq.-(i+1000000)).or.
     1        (irray0.eq.-(i+2000000)).or.irray0.eq.-(i+2000).or.
     2        (irray0.lt.-100000000.and.irray0.gt.-200000000)) 
     3        then
c i is a non-OMR boundary node, a well-capture node on boundary
c or a spring node on a boundaryor a cliff node.
c Set gotcord = cord of i
            gotcord=cord(i,l)
            kb=0

c......s kelkar 1/27/04 3-D stuff................
         elseif(irray(i,0).eq.-(i+1000).or.
     1          (irray(i,0).lt.-200000000)) then
c handle OMR nodes on exterior boundaries, including omr-cliff nodes
            itempf = 0
            do ibou=1,6
               if(iboulist(i,ibou).eq.la) then
c     node i is an OMR node on an exterior boundary, with the exterior in 
c     the la direction. set gotcord =cord of i
                  if (irray(i, la) .lt. 0) then
                     kb = irray(i, la)
                     gotcord = cord (abs(kb), l)
                  else
                     gotcord=cord(i,l)
                     kb=0
                   end if
                   itempf = 1
              endif
            enddo

         endif
c.................................................

         if(irray(i,0).eq.0.or.((irray(i,0).eq.(-i-1000).or.
     1          (irray(i,0).lt.-200000000)).and.
     2           itempf.eq.0)) then

c we have either an non-boundary omr node 
c , or a boundary OMR node but with boundary oriented in a plane 
c different from that given by la. Define a 
c fictious control volume face
c find the node furtherest from i in the la direction
c form the list of nodes with a bad connection with i in the
c la direction (these are saved in isave_omr)            
            dmax=0.
            jtemp_count=0
            do j=0,komrmax,1
               if(j.eq.0) then
                  jkb=abs(irray(i,la))
               else
                  jkb=0
                  jk=isave_omr(la,j)
                  if(jk.gt.0) jkb=nelm(jk)
               endif
               if(jkb.gt.0) then
                  d=abs(cord(i,l)-cord(jkb,l))
                  if(d.gt.dmax) then
c saving jkb in jtemp if a search is needed below over 
c the neighbours of jkb to avoind creating holes in the mesh.
                     jtemp_count=jtemp_count+1
                     jtemp(jtemp_count)=jkb
                     dmax=d
                     jkbmax=jkb
                  endif
               endif
            enddo

            if(jkbmax.gt.0 ) then
               gotcord=cord(jkbmax,l)
               kb=jkbmax
            else
               write(ierr,*)'STOP. getcord found jkbmax=0'
               write(ierr,*)'i,kb,la=',i,kb,la
               if (iptty .ne. 0) then
                  write(iptty,*)'STOP. getcord found jkbmax=0'
                  write(iptty,*)'i,kb,la=',i,kb,la
               end if
               call exit_ptrac1
            endif

c now check for the rare situation when
c the node i lies on the side of a rectangle which has the central 
c node missing due to refinement on all sides,
c and if so, getcord is set equal to the node on the other side
c of the squarerectangle, thus creating cc's that overlap in the
c middle of this rectangle, but avoid creating holes in the model
c this situation can only arise if the grid on the 'la' side is
c one level coarse compared to the grid on the '-la' side.
c the search for a node needs to be only over the neighbours of the 
c nodes stored in jtemp(1:jtemp_count)
            jkc=jkbmax
            j1=nelm(jkc)+1
            j2=nelm(jkc+1)
            do j=j1,j2
               jkb=nelm(j)
               if(jkb.ne.i) then
                  djkb=abs(cord(i,l)-cord(jkb,l))
                  if(djkb.eq.dmax*2.) then
c     the node jkb is at the expected distance from i along  la axix
c     check if the other 2 coordinates match
                     do jjjj=1,3
                        if(jjjj.ne.l) then
                          if(cord(i,jjjj).ne.cord(jkb,jjjj)) goto 91911
                        endif
                     enddo
c     the coordinates match, 
c     find the normal axis to the coordinate plane formed the nodes 
c     i,jkc. note 
c     that i and jkb lie along the la coordinate axis
                     do jjjj=1,3
                        if(jjjj.ne.l) then
                           if(cord(i,jjjj).eq.cord(jkc,jjjj)) then
                              i_dir=jjjj
                              goto 91913
                           endif
                        endif
                     enddo
c     i and jkc not in a coordinate plane, continue with other 
c     neighbours of jkc
                     goto 91911
                     
91913                continue
c     now check if the nodes i, jkc and jkb form
c     a closed figure with another neighbour of i in the plane 
c     corrosponding to the coordinate i_dir. need to search only
c     only the nodes on the la side of i, ie isave_omr(la,iii)
                     do iii=1,komrmax
                        jk=isave_omr(la,iii)
                        if(jk.gt.0) then
                           jkd=nelm(jk)
                           if(jkd.ne.jkc) then
c     jkd is already ne jkb, and ne i 
                              if(cord(i,i_dir).eq.cord(jkd,i_dir)) then
c     found a neighbour of i in the plane of i-jkc-jkb. See if it is
c     connected to jkb, forming a closed figure
                                 jkb1=nelm(jkb)+1
                                 jkb2=nelm(jkb+1)
                                 do jkbb=jkb1,jkb2
                                    if(jkd.eq.nelm(jkbb)) then
c     closed figure is formed. By construction, we know that jkb and 
c     i are on the opposite sides of the line jkc-jkd: this is because 
c     jkc is the neighbour of i that is furthest from i in the la 
c     direction, and jkb is twice as far. Now check if 
c     jkc and jkd are on the opposite sides of the line i-jkb. Note that 
c     the line i-jkb is parallel to the la axis. First find l_perp, 
c     the axis normal to la(and l) and i_dir
                                       if(l.eq.1) then
                                          if(i_dir.eq.2) then
                                             l_perp=3
                                          elseif(i_dir.eq.3) then
                                             l_perp=2
                                          else
                     write(iptty,*)' Error in gotcord. l=i_dir=1'
                                             stop
                                          endif
                                       elseif(l.eq.2) then
                                          if(i_dir.eq.3) then
                                             l_perp=1
                                          elseif(i_dir.eq.1) then
                                             l_perp=3
                                          else
                     write(iptty,*)' Error in gotcord. l=i_dir=2'
                                             stop
                                          endif
                                       elseif(l.eq.3) then
                                          if(i_dir.eq.2) then
                                             l_perp=1
                                          elseif(i_dir.eq.1) then
                                             l_perp=2
                                          else
                     write(iptty,*)' Error in gotcord. l=i_dir=3'
                                             stop
                                          endif
                                       endif
                     dtemp1=cord(i,l_perp)-cord(jkd,l_perp)
                     dtemp2=cord(i,l_perp)-cord(jkc,l_perp)
                                       if((dtemp1*dtemp2).lt.0.) then

c     set getchord equal to cord of  node jkb
c     and exit the search loop. jkb is the neighbour of a neighbour, 
c     and not in the original neighbour list for i in nelm. So 
c     save jkb as the la-neighbour of i in irray(i,la)
c     irray(i,la)=jkb
c     irray(jkb,-la)=i
                                          gotcord=cord(jkb,l)
                                          kb=jkb
                                          goto 91912
                                       endif
                                    endif
                                 enddo
                              endif
                           endif
                        endif
                     enddo
                  endif
               endif              
91911          continue
            enddo
91912       continue
         endif
      endif
      
      return
      
      end
      
c..................................................................

      subroutine flag_boundry(i,i1,i2,i3,iflag_boundry)

      use combi
      use comdi, only : ps
      use comsk
      
      implicit none

      integer i,i1,i2,i3,iflag_boundry,i3ab,i3sign,kb,ksign
      integer k

      real*8 dist

c....s kelkar  feb 3, 04, 3D ORM stuff............
c flag OMR nodes that are on the exterior boundary. The subroutine
c flag_boundry checks orientations for missing nodes, and if they are
c present, then checks if any other connetcions exist on the same side 
c of the axis. If not, then it is a boundry node.

      i3ab=iabs(i3)
      i3sign=isign(1,i3)
      do k=i1,i2
         kb=nelm(k)
         if (ps(kb) .gt. 0.d0) then
            dist=cord(kb,i3ab)-cord(i,i3ab)
            ksign=dsign(1.d0,dist)
            if(ksign.eq.i3sign .and. abs(dist).gt.1.e-20) then
c found a neighbour node on the same side as the missing
c node. So i is not a boundry node. set flag and return
               iflag_boundry = 0
               goto 9999
            endif
         end if
      enddo
c did not find any neighbour nodes on the same side as the 
c missing node, so node i must be on an exterior boundry.
c set flag 
      iflag_boundry= +1

 9999 continue

      return

      end subroutine flag_boundry

c...................................................................

      
      subroutine  exit_ptrac1

      stop

      return

      end subroutine  exit_ptrac1

c...........................................................

      subroutine ddx_corn_array
c set up ddx, ddy,ddz and corn arrays
c NOTE: using ggg as scratch storage for saving del+ and del-
c temperorl

      use comai, only : neq, isptr9
      use combi, only : cord
      use comsptr
      use comsk

      implicit none

      integer i,j

      real*8 del_plus,del_minus,dtemp(3)

      do i=1,neq
         do j=1,3
            del_plus=ggg(i,j)
            del_minus=ggg(i,-j)
            dtemp(j)=del_plus+del_minus
            corn(i,j)=cord(i,j)-del_minus
         enddo
         ddx(i)=dtemp(1)
         ddy(i)=dtemp(2)
         ddz(i)=dtemp(3)
      enddo
c s kelkar sep 12 05 volume output for plumecalc
      if(sptrx_flag) then
         call sptr_volume_out
      endif

      return

      end subroutine ddx_corn_array

c....................................................................

      subroutine boundary_nodes(i,omrflag)
c....s kelkar  April 4, 2005
c flag nodes that are on the exterior boundary. The subroutine
c flag_boundry checks orientations for missing nodes, and if they are
c present, then checks if any other connetcions exist on the same side 
c of the axis. If not, then it is a boundry node.
c
c iboulist(i,7  )=# of boundary faces for node i (max 6)
c iboulist(i,1:6)= codes for  boundary faces of node i
c     irray(i,0) = +i : regular node, not a source/sink
c     irray(i,0) = -i-2000 : regular node on a external boundary
c     irray(i,0) = -i : regular interior node that is a sink/source
c                        but not explicitly specified in sptr macro
c     -100000000 < irray(i,0) < -10000000 : regular interior node that is 
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c                = -(i+2000000) : spring node
c                = -(i+1000000) : well-capture node on external bound
c                   simillar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c     -200000000 < irray(i,0) < -100000000 : non-OMR cliff node 
c     irray(i,0) < -200000000 : OMR cliff node 
c

      use comai, only : ierr, iptty
      use combi, only : nelm
      use comsk
      use comsptr

      implicit none

      integer i,ii1,ii2,i3,iflag_boundry,omrflag
      integer ibou,upper_limit

      ibou=0
      upper_limit=3
      ii1=nelm(i)+1
      ii2=nelm(i+1)

      do i3=1,upper_limit
         iflag_boundry=0
         if(irray(i,+i3).le.0) then
            call flag_boundry(i,ii1,ii2,+i3,iflag_boundry)
         endif
         if(iflag_boundry.eq.1) then
            if(irray(i,0).gt.-100000000) then
               if(omrflag.eq.1) then
                  irray(i,0)=-i-1000
               elseif(irray(i,0).ne.-(i+1000000).and.irray(i,0).ne.
     $                 -(i+2000000)) then
                  irray(i,0)=-i-2000
               endif
            endif
            ibou=ibou+1
            if(ibou.gt.6) then
               write(ierr, 1002)
               if (iptty .ne. 0) write(iptty, 1002)
               call exit_ptrac1
            endif
            if (omr_flag) iboulist(i,ibou)=i3
         endif
         iflag_boundry=0
         if(irray(i,-i3).le.0) then
            call flag_boundry(i,ii1,ii2,-i3,iflag_boundry)
         endif
         if(iflag_boundry.eq.1) then
            if(irray(i,0).gt.-100000000) then
               if(omrflag.eq.1) then
                  irray(i,0)=-i-1000
               elseif(irray(i,0).ne.-(i+1000000).and.irray(i,0).ne.
     $                 -(i+2000000)) then
                  irray(i,0)=-i-2000
               endif
            endif
            ibou=ibou+1
            if(ibou.gt.6) then
               write(ierr, 1002)
               if (iptty .ne. 0) write(iptty, 1002)
               call exit_ptrac1
            endif
            if (omr_flag) iboulist(i,ibou)=-i3
         endif
      enddo
 1002 format ('Error in stuc_conn_array, ibou > 6: STOPPING')
c     save the number of faces on the boundary
      if (omr_flag) iboulist(i,7)=ibou
      
      return
      
      end subroutine boundary_nodes
c...................................................


      subroutine subomr(i)

c s kelkar 1 onv 0
c determin which faces the omr nodes lie

      use comai
      use combi
      use comci
      use comdi
      use comflow
      use comsptr
      use comsk
      use compart

      implicit none

      integer i,i1,k,ia,kb,ikb,ipos

      real*8 epsilon,xkb,xia,d

      epsilon=1.e-8

      do i1=-3,3
         komr_count(i1)=0
      enddo

      do i1=1,iomr_count
         
c     NOTE: in iomr_neighbour,saved ikb rather than the node # kb 
c     because ikb can be directly used as a pointer to sx and a_axy 
c     arrays and kb can be retrieved from nelm(ikb)
         ikb=iomr_neighbour(i1)
         kb=nelm(ikb)
            
c     now figure out the geometry stuff

c     Look for type II connection, ie, kb is on a face defined
c     without a central node, going from higher to lower level omr.
c     we need to look at only those values of l for which irray(i,l)=0
c     a connection node is missing. Use the sign of the difference
c     in the l coordinate as an indicator
            
            do k=1,3            
               ia=abs(irray(i,k))
               if(ia.eq.0) then
                  d=cord(kb,k)-cord(i,k)
                  if(d.gt.epsilon) then
                     komr_count(k)=komr_count(k)+1
                     if(komr_count(k).gt.komrmax) then
                        write(ierr,*)'komr_count(k).gt.komrmax in ',
     &                       'subomr'
                        write(ierr,*)'change dimension in comomr.'
                        write(ierr,*)'i,iomr_count',i,iomr_count
                        write(ierr,*)'STOP'
                        if (iptty .ne. 0) then
                           write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                          'subomr'
                           write(iptty,*)'change dimension in comomr.'
                           write(iptty,*)'i,iomr_count',i,iomr_count
                           write(iptty,*)'STOP'
                        end if
                         call exit_ptrac1
                     endif
                     isave_omr(k,komr_count(k))=ikb
! commenting out the next goto 9191 to allow node kb to be included
! as a neighbour on multiple sides of i
c                     goto 9191
                  endif
               endif
               
               ia=abs(irray(i,-k))
               if(ia.eq.0) then
                  d=cord(kb,k)-cord(i,k)
                  if(d.lt.-(epsilon)) then
                     komr_count(-k)=komr_count(-k)+1
                     if(komr_count(-k).gt.komrmax) then
                        write(ierr,*)'komr_count(k).gt.komrmax in ',
     &                       'subomr'
                        write(ierr,*)'change dimension in comomr.'
                        write(ierr,*)'STOP'
                        if (iptty .ne. 0) then
                           write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                          'subomr'
                           write(iptty,*)'change dimension in comomr.'
                           write(iptty,*)'STOP'
                        end if
                         call exit_ptrac1
                     endif
                     isave_omr(-k,komr_count(-k))=ikb
c                     goto 9191
                  endif
               endif            
            enddo
            
c     check if its a type-I connection, ie if kb is on a face defined 
c     by a central node, going from lower to higher level omr.
c     check distance from each normal plane to kb to see if it is in the
c     plane
c     
            do k=1,3
               xkb=cord(kb,k)
               ia=abs(irray(i,k))
               if(ia.gt.0) then
                  xia=cord(ia,k)
                  d=abs(xkb-xia)
                  if(d.lt.epsilon) then
                     komr_count(k)=komr_count(k)+1
                     if(komr_count(k).gt.komrmax) then
                        write(ierr,*)'komr_count(k).gt.komrmax in ',
     &                       'subomr'
                        write(ierr,*)'change dimension in comomr.'
                        write(ierr,*)'STOP'
                        if (iptty .ne. 0) then
                           write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                          'subomr'
                           write(iptty,*)'change dimension in comomr.'
                           write(iptty,*)'STOP'
                        end if
                         call exit_ptrac1
                     endif
                     isave_omr(k,komr_count(k))=ikb
c                     goto 9191
                  endif
               endif
               
               ia=abs(irray(i,-k))
               if(ia.gt.0) then
                  xia=cord(ia,k)
                  d=abs(xkb-xia)
                  if(d.lt.-(epsilon)) then
                     komr_count(-k)=komr_count(-k)+1
                     if(komr_count(-k).gt.komrmax) then
                        write(ierr,*)'komr_count(-k).gt.komrmax in ',
     &                       'subomr'
                        write(ierr,*)'change dimension in comomr.'
                        write(ierr,*)'STOP'
                        if (iptty .ne. 0) then
                           write(iptty,*)'komr_count(k).gt.komrmax in ',
     &                          'subomr'
                           write(iptty,*)'change dimension in comomr.'
                           write(iptty,*)'STOP'
                        end if
                         call exit_ptrac1
                     endif
                     isave_omr(-k,komr_count(-k))=ikb
c                     goto 9191
                  endif
               endif
               
            enddo
                        
 9191       continue
            
      enddo
      
      return
      
      end
      
c.......................................................................

      subroutine wtsi_ptrac1_init
c     s kelkar aug 29, 05
c     for water table nodes, if S(ijkv(np1))>Smin then set the initial
c     position below ddz*S, if S<Smin then move the particle vertically
c     downward until a node with S>Smin is encountered.

      use comai, only : days
      use comdi, only : izone_free_nodes,s
      use comsptr
      use comsk

      implicit none

      integer i,inp1,newnode

      real*8 xp,yp,zp,dumm,zwt

      do i=1,num_part
         inp1=ijkv(i)
c     move only if it is time for the particle to move
         if (izone_free_nodes(inp1).gt.1 .and. ttp1(i) .le. days) then
            zp=zo(i)-corn(inp1,3)
            zwt=ddz(inp1)*s(inp1)
            if(zp.gt.zwt) then
               xp=xo(i)-corn(inp1,1)
               yp=yo(i)-corn(inp1,2)
               newnode=inp1
               call wtsi_find_water(inp1,i,xp,yp,zp,newnode)
               if (newnode .ne. 0) then
                  call wtsi_displace_node(inp1,i,xp,yp,zp,newnode)
                  ijkv(i)=newnode
                  xo(i)=xp+corn(newnode,1)
                  yo(i)=yp+corn(newnode,2)
                  zo(i)=zp+corn(newnode,3)
               end if
            endif
         endif
      enddo

      end subroutine wtsi_ptrac1_init

c....................................................................


      subroutine wtsi_find_water(inp1,np1,xp,yp,zp,newnode)
      
c     s kelkar  aug 30, 05
c     If newnode has irreducible water
c     saturation, flagged by izone_free_nodes(inp1).gt.1 then
c     search vertically downward to find a node with flowing water.

      use comai, only : ierr, iptty 
      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,newnode,j,nextnode, node_flag, nodetemp
      integer node_previous,ibou,i1
      
      real*8 xp,yp,zp
      real*8  epsilon

      epsilon= 1.e-10
      
      if (izone_free_nodes(newnode).ge.3.or.s(newnode).lt.smin) then
         node_previous=newnode
         do j=1,1000000
            node_flag=irray(node_previous,0)
            if(node_flag.eq.0) then
c     nextnode is interior OMR, do OMR check
c     xc and yc are updated wrt nodetemp in wtsi_neighbour and also
c     zc is set on the +3 boundary of nodetemp.
               call wtsi_neighbour2(inp1,np1,node_previous,
     1              xp,yp,zp, nodetemp) 
               if(izone_free_nodes(nodetemp).le.2.and.s(nodetemp)
     1              .ge.smin) then
c     found a valid flowing node,  return
                  newnode=nodetemp
                  goto 99999
               else
                  node_previous = nodetemp
               endif
            elseif(node_flag.eq.-(1000+node_previous)) then
c     node_previous is a boundary OMR node. check boundary faces,
c     stored in iboulist(), to see if the -3 plane is a
c     boundary plane. if so, particle is exiting the model- return.
c     If it is not, find the new nearest node.
               ibou=iboulist(node_previous,7)
               do i1=1,ibou
                  if(-3.eq.iboulist(node_previous,i1)) then
c particle exited the model
                     newnode=0
                     write(ierr,*)'ptrac1:find_water_table. cant find'
                     write(ierr,*)'water table. Check initial particle'
                     write(ierr,*)'locations. np1=',np1
                     write(ierr,*)'exit_ptrac1. stop.'
                     call exit_ptrac1
c                     goto 99999
                  endif
               enddo
c     node_previous not on -3 boundary. In wtsi_neighbour,
c     xc and yc are updated wrt nodetemp in wtsi_neighbour and also
c     zc is set on the +3 boundary of nodetemp.
               call wtsi_neighbour2(inp1,np1,node_previous,
     1              xp,yp,zp, nodetemp) 
               if(izone_free_nodes(nodetemp).le.2.and.s(nodetemp)
     1              .ge.smin) then
c     found a valid flowing node,  return
                  newnode=nodetemp
                  goto 99999
               else
                  node_previous = nodetemp
               endif
            else
c     regular interior node
               nextnode=irray(node_previous,-3)
               if(nextnode.gt.0) then
c     wtsi_neighbour has not been called, so need to update xc,yc,zc
c     change xc,yc to refere to nextnode and zc slightly below
c     + 3 face of nextnode. 
                  xp=xp+corn(node_previous,1)-corn(nextnode,1)
                  yp=yp+corn(node_previous,2)-corn(nextnode,2)
                  zp=(1.-epsilon)*ddz(nextnode)
                  If(izone_free_nodes(nextnode).le.2.and.s(nextnode)
     1              .ge.smin) then
c     found a valid flowing node,  return
                     newnode=nextnode
                     goto 99999
                  else
c     nextnode doesnt have flowing water,continue checking -3 neighbour
                     node_previous = nextnode
                  endif
                  
               elseif(nextnode.eq.0) then
c     bottom node- particle has to exit
                  newnode=0
               else
                  write(iptty,*)'error wtsi_newnode.nextnode=',nextnode
                  write(iptty,*)'STOP'
                  write(ierr,*)'error wtsi_newnode.nextnode=',nextnode
                  write(ierr,*)'STOP'
                  stop
               endif
            endif
   
         enddo
         
      else
c     found a valid flowing node,  return

      endif
               
99999 continue
      
      return
      
      end subroutine wtsi_find_water

c...........................................................
      
      subroutine wtsi_neighbour2(inp1,np1,node_previous,xc,yc,zc,
     1     nodetemp)    

      use comai, only : ierr, iptty 
      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,j, nodetemp, node_previous
      
      real*8 xc,yc,zc,epsilon,zcc

      epsilon=1.e-8

c     an Internal OMR node, check for neighbours on -3 side 
c     to begin ,place zcc slightly below -3 face of node_previous

      zcc=-epsilon*ddz(node_previous)

      call nearest_node(node_previous,np1,-3,xc,yc,zcc, nodetemp)

      if(nodetemp.le.0) then
         write(ierr,*)'Hole in the model. node=',node_previous
         write(ierr,*)'stop in wtsi_neighbour2'
         if (iptty .ne. 0) then
            write(iptty,*)'Hole in the model. node=',node_previous
            write(iptty,*)'stop in wtsi_neighbour2'
         endif
         call update_exit(-inp1,np1,-100,nodetemp,
     $     0.,xc,yc,zc)
c         stop
      end if   

c     change xc,yc to refere to nodetemp 
c     place zc at +3 face of nodetemp
      xc=xc+corn(node_previous,1)-corn(nodetemp,1)
      yc=yc+corn(node_previous,2)-corn(nodetemp,2)
      zc=(1.-epsilon)*ddz(nodetemp)

      return
      
      end subroutine wtsi_neighbour2

c...........................................................

      subroutine wtsi_newlocation2(nodeabove,np1,nodebelow,zc)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer nodeabove,np1,nodebelow,j
      
      real*8 dzw,zc, epsilon

      epsilon=1.e-12

      dzw=s(nodebelow)*ddz(nodebelow)
      if((dzw+corn(nodebelow,3)).lt.(ddz(nodeabove)+corn(nodeabove,3)))
     1     then
c     at this point zc is assumed to be wrt nodebelow already.
         if(zc.gt.dzw) then
            zc=dzw*(1.-deltawt)
         endif
      endif

      return
      
      end subroutine wtsi_newlocation2
c...........................................................

      subroutine wtsi_displace_node(inp1,np1,xp,yp,zp,newnode)
      
c     s kelkar  aug 30, 05
c     displace particle location vertically downward to deltawt
c     meters below the water table and find the new node 
      use comai, only : ierr, iptty 
      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,newnode,j,nextnode, node_flag, nodetemp
      integer node_previous,ibou,i1
      
      real*8 xp,yp,zp,zwt,zptemp
      real*8  epsilon

      epsilon= 1.e-10

c     place the particle deltawt(m) below the water table. note that
c     xp,yp,zp are wrt newnode
      zwt=ddz(newnode)*s(newnode)
      zp=zwt*(1-epsilon)-abs(deltawt)
      
      if (zp.lt.0.) then
         node_previous=newnode
         do j=1,1000000
            node_flag=irray(node_previous,0)
            if(node_flag.eq.0) then
c     nextnode is interior OMR, do OMR check
c     xp,yp,zp are updated wrt nodetemp in wtsi_neighbour3
c     zp might be far below corn of node_previous, in which
c     case nearest_node  might not find nodetemp in level_max
c     iterations, so use zptemp to place the particle 
c     just below corn(node_previous,3)
               zptemp=-epsilon
               call wtsi_neighbour3(inp1,np1,node_previous,
     1              xp,yp,zptemp, nodetemp)
c     recalculate zp wrt nodetemp
               zp=zp+corn(node_previous,3)-corn(nodetemp,3)
               if(zp.gt.0.) then
c     found the new CC,  return
                  newnode=nodetemp
                  goto 99999
               else
                  node_previous = nodetemp
               endif
            elseif(node_flag.eq.-(1000+node_previous)) then
c     node_previous is a boundary OMR node. check boundary faces,
c     stored in iboulist(), to see if the -3 plane is a
c     boundary plane. if so, particle is exiting the model- return.
c     If it is not, find the new nearest node.
               ibou=iboulist(node_previous,7)
               do i1=1,ibou
                  if(-3.eq.iboulist(node_previous,i1)) then
c particle exited the model
                     newnode=0
                     write(ierr,*)'ptrac1:find_water_table. cant find'
                     write(ierr,*)'water table. Check initial particle'
                     write(ierr,*)'locations. np1=',np1
                     write(ierr,*)'exit_ptrac1. stop.'
                     call exit_ptrac1
c                     goto 99999
                  endif
               enddo
c     node_previous not on -3 boundary. In wtsi_neighbour,
c     xc and yc are updated wrt nodetemp  and also
c     zc is set on the +3 boundary of nodetemp.
               zptemp=-epsilon
               call wtsi_neighbour3(inp1,np1,node_previous,
     1              xp,yp,zptemp, nodetemp) 
c     recalculate zp wrt nodetemp
               zp=zp+corn(node_previous,3)-corn(nodetemp,3)
               if(zp.gt.0.) then
c     found the new CC,  return
                  newnode=nodetemp
                  goto 99999
               else
                  node_previous = nodetemp
               endif
            else
c     regular interior node
               nextnode=irray(node_previous,-3)
               if(nextnode.gt.0) then
c     wtsi_neighbour has not been called, so need to update xc,yc,zc
c     change xc,yc to refere to nextnode and zc slightly below
c     + 3 face of nextnode. 
                  xp=xp+corn(node_previous,1)-corn(nextnode,1)
                  yp=yp+corn(node_previous,2)-corn(nextnode,2)
                  zp=zp+corn(node_previous,3)-corn(nextnode,3)
                  If(zp.gt.0.) then
c     found the new CC,  return
                     newnode=nextnode
                     goto 99999
                  else
c     not in CC of nextnode  ,continue checking -3 neighbour
                     node_previous = nextnode
                  endif
                  
               elseif(nextnode.eq.0) then
c     bottom node- particle displacement outside the model
c     place the particle in the last CC
                  newnode=node_previous
               else
                  write(iptty,*)'error wtsi_displace_node.nextnode=',
     1                 nextnode
                  write(iptty,*)'STOP'
                  write(ierr,*)'error wtsi_displace_node.nextnode=',
     1                 nextnode
                  write(ierr,*)'STOP'
                  stop
               endif
            endif
   
         enddo
         
      else
c     found the new CC,  return

      endif
               
99999 continue
      
      return
      
      end subroutine wtsi_displace_node

c...........................................................
      
      subroutine wtsi_neighbour3(inp1,np1,node_previous,xc,yc,zc,
     1     nodetemp)    

      use comai, only : ierr, iptty 
      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,j, nodetemp, node_previous
      
      real*8 xc,yc,zc,epsilon,zcc

      epsilon=1.e-8

c     an Internal OMR node, check for neighbours on -3 side 

      call nearest_node(node_previous,np1,-3,xc,yc,zc, nodetemp)

      if(nodetemp.le.0) then
         write(ierr,*)'Hole in the model. node=',node_previous
         write(ierr,*)'stop in wtsi_neighbour'
         write(iptty,*)'Hole in the model. node=',node_previous
         write(iptty,*)'stop in wtsi_neighbour'
         stop
      endif   

c     change xc,yc,zc to refere to nodetemp 

      xc=xc+corn(node_previous,1)-corn(nodetemp,1)
      yc=yc+corn(node_previous,2)-corn(nodetemp,2)
      zc=zc+corn(node_previous,3)-corn(nodetemp,3)

      return
      
      end subroutine wtsi_neighbour3

c...........................................................

      subroutine sptr_volume_out
c sep 13 05 s kelkar
c store volumes for plumecalc

      use comai, only : iw, neq, neq_primary, isptr9
      use combi
      use comsptr
      use comsk
      use comxi, only : cform
      implicit none

      logical opened

      integer i1, i1flag, i, ncoef, max_con

      real*8, allocatable :: sx1temp(:)

      allocate (sx1temp(neq))

      do i1=1,neq
         i1flag=irray(i1,0)
         if(i1flag.eq.0.or.i1flag.eq.-(i1+1000).or.
     1        (i1flag.lt.-200000000)) then
c omr node- compute the volume of the approximate brick shape
            sx1temp(i1)=ddx(i1)*ddy(i1)*ddz(i1)
         else
            sx1temp(i1)=sx1(i1)
         endif
      enddo

      ncoef=0
      max_con=0
      
      if(cform(26).eq.'formatted') then
c   formatted
         write(isptr9, '(5i10)' ) iw, neq_primary, 
     &        nelm(neq_primary+1), ncoef, max_con
         write(isptr9, '(5(1pe20.10))') 
     &        (sx1temp(i), i = 1, neq_primary)
         write(isptr9, '(5i10)'  ) 
     &        (nelm(i), i = 1, nelm(neq_primary+1))
      elseif(cform(26).eq.'unformatted') then
c unformatted
         write(isptr9) iw, neq_primary, nelm(neq_primary+1), 
     &        ncoef, max_con
         write(isptr9) 
     &        (sx1temp(i), i = 1, neq_primary)
         write(isptr9) 
     &        (nelm(i), i = 1, nelm(neq_primary+1) )

      endif

      if(allocated(sx1temp)) deallocate(sx1temp)      

      close(isptr9)

      return
      
      end subroutine sptr_volume_out

******************************************************************

      subroutine init_sptr_params

      use comai, only : neq
      use combi, only : izonef
      use comci, only : rolf
      use comdi, only : denr, diffmfl, ifree, itrc, ps_trac, s
      use compart, only : aperture, kd, matrix_por, secondary
      use comsptr
      use davidi, only : irdof
      implicit none

      integer i, np1
      real*8 denominator, rprime, spacing
     
c     Subroutine to initialize particle tracking parameters

      if(nzbtc.gt.0) then
c zvd - initialized in allocmem and used in insptr, don't reset here
c         izonebtc = 0
      end if

cHari 01-Nov-06 include colloid diversity model (tprpflag=11)
      do i = 1, neq
         if(tprpflag(itrc(i)).eq.1.or.tprpflag(itrc(i)).eq.2.or.
     2        tprpflag(itrc(i)).eq.11) then

c     Compute denominator, make sure no divide by 0

            if (irdof .ne. 13 .or. ifree .ne. 0) then
               denominator = rolf(i)*s(i)*ps_trac(i)
            else
               denominator = rolf(i)*ps_trac(i)
            endif
            denominator = max(1.d-30, denominator)
            omega_partial(i) = 1.+kd(itrc(i),1)*
     2           denr(i)/denominator

         elseif(tprpflag(itrc(i)).eq.3.or.tprpflag(itrc(i)).eq.4) then
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               denominator = rolf(i)*s(i)*matrix_por(itrc(i))
            else
               denominator = rolf(i)*matrix_por(itrc(i))
            endif
            denominator = max(1.d-30, denominator)
            rprime = 1.+kd(itrc(i),1)*
     2           denr(i)/denominator

c     aperture = 2 * b in Sudicky and Frind solution
c     spacing = 2 * B in Sudicky and Frind solution

            if(aperture(itrc(i)).lt.0.) then
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  omega_partial(i) = (matrix_por(itrc(i))*s(i))**2*
     2                 diffmfl(1,itrc(i))*rprime
               else
                  omega_partial(i) = (matrix_por(itrc(i)))**2*
     2                 diffmfl(1,itrc(i))*rprime
               endif
               sigma_partial(i) = aperture(itrc(i))
            else
               if (secondary(itrc(i)) .ne. 0.) then
                  spacing = secondary(itrc(i))
               else
                  spacing = abs(aperture(itrc(i)))/max(1.d-30,
     &                 ps_trac(i))
               end if
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  omega_partial(i) = s(i)*matrix_por(itrc(i))*
     2                 sqrt(rprime*diffmfl(1,itrc(i)))/(0.5*
     3                 abs(aperture(itrc(i))))
               else
                  omega_partial(i) = matrix_por(itrc(i))*
     2                 sqrt(rprime*diffmfl(1,itrc(i)))/(0.5*
     3                 abs(aperture(itrc(i))))
               endif
               sigma_partial(i) = sqrt(rprime/diffmfl(1,itrc(i)))*0.5*
     3              (spacing-abs(aperture(itrc(i))))
            end if
         end if
      end do

      return
      end subroutine init_sptr_params

******************************************************************
