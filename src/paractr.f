      subroutine paractr(iflg)
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
!D1 Partition the zones into separate pieces for serial processing.
!D1 Eventually the zones will be shared between parallel processors.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 21-MAY-03, Programmer: P. Fasel
!D2
!D2 $Log:   /pvcs.config/fehm90/src/paractr.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
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

      use davidi
      use comai
      use combi
      use comdi
      use comdti
      use comgi
      use com_part
      implicit none

      integer i,i1,i2,j,ii,kb,neqp1,it, jpara
      integer nsizea_tmp,nsizeb_tmp
      integer zone, node, new_node
      integer last_conn, first_conn
      integer pos(4)
      real*8 strd_para, aiped_mult
      integer old_node
      integer nindex, neighbor ,mindex
      integer first_conn_new,last_conn_new
      integer neighbor_new,neighbor_new_old
      integer iflg, neqp1_part, igr
      integer max_zone_nodes    
      integer ic, inode, izone, nzone_para_dum, nzone_dum     
      integer maxcon, nzone_max             
      parameter (maxcon = 10, nzone_max = 100)
      parameter (strd_para=1.000, aiped_mult=1.d08)

c================================================================
      if(ipara.eq.0) return
c================================================================
      

      if(iflg.eq.0) then
c     
c     read input and open output file for zone fluxes      
c     
         if(nzone_para.eq.0) then
            open(101,file='zone_part_flux',status='unknown')
            write(101,601) 
 601        format (1x,'Flux budget for zone partitions')

            ipara = 1
            igrow = 0
            npall = n0

            allocate (bound_zone(neq))
            bound_zone = 0

            read(inpt,*) nzone_para
            if(.not.allocated(izone_para)) then
               allocate(izone_para(max(1,nzone_para)))
               allocate(igrow_part(max(1,nzone_para)))
               allocate(ipart_ini(max(1,nzone_para)))
               allocate(decomp_part(npall))
               allocate(zone_level(nzone_max))
            end if
            backspace inpt 
            read(inpt,*) nzone_para, (izone_para(i),
     &           igrow_part(i),ipart_ini(i),i=1,nzone_para)

c     
c     Loop over each zone for determining izone_para array
c     
            allocate (neq_part(nzone_max))
            max_zone_nodes = 0
            decomp_part = 0
            
            do izone = 1, nzone_para
               ic = 0
               zone_level(izone) = ipara
               if(igrow_part(izone).ne.0) igrow = 1
               do inode = 1, n0
                  if(izonef(inode).eq.izone_para(izone)) then
                     decomp_part(inode) = izone_para(izone)
                     ic = ic +1
                  end if
               end do
               neq_part(izone) = ic
               max_zone_nodes = max(max_zone_nodes,ic) 
            end do
         else
            ipara=ipara+1
            nzone_para0 = nzone_para
            npall=npall + n0
            read(inpt,*) nzone_para
            nzone_para = nzone_para+nzone_para0
            allocate(izone_para_temp(max(1,nzone_para)))
            allocate(igrow_part_temp(max(1,nzone_para)))
            allocate(ipart_ini_temp(max(1,nzone_para)))
            allocate(decomp_part_temp(npall))
            izone_para_temp(1:nzone_para0)= izone_para(1:nzone_para0)
            decomp_part_temp(1:npall-n0)= decomp_part(1:npall-n0)
            igrow_part_temp(1:nzone_para0)= igrow_part(1:nzone_para0)
            ipart_ini_temp(1:nzone_para0)= ipart_ini(1:nzone_para0)
            deallocate(izone_para,igrow_part,decomp_part)
            deallocate(ipart_ini)
            allocate(izone_para(nzone_para))
            allocate(igrow_part(nzone_para))
            allocate(ipart_ini(nzone_para))
            allocate(decomp_part(npall))
            izone_para = izone_para_temp
            igrow_part = igrow_part_temp
            ipart_ini = ipart_ini_temp
            decomp_part = decomp_part_temp
            deallocate(izone_para_temp,igrow_part_temp)
            deallocate (decomp_part_temp,ipart_ini_temp)
            backspace inpt
            read(inpt,*) nzone_dum,
     &           (izone_para(i),igrow_part(i),ipart_ini(i),
     &           i=nzone_para0+1,nzone_para)
c     
            decomp_part(n0*(ipara-1)+1:n0*ipara) = 0
c     
c     Loop over each zone for determining izone_conv array
c     
            do izone = nzone_para0+1, nzone_para
               ic = 0
               zone_level(izone) = ipara
               if(igrow_part(izone).ne.0) igrow = 1
               do inode = 1, n0
                  if(izonef(inode).eq.izone_para(izone)) then
                     decomp_part(inode+n0*(ipara-1)) = 
     &                    izone_para(izone)
                     ic = ic +1
                  end if
               end do
               neq_part(izone) = ic
               max_zone_nodes = max(max_zone_nodes,ic) 
            end do
         endif

      else if(iflg.eq.1) then

c     
c     first allocate maximum space for index_part
c**** initialize the reverse index because 0 will mean the old node
c**** number is in another zone and shouldn't be included
c     
         allocate (index_part(nzone_para,n0))
         allocate (rindex_part(nzone_para,n0))

         
         rindex_part = 0
         neq_part = 0

         if(igrow.ne.0) then

            allocate (decomp_part_temp(npall))
            allocate (decomp_part_temp1(npall))
            decomp_part_temp = 0
            decomp_part_temp1 = 0

c**** pass 1 builds the nelm for this partition by allowing space
c**** and filling in index and rindex references to all nodes
            do zone = 1, nzone_para
               ipara = zone_level(zone)
               decomp_part_temp(n0*(ipara-1)+1:n0*ipara) = 
     &              decomp_part(n0*(ipara-1)+1:n0*ipara)
               decomp_part_temp1(n0*(ipara-1)+1:n0*ipara) = 
     &              decomp_part(n0*(ipara-1)+1:n0*ipara)
               izone = izone_para(zone)
               if(igrow_part(zone).ne.0) then
                  do igr = 1,igrow_part(zone)
                     do node = 1, neq
                        if (decomp_part_temp(node+n0*(ipara-1)) == 
     &                       izone) then
c     grow zone neighbor list larger
                           i1 = nelm(node)+1
                           i2 = nelm(node+1)
                           do ii = i1,i2
                              kb = nelm(ii)
                              decomp_part_temp1(kb+n0*(ipara-1)) = izone
                           enddo
                        endif
                     enddo
                     do ii = 1,neq
                        decomp_part_temp(ii+n0*(ipara-1)) = 
     &                       decomp_part_temp1(ii+n0*(ipara-1))
                     enddo
                  enddo
                  do node = 1, neq 
                     if (decomp_part_temp(node+n0*(ipara-1)) == izone) 
     &                    then
                        neq_part(zone) = neq_part(zone) + 1
                        index_part(zone, neq_part(zone)) = node
                        rindex_part(zone, node) = neq_part(zone)
                     endif
                  enddo
               else
                  do node = 1, neq 
                     if (decomp_part_temp(node+n0*(ipara-1)) == izone) 
     &                    then
                        neq_part(zone) = neq_part(zone) + 1
                        index_part(zone, neq_part(zone)) = node
                        rindex_part(zone, node) = neq_part(zone)
                     endif
                  enddo
               endif
            enddo

            deallocate(decomp_part_temp, decomp_part_temp1)

c     endif on growth of zones
         else

c**** pass 1 builds the nelm for this partition by allowing space
c**** and filling in index and rindex references to all nodes
            do zone = 1, nzone_para
               izone = izone_para(zone)
               ipara = zone_level(zone)
               do node = 1, neq
                  if (decomp_part(node+n0*(ipara-1)) == izone) then
                     neq_part(zone) = neq_part(zone) + 1
                     index_part(zone, neq_part(zone)) = node
                     rindex_part(zone, node) = neq_part(zone)
                  endif
               enddo
            enddo

         endif

c     
c     find maximum sizes of zones, resize some arrays
c     
         max_zone_nodes =0
         do zone = 1, nzone_para
            max_zone_nodes = max(max_zone_nodes,neq_part(zone))
         enddo
c     
c     note repeated nodes and for nodes with no zone partition
c     
         allocate (repeat_part(n0))
	 repeat_part = 0.0
         do zone = 1,nzone_para
            do i = 1,neq_part(zone)
               kb = index_part(zone,i)
               repeat_part(kb) = repeat_part(kb) + 1.0
            enddo
         enddo
         ic =0
         do i = 1,neq
            if(repeat_part(i).gt.0.0) then
               repeat_part(i) = 1.0/repeat_part(i)*strd_para
            else
               ic = ic + 1
            endif
         enddo
         nzone_para0 = nzone_para

         if(ic.ne.0) then
            nzone_para = nzone_para + 1 
            max_zone_nodes = max(max_zone_nodes,ic)
            allocate(index_part_temp(nzone_para,max_zone_nodes))
            allocate (rindex_part_temp(nzone_para,n0))      
            do zone = 1, nzone_para0
               index_part_temp(zone,1:max_zone_nodes) = 
     &              index_part(zone,1:max_zone_nodes)
               rindex_part_temp(zone,1:n0) = rindex_part(zone,1:n0)
            enddo
c     must resize other arrays
            
            ic = 0
            izone = nzone_para
            do i = 1,neq
               if(repeat_part(i).eq.0.0) then
                  repeat_part(i) = 1.0
                  ic = ic +1
                  index_part_temp(izone,ic) = i
                  rindex_part_temp(izone, i) = ic
               endif
            enddo
            neq_part(izone) = ic
            deallocate(index_part)
            allocate(index_part(nzone_para,max_zone_nodes))
            index_part = index_part_temp
            deallocate(index_part_temp)
            allocate(izone_para_temp(nzone_para))
            izone_para_temp(1:nzone_para-1) = izone_para(1:nzone_para-1)
            deallocate(izone_para)
            allocate(izone_para(nzone_para))
            izone_para = izone_para_temp
            izone_para(nzone_para) = -1
            deallocate(izone_para_temp)
            deallocate(rindex_part)
            allocate (rindex_part(nzone_para,n0))
            rindex_part = rindex_part_temp
            deallocate(rindex_part_temp)
         endif

c     
c     Allocate some more arrays
c     

         allocate (iad_part(nzone_para))
         allocate (itert_part(nzone_para))
         allocate (itotal_part(nzone_para))
         allocate (itotals_part(nzone_para))
         allocate (timing_part(nzone_para))
         allocate (fdum_part(nzone_para))
         allocate (f0_part(nzone_para))
         allocate (bound_pos(nzone_para))
         allocate (nmat_part(nzone_para,36))
         allocate (nb_part(nzone_para,36))
         allocate (nrhs_part(nzone_para,36))

         allocate (nelmdg_part(nzone_para,max_zone_nodes))
         allocate (npvt_part(nzone_para,max_zone_nodes))
         allocate (bp_part(nzone_para,2*max_zone_nodes))
         allocate (nelm_part(nzone_para, max_zone_nodes*maxcon))
         allocate (bound_part(nzone_para, max_zone_nodes*maxcon))
         allocate (flux_l_part(nzone_para, max_zone_nodes*maxcon))
         allocate (flux_v_part(nzone_para, max_zone_nodes*maxcon))
         allocate (bound_istrw(nzone_para, max_zone_nodes*maxcon))
         allocate (istrw_part(nzone_para, max_zone_nodes*maxcon))

         itotal_part = 0
         itotals_part = 0

c**** pass 2
         do zone = 1, nzone_para
            pos(zone) = neq_part(zone) + 1
            bound_pos(zone) = neq_part(zone) + 1
            nelm_part(zone,1) = pos(zone)
            bound_part(zone,1) = bound_pos(zone)
            neqp1_part = neq_part(zone)+1
            do new_node = 1, neq_part(zone)
               
c**** get the index of the first and last connections
               old_node = index_part(zone, new_node)
               first_conn = nelm(old_node) + 1
               last_conn = nelm(old_node + 1)

               do nindex = first_conn, last_conn

                  neighbor = nelm(nindex)
                  
c**** is the connection in the new node list
                  if (rindex_part(zone, neighbor) /= 0) then
                     pos(zone) = pos(zone) + 1
                     nelm_part(zone, pos(zone)) = 
     &                    rindex_part(zone, neighbor)

                  else

c**** this connection is a boundary
                     bound_pos(zone) = bound_pos(zone) + 1
                     bound_part(zone, bound_pos(zone)) = neighbor
                     bound_istrw(zone, bound_pos(zone)-neqp1_part)
     &                    = istrw(nindex-(neq+1))

                  endif
               enddo

c**** fill in the index of the last connection
               nelm_part(zone, new_node + 1) = pos(zone)
               bound_part(zone, new_node + 1) = bound_pos(zone)

            enddo

         enddo

         do zone=1,nzone_para
            neqp1=neq_part(zone)+1
            nsizea_tmp = nelm_part(zone,neqp1) - neqp1
            nsizeb_tmp = nelm_part(zone,neqp1) - neqp1

            nmat_part(zone,1) = 0
            nb_part(zone,1) = 0
            
            do i = 2, 36
               nb_part(zone,i) = nb_part(zone,i-1) + nsizeb_tmp
               nmat_part(zone,i) = nmat_part(zone,i-1) + nsizea_tmp
            enddo

            do i = 1, 6
               nrhs_part(zone,i) = (i-1) * neq_part(zone)
            enddo
         enddo

c**** calculate the nelmdg_part for each zone
         do zone = 1, nzone_para
            do node = 1, neq_part(zone)
               loop: do j = nelm_part(zone,node)+1, 
     &              nelm_part(zone,node+1)
               if (nelm_part(zone,j) .eq. node) then
                  nelmdg_part(zone,node) = j
                  exit loop
               endif
            enddo loop
         enddo
      enddo

c**** calculate npvt_part for each zone
      do zone = 1, nzone_para
         do node = 1, neq_part(zone)
            i1 = nelm_part(zone,node) + 1
            i2 = nelm_part(zone,node+1)
            do i = i1, i2
               if (nelm_part(zone,i) .eq. node) then
                  npvt_part(zone,node) = i
                  goto 100
               endif
            enddo
 100        continue
         enddo
      enddo

c**** calculate istrw_part for each zone
      istrw_part = 0
      do zone = 1, nzone_para
         do new_node = 1, neq_part(zone)

            first_conn_new = nelm_part(zone,new_node) + 1
            last_conn_new = nelm_part(zone,new_node + 1)
            
            old_node = index_part(zone, new_node)
            first_conn = nelm(old_node) + 1
            last_conn = nelm(old_node + 1)

            do mindex = first_conn_new, last_conn_new
               neighbor_new = nelm_part(zone,mindex)
               neighbor_new_old = index_part(zone, neighbor_new)

               do nindex = first_conn, last_conn
                  neighbor = nelm(nindex)
                  if(neighbor.eq.neighbor_new_old) then
                     iw = istrw(nindex-(neq+1))
                     istrw_part(zone,mindex-
     &                    (neq_part(zone)+1)) = iw
                     go to 200
                  endif
               enddo
 200           continue
            enddo
         enddo
      enddo
      else if(iflg.eq.2) then
         do zone = 1, nzone_para
            do new_node = 1, neq_part(zone)
               old_node = index_part(zone,new_node) 
               bp(old_node+nrhs(1)) = 
     &              bp_part(zone,new_node+nrhs_part(zone,1))
               bp(old_node+nrhs(2)) = 
     &              bp_part(zone,new_node+nrhs_part(zone,2))
            enddo
         enddo
      else if(iflg.eq.3) then
c     
c     printout partition information
c     
         write(iout,501)  nzone_para
         if(iatty.gt.0) then
            write(iatty,501)  nzone_para
         endif
 501     format(/,1x,'Partition statistics, partitions = ',i6)
         
         do zone = 1, nzone_para
            izone = izone_para(zone)
            write(iout,502) izone, neq_part(zone), iad_part(zone),
     &           itotal_part(zone),itert_part(zone),itotals_part(zone),
     &           timing_part(zone)
            if(iatty.gt.0) then
               write(iatty,502) izone, neq_part(zone), iad_part(zone),
     &              itotal_part(zone),itert_part(zone),
     &              itotals_part(zone),timing_part(zone)
            endif
         enddo
 502     format(1x,'zone, nodes',1x,',N-R iters,Total N-R iters,',
     &        ' Solver iters, Total solver iters, CPU (time step)',/,
     &        i5,1x,i6,6x,i5,10x,i6,9x,i5,15x,i5,6x,g10.2)

      else if(iflg.eq.4) then
c     
c     printout partition flux info
c     
         write(101,701) l,days
 701     format (/,'time step ',i7,' time(days) ',g13.5)
	 do zone = 1, nzone_para
            izone = izone_para(zone)
            do node = 1, neq_part(zone)
               i1=bound_part(zone,node)+1
               i2=bound_part(zone,node+1)
               if(i2.ge.i1) then
                  do ii=i1,i2
                     kb = bound_part(zone,ii)
                     if(kb.ne.0) then
                        write (101,503) izone, index_part(zone,node), 
     &                       kb,flux_l_part(zone, node),
     &                       flux_v_part(zone,node)
                     endif
                  enddo
               endif
            enddo
	 enddo
         
 503     format(1x,'zone ',i4, ' zone node ', i7,' outside node ',i7,
     &        ' flux liquid(kg/s) ',g13.5,' flux vapor(kg/s) ',g13.5)

      else if(iflg.eq.5) then
c     
c     apply  head BC on child grid
c     
         i2 = 0
         do zone = 1, nzone_para
            jpara = ipart_ini(zone)
            if(jpara.ne.0) then
               ipart_ini(zone) = 1
               do node = 1, neq_part(zone)
                  
                  kb = index_part(zone,node)
                  if(ka(kb).eq.-6.or.bound_zone(kb).eq.-6) then
                     i2 = i2 +1 
                     ka(kb) = -1
                     wellim(kb) = sx1(kb)*aiped_mult
                     esk(kb) = 1.0
                     bound_zone(kb) = -6
                  endif
               enddo
            endif
         enddo
         i2 = 0
         do i = 1,neq
            if(bound_zone(i).eq.-6) then	
               i2 = i2 +1 
            endif
         enddo 

c     
         if (iptty .ne. 0) write (iptty,*) 'couple BC changes', i2
         if (iout .ne. 0) write (iout,*) 'couple BC changes', i2  
      else if(iflg.eq.6) then
c     remove head BC on child grid
         
c     write (*,*) 'iad ', iad,' fdum ', fdum
c     write(*,*) 'well',  nskw(1), phi(nskw(1))
c     write (iout,*) 'iad ', iad,' fdum ', fdum
c     write(iout,*) 'well',  nskw(1), phi(nskw(1))
c     return
c     write(iout,*) iad, fdum
c     write(iout,*) 'node 22091 ',phi(22091), bp(22091)
c     write(iout,*) 'node 22005 ',phi(22005), bp(22005)
c     write(iout,*) 'node 22006 ',phi(22006), bp(22006)
c     write(iout,*) 'node 27121 ',phi(27121), bp(27121)
c     write(iout,*) 'node 3114 ',phi(3114), bp(3114)
         do zone = 1, nzone_para
            jpara = ipart_ini(zone)
            if(jpara.gt.0) then
               ipart_ini(zone) = -1
               do node = 1, neq_part(zone)
                  kb = index_part(zone,node)
                  if(bound_zone(kb).eq.-6) then
                     ka(kb) = 0
                  endif
               enddo
            endif
         enddo
         
      endif

      end
