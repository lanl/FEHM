      subroutine user(k)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine can be invoked at the beginning of every time 
CD1  step to change or write-out any variable in module. 
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY (ECD Number 22)
CD2
CD2 $Log:   /pvcs.config/fehmn/src/user_ymp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
CD2 
CD2    Rev 1.6   Fri Feb 09 09:06:36 1996   robinson
CD2 Added option for cannister heat and episodic infil.
CD2 
CD2    Rev 1.5   Thu Jan 11 12:11:00 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   12/13/95 10:31:50   robinson
CD2 Added radioactive decay heat curve
CD2 
CD2    Rev 1.3   05/16/95 09:35:16   robinson
CD2 Added ngas capability to the zone infiltration option
CD2 
CD2    Rev 1.2   01/27/95 15:56:48   llt
CD2 implemented new infiltration boundary condition
CD2 
CD2    Rev 1.1   03/18/94 16:19:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:04   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.7 Sources and sinks
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use comai
      use combi
      use comci
      use comco2
      use comdi
      use comei
      use comfi
      use comhi
      use comii
      use comdti
      use comflow
      use compart, only : pcnsk
      implicit none

      integer i, i1, i2, i3, i4, idiff, ii, iii, inode, inneq, ix, izone
      integer j, jj, jk, jx, k, kb, kb1, kb_max, addnode
      integer i_neigh, max_neigh, neqp1, nx,  nx1,  nx2, ny, ny1, ny2
      integer idummy, iconn, indexa_axy, ifl2, imat, mod_type
      integer open_file, ifile1, ifile2, ifile3, num_nodes, nawt
      integer counter, dumcyc, incf, ntimes, nfix
      integer :: readflag = 0
      integer, allocatable :: node_a(:)
      integer, allocatable :: node_b(:)
      integer, allocatable :: idum(:)
      integer, allocatable :: idum_a(:,:)
      integer, allocatable :: idum_b(:)
      integer, allocatable :: idum_c(:)
      integer, allocatable :: nel_node(:)
      integer, allocatable :: izone_awt(:)
c      integer, allocatable :: elem_temp(:,:)
      real*8, allocatable :: cons(:,:)
      real*8, allocatable :: cord_add(:,:)
      real*8, allocatable :: times(:)
      real*8 tdum,cstar,c1star,c2star,astar,bstar
      real*8 delx, delx1, delx2, dely, dely1, dely2
      real*8 dumx, dumy, dis, dis2, dis_min, discal, discal2
      real*8 x0, x01, x02, xa1, y0, y01, y02, ya1, xx, yy, zz 
      real*8 xbmin, xbmax, ybmin, ybmax, xcal, ycal, xmin, ymin
      real*8 x_level_1, x_level_2, y_level_1, y_level_2, y_split
      real*8 x_cutoff_E, x_cutoff_W, y_cutoff_N, y_cutoff_S
      real*8 xw1, xe, ys, yn, deldist, x_center, y_center, radius
      real*8 head_val, hgt_max, perm_fac, perm_tol, perm_set 
      real*8 recharge_new, recharge_orig
      real*8 sum_volj, sum_volj100, tempa, temp_grad
      real*8 total_recharge, total_recharge_in, total_recharge1
      real*8 total_recharge2, tol, fluxd, vol, frac, shut_time
      real*8 rzw, sww, scc
      real*8 t1, t0, avap1, avap2, ms, xf, and, dand , dela 
      real*8 pvwn,pvwn1,pvwn2,pvwn3, pvwd,pvwd1,pvwd2,pvwd3
      real*8 pvw,pvw0, bcoef1,bcoef2
      real*8 dla0, dlpa1, dlpa2, dlpa3, dlta1, dlta2, dlta3
      real*8 dlpta, dlpt2a, dlp2ta
      real*8 dlb0, dlpb1, dlpb2, dlpb3, dltb1, dltb2, dltb3
      real*8 dlptb, dlpt2b, dlp2tb
      real*8 tl,tl2,tl3,xm,xm2,xm3,tlxm,tl2xm,tlxm2
      integer iuser1,iuser2,icount
      integer i_outer_elem,ie_1,ie_2,ie_3,ie_4,ie_5,ie_6,ie_7,ie_8
      integer i_edge, nei_old, neq_add, i_new, i_elem_new, kk, ns_old
      integer i_level, i_level_max
      real*8 permxd,x_old,y_old,rad_fac, rad_max, rad_max2
      real*8 rad_old, rad_old2, rad_new, rad_new2, rad_use
      real*8 rad_fac2,rad_use2, x_max_new, y_max_new
      real*8, allocatable :: sk_save(:)
      real*8, allocatable :: sk_save1(:)
      real*8, allocatable :: sk_save2(:)
      real*8, allocatable :: temp1(:), temp2(:,:)
      real*8, allocatable :: ptemp1(:),ptemp2(:),ptemp12(:)
      real*8, allocatable :: water(:)
      real*8, allocatable :: x_wt(:),y_wt(:)
      character*80 dum_user
c gaz 041219 - 043019     
      integer nblk, izone_recharge, izone_num1, izone_num2, idir_zone
      integer ifile4, max_lines, izone_num0, num_rech_blks, kb_model
      integer, allocatable :: izone_node(:)
      real*8, allocatable :: area_surface(:)
      integer, allocatable :: node_surface(:)
      real*8 recharge1, recharge2, area_total, volume_total
      real*8 por_air, por_tbu, por_wlt, por_cddl
      real*8 area_tbu, area_wlt, area_cddl
      character*200 wdd1_t,file_zone1,file_area1,file_area2
      character*200 file_node1, file_node2, file_recharge_out
      character*8 macro
      logical approx
      
      character*80 water_file, above_wt_file, outside_zone_file
      logical matrix_node, null1, xy, ex

      save ntimes, times, cons, counter, dumcyc, nfix, readflag,
     &     shut_time, num_nodes, idum

      select case (k)
c gaz added for identifing edge nodels for extending the grid in radial direction
      case (2001)
          write(*,*) '**********************************'
          write(*,*) '  '
          write(*,*)'enter number of levels,rad_fac, rad_max, NOW'
          read(*,*) i_level_max, rad_fac, rad_max
          rad_fac2 = rad_fac*rad_fac
          rad_max2 = rad_max*rad_max
          rad_new2 = 0.0
c perform at least one level
          if(i_level_max.lt.1) then
           write(*,*)'setting number of levels to 1'
            i_level_max = 1
          endif          
          i_level_max = max(i_level_max,1)
          i_level = 0
2002      continue
          i_level = i_level + 1
          if(i_level.gt.i_level_max.or.rad_new2.ge.rad_max2) go to 2003
          if(i_level.eq.1) then
            incf = open_file('zone_299.txt','unknown')
          else
            incf = open_file('zone_edge.txt','unknown')
          endif
           read(incf,*) i1
           allocate(node_a(i1))   
           read(incf,*)(node_a(i),i = 1,i1)
           i_edge = i1
          close(incf)
          if(i_level.eq.1) then
           incf = open_file('Coarse_Grid_temp.grid','unknown')
          else
           incf = open_file('Coarse_Grid_new.grid','unknown')
          endif
           read(incf,*) dum_user (1:10)
           read(incf,*) ii
            allocate(cord_add(ii+i_edge,3))
            idum = 0
            neq_add = ii
           do i = 1,ii
            read(incf,*) jj,(cord_add(jj,k),k =1,3)
           enddo
           read(incf,*) dum_user (1:10)
           read(incf,*) jj, ii
           ns_old = jj
           allocate(idum_a(neq_add,ns_old))
c storage for elements need to be thought out  (&&&&&&&&&&&&&&&&&&&)         
           allocate(elem_temp(ii+ii,jj))
           allocate(idum_b(ii))
           allocate(nel_node(neq_add))
            nel_node = 0
            idum_b = 0
            nei_old = ii
           do i = 1,ii
            read(incf,*) kk,(elem_temp(kk,j),j = 1,jj)
c id nodes in this element set            
            do j = 1,jj
              i2 = elem_temp(kk,j)  
              nel_node(i2) = nel_node(i2)+1
              idum_a(i2,nel_node(i2)) = i
            enddo
           enddo
           close(incf)
c now find outside edge elements
          do i = 1, i_edge
            j = node_a(i)
           do i3 = 1,nel_node(j)
           i4 = idum_a(j,i3)
            if(i4.ne.0) then
             idum_b(i4) = 1
            endif
           enddo
          enddo
c count elements
         jk = 0
         do i = 1, ii
           if(idum_b(i).ne.0) then
            jk = jk + 1
           endif       
         enddo
         allocate (idum_c(jk))
         write(*,*) 'i_level = ',i_level,' outer edge nodes = ', i_edge
         write(*,*) 'i_level = ',i_level,' outer edge elements = ', jk
         jk = 0
         do i = 1, ii
           if(idum_b(i).ne.0) then
            jk = jk + 1
            idum_c(jk) = i
           endif       
         enddo  
         i_outer_elem = jk
          allocate(node_b(neq_add))
          i_new = neq_add
          x_max_new = 0.0
          y_max_new = 0.0
         do i = 1,i_edge
           i_new = i_new +1
           node_b(node_a(i)) = i_new
c assign coordinates
           x_old = cord_add(node_a(i),1)
           y_old = cord_add(node_a(i),2)  
           rad_old2 = x_old**2 + y_old**2
           rad_new2 = rad_fac2*rad_old2
           if(rad_new2.gt.rad_max2)then
            rad_use2 = rad_max2/rad_old2
            rad_use = sqrt(rad_use2)
            rad_new2 = rad_max2
           else
            rad_use = rad_fac
           endif
           cord_add(i_new,1)=rad_use*x_old
           cord_add(i_new,2)=rad_use*y_old
           cord_add(i_new,3)=cord_add(node_a(i),3)
           x_max_new = max(rad_use*x_old,x_max_new)
           y_max_new = max(rad_use*y_old,y_max_new)
         enddo
         write(*,*) 'i_level = ',i_level,' radius = ', sqrt(rad_new2)
         write(*,*) 'i_level = ',i_level,' xmax = ', x_max_new
         write(*,*) 'i_level = ',i_level,' ymax = ', y_max_new
         i_elem_new = nei_old
         do i = 1, i_outer_elem
          i1 = idum_c(i)
          ie_1 = elem_temp(i1,1)
          ie_2 = elem_temp(i1,2)
          ie_3 = elem_temp(i1,3)
          ie_4 = elem_temp(i1,4)
          ie_5 = elem_temp(i1,5)
          ie_6 = elem_temp(i1,6)
          ie_7 = elem_temp(i1,7)
          ie_8 = elem_temp(i1,8)
          i_elem_new = i_elem_new + 1
           elem_temp(i_elem_new,1) = ie_2
           elem_temp(i_elem_new,4) = ie_3
           elem_temp(i_elem_new,5) = ie_6
           elem_temp(i_elem_new,8) = ie_7
           elem_temp(i_elem_new,2) = node_b(ie_2)
           elem_temp(i_elem_new,3) = node_b(ie_3)
           elem_temp(i_elem_new,6) = node_b(ie_6)
           elem_temp(i_elem_new,7) = node_b(ie_7)
         enddo
         write(*,*)'New grid: nodes ',i_new,' elements ', i_elem_new
         incf = open_file('zone_edge.txt','unknown')
           write(incf,*) i_edge
           write(incf,*)(node_b(node_a(i)),i = 1,i_edge)
         close(incf)
         incf = open_file('Coarse_Grid_new.grid','unknown')
          write(incf,'(a4)') 'coor'
          write(incf,*) i_new
          do i = 1, i_new
           write(incf,'(i8,1x,3(f20.7))') i,(cord_add(i,jj),jj=1,3)
          enddo
          write(incf,'(a4)') '    '
          write(incf,'(a4)') 'elem'
          write(incf,*) ns_old,  i_elem_new, 0
          do i = 1, i_elem_new
           write(incf,'(i8,8(1x,i8))') i,(elem_temp(i,jj),jj=1,ns_old)  
          enddo
          write(incf,'(a4)') '    '
          write(incf,'(a4)') 'stop'
         close(incf)
         deallocate(idum_a,idum_b,idum_c,
     &          node_a,node_b,elem_temp,cord_add,nel_node)
         go to 2002
2003     continue  
         write(*,*) 'Generated new grid: Coarse_Grid_new.grid '
         write(*,*) 'nodes = ',i_new,' elements = ',i_elem_new
         write(*,*) '  '
         write(*,*) '**********************************'    
          stop

c gaz added for turning off CO2 identified CO2 injection wells injection wells
      case (1001)
         if(readflag.eq.0) then
            incf = open_file('co2_inj.txt','unknown')
            read(incf,*) shut_time
            read(incf,*)  num_nodes
            allocate(idum(num_nodes))
            read(incf,*)  (idum(i),i= 1,num_nodes)
            close (incf)
            readflag = 1
         end if	
         if(days.ge.shut_time*365.25.and.readflag.eq.1) then
          readflag = 2
          write(*,*) '*************WEllS SHUT IN***********************'
          write(iout,*) '*************WEllS SHUT IN********************'
            do j = 1, num_nodes
              i = idum(j)
               if(kaco2(i).le.-1) then
                  wellco2(i) = 0.d0
               elseif(skco2(i).lt.0.d0) then
                  skco2(i) = 0
               endif
           enddo
         endif          
c RJP added for turning off CO2 injection wells
      case (999)
         if(readflag.NE.1) then
            incf = open_file('co2_inj.txt','unknown')
            read(incf,*) shut_time
            close (incf)
            readflag = 1
         end if	
         if(days.ge.shut_time*365.25.and.readflag.eq.1) then
          readflag = 2
          write(*,*) '*************WEllS SHUT IN***********************'
          write(iout,*) '*************WEllS SHUT IN********************'
            do i = 1, neq
               if(kaco2(i).le.-1) then
                  wellco2(i) = 0.d0
               elseif(kaco2(i).eq.-4) then
                  wellco2(i) = 0.d0
               elseif(skco2(i).lt.0.d0) then
                  skco2(i) = 0
               endif
           enddo
         endif
c RJP added 10/17/07 for changing the saturations for goldsim
      case(998)
         if(readflag.NE.1) then
            incf = open_file('time_dep_sat.txt','unknown')
            read(incf,*) rzw,sww,scc
            close(incf)
            readflag = 1
            fow(1) = sww
            fol(1) = scc
            fw(1) = sww
            fl(1) = scc
            fog(1) = 1.d0-fow(1)-fol(1)
            fg(1) = 1.d0-fw(1)-fl(1)
            if(rzw.ne.0) then
               fow(2) = sww
               fol(2) = scc
               fw(2) = sww
               fl(2) = scc
               fog(2) = 1.d0-fow(2)-fol(2)
               fg(2) = 1.d0-fw(2)-fl(2)
            endif
         endif
c      PHS  10/05/05   Adding user option 666 to allow D. Neeper to 
c      test the analytical solution
c      PHS   4/24/08   Updating to allow nfix nodes to be fixed.
c-----------------------------------------------------------------------
      case (666)
         if(readflag.NE.1) then
            incf = open_file('nodes_conc_fixed.macro','unknown')
            read(incf,*) ntimes, nfix
            allocate(times(ntimes),cons(ntimes,nfix),node_a(nfix))
            read(incf,*) (node_a(iii), iii = 1,nfix)
            do iii = 1,ntimes
               read(incf,*) times(iii), (cons(iii,j), j = 1,nfix)
            end do
            close (incf)
            counter = 1
            dumcyc = 1
            readflag = 1
            do iii = 1,nfix
               pcnsk(node_a(iii)) = -1
            end do
         end if
         tdum = days - (dumcyc -1)*times(ntimes)
         if(tdum.GE.times(counter+1)) counter = counter + 1
         do iii = 1,nfix
            cnsk(node_a(iii)) = -cons(counter,iii)
         end do
         if(counter.GE.ntimes) then
            counter = 1
            dumcyc = dumcyc + 1
         end if
      case (-291)
         if(l.eq.1) then
c     
c     edit wt data for boundaries                  
c     
!           open(unit=98,file='water_table_on_boun',status='unknown')
!           open(unit=99,file='wt.dat',status='unknown')
            write(*,*) 'user sub (-291), ',
     &           'output water table on boundaries'
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile1 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile1 = open_file(dum_user,'old')
            end if
            read(ifile1,'(a80)') water_file
            read (ifile1,'(a80)') dum_user
!            write(*,*) 'Enter coordinates for model boundaries, ',
!     &           'xw1,xe,ys,yn'
            read(ifile1,*) xw1,xe,ys,yn
!            write(*,*) 'Enter distance criteria for boundary search'
            read(ifile1,*) deldist
            close (ifile1)
            ifile1 = open_file(dum_user,'unknown')
            ifile2 = open_file(water_file,'old')
            write(*,*) 'BCs from ', trim(water_file)
            write(ifile1,*) 'BCs from ', trim(water_file)
            read(ifile2,*)
            read(ifile2,*)
            read(ifile2,*)
            read(ifile2,*)
            i1 = 0
            do i=1,10000000
c     top boundary side node with new coords reflecting new model
               read(ifile2,*, end = 987) x0,y0, dumx
               if(abs(x0-xw1).lt.deldist) then 
                  i1 = i1 + 1
                  write(ifile1,'(3g18.9)') x0,y0,dumx
               else if(abs(x0-xe).lt.deldist) then
                  i1 = i1 + 1
                  write(ifile1,'(3g18.9)') x0,y0,dumx            
               else if(abs(y0-ys).lt.deldist) then
                  i1 = i1 + 1
                  write(ifile1,'(3g18.9)') x0,y0,dumx
               else if(abs(y0-yn).lt.deldist) then
                  i1 = i1 + 1
                  write(ifile1,'(3g18.9)') x0,y0,dumx
               endif
            enddo
 987        continue
            write(*,*) 'user sub completed, node found = ', i1
            write(*,*) 'file "', trim(dum_user), '" created '
            close(ifile1)
            close(ifile2)
            stop
         endif
      case (-292)
         if(l.eq.1) then
c     
c     print out new flow boundaries for new model                 
c     
!           open(unit=99,file='boun_head_base_case',status='unknown')
!           open(unit=98,file='boun_head_new.macro',status='unknown')
!           open(unit=97,file='supicious_nodes_boun',status='unknown')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile3 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile3 = open_file(dum_user,'old')
            end if
! Name of macro file to be created
            read (ifile3, '(a80)') water_file
            ifile1 = open_file(water_file,'unknown')
! Name of water table boundary file (created with -291)
            read (ifile3, '(a80)') dum_user
            ifile2 = open_file(dum_user,'old')
! Name of file to write questionable boundary nodes to
            read (ifile3, '(a80)') dum_user
! Distance criteria for screening suspicious nodes
            read (ifile3,*) dis
            close (ifile3)
            ifile3 = open_file(dum_user,'unknown')
!            dis = 250.
            write(ifile3,*) 
     &           'Suspicious nodes with distances greater than ', dis
            read(ifile2,*)
            do i2 = 1,neq
               read(ifile2,*,end = 299) pnx(i2),pny(i2),pnz(i2)
            enddo
 299        continue
            i2 = i2 -1
            i1 = 0
            ii = 0
            write(*,*) 'user sub (-292)'
            write(ifile1,'(a4)') 'flow'
            do i=1,neq
               if(ka(i).eq.-1) then
c     top boundary side node with new coords reflecting new model

                  x0 = cord(i,1)
                  y0 = cord(i,2)
                  dumx = 100000000.
                  do j = 1,i2
                     dumy = (x0-pnx(j))**2+(y0-pny(j))**2
                     if(dumy.lt.dumx) then
                        i3 = j
                        dumx = dumy
                     endif
                  enddo
                  if(dumx.gt.dis*dis) then
                     ii = ii +1
                     write(ifile3,876) i3,i,x0,y0,cord(i,3),sqrt(dumx)  
                  endif
                  i1 = i1 +1
                  call headctr(4,i,pflow(i),head_val)
                  
                  write(ifile1,877) i,i,1, pnz(i3), 1., 1.e06,
     &                 head_val-pnz(i3),sqrt(dumx) 
               endif
            enddo
 876        format('base',i6,' new node ',i8,' x ',g15.7,' y ',g13.7,
     &           ' z ',g13.7,' diff ',g13.5)      
 877        format(3(1x,i7),3(1x,g14.7),' # diff ',2g12.4)
            write(ifile1,'(a4)') '    '
            write(ifile1,'(a4)') 'stop'
            write(*,*) 'user sub completed, nodes written = ', i1
            write(*,*) 'file "', trim(water_file),'" created '
            write(*,*) 'number of suspicious nodes(dis > ',dis,' ) = ', 
     &           ii
            write(*,*) 'nodes are in file: ', trim(dum_user)
            close(ifile1)
            close(ifile2)
            close(ifile3)
            stop
         endif
      case (-303)
         if(l.eq.1) then
c     
c     read recharge map from bill arnold              
c     write flow macro for fehm
c     
            write(*,*) 
     &           '>>> User sub for recharge(-303) invoked <<< '
            write(*,*) '   '
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile2 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile2 = open_file(dum_user,'old')
            end if
! Recharge data file
            read(ifile2,'(a80)') dum_user
            ifile1 = open_file(dum_user,'old')
! Output macro name
            read(ifile2,'(a80)') water_file
! Water table zone file
            read(ifile2,'(a80)') dum_user
            read(ifile2,*) x0,y0,delx,dely,nx,ny
            close (ifile2)
            allocate(sk_save(n0))
            allocate(idum(n0))
            sk_save=0.0d00
            total_recharge = 0.0d00
            do
               read(ifile1,*, end=555, err=555) ii, i1, i2, sk_save(ii)
               total_recharge = total_recharge + sk_save(ii)
            enddo
 555        continue
            close(ifile1)
            total_recharge_in = total_recharge 
!           open(unit=98,file='flow_recharge.macro',status='unknown')
            ifile2 = open_file(water_file,'unknown')
            write(ifile2,'(a20)') 'flow      # recharge'
!            write(*,*) 'Enter name of water table zone file'
!            read(*,'(a80)') dum_user      
!           open(unit=99,file=dum_user,status='old')
            ifile3 = open_file(dum_user,'old')
            read(ifile3,'(a4)') dum_user(1:4)
            read(ifile3,'(a4)') dum_user(1:4)
            read(ifile3,'(a4)') dum_user(1:4)
            read(ifile3,*) num_nodes,(irb(i),i=1,num_nodes)      
            close(ifile3)
!           open(unit=96,file='recharge.plt',status='unknown')
            dum_user = trim(water_file) // ".plt"
            ifile1 = open_file(dum_user,'unknown')
            write(ifile1,*) 
     &           'recharge file for tecplot-delete before use'
            write(ifile1,557)
 557        format(13x,'x',13x,'y',10x,'z',3x,'recharge',
     &           3x,'node in',1x,'node fehm')
            total_recharge = 0.0d00
            iii = nx*ny
            do ii=1,num_nodes
               kb = irb(ii)                   
c     check for unconnected nodes
               if(nelm(kb)+1.ne.nelm(kb+1)) then
                  xx = cord(kb,1)
                  yy = cord(kb,2)
                  ix = nint((xx-x0)/delx) +1
                  jx = nint((yy-y0)/dely) +1
                  i = (jx-1)*nx +ix
                  zz = cord(kb,3)
                  if (i.gt.0.and.i.le.iii) then
                     write(ifile2,556) kb,kb,1,sk_save(i),1.0d00,0.0d00
                     write(ifile1,'(2(1x,e13.6),2(1x,f10.4),2(2x,i8))') 
     &                    xx,yy,zz,sk_save(i),i,kb
                     total_recharge = total_recharge +sk_save(i)
                     sk_save(i) = 0.0d00
                  endif
               endif
            enddo
 556        format(3(1x,i8),3(1x,e12.4))
            write(ifile2,'(a4)') '    '
            write(ifile2,'(a4)') 'stop'
            write(*,558) total_recharge_in, total_recharge
 558        format('total_recharge(input) = ',f12.4,' output = ',f12.4) 
            write(*,*) 'created new file : ', trim(water_file)
            close(ifile1)
            close(ifile2)
            stop
         endif
      case (-304)
         if(l.eq.1) then
c     
c     modify recharge to remove low infiltration flux
c     
            write(*,*) 'User sub(-304) invoked '
            write(*,*) 'Procedure to modify distribution of recharge'
            write(*,*) 'Need top nodes in file : top.zone '
            write(*,*) 'top.zone should not contain side BC nodes '
            write(*,*) 'or the pure top nodes can be labeled with '
            write(*,*)'porosity < 0.5 with side nodes porosity > 0.5'
            write(*,*) '   '
!           open(unit=95,file='top.zone',status='old')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile2 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile2 = open_file(dum_user,'old')
            end if
! Top zone file
            read(ifile2,'(a80)') dum_user
            ifile1 = open_file(dum_user,'old')
! Permeability cutoff macro
            read(ifile2,'(a80)') above_wt_file
! Recharge plot file
            read(ifile2,'(a80)') water_file
! Permeability cutoff (suggested 1.d-17)
            read(ifile2,*) perm_tol
! Level of neighbor search (suggested 2) 
            read(ifile2,*) max_neigh 
            close (ifile2)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,*) num_nodes,(irb(i),i=1,num_nodes)
            close(ifile1)
c     make sure water table nodes are tagged
c     sum recharge flux                       
            allocate(sk_save(n0))
            recharge_orig = 0.0d00
            iirb = 0
            sk_save = 0.0d00
            do i=1,num_nodes
               kb = irb(i)
               if(ps(kb).lt.0.5d00) then
                  iirb(kb) = 1
                  recharge_orig = recharge_orig +sk(kb)
                  sk_save(kb) = sk(kb)
               endif
            enddo
!           open(unit=99,file='recharge_tec.plt',status='unknown')
!  open(unit=97,file='recharge_perm_cutoff_new.macro',status='unknown')
            ifile1 = open_file(water_file,'unknown')
            ifile2 = open_file(above_wt_file,'unknown')
            write(ifile2,'(a4)') 'flow'
            write(*,*) 'removing infitration for low perms'
!            write(*,*) 'Enter permeability cutoff(suggest 1.d-17)'
!            read(*,*) perm_tol     
!            write(*,*) 'Enter level of neighbor search(suggest 2)'
!            read(*,*) max_neigh    
            perm_tol=1.d+6*perm_tol
            tol = 1.0d-5
            neqp1 = neq+1
            do i_neigh = 1,max_neigh
               do i=1,neq
                  if(iirb(i).ne.0.and.sk(i).ne.0.0d00) then  
                     if(pnx(i).lt.perm_tol) then    
                        i1 = nelm(i)+1
                        i2 = nelm(i+1)
                        hgt_max = 0.0d00   
                        kb_max = 0
                        do jj=i1,i2
                           kb = nelm(jj)
                           if(iirb(kb).ne.0.and.kb.ne.i) then
                              if(cord(i,3)-cord(kb,3).gt.hgt_max+tol)
     &                             then
                                 kb_max=kb
                                 hgt_max = cord(i,3)-cord(kb,3)
                              endif
                           endif
                        enddo
                        if(kb_max.eq.0) then
                           do jj=i1,i2
                              kb = nelm(jj)
                              if(iirb(kb).ne.0.and.kb.ne.i) then
                                 if(cord(i,3)-cord(kb,3).le.
     &                                hgt_max+tol .and. cord(i,3)-
     &                                cord(kb,3).gt.hgt_max-tol) then
                                    if(pnx(kb).gt.pnx(kb_max)) then
                                       kb_max=kb
                                    endif
                                 endif
                              endif
                           enddo
                        endif
                        if(kb_max.ne.0.and.kb_max.ne.i) then
                           sk(kb_max) = sk(kb_max) + sk(i)
                           sk(i) = 0.0d00
                        else
                           do jj=i1,i2
                              kb = nelm(jj)
                              dis_min = 1.d+30
                              if(iirb(kb).ne.0.and.kb.ne.i) then
                                 dis = (cord(i,1)-cord(kb,1))**2 +
     &                                (cord(i,2)-cord(kb,2))**2 
                                 if(dis.lt.dis_min) then
                                    kb_max = kb
                                 endif
                              endif
                           enddo
                           if(kb_max.ne.i) then
                              sk(kb_max) = sk(kb_max) + sk(i)
                              sk(i) = 0.0d00
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
            recharge_new = 0.0d00
            write(ifile1,681)
 681        format(12x,'x',12x,'y',5x,'flux old',5x,'flux new',
     &           4x,'flux diff',4x,'node')
            idiff = 0 
            do i= 1,neq
               if(iirb(i).ne.0) then
                  write(ifile2,'(3i7,3g12.4)') 
     &                 i,i,1,sk(i),1.0d00,0.0d00
                  if(abs(sk_save(i)-sk(i)).gt.1.d-20) idiff = idiff+1
                  recharge_new = recharge_new +sk(i)
                  write(ifile1,'(5(1x,g12.4),1x,i7)') cord(i,1),
     &                 cord(i,2),sk_save(i),sk(i),sk_save(i)-sk(i),i
               endif
            enddo
            write(ifile2,'(a4)') '    '
            write(ifile2,'(a4)') 'stop'
            write(*,682) recharge_orig,recharge_new, idiff
 682        format('original recharge : ',f9.4,
     &           ' new recharge : ',f9.4,' sources modified ',i6)
            write(*,*) 'file written: ', trim(above_wt_file)
            write(*,*) 'file written: ', trim(water_file) 
            close(ifile2)
            close(ifile1)
            stop
         endif
      case (-305)
         if(l.eq.1) then
c     
c     modify permeability if necessary in recharge zones   
c     
            write(*,*) 'User sub(-305) invoked '
            write(*,*) 
     &           'Procedure to modify permeability in recharge nodes'
            write(*,*) 'Only changes y-z permeabilities'
            write(*,*) 'Need top nodes in file : top.zone'
            write(*,*) 'top.zone should not contain side BC nodes'
            write(*,*) 'or the pure top nodes can be labeled with'
            write(*,*) 
     &           'porosity < 0.5 with side nodes porosity > 0.5'
            write(*,*) '   '
!           open(unit=95,file='top.zone',status='old')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile2 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile2 = open_file(dum_user,'old')
            end if
! Top zone file
            read(ifile2,'(a80)') dum_user
            ifile1 = open_file(dum_user,'old')
! Permeability cutoff macro
            read(ifile2,'(a80)') water_file
! Permeability cutoff (suggested 1.d-17)
            read(ifile2,*) perm_tol     
! Level of neighbor search(0, 1, 2, 3)
            read(ifile2,*) max_neigh 
! Bounding box for changes (xmin, xmax, ymin, ymax)  
            read(ifile2,*) x_level_1, x_level_2, y_level_1, y_level_2
            close (ifile2)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,'(a4)') dum_user(1:4)
            read(ifile1,*) num_nodes,(irb(i),i=1,num_nodes)  
            close (ifile1)
c     make sure water table nodes are tagged
c     sum recharge flux                       
            allocate(sk_save(n0))
            recharge_orig = 0.0d00
            iirb = 0
            sk_save = 0.0d00
            do i=1,num_nodes
               kb = irb(i)
               if(ps(kb).lt.0.5d00) then
                  iirb(kb) = 1
                  recharge_orig = recharge_orig +sk(kb)
                  sk_save(kb) = sk(kb)
               endif
            enddo
            write(*,*)'total recharge top.zone = ',recharge_orig
            do i = 1,n
               pnx(i)=pnx(i)*1.d-6
               pny(i)=pny(i)*1.d-6
               pnz(i)=pnz(i)*1.d-6
            enddo
            write(*,*) 'Creating new file: ', trim(water_file)
!           open(unit=97,file='perm_cutoff_recharge.macro',status='new')
            ifile1 = open_file(water_file,'unknown')
            write(ifile1,'(a4)') 'perm'
            idiff = 0
            dum_user(1:1) = 'n'
 977        continue
            write(*,*) 'modifying perms in recharge nodes'
            if(dum_user(1:1).eq.'y'.or.dum_user(1:1).eq.'Y') then
               write(*,*) 'Enter permeability cutoff(suggest 1.d-17)',
     &              'previous value ', perm_tol
               read(*,*) perm_tol     
!            write(*,*) 'Enter level of neighbor search(0, 1, 2, 3)'
!            read(*,*) max_neigh    
!            write(*,*) 'Enter y_level_1: only changes to the North'
!            read(*,*) y_level_1      
!            write(*,*) 'Enter y_level_2,x_level_1: ',
!     &           'changes North_and East'
!            read(*,*) y_level_2,x_level_1
            end if      
            if(max_neigh.eq.0) then
               write(*,*) 'set all perms in Northern region'
            else if(max_neigh.lt.1) then
               write(*,*) 'neighbor search entered < 1,set =1'
               max_neigh=1
            else if(max_neigh.gt.3) then
               write(*,*) 'neighbor search entered > 3,set =3'
               max_neigh=3
            endif
            tol = 1.0d00
            irb = 0
            neqp1 = neq+1
            do i=1,neq
               if(iirb(i).ne.0.and.sk(i).lt.-1.d-12) then  
                  if(cord(i,1).ge.x_level_1 .and. cord(i,1).le.x_level_2
     &                 .and. cord(i,2).ge.y_level_1 .and.
     &                 cord(i,2).le.y_level_2) then  
                     if(pnx(i).lt.perm_tol) then    
                        pnx(i)= perm_tol   
                        pny(i)= perm_tol   
                        pnz(i)= perm_tol   
                        if(irb(i).eq.0) then
                           idiff = idiff +1
                           irb(i) = 1
                        endif
                        if(max_neigh.gt.1) then
                           i1 = nelm(i)+1
                           i2 = nelm(i+1)
                           do jj=i1,i2
                              kb = nelm(jj)
                              if(iirb(kb).ne.0) then
                                 perm_set = max(perm_tol,pnx(kb),
     &                                pny(kb),pnz(kb))
                                 pnx(kb)= perm_set
                                 pny(kb)= perm_set
                                 pnz(kb)= perm_set
                                 if(irb(kb).eq.0) then
                                    idiff = idiff +1
                                    irb(kb) = 1
                                 endif
                                 if(max_neigh.gt.2) then
                                    i3 = nelm(kb)+1
                                    i4 = nelm(kb+1)
                                    do jk=i3,i4
                                       kb1 = nelm(jk)
                                       if(iirb(kb1).ne.0) then
                                          perm_set = max(perm_tol,
     &                                         pnx(kb1),pny(kb1),
     &                                         pnz(kb1))
                                          pnx(kb1)= perm_set
                                          pny(kb1)= perm_set
                                          pnz(kb1)= perm_set
                                          if(irb(kb1).eq.0) then
                                             idiff = idiff +1
                                             irb(kb1) = 1
                                          endif
                                       endif
                                    enddo
                                 endif
                              endif
                           enddo
                        endif
                     endif
                  endif
               endif
            enddo
            write(*,*) '>>>> number of permeability changes ', idiff 
            write(*,*) 'Try again with new perms (y/n) ?'
            read(*,'(a4)') dum_user
            if(dum_user(1:1).eq.'y'.or.dum_user(1:1).eq.'Y') then
               go to 977
            endif
            do i= 1,neq
               if(irb(i).ne.0) then
                  write(ifile1,'(3(1x,i7),3(1x,1p,g12.4))') 
     &                 i,i,1,pnx(i),pny(i),pnz(i)
               endif
            enddo
            write(ifile1,'(a4)') '    '
            write(ifile1,'(a4)') 'stop'
            write(*,*) ' file written: ', trim(water_file) 
            close(ifile1)
            stop
         endif
      case (-401)

         if(l.eq.1) then
c     
c     set flux zones                       
c     
            write(*,*)
     &           ' >>>> User sub for extracting fluxes <<<<  '
            write(*,*)
     &           ' >>>> Uses zones defined in flxz macro <<<< '
            write(*,*)
     &           ' >>>> Also has material number and volume <<<< '
            write(*,*)
     &           ' >>>> enter file name for mapping fluxes <<<< '
            read(*,'(a80)') dum_user 
!           open(unit=97,file=dum_user,status='unknown')
            ifile1 = open_file(dum_user,'unknown')
            write(*,*)
     &           ' >>>> enter 2nd file name for mapping fluxes <<<< '
            write(*,*)
     &           ' >>>> return means no mapping of fluxes <<<< '
            read(*,'(a80)') dum_user 
            write(*,*) 'running model ...'
         endif
         if(ifinsh.ne.0.or.l.eq.nstep-1) then
c     
c     print out flow boundaries        
c     
            
            allocate(sk_save1(n0),sk_save2(n0))
            allocate(idum(n0),idum_b(n0))
            allocate(water(n0))
            allocate(x_wt(n0),y_wt(n0))

            sk_save1 = 0.0
            sk_save2 = 0.0
            water = 0.0
            do izone = 1, nflxz
               idum = 0
c     Loop over all nodes
               do inode = 1, n0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
                  if(inode.gt.neq) then
                     inneq = inode-neq
                     matrix_node = .true.
                     addnode = nelm(neq+1)-neq-1
                     idummy = neq
                  else
                     matrix_node = .false.
                     inneq = inode
                     addnode = 0
                     idummy = 0
                  end if
c     Determine if node is part of the zone being summed
                  if(izoneflxz(inode).eq.izone) then
c     identify node belonging to this zone
                     idum(inode) = izone
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
                     i1 = nelm(inneq)+1
                     i2 = nelm(inneq+1)
c     loop over every connecting node
                     do iconn = i1, i2
                        indexa_axy = iconn-neq-1+addnode
                        if(a_axy(indexa_axy).gt.0.) then
                           if(izoneflxz(idummy+nelm(iconn))
     2                          .ne.izone.or.nelm(iconn)
     3                          .eq.inneq) then
                              sk_save1(inode) = 
     4                             sk_save1(inode) + a_axy(indexa_axy)
                           end if
                        end if
                     end do
                  end if
               end do
c     Write results
               do jj = 1,n0
                  if(idum(jj).eq.izone) then
                     write(ifile1,1045) jj,iflxz(izone),
     2                    sk_save1(jj),irlp(jj),sx1(jj)
                     water(irlp(jj)) = water(irlp(jj)) +
     &                    sk_save1(jj) 
                  endif
               enddo
 1045          format('node ',i7,' zone ',i5,
     &              ' flux ',g12.5,1x,'mat ',i5,' volume ',g12.5)
            enddo

            close(ifile1)
            if(dum_user(1:10).ne.'          ') then
!              open(unit=98,file=dum_user,status='unknown')
               ifile1 = open_file(dum_user,'old')
               x_wt = 0.0
               y_wt = 0.0
               idum_b = 0
               ifl2 =1
 1040          read(ifile1,*,end=935) dum_user(1:5),idum_b(ifl2),
     &              dum_user(1:6),izone,dum_user(1:6),fluxd,
     &              dum_user(1:4),imat,dum_user(1:7),vol
               idum(ifl2) = imat
               y_wt(ifl2) = vol
               x_wt(imat) = x_wt(imat) +vol
               ifl2 = ifl2 +1
               go to 1040
 935           continue
               ifl2 = ifl2 - 1
               close(ifile1)
 1047          format(a5,i7,a30,a4,i5,a8,f12.0)
!              open(unit=99,file='flux_boun.macro',status='unknown')
               ifile1 = open_file('flux_boun.macro','unknown')
               call headctr(2,0,0.0d00,0.0d00)        
               write(*,*)'>>>> Printing new flow macro: flux_boun.macro'
               write(ifile1,1048) 
 1048          format('flow   [user sub = -401]')
               do i=1,ifl2
                  frac = y_wt(i)/x_wt(idum(i))
                  ii = idum_b(i)
                  write(ifile1,1046) ii,ii,1,-frac*water(idum(i)),1.,0.,
     &                 idum(i)
               enddo
 1046          format(3(1x,i6),1x,g12.4,2(1x,f6.2),i6)
               write(ifile1,'(a4)') '    '
               write(ifile1,'(a4)') 'stop'
               close(ifile1)
            endif
            stop
         endif
      case (-402)
         if(l.eq.1) then
c     
c     print out volume fractions and coordinates(split zones)
c     print new zones 
c     
            write(*,*) 'k=-402, calculate new zones in LHG area'
            write(*,*) 'HFM units must be each given unique rlp model'
!           open(unit=99,file='volume_fractions.txt',status='unknown')
!           open(unit=89,file='zones_north_split.macro',status='unknown')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile1 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile1 = open_file(dum_user,'old')
            end if
! Volume fraction file
            read(ifile1,'(a80)') dum_user
! Split zones files
            read(ifile1,'(a80)') water_file
            ifile2 = open_file(water_file,'unknown')
!           write(*,*) ' ENTER model type : 1 = min y value' 
!           write(*,*) '                    2 = dis from caldera'
!           write(*,*) '                    3 = xy box'
            read(ifile1,*) mod_type
            if(mod_type.eq.1) then
!               write(*,*) ' ENTER y value(UTM) for split(N-S)'
               read(ifile1,*) y_split
               write(ifile2,'(a26,f12.3)') 'zonn    # northing value: ',
     &              y_split
            else if(mod_type.eq.2) then
!               write(*,*) ' ENTER r value and x,y for caldera center'
               read(ifile1,*) discal,xcal,ycal
               discal2 = discal*discal
               write(ifile2,'(a26,3f12.3)') 'zonn   # dis from cadera: '
     &              , discal , xcal, ycal
            else if(mod_type.eq.3) then
!               write(*,*) 
!     &              ' ENTER coordinates for box boundaries xw,xe,ys,yn'
               read(ifile1,*) xw1,xe,ys,yn
               write(ifile2,'(a26,4f12.3)') 'zonn   # xw, xe, ys, yn: ',
     &             xw1, xe, ys, yn 
            endif
            close (ifile1)
            write(*,*) ' calculating volume percentages....'
            ifile1 = open_file(dum_user,'unknown')
             do  j=1,nrlp
               sum_volj = 0.0d00
               sum_volj100 = 0.0d00
               ii = 0
               iirb = 0
               if(mod_type.eq.1) then
                  do i=1,neq
                     if(irlp(i).eq.j) then    
                        if(cord(i,2).gt.y_split) then
                           ii = ii+1
                           iirb(ii) = i
                           sum_volj100 = sx1(i)+sum_volj100
                        else
                           sum_volj = sx1(i)+sum_volj
                        endif
                     endif
                  enddo
               else if(mod_type.eq.2) then
                  do i=1,neq
                     dis2 = (cord(i,1)-xcal)**2 + (cord(i,2)-ycal)**2
                     if(irlp(i).eq.j) then    
                        if(dis2.le.discal2) then
                           ii = ii+1
                           iirb(ii) = i
                           sum_volj100 = sx1(i)+sum_volj100
                        else
                           sum_volj = sx1(i)+sum_volj
                        endif
                     endif
                  enddo
               else if(mod_type.eq.3) then
                  do i=1,neq
                     if(irlp(i).eq.j) then 
                        if (cord(i,1) .ge. xw1 .and. cord(i,1) .le. xe
     &                       .and. cord(i,2) .ge. ys .and. cord(i,2)
     &                       .le. yn) then
                           ii = ii+1
                           iirb(ii) = i
                           sum_volj100 = sx1(i)+sum_volj100
                        else
                           sum_volj = sx1(i)+sum_volj
                        endif
                     endif
                  enddo                  
               endif
               if(ii.ge.1) then
                  write(ifile2,'(i6)') j+100
                  write(ifile2,'(a4)') 'nnum'
                  write(ifile2,'(i6)') ii
                  write(ifile2,'(8i10)') (iirb(jj),jj=1,ii)
               endif
               write(*,*) 'unit ',j,' volume percent  ',
     &              sum_volj/vtot*100.0d00
               write(*,*) 'unit ',j+100,' volume percent  ',
     &              sum_volj100/vtot*100.0d00
               write(ifile1,*) 'unit ',j,' volume percent  ',
     &              sum_volj/vtot*100.0d00
               write(ifile1,*) 'unit ',j+100,' volume percent  ',
     &              sum_volj100/vtot*100.0d00
            enddo
            write(ifile2,*) ' '
            write(ifile2,'(a4)') 'stop'
            xmin=1.d30
            ymin=1.d30
            do i=1,neq
               xmin = min(xmin,cord(i,1))
               ymin = min(ymin,cord(i,2))
            enddo
            write(ifile1,*) ' '
            write(ifile1,*) 'x coordinates for 1000m spacing'
            write(ifile1,'(5(1x,f12.0))') (xmin+1000.*(i-1),i=1,31)
            write(ifile1,*) 'y coordinates for 1000m spacing'
            write(ifile1,'(5(1x,f12.0))') (ymin+1000.*(i-1),i=1,46)
            write(ifile1,*) 'x coordinates for 500m spacing'
            write(ifile1,'(5(1x,f12.0))') (xmin+500.*(i-1),i=1,61)
            write(ifile1,*) 'y coordinates for 500m spacing'
            write(ifile1,'(5(1x,f12.0))') (ymin+500.*(i-1),i=1,91)
            close(ifile1)
            close(ifile2)
            stop
         endif
      case (-408)
c
c create temperature profile for isothermal model
c
         if(l.eq.1) then
            write(*,*) ' >>>> User sub for setting temperature ',
     &           'invoked(-408) <<<< '
            write(*,*)
     &           ' >>>> Procedure to set initial conditions  <<<< '
            allocate(izone_awt(n0))
            allocate(temp1(n0))
            write(*,*) '**************************'
            write(*,*) 
            write(*,*) 'Temperature control file'
            write(*,*) 'Contains:'
            write(*,*) 'Land surface zone file'
            write(*,*) 'Ambient Temperature (C),Temperature gradient ',
     &           '(deg C/km) '
            write(*,*) 'delta x, delta y'
            write(*,*) 'Initial head'
            write(*,*) '"pres" macro file name '
!            write(*,*) 'Enter Temperature control file name NOW'
!            read(*,'(a80)') water_file
!           open(89,file=water_file,status='unknown')
!            ifile1 = open_file(water_file,'old')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile1 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile1 = open_file(dum_user,'old')
            end if
            read(ifile1,'(a80)') above_wt_file
            read(ifile1,*) tempa,temp_grad             
            read(ifile1,*) delx,dely
            read(ifile1,*) head_val
c     convert C/km/ to C/m
            temp_grad=temp_grad/1000.
            read(ifile1,'(a80)') outside_zone_file
            close(ifile1)
!           open(90,file=above_wt_file,status='unknown')
            ifile1 = open_file(above_wt_file,'old')
            read(ifile1,'(a80)') wdd
            read(ifile1,'(a80)') wdd
            read(ifile1,'(a80)') wdd
            read(ifile1,*)  nawt
            read(ifile1,*) (izone_awt(i),i=1,nawt)
            close(ifile1)
!           open(91,file=outside_zone_file,status='unknown')
            ifile1 = open_file(outside_zone_file,'unknown')
            write(ifile1,'(a41)') 
     &           'pres           #with temperature gradient'
            write(*,*) 'setting temperatures.....'
            x0 = 1.d20
            xa1 = -1.d20
            y0 = 1.d20
            ya1 = -1.d20
            do i=1,n0
               x0=min(x0,cord(i,1))
               xa1=max(xa1,cord(i,1))
               y0=min(y0,cord(i,2))
               ya1=max(ya1,cord(i,2))
            enddo
            ix = int((xa1 - x0)/delx) + 1
            jx = int((ya1 - y0)/dely) + 1
            allocate(temp2(ix,jx))
            temp2=0.
            do k=1,nawt
               kb=izone_awt(k)
               xa1=cord(kb,1)
               ya1=cord(kb,2)
               ix= int((xa1-x0)/delx) + 1
               jx= int((ya1-y0)/dely) + 1
               temp2(ix,jx) = cord(kb,3)
            enddo
            
            do kb=1,n0
               xa1=cord(kb,1)
               ya1=cord(kb,2)
               ix= int((xa1-x0)/delx) + 1
               jx= int((ya1-y0)/dely) + 1
               hgt_max=temp2(ix,jx)-cord(kb,3)
               temp1(kb) = tempa+temp_grad*hgt_max
            enddo

!           open(92,file='temps_not_set.err',status='unknown')
            ifile2 = open_file('temps_not_set.err','unknown')
            write(ifile2,'(a40)') 'node, x,y,z for node with no ',
     &           'temperature'
            iii=0
            do i=1,n0
               if(temp1(i).eq.0.0) then
                  iii=iii+1
                  write(ifile2,'(i9,1x,1p,3g16.7)') i,
     &                 (cord(i,jk),jk=1,3)
               else
                  write(ifile1,'(3i9,1x,1p,2g16.7,1x,i9)') 
     &                 i,i,1,head_val,temp1(i),1
               endif
            enddo
            write(ifile1,'(a4)') '    '
            write(ifile1,'(a4)') 'stop'
            if(iii.ne.0) then
               write(*,*) '******************'
               write(*,*) iii,' nodes found with no set temperature'
               write(*,*) 'see file temps_not_set.err'
               write(*,*) '******************'
            endif
            close(ifile1)
            close(ifile2)
            stop
         endif
c
      case (-439)
         if(l.eq.1) then
c     
c     print out volume fractions and coordinates(split zones)
c     print new zones 
c     
            write(*,*) 'k=-439, calculate new zones in YM area'
            write(*,*) 'HFM units must each be given  a unique rlp ',
     &           'model unit'
!           open(unit=99,file='volume_fractions.txt',status='unknown')
!           open(unit=89,file='zones_near_YM.macro',status='unknown')
            inquire (file = 'user_ymp.ctl', exist = ex)
            if (ex) then
               ifile1 = open_file('user_ymp.ctl','old')
            else
               write(*,*) 'Enter name of control file'
               read(*,'(a80)') dum_user
               ifile1 = open_file(dum_user,'old')
            end if
            allocate(idum(200))
! Volume fractions file
            read (ifile1,'(a80)') dum_user
! Zone macro file
            read (ifile1,'(a80)') water_file
            ifile2 = open_file(water_file,'unknown')
!            write(*,*) ' ENTER x1,x2 and y1,y2 to define bounding box'
            read(ifile1,*) x01,x02,y01,y02

            write(ifile2,'(a28,4f12.3)') 'zonn   # x01, x02, y01, y02 ',
     &           x01,y01,x02,y02
!            write(*,*) ' ENTER the number of units followed by ', 
!     &           'the desired unit ID numbers'
!            write(*,*) ' FORM: total_number of units, then ',
!     &           'ID1, ID2, ID3 .... IDn'
            read(ifile1,*) iii,(idum(i),i=1,iii)
            close (ifile1)
            ifile1 = open_file(dum_user,'unknown')

            write(*,*) ' calculating volume percentages....'
            do  kb=1,iii
               j = idum(kb)
               sum_volj = 0.0d00
               sum_volj100 = 0.0d00
               ii = 0
               iirb = 0

               do i=1,neq
                  if(irlp(i).eq.j) then  
                     if(cord(i,1).gt.x01.and.cord(i,1).lt.x02) then  
                        if(cord(i,2).gt.y01.and.cord(i,2).lt.y02) then
                           ii = ii+1
                           iirb(ii) = i
                           sum_volj100 = sx1(i)+sum_volj100
                           go to 346
                        endif
                     endif
                     sum_volj = sx1(i)+sum_volj
 346                 continue             
                  endif
               enddo
               
               if(ii.ge.1) then
                  write(ifile2,'(i6)') j+200
                  write(ifile2,'(a4)') 'nnum'
                  write(ifile2,'(i6)') ii
                  write(ifile2,'(8i10)') (iirb(jj),jj=1,ii)
               endif
               write(*,*) 'unit ',j,' volume percent  ',
     &              sum_volj/vtot*100.0d00
               write(*,*) 'unit ',j+200,' volume percent  ', 
     &              sum_volj100/vtot*100.0d00
               write(ifile1,*) 'unit ',j,' volume percent  ', 
     &              sum_volj/vtot*100.0d00
               write(ifile1,*) 'unit ',j+200,' volume percent  ',
     &              sum_volj100/vtot*100.0d00
            enddo
         endif
         xmin=1.d30
         ymin=1.d30
         do i=1,neq
            xmin = min(xmin,cord(i,1))
            ymin = min(ymin,cord(i,2))
         enddo
         write(ifile2,*) ' '
         write(ifile2,'(a4)') 'stop'
         
         close(ifile1)
         close(ifile2)
         stop
      case(455)
c read permeability from daniil and out in permeability file   
        write(*,*) '>>>> enter raw permeability file  <<<<<<'
        read(*,'(a80)') dum_user  
        write(*,*) '>>>> enter nx,ny  <<<<<<'
        read(*,*) i1,i2          
        iuser1 = open_file(dum_user,'old')
        iuser2 = open_file('perm_fehm_ner.macro','unknown')
        icount = 0
        i3 = i1*i2 + 1
        write(iuser2,'(a4)') 'perm'
455     continue   
        read(iuser1,*,end = 456) permxd
        icount = icount + 1
        i3 = i3 -1
        write(iuser2,457) i3, i3, 1, permxd,permxd,permxd
        go to 455
456     continue  
        write(iuser2,*) ' '
        write(*,*) 'perm data read and transformed, n = ',icount
        write(*,*) 'file perm_fehm_ner.macro created'
        pause
        stop
457     format(3(1x,i8),1p,3(1x,g14.6))   
      case(-900)
      open(unit = 97,file='fehm_volumes',status='unknown')  
      sum_volj = 0.0
      do i = 1,n
       sum_volj = sum_volj + sx1(i)
       write(97,'(i8,1x,f15.6)') i,sx1(i)
      enddo 
       write(97,*) 'volume total ', sum_volj
      stop
      case(-901)
       if(l.eq.1) then
        isalt = 1
        open(unit = 98,file='pres_vap.dat',status='unknown')  
        open(unit = 99,file='pres_vap.out',status='unknown')  
        read(98,'(a80)') wdd(1:80)
        read(98,*) iii
c        read(98,*) t1, t0, avap1, avap2
        write (99,'(a80)') wdd(1:80)
        write (99,*) iii
        allocate(temp1(iii),ptemp1(iii),ptemp2(iii))
        allocate(ptemp12(iii))
        if(allocated(an)) deallocate(an)
        allocate(an(iii))
        deallocate(phi)
        allocate(phi(iii))
        deallocate(pcp)
        allocate(pcp(iii))
        pcp = 0.
        deallocate(pci)
        allocate(pci(iii))
        deallocate(t)
        allocate(t(iii))
        deallocate(ieos)
        allocate(ieos(iii))
        ieos = 3
        deallocate(dpcef)
        allocate(dpcef(iii))
         read(98,*) (temp1(i), i = 1,iii)

c loop on an
         dela = 0.1
         dand = 2.  
         and = -dand 
         do j = 1, 11
          and = and + dand
          an = and
                     
          ms = an(1) * 58.55 / 1000.
          xf = ms / ( ms + 1.0 ) 
          write(99,*)    
          write(99,*) 'conc',  an(1), ms, xf  
         do i = 1, iii
          phi(i) = 1.0
          t(i) = temp1(i)
          ivaprsalt = 2
          call saltctr(1,i,tdum,tdum) 
          ptemp1(i) = phi(i) - pci(i)
          ivaprsalt = 2
          an(i) = an(i)+dela
          call saltctr(1,i,tdum,tdum) 
          ptemp2(i) = phi(i) - pci(i)

          write(99,984) t(i), ptemp1(i),(ptemp2(i)-ptemp1(i))/dela
          go to 371
          tdum = (t(i)-t0+1.e-10)/(t1-t0)
c          ptemp12(i) = an(i)*avap2*exp(avap1*tdum)
c          ptemp12(i) = an(i)*(avap1*tdum**4 + avap2)
c simple inverse t model
           cstar = 5
           c1star = 1./(1.+cstar)
           c2star = 1./cstar
           astar = (avap1-avap2)/(c2star-c1star)
           bstar = avap1-astar*c2star
           ms = an(i) * 58.55 / 1000.
	     xf = ms / ( ms + 1.0 )
           ptemp12(i) = 
     &      xf**0.95*(astar*(1./(1.-tdum**3+cstar)) + bstar)
c
          write(99,984) 
     &    t(i), ptemp1(i), ptemp2(i), xf, ptemp2(i)-ptemp12(i),
     &    ptemp1(i)-ptemp2(i)+ptemp12(i),ptemp12(i)
984       format(7(1x,f15.5))     
371      continue
         enddo
         enddo
         close (98) 
         close (99)         
         stop
       endif
      case(-909)   
c
c compare temperatures for nts_thermal verification
c 
      if(l.eq.1) then   
        ex = .false.
        inquire(file = 'temp_test.txt', exist = ex) 
        if(.not.ex) then
          if(iout.ne.0) 
     &      write (iout,*)'file 1 for nts_thermal comp missing'            
          if(iptty.ne.0) 
     &      write (iptty,*)'file 1 for nts_thermal comp missing'
          stop
        endif
        open(unit = 98,file='temp_test.txt',status='unknown')  
        ex = .false.
        inquire(file = 'temp_test_LANL.txt', exist = ex) 
        if(.not.ex) then
          if(iout.ne.0) 
     &      write (iout,*)'file 2 for nts_thermal comp missing'            
          if(iptty.ne.0) 
     &      write (iptty,*)'file 2 for nts_thermal comp missing'
          stop
        endif
        open(unit = 99,file='temp_test_LANL.txt',status='unknown')  
c
        allocate(temp1(7))
        read(98,*) 
        read(98,*) 
        read(98,*) 
        read(98,*) 
        read(98,*) 
        read(99,*) 
        read(99,*) 
        read(99,*) 
        read(99,*) 
        read(99,*) 
c
        temp1(4) = 0.0
        ntimes = 0
        kb_max = 0
985     continue        
        read (98,*, end = 986) i, dumx, temp1(1)
        read (99,*, end = 986) j, dumx, temp1(2)
        ntimes = ntimes + 1
        temp1(3) = abs(temp1(1)-temp1(2))
        temp1(7) = temp1(7) + temp1(3)
        if(temp1(3).gt.temp1(4)) then
            temp1(4) = temp1(3)
            temp1(5) = temp1(1)
            temp1(6) = temp1(2)
            kb_max = i
        endif
        go to 985
986     continue   
        
          
           if(iout.ne.0) write(iout,*) 
     &        'max temp diff node ',kb_max,' diff = ', temp1(4)
           if(iptty.ne.0) write(iptty,*) 
     &       'max temp diff node ',kb_max,' diff = ', temp1(4)            
           write (iout,*) ntimes,' lines , end of data,stopping'            
          if(iptty.ne.0) 
     &      write (iptty,*) ntimes,' lines , end of data,stopping' 
          stop        
      endif
      case(-910)
c
c vap_press fit
c
      if(l.eq.1) then
c
        open(unit = 98,file='pres_vap3.dat',status='unknown')  
        open(unit = 99,file='pres_vap2.out',status='unknown')  
c
      read(98,'(a80)') wdd(1:80)
      read(98,*) t1
      read(98,'(a80)') wdd(1:80)
      read(98,*) iii
c read temperature data
      allocate(temp1(iii),ptemp1(iii),ptemp2(iii),ptemp12(iii))
      deallocate(an)
      allocate(an(iii))
      read(98,*) (temp1(i), i = 1,iii)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,*) dla0, dlpa1, dlpa2, dlpa3, dlta1, dlta2, dlta3,
     &   dlpta, dlpt2a, dlp2ta
      read(98,*) dlb0, dlpb1, dlpb2, dlpb3, dltb1, dltb2, dltb3,
     &   dlptb, dlpt2b, dlp2tb
      read(98,*) bcoef1

c     loop on tracer
         dela = 0.1
         dand = 2.  
         and = -dand 
         do j = 1, 11
          and = and + dand    
          an = and         
          ms = and*58.55 / 1000.
          xf = ms / ( ms + 1.0 ) 
          xm = xf 
          xm2 = xm*xm
          xm3 = xm2*xm
          write(99,786) and,ms,xf
786       format('concentration = ',5(1x,f15.8))
787       format(3(1x,g20.10)) 
c     loop on temperture
         isalt = 1
         do i = 1, iii   
          tdum = temp1(i)
          t(i) = tdum
          ivaprsalt = 2
          call saltctr(1,i,tdum,tdum) 
          ptemp1(i) = phi(i) - pci(i)
c get sat pressure
          ivaprsalt = 4
          call saltctr(1,i,tdum,tdum) 
          ptemp2(i) = phi(i) - pci(i)
          tl = tdum/t1
          tl2 = tl*tl 
          tl3 = tl2*tl
          tlxm = tl*xm
          tl2xm = tl2*xm
          tlxm2 = tl*xm2
          pvwn1=dla0+dlpa1*xm+dlpa2*xm2+dlpa3*xm3
          pvwn2=dlta1*tl+dlta2*tl2+dlta3*tl3
          pvwn3=dlpta*tlxm+dlpt2a*tl2xm+dlp2ta*tlxm2
          pvwn=pvwn1+pvwn2+pvwn3
          pvwd1=dlb0+dlpb1*xm+dlpb2*xm2+dlpb3*xm3
          pvwd2=dltb1*tl+dltb2*tl2+dltb3*tl3
          pvwd3=dlptb*tlxm+dlpt2b*tl2xm+dlp2tb*tlxm2
          pvwd=pvwd1+pvwd2+pvwd3
c          pvw0=pvwn/pvwd
c
c power law
c
          

          ptemp12(i) = pvw0*(xm**bcoef1)
        enddo
         write(99, 787) (temp1(k),ptemp2(k)-ptemp12(k),
     &      ptemp1(k), k = 1,iii)
         write(99,*)
        enddo
        close(98)
        close(99)
        stop
      endif
      case(-911)
c
c vap_press fit
c
      if(l.eq.1) then
c
        open(unit = 98,file='pres_vap.dat',status='unknown')  
        open(unit = 99,file='pres_vap.out',status='unknown')  
c
      read(98,'(a80)') wdd(1:80)
      read(98,*) iii
c read temperature data
      allocate(temp1(iii),ptemp1(iii),ptemp2(iii),ptemp12(iii))
      deallocate(an)
      allocate(an(iii))
      read(98,*) (temp1(i), i = 1,iii)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,*) t1,cstar,avap1,avap2
      read(98,*) bcoef1,bcoef2

c     loop on tracer
         dela = 0.1
         dand = 2.  
         and = -dand 
         do j = 1, 11
          and = and + dand    
          an = and         
          ms = and*58.55 / 1000.
          xf = ms / ( ms + 1.0 ) 
          xm = xf 
          xm2 = xm*xm
          xm3 = xm2*xm
          write(99,786) and,ms,xf
c786       format('concentration = ',5(1x,g15.10))
c787       format(2(1x,g20.10)) 
c     loop on temperture
         isalt = 1
         do i = 1, iii   
          tl = temp1(i)
          t(i) = tl
          ivaprsalt = 2
          call saltctr(1,i,tdum,tdum) 
          ptemp1(i) = phi(i) - pci(i)
c get sat pressure
          ivaprsalt = 4
          call saltctr(1,i,tdum,tdum) 
          ptemp2(i) = phi(i) - pci(i)
           tdum = tl/t1
           c1star = 1./(1.+cstar)
           c2star = 1./cstar
           astar = (avap1-avap2)/(c2star-c1star)
           bstar = avap1-astar*c2star
 
           pvw = 
     &      xm**bcoef1*(astar*(1./(1.-tdum**bcoef2+cstar)) + bstar)

          ptemp12(i) = pvw
        enddo
         write(99, 787) (temp1(k),ptemp2(k)-ptemp12(k),
     &      ptemp1(k), k = 1,iii)
         write(99,*)
        enddo
        close(98)
        close(99)
        stop
      endif
      case(-912)
c
c vap_press fit
c
      if(l.eq.1) then
c
        open(unit = 98,file='pres_vap3.dat',status='unknown')  
        open(unit = 99,file='pres_vap3.out',status='unknown')  
c
      read(98,'(a80)') wdd(1:80)
      read(98,*) t1
      read(98,'(a80)') wdd(1:80)
      read(98,*) iii
c read temperature data
      allocate(temp1(iii),ptemp1(iii),ptemp2(iii),ptemp12(iii))
      deallocate(an)
      allocate(an(iii))
      read(98,*) (temp1(i), i = 1,iii)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,'(a80)') wdd(1:80)
      read(98,*) dla0, dlpa1, dlpa2, dlpa3, dlta1, dlta2, dlta3,
     &   dlpta, dlpt2a, dlp2ta
      read(98,*) dlb0, dlpb1, dlpb2, dlpb3, dltb1, dltb2, dltb3,
     &   dlptb, dlpt2b, dlp2tb
      read(98,*) bcoef1

c     loop on tracer
         dela = 0.1
         dand = 2.  
         and = -dand 
         do j = 1, 11
          and = and + dand    
          an = and         
          ms = and*58.55 / 1000.
          xf = ms / ( ms + 1.0 ) 
          xm = xf 
          xm2 = xm*xm
          xm3 = xm2*xm
          write(99,786) and,ms,xf
c786       format('concentration = ',5(1x,g15.10))
c787       format(3(1x,g20.10)) 
c     loop on temperture
         isalt = 1
         do i = 1, iii   
          tdum = temp1(i)
          t(i) = tdum
          tl = tdum/t1
          ivaprsalt = 2
          call saltctr(1,i,tdum,tdum) 
          ptemp1(i) = phi(i) - pci(i)
c get sat pressure
          ivaprsalt = 4
          call saltctr(1,i,tdum,tdum) 
          ptemp2(i) = phi(i) - pci(i)
c          pvwn1=dla0+dlpa1*xm**dlpa2
          pvwn1=dla0
          pvwn2=dlta1*tl**dlta2
          pvwn3=0.0
          pvwn=pvwn1+pvwn2+pvwn3
          pvwd1=dlb0+dlpb1*xm**dlpb2
          pvwd2=dltb1*tl**dltb2
          pvwd3=0.0
          pvwd=pvwd1+pvwd2+pvwd3
          pvw0=pvwn/pvwd 
c
c power law
c
          
          ptemp12(i) = dlpa1*xm**dlpa2*(pvw0+bcoef1)
        enddo
         write(99, 787) (temp1(k),ptemp2(k)-ptemp12(k),
     &      ptemp1(k), k = 1,iii)
         write(99,*)
        enddo
        close(98)
        close(99)
        stop
      endif
      case(-913)
c george zyvoloski 041319 (NEPTUNE PROJECT)          
c will add read voronoi areas - apply flowrates  
c make sure to close files          
       if(l.eq.1) then
c open  files for recharge (recharge.files)
          ex = .false.
          inquire(file = 'recharge.files', exist = ex) 
          if(ex) then
            ifile1 = open_file('recharge.files','old')  
699         read (ifile1, '(a200)',end=988) wdd1_t          
            if (wdd1_t(1:1) .eq. '#') go to 699
            read (wdd1_t, '(a4)') macro
            if (macro(1:4) .eq. 'stop') then
              go to 788
            elseif (macro(1:4) .eq. 'nblk') then
               read (wdd1_t(5:200), *) nblk
               allocate(node_a(nblk))
            elseif (macro(1:4) .eq. 'rech') then
               read(ifile1,*) recharge1, file_zone1
               ex = .false.
               inquire(file = file_zone1, exist = ex)   
               if(ex) then
                ifile2= open_file(file_zone1,'old')
               else
                 write(iout,239) 
                 stop
               endif
698            read (ifile2, '(a200)',end=988) wdd1_t 
               if (wdd1_t(1:1) .eq. '#') go to 698
               if (wdd1_t(1:3) .eq. 'zon') then
                read(ifile2,*) izone_recharge
                read(ifile2,'(a)',end=988) wdd1_t(1:4)
                if(wdd1_t(1:4).ne.'nnum') then
                 write(iout,240)
                 stop
                endif
               endif
               read(ifile2,*) num_rech_blks
               allocate(izone_node(num_rech_blks))
               read(ifile2,*) (izone_node(i), i =1, num_rech_blks)
            elseif (macro(1:4) .eq. 'area') then  
c just read "outside top  areas and nodes
c file_area1 is usually the LaGrit outside zone file 
c izone_num1 is the zone number in  file_area1 that is the top  
             read(ifile1,'(a)') wdd1_t(1:6)   
             if(wdd1_t(1:6).eq.'approx') then
              approx = .true.
              go to 699
             else 
              approx = .false.
              backspace ifile1
             endif
             read (ifile1,'(a)') file_area1
             ex = .false.
             inquire(file = file_area1, exist = ex)
             if(ex) then
              ifile2= open_file(file_area1,'old')
             else
              if(iout.ne.0) write(iout,239) 
             endif
             read (ifile1,*) izone_num1
             max_lines = 2*nblk
             do i = 1, max_lines
              read(ifile2,'(a4)') dum_user(1:4)
              if(dum_user(1:4).eq.'nnum') then
               backspace ifile2 
               backspace ifile2
               read(ifile2,*) izone_num0
               if(izone_num0.eq.izone_num1) then
                read(ifile2,*) dum_user(1:4)
                read(ifile2,*) idummy
                  read(ifile2,*)(node_a(j), j=1, idummy) 
                go to 766
               else 
                read(ifile2,*) dum_user(1:4)
                read(ifile2,*) dum_user(1:4)
                read(ifile2,*) dum_user(1:4)                  
               endif
              endif
             enddo
766    continue             
c file_area2 is usually the LaGrit outside voronoi area file 
c izone_num1 is the zone number in  file_area2 that is the top 
c idir_zone in the direction of recharge for file_area2. usually=3    
c
             read (ifile1,'(a)') file_area2
             ex = .false.
             inquire(file = file_area2, exist = ex)
             if(ex) then
              ifile3= open_file(file_area2,'old')
             else
              if(iout.ne.0) then
               write(iout,239) 
              endif             
             endif   
            read (ifile1,*) izone_num2, idir_zone 
             allocate(temp2(nblk,3))
             do i = 1, max_lines
              read(ifile3,'(a4)') dum_user(1:4)
              if(dum_user(1:4).eq.'nnum') then
               backspace ifile3 
               backspace ifile3
               read(ifile3,*) izone_num0
               if(izone_num0.eq.izone_num2) then
                read(ifile3,*) dum_user(1:4)
                read(ifile3,*) idummy
                 read(ifile3,*)
     &            (temp2(j,1),temp2(j,2),temp2(j,3), j=1, idummy)
                go to 767
               else 
                read(ifile3,*) dum_user(1:4)
                read(ifile3,*) dum_user(1:4)
                read(ifile3,*) dum_user(1:4)                  
               endif
              endif
             enddo
767          continue  
   
c             
       else if(macro(1:4).eq.'outp') then
c
c merge data and write as FEHM file
c assume recharge (recharge1) is m/yr converted to kg/sec
c yr/sec= 3.170979e-8,recharge2= (kg/sec)/area 
        recharge2 = recharge1*3.170979e-8*997.
        read (ifile1,'(a)') file_recharge_out
        ifile4= open_file(file_recharge_out,'unknown')
c first fill in array (nblk-sized) of surface grid blocks and surface areas  
c    
        allocate (area_surface(nblk))
        area_surface = 0.0d0
        allocate (sk_save(nblk))
        sk_save = 0.0d0
        area_total = 0.0
        if(approx.eqv..false.) then
c area obtained from outside zone list and outside_vor.area            
        write(ierr,*) '>>>>> top outside zone >>>>> ',izone_num2
        do j = 1, idummy
          i = node_a(j)
          area_surface(i) = temp2(j,idir_zone)
          area_total = area_total + area_surface(i)
         write(ierr,'(a7,i8,a2,g12.3,a6,3g16.6)')'node_o ', i ,' area ',
     &   area_surface(i),' x y z = ', (cord(i,i1), i1 = 1,3) 
        enddo
        write(ierr,242) izone_num2, area_total
        write(ierr,*)
        write(ierr,*)  'analysis of zone ', izone_num0
        else
c area obtained from gridblock volume and gridblock thickness
         do j = 1, num_rech_blks
          i = izone_node(j)
          area_surface(i) = volume(i)/dzrg(i)
          area_total = area_total + area_surface(i)
         enddo
        endif

       iii = 0
       area_total =0.0
       do j = 1, num_rech_blks
        i = izone_node(j)
        sk_save(i) = -recharge2*area_surface(i)
        if(sk_save(i).ne.0.0) iii = iii+1
        write(ierr,'(a7,i8,a2,g12.3,a6,3g16.6)') 'node_r ', i ,' area ',
     &   area_surface(i),' x y z = ', (cord(i,i1), i1 = 1,3)   
         area_total = area_total + area_surface(i)
       enddo
       write (ierr,*) 'non zero sources ', iii,' area_tot',
     &  area_total
c
c  write FEHM flow macro

c 
       write(ifile4,*) 'flow  ',
     &   '    # created from Aaron Bandler analysis'
       do j = 1, num_rech_blks
        i = izone_node(j)
        write(ifile4,241) i,i,1,sk_save(i),1.0,0.0
       enddo
       write(ifile4,'(a4)') '    '
       close(ifile4)
       endif
       go to 699
788    if(iout.ne.0) then
        write(iout,238) 
c        close(ifile1) 
c        close(ifile2)
c        close(ifile3)
c        close(ifile4)
        stop
       endif 
       
988    continue
       if(iout.ne.0) write(iout,*) 
     &     'EOF found recharge.files (user -913)'
       stop
238    format ('stop record found in recharge.files',
     &          ', ending recharge calcs')
239   format('recharge.files is missing, stopping')
240   format('zone recharge file wrong format(need nnum), stopping')
241   format(1x,i8,1x,i8,1x,i4,1x,g14.5,1x,f8.2,1x,g10.2)    
242   format('voronoi zone number ', i5, ' total_area ', 1p, g14.5)      
      endif
      endif
      case(-914)
c george zyvoloski 043019 (NEPTUNE PROJECT)          
c create new air zone that 
c make sure to close files 
       if(l.eq.1) then
       if(iout.ne.0) write(iout,*) '>>>>>> starting  user -914 '
       if(iptty.ne.0) write(iptty,*) '>>>>> starting user -914 '           
       por_air = 5.55d-5
       iii = 0
       ii = 0
       i3 = 0      
       allocate(node_a(n))
       allocate(node_b(n))
       allocate(izone_node(n))
       node_a = 0
       node_b = 0
       izone_node = 0
        do i = 1,n
          if(abs(ps(i)-por_air).le.1.d-8) then
           if(ka(i).ne.0) then
            i3 = i3 +1    
            node_b(i3) = i
            write(ierr,'(a17,1x,i8,1p,3(1x,g14.4))') 
     &         'flow bc in air nodes', i, sk(i), pflow(i), wellim(i)
           else
             iii = iii+ 1   
             izone_node(iii) = i
              i1 = nelm(i)+1
              i2 = nelm(i+1)
              do j = i1,i2
               kb = nelm(j)
               if(ps(kb).ne.por_air) then
                ii = ii +1
                node_a(ii) = i
                go to 261
               endif
              enddo
           endif
          endif
261     continue        
        enddo
c write files for air zones 
       ifile1 = open_file('air_total_401_gaz.zonn','unknown')   
       ifile2 = open_file('air_near_402_gaz.zonn','unknown') 
       ifile3 = open_file('air_src_403_gaz.zonn','unknown')   
       write(ifile1,*) 'zonn'
       write(ifile1,*) ' 401'
       write(ifile1,*) 'nnum'
       write(ifile1,'(i9)') iii
       write(ifile1,'(10(1x,i9))') (izone_node(i), i = 1, iii)
       write(ifile1,*) '    '
       close(ifile1)
       write(ifile2,*) 'zonn'
       write(ifile2,*) ' 402'
       write(ifile2,*) 'nnum'
       write(ifile2,'(i9)') ii
       write(ifile2,'(10(1x,i9))') (node_a(i), i = 1, ii)
       write(ifile2,*) '    '
       close(ifile2)
       write(ifile3,*) 'zonn'
       write(ifile3,*) ' 403'
       write(ifile3,*) 'nnum'
       write(ifile3,'(i9)') i3
       write(ifile3,'(10(1x,i9))') (node_b(i), i = 1, i3)
       write(ifile3,*) '    '
       close(ifile3)
c calculate nodes in air-model interface
       node_a = 0
       node_b = 0
       izone_node = 0
       if(.not.allocated(izone_awt)) allocate(izone_awt(n))
       izone_awt =0
       ii = 0
        do i = 1,n
          if(abs(ps(i)-por_air).le.1.d-8) then
              i1 = nelm(i)+1
              i2 = nelm(i+1)
              do j = i1,i2
               kb = nelm(j)
               if(abs(ps(kb)-por_air).gt.1.d-8) then
                izone_awt(kb) = izone_awt(kb) + 1
                write(ierr,'(a17,1x,i7,1x,i7)') 'found model node ', kb,
     &            izone_awt(kb)
               endif
              enddo
           endif
271        continue        
        enddo
c node_a() contains the rind (air-model interface)
      write(ierr,*) ' air model interface contains nodes = ', ii
       por_tbu = 0.11d0
       por_wlt = 0.12d0
       por_cddl = 0.223d0
       i1 = 0
       i2 = 0
       i3 = 0
       area_tbu = 0.0
       area_wlt = 0.0
       area_cddl = 0.0
       do i = 1, n
        j = izone_awt(i)
        if(j.ne.0) then
         if(abs(ps(i)-por_tbu).le.1.d-8) then
          i1 = i1 +1
          node_a(i1) = i
          area_tbu =  area_tbu + volume(i)/dzrg(i)
         else if(abs(ps(i)-por_wlt).le.1.d-8) then
          i2 = i2 +1
          node_b(i2) = i
          area_wlt =  area_wlt + volume(i)/dzrg(i)
         else if(abs(ps(i)-por_cddl).le.1.d-8) then
          i3 = i3 +1
          izone_node(i3) = i  
          area_cddl =  area_cddl + volume(i)/dzrg(i)
         endif
        endif
       enddo
       write(ierr,*) 'node_tbu ', i1,' node_wlt ', i2,
     &   ' node_cddl ', i3       
       write(ierr,*) 'area_tbu ', area_tbu,' area_wlt ', area_wlt,
     &   ' area_cddl ', area_cddl
c write files for air zones 
       ifile1 = open_file('tbu_rechg_gaz_222.zonn','unknown')   
       ifile2 = open_file('wlt_rechg_gaz_221.zonn','unknown') 
       ifile3 = open_file('cddl_rechg_gaz_223.zonn','unknown')   
       write(ifile1,'(a)') 'zonn'
       write(ifile1,'(a)') ' 222'
       write(ifile1,'(a)') 'nnum'
       write(ifile1,'(i9)') i1
       write(ifile1,'(10(1x,i9))') (node_a(i), i = 1, i1)
       write(ifile1,'(a)') '    '
       close(ifile1)
       write(ifile2,'(a)') 'zonn'
       write(ifile2,'(a)') ' 221'
       write(ifile2,'(a)') 'nnum'
       write(ifile2,'(i9)') i2
       write(ifile2,'(10(1x,i9))') (node_b(i), i = 1, i2)
       write(ifile2,'(a)') '    '
       close(ifile2)
       write(ifile3,'(a)') 'zonn'
       write(ifile3,'(a)') ' 223'
       write(ifile3,'(a)') 'nnum'
       write(ifile3,'(i9)') i3
       write(ifile3,'(10(1x,i9))') (izone_node(i), i = 1, i3)
       write(ifile3,'(a)') '    '
       close(ifile3)       
       endif
c l = 0 ifblock ends
       if(iout.ne.0) write(iout,*) '>>>>>> stopping user -914 '
       if(iptty.ne.0) write(iptty,*) '>>>>> stopping user -914 '
       stop
      end select
      end
