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
      integer counter, dumcyc, readflag, incf, ntimes, node_a, node_b
      integer, allocatable :: idum(:)
      integer, allocatable :: idum_b(:)
      integer, allocatable :: izone_awt(:)
      real*8, allocatable :: cons(:,:)
      real*8, allocatable :: times(:)
      real*8 tdum
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
      real*8, allocatable :: sk_save(:)
      real*8, allocatable :: sk_save1(:)
      real*8, allocatable :: sk_save2(:)
      real*8, allocatable :: temp1(:), temp2(:,:)
      real*8, allocatable :: water(:)
      real*8, allocatable :: x_wt(:),y_wt(:)
      character*80 dum_user
      character*80 water_file, above_wt_file, outside_zone_file
      logical matrix_node, null1, xy, ex

      save ntimes, times, cons, counter, dumcyc

      select case (k)
c RJP added for turning off CO2 injection wells
      case (999)
         if(readflag.NE.1) then
            incf = open_file('co2_inj.txt','unknown')
            read(incf,*) shut_time
            close (incf)
            readflag = 1
         end if	
         if(days.ge.shut_time*365.25) then
            do i = 1, neq
               if(skco2(i).lt.0.d0) skco2(i) = 0
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
c-----------------------------------------------------------------------
c      test the analytical solution
      case (666)
         if(readflag.NE.1) then
            incf = open_file('nodes_conc_fixed.macro','unknown')
            read(incf,*) ntimes, node_a, node_b
            allocate(times(ntimes),cons(ntimes,2))
            do iii = 1,ntimes
               read(incf,*) times(iii), cons(iii,1), cons(iii,2)
            end do
            close (incf)
            counter = 1
            dumcyc = 1
            readflag = 1
            pcnsk(node_a) = -1
            pcnsk(node_b) = -1
         end if
         tdum = days - (dumcyc -1)*times(ntimes)
         if(tdum.EQ.times(counter+1)) counter = counter + 1
         cnsk(node_a) = -cons(counter,1)
         cnsk(node_b) = -cons(counter,2)
         if(counter.EQ.ntimes) then
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
!     &           'xw,xe,ys,yn'
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
 876        format('base',i6,' new node ',i8,' x ',g13.7,' y ',g13.7,
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
               call headctr(2,0,0.0,0.0)        
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

      end select
      end
