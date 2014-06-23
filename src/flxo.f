      subroutine flxo(iz)
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
CD1  This subroutine handles the reading and writing of internode 
CD1  fluxes, either from the input file or for output during a run.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/flxo.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:06   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.14   Mon Mar 31 08:35:40 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.13   Wed Jun 12 15:21:00 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.12   Mon Jun 03 11:17:48 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.11   Fri May 31 10:55:50 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.10   Wed May 08 14:15:12 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.9   Thu Apr 04 10:36:42 1996   hend
CD2 Increased Precision of velocity output
CD2 
CD2    Rev 1.8   Tue Apr 02 10:55:54 1996   hend
CD2 Corrected Sign Error in Output
CD2 
CD2    Rev 1.7   Wed Feb 07 12:05:40 1996   gaz
CD2 changes to protect from divide by zero
CD2 
CD2    Rev 1.6   Wed Jan 10 11:15:24 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.5   Wed Jan 10 09:26:56 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   09/11/95 17:29:52   awolf
CD2 seh fixed indexing of matrix nodes for avs output
CD2 
CD2    Rev 1.3   08/07/95 11:15:28   awolf
CD2 Now writes internode velocities and fluxes for DPDP cases
CD2 too. Also write fluxes or velocities between matrix and fracture
CD2 nodes. Fluxes based on a_axy term
CD2 
CD2    Rev 1.2   03/10/95 10:38:50   llt
CD2 changed extraction of fluxes or velocities - gaz
CD2 
CD2    Rev 1.1   03/18/94 15:43:06   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:56   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
c Henderson --> 11 Sep 1995
c This code currently computes a darcy velocity. If you want the 
c velocity particles will travel at, for instance,
c 
c FRACTURES: (Computed Velocity) / (porosity*sat*fracture volume frac)
c
c MATRIX: (Computed Velocity) / (porosity*sat)
c

      use combi
      use comci
      use comdi
      use comei
      use comflow
      use davidi
      use comgi
      use comfi
      use comdti
      use comxi
      use comai
      use comco2, only : icarb, c_axy, c_vxy
      implicit none

      integer chngsign, iroot, ifile(2), open_file
      integer i, i1, i2, if, iff1, ii, iii, iz, j, jj, jjj, kb
      integer neqp1, nflx_old, nmatavw, nmat2avw,  node1, node2
      real*8  area_t, axyf
      real*8  coordx1, coordy1, coordz1, coordx2, coordy2, coordz2
      real*8  cosx, cosy, cosz
      real*8  dili, dilkb, divi, divkb, dis, dis_max
      real*8  fid, fid1, flx12vd, flx12ld, tol, tolf, vxyf
      real*8  x1, x2, y1, y2, z00, z1, z2
      real*8  x_flx, y_flx, z_flx
      real(8) :: axy = 0., vxy = 0.
      character(200) :: internode_file_root, internode_file
      logical null1
      parameter(z00=0.0d00, tolf = 1.d-50)
      save ifile

c next 3 lines added for andy flow rate arrays
      neqp1=neq+1
      nmatavw=nelm(neqp1)-neqp1
      nmat2avw=2*nmatavw

      if(nflxt.eq.0) return

      if(iz.eq.0) then
c
c        input section
c

         iflxc = iflxc +1
         if(iflxc.eq.1) then
            if (null1(root_name)) then
               if (nmfil(5)  .ne. nmfily(3) .and. nmfil(5) .ne. ' ') 
     &              then
! Use output file prefix to name internode flux output file if it exists
                  call file_prefix(nmfil(5), iroot)
                  internode_file_root = nmfil(5)(1:iroot)
               else 
! Use input file prefix to name internode flux output file
                  call file_prefix(nmfil(2), iroot)
                  internode_file_root = nmfil(5)(1:iroot)
              endif
            else
               iroot = len_trim (root_name)
               internode_file_root = root_name(1:iroot)
            end if
            internode_file = internode_file_root(1:iroot) //
     &           '.internode_fluxes.out'
            ifile(1) = open_file (internode_file, 'unknown')
            if (icarb .ne. 0) then
               internode_file = internode_file_root(1:iroot) //
     &               '.internode_co2_fluxes.out'
               ifile(2) = open_file (internode_file, 'unknown')
            end if
            nflx_old = 0
         else
            nflx_old = nflxc(iflxc-1)
         endif
         read(inpt,*) nflx
         nflxc(iflxc) = nflx + nflx_old
         read(inpt,*) (iflx1(nflx_old+i)
     &                ,iflx2(nflx_old+i),i=1,nflx)
         do j=1,nflx
            if(iflx1(nflx_old+j).lt.0) then
               read(inpt,*) x_flx, y_flx, z_flx
               call near3(x_flx,y_flx,z_flx,node1,0)
               iflx1(nflx_old+j)=node1
            else if (iflx1(nflx_old+j).gt.n0) then
               write (ierr, 6008)
               write (ierr, 6009) iflx1(nflx_old+j), n0
               if (iout .ne. 0) write(iout, 6009) iflx1(nflx_old+j), n0
               if (iptty .gt. 0) 
     .              write (iptty, 6009) iflx1(nflx_old+j), n0
               stop
            endif
            if(iflx2(nflx_old+j).lt.0) then
               read(inpt,*) x_flx, y_flx, z_flx
               call near3(x_flx,y_flx,z_flx,node1,0)
               iflx2(nflx_old+j)=node1
            else if (iflx2(nflx_old+j).gt.n0) then
               if (ivelo .lt. 0) then
                  write (ierr, 6008) 'dvel'
               else
                  write (ierr, 6008) 'flxo'
               end if
               write (ierr, 6009) iflx2(nflx_old+j), n0
               if (iout .ne. 0) write (iout, 6009) iflx2(nflx_old+j), n0
               if (iptty .gt. 0) 
     .              write (iptty, 6009) iflx2(nflx_old+j), n0
               stop
            endif
         enddo
 6008    format (' **** Invalid input: macro ', a4, ' ****')
 6009    format (' **** Invalid node specified, ', i8, 
     .        ' is greater than ', 'n0 (', i8, ' ): stopping ****')

      else if(iz.eq.1) then
         do if=1,nflxt
            if(iflx2(if).eq.0) then
               node1=iflx1(if)
               i1=nelmdg(node1)+1
               i2=nelm(node1+1)
               dis_max=-1.d50
                do j=i1,i2
                 kb=nelm(j) 
                     if(icnl.eq.0) then
                        dis=(cord(node1,1)-cord(kb,1))**2
     &                       + (cord(node1,2)-cord(kb,2))**2 
     &                       + (cord(node1,3)-cord(kb,3))**2
                     else
                        dis=(cord(node1,1)-cord(kb,1))**2
     &                       + (cord(node1,2)-cord(kb,2))**2 
                     endif                
                     if(dis.gt.dis_max) then
                      dis_max=dis
                      iflx2(if)=kb
                     endif
                enddo
            endif
         enddo
      else if(iz.eq.2.or.iz.eq.3) then
c     
c     extract fluxes or velocities
c     
         neqp1=neq+1
         tol=1.d-20
         flx12v = 0.0d00
         flx12l = 0.0d00
         do if=1,nflxt
            i=iflx1(if)
            j=iflx2(if)
            chngsign=0
            if(i.gt.j) then
               ii=i
               i=j
               j=ii
               chngsign=1
            endif
            iii=i
            jjj=j
            if(i.gt.neq) iii=i-neq
            if(j.gt.neq) jjj=j-neq
            do jj=nelmdg(iii)+1,nelm(iii+1)
               kb=nelm(jj)
               if(kb.eq.jjj.or.(i.le.neq.and.j.gt.neq)) then
                  iw=istrw(jj-neqp1)
c     use internode fluxes already stored         
c     axy=a_axy(jj-neqp1)
c     vxy=a_vxy(jj-neqp1)

                  if(iz.eq.2) then
                     if(i.gt.neq.and.j.gt.neq) then
                        axy=a_axy(jj-neqp1+nmatavw)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=a_vxy(jj-neqp1+nmatavw)
                     else if (i.le.neq.and.j.le.neq) then
                        axy=a_axy(jj-neqp1)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=a_vxy(jj-neqp1)
                     else
                        axy=a_axy(nmat2avw+i)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=a_vxy(nmat2avw+i)
                     endif
                  elseif(iz.eq.3) then
                     if(i.gt.neq.and.j.gt.neq) then
                        axy=c_axy(jj-neqp1+nmatavw)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=c_vxy(jj-neqp1+nmatavw)
                     else if (i.le.neq.and.j.le.neq) then
                        axy=c_axy(jj-neqp1)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=c_vxy(jj-neqp1)
                     else
                        axy=c_axy(nmat2avw+i)
                        if (irdof .ne. 13 .and. jswitch .eq. 0) 
     &                       vxy=c_vxy(nmat2avw+i)
                     endif
                  endif
                  
                  if(ivelo.lt.0) then
c     velocities         
                     if(icnl.eq.0) then
                        dis=sqrt((cord(iii,1)-cord(kb,1))**2
     &                       + (cord(iii,2)-cord(kb,2))**2 
     &                       + (cord(iii,3)-cord(kb,3))**2)
                     area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
                     else
                        dis=sqrt((cord(iii,1)-cord(kb,1))**2
     &                       + (cord(iii,2)-cord(kb,2))**2) 
                     area_t=-(sx(iw,isox)+sx(iw,isoy))*dis        
                     endif                
                     flx12l(if) = 0.0d00
                     if(abs(axy).gt.tol) then
                     if(axy.gt.0.0) then
                        fid=upwgt
                        fid1=1.0-fid
                     else if(axy.lt.0.0) then
                        fid1=upwgt
                        fid=1.0-fid1
                     else
                        fid=0.5
                        fid1=0.5
                     endif
                     dili=dil(i)
                     dilkb=dil(j)
                     axyf=max(tol,(fid*dilkb+fid1*dili))
                     if(abs(axy).gt.tol) then
                        flx12l(if)=axy/axyf/area_t
     &                       *(fid1*dil(i)/max(rolf(i),tol)+fid*
     &                       dil(j)/max(rolf(j),tol))
                     endif
                     endif
                     flx12v(if) = 0.0d00
                     if(abs(vxy).gt.tol) then
                     if(vxy.gt.0.0) then
                        fid=upwgt
                        fid1=1.0-fid
                     else if(vxy.lt.0.0) then
                        fid1=upwgt
                        fid=1.0-fid1
                     else
                        fid=0.5
                        fid1=0.5
                     endif
                     divi=div(i)
                     divkb=div(j)
                     vxyf=max(tol,(fid*divkb+fid1*divi))
                     if(abs(vxy).gt.tol) then
                        flx12v(if)=vxy/vxyf/area_t
     &                       *(fid1*div(i)/max(rovf(i),tol)+fid*
     &                       div(j)/max(rovf(j),tol))
                     endif
                     endif
                  else
c     mass fluxes         
                     flx12l(if)=axy  
                     flx12v(if)=vxy         
                  endif
                  if (chngsign.eq.1) then
                     flx12l(if)=-flx12l(if)
                     flx12v(if)=-flx12v(if)
                  endif
               endif
            enddo
         enddo  
c     
c     write out fluxes or velocities
c     
	if(iz.eq.2) then
         if(ivelo.ge.0) then
	      if (iptty.gt.0) write(iptty,20)
            if (iptty.gt.0) write(iptty,30)
            if (iout .ne. 0) write(iout,20)
            if (iout .ne. 0) write(iout,30)
         else
            if (iptty.gt.0) write(iptty,21)
            if (iptty.gt.0) write(iptty,35)
            if (iout .ne. 0) write(iout,21)
            if (iout .ne. 0) write(iout,35)
         endif
	elseif(iz.eq.3) then
         if(ivelo.ge.0) then
	      if (iptty.gt.0) write(iptty,22)
            if (iptty.gt.0) write(iptty,30)
            if (iout .ne. 0) write(iout,22)
            if (iout .ne. 0) write(iout,30)
         else
            if (iptty.gt.0) write(iptty,23)
            if (iptty.gt.0) write(iptty,35)
            if (iout .ne. 0) write(iout,23)
            if (iout .ne. 0) write(iout,35)
         endif
	endif


 20      format(/,20x,'Internode Fluxes')
 21      format(/,20x,'Internode Velocities')
 22      format(/,20x,'Internode Fluxes for CO2')
 23      format(/,20x,'Internode Velocities for CO2')	
 30      format(1x,'Node1  Node2',11x,'X,Y,Z Vapor Flux(kg/sec)  ',
     &        '        X,Y,Z Liquid Flux (kg/sec)')
 35      format(1x,'Node1  Node2',11x,'X,Y,Z Vapor Vel.(m/sec)  ',
     &        '        X,Y,Z Liquid Vel.(m/sec)')
         do jj =1,iflxc
         if(jj.eq.1) then
          iff1 = 1 
         else
          iff1 = nflxc(jj-1) + 1
         endif
         if (iout .ne. 0) write(iout,355) jj,nflxc(jj)-iff1+1
         if (iptty .gt. 0) write(iptty,355) jj,nflxc(jj)-iff1+1
         write(ifile(iz-1),356) jj, nflxc(jj)-iff1+1, l, days
355      format(10x,2i10,'   (call number,number of pairs)')
356      format(3i6,1x,e14.6, 
     &        ' (call number,number of pairs,time step,days)')
         do if=iff1,nflxc(jj)
                 kb=iflx1(if)
                 iii=iflx2(if)
                 if(iii.gt.neq) iii=iii-neq
                 if(kb.gt.neq) kb=kb-neq
                     if(icnl.eq.0) then
                        dis=sqrt((cord(iii,1)-cord(kb,1))**2
     &                       + (cord(iii,2)-cord(kb,2))**2 
     &                       + (cord(iii,3)-cord(kb,3))**2)
                        cosx=(cord(iii,1)-cord(kb,1))/dis
                        cosy=(cord(iii,2)-cord(kb,2))/dis
                        cosz=(cord(iii,3)-cord(kb,3))/dis
                     else
                        dis=sqrt((cord(iii,1)-cord(kb,1))**2
     &                       + (cord(iii,2)-cord(kb,2))**2) 
                        cosx=(cord(iii,1)-cord(kb,1))/dis
                        cosy=(cord(iii,2)-cord(kb,2))/dis
                     endif                
                 flx12vd=flx12v(if)
                 flx12ld=flx12l(if)
                 if(abs(flx12vd).lt.tolf) flx12vd = 0.0d00
                 if(abs(flx12ld).lt.tolf) flx12ld = 0.0d00
            write(ifile(iz-1),40) 
     &           iflx1(if),iflx2(if),flx12vd,flx12ld             
            if(icnl.eq.0) then
              if(iptty.gt.0) write(iptty,40) iflx1(if),iflx2(if),
     &           cosx*flx12v(if),cosy*flx12v(if),cosz*flx12v(if),
     &           cosx*flx12l(if),cosy*flx12l(if),cosz*flx12l(if)
              if (iout .ne. 0) write(iout,40) iflx1(if),iflx2(if),
     &           cosx*flx12v(if),cosy*flx12v(if),cosz*flx12v(if),
     &           cosx*flx12l(if),cosy*flx12l(if),cosz*flx12l(if)
            else if(icnl.ne.0) then
              if(iptty.gt.0) write(iptty,40) iflx1(if),iflx2(if),
     &           cosx*flx12v(if),cosy*flx12v(if),z00,
     &           cosx*flx12l(if),cosy*flx12l(if),z00
              if (iout .ne. 0) write(iout,40) iflx1(if),iflx2(if),
     &           cosx*flx12v(if),cosy*flx12v(if),z00,
     &           cosx*flx12l(if),cosy*flx12l(if),z00
            endif
         enddo
         enddo
c
c write out coordinate information 
            if (iptty.gt.0) write(iptty,45)
            if (iout .ne. 0) write(iout,45)
 45      format(1x,'Node1 ',11x,'X,Y,Z      Coordinates  ',
     &        'Node2          X,Y,Z        Coordinates')
         do jj =1,iflxc
         if(jj.eq.1) then
          iff1 = 1 
         else
          iff1 = nflxc(jj-1) + 1
         endif
         if (iout .ne. 0) write(iout,355) jj,nflxc(jj)-iff1+1
         if (iptty .gt. 0) write(iptty,355) jj,nflxc(jj)-iff1+1
         do if=iff1,nflxc(jj)
                 kb=iflx1(if)
                 iii=iflx2(if)
                 if(iii.gt.neq) iii=iii-neq
                 if(kb.gt.neq) kb=kb-neq
                     if(icnl.eq.0) then
                        coordx1=cord(kb,1)
                        coordy1=cord(kb,2)
                        coordz1=cord(kb,3)
                        coordx2=cord(iii,1)
                        coordy2=cord(iii,2)
                        coordz2=cord(iii,3)
                     else
                        coordx1=cord(kb,1)
                        coordy1=cord(kb,2)
                        coordz1=0.0       
                        coordx2=cord(iii,1)
                        coordy2=cord(iii,2)
                        coordz2=0.0       
                     endif                
              if(iptty.gt.0) write(iptty,50) iflx1(if),coordx1,
     &         coordy1,coordz1,iflx2(if),coordx2,coordy2,coordz2
              if (iout .ne. 0) write(iout,50) iflx1(if),coordx1,
     &         coordy1,coordz1,iflx2(if),coordx2,coordy2,coordz2
         enddo
         enddo
      endif
 40   format(i6,1x,i6,1x,6(e12.6,1x))
 50   format(i6,1x,3(f10.1,1x),i6,3(1x,f10.1))
      return
      end
      

