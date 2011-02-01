      subroutine load_omr_flux_array
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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

!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To .
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 14-Jan-02, Programmer: S Kelkar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/load_omr_flux_array.f_a  $
!D2    Rev 2.4   29 Jan 2003 09:09:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
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

      use comai
      use combi
      use comci
      use comdti
      use comdi
      use comflow
      use comsptr
      use comsk
      use davidi
 
      implicit none

      integer k,ikb,kb,ipos,iwsk,ij,ij2,i3,i,list_max
      integer icnl_subst,lower_limit,upper_limit,ism
      integer find_index, numomr, open_file, ifile, iir
      integer irray0, flag_por

c FOR WTSI problems
c       flag_sat_wtsi=1, use average saturation for area calculations
c       flag_sat_wtsi=2, use upstream saturation for area calculations
      integer flag_sat_wtsi
      parameter(flag_sat_wtsi=2)

      real*8 xface1,xface2,xface3,face1,face2,face3,xxx
      real*8 xcos,ycos,zcos

      real*8 gggin,gggout,facin,facout,axy,sxkb,d,area
      real*8 axytotal,masstotal
      real*8 epsilon
      real*8 fxl,fxr,fyf,fyb,fzt,fzb

      real*8 sum_k_mass,epsilong,gotcord,gotcord2
c...s kelkar 3/4/04 3d omr.....
      integer ibou,n6lsq
      real*8 b6(6),b5(5),b4(4),b3(3),dmass_6lsq,taxy,sqaxy

c s kelkar may 13 2009
      real*8 facex1,facex2,facey1,facey2,facez1,facez2

c............................

c....s kelkar  march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular interior node, not a source/sink
c     irray(i,0) = -i : regular node, producing well, particle capture
c                        but not explicitly specified in sptr macro
c     irray(i,0) < -10000000 : regular interior node that is 
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c                = -(i+2000000) : spring node
c                = -(i+1000000) : well-capture node on extrernal bound
c                   simillar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c.........................................................

c      ifile = open_file('axy.dbg','unknown')

      list_max=200

      epsilong=1.e-20
      epsilon=1.e-10
c.......................................................

      icnl_subst = icnl
      
c     ************ put sinks on bdry ********
      
      if(icnl.eq.0) then
         lower_limit = -3
         upper_limit = 3
      else
         lower_limit = -2
         upper_limit = 2
      end if
      do i=1,neq
         if (irdof .ne. 13 .and. ifree .eq. 0) then
            xxx=+ps_trac(i)*rolf(i)*s(i)
         else
            xxx=+ps_trac(i)*rolf(i)
         endif

c....s kelkar  march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular node, not a source/sink
c     irray(i,0) = -i : regular node, producing well, particle capture
c                        but not explicitly specified in sptr macro
c     -100000000< irray(i,0) < -10000000 : regular interior node that is 
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c                = -(i+2000000) : spring node 
c                = -(i+1000000) : well-capture node on extrernal bound
c                   simillar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     -200000000< irray(i,0) <-100000000  : regular cliff node
c      irray(i,0) <-200000000  : OMR cliff node
c     
c.........................................................

         irray0=irray(i,0)
         if(abs(irray0).eq.i.or.irray0.eq.-(i+2000000).or.
     1        irray0.eq.-(i+1000000).or.irray0.eq.-(i+2000).or.
     2        (irray0.gt.-200000000.and.irray0.lt.-100000000) 
     $        )then 
c 3/10/04 s kelkar  3DOMR calling this routine only for regular nodes
c the -i covers the cases of  capture or spring nodes.
            call stream_tube(i,1,icnl_subst,n0,neq,nelm,a_axy,
     $           a_vxy,cord)
c     ******* count # rectangle faces on bdry ***

c########debugging output s kelkar sep 2 05
c       if(cord(i,2).eq.1.25.and.ggg(i,1).gt.0..and.ggg(i,3).lt.0.) then
c          write(ifile,9797)i,cord(i,1),cord(i,3),ggg(i,+1),ggg(i,+3)
c            write(ifile,9797)i,cord(i,1),cord(i,2),cord(i,3),
c     1           ggg(i,-1),ggg(i,-2),ggg(i,-3),ggg(i,+1),
c     2                      ggg(i,+2),ggg(i,+3),s(i)
c       endif
c9797       format(i10,10(g12.5,2x))
c#######################################


            ism=0
            do i3=lower_limit, upper_limit
               if(irray(i,i3).eq.0.and.i3.ne.0) ism=ism+1
            enddo

            if(ism.eq.0) then
c s kelkar Jan 27 05 Internal node..........................
c default, unless otherwise specified in the 'sptr' macro in the 
c keyword 'captur' and 'spring', an internal producing node is 
c taken to be producing well and particle capture routines are
c implemented, where as a producing node on an external boundary
c is taken to be a spring node and the flux is uniformly 
c distributed on the +z (ie ipc=+3) face. Note that even an internal
c node can be specified in the keyword 'spring' as a spring node.
c also NOTE that a boundary node can be specified as a 'well capture'
c node in the sptr macro.
c The default is superseded by specification in the sptr macro 

               if(sk(i).gt.1.e-20) then
                  if(reverse_flow) then
                     irray(i,0) = -i-2000000                        
                  else
                     if(irray(i,0).ne.-i) irray(i,0) = -(i+2000000)
                  endif
               endif
               if(a_axy(nelmdg(i)-1-neq).gt.1.e-20)then
                  if(reverse_flow) then
                     irray(i,0) = -i-2000000                        
                  else
                     if(irray(i,0).ne.-i) irray(i,0) = -(i+2000000)
                  endif
               endif

                  
            else
c node on external boundary               
c if the node is a cliff node, but has a specified boundary
c outflow at it, remove cliff tag and mark as a regular
c boundry node
               if(irray0.lt.-100000000) then
                  if(a_axy(nelmdg(i)-1-neq).gt.0.) then
                     irray(i,0)=-i-2000
                     irray0=-i-2000
                  endif
               endif

               if(irray0.ne.-i.and.irray0.ne.-(i+2000000).and.
     1            irray0.gt. -10000000  ) then
c  for nodes not specified as well-capture or spring nodes 
c  in the sptr macro, partition total to faces on bdry *****
                  do i3=lower_limit, upper_limit
c     if(irray(i,i3).eq.0.and.i3.ne.0) ggg(i,i3)=sk(i)/ism
                     if(irray(i,i3).eq.0.and.i3.ne.0) 
     #                    ggg(i,i3)=a_axy(nelmdg(i)-1-neq)/ism
                  enddo
               else
c node on extenal boundary specified as a well-capture node.
c flag it for special treatment if not a cliff node and 
c not a spring node
                  if(irray0.gt.-100000000.and.irray0.ne.-(i+2000000))
     1                 irray(i,0)=-(i+1000000)
               endif
            end if
            if(irray0.eq.-(i+2000000)) then
c node specified as a spring node. distribute the flux uniformly
c on the +z (ipc=+3) boundary
               ggg(i,+3)=a_axy(nelmdg(i)-1-neq) 
            endif              

c     ************done putting sinks on bdry ********
      
            face1=ddy(i)*ddz(i)
            face2=ddz(i)*ddx(i)
            face3=ddx(i)*ddy(i)
c            facex1=ddy(i)*ddz(i)
c            facey1=ddz(i)*ddx(i)
c            facez1=ddx(i)*ddy(i)
c            facex2=0.
c            facey2=0.
c            facez2=0.
c--------------------------------------------------------------------
c s kelkar 5/13/09 
c zero cross-areas can erroneously result from 
c the ddx,y,z calculated in struct_geom_array(ptrac1) being zero
c when both + and - neighbours on an axis have been removed as 
c porosity <=0 nodes. In this case use area-values stored in sx()
c for OMR nodes, this may give wrong areas
c            if(ism.lt.6) then
c               if(ddx(i).le.0.) then
c                  call area_normal(i,2,facey1)
c                  call area_normal(i,3,facez1)
c               endif
c               if(ddy(i).le.0.) then
c                  call area_normal(i,3,facez2)
c                  call area_normal(i,1,facex1)
c               endif
c               if(ddz(i).le.0.) then
c                  call area_normal(i,1,facex2)
c                  call area_normal(i,2,facey2)
c               endif
c               face1=dmax(facex1,facex2)
c               face2=dmax(facey1,facey2)
c               face3=dmax(facez1,facez2)
c            endif
c--------------------------------------------------------------------

c s kelkar dec 15, 05
c define saturation factors for WTSI problems
            fxr=1.
            fxl=1.
            fyf=1.
            fyb=1.
            fzt=1.
            fzb=1.
            if(ifree.ne.0) then
               if(flag_sat_wtsi.eq.1) then
c upstream saturations in area factors
                  iir=irray(i,+1)
                  if(iir.gt.0) then
                     fxr=0.5*abs(s(i)+s(iir))
                  else
                     fxr=s(i)
                  endif
                  iir=irray(i,-1)
                  if(iir.gt.0) then
                     fxl=0.5*abs(s(i)+s(iir))
                  else
                     fxl=s(i)
                  endif
                  iir=irray(i,+2)
                  if(iir.gt.0) then
                     fyf=0.5*abs(s(i)+s(iir))
                  else
                     fyf=s(i)
                  endif
                  iir=irray(i,-2)
                  if(iir.gt.0) then
                     fyb=0.5*abs(s(i)+s(iir))
                  else
                     fyb=s(i)
                  endif
                  iir=irray(i,-3)
                  if(iir.gt.0) then
                     fzb=0.5*(rlzf(iir)+rlzf(i))
                  else
                     fzb=rlzf(i)
                  endif
               elseif(flag_sat_wtsi.eq.2) then
c upstream saturations in area factors
                  iir=irray(i,+1)
                  if(iir.gt.0) then
                     if(+ggg(i,+1).gt.0.) then
                        fxr=abs(s(iir))
                     else
                        fxr=abs(s(iir))
                     endif
                  else
                     fxr=s(i)
                  endif
                  iir=irray(i,-1)
                  if(iir.gt.0) then
                     if(-ggg(i,-1).gt.0.) then
                        fxl=abs(s(iir))
                     else
                        fxl=abs(s(i))
                     endif
                  else
                     fxl=s(i)
                  endif
                  iir=irray(i,+2)
                  if(iir.gt.0) then
                     if(+ggg(i,+2).gt.0.) then
                        fyf=abs(s(i))
                     else
                        fyf=abs(s(iir))
                     endif
                  else
                     fyf=s(i)
                  endif
                  iir=irray(i,-2)
                  if(iir.gt.0) then
                     if(-ggg(i,-2).gt.0.) then
                        fyb=abs(s(iir))
                     else
                        fyb=abs(s(i))
                     endif
                  else
                     fyb=s(i)
                  endif
c               iir=irray(i,+3)
c               if(iir.gt.0) then
c                  fzt=0.5*abs(s(i)+s(iir))
c               else
c                  fzt=s(i)
c               endif
                  iir=irray(i,-3)
                  if(iir.gt.0) then
                     if(-ggg(i,-3).gt.0.) then
                        fzb=rlzf(iir)
                     else
                        fzb=rlzf(i)
                     endif
                  else
                     fzb=rlzf(i)
                  endif

               endif
            endif
            if(fxr.le.0.)  fxr=1.
            if(fxl.le.0.)  fxl=1.
            if(fyf.le.0.)  fyf=1.
            if(fyb.le.0.)  fyb=1.
c zvd 019-Sep-07 added back in check for 0
c            if(fzt.le.0.)  fzt=1.
c            if(fzb.le.0.)  fzb=1.
            if(fzt.le.0.)  fzt=1.
            if(fzb.le.0.)  fzb=1.
            
c     handle the case of ddx, ddy, ddy, = 0
c     or negative porosity
         
            xface1 = xxx*face1
            xface2 = xxx*face2
            xface3 = xxx*face3
            
            if(xface1.gt.0.) then
               ggg(i, 1)= ggg(i, 1)/(xxx*face1*fxr)
               ggg(i,-1)= ggg(i,-1)/(xxx*face1*fxl)
            else
               ggg(i, 1)= 0.
               ggg(i,-1)= 0.
            end if
            
            if(xface2.gt.0.) then
               ggg(i, 2)= ggg(i, 2)/(xxx*face2*fyf)
               ggg(i,-2)= ggg(i,-2)/(xxx*face2*fyb)
            else
               ggg(i, 2)= 0.
               ggg(i,-2)= 0.
            end if
            
            if(xface3.gt.0..and.icnl.eq.0) then
               ggg(i, 3)= ggg(i, 3)/(xxx*face3*fzt)
               ggg(i,-3)= ggg(i,-3)/(xxx*face3*fzb)
            else
               ggg(i, 3)= 0.
               ggg(i,-3)= 0.
            end if
c now check if it is constant pressure boundary
c and if so,make sure there is nonzero exit velocity to guarantee that
c particles dont enter the cell and then get hung up
            if(ka(i).lt.0) call boundry_vel(i)
            
         elseif(irray0.eq.0 .or. irray0.eq.(-i-1000).or.
     1           irray0.lt.-200000000) then
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c     irray(i,0) < -200000000: OMR cliff node
c s kelkar 3/4/04 3D OMR
c find_v_6lsq does a least squares fit to the mass balance for
c the 6 face velocities needed for pollock interpolation
c vx1,vx2,vy1,vy2,vz1,vz2.
cc***** s kelkar Sep 8 05
c if wtsi option is being used, for nodes with Sw>Sw-min, call find_v
c but for nodes with Sw<Smin, set ggg=0.
            if(ifree.ne.0) then
               if(s(i).ge.smin) then
                  if (ps(i) .gt. 0.d0) call find_v_lsq(i,xxx)
               else
                  do i3=-3,3
                     if(i3.ne.0) ggg(i,i3)=0.
                  enddo
               endif
            else
               if (ps(i) .gt. 0.d0) then
                  flag_por=0
                  do i3=-3,3
                     if(i3.ne.0.and.irray(i,i3).gt.0) 
     1                    flag_por=flag_por+1
                  enddo
                  if(flag_por.gt.0) call find_v_lsq(i,xxx)
               endif
            endif

c for debugging, output mass balance error
!            write(ifile,9511)i,n6lsq,dmass_6lsq,taxy,sqaxy
 9511       format(2(i8,1x),3(2x,g12.5))

         elseif(irray0.eq.-(i+2000000)) then
c node specified as a spring node. distribute the flux uniformly
c on the +z (ipc=+3) boundary
            face3=ddx(i)*ddy(i)
            ggg(i,+3)=a_axy(nelmdg(i)-1-neq)/(xxx*face3) 
            
         endif
         
      enddo

c s kelkar july 9 2004 temporarily fix ggg for corner nodes
!      do i=1,neq      
!         ibou=iboulist(i,7)
!         if(ibou.eq.3) call find_v_3lsq(i,xxx,b3)
!      enddo
c....................................................
c s kelkar feb 8 2005. Modify ggg if velocities at a face from the 
c two sides are pointing into each-other's control volume
c      open(unit=89,file='ggg_debug.dbg')
      call ggg_modify(neq)
c      close(89)
c      stop

c.............................................................


c........kelkar  march 27,2001, for reverse particle tracing
      if (irevers.eq.1) then
         do i=1,neq
            do i3=-3,3
               ggg(i,i3) = -ggg(i,i3)
            enddo
         enddo
      endif
c.....................

c      close(ifile)

c##########debugging ggg
c      ifile = open_file('ggg.dbg','unknown')
c      do i=1,neq
c         write(ifile,9797)i,cord(i,1),cord(i,2),cord(i,3),
c     1           ggg(i,-1),ggg(i,-2),ggg(i,-3),ggg(i,+1),
c     2                      ggg(i,+2),ggg(i,+3),s(i)
c       if(cord(i,2).eq.1.25.and.ggg(i,1).gt.0..and.ggg(i,3).lt.0.) then
c          write(ifile,9797)i,cord(i,1),cord(i,3),ggg(i,+1),ggg(i,+3)
c       endif
c      enddo
c      close(ifile)
c#####################################

      return 

      end 

******************************************************************

      function find_index(i,kb)

      use comai, only : ierr, iptty
      use combi

      implicit none

      integer i,kb,find_index,i1,i2,j

      i1=nelm(i)+1
      i2=nelm(i+1)

      do j=i1,i2
         if(nelm(j).eq.kb) then
            find_index=j
            goto 9191
         endif
      enddo

      write(ierr,*)'Error in find_index in load_struct_flux_array.'
      write(ierr,*)'Cant find index. i,kb=',i,kb
      write(ierr,*)'STOP'
      if (iptty .ne. 0) then
         write(iptty,*)'Error in find_index in load_struct_flux_array.'
         write(iptty,*)'Cant find index. i,kb=',i,kb
         write(iptty,*)'STOP'
      end if
      call exit_ptrac_omr

 9191 continue

      return

      end

c...............................................................

      subroutine exit_ptrac_omr

      stop

      return

      end

c........................................

      subroutine boundry_vel(i)
c s kelkar 1/5/06
c s kelkar  3/11/04  3domr
c calculate velocities between connections between constant pressure
c boundary nodes. FEHM sets these to 0.

      use comai
      use combi
      use comci
      use comdi
      use comsk
      use comsptr

      implicit none

      integer i,i3,i3a,ij

      real*8 perme,visrlp,vel, dist

c      do i3=-3,+3
c         if(i3.ne.0) then
c            ij=irray(i,i3)
c            if(ij.gt.0) then
c               if(ka(ij).lt.0) then
c                  i3a=abs(i3)
c!                  if (ggg(i,i3) .eq. 0 .and. 
cc harmonic averaging of permeability along the connection
c                  if(i3a.eq.1) perme=2.*(pnx(i)*pnx(ij))/
c     &                 (pnx(i)+pnx(ij))
c                  if(i3a.eq.2) perme=2.*(pny(i)*pny(ij))/
c     &                 (pny(i)+pny(ij))
c                  if(i3a.eq.3) perme=2.*(pnz(i)*pnz(ij))/
c     &                 (pnz(i)+pnz(ij))
cc in thermw, dil =relperm*density/viscosity and rolf=density
c                  visrlp=dil(i)/rolf(i)
cc distance along the connection
c                  dist=cord(ij,i3a)-cord(i,i3a)
cc Dacry's law
c                  vel=-perme*visrlp*(phi(ij)-phi(i))/dist
cc save in ggg
c                  ggg(i,i3)=vel
c               endif
c            endif
c         endif
c      enddo

      do i3=-3,+3
         if(i3.ne.0) then
            ij=irray(i,i3)
            if(ij.eq.0) then
               if(abs(ggg(i,i3)).eq.0.) then
c     modify ggg(i,i3) only if it is zero
                  ggg(i,i3) = -ggg(i,-i3)
               endif
            endif
         endif
      enddo
      
      return

      end

c.................................................................
      subroutine ggg_modify(neq)
c check the grid for cases where velocities at two sides of a face
c point into eachother's control volume.

      use comsptr

      implicit none

      integer i,neq,kb,j
      integer, allocatable :: flag_change_ggg(:,:)

      real*8 vi,vkb

      allocate(flag_change_ggg(1:neq,-3:3))
      flag_change_ggg=0

      do i=1,neq
         do j=1,3
            vi=ggg(i,j)
            if(vi.gt.0.) then
               kb=irray(i,j)
               if(kb.eq.0) then
                  if(irray(i,0).eq.0) 
     $                 call ggg_modify_omr(i,j,neq,flag_change_ggg) 
                  if(irray(i,0).eq.-(i+1000)) 
     $                 call ggg_modify_omr_boun(i,j,neq,flag_change_ggg)
               elseif(kb.gt.0) then
                  vkb=-ggg(kb,-j)
                  if(vkb.lt.0.) then
                     if(-ggg(i,-j).gt.0.) then
c case 1
                        if(+ggg(kb,+j).lt.0.) then
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                           if(abs(ggg(i,j)).le.abs(ggg(kb,-j))) then
ccc                              ggg(i,j)=-ggg(kb,-j)
ccc                           else
ccc                              ggg(kb,-j)=-ggg(i,j)
ccc                           endif
                        else
c case 2
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                           ggg(kb,-j)=-ggg(i,+j)
                        endif
                     elseif(-ggg(i,-j).lt.0.) then
c case 3
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                        if(+ggg(kb,+j).lt.0.) then
                           ggg(i,j)=-ggg(kb,-j)
                        else
c case 4
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                           ggg(i,j)=-ggg(i,j)
ccc                           ggg(kb,-j)=-ggg(kb,-j)
                        endif
                     endif
                  endif
               endif
            endif
         enddo

         do j=-3,-1
            vi=-ggg(i,j)
            if(vi.lt.0.) then
               kb=irray(i,j)
               if(kb.eq.0) then
                  if(irray(i,0).eq.0) 
     $                 call ggg_modify_omr(i,j,neq,flag_change_ggg)
                if(irray(i,0).eq.-(i+1000).or.irray(i,0).lt.-200000000) 
     $               call ggg_modify_omr_boun(i,j,neq,flag_change_ggg) 
               elseif(kb.gt.0) then
                  vkb=+ggg(kb,-j)
                  if(vkb.gt.0.) then
                     if(+ggg(i,-j).lt.0.) then
c case 1
                        if(-ggg(kb,j).gt.0.) then
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                           if(abs(ggg(i,j)).le.abs(ggg(kb,-j))) then
ccc                              ggg(i,j)=-ggg(kb,-j)
ccc                           else
ccc                              ggg(kb,-j)=-ggg(i,j)
ccc                           endif
                        else
c case 2
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                           ggg(kb,-j)=-ggg(i,j)
                        endif
                     elseif(+ggg(i,-j).gt.0.) then
c case 3
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                        if(-ggg(kb,j).gt.0.) then
                           ggg(i,j)=-ggg(kb,-j)
                        else
c case 4
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                           ggg(i,j)=-ggg(i,j)
ccc                           ggg(kb,-j)=-ggg(kb,-j)
                        endif
                     endif
                  endif
               endif
            endif
         enddo

      enddo

      deallocate(flag_change_ggg)

      return

      end

c.....................................................................
                
      subroutine ggg_modify_omr(i,j,neq,flag_change_ggg) 
c this an an OMR node not on an external boundary. in this case 
c newnode=0 indicates  non-unique neighbours on the ipc side.
c search all connecting nodes on the ipc side for the nearest
c neighbour to xc,yc,zc

      use combi, only : cord, nelm
      use comsptr
      use comsk

      implicit none

      integer i,j,jsign,i1,i2,jab,kb,ik,itemp(200)
      integer index_temp,aflag
      integer neq,flag_change_ggg(1:neq,-3:3)

      real*8 epsilon,diff,v1bi,v2bi,v1bj,v2bj

      epsilon=1.e-18

      jsign=isign(1,j)
      jab=iabs(j)

c     make a list neighbours of i on the j side
      i1=nelm(i)+1
      i2=nelm(i+1)
      index_temp=0
      do ik=i1,i2
         kb=nelm(ik)
         diff=jsign*(cord(kb,jab)-cord(i,jab))
         if(diff.gt.epsilon) then
c     the node is on the correct side of the face, save it 
            index_temp=index_temp+1
            itemp(index_temp)=kb
         endif
      enddo
      
      v1bi = -jsign*ggg(i,-j)
      if(jsign*v1bi.lt.0.) then

         aflag=0
         do ik=1,index_temp
            kb=itemp(ik)
            v1bj=-jsign*ggg(kb,-j)
            if(jsign*v1bj.lt.0.) then
               v2bj=jsign*ggg(kb,j)
               if(jsign*v2bj.lt.0.) then
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                  if(flag_change_ggg(i,j).eq.0) then
                     ggg(i,j)=jsign*v1bj
                     flag_change_ggg(i,j)=kb
                     aflag=+1
                     goto 99111
                  else
                     call error_ggg_modify_omr
     $                       (i,j,neq,flag_change_ggg)
                  endif
               endif
            endif
         enddo
         do ik=1,index_temp
            kb=itemp(ik)
            v1bj=-jsign*ggg(kb,-j)
            if(jsign*v1bj.lt.0.) then
               v2bj=jsign*ggg(kb,j)
               if(jsign*v2bj.gt.0.) then
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                  if(flag_change_ggg(kb,-j).eq.0) then
ccc                     ggg(kb,-j)=-ggg(i,j)
ccc                     flag_change_ggg(kb,-j)=i
ccc                  else
ccc                     call error_ggg_modify_omr
ccc     $                       (i,j,neq,flag_change_ggg)
ccc                  endif
               endif
            endif
         enddo
99111    continue
         
      elseif(jsign*v1bi.gt.0.) then
         aflag=0
         do ik=1,index_temp
            kb=itemp(ik)
            v1bj=-jsign*ggg(kb,-j)
            if(jsign*v1bj.lt.0.) then
               v2bj=jsign*ggg(kb,j)
               if(jsign*v2bj.lt.0.) then
c out of the 4 velocities, there are two pointing one way
c and two the other way- this ambiguous. Do no tmodify ggg,
c v_newnode will be checked in compute_exit_new  
ccc                  if(abs(ggg(i,j)).le.abs(ggg(kb,-j))) then
ccc                     if(flag_change_ggg(i,j).eq.0) then
ccc                        ggg(i,j)=-ggg(kb,-j)
ccc                        flag_change_ggg(i,j)=kb
ccc                        aflag=+1
ccc                        goto 99112
ccc                     else
ccc                        call error_ggg_modify_omr
ccc     $                       (i,j,neq,flag_change_ggg)
ccc                     endif
ccc                  else
ccc                     if(flag_change_ggg(kb,-j).eq.0) then
ccc                        ggg(kb,-j)=-ggg(i,j)
ccc                        flag_change_ggg(kb,-j)=i
ccc                     else
ccc                        call error_ggg_modify_omr
ccc     $                       (i,j,neq,flag_change_ggg)
ccc                     endif
ccc                  endif
               endif
            endif
         enddo
         if(aflag.eq.0) then
            do ik=1,index_temp
               kb=itemp(ik)
               v1bj=-jsign*ggg(kb,-j)
               if(jsign*v1bj.lt.0.) then
                  v2bj=jsign*ggg(kb,j)
                  if(jsign*v2bj.gt.0.) then
c out of the 4 velocities, 3 are pointing one way and 1 is
c pointing the other way- Fix this one. 
                     if(flag_change_ggg(kb,-j).eq.0) then
                        ggg(kb,-j)=-ggg(i,j)
                        flag_change_ggg(kb,-j)=i
                        goto 99112
                     else
                        call error_ggg_modify_omr
     $                       (i,j,neq,flag_change_ggg)
                     endif
                  endif
               endif
            enddo
99112       continue
         endif
         
      endif
      
      return
      
      end
      
c......................................................................

      subroutine ggg_modify_omr_boun(i,j,neq,flag_change_ggg) 
c this is an OMR node on an external boundary. check its boundary
c faces, stired in iboulist(), to see if the ipc plane is a
c boundary plane. If it is not, find the new nearest node.

      use comsptr
      use comsk

      implicit none

      integer neq,flag_change_ggg(1:neq,-3:3)
      integer i,j,ibou,i1

      ibou=iboulist(i,7)
      do i1=1,ibou
         if(j.eq.iboulist(i,i1)) then
            goto 9999
         endif
      enddo

      call ggg_modify_omr(i,j,neq,flag_change_ggg)
      
 9999 return
      
      end

c...........................................................................

      subroutine error_ggg_modify_omr (i,j,neq,flag_change_ggg)

      use comai, only : ierr, iptty 
      implicit none

      integer neq,flag_change_ggg(1:neq,-3:3)
      integer i,j

      write(ierr,*)'Error in ggg_modify_omr. ggg(i,j) has already'
      write(ierr,*)' been changed wrt neighbour kb.'
      write(ierr,*)'i,j,neq,kb=',i,j,neq,flag_change_ggg(i,j)
      write(ierr,*)'STOP'
      if (iptty .ne. 0) then
         write(iptty,*)'Error in ggg_modify_omr. ggg(i,j) has already'
         write(iptty,*)' been changed wrt neighbour kb.'
         write(iptty,*)'i,j,neq,kb=',i,j,neq,flag_change_ggg(i,j)
         write(iptty,*)'STOP'
      end if

      stop

      return

      end

c.......................................................................

      subroutine area_normal(i,idir,face)

c s kelkar 5/13/09 
c Zero cross-areas can erroneously result from 
c the ddx,y,z calculated in struct_geom_array(ptrac1) being zero
c when both + and - neighbours on an axis have been removed as 
c porosity <=0 nodes. In this case use area-values stored in sx()
c for OMR nodes, this may give wrong areas
      
      use comai, only : ierr, iptty, neq
      use combi, only : cord, isox, isoy, isoz, istrw, nelm, sx
      use comsptr, only: irray

      implicit none

      integer i,i1,i2,ii1,ii2,j,iwskj, idir, k

      real*8 face,face1,face2

      face1=0.
      face2=0.
      i1=irray(i,-idir)
      i2=irray(i,+idir)
      ii1=nelm(i)+1
      ii2=nelm(i+1)
c     area normal to axis
      if(i1.gt.0) then
         do k=ii1,ii2 
            if(nelm(k).eq.i1) goto 9991
         enddo 
         if (iptty .ne. 0) then
            write(iptty,*)'Error in area_normal. Stop.'
            write(iptty,*)' Cant find nelm(k)=i1'
            write(iptty,*)'i,ii1,ii2,i1,idir=',i,ii1,ii2,i1,idir
         end if
         stop
 9991    iwskj=istrw(k-neq-1)
         face1=sx(iwskj,isox)+sx(iwskj,isoy)+sx(iwskj,isoz)
         face1=abs(face1*(cord(i,-idir)-cord(i1,-idir)))
      endif
      if(i2.gt.0) then
         do k=ii1,ii2 
            if(nelm(k).eq.i2) goto 9992
         enddo 
         if (iptty .ne. 0) then
            write(iptty,*)'Error in area_normal. Stop.'
            write(iptty,*)' Cant find nelm(k)=i2'
            write(iptty,*)'i,ii1,ii2,i2,idir=',i,ii1,ii2,i2,idir
         end if
         stop
 9992    iwskj=istrw(k-neq-1)
         face2=sx(iwskj,isox)+sx(iwskj,isoy)+sx(iwskj,isoz)
         face2=abs(face2*(cord(i,+idir)-cord(i1,+idir)))
      endif
      face=max(face1,face2)
      
      end subroutine area_normal
c---------------------------------------------------------------------
