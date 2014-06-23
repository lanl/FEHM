      subroutine compute_exit_new(np1,ep5,ep,small_slope,
     $     idum_count,time_step_flag,inp_adv,
     1     vx_zheng,vy_zheng,vz_zheng)
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
!D1 Compute new location of particle.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.22
!D2 
!D2 Initial implementation: 31-Mar-04, Programmer: S Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/compute_exit_new_a  $
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
c modified compute_exit_struct, S kelkar, March 25 04

      use comai
      use combi
      use comflow
      use comdi
      use comsptr
      use comsk
      
      implicit none
      
      integer i,inp1,icapture,np1,newnode,ipc,ipcsave,oldinp
      integer iplane(3),nplanes,ipsign,ipab,nplanes2
      integer idum_count,j,ibou,ifix_flag,time_step_flag
      integer flag_x1,flag_x2,flag_y1,flag_y2,flag_z1,flag_z2
      integer flag_slow,flag_fixggg,ipcsign,inp_adv
      integer flag_divd,current_flag,current_model
      integer indexrw,inp_exit_zheng
      integer irray0,newray0
      
      real*8 ep5,ep,small_slope,x2dim,y2dim,z2dim
      real*8 dx,dy,dz,x,y,z,xx,yy,zz,pi
      real*8 vx_zheng,vy_zheng,vz_zheng,rw,rwa
      real*8 xw1,xw2,yw1,yw2,zw1,zw2
      real*8 epsilon,xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc,cl,ccm,cn,param
      real*8 xpnew, ypnew, zpnew,epsk,vomrx,vomry,vomrz
      real*8 aomrx,aomry,aomrz,axydg,epomr,epvel
      real*8 abv_current,abv_new,abrat,eprat,dtc,dtnew
      real*8 vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,vexit,vother,venter
      real*8 thickness_factor,v_newnode
      real*8 vel_divd

      epsilon=1.e-10
      epsk=1.e-20
      epomr=1.e-22
      eprat=1.e-8
      epvel=1.e-8
      pi=3.14159

c if  particle has made a random jump on the
c  previous call to random_walk, oldnode and oldnode2 
c were reset to ijkv(np1) in new_part _loc_struct
c     to avoid confusing this with a bouncing particle
      
      inp1=ijkv(np1)
      oldinp=oldnode(np1)
      irray0=irray(inp1,0)
      inp_adv=inp1
	
c     at this point ipc is referenced to the node inp1, and denotes
c     the face thru which the particle entered the cc of inp1
      ipcsave=ipc_save(np1)
      ipc=ipcsave
      
      if(inp1.gt.0) then

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c     set flag to add the grad-dispersion term here. 
c     gradd was calculated in 
c     subroutine dispersion_divergence_omr which was called 
c     from ptrac3 
c     The sign of itensor is the flag that says if we are to
c     use the gradd term or not
            
         flag_divd=0
         current_model= itrc(inp1) 
         current_flag = tprpflag(current_model)
         if(current_flag.eq.2.or.current_flag.eq.4) then
            if(omr_flag) then
               if(itensor.gt.0) then
                  flag_divd=1
               endif
            endif
         endif               
                        
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

         rw=0.
         indexrw=0
         if(irray0.gt.-100000000) indexrw=-(irray0+10000000)
         if(irray(inp1,0).eq.-(inp1+1000000).or.
     1        irray(inp1,0).eq.-inp1.or.indexrw.gt.0) then
            if(indexrw.gt.0.and.indexrw.le.neq) then
               rw=well_radius(indexrw)
            else
c     rw = 0.1*min(ddx(inp1),ddy(inp1))
               rw=0.1
            endif
         endif
         
         if(rw.gt.0.) then
            call zheng2(inp1,np1,rw,vx_zheng,vy_zheng,vz_zheng)
         else
            
            ic_new(np1) = 0
            
c     initialize the particle position and velocities and derivatives
c     note that this can be a complicated matter in some omr cc's
            call initialize_particle(flag_divd,inp1,oldinp,np1,
     $           ipcsave,xp1,yp1,
     1           zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2           vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
            
c...  june 10 04, s kelkar..........................................
c     next is a quick fix to deal with andy's recharge nodes on the top 
c     surface, sometimes they get spurious upstream velocities
c     first check if the current node is an OMR boundary node,
c     or an internal, non-OMR  node, or a boundary-nonOMR node
c     Then check if it is an injection or 
c     a recharge node and if the node has a missing neighbour in +z
            irray0=irray(inp1,0)
            if(irray0.eq.-(inp1+1000).or.irray0.eq.inp1.or.
     &           irray0.eq.-(inp1+2000).or.irray0.lt.-100000000) then

               axydg=a_axy(nelmdg(inp1)-1-neq)
               if(abs(axydg).gt.epomr) then
                  if(irray(inp1,3).eq.0) then
                     call fixggg_recharge(flag_divd,inp1,oldnode(np1),
     $                    np1,ipc,
     $                    xp1,yp1,zp1,vomrx,vomry,vomrz,aomrx,aomry,
     2                    aomrz,vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
                  endif
               endif
            endif
c.................................................................
            
c     calculate the timestep that puts the particle at the exit 
c     boundary from the cc. If the input sptr file has a -ve
c     value of dtmx, then the abs of that number is time step. 
            if(time_step_flag.ne.1) then
c               call newdt3(inp1,np1,xp1,yp1,zp1,vomrx,
c     1              vomry,vomrz,aomrx,aomry,aomrz,
c     2              vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,dtnew)
c s kelkar april 12 2005, to calculate time to cross the whole
c cell and use a currant_fraction multiplier.
c In case of random walk, 
c also include the time step based on Div(D) term 
               call newdt_courant(inp1,np1,aomrx,aomry,aomrz,
     1              vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,
     2              vomrx,vomry,vomrz,xp1,yp1,zp1,
     3              dtnew)
            else
               dtnew=abs(dtmx)
            endif
            
c....................................................
c     calculate preliminary new particle location xp2,yp2,zp2
            
            call newlocation(inp1,np1,small_slope,ep5,xp1,yp1,zp1,
     1           dtnew,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2           vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,
     3           xp2,yp2,zp2)
            
            
c.... check for well-capture
            icapture=0
            if(rw.lt.0.) then
c     do a particle capture without invoking the semi-analytical 
c     modification to the flow field
               rwa=-rw
               xw1=xp1-0.5*ddx(inp1)
               yw1=yp1-0.5*ddy(inp1)
               xw2=xp2-0.5*ddx(inp1)
               yw2=yp2-0.5*ddy(inp1)
               call check_well_capture(inp1,xw1,yw1,xw2,yw2,
     1              icapture,vx_zheng,vy_zheng,vz_zheng,rwa)
            endif
            if(icapture.eq.1) then
               ipc=-10
               inp_exit_zheng=-inp1
               call update_exit(inp1,np1,ipc,inp_exit_zheng,dtc,
     1              xc,yc,zc)
               call seek_btczone(np1, inp1)
c     call remove_particle(np1,5,exit_node(np1),ipc)
               exit_node(np1)=inp1
               call capture_particle(np1)
            else
               dtc=dtnew
               xc=xp2
               yc=yp2
               zc=zp2
               
c     now find where the new location is........................
               
c     s kelkar nov 24 04..............................     
c     nplanes=1
c     nplanes2=nplanes
c     call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
c     if the point lands on a face,
c     fix xc,yc,zc, to be slightly inside the new control volume 
c     call shift_xc_inside(inp1,newnode,np1,ipc,xc,yc,zc)        
c     c s kelkar 11/18/04 ipc is calculated in newdt3 and new location 
c     c is calculated in newlocation assuming that the particle crosses
c     c the entire control voulme.
c     c  until RANDOM_WALK is activated, no need
c     c to use the more general routines find_exit_planes and 
c     c find_exit_point.
c..............................................................
               
               call find_exit_planes(inp1,np1,xp2,yp2,zp2,nplanes,
     1              iplane)
               nplanes2=nplanes
               
               if(nplanes.eq.0) then
c     if nplanes=0 then particle still in the same cc
                  
                  count_steps(np1)=count_steps(np1)+1
                  if(count_steps(np1).gt.count_steps_max) then
                     call remove_particle(np1,1,inp1,ipc)
                  else
                     call update_samenode(inp1,np1,ipc,newnode,
     $                    dtc,xp2,yp2,zp2)
                  endif
                  
               else
c     find new control volume. nplanes2 is used later to check bouncing
c     ipc gets changed to the face at which the particle will exit
c     the cc of inp1, it is still wrt the current node inp1
                  
                  call find_exit_point(inp1,np1,nplanes,iplane,xp1,
     1                 yp1,zp1,xp2,yp2,zp2,
     2                 vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3                 ipc,xc,yc,zc,dtc,newnode,vexit)
	

! If ioutt has been set to 1 to signal end of timestep for this particle
! but the time was reduced because we exited a control volume, reset
! ioutt so calculations will continue
                  if (dtc .lt. dtnew .and. ioutt(np1) .eq. 1) 
     &                 ioutt(np1) = 0
c  if newnode is wtsi (Water table) node 
c   then reset z-velocities and recalculate zp2
c  to be within dzw of newnode. I
                  if(newnode.gt.0.and.ifree.ne.0) then
                     if(izone_free_nodes(newnode).gt.1) then
c      partially saturated water table node
                        if(a_axy(nelmdg(newnode)-1-neq).le.0.) then
c	 not a sink node

c                        if(abs(ipc).ne.3) then
                        call wtsi_newnode(inp1,np1,newnode,
     1                       xc,yc,zc,vomrz,aomrz,vz1b,vz2b)
                        goto 99999
c                        endif
					endif
                     endif
                  endif
                  
c     check if a spring node, if the particle has exited from top, 
c     remove it from the model. Note spring nodes can be internal
                  if(irray(inp1,0).eq.-(inp1+2000000) ) then
                     if(ipc.eq.+3) then
                        call exit_spring(inp1,np1,ipc,newnode,
     1                       dtc,xc,yc,zc)
                        goto 9999
                     endif
                  endif
c>>>
                     if(omr_flag) then
                     if(newnode.gt.0) then
                        
                        if(irray(newnode,0).ne.-newnode) then

c     1                    . and.irray(newnode,0).gt.-(10000000)) then
c     if there is a legit newnode, and it is not a capture node, then
c..   s kelkar  june 25 04 .. if the particle jumps into a slow moving
c     node, remove it from the computations. But the tollarance, eprat
c     is set quite low, because this situation can sometimes be legit 
c     in regular nodes. A special chack with higher tolarance, epvel,
c     is done later for OMR nodes
                           
c     define vexit, the velocity at the exit face
                           if(ipc.eq.+1) then
                              vexit=vx2b
                           elseif(ipc.eq.-1) then
                              vexit=vx1b
                           elseif(ipc.eq.+2) then
                              vexit=vy2b
                           elseif(ipc.eq.-2) then
                              vexit=vy1b
                           elseif(ipc.eq.+3) then
                              vexit=vz2b
                           elseif(ipc.eq.-3) then
                              vexit=vz1b
                           endif
                           
                           flag_fixggg=0                  
                           call slow_particle(flag_divd,inp1,ipc,
     1                          np1,newnode,vexit,eprat,abrat,
     2                          flag_slow)
                           
                           if (flag_slow.eq.1) then
                              flag_fixggg=+1                  
                              call fixggg_3(flag_divd,inp1,oldnode(np1),
     $                             np1,ipc, xp1,yp1, zp1,vomrx,
     $                             vomry,vomrz,aomrx,aomry, aomrz,
     2                             vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
c     s kelkar.......april 8 04  check for bouncing between two cells
                              
                           elseif((newnode.eq.oldnode(np1)).and.
     &                             (oldnode(np1).ne.inp1)) then
                              flag_fixggg=+1                  
                              call fixggg_1(flag_divd,inp1,oldnode(np1),
     $                             np1,ipc,
     1                             vomrx,vomry,vomrz,aomrx,aomry,aomrz)      
                              idum_count=idum_count+1
                           elseif((newnode.eq.oldnode2(np1)).and.
     &                             (oldnode2(np1).ne.inp1)) then
                              flag_fixggg=+1                  
                              call fixggg_1(flag_divd,inp1,
     $                             oldnode2(np1),np1,ipc, vomrx,
     1                            vomry,vomrz,aomrx,aomry,aomrz) 
                              idum_count=idum_count+1
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     s kelkar march 10 05 ggg has been preprossed to remove opposing 
c     velocities on faces for tha cases where 3 of the 4 velocities
c     for the pair cc's in the abs(ipc) direction are pointing the same
c     way. But there are many instances where 2 velocities point one way 
c     and two point the other way- in these cases zero out the velocity 
c     within inp1 pointing to newnode and recalculate the trajectory
                              
                           else

                              irray0=irray(inp1,0)
                              newray0=irray(newnode,0)
                              if((newray0.ne.newnode.and.
     $                            newray0.ne.-(newnode+2000).and.
     $               (newray0.lt.-200000000.or.newray0.gt.-100000000)).
     $                             or.(irray0.ne.inp1.and.
     $                             irray0.ne.-(inp1+2000).and.
     $                 (irray0.lt.-200000000.or.irray0.gt.-100000000))
     $                            ) then
                                 ipcsign=isign(1,ipc)
c^^^  
                              call sub_vel_divd(flag_divd,newnode,-ipc,
     $                                vel_divd)
c^^^  
                                 v_newnode=vel_divd
                                 if((v_newnode*vexit).le.0.)  then
                                    flag_fixggg=+1                  
                                    call fixggg_3(flag_divd,inp1,
     $                                   oldnode(np1),np1,ipc,
     $                                   xp1,yp1,zp1,vomrx,vomry,vomrz,
     2                                aomrx,aomry,aomrz,vx1b,vy1b,vz1b,
     3                                   vx2b,vy2b,vz2b)
                                    idum_count=idum_count+1
                                 endif
                              endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                           endif
                           
                           if(flag_fixggg.eq.1) then
c     re-calculate the timestep that puts the particle at the exit 
c     boundary from the cc. If the input sptr file has a -ve
c     value of dtmx, then the abs of that number is time step. 
                              if(time_step_flag.ne.1) then
c     s kelkar april 12 2005, to calculate time to cross the whole
c     cell and use a currant_fraction multiplier
                                call newdt_courant(inp1,np1,aomrx,aomry,
     1                                aomrz,vx1b,vy1b,vz1b,vx2b,vy2b,
     2                                vz2b,
     3                                vomrx,vomry,vomrz,xp1,yp1,zp1,
     4                                dtnew)
                              else
                                 dtnew=abs(dtmx)
                              endif
                              
c     re-calculate new location
                              call newlocation(inp1,np1,small_slope,ep5,
     1                             xp1,yp1,zp1,dtnew,
     2                             vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3                             vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,
     4                             xp2,yp2,zp2)
                              dtc=dtnew
                              xc=xp2
                              yc=yp2
                              zc=zp2
                              call find_exit_planes(inp1,np1,
     $                             xp2,yp2,zp2,nplanes2,iplane)
                              if(nplanes2.eq.0 ) then
c     c     if nplanes2=0 then modified position still in the same cc
                                 count_steps(np1)=count_steps(np1)+1
                                 if(count_steps(np1).gt.count_steps_max)
     &                                then
                                    call remove_particle(np1,1,inp1,ipc)
                                 else
                                    call update_samenode(inp1,np1,ipc,
     $                                   newnode,dtc,xp2,yp2,zp2)
                                 endif
                              else
                                 
c     c     new control volume
                                 call find_exit_point(inp1,np1,nplanes2,
     $                                iplane,xp1,yp1,zp1,xp2,yp2,zp2,
     2                              vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3                                ipc,xc,yc,zc,dtc,newnode,vexit)
! If ioutt has been set to 1 to signal end of timestep for this particle
! but the time was reduced because we exited a control volume, reset
! ioutt so calculations will continue
                                 if (dtc .lt. dtnew .and. 
     &                                ioutt(np1) .eq. 1) ioutt(np1) = 0
                                 
c     check if a spring node, if the particle has exited from top, 
c     remove it from the model. Note spring nodes can be internal
                                if(irray(inp1,0).eq.-(inp1+2000000))then
                                    if(ipc.eq.+3) then
                                       call exit_spring(inp1,np1,ipc,
     1                                      newnode,dtc,xc,yc,zc)
                                       goto 9999
                                    endif
                                 endif
                                 
                              endif
                           endif
                        endif
                     endif
                  endif
c>>>
                  
               endif

99999          continue

               if(nplanes2.gt.0) then            
                  if(newnode.gt.0) then
                     
c     found a legitimate new neighbour, proceed with updating
c     location and indices
                     
                     count_steps(np1)=0
                     idum_count=0
                     call update_newnode(inp1,np1,ipc,newnode,
     $                    dtc,xc,yc,zc)                     
                     call seek_btczone(np1, newnode)
                  elseif(newnode .le. 0) then
c     if newnode=0, particle has exited the model domain
                     call update_exit(inp1,np1,ipc,newnode,
     $                    dtc,xc,yc,zc)
                     if(newnode.lt.0) then
c particle has entered a -ve or 0 porosity node. Not allowed
c Flag and remove particle
                        istop(np1) = 10
                        write(ierr, 9998) inp1,np1,ipc,newnode
                        if (iptty .ne. 0) 
     &                       write(iptty, 9998) inp1,np1,ipc,newnode
                     end if
                  endif         
               endif
               
            endif
            
         endif
         
 9999    continue
 9998    format ('STOP In compute_exit_new.', 
     &        ' Particle has entered a -ve porosity node', /,
     &        'inp1, np1, ipc, newnode =', 3(1x, i8), /)
      else
         write(ierr,*)'stop. ijkv(np1) = 0 on entry to compute_exit_new'
        write(iptty,*)'stop. ijkv(np1) = 0 on entry to compute_exit_new'
c     stop
         call remove_particle(np1,3,inp1,ipc)
      endif
      
      return
      
      end subroutine compute_exit_new
      
c..........................................................

      subroutine get_velocities(inp1,np1,xp1,yp1,zp1,
     1           aomrx,aomry,aomrz,vx1b,vy1b,vz1b,vomrx,vomry,vomrz)

      use comsptr
      use comsk

      implicit none

      integer inp1,np1

      real*8 xp1,yp1,zp1,aomrx,aomry,aomrz,vomrx,vomry,vomrz
      real*8 vx1b,vy1b,vz1b

      vomrx=aomrx*xp1+vx1b
      vomry=aomry*yp1+vy1b
      vomrz=aomrz*zp1+vz1b

      return

      end subroutine get_velocities

c...............................................................


      subroutine newlocation(inp1,np1,small_slope,ep5,xp1,yp1,zp1,
     1     dtnew,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,
     4     xp2,yp2,zp2)

      use comsk
      use comsptr

      implicit none

      integer np1,inp1
      integer flag_x1,flag_x2,flag_y1,flag_y2,flag_z1,flag_z2

      real*8 ep5,xp1,yp1,zp1,xp2,yp2,zp2,vomrx,vomry,vomrz,dtnew
      real*8 vx2b, vx1b,vy2b,vy1b,vz2b,vz1b,ep4
      real*8 small_slope,aomrx,aomry,aomrz,ddxinp,ddyinp,ddzinp

      ep5=1.e-10
      ep4=1.e-4
      xp2=xp1
      yp2=yp1
      zp2=zp1
      ddxinp=ddx(inp1)
      ddyinp=ddy(inp1)
      ddzinp=ddz(inp1)
      
c recall that at this point ipc_save has the entry face for
c the cc of inp1 wrt inp1.

      if(vomrx.gt.0.and.ipc_save(np1).ne.+1) then
c in case of o  velocity on the right face, make
c sure particle doesnt get pushed into it, or into
c next cc by numerical error
         if(vx2b.le.0.) then
            if(abs((ddxinp-xp1)/ddxinp).gt.ep4) then
               call x_displacement(dtnew,xp1,vomrx,aomrx,
     $              small_slope,xp2)
            endif
            if(((ddxinp-xp2)/ddxinp).lt.ep4) 
     $           xp2=ddxinp*(1.-ep4)
         else
            call x_displacement(dtnew,xp1,vomrx,aomrx,
     $           small_slope,xp2)
         endif
      elseif(vomrx.lt.-0.and.ipc_save(np1).ne.-1)then
c in case of o  velocity on the left face, make
c sure particle doesnt get pushed into it, or into the 
c the next cc  by numerical error
         if(vx1b.ge.0.) then
            if(abs(xp1/ddxinp).gt.ep4) then
               call x_displacement(dtnew,xp1,vomrx,aomrx,
     $              small_slope,xp2)
            endif
            if((xp2/ddxinp).lt.ep4) xp2=ddxinp*ep4
         else
            call x_displacement(dtnew,xp1,vomrx,aomrx,
     $           small_slope,xp2)
         endif
      endif
      xp2=(1.+ep5)*xp2-ep5*xp1
      
c     Set y2 and y3
      if(vomry.gt.0.and.ipc_save(np1).ne.+2) then
c in case of o  velocity on the back face, make
c sure particle doesnt get pushed into it, or into next
c cc by numerical error
         if(vy2b.le.0.) then
            if(abs((ddyinp-yp1)/ddyinp).gt.ep4) then
               call y_displacement(dtnew,yp1,vomry,aomry,
     $              small_slope,yp2)
            endif
            if(((ddyinp-yp2)/ddyinp).lt.ep4) 
     $           yp2=ddyinp*(1.-ep4)
         else
            call y_displacement(dtnew,yp1,vomry,aomry,
     $           small_slope,yp2)
         endif
      elseif(vomry.lt.-0.and.ipc_save(np1).ne.-2)then
c in case of o  velocity on the front face, make
c sure particle doesnt get pushed into it, or into
c the next cc  by numerical error
         if(vy1b.ge.0.) then
            if(abs(yp1/ddyinp).gt.ep4) then
               call y_displacement(dtnew,yp1,vomry,aomry,
     $              small_slope,yp2)
            endif
            if((yp2/ddyinp).lt.ep4) yp2=ddyinp*ep4
         else
            call y_displacement(dtnew,yp1,vomry,aomry,
     $           small_slope,yp2)
         endif
      endif
      yp2=(1.+ep5)*yp2-ep5*yp1
      
c     Set z2 and z3
      if(vomrz.gt.0.and.ipc_save(np1).ne.+3) then
c in case of o  velocity on the top face, make
c sure particle doesnt get pushed into it by numerical error
         if(vz2b.le.0.) then
            if(abs((ddzinp-zp1)/ddzinp).gt.ep4) then
               call z_displacement(dtnew,zp1,vomrz,aomrz,
     $              small_slope,zp2)
            endif
            if(((ddzinp-zp2)/ddzinp).lt.ep4)
     $           zp2=ddzinp*(1.-ep4)
         else
            call z_displacement(dtnew,zp1,vomrz,aomrz,
     $           small_slope,zp2)
         endif
      elseif(vomrz.lt.-0.and.ipc_save(np1).ne.-3)then
c in case of o  velocity on the bottom face, make
c sure particle doesnt get pushed into it by numerical error
         if(vz1b.ge.0.) then
            if(abs(zp1/ddzinp).gt.ep4) then
               call z_displacement(dtnew,zp1,vomrz,aomrz,
     $              small_slope,zp2)
            endif
            if((zp2/ddzinp).lt.ep4) zp2=ddzinp*ep4
         else
            call z_displacement(dtnew,zp1,vomrz,aomrz,
     $           small_slope,zp2)
         endif
      endif
      zp2=(1.+ep5)*zp2-ep5*zp1
      
      return
       
      end subroutine newlocation
      
c.......................................................................

      subroutine x_displacement(dtnew,xp1,vomrx,aomrx,small_slope,xp2)

      implicit none

      real*8 dtnew,xp1,vomrx,aomrx,xp2,small_slope

      if(abs(aomrx*dtnew).le.small_slope) then
         xp2=xp1+vomrx*dtnew
      else
         xp2=xp1+(dexp(aomrx*dtnew)-1.)*vomrx/aomrx
      end if

      return

      end subroutine x_displacement

c..........................................................................


      subroutine y_displacement(dtnew,yp1,vomry,aomry,small_slope,yp2)

      implicit none

      real*8 dtnew,yp1,vomry,aomry,yp2,small_slope

      if(abs(aomry*dtnew).le.small_slope) then
         yp2=yp1+vomry*dtnew
      else
         yp2=yp1+(dexp(aomry*dtnew)-1.)*vomry/aomry
      end if

      return

      end subroutine y_displacement 

c.................................................................

      subroutine z_displacement(dtnew,zp1,vomrz,aomrz,small_slope,zp2)

      implicit none

      real*8 dtnew,zp1,vomrz,aomrz,zp2,small_slope

      if(abs(aomrz*dtnew).le.small_slope) then
         zp2=zp1+vomrz*dtnew
      else
         zp2=zp1+(dexp(aomrz*dtnew)-1.)*vomrz/aomrz
      end if

      return

      end subroutine z_displacement  

c.........................................................


             
      subroutine fixggg_recharge(flag_divd,inp1,inp1old,np1,ipc,
     $     xp1,yp1,zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
      
      use comsk
      use comsptr

      implicit none

      integer inp1,np1,inp1old,ipc
      integer flag_divd

      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vvinp,vvold,xp1,yp1,zp1
      real*8 vx1b,vy1b,vz1b,vx2b,vy2b,vz2b
      real*8 vel_divd
      
c ipc_save(np1)=face particle entered the cc of inp1, wrt the inp1
      if(ipc_save(np1).eq.+1) then
         vvinp=vx2b
         call sub_vel_divd(flag_divd,inp1old,-1,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vx2b=vvold
            aomrx=(vx2b-vx1b)/ddx(inp1)
            vomrx=vx1b+aomrx*xp1
         endif
      elseif(ipc_save(np1).eq.-1) then
         vvinp=vx1b
         call sub_vel_divd(flag_divd,inp1old,+1,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vx1b=vvold
            aomrx=(vx2b-vx1b)/ddx(inp1)
            vomrx=vx1b+aomrx*xp1
         endif
      elseif(ipc_save(np1).eq.+2) then
         vvinp=vy2b
         call sub_vel_divd(flag_divd,inp1old,-2,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vy2b=vvold
            aomry=(vy2b-vy1b)/ddy(inp1)
            vomry=vy1b+aomry*yp1
         endif
      elseif(ipc_save(np1).eq.-2) then
         vvinp=vy1b
         call sub_vel_divd(flag_divd,inp1old,+2,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vy1b=vvold
            aomry=(vy2b-vy1b)/ddy(inp1)
            vomry=vy1b+aomry*yp1
         endif
      elseif(ipc_save(np1).eq.+3) then
         vvinp=vz2b
         call sub_vel_divd(flag_divd,inp1old,-3,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vz2b=vvold
            aomrz=(vz2b-vz1b)/ddz(inp1)
            vomrz=vz1b+aomrz*zp1
         endif
      elseif(ipc_save(np1).eq.-3) then
         vvinp=vz1b
         call sub_vel_divd(flag_divd,inp1old,+3,vel_divd)
         vvold=vel_divd
         if((vvold*vvinp).lt.0.) then
            vz1b=vvold
            aomrz=(vz2b-vz1b)/ddz(inp1)
            vomrz=vz1b+aomrz*zp1
         endif
      endif

      return

      end
c.......................................................................

      subroutine fixggg_1(flag_divd,inp1,inp1old,np1,ipc,
     $     vomrx,vomry,vomrz,aomrx,aomry,aomrz)
      
      use comsk
      use comsptr

      implicit none

      integer inp1,np1,ipc,inp1old,flag_divd

      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vvinp,vvold, vvopp

c     in case the velocity has opposite signs on the same face
c     on the two different sides, use the value from the upstream
c     face, ie value out of the control volume of the oldnode.
c     NOTE: ipc is wrt the current node, and the face has the 
c     opposite sign wrt oldnode

      if(ipc.eq.+1) then
         call sub_vel_divd(flag_divd,inp1,+1,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,-1,vvold) 
         call sub_vel_divd(flag_divd,inp1,-1,vvopp) 
c         vvinp=+ggg(inp1,+1)
c         vvold=-ggg(inp1old,-1)
c         if(vvinp*(-ggg(inp1,-1)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if x velocity on both the left and right face are pointing
c     to right, zero out the normal component of the velocity
            aomrx=0.
            vomrx=0.
         else
c     but if the velocity on the left face is pointing left
c     , ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomrx=vvold
            aomrx=(vvold-vvopp)/ddx(inp1)
         endif

      elseif(ipc.eq.-1) then
         call sub_vel_divd(flag_divd,inp1,-1,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,+1,vvold) 
         call sub_vel_divd(flag_divd,inp1,+1,vvopp) 
c         vvinp=-ggg(inp1,-1)
c         vvold=+ggg(inp1old,+1)
c         if(vvinp*(+ggg(inp1,+1)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if x velocity on both the left and right face are pointing
c     left, zero out the normal component of the velocity
            aomrx=0.
            vomrx=0.
         else
c     but if the velocity on the right face is pointing to right
c     , ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomrx=vvold
            aomrx=(vvopp-vvold)/ddx(inp1)
         endif

      elseif(ipc.eq.+2) then
         call sub_vel_divd(flag_divd,inp1,+2,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,-2,vvold) 
         call sub_vel_divd(flag_divd,inp1,-2,vvopp) 
c         vvinp=+ggg(inp1,+2)
c         vvold=-ggg(inp1old,-2)
c         if(vvinp*(-ggg(inp1,-2)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if y velocity on both the front and back face are pointing
c     the wrong way, zero out the normal component of the velocity
            aomry=0.
            vomry=0.
         else
c     but if the velocity on the back face is pointing the correct
c     way, ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomry=vvold
            aomry=(vvold-vvopp)/ddy(inp1)
         endif

      elseif(ipc.eq.-2) then
         call sub_vel_divd(flag_divd,inp1,-2,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,+2,vvold) 
         call sub_vel_divd(flag_divd,inp1,+2,vvopp) 
c         vvinp=-ggg(inp1,-2)
c         vvold=+ggg(inp1old,+2)
c         if(vvinp*(+ggg(inp1,+2)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if y velocity on both the front and back face are pointing
c     the wrong way, zero out the normal component of the velocity
            aomry=0.
            vomry=0.
         else
c     but if the velocity on the front face is pointing the correct
c     way, ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomry=vvold
            aomry=(vvopp-vvold)/ddy(inp1)
         endif

      elseif(ipc.eq.+3) then
         call sub_vel_divd(flag_divd,inp1,+3,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,-3,vvold) 
         call sub_vel_divd(flag_divd,inp1,-3,vvopp) 
c         vvinp=+ggg(inp1,+3)
c         vvold=-ggg(inp1old,-3)
c         if(vvinp*(-ggg(inp1,-3)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if z velocity on both the upper and lower face are pointing
c     the wrong way, zero out the normal component of the velocity
            aomrz=0.
            vomrz=0.
         else
c     but if the velocity on the lower face is pointing the correct
c     way, ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomrz=vvold
            aomrz=(vvold-vvopp)/ddz(inp1)
         endif

      elseif(ipc.eq.-3) then
         call sub_vel_divd(flag_divd,inp1,-3,vvinp) 
         call sub_vel_divd(flag_divd,inp1old,+3,vvold) 
         call sub_vel_divd(flag_divd,inp1,+3,vvopp) 
c         vvinp=-ggg(inp1,-3)
c         vvold=+ggg(inp1old,+3)
c         if(vvinp*(+ggg(inp1,+3)).gt.0.) then
         if(vvinp*vvopp.gt.0.) then
c     if z velocity on both the upper and lower face are pointing
c     the wrong way, zero out the normal component of the velocity
            aomrz=0.
            vomrz=0.
         else
c     but if the velocity on the upper face is pointing the correct
c     way, ie there is a flow devide, use the velocity and accelaration
c     from the oldnode
            vomrz=vvold
            aomrz=(vvopp-vvold)/ddz(inp1)
         endif
      endif

      return

      end subroutine fixggg_1
c.......................................................................

      subroutine fixggg_2(flag_divd,inp1,np1,ipc,
     $     vomrx,vomry,vomrz,aomrx,aomry,aomrz)
      
      use comsk
      use comsptr

      implicit none

      integer inp1,np1,ipc,flag_divd

      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vvinp,vvold,vvopp

c     in case the velocity has opposite signs on the same face
c     on the two different sides, use the value from the upstream
c     face, ie value out of the control volume of the oldnode2.
c     NOTE: ipc is wrt the current node, and the face has the 
c     opposite sign wrt oldnode2

      if(ipc.eq.+1) then
         call sub_vel_divd(flag_divd,inp1,+1,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),-1,vvold)
         call sub_vel_divd(flag_divd,inp1,-1,vvopp)
c         vvinp=+ggg(inp1,+1)
c         vvold=-ggg(oldnode2(np1),-1)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(-ggg(inp1,-1)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if x velocity on both the left and right face are pointing
c to right, zero out the normal component of the velocity
               aomrx=0.
               vomrx=0.
            else
c but if the velocity on the left face is pointing left
c , ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomrx=vvold
               aomrx=(vvold-vvopp)/ddx(inp1)
            endif
         endif
      elseif(ipc.eq.-1) then
         call sub_vel_divd(flag_divd,inp1,-1,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),+1,vvold)
         call sub_vel_divd(flag_divd,inp1,+1,vvopp)
c         vvinp=-ggg(inp1,-1)
c         vvold=+ggg(oldnode2(np1),+1)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(+ggg(inp1,+1)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if x velocity on both the left and right face are pointing
c left, zero out the normal component of the velocity
               aomrx=0.
               vomrx=0.
            else
c but if the velocity on the right face is pointing to right
c , ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomrx=vvold
               aomrx=(vvopp-vvold)/ddx(inp1)
            endif
         endif
      elseif(ipc.eq.+2) then
         call sub_vel_divd(flag_divd,inp1,+2,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),-2,vvold)
         call sub_vel_divd(flag_divd,inp1,-2,vvopp)
c         vvinp=+ggg(inp1,+2)
c         vvold=-ggg(oldnode2(np1),-2)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(-ggg(inp1,-2)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if y velocity on both the front and back face are pointing
c the wrong way, zero out the normal component of the velocity
               aomry=0.
               vomry=0.
            else
c but if the velocity on the back face is pointing the correct
c way, ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomry=vvold
               aomry=(vvold-vvopp)/ddy(inp1)
            endif
         endif
      elseif(ipc.eq.-2) then
         call sub_vel_divd(flag_divd,inp1,-2,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),+2,vvold)
         call sub_vel_divd(flag_divd,inp1,+2,vvopp)
c         vvinp=-ggg(inp1,-2)
c         vvold=+ggg(oldnode2(np1),+2)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(+ggg(inp1,+2)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if y velocity on both the front and back face are pointing
c the wrong way, zero out the normal component of the velocity
               aomry=0.
               vomry=0.
            else
c but if the velocity on the front face is pointing the correct
c way, ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomry=vvold
               aomry=(vvopp-vvold)/ddy(inp1)
            endif
         endif
      elseif(ipc.eq.+3) then
         call sub_vel_divd(flag_divd,inp1,+3,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),-3,vvold)
         call sub_vel_divd(flag_divd,inp1,-3,vvopp)
c         vvinp=+ggg(inp1,+3)
c         vvold=-ggg(oldnode2(np1),-3)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(-ggg(inp1,-3)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if z velocity on both the upper and lower face are pointing
c the wrong way, zero out the normal component of the velocity
               aomrz=0.
               vomrz=0.
            else
c but if the velocity on the lower face is pointing the correct
c way, ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomrz=vvold
               aomrz=(vvold-vvopp)/ddz(inp1)
            endif
         endif
      elseif(ipc.eq.-3) then
         call sub_vel_divd(flag_divd,inp1,-3,vvinp)
         call sub_vel_divd(flag_divd,oldnode2(np1),+3,vvold)
         call sub_vel_divd(flag_divd,inp1,+3,vvopp)
c         vvinp=-ggg(inp1,-3)
c         vvold=+ggg(oldnode2(np1),+3)
         if((vvold*vvinp).lt.0.) then
c            if(vvinp*(+ggg(inp1,+3)).gt.0.) then
            if(vvinp*vvopp.gt.0.) then
c if z velocity on both the upper and lower face are pointing
c the wrong way, zero out the normal component of the velocity
               aomrz=0.
               vomrz=0.
            else
c but if the velocity on the upper face is pointing the correct
c way, ie there is a flow devide, use the velocity and accelaration
c from the oldnode2
               vomrz=vvold
               aomrz=(vvopp-vvold)/ddz(inp1)
            endif
         endif
      endif

      return

      end subroutine fixggg_2
c......................................................................

      subroutine fixggg_3(flag_divd,inp1,inp1old,np1,ipc,xp1,yp1,zp1,
     $     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
      
c     if the particle is being artifically
c     pushed into a stagnant node due to ficticious velocities 
c     resulting from interpolation at omr nodes. If so, in the cc of inp1
c     zero out velocity in that direction and recalculate the trajectory

      use comsk
      use comsptr

      implicit none

      integer inp1,np1,ipc,inp1old,flag_divd

      real*8 xp1,yp1,zp1
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vvinp,vvold
      real*8 vx1b,vy1b,vz1b,vx2b,vy2b,vz2b

c     in case the velocity has opposite signs on the same face
c     on the two different sides, use the value from the upstream
c     face, ie value out of the control volume of the oldnode.
c     NOTE: ipc is wrt the current node, and the face has the 
c     opposite sign wrt oldnode

      if(ipc.eq.+1) then
         vx2b=0.
         aomrx=(vx2b-vx1b)/ddx(inp1)
         vomrx=vx1b+aomrx*xp1         

      elseif(ipc.eq.-1) then
         vx1b=0.
         aomrx=(vx2b-vx1b)/ddx(inp1)
         vomrx=vx1b+aomrx*xp1         

      elseif(ipc.eq.+2) then
         vy2b=0.
         aomry=(vy2b-vy1b)/ddy(inp1)
         vomry=vy1b+aomry*yp1

      elseif(ipc.eq.-2) then
         vy1b=0.
         aomry=(vy2b-vy1b)/ddy(inp1)
         vomry=vy1b+aomry*yp1

      elseif(ipc.eq.+3) then
         vz2b=0.
         aomrz=(vz2b-vz1b)/ddz(inp1)
         vomrz=vz1b+aomrz*zp1

      elseif(ipc.eq.-3) then
         vz1b=0.
         aomrz=(vz2b-vz1b)/ddz(inp1)
         vomrz=vz1b+aomrz*zp1

      endif

      return

      end subroutine fixggg_3
c.......................................................................

      subroutine pollock_intersect_plane(i,np1,ip,xp1,yp1,zp1,xp2,yp2,      
     &     zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     xc,yc,zc,dtc)
c since this routine is called only when we already know that the particle has 
c exited the CV of node i, the particle is placed slightly outside ipc-plane 

      use comsptr

      implicit none

      integer i,np1,ip
      integer debug_pollock

      real*8 epsilon,xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc
      real*8 alamda,amue,anue,ap,d,cl,ccm,cn,fac1,fac2
      real*8 epomr,dtc,small_slope,epsilon2
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vx1b,vx2b,vy1b,vy2b,vz1b,vz2b

      small_slope=1.e-20
      epsilon=1.e-10
	epsilon2=1.e-12
      epomr=1.e-4
      debug_pollock = 1

c     particle crossed the y-z plane in the x-direction
         
      if(ip.eq.+1) then

c calculating vx2b from vomrx and aomrx because these may have
c been modified in fixggg and no longer equal vx2bv(np1)
         vx2b = vomrx + aomrx*(ddx(i)-xp1)
         if((vomrx.le.0..or.vx2b.le.0.).and.debug_pollock.eq.1) call 
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

c     set the particle at the cell face and calculate the actual time 
c     needed to get there
         xc=ddx(i)*(1.+epsilon2)
         x_new(np1)=vx2b/vomrx
c     s kelkar  april 13 04 3DOMR modified the next conditioner
c     to exclude the case when particle starts on a face. this is a 
c     quick fix. 
         if(abs(aomrx).le.small_slope.or.
     $        abs(x_new(np1)-1.).le.epomr) then
            dtc=(ddxv(np1)-x1(np1))/vomrx
         else
            dtc=dlog(x_new(np1))/aomrx
         end if
c     now that we have delta-time, calculate new yc and zc corrosponding 
c     to the xc on the cell face.
         if(abs(aomry).le.small_slope) then
            yc=y1(np1)+vomry*dtc
         else
            yc=y1(np1)+vomry*(dexp(aomry*dtc)-1.)/aomry
         endif
c make sure yc is between yp1 and yp2
         if((yc-yp2)*(yc-yp1).gt.0.) yc=yp2
         if(abs(aomrz).le.small_slope) then
            zc=z1(np1)+vomrz*dtc
         else
            zc=z1(np1)+vomrz*(dexp(aomrz*dtc)-1.)/aomrz
         endif
         if((zc-zp2)*(zc-zp1).gt.0.) zc=zp2
         
      elseif(ip.eq.-1) then

c calculating vx1b from vomrx and aomrx because these may have
c been modified in fixggg and no longer equal vx2bv(np1)
         vx1b = vomrx - aomrx*xp1
         if((vomrx.ge.0..or.vx1b.ge.0.).and.debug_pollock.eq.1) call  
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

         xc=0.-ddx(i)*epsilon2
         x_new(np1)=vx1b/vomrx
         if(abs(aomrx).le.small_slope.or.
     $        abs(x_new(np1)-1.).le.epomr) then
            dtc=-x1(np1)/vomrx
         else
            dtc=dlog(x_new(np1))/aomrx
         end if
         if(abs(aomry).le.small_slope) then
            yc=y1(np1)+vomry*dtc
         else
            yc=y1(np1)+vomry*(dexp(aomry*dtc)-1.)/aomry
         endif
         if((yc-yp2)*(yc-yp1).gt.0.) yc=yp2
         if(abs(aomrz).le.small_slope) then
            zc=z1(np1)+vomrz*dtc
         else
            zc=z1(np1)+vomrz*(dexp(aomrz*dtc)-1.)/aomrz
         endif
         if((zc-zp2)*(zc-zp1).gt.0.) zc=zp2
         
      elseif(ip.eq.+2) then
c     particle crossed the x-z plane in the y-direction

c calculating vy2b from vomry and aomry because these may have
c been modified in fixggg and no longer equal vy2bv(np1)
         vy2b = vomry + aomry*(ddy(i)-yp1)
         if((vomry.le.0..or.vy2b.le.0.).and.debug_pollock.eq.1) call  
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

         yc=ddy(i)*(1.+epsilon2)
         y_new(np1)=vy2b/vomry
         if(abs(aomry).le.small_slope.or.
     $        abs(y_new(np1)-1.).le.epomr) then
            dtc=(ddyv(np1)-y1(np1))/vomry
         else
            dtc=dlog(y_new(np1))/aomry
         end if
         if(abs(aomrx).le.small_slope) then
            xc=x1(np1)+vomrx*dtc
         else
            xc=x1(np1)+vomrx*(dexp(aomrx*dtc)-1.)/aomrx
         endif
         if((xc-xp2)*(xc-xp1).gt.0.) xc=xp2
         if(abs(aomrz).le.small_slope) then
            zc=z1(np1)+vomrz*dtc
         else
            zc=z1(np1)+vomrz*(dexp(aomrz*dtc)-1.)/aomrz
         endif
         if((zc-zp2)*(zc-zp1).gt.0.) zc=zp2
         
      elseif(ip.eq.-2) then

c calculating vy1b from vomry and aomry because these may have
c been modified in fixggg and no longer equal vy1bv(np1)
         vy1b = vomry - aomry*yp1
         if((vomry.ge.0..or.vy1b.ge.0.).and.debug_pollock.eq.1) call  
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

         yc=0.-ddy(i)*epsilon2
         y_new(np1)=vy1b/vomry
         if(abs(aomry).le.small_slope.or.
     $        abs(y_new(np1)-1.).le.epomr) then
            dtc=-y1(np1)/vomry
         else
            dtc=dlog(y_new(np1))/aomry
         end if
         if(abs(aomrx).le.small_slope) then
            xc=x1(np1)+vomrx*dtc
         else
            xc=x1(np1)+vomrx*(dexp(aomrx*dtc)-1.)/aomrx
         endif
         if((xc-xp2)*(xc-xp1).gt.0.) xc=xp2
         if(abs(aomrz).le.small_slope) then
            zc=z1(np1)+vomrz*dtc
         else
            zc=z1(np1)+vomrz*(dexp(aomrz*dtc)-1.)/aomrz
         endif
         if((zc-zp2)*(zc-zp1).gt.0.) zc=zp2
         
      elseif(ip.eq.+3) then
c     particle crossed the x-y plane in the z-direction

c calculating vz2b from vomrz and aomrz because these may have
c been modified in fixggg and no longer equal vz2bv(np1)
         vz2b = vomrz + aomrz*(ddz(i)-zp1)
         if((vomrz.le.0..or.vz2b.le.0.).and.debug_pollock.eq.1) call  
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

         zc=ddz(i)*(1.+epsilon2)
         z_new(np1)=vz2b/vomrz
         if(abs(aomrz).le.small_slope.or.
     $        abs(z_new(np1)-1.).le.epomr) then
            dtc=(ddzv(np1)-z1(np1))/vomrz
         else
            dtc=dlog(z_new(np1))/aomrz
         end if
         if(abs(aomrx).le.small_slope) then
            xc=x1(np1)+vomrx*dtc
         else
            xc=x1(np1)+vomrx*(dexp(aomrx*dtc)-1.)/aomrx
         endif
         if((xc-xp2)*(xc-xp1).gt.0.) xc=xp2
         if(abs(aomry).le.small_slope) then
            yc=y1(np1)+vomry*dtc
         else
            yc=y1(np1)+vomry*(dexp(aomry*dtc)-1.)/aomry
         endif
         if((yc-yp2)*(yc-yp1).gt.0.) yc=yp2
         
      elseif(ip.eq.-3) then

c calculating vz1b from vomrz and aomrz because these may have
c been modified in fixggg and no longer equal vz1bv(np1)
         vz1b = vomrz - aomrz*zp1
         if((vomrz.ge.0..or.vz1b.ge.0.).and.debug_pollock.eq.1) call  
     2    exit_pollock(i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomrx)

         zc=0.-ddz(i)*epsilon2
         z_new(np1)=vz1b/vomrz
         if(abs(aomrz).le.small_slope.or.
     $        abs(z_new(np1)-1.).le.epomr) then
            dtc=-z1(np1)/vomrz
         else
            dtc=dlog(z_new(np1))/aomrz
         end if
         if(abs(aomrx).le.small_slope) then
            xc=x1(np1)+vomrx*dtc
         else
            xc=x1(np1)+vomrx*(dexp(aomrx*dtc)-1.)/aomrx
         endif
         if((xc-xp2)*(xc-xp1).gt.0.) xc=xp2
         if(abs(aomry).le.small_slope) then
            yc=y1(np1)+vomry*dtc
         else
            yc=y1(np1)+vomry*(dexp(aomry*dtc)-1.)/aomry
         endif
         if((yc-yp2)*(yc-yp1).gt.0.) yc=yp2
         
      End if
      
      return
      
      end subroutine pollock_intersect_plane

c..........................................................................   

      subroutine shift_xc_inside(inp1,newnode,np1,
     2     vomrx,vomry,vomrz,ipc,xc,yc,zc)

      use comsk
      use comsptr

      implicit none

      integer inp1,newnode,ip1,np1,ipc

      real*8 epsilon,xc,yc,zc,ddtc,dxc,dyc,dzc
      real*8 vomrx,vomry,vomrz

      epsilon=1.e-13

      ip1=ipc

      if(ip1.eq.-1) then
c place the point slightly to the left of the plane x=-1
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dxc=-epsilon*ddx(newnode)
         ddtc=abs(dxc/vomrx)
         xc=xc+dxc
         yc=yc+vomry*ddtc
         zc=zc+vomrz*ddtc
c change reference point to corn(newnode)
         xc=xc+corn(inp1,1)-corn(newnode,1)
         yc=yc+corn(inp1,2)-corn(newnode,2)
         zc=zc+corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move it slightly to place it inside
         if(yc.le.0.) then 
            yc=+epsilon*ddy(newnode)
         elseif(yc.ge.ddy(newnode)) then
            yc=ddy(newnode)*(1.-epsilon)
         endif
         if(zc.le.0.) then 
            zc=+epsilon*ddz(newnode)
         elseif(zc.ge.ddz(newnode)) then
            zc=ddz(newnode)*(1.-epsilon)
         endif
      elseif(ip1.eq.+1) then
c place the point slightly to the right of the plane x=+1 (wrt i)
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dxc=epsilon*ddx(newnode)
         ddtc=abs(dxc/vomrx)
         xc=xc+dxc
         yc=yc+vomry*ddtc
         zc=zc+vomrz*ddtc
c change reference point to corn(newnode)
         xc   =xc   +corn(inp1,1)-corn(newnode,1)
         yc   =yc   +corn(inp1,2)-corn(newnode,2)
         zc   =zc   +corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move it slightly to place it inside newnode
         if(yc.le.0.) then 
            yc=+epsilon*ddy(newnode)
         elseif(yc.ge.ddy(newnode)) then
            yc=ddy(newnode)*(1.-epsilon)
         endif
         if(zc.le.0.) then 
            zc=+epsilon*ddz(newnode)
         elseif(zc.ge.ddz(newnode)) then
            zc=ddz(newnode)*(1.-epsilon)
         endif
      elseif(ip1.eq.-2) then
c place the point slightly below  the plane y=-2
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dyc=-epsilon*ddy(newnode)
         ddtc=abs(dyc/vomry)
         yc=yc+dyc
         xc=xc+vomrx*ddtc
         zc=zc+vomrz*ddtc
c change reference point to corn(newnode)
         xc=xc+corn(inp1,1)-corn(newnode,1)
         yc=yc+corn(inp1,2)-corn(newnode,2)
         zc=zc+corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move itr slightly to place it inside
         if(xc.le.0.) then 
            xc=+epsilon*ddx(newnode)
         elseif(xc.ge.ddx(newnode)) then
            xc=ddx(newnode)*(1.-epsilon)
         endif
         if(zc.le.0.) then 
            zc=+epsilon*ddz(newnode)
         elseif(zc.ge.ddz(newnode)) then
            zc=ddz(newnode)*(1.-epsilon)
         endif
      elseif(ip1.eq.+2) then
c place the point slightly above the plane +y
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dyc=epsilon*ddy(newnode)
         ddtc=abs(dyc/vomry)
         yc=yc+dyc
         xc=xc+vomrx*ddtc
         zc=zc+vomrz*ddtc
c change reference point to corn(newnode)
         xc   =xc   +corn(inp1,1)-corn(newnode,1)
         yc   =yc   +corn(inp1,2)-corn(newnode,2)
         zc   =zc   +corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move it slightly to place it inside newnode
         if(xc.le.0.) then 
            xc=+epsilon*ddx(newnode)
         elseif(xc.ge.ddx(newnode)) then
            xc=ddx(newnode)*(1.-epsilon)
         endif
         if(zc.le.0.) then 
            zc=+epsilon*ddz(newnode)
         elseif(zc.ge.ddz(newnode)) then
            zc=ddz(newnode)*(1.-epsilon)
         endif
      elseif(ip1.eq.-3) then
c place the point slightly behind  the plane -z
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dzc=-epsilon*ddz(newnode)
         ddtc=abs(dzc/vomrz)
         zc=zc+dzc
         xc=xc+vomrx*ddtc
         yc=yc+vomry*ddtc
c change reference point to corn(newnode)
         xc=xc+corn(inp1,1)-corn(newnode,1)
         yc=yc+corn(inp1,2)-corn(newnode,2)
         zc=zc+corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move itr slightly to place it inside
         if(xc.le.0.) then 
            xc=+epsilon*ddx(newnode)
         elseif(xc.ge.ddx(newnode)) then
            xc=ddx(newnode)*(1.-epsilon)
         endif
         if(yc.le.0.) then 
            yc=+epsilon*ddy(newnode)
         elseif(yc.ge.ddy(newnode)) then
            yc=ddy(newnode)*(1.-epsilon)
         endif
      elseif(ip1.eq.+3) then
c place the point slightly in front of the plane +z
c note that xc,yc,zc,ipc_save are still wrt node i and not newnode
         dzc=epsilon*ddz(newnode)
         ddtc=abs(dzc/vomrz)
         zc=zc+dzc
         xc=xc+vomrx*ddtc
         yc=yc+vomry*ddtc
c change reference point to corn(newnode)
         xc   =xc   +corn(inp1,1)-corn(newnode,1)
         yc   =yc   +corn(inp1,2)-corn(newnode,2)
         zc   =zc   +corn(inp1,3)-corn(newnode,3)
c check if the point is actually within the cc of newnode
c if not, move it slightly to place it inside newnode
         if(xc.le.0.) then 
            xc=+epsilon*ddx(newnode)
         elseif(xc.ge.ddx(newnode)) then
            xc=ddx(newnode)*(1.-epsilon)
         endif
         if(yc.le.0.) then 
            yc=+epsilon*ddy(newnode)
         elseif(yc.ge.ddy(newnode)) then
            yc=ddy(newnode)*(1.-epsilon)
         endif

      endif
 
      return

      end subroutine shift_xc_inside

c.................................................................... 

      subroutine  update_newnode(inp1,np1, ipc,newnode,
     $  dtc,xpnew, ypnew, zpnew)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer i,inp1,icapture,np1,icnl,newnode,ipc
      integer iplane(3),nplanes,ipsign,ipab

      real*8 xpnew, ypnew, zpnew, dtc
      
c     update arrays
c     set the x2 arrays- they have the local coords of the point 
c     with reference to the corn of the new node. xpnew is already 
c     corrected for this

      dt(np1)=dtc
         
      x2(np1)=xpnew
      y2(np1)=ypnew
      z2(np1)=zpnew
c up to this point, ipc is wrt inp1. But all the updated arrays are 
c wrt noewnode. Hence switch the sign of ipc so that it becomes
c referenced to the new control voulme.
      ipc_save(np1)=-ipc
      oldnode2(np1)=oldnode(np1)
      oldnode(np1)=inp1
      ijkvs(np1)=inp1
      ijkv(np1)=newnode
      
      x3(np1) = x2(np1) + corn(ijkv(np1), 1)
      y3(np1) = y2(np1) + corn(ijkv(np1), 2)
      z3(np1) = z2(np1) + corn(ijkv(np1), 3)

c ddxv's are needed in new_part_loc_struct & move_particle
      ddxv(np1)=ddx(ijkv(np1))
      ddyv(np1)=ddy(ijkv(np1))
      ddzv(np1)=ddz(ijkv(np1))
      
      return
      
      end subroutine  update_newnode
      
c........................................................................

      subroutine update_samenode(inp1,np1,ipc,newnode,
     $     dtc,xp2,yp2,zp2)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer i,inp1,np1,icnl,newnode,ipc

      real*8 xp2, yp2, zp2, dtc
      
      dt(np1)=dtc

      newnode=ijkv(np1)
      ijkvs(np1) = inp1
      x2(np1)=xp2
      y2(np1)=yp2
      z2(np1)=zp2
      x3(np1) = x2(np1) + corn(ijkv(np1), 1)
      y3(np1) = y2(np1) + corn(ijkv(np1), 2)
      z3(np1) = z2(np1) + corn(ijkv(np1), 3)
      
      return

      end subroutine update_samenode

c........................................................................

      subroutine update_exit(inpin,np1,ipc,newnode,
     $     dtc,xc,yc,zc)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer i,inp1,np1,icnl,newnode,ipc,inpin

      real*8 xc, yc, zc, dtc
      
      inp1=abs(inpin)
      dt(np1)=dtc

      x2(np1)=xc
      y2(np1)=yc
      z2(np1)=zc
      x3(np1) = x2(np1) + corn(inp1, 1)
      y3(np1) = y2(np1) + corn(inp1, 2)
      z3(np1) = z2(np1) + corn(inp1, 3)
      
      exit_node(np1)=inp1
      ijkvs(np1)=inpin
      ijkv(np1)=newnode
      istop(np1)=1
      
      return

      end subroutine update_exit

c....................................................................

      subroutine remove_particle(np1,k,inp1,ipc)

      use comai, only : ierr, iptty 
      use comsptr
      use comsk

      implicit none

      integer np1,k,inp1,ipc

c      call pause_compute_exit

      istop(np1)=3

      if(k.eq.1) then
!         write(ierr,*)'WARNING: REMOVING particle #,inp1:',np1,inp1
!         write(ierr,*)'count_steps.gt.count_steps_max = ', 
!     &        count_steps_max
!         if (iptty .ne. 0) then
!            write(iptty,*)'WARNING: REMOVING particle #,inp1:',np1,inp1
!            write(iptty,*)'count_steps.gt.count_steps_max = ', 
!     &           count_steps_max
!         end if
!      elseif(k.eq.2) then
!         write(ierr,*)'WARNING: REMOVING particle #',np1
!         write(ierr,*)'abrat.le.eprat'
!         if (iptty .ne. 0) then
!            write(iptty,*)'WARNING: REMOVING particle #',np1
!            write(iptty,*)'abrat.le.eprat'
!         end if
      elseif(k.eq.5) then
         istop(np1)=2
!      elseif(k.eq.8) then
!         write(ierr,*)'WARNING: REMOVING particle #',np1
!         write(ierr,*)'Error in newdt2(compute_exit_new), dtnew=1.e30'
!         write(ierr,*)'np1,inp1,ipc'
!         write(ierr,*) np1,inp1,ipc
!         if (iptty .ne. 0) then
!            write(iptty,*)'WARNING: REMOVING particle #',np1
!            write(iptty,*)'Error in newdt2(compute_exit_new),',
!     &           ' dtnew=1.e30'
!            write(iptty,*)'np1,inp1,ipc'
!            write(iptty,*) np1,inp1,ipc
!         end if
      endif

      return

      end subroutine remove_particle

c.......................................................

      subroutine slow_particle(flag_divd,inp1,ipc,np1,newnode,vexit,
     1     eprat,abrat,flag_slow)

      use comsptr
      use comsk

      implicit none

      integer np1,inp1,ipc, newnode,j,flag_slow,ipsign,ipab
      integer flag_divd

      real*8 eprat,abrat,abv_current,abv_new,agggmax,aggg
      real*8 agggrat,v_current,v_new,v_opposite,vexit

      flag_slow=0

      ipsign=isign(1,ipc)
      ipab=iabs(ipc)
c note: sign if ipc in newnode is reversed wrt inp1
      call sub_vel_divd(flag_divd,newnode,-ipc,v_new) 
      call sub_vel_divd(flag_divd,newnode,+ipc,v_opposite) 
c      v_new    = -ipsign*ggg(newnode,-ipc)
c      v_opposite=+ipsign*ggg(newnode,+ipc)

      abv_current=abs(vexit)
      abv_new=abs(v_new)
      if(abv_current.gt.0.) then
         abrat=abv_new/abv_current
      else
         abrat=1.
      endif

c if particle has jumped into a slow moving place, check if the 
c new cc has all velocities stagnant by looking at max-ggg
      if(abrat.lt.eprat) then
         agggmax=-1.e20
         do j=-3,3
            if (j.ne.0) then
               callsub_vel_divd(flag_divd,newnode,j,aggg)
               aggg=abs(aggg)
c               aggg=abs(ggg(newnode,j))
               if(agggmax.lt.aggg)  agggmax=aggg
            endif
         enddo
         if(abv_current.gt.0.) then
            agggrat=agggmax/abv_current
         else
            agggrat=1.
         endif
c if max-norm of ggg in newnode is small, set flag to reset velocities
c in inp1. If not, do another check on the signs of velotity
c in newnode corrosponding to ipc: if they are both opposit
c of the velocity from inp1, set flag top reset velocities in newnode
         if(agggrat.lt.eprat) then
            flag_slow=1
         else
            if (vexit*v_new.lt.0.) then
               if (vexit*v_opposite.lt.0.) then
                  flag_slow=1
               endif
            endif
         endif
      endif
      
      return
      
      end subroutine slow_particle

c........................................................

      subroutine  exit_pollock
     2 (i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomr)

      use comai, only : ierr, iptty 

      integer i,np1,ip
      real*8 xp1,yp1,zp1,vomrx,vomry,vomrz,vomr

      write(ierr,*)'pollock_intersect_plane. Sign of velocity wrong.'
      write(ierr,*)'i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomr-ip'
      write(ierr,*)i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomr

      if (iptty .ne. 0) then
      write(iptty,*)'pollock_intersect_plane. Sign of velocity wrong.'
      write(iptty,*)'i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomr-ip'
      write(iptty,*)i,np1,ip,xp1,yp1,zp1,vomrx,vomry,vomrz,vomr
      end if

c      stop
      call pause_compute_exit
      call remove_particle(np1,1,i,1)

      return

      end subroutine exit_pollock

c............................................................

      subroutine newdt3(inp1,np1,xp1,yp1,zp1,xvel,yvel,zvel,
     1     aomrx,aomry,aomrz,
     2     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,dtnew)
c do a simple estimate of dt for estimating the direction of
c the next step. if the particle reaches the boundary,
c then the actual time step is calculated in 
c pollock-intersect_plane.

      use comsk
      use comsptr

      implicit none

      integer inp1,np1,newnode,ipab,ipc,ipcx,ipcy,ipcz

      real*8 xvel,yvel,zvel,del,epsilon,dtnew,xdt,ydt,zdt
      real*8 xp1,yp1,zp1,dt_factor,aomrx,aomry,aomrz
      real*8 vx2b, vx1b,vy2b,vy1b,vz2b,vz1b,epfac,fac1

      dt_factor=1./3.

      epsilon=1.e-22
      epfac=1.e-8
      xdt=1.e+30
      ydt=1.e+30
      zdt=1.e+30
      
      if(xvel.gt.+epsilon.and.vx2b.gt.+epsilon.and.
     $     ipc_save(np1).ne.+1) then
c the first part of this conditional is checking for a situation where
c there is a flow devide in the middle of the control volume. The particle  
c shopuld not go there, and Pollock gives wrong answers. 3-rd conditional:
c this is the case where the particle might exit on the right hand
c x boundary. In that case we want to make sure that it did not enter the 
c current cc from that face in order to avoid bouncing in case of a flow
c devide situation. at this stage this face is given by ipc_save(np1),
c Ipc_save is stored in the last step in "update_position" wrt the 
c "newnode" for the last step, which is the current node for this step.
         del=(ddx(inp1)-xp1)*(1.+epfac)
         fac1=1.+aomrx*del/xvel
         if (abs(aomrx).gt.epsilon.and.fac1.gt.epsilon) then
            xdt=(dlog(fac1))/aomrx
         else
            xdt=(del/xvel)
         endif
         ipcx=+1
      elseif(xvel.lt.-epsilon.and.vx1b.lt.-epsilon.and.
     $        ipc_save(np1).ne.-1)then
         del=-xp1*(1.+epfac)
         fac1=1.+aomrx*del/xvel
         if (abs(aomrx).gt.epsilon.and.fac1.gt.epsilon) then
            xdt=(dlog(fac1))/aomrx
         else
            xdt=(del/xvel)
         endif
         ipcx=-1
      endif
      if(yvel.gt.+epsilon.and.vy2b.gt.+epsilon.and.
     $     ipc_save(np1).ne.+2) then
         del=(ddy(inp1)-yp1)*(1.+epfac)
         fac1=1.+aomry*del/yvel
         if (abs(aomry).gt.epsilon.and.fac1.gt.epsilon) then
            ydt=(dlog(fac1))/aomry
         else
            ydt=+(del/yvel)
         endif
         ipcy=+2
      elseif(yvel.lt.-epsilon.and.vy1b.lt.-epsilon.and.
     $        ipc_save(np1).ne.-2)then
         del=-yp1*(1.+epfac)
         fac1=1.+aomry*del/yvel
         if(abs(aomry).gt.epsilon.and.fac1.gt.epsilon) then
            ydt=(dlog(fac1))/aomry
         else
            ydt=(del/yvel)
         endif
         ipcy=-2
      endif
      if(zvel.gt.+epsilon.and.vz2b.gt.+epsilon.and.
     $     ipc_save(np1).ne.+3) then
         del=(ddz(inp1)-zp1)*(1.+epfac)
         fac1=1.+aomrz*del/zvel
         if(abs(aomrz).gt.epsilon.and.fac1.gt.epsilon) then
            zdt=(dlog(fac1))/aomrz
         else
            zdt=+(del/zvel)
         endif
         ipcz=+3
      elseif(zvel.lt.-epsilon.and.vz1b.lt.-epsilon.and.
     $        ipc_save(np1).ne.-3)then
         del=-zp1*(1.+epfac)
         fac1=1.+aomrz*del/zvel
         if(abs(aomrz).gt.epsilon.and.fac1.gt.epsilon) then
            zdt=(dlog(fac1))/aomrz
         else
            zdt=-(del/zvel)
         endif
         ipcz=-3
      endif
      
c find the min of the nonzero times

      dtnew=1.e+30

      if((xdt.lt.ydt.and.xdt.gt.0.).or.ydt.le.0.) then
         if(xdt.lt.zdt.or.zdt.le.0.) then
            dtnew=xdt
c            ipc=ipcx
         else
            if(zdt.gt.0.) then
               dtnew=zdt
c               ipc=ipcz
            endif
         endif
      else
         if((ydt.lt.zdt.and.ydt.gt.0.).or.zdt.le.0.) then
            dtnew=ydt
c            ipc=ipcy
         else
            if(zdt.gt.0.) then
               dtnew=zdt
c               ipc=ipcz
            endif
         endif
      endif

      dtnew=courant_factor*dtnew
      dtnew = min(x61(np1),dtnew)

c     Reset dtnew to hit the end of the time step exactly if it 
c     will be in the current cell at that time. Set flag ioutt
c     so record this
      
      if(ioutt(np1).eq.0) then
         if(dtnew/86400..gt.tt1-ttp1(np1)) then
            dtnew=86400.*(tt1-ttp1(np1))
            ioutt(np1) = 1
         end if
      end if
      
c     Add check to make sure dt is positive
      
      if(dtnew.lt.0.) call pause_compute_exit
      dtnew = max(1.d-10,dtnew)
      
      if(dtnew.ge.1.e+30) then
         call remove_particle(np1,8,inp1,ipc_save(np1))
      endif
      
      return
      
      end subroutine newdt3

c......................................................................

      subroutine pause_compute_exit

      return

      end subroutine pause_compute_exit

c.....................................................................

      subroutine  time_zheng_new(icnl,np1,
     1     xp1,yp1,zp1,rw,vx_zheng,vy_zheng,vz_zheng,dtnew) 

      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,icnl

      real*8 tfac,epsilon
      real*8 dx,dy,dz,x,y,z,xx,yy,zz,skz,rr
      real*8 qwp,qqw,a,fac,xyw,pi,xp1,yp1,zp1

      real*8 vx_zheng,vy_zheng,vz_zheng,rw,dtnew

      pi=3.14159
      epsilon=1.e-20
      tfac = 0.1

c use equations from Zheng, 1994

      inp1=ijkv(np1)
      dx=ddx(inp1)
      dy=ddy(inp1)
      dz=ddz(inp1)
      x=xp1-0.5*dx
      y=yp1-0.5*dy
      z=zp1-0.5*dz

c distance to the particle from the welbore, for timestep size
      rr=sqrt(x*x+y*y)

      call vel_zheng(inp1,icnl,x,y,z,dx,dy,dz,
     1     vx_zheng,vy_zheng,vz_zheng)

      if(abs(vx_zheng).gt.epsilon) then
c using the distance to the wellbore, x, as the  sale
         dtx(np1)=abs(tfac*rr/vx_zheng)
      else
         dtx(np1)=1.e+30
      endif

      if(abs(vy_zheng).gt.epsilon) then
c using the distance to the wellbore, x, as the  sale
         dty(np1)=abs(tfac*rr/vy_zheng)
      else
         dty(np1)=1.e+30
      endif
      
      if(icnl.eq.0) then
         if(abs(vz_zheng).gt.epsilon) then
            dtz(np1)=abs(tfac*dz/vz_zheng)
         else
            dtz(np1)=1.e+30
         endif
      endif

      dtnew = min(x61(np1),dtx(np1),dty(np1),dtz(np1))

c     Reset dtnew to hit the end of the time step exactly if it 
c     will be in the current cell at that time. Set flag ioutt
c     so record this
      
      if(ioutt(np1).eq.0) then
         if(dtnew/86400..gt.tt1-ttp1(np1)) then
            dtnew=86400.*(tt1-ttp1(np1))
            ioutt(np1) = 1
         end if
      end if
      
c     Add check to make sure dt is positive
      
      if(dtnew.lt.0.) call pause_compute_exit
      dtnew = max(1.d-10,dtnew)
      
      if(dtnew.ge.1.e+30) then
         call remove_particle(np1,8,inp1,ipc_save(np1))
      endif
      
      return
      
      end  subroutine  time_zheng_new

c.........................................................

      subroutine exit_zheng_new(icapture, icnl,np1,
     1     xp1,yp1,zp1,vx_zheng,vy_zheng,vz_zheng,rw,dtnew,
     2     xp2,yp2,zp2) 

c.....1/22/01 s kelkar particle capture
c     use equations from Zheng, 1994. velocities are already calculated
c     in time_zheng

      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,icapture,np1,icnl,newnode,ipc
      integer iplane(3),nplanes

      real*8 ep5,ep,small_slope,x2dim,y2dim,z2dim
      real*8 dx,dy,dz,x,y,z,xx,yy,zz,pi
      real*8 vx_zheng,vy_zheng,vz_zheng,rw,dtnew

      real*8 epsilon,xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc

      epsilon=1.e-10
      ep5 = epsilon
      pi=3.14159

      inp1=ijkv(np1)
      dx=ddx(inp1)
      dy=ddy(inp1)
      dz=ddz(inp1)
      x=xp1-0.5*dx
      y=yp1-0.5*dy
      z=zp1-0.5*dz
      
      if(irevers.eq.1) then
         xx=x-vx_zheng*dtnew
         yy=y-vy_zheng*dtnew
      else
         xx=x+vx_zheng*dtnew
         yy=y+vy_zheng*dtnew
      endif
      
c     check if particle has entered the well.
      if(irevers.ne.1) call check_well_capture(inp1,x,y,xx,yy,icapture,
     1     vx_zheng,vy_zheng,vz_zheng,rw)
      
      if(icapture.ne.1) then
c     particle has not entered the well, compute new location
         
         if(icnl.eq.0) then
            
            if(irevers.eq.1) then
               z=z-vz_zheng*dtnew
            else
               z=z+vz_zheng*dtnew
            endif
            
         endif
         
         xp2=xx+0.5*dx
         yp2=yy+0.5*dy
         xp2=(1.+ep5)*xp2-ep5*xp1
         yp2=(1.+ep5)*yp2-ep5*yp1
         if(icnl.eq.0) then
            zp2=z+0.5*dz
            zp2=(1.+ep5)*zp2-ep5*zp1
         endif
         
      endif
      
      return

      end subroutine exit_zheng_new

c.........................................................

      subroutine exit_point_zheng(inp1,np1,xp1,yp1,zp1,
     1     xp2,yp2,zp2,nplanes,iplane,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     3     ipc,xc,yc,zc,dtc,newnode)

c     find exit point(xc,yc,zc with reference to corn(inp1,)) 
c     and exit plane(ipc), and the new neighbour
c     xc is taken to be the intersection of the STRAIGHT line xp1-xp2
c     with the plane ipc- this is an approximation of the 
c     Pollock trajectory
      
      use comsptr
      use comsk
      
      implicit none
      
      integer inp1,np1,icnl,newnode
      integer iplane(3),nplanes,ipc
      
      real*8 xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc,dtc
      real*8 vx_zheng,vy_zheng,vz_zheng,rw
      
      newnode=0
      ipc=0
      
c     if nplanes > 0 then particle has entered a new cc
      if(nplanes.eq.1) then
         call find_newpoint1_zheng(inp1,np1,iplane(1),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     2        vx_zheng,vy_zheng,vz_zheng,rw,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.2) then
         call find_newpoint2_zheng(inp1,np1,iplane(1),iplane(2),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     2        vx_zheng,vy_zheng,vz_zheng,rw,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.3) then
         call find_newpoint3_zheng(inp1,np1,iplane(1),iplane(2),
     1        iplane(3), xp1,yp1,zp1,xp2,yp2,zp2,
     2        vx_zheng,vy_zheng,vz_zheng,rw,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      endif
      
c     if the point lands on a face,
c     fix xc,yc,zc, to be slightly inside the new control volume 
      if(newnode.gt.0) call shift_xc_inside(inp1,newnode,np1,
     2     vx_zheng,vy_zheng,vz_zheng,
     3     ipc,xc,yc,zc)
      
      return
      
      end subroutine exit_point_zheng

c...........................................................
      
      subroutine find_newpoint1_zheng(i,np1,ip1,
     1     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     3     ipc,xc,yc,zc,dtc)
c     s kelkar march 29 04
c     find exit point accross the plane given by ip
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip,ip1,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 dtc,dtc1,xc1,yc1,zc1
      real*8 vx_zheng,vy_zheng,vz_zheng,rw
      
      call zheng_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc1,yc1,zc1,dtc1)
      
      xc=xc1
      yc=yc1
      zc=zc1
      dtc=dtc1
      ipc=ip1
      
      return
      
      end subroutine find_newpoint1_zheng
      
c...........................................................
      
      subroutine find_newpoint2_zheng(i,np1,ip1,ip2,
     &     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     3     ipc,xc,yc,zc,dtc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1, and ip2 and taking the closest of the two
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,d1,d2,dc,dtc,dtc1,dtc2
      real*8 vx_zheng,vy_zheng,vz_zheng,rw
      
      call zheng_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc1,yc1,zc1,dtc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call zheng_intersect_plane(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc2,yc2,zc2,dtc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
         dtc=dtc1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
         dtc=dtc2
      endif
      
      return
      
      end subroutine find_newpoint2_zheng
      
c............................................................
      
      subroutine find_newpoint3_zheng(i,np1,ip1,ip2,ip3,
     &     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     3     ipc,xc,yc,zc,dtc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1,ip2 and ip3 and taking the closest of the three
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,ip3,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,d1,d2,d3,dc
      real*8 dtc,dtc1,dtc2,dtc3
      real*8 vx_zheng,vy_zheng,vz_zheng,rw
      
      call zheng_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc1,yc1,zc1,dtc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call zheng_intersect_plane(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc2,yc2,zc2,dtc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      call zheng_intersect_plane(i,np1,ip3,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     &     xc3,yc3,zc3,dtc3)
      d3=sqrt((xp1-xc3)**2.+(yp1-yc3)**2.+(zp1-zc3)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
         dtc=dtc1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
         dtc=dtc2
      endif
      if(dc.gt.d3) then
         xc=xc3
         yc=yc3
         zc=zc3
         ipc=ip3
         dc=d3
         dtc=dtc3
      endif
      
      return
      
      end subroutine find_newpoint3_zheng
c......................................................................

      subroutine zheng_intersect_plane(i,np1,ip,xp1,yp1,zp1,xp2,yp2,      
     &     zp2,
     2     vx_zheng,vy_zheng,vz_zheng,rw,
     3     xc,yc,zc,dtc)

c this is a xp1-xp2 straight line intersecting the plane at xc etc

      use comsptr

      implicit none

      integer i,np1,ip

      real*8 epsilon,xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc
      real*8 alamda,amue,anue,ap,d,cl,ccm,cn,fac1,fac2
      real*8 epomr,dtc,small_slope
      real*8 vx_zheng,vy_zheng,vz_zheng,rw

      small_slope=1.e-20
      epsilon=1.e-10
      epomr=1.e-4

c     particle crossed the y-z plane in the x-direction
         
      if(ip.eq.+1) then

c     set the particle at the cell face and calculate the actual time 
c     needed to get there
         xc=ddx(i)
         dtc=(xc-x1(np1))/vx_zheng
c     now that we have delta-time, calculate new yc and zc corrosponding 
c     to the xc on the cell face.
         yc=y1(np1)+vy_zheng*dtc
         zc=z1(np1)+vz_zheng*dtc
         
      elseif(ip.eq.-1) then

         xc=0.
         dtc=-x1(np1)/vx_zheng
         yc=y1(np1)+vy_zheng*dtc
         zc=z1(np1)+vz_zheng*dtc
         
      elseif(ip.eq.+2) then
c     particle crossed the x-z plane in the y-direction

         yc=ddy(i)
         dtc=(yc-y1(np1))/vy_zheng
         xc=x1(np1)+vx_zheng*dtc
         zc=z1(np1)+vz_zheng*dtc
         
      elseif(ip.eq.-2) then

         yc=0.
         dtc=-y1(np1)/vy_zheng
         xc=x1(np1)+vx_zheng*dtc
         zc=z1(np1)+vz_zheng*dtc
         
      elseif(ip.eq.+3) then
c     particle crossed the x-y plane in the z-direction

         zc=ddz(i)
         dtc=(zc-z1(np1))/vz_zheng
         xc=x1(np1)+vx_zheng*dtc
         yc=y1(np1)+vy_zheng*dtc
         
      elseif(ip.eq.-3) then

         zc=0.
         dtc=-z1(np1)/vz_zheng
         xc=x1(np1)+vx_zheng*dtc
         yc=y1(np1)+vy_zheng*dtc
         
      End if
      
      return
      
      end subroutine zheng_intersect_plane

c......................................................................

      subroutine initialize_particle(flag_divd,inp1,oldinp,np1,
     $     ipcsave,xp1,yp1,
     1     zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)

c s kelkar  Dec 14, 04986 9108
c initialize the particle position and velocities and derivatives
c note that this can be a complicated matter in some omr cc's

      use comai, only : day
      use comdi
      use comsptr
      use comsk

      implicit none

      integer i,j,k,inp1,np1,ipc, ipcsave,oldinp
      integer current_model,current_flag
      integer flag_divd, ixy, h_flag

      real*8 xp1,yp1,zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 vx1b,vy1b,vz1b,vx2b,vy2b,vz2b
      real*8 vdisp(3),vdphi(3),voldinp,dzw,vwt,snoise
      real*8 epsilon_ex

      epsilon_ex = 1.e-12
      snoise=0.0002

      xp1=x1(np1)
      yp1=y1(np1)
      zp1=z1(np1)

c s kelkar 2006 Jan 31 fix position for corner nodes on the top 
c surface of the model so that the particle doesnt get caught
c in the low velocity zone. Move the particle epsilon_corner 
c away from the boundary face if it is closer than that.
      if(corner_flag) then
c         call fix_corner_position(inp1,xp1,yp1,zp1)
c fix_corner_position uses iboulist which is not defined for a non-omr problem
      endif
c

c      vx1b=-ggg(inp1,-1)
c      vy1b=-ggg(inp1,-2)
c      vz1b=-ggg(inp1,-3)
c      vx2b=+ggg(inp1,+1)
c      vy2b=+ggg(inp1,+2)
c      vz2b=+ggg(inp1,+3)
c actual velocities are modified by adding Div(D) and 
c grad(porosity) terms in sub_vel_divd()
      call sub_vel_divd(flag_divd,inp1,-1,vx1b)
      call sub_vel_divd(flag_divd,inp1,+1,vx2b)
      call sub_vel_divd(flag_divd,inp1,-2,vy1b)
      call sub_vel_divd(flag_divd,inp1,+2,vy2b)
      if(ifree.ne.0) then
         if(izone_free_nodes(inp1).le.1) then
c     fully saturated or not a water table node,
            call sub_vel_divd(flag_divd,inp1,-3,vz1b)
            call sub_vel_divd(flag_divd,inp1,+3,vz2b)
         else
c     partially saturated water table node, set aomrz and vz2b 
c     so that Vz=ddz*(dSw/dt) at the water table. do not include
c     Del(D)-z terms
c*****Aug 22, 05 next line a quick fix hard wired to smooth out 
c*****the effect of water table flucutaions
            if(abs(s(inp1)-so(inp1)).gt.snoise) then
               vwt= ddz(inp1)*((s(inp1)-so(inp1))/(day*86400.))
            else
               vwt= 0.
            endif
            vz1b=-ggg(inp1,-3)
            dzw=ddz(inp1)*s(inp1)
            vz2b=vz1b+((vwt-vz1b)/dzw)*ddz(inp1)

c****
c 12/21/05 particles very near the water table can get stuck 
c because they have ~0 vertical velocity. 
c check if the particle has an exit in the horizontal
c if not, set the z velocity of the particle to that on the 
c exit face.
            h_flag=0
            if(vx1b.lt.-epsilon_ex) h_flag=1
            if(vx2b.gt.+epsilon_ex) h_flag=1
            if(vy1b.lt.-epsilon_ex) h_flag=1
            if(vy2b.gt.+epsilon_ex) h_flag=1
            if(h_flag.eq.0) then
               if(vz1b.lt.-epsilon_ex) then
                  vz2b=vz1b
               else
                  if(vz2b.gt.+epsilon_ex) vz1b=vz2b
               endif
            endif
c***************************

c****
c     in case there are neighbours with lower water saturation in 
c     the x-y plane, set velocity=opposite face so that the particle
c     does not slow down
            ixy=irray(inp1,+1)
            if(ixy.gt.0) then
               if((s(ixy)-s(inp1)).le.-snoise) vx2b=vx1b
            endif
            ixy=irray(inp1,-1)
            if(ixy.gt.0) then
               if((s(ixy)-s(inp1)).le.-snoise) vx1b=vx2b
            endif
            ixy=irray(inp1,+2)
            if(ixy.gt.0) then
               if((s(ixy)-s(inp1)).le.-snoise) vy2b=vy1b
            endif
            ixy=irray(inp1,-2)
            if(ixy.gt.0) then
               if((s(ixy)-s(inp1)).le.-snoise) vy1b=vy2b
            endif
c****
         endif  
      else        
c     not a water table node
            call sub_vel_divd(flag_divd,inp1,-3,vz1b)
            call sub_vel_divd(flag_divd,inp1,+3,vz2b)
      endif

c compute velocity gradients
      aomrx=(vx2b-vx1b)/ddx(inp1)
      aomry=(vy2b-vy1b)/ddy(inp1)
      aomrz=(vz2b-vz1b)/ddz(inp1)

c if inp1 happens to be a cc with all non-ipc_save velocites
c flowing 'inwards' (at this point, ipcsave is wrt inp1), then
c check if the face accross which the particle entered the cc
c has velocities in opposite directions accross
c from it, and if so, use the velocities from the upstream node      
      if(oldinp.ne.0.and.oldinp.ne.inp1) then
         if(ipcsave.eq.+1) then
            if(vy1b.gt.0..and.vy2b.lt.0.) then
               if(vz1b.gt.0..and.vz2b.lt.0.) then
                  call sub_vel_divd(flag_divd,oldinp,-1,voldinp)
                  if(vx2b*voldinp.lt.0.) then
                     vx2b=voldinp
                     aomrx=(vx2b-vx1b)/ddx(inp1)
                  endif
               endif
            endif
         elseif(ipcsave.eq.-1) then
            if(vy1b.gt.0..and.vy2b.lt.0.) then
               if(vz1b.gt.0..and.vz2b.lt.0.) then
c                  if(vx1b*ggg(oldinp,+1).lt.0.) then
c                     vx1b=+ggg(oldinp,+1)
                  call sub_vel_divd(flag_divd,oldinp,+1,voldinp)
                  if(vx1b*voldinp.lt.0.) then
                     vx1b=voldinp
                     aomrx=(vx2b-vx1b)/ddx(inp1)
                  endif
               endif
            endif
         elseif(ipcsave.eq.+2) then
            if(vx1b.gt.0..and.vx2b.lt.0.) then
               if(vz1b.gt.0..and.vz2b.lt.0.) then
c                  if(vy2b*(-ggg(oldinp,-2)).lt.0.) then
c                     vy2b=-ggg(oldinp,-2)
                  call sub_vel_divd(flag_divd,oldinp,-2,voldinp)
                  if(vy2b*voldinp.lt.0.) then
                     vy2b=voldinp
                     aomry=(vy2b-vy1b)/ddy(inp1)
                  endif
               endif
            endif
         elseif(ipcsave.eq.-2) then
            if(vx1b.gt.0..and.vx2b.lt.0.) then
               if(vz1b.gt.0..and.vz2b.lt.0.) then
c                  if(vy1b*(+ggg(oldinp,+2)).lt.0.) then
c                     vy1b=+ggg(oldinp,+2)
                  call sub_vel_divd(flag_divd,oldinp,+2,voldinp)
                  if(vy1b*voldinp.lt.0.) then
                     vy1b=voldinp
                     aomry=(vy2b-vy1b)/ddy(inp1)
                  endif
               endif
            endif
         elseif(ipcsave.eq.+3) then
            if(vx1b.gt.0..and.vx2b.lt.0.) then
               if(vy1b.gt.0..and.vy2b.lt.0.) then
c                  if(vz2b*(-ggg(oldinp,-3)).lt.0.) then
c                     vz2b=-ggg(oldinp,-3)
                  if(ifree.ne.0) then
                     if(izone_free_nodes(oldinp).le.1) then
c     fully saturated or not a water table node
                        call sub_vel_divd(flag_divd,oldinp,-3,voldinp)
                     else
c     oldinp is partially saturated water table node,  do not include
c     Del(D)-z terms
                        voldinp=-ggg(inp1,-3)
                     endif
                  else
c     not a water table node
                        call sub_vel_divd(flag_divd,oldinp,-3,voldinp)
                  endif
                  if(vz2b*voldinp.lt.0.) then
                     vz2b=voldinp
                     aomrz=(vz2b-vz1b)/ddz(inp1)
                  endif
               endif
            endif
         elseif(ipcsave.eq.-3) then
            if(vx1b.gt.0..and.vx2b.lt.0.) then
               if(vy1b.gt.0..and.vy2b.lt.0.) then
c                  if(vz1b*(+ggg(oldinp,+3)).lt.0.) then
c                     vz1b=+ggg(oldinp,+3)
                  if(ifree.ne.0) then
                     if(izone_free_nodes(inp1).le.1) then
c     fully saturated or not a water table node
                        call sub_vel_divd(flag_divd,oldinp,+3,voldinp)
                     else
c     inp1 is a partially saturated water table node. Do not include 
c     Del(D)-z term in the
c     calculation of voldinp even though oldinp is not a water table node
                        voldinp=+ggg(oldinp,+3)
                     endif
                  else
c     not a water table problem
                        call sub_vel_divd(flag_divd,oldinp,+3,voldinp)
                  endif
                  if(vz1b*voldinp.lt.0.) then
                     vz1b=voldinp
                     aomrz=(vz2b-vz1b)/ddz(inp1)
                  endif
               endif
            endif
         endif
      endif

      call get_velocities(inp1,np1,xp1,yp1,zp1,
     1     aomrx,aomry,aomrz,vx1b,vy1b,vz1b,vomrx,vomry,vomrz)

c-deal with cliff nodes s kelkar jan 9 0-----------------------
      if(cliff_flag) then

c     if inp1 is a surface node on a stair-step on an internal corner
c     then it does not have a brickshape cc. create reflecting
c     boundaries for this. Such nodes are tagged in Tag_stairstep
c     as a part of the mesh generating. 
c     n_stairstep= # of such nodes. For these nodes,
c     irray(i,0)=-100000000-icounter for non-OMR nodes and
c     =-200000000-icounter for OMR nodes
c     where icounter points to the
c     slot in the array istep_cases. istep_cases(icounter) is in turn
c     another pointer that points to the slot in the array istep_cases
c     where information for the node inp1 is saved.
c     ipointer1=irray(i,0)-100000000
c     ipointer=istep_cases(ipointer1)
c     ncases=istep(ipointer)
c     istep_cases(ipointer:1 thru +ncases)= case # given below.
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
         
         if(irray(inp1,0).lt.-100000000) then
            call missing_corner_velocities(inp1,xp1,yp1,
     1           zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2           vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)
         endif

      endif
c------------------------------------------------------------------

c interpolate velocities to particle position
      call get_velocities(inp1,np1,xp1,yp1,zp1,
     1     aomrx,aomry,aomrz,vx1b,vy1b,vz1b,vomrx,vomry,vomrz)

      vx1(np1)=vomrx
      vy1(np1)=vomry
      vz1(np1)=vomrz
      
      return
      
      end
      
c.....................................................................

      subroutine  zheng2(inp1,np1,rw,vx_zheng,vy_zheng,vz_zheng)

      use comai
      use combi
      use comflow
      use comdi
      use comsptr
      use comsk
      
      implicit none

      integer np1,icapture,nplanes,iplane(3),inp1,ipc,newnode
      integer inp_exit_zheng,indexrw

      real*8 xp1,yp1,zp1,rw,vx_zheng,vy_zheng,vz_zheng,dtnew
      real*8 xp2,yp2,zp2,xc,yc,zc,dtc,rwa
            
c.... well radius
c     indexrw=-(irray(inp1,0)+10000000)
c     if(indexrw.gt.0.and.indexrw.le.neq) then
c     rw=well_radius(indexrw)
c     else
c     rw = 0.1*min(ddx(inp1),ddy(inp1))
c     endif
      
c     Modify the flow field using
c     the semi-analytical capture solution from Zheng et al.
      xp1=x1(np1)
      yp1=y1(np1)
      zp1=z1(np1)
      
      call time_zheng_new(icnl,np1,
     1     xp1,yp1,zp1,rw,vx_zheng,vy_zheng,vz_zheng,dtnew) 
      
      call exit_zheng_new(icapture, icnl,np1,
     1     xp1,yp1,zp1,vx_zheng,vy_zheng,vz_zheng,rw,dtnew,
     2     xp2,yp2,zp2) 
      dtc=dtnew
      xc=xp2
      yc=yp2
      zc=zp2
      
      if(icapture.eq.1) then
         ipc=-10
         inp_exit_zheng=-inp1
         call update_exit(inp1,np1,ipc,inp_exit_zheng,dtc,
     1        xc,yc,zc)
c     call remove_particle(np1,5,exit_node(np1),ipc)
         exit_node(np1)=inp1
         call capture_particle(np1)
      else
         
         call find_exit_planes(inp1,np1,xp2,yp2,zp2,
     1        nplanes,iplane)
         
         if(nplanes.eq.0) then
            call update_samenode(inp1,np1,ipc,newnode,
     $           dtc,xc,yc,zc)
         else
            call exit_point_zheng(inp1,np1,xp1,yp1,zp1,
     1           xp2,yp2,zp2,nplanes,iplane,
     2           vx_zheng,vy_zheng,vz_zheng,rw,
     3           ipc,xc,yc,zc,dtc,newnode)
            
            if(newnode.gt.0) then
               
c     found a legitimate new neighbour, proceed with updating
c     location and indices
               
               call update_newnode(inp1,np1,ipc,newnode,
     $              dtc,xc,yc,zc)
               
            else
c     if newnode=0, particle has exited the model domain
               call update_exit(inp1,np1,ipc,newnode,
     $              dtc,xc,yc,zc)
            endif         
         endif
      endif
      
      return
      
      end
      
c...................................................................
      
      subroutine exit_spring(inp1,np1,ipc,newnode,dtc,xc,yc,zc)

      implicit none

      integer inp1,np1,ipc,newnode,inp_exit_spring

      real*8 dtc,xc,yc,zc

      ipc=-30
      inp_exit_spring=+inp1
      call update_exit(inp1,np1,ipc,inp_exit_spring,dtc,xc,yc,zc)
      call remove_particle(np1,5,inp1,ipc)

      return
      
      end

c....................................................................

      subroutine nearest_node_2(i,np1,ipc,level_max,xc,yc,zc, newnode)
c search all connecting nodes on the ipc side for the nearest
c neighbour to xc,yc,zc, and if cc is not found then search the
c neighbours of the neighbours. 
c This subroutine is modified
c from tree_search by adding the check if the node kb is on the
c correct side of i

      use comai, only : ierr, iptty 
      use combi, only : cord, nelm
      use comsptr
      use comsk
      
      implicit none

      integer i, newnode,i1,i2,ik,ipsign,ipab,kb,np1,ipc

      real*8 epsilon,diff,d,dmin,xc,yc,zc,xp,yp,zp
      real*8 xmin,xmax,ymin,ymax,zmin,zmax

      integer n_oldlist,n_newlist,n1,n2,index
      integer n_level,level_max, maxdim,flag_box,n_total

      epsilon = +1.e-18

c find global coordinates
      xp=xc+corn(i,1)
      yp=yc+corn(i,2)
      zp=zc+corn(i,3)
      ipsign=isign(1,ipc)
      ipab=iabs(ipc)

c in compress_list, n_level=0 if the node itself, =1 is the 
c connected neighbours and =2 is the neighbours of beighbours
c      level_max=2
      flag_box=0
      maxdim=10000
      if (.not. allocated(oldlist)) allocate(oldlist(maxdim))
      if (.not. allocated(newlist)) allocate(newlist(maxdim))
      
      n_oldlist=0
      n_newlist=1
      oldlist(1)=i
      flag_box=0
      
      do n_level=1,level_max
         
c     going to compress_list, n_oldlist is the #nodes 2 levels, back
c     and n_newlist is the # of nodes 1 level back. On return,
c     n_oldlist is the #nodes 1 level back, and n_newlist is the 
c     #nodes at the current level
         call compress_list(maxdim,n_oldlist,oldlist,n_newlist,
     1        newlist,n_level)
         
         do n1=1,n_newlist
            kb=newlist(n1)
            flag_box=0
            diff=ipsign*(cord(kb,ipab)-cord(i,ipab))
            if(diff.gt.epsilon) then
               call inside_box_2(kb,xp,yp,zp,flag_box)
               if(flag_box.gt.0) then
                  newnode=kb
                  goto 99999
               endif
            endif
         enddo
         
c     augment oldlist with newlist. Keeping the previous oldlist
c     to go 2 levels back in compress_list
         n_total=n_oldlist+n_newlist
         if(n_total.gt.maxdim) goto 99998
         do n1=1,n_newlist
            oldlist(n1+n_oldlist)=newlist(n1)
         enddo
         
      enddo
c xp,yp,zp not in the brick shaped approximate control volume
c of a neighbour of i up to 2 levels. Assume that the particle
c has exited the domain.

      newnode=0
      goto 99999
      
99998 write(ierr,*)'error in nearest_node_2.maxdim,n_total= ',
     1     maxdim,n_total
      write(ierr,*)'i,np1,ipc,xc,yc,zc',i,np1,ipc,xc,yc,zc
      write(ierr,*)' level_max=',level_max
      write(ierr,*)'STOP'
      write(iptty,*)'error in nearest_node_2.maxdim,n_total= ',
     1     maxdim,n_total
      write(iptty,*)'i,np1,ipc,xc,yc,zc',i,np1,ipc,xc,yc,zc
      write(iptty,*)' level_max=',level_max
      write(iptty,*)'STOP'
      stop
      
99999 continue

      return

      end subroutine nearest_node_2
      
c......................................................

      subroutine inside_box_2(kb,xp,yp,zp,flag_box)

c s kelkar Jan 21, 05
c check if the point xp,yp,zp is within the brick shape around the 
c node i. If so, flag_box=i, if not, flag_box=0

      use comsk
      use comsptr

      implicit none
      
      integer kb,flag_box

      real*8 xp,yp,zp

      flag_box=0

      if(corn(kb,1).le.xp) then
         if ((corn(kb,1)+ddx(kb)).ge.xp) then
            if(corn(kb,2).le.yp) then
               if ((corn(kb,2)+ddy(kb)).ge.yp) then
                  if(corn(kb,3).le.zp) then
                     if ((corn(kb,3)+ddz(kb)).ge.zp) then
                        flag_box=kb
                     endif
                  endif
               endif
            endif
         endif
      endif

      return

      end

c...................................................................
      subroutine newdt_courant(inp1,np1,aomrx,aomry,aomrz,
     1     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b,
     2     vomrx,vomry,vomrz,xp1,yp1,zp1,
     3     dtnew)
c do a simple estimate of dt for estimating the direction of
c the next step. if the particle reaches the boundary,
c then the actual time step is calculated in 
c pollock-intersect_plane.

      use comsk
      use comsptr

      implicit none

      integer inp1,np1,newnode,ipab,ipc,ipcx,ipcy,ipcz

      real*8 xvel,yvel,zvel,del,epsilon,dtnew,xdt,ydt,zdt
      real*8 dt_factor,aomrx,aomry,aomrz
      real*8 vx2b, vx1b,vy2b,vy1b,vz2b,vz1b,epfac,fac1
      real*8 dt_rand, zpp
      real*8 vomrx,vomry,vomrz,xp1,yp1,zp1,dzw

      dt_factor=1./3.

      epsilon=1.e-22
      epfac=1.e-8
      xdt=1.e+30
      ydt=1.e+30
      zdt=1.e+30
      
      if(vx1b.gt.+epsilon.and.vx2b.gt.+epsilon) then
         del=(ddx(inp1))*(1.+epfac)
         xvel=vx1b
         fac1=1.+aomrx*del/xvel
         if (abs(aomrx).gt.epsilon.and.fac1.gt.epsilon) then
            xdt=(dlog(fac1))/aomrx
         else
            xdt=(del/xvel)
         endif
      elseif(vx1b.lt.-epsilon.and.vx2b.lt.-epsilon) then
         del=-ddx(inp1)*(1.+epfac)
         xvel=vx2b
         fac1=1.+aomrx*del/xvel
         if (abs(aomrx).gt.epsilon.and.fac1.gt.epsilon) then
            xdt=(dlog(fac1))/aomrx
         else
            xdt=(del/xvel)
         endif
      else
c     velocity sign reversed on opposite faces. To be correct, we
c     should find the v=0 lplane, but it doesnt matter because this
c     time estimate is used only as a first guess, and correct time
c     is calculated as dtc later. So use the entire cell length for
c     both +ve and -ve V
         if(abs(vomrx).gt.epsilon) then
            xdt=abs(ddx(inp1)/vomrx)
         endif
      endif

      if(vy1b.gt.+epsilon.and.vy2b.gt.+epsilon)then
         del=(ddy(inp1))*(1.+epfac)
         yvel=vy1b
         fac1=1.+aomry*del/yvel
         if (abs(aomry).gt.epsilon.and.fac1.gt.epsilon) then
            ydt=(dlog(fac1))/aomry
         else
            ydt=+(del/yvel)
         endif
      elseif(vy2b.lt.-epsilon.and.vy1b.lt.-epsilon)then
         del=-ddy(inp1)*(1.+epfac)
         yvel=vy2b
         fac1=1.+aomry*del/yvel
         if(abs(aomry).gt.epsilon.and.fac1.gt.epsilon) then
            ydt=(dlog(fac1))/aomry
         else
            ydt=(del/yvel)
         endif
      else
         if(abs(vomry).gt.epsilon) then
            ydt=abs(ddy(inp1)/vomry)
         endif
      endif

      if(vz1b.gt.+epsilon.and.vz2b.gt.+epsilon) then
         del=(ddz(inp1))*(1.+epfac)
         zvel=vz1b
         fac1=1.+aomrz*del/zvel
         if(abs(aomrz).gt.epsilon.and.fac1.gt.epsilon) then
            zdt=(dlog(fac1))/aomrz
         else
            zdt=+(del/zvel)
         endif
      elseif(vz2b.lt.-epsilon.and.vz1b.lt.-epsilon) then
         del=-ddz(inp1)*(1.+epfac)
         zvel=vz2b
         fac1=1.+aomrz*del/zvel
         if(abs(aomrz).gt.epsilon.and.fac1.gt.epsilon) then
            zdt=(dlog(fac1))/aomrz
         else
            zdt=(del/zvel)
         endif
      else
         if(vomrz.lt.-epsilon.and.vz1b.lt.-epsilon) then
c     multiplying the distance from the Vz=0 level to the
c     lower bondary (index=-3) in order to avoid 0's and
c     infinities, because strictly speaking, a particle cant
c     leave a v=0 point. Then selecting 
c     dzw=0.9*abs(vz1b/aomrz), and calculating zdt from
c     zdt=abs(dlog(abs(aomrz*dzw/vz1b+1))/aomrz). 2.3 = ln(10.)
            if(abs(aomrz).gt.epsilon) then
               zdt=abs(2.3/aomrz)
            else
               zdt=1./epsilon
            endif
         elseif(vomrz.gt.+epsilon.and.vz2b.gt.+epsilon) then
c     multiplying the distance from the Vz=0 level to the
c     upper bondary (index=+3) in order to avoid 0's and
c     infinities, because strictly speaking, a particle cant
c     leave a v=0 point. Then selecting (arbitrarily)
c     z(particle)=zpp=dzw+.05*ddz, where dzw=abs(ddz-vz2b/aomrz),
c     and calculating zdt from
c     zdt=abs(dlog(abs(vz2b/(vz2b-aomrz*(ddz-zpp))))/aomrz)
            if(abs(aomrz).gt.epsilon) then
               zpp=abs(ddz(inp1)-vz2b/aomrz)+0.05*ddz(inp1)
               zdt=abs(dlog(abs(vz2b/(vz2b-aomrz*(ddz(inp1)-zpp))))/
     &              aomrz)
            else
               zdt=1./epsilon
            endif
         endif
      endif

c find the min of the nonzero times

      dtnew=1.e+30

      if((xdt.lt.ydt.and.xdt.gt.0.).or.ydt.le.0.) then
         if(xdt.lt.zdt.or.zdt.le.0.) then
            dtnew=xdt
c            ipc=ipcx
         else
            if(zdt.gt.0.) then
               dtnew=zdt
c               ipc=ipcz
            endif
         endif
      else
         if((ydt.lt.zdt.and.ydt.gt.0.).or.zdt.le.0.) then
            dtnew=ydt
c            ipc=ipcy
         else
            if(zdt.gt.0.) then
               dtnew=zdt
c               ipc=ipcz
            endif
         endif
      endif

c for random walk get time step based on Div(D)
      dt_rand=1.e+30
c      call dt_divd(inp1,dt_rand)

      dtnew = courant_factor*min(x61(np1),dtnew,dt_rand)

c     Reset dtnew to hit the end of the time step exactly if it 
c     will be in the current cell at that time. Set flag ioutt
c     so record this
      
      if(ioutt(np1).eq.0) then
         if(dtnew/86400..gt.tt1-ttp1(np1)) then
            dtnew=86400.*(tt1-ttp1(np1))
            ioutt(np1) = 1
         end if
      end if
      
c     Add check to make sure dt is positive
      
      if(dtnew.lt.0.) call pause_compute_exit
      dtnew = max(1.d-10,dtnew)
      
      if(dtnew.ge.1.e+30) then
         call remove_particle(np1,8,inp1,ipc_save(np1))
      endif
      
      return
      
      end subroutine newdt_courant

c......................................................................

      subroutine dt_divd(inp1,dt_rand)
c s kelkar June 15, 05
c calculate time step based on Div(D) for random walk problems

      use comdi
      use comsptr

      implicit none
     
      integer current_model,inp1,i

      real*8 dt_rand,dt_r(3),avdivd,epsilon

      epsilon=1.e-22
      dt_rand=1.e+30

      if(inp1.ne.0) then
         current_model = itrc(inp1)
         if(tprpflag(current_model).eq.2.or.
     $        tprpflag(current_model).eq.4) then
            do i=1,3
               dt_r(i)=1.e+30
               avdivd=abs(divd_omr(inp1,i)+
     $              dispersivity1(current_model)*dpordx_omr(inp1,i))
               if(avdivd.gt.epsilon) dt_r(i)=ddx(inp1)/avdivd
            enddo            
            dt_rand=min(dt_r(1),dt_r(2),dt_r(3))
            
         endif
      endif

      return

      end

c.....................................................................

      subroutine sub_vel_divd(flag_divd,i,id,vel_divd)
c  s kelkar June 21 05
c This sub corrects the velocities in ggg for Div(D) and grad(porosity)
c terms. 
c     we are using the bilinear interpolation of Darcy velocities in
c     calculating the Divergence of Dispersion tensor term
c     the gradient of diffusivity is being egnored for the present.

      use comdi
      use comsptr

      implicit none

      integer i,id,idir,flag_divd,idir_sign,current_model

      real*8 vel_divd

      idir_sign=isign(1,id)
      idir=abs(id)
      vel_divd=idir_sign*ggg(i,id)

      if(flag_divd.eq.1) then
         current_model= itrc(i) 
         if(divd_omr(i,1).le.-1.e+20) then
            call dispersion_divergence_omr(i)
         endif
         vel_divd=vel_divd+divd_weight*(divd_omr(i,idir))
      endif
      
      return

      end

c.....................................................................

      subroutine new_neighbour(i,np1,ipc,xc,yc,zc, newnode)
      
      use comsptr
      use comsk
      
      implicit none

      integer i, newnode,i1,i2,ik,ipsign,ipab,ibou,np1,ipc
      integer irray0

      real*8 xc,yc,zc
      real*8 epsilon,diff,d,dmin

      epsilon = +1.e-18

c....s kelkar  Jan 27march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular interior node, not a source/sink
c     irray(i,0) = -i-2000 : regular  node on a external boundary
c     irray(i,0) = -i : regular interior node that is a sink/source
c                        but not explicitly specified in sptr macro
c     -100000000 < irray(i,0) < -10000000 : regular interior node that is
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c     -200000000 < irray(i,0) < -100000000 : non-OMR cliff node 
c     irray(i,0) < -200000000 : OMR cliff node 
c                = -(i+2000000) : spring node
c                = -(i+1000000) : well-capture node on extrernal bound
c                   simillar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c
      irray0=irray(i,0)
      if(irray0.eq.0) then
c this an an OMR node not on an external boundary. in this case 
c newnode=0 indicates  non-unique neighbours on the ipc side.
c search all connecting nodes on the ipc side for the nearest
c neighbour to xc,yc,zc
         call nearest_node_3(i,np1,ipc,2,xc,yc,zc, newnode)
      elseif(irray0.eq.-(i+1000).or.
     3        irray0.lt.-200000000) then
c this is an OMR node on an external boundary, and perhaps
c a cliff node. check its boundary
c faces, stored in iboulist(), to see if the ipc plane is a
c boundary plane. If it is not, find the new nearest node.
         ibou=iboulist(i,7)
         do i1=1,ibou
            if(ipc.eq.iboulist(i,i1)) then
               newnode=0
               goto 9999
            endif
         enddo
         call nearest_node_3(i,np1,ipc,2,xc,yc,zc, newnode)
      else
c non-omr node
         if(irray0.eq.i.or.irray0.eq.-i.or.
     1        irray0.eq.-(i+2000000).or.
     2        irray0.eq.-(i+1000000).or.irray0.eq.-(i+2000).or.
     3        (irray0.lt.-(10000000).and.irray0.gt.-200000000)) 
     $        then
            newnode=irray(i,ipc)
         endif
      endif
      
 9999 return
      
      end subroutine new_neighbour

c...........................................................

      subroutine nearest_node(i,np1,ipc,xc,yc,zc, newnode)
c search all connecting nodes on the ipc side for the nearest
c neighbour to xc,yc,zc

      use comai, only : ierr, iptty
      use combi, only : cord, nelm
      use comdi, only : ps
      use comsptr
      use comsk
      
      implicit none

      integer i,newnode,i1,i2,ik,ipsign,ipab,kb,np1,ipc

      real*8 epsilon,diff,d,dmin,xc,yc,zc,xp,yp,zp
      real*8 xmin,xmax,ymin,ymax,zmin,zmax

      epsilon = +1.e-18

c find global coordinates
      xp=xc+corn(i,1)
      yp=yc+corn(i,2)
      zp=zc+corn(i,3)

      i1=nelm(i)+1
      i2=nelm(i+1)
      ipsign=isign(1,ipc)
      ipab=iabs(ipc)
      dmin=1.e+20
      do ik=i1,i2
         kb=nelm(ik)
         diff=ipsign*(cord(kb,ipab)-cord(i,ipab))

         if(diff.gt.epsilon) then
c the node is on the correct side of the face, check if the particle 
c in within its approximate brickshaped control volume. The distance 
c criteria doesnt always work- need to check actual cord
c            d=(xp-cord(kb,1))**2.+
c     &        (yp-cord(kb,2))**2.+
c     &        (zp-cord(kb,3))**2.
            xmin= corn(kb,1)
            ymin= corn(kb,2)
            zmin= corn(kb,3)
            xmax= xmin+ddx(kb)
            ymax= ymin+ddy(kb)
            zmax= zmin+ddz(kb)
            if (xmin.le.xp.and.xmax.ge.xp) then
               if (ymin.le.yp.and.ymax.ge.yp) then
                  if (zmin.le.zp.and.zmax.ge.zp) then
                     newnode=kb
                     if(ps(kb).le.0.) newnode=-kb
                     if(newnode.lt.0) then
c particle has entered a -ve or 0 porosity node. Not allowed
                        if (iptty .ne. 0)
     &                       write(iptty, 9999) i,np1,ipc,newnode
                        write(ierr, 9999) i,np1,ipc,newnode
                        istop(np1) = 10
                     endif
                  endif
               endif
            endif
         end if
      enddo
 9999 format ('STOP In nearest_node(compute_exit_new).', /,
     &     'particle has entered a -ve porosity node', /,
     &     'inp1,np1,ipc,newnode=', 4(1x, i8))
      return

      end subroutine nearest_node
      
c......................................................

      subroutine find_newpoint1(i,np1,ip1,
     &     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     ipc,xc,yc,zc,dtc)
c     s kelkar march 29 04
c     find exit point accross the plane given by ip
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip,ip1,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 dtc,dtc1,xc1,yc1,zc1
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      
      call pollock_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc1,yc1,zc1,dtc1)
      
      xc=xc1
      yc=yc1
      zc=zc1
      dtc=dtc1
      ipc=ip1
      
      return
      
      end subroutine find_newpoint1
      
c...........................................................
      
      subroutine find_newpoint2(i,np1,ip1,ip2,
     &     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     ipc,xc,yc,zc,dtc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1, and ip2 and taking the closest of the two
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,d1,d2,dc,dtc,dtc1,dtc2
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      
      call pollock_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc1,yc1,zc1,dtc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call pollock_intersect_plane(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc2,yc2,zc2,dtc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
         dtc=dtc1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
         dtc=dtc2
      endif
      
      return
      
      end subroutine find_newpoint2
      
c............................................................
      
      subroutine find_newpoint3(i,np1,ip1,ip2,ip3,
     &     xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     ipc,xc,yc,zc,dtc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1,ip2 and ip3 and taking the closest of the three
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,ip3,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,d1,d2,d3,dc
      real*8 dtc,dtc1,dtc2,dtc3
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      
      call pollock_intersect_plane(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc1,yc1,zc1,dtc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call pollock_intersect_plane(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc2,yc2,zc2,dtc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      call pollock_intersect_plane(i,np1,ip3,xp1,yp1,zp1,xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     &     xc3,yc3,zc3,dtc3)
      d3=sqrt((xp1-xc3)**2.+(yp1-yc3)**2.+(zp1-zc3)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
         dtc=dtc1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
         dtc=dtc2
      endif
      if(dc.gt.d3) then
         xc=xc3
         yc=yc3
         zc=zc3
         ipc=ip3
         dc=d3
         dtc=dtc3
      endif
      
      return
      
      end subroutine find_newpoint3
c......................................................................

      subroutine find_exit_point(inp1,np1,nplanes,iplane,xp1,yp1,zp1,
     1     xp2,yp2,zp2,
     2     vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3     ipc,xc,yc,zc,dtc,newnode,vexit)
      
c     find exit point(xc,yc,zc with reference to corn(inp1,)) 
c     and exit plane(ipc), and the new neighbour
c     xc is taken to be the intersection of the STRAIGHT line xp1-xp2
c     with the plane ipc- this is an approximation of the 
c     Pollock trajectory      

      use comai, only: ierr, iptty
      use comsptr
      use comsk
      
      implicit none
      
      integer inp1,np1,icnl,newnode
      integer iplane(3),nplanes,ipc
      
      real*8 xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc,dtc,vexit
      real*8 vomrx,vomry,vomrz,aomrx,aomry,aomrz
      
      newnode=0
      ipc=0
      
c     if nplanes > 0 then particle has entered a new cc
      if(nplanes.eq.1) then
         call find_newpoint1(inp1,np1,iplane(1),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     2        vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.2) then
         call find_newpoint2(inp1,np1,iplane(1),iplane(2),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     2        vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.3) then
         call find_newpoint3(inp1,np1,iplane(1),iplane(2),
     1        iplane(3), xp1,yp1,zp1,xp2,yp2,zp2,
     2        vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     3        ipc,xc,yc,zc,dtc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      endif
      
c     if there is a nonzero newnode and 
c     if the point lands on a face,
c     fix xc,yc,zc, to be slightly inside the new control volume 
      if(newnode.gt.0) then
         call shift_xc_inside(inp1,newnode,np1,
     2        vomrx,vomry,vomrz,ipc,xc,yc,zc)
      elseif(newnode.lt.0) then
c     particle has entered a -ve or 0 porosity node. Not allowed
         if (iptty .ne. 0) write(iptty, 9999) inp1,np1,ipc,newnode
         write(ierr, 9999)inp1,np1,ipc,newnode
         istop(np1) = 10
      endif
      
      return
      
 9999 format ('STOP In find_exit_point(compute_exit_new).', /,
     &     'particle has entered a -ve porosity node', /,
     &     'inp1,np1,ipc,newnode=', 4(1x, i8))
      end subroutine find_exit_point

c...........................................................


      subroutine wtsi_newnode(inp1,np1,newnode,xc,yc,zc,
     1     vomrz,aomrz,vz1b,vz2b)
      
c     s kelkar  Jul 25, 05
c     if the newnode is a water table node with a lower
c     saturation, drop the particle down so that zp2=
c     ddz(newnode)*s(newnode). If newnode has irreducible water
c     saturation, flagged by izone_free_nodes(inp1) = 3 
c     search vertically downward to find a node with 
c     flowing water.

      use comai, only : ierr, iptty 
      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,newnode,j,nextnode, node_flag, nodetemp
      integer node_previous,ibou,i1
      
      real*8vomrz,aomrz,xc,yc
      real*8 vz1b,vz2b
      real*8 dzw,zc, epsilon

      epsilon= 1.e-10
      
      if (izone_free_nodes(newnode).ge.3.or.s(newnode).lt.sirr) then
         node_previous=newnode
         do j=1,1000000
            node_flag=irray(node_previous,0)
            if(node_flag.eq.0) then
c     nextnode is interior OMR, do OMR check
c     xc and yc are updated wrt nodetemp in wtsi_neighbour and also
c     zc is set on the +3 boundary of nodetemp.
               call wtsi_neighbour(inp1,np1,node_previous,
     1              xc,yc,zc, nodetemp) 
               if(izone_free_nodes(nodetemp).le.2.and.s(nodetemp)
     1              .ge.sirr) then
c     found a valid flowing node, jump out of the do loop  and
c     update zc and z-velocity and return
                  newnode=nodetemp
                  call wtsi_newlocation(node_previous,np1,nodetemp,
     1                 vomrz,aomrz,vz1b,vz2b,zc)
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
                     goto 99999
                  endif
               enddo
c     node_previous not on -3 boundary. In wtsi_neighbour,
c     xc and yc are updated wrt nodetemp in wtsi_neighbour and also
c     zc is set on the +3 boundary of nodetemp.
               call wtsi_neighbour(inp1,np1,node_previous,
     1              xc,yc,zc, nodetemp) 
               if(izone_free_nodes(nodetemp).le.2.and.s(nodetemp)
     1              .ge.sirr) then
c     found a valid flowing node, jump out of the do loop  and
c     update particle location and velocity and return
                  newnode=nodetemp
                  call wtsi_newlocation(node_previous,np1,nodetemp,
     1                 vomrz,aomrz,vz1b,vz2b,zc)
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
                  xc=xc+corn(node_previous,1)-corn(nextnode,1)
                  yc=yc+corn(node_previous,2)-corn(nextnode,2)
                  zc=(1.-epsilon)*ddz(nextnode)
                  If(izone_free_nodes(nextnode).le.2.and.s(nextnode)
     1              .ge.sirr) then
c     found a valid flowing node, jump out of the do loop  and
c     update particle location and velocity and return
                     newnode=nextnode
                     call wtsi_newlocation(node_previous,np1,nextnode,
     1                    vomrz,aomrz,vz1b,vz2b,zc)
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
c     found a valid flowing node, jump out of the do loop  and
c     update particle location and velocity and return
         call wtsi_newlocation(inp1,np1,newnode,
     1        vomrz,aomrz,vz1b,vz2b,zc)
      endif
               
99999 continue
      
      return
      
      end subroutine wtsi_newnode

c...........................................................
      
      subroutine wtsi_neighbour(inp1,np1,node_previous,xc,yc,zc,
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
c     place zcc at -3 face of nodetemp.

      zcc=-epsilon*ddz(node_previous)

      call nearest_node(node_previous,np1,-3,xc,yc,zcc, nodetemp)

      if(nodetemp.le.0) then
         write(ierr,*)'Hole in the model. node=',node_previous
         write(ierr,*)'stop in wtsi_neighbour'
         write(iptty,*)'Hole in the model. node=',node_previous
         write(iptty,*)'stop in wtsi_neighbour'
         call update_exit(-inp1,np1,-100,nodetemp,
     $     0.d0,xc,yc,zc)
c         stop
      endif   

c     change xc,yc to refere to nodetemp 
c     place zc at +3 face of nodetemp
      xc=xc+corn(node_previous,1)-corn(nodetemp,1)
      yc=yc+corn(node_previous,2)-corn(nodetemp,2)
      zc=(1.-epsilon)*ddz(nodetemp)*s(nodetemp)

      return
      
      end subroutine wtsi_neighbour

c...........................................................

      subroutine wtsi_newlocation(nodeabove,np1,nodebelow,
     1     vomrz,aomrz,vz1b,vz2b,zc)

      use comai, only : day
      use comdi
      use comsptr
      use comsk

      implicit none

      integer nodeabove,np1,nodebelow,j
      
      real*8 vomrz,aomrz
      real*8 vz1b,vz2b
      real*8 dzw,zc, epsilon, vwt,snoise,ds

      epsilon=1.e-6
      snoise=0.0002

      dzw=s(nodebelow)*ddz(nodebelow)
      if((dzw+corn(nodebelow,3)).lt.(ddz(nodeabove)+corn(nodeabove,3)))
     1     then
c     at this point zc is assumed to be wrt nodebelow already.
         if(zc.gt.dzw) then
            zc=dzw*(1.-epsilon)
c     partially saturated water table node, set aomrz and vz2b 
c     so that Vz=ddz*(dSw/dt) at the water table. do not include
c     Del(D)-z terms
c*****Aug 22, 05 next line a quick fix hard wired to smooth out 
c*****the effect of water table flucutaions
            ds=s(nodebelow)-so(nodebelow)
            if(abs(ds).gt.snoise) then
               vwt= (ddz(nodebelow)*ds)/(day*86400.)
            else
               vwt= 0.
            endif
            vz1b=-ggg(nodebelow,-3)
            aomrz=(vwt-vz1b)/dzw
            vz2b=vz1b+(aomrz)*ddz(nodebelow)
            vomrz=vwt
         endif
      endif

      return
      
      end subroutine wtsi_newlocation
c...........................................................

      subroutine vel_zheng(inp1,icnl,x,y,z,dx,dy,dz,
     1     vx_zheng,vy_zheng,vz_zheng)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,icnl

      real*8 dx,dy,dz,skz,pi,x,y,z
      real*8 qwp,qqw,a,fac,xyw,xp1,yp1,zp1
      real*8 vx_zheng,vy_zheng,vz_zheng

      pi=3.14159
c convert sk from kg/sec to m^3/sec assuming rhow=1000kg/m^3
c To be consistant with Zheng's conventions, devide sk by dx*dy*dz and 
c change sign because fehm convention is -ve for injection, but
c Zheng convention is +ve for injection
      skz = -1.*sk(inp1)/(1000.)/(dx*dy*dz)/ps_trac(inp1)

c.......s kelkar march 30 2001 
c if irevers.eq.1 then the sign of ggg is changed in ptrac3. 
c hence it has to be changed again

      if(irevers.eq.1) then
         qwp=skz-(-ggg(inp1,+3)-ggg(inp1,-3))/dz
      else
         qwp=skz-(+ggg(inp1,+3)+ggg(inp1,-3))/dz
      endif

      qqw=qwp*dx*dy*dz
      a=1.
c      a=(permy/permx)
      fac=qqw*sqrt(a)/(2.*pi)/dz
      xyw=x*x/a+y*y
      
      if(irevers.eq.1) then
         vx_zheng=fac*x/xyw+0.5*(-ggg(inp1,1)+ggg(inp1,-1))
         vy_zheng=fac*y/xyw+0.5*(-ggg(inp1,2)+ggg(inp1,-2))
         if(icnl.eq.0) then
            vz_zheng=((-ggg(inp1,3)-ggg(inp1,-3))/dz)*(z-0.5*dz)
     1           +ggg(inp1,-3)
         endif
      else
         vx_zheng=fac*x/xyw+0.5*(+ggg(inp1,1)-ggg(inp1,-1))
         vy_zheng=fac*y/xyw+0.5*(+ggg(inp1,2)-ggg(inp1,-2))
         if(icnl.eq.0) then
            vz_zheng=((+ggg(inp1,3)+ggg(inp1,-3))/dz)*(z-0.5*dz)
     1           -ggg(inp1,-3)
         endif
      endif
         
      return

      end

c..................................................................

      subroutine missing_corner_velocities(i,xp1,yp1,
     1     zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz,
     2     vx1b,vy1b,vz1b,vx2b,vy2b,vz2b)

      use comai, only : ierr, iptty 
      use combi, only : cord, nelm
      use comsptr
      use comsk
      
      implicit none

      integer i,isf,ncases,ip,i1,i2,icase
      integer cliff_pointer_zero , irray0

      real*8 xp1,yp1,xpg,ypg,zpg
      real*8 zp1,vomrx,vomry,vomrz,aomrx,aomry,aomrz
      real*8 aomrcx,aomrcy,aomrcz
      real*8 vx1b,vy1b,vz1b,vx2b,vy2b,vz2b
      real*8 vxc1b,vyc1b,vzc1b,vxc2b,vyc2b,vzc2b

      cliff_pointer_zero = 100000000

c     -200000000 < irray(i,0) < -100000000 : non-OMR cliff node 
c     irray(i,0) < -200000000 : OMR cliff node
c if the node is a cliff node, but has a specified boundary
c outflow at it, remove cliff tag and mark as a regular
c boundry node in load_omr_flux_array

      vxc1b=vx1b
      vyc1b=vy1b
      vzc1b=vz1b
      vxc2b=vx2b
      vyc2b=vy2b
      vzc2b=vz2b

      xpg=xp1+corn(i,1)
      ypg=yp1+corn(i,2)
      zpg=zp1+corn(i,3)

      irray0=irray(i,0)
      if(irray0.lt.-2*cliff_pointer_zero) then
         ip=-irray0-2*cliff_pointer_zero
      elseif(irray0.lt.-cliff_pointer_zero) then
         ip=-irray0-cliff_pointer_zero 
      endif
      i1=istep_cases(ip)
      i2=istep_cases(ip+1)
      ncases=i2-i1

      do isf=i1,i2-1
         icase=istep_cases(isf)
         if(icase.eq.1) then
c     missing faces normal1=+x, normal2=+y
            if(ypg.gt.cord(i,2)) then
               vxc2b=-vxc1b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.gt.cord(i,1)) then
               vyc2b=-vyc1b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            endif
         elseif(icase.eq.2) then
c     missing faces normal1=+x, normal2=-y
            if(ypg.lt.cord(i,2)) then
               vxc2b=-vxc1b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.gt.cord(i,1)) then
               vyc1b=-vyc2b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            endif
         elseif(icase.eq.3) then
c     missing faces normal1=+x, normal2=+z
            if(zpg.gt.cord(i,3)) then
               vxc2b=-vxc1b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.gt.cord(i,1)) then
               vzc2b=-vzc1b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.4) then
c     missing faces normal1=+x, normal2=-z
            if(zpg.lt.cord(i,3)) then
               vxc2b=-vxc1b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.gt.cord(i,1)) then
               vzc1b=-vzc2b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.5) then
c     missing faces normal1=-x, normal2=+y
            if(ypg.gt.cord(i,2)) then
               vxc1b=-vxc2b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.lt.cord(i,1)) then
               vyc2b=-vyc1b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            endif
         elseif(icase.eq.6) then
c     missing faces normal1=-x, normal2=-y
            if(ypg.lt.cord(i,2)) then
               vxc1b=-vxc2b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.lt.cord(i,1)) then
               vyc1b=-vyc2b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            endif
         elseif(icase.eq.7) then
c     missing faces normal1=-x, normal2=+z
            if(zpg.gt.cord(i,3)) then
               vxc1b=-vxc2b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.lt.cord(i,1)) then
               vzc2b=-vzc1b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.8) then
c     missing faces normal1=-x, normal2=-z
            if(zpg.lt.cord(i,3)) then
               vxc1b=-vxc2b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
            elseif(xpg.lt.cord(i,1)) then
               vzc1b=-vzc2b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.9) then
c     missing faces normal1=+y, normal2=+z
            if(zpg.gt.cord(i,3)) then
               vyc2b=-vyc1b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            elseif(ypg.gt.cord(i,2)) then
               vzc2b=-vzc1b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.10) then
c     missing faces normal1=+y, normal2=-z
            if(zpg.lt.cord(i,3)) then
               vyc2b=-vyc1b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            elseif(ypg.gt.cord(i,2)) then
               vzc1b=-vzc2b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.11) then
c     missing faces normal1=-y, normal2=+z
            if(zpg.gt.cord(i,3)) then
               vyc1b=-vyc2b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            elseif(ypg.lt.cord(i,2)) then
               vzc2b=-vzc1b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         elseif(icase.eq.12) then
c     missing faces normal1=-y, normal2=-z
            if(zpg.lt.cord(i,3)) then
               vyc1b=-vyc2b
               aomrcy=(vyc2b-vyc1b)/ddy(i)
            elseif(ypg.lt.cord(i,2)) then
               vzc1b=-vzc2b
               aomrcz=(vzc2b-vzc1b)/ddz(i)
            endif
         endif
      enddo

c if the particle is being driven against a cliffface, move it
c into the lower half of the CC and reset velocities and
c accelarations

      if(vxc1b.gt.0..and.vxc2b.lt.0.) then
         if(vyc1b.gt.0..and.vyc2b.lt.0.) then
            if(zpg.gt.cord(i,3)) then
               vxc1b=vx1b
               vyc1b=vy1b
               vzc1b=vz1b
               vxc2b=vx2b
               vyc2b=vy2b
               vzc2b=vz2b
               aomrcx=(vxc2b-vxc1b)/ddx(i)
               aomrcy=(vyc2b-vyc1b)/ddy(i)
               aomrcz=(vzc2b-vzc1b)/ddz(i)
               zp1=0.25*ddz(i)
            endif
         endif
      endif

      vx1b=vxc1b
      vy1b=vyc1b
      vz1b=vzc1b
      vx2b=vxc2b
      vy2b=vyc2b
      vz2b=vzc2b
      aomrx=aomrcx
      aomry=aomrcy
      aomrz=aomrcz

      return

      end subroutine missing_corner_velocities
      
c......................................................

      subroutine fix_corner_position(i,xp1,yp1,zp1)

      use comsk
      use comsptr
      
      implicit none

      integer i,iii,j,ibou

      real*8 xp1,yp1,zp1

      ibou=iboulist(i,7)
      if(ibou.ge.2) then
         do iii=1,ibou
            if(iboulist(i,iii).eq.+3) then
c if the node i is on the top surface of the model 
c first fix the z-position
               if((ddz(i)-zp1).lt.epsilon_corner) then
                  zp1 = ddz(i)-epsilon_corner
               endif
c now find out which other faces make the corner and fix those
c positions
               do j=1,ibou
                  if(iboulist(i,j).eq.+1) then
                     if((ddx(i)-xp1).lt.epsilon_corner) then
                        xp1=ddx(i)-epsilon_corner
                     endif
                  elseif(iboulist(i,j).eq.-1) then
                     if(xp1.lt.epsilon_corner) then
                        xp1=epsilon_corner
                     endif
                  elseif(iboulist(i,j).eq.+2) then
                     if((ddy(i)-yp1).lt.epsilon_corner) then
                        yp1=ddy(i)-epsilon_corner
                     endif
                  elseif(iboulist(i,j).eq.-2) then
                     if(yp1.lt.epsilon_corner) then
                        yp1=epsilon_corner
                     endif
                  endif

                  goto 91919
               enddo   
            endif
         enddo
      endif
      
91919 continue

      return
      
      end

c.........................................................


      subroutine nearest_node_3(i,np1,ipc,level_max,xc,yc,zc, newnode)
c search all connecting nodes on the ipc side for the nearest
c neighbour to xc,yc,zc, and if cc is not found then search the
c neighbours of the neighbours. 
c This subroutine is modified
c from tree_search by adding the check if the node kb is on the
c correct side of i- this has changed- see below
c 3/14/06 s kelkar- doing a search in all directions, not just the
c half space, since in some rare cases with overlapping volumes,can
c get into the cc of a node which is on the a different side of i
c 3/15/06 since there can be multiple overlapping volumes, when
c possible we want to choose a node in the forward half space
c hence that search is completed first, the nodes from the backward
c half space are saved in array templist() and used only if the 
c forward search has failed

      use comai, only : ierr, iptty 
      use combi, only : cord, nelm
      use comsptr
      use comsk
      
      implicit none

      integer i, newnode,i1,i2,ik,ipsign,ipab,kb,np1,ipc
      integer n_oldlist,n_newlist,n1,n2,index
      integer n_level,level_max, maxdim,flag_box,n_total
      integer n_templist

      real*8 epsilon,diff,d,dmin,xc,yc,zc,xp,yp,zp
      real*8 xmin,xmax,ymin,ymax,zmin,zmax


      epsilon = +1.e-18

c find global coordinates
      xp=xc+corn(i,1)
      yp=yc+corn(i,2)
      zp=zc+corn(i,3)
      ipsign=isign(1,ipc)
      ipab=iabs(ipc)

c in compress_list, n_level=0 if the node itself, =1 is the 
c connected neighbours and =2 is the neighbours of beighbours
c      level_max=2
      flag_box=0
      maxdim=10000
      if (.not. allocated(templist)) allocate(templist(maxdim))
      if (.not. allocated(oldlist)) allocate(oldlist(maxdim))
      if (.not. allocated(newlist)) allocate(newlist(maxdim))
      
      n_oldlist=0
      n_newlist=1
      oldlist(1)=i
      flag_box=0
      n_templist=0

      do n_level=1,level_max
         
c     going to compress_list, n_oldlist is the #nodes 2 levels, back
c     and n_newlist is the # of nodes 1 level back. On return,
c     n_oldlist is the #nodes 1 level back, and n_newlist is the 
c     #nodes at the current level
         call compress_list(maxdim,n_oldlist,oldlist,n_newlist,
     1        newlist,n_level)
         
         do n1=1,n_newlist
            kb=newlist(n1)
            flag_box=0
            diff=ipsign*(cord(kb,ipab)-cord(i,ipab))
            if(diff.gt.epsilon) then
c kb in the forward half
               call inside_box_2(kb,xp,yp,zp,flag_box)
               if(flag_box.gt.0) then
                  newnode=kb
                  goto 99999
               endif
            else
c     kb in or at the backward halfspace,save for later use
c     if the forward search fails 
               n_templist=n_templist+1
               templist(n_templist)=kb
            endif
         enddo
         
c     augment oldlist with newlist. Keeping the previous oldlist
c     to go 2 levels back in compress_list
         n_total=n_oldlist+n_newlist
         if(n_total.gt.maxdim) goto 99998
         do n1=1,n_newlist
            oldlist(n1+n_oldlist)=newlist(n1)
         enddo
      enddo

c if we get to this point, the forward search has failed, 
c do a search on the backward half space using the saved nodes
c in templist()
      do n1=1,n_templist
         kb=templist(n1)
         flag_box=0
         call inside_box_2(kb,xp,yp,zp,flag_box)
         if(flag_box.gt.0) then
            newnode=kb
            goto 99999
         endif
      enddo

c if we get to this point both forward and backward searches have 
c failed. 
c xp,yp,zp not in the brick shaped approximate control volume
c of a neighbour of i up to 2 levels. Assume that the particle
c has exited the domain.

      newnode=0
      goto 99999
      
99998 write(ierr,*)'error in nearest_node_3.maxdim,n_total= ',
     1     maxdim,n_total
      write(ierr,*)'i,np1,ipc,xc,yc,zc',i,np1,ipc,xc,yc,zc
      write(ierr,*)' level_max=',level_max
      write(ierr,*)'STOP'
      write(iptty,*)'error in nearest_node_3.maxdim,n_total= ',
     1     maxdim,n_total
      write(iptty,*)'i,np1,ipc,xc,yc,zc',i,np1,ipc,xc,yc,zc
      write(iptty,*)' level_max=',level_max
      write(iptty,*)'STOP'
      stop
      
99999 continue

      return

      end subroutine nearest_node_3
      
c......................................................
