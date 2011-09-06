      subroutine ptrac3
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
!D1 Perform streamline particle tracking simulations at the current
!D1 time step.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/ptrac3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:18   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:30 1999   pvcs
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

      use comflow
      use davidi
      use comji
      use comhi
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsptr
      use comsk
      use compart

      implicit none

      real*8 small_slope
      parameter(small_slope = 1.e-12)
      integer nelmm,n50
      integer i3,i,ism,idone
      integer itime,istep,np1,ip,itloop
      real*8 xcoordw, ycoordw, zcoordw
      real*8 xface1,xface2,xface3
      real*8 ddt,edt,face1,face2,face3
      real*8 xxx,tto
      real*8 ep,ep5,x60
      real*8 days_save, tt1_old
      real*8 x2nondim, y2nondim, z2nondim
      real*8 x2dim, y2dim, z2dim
      integer walk_counter
      real*8 dt_delay, denom
      real*8 sigmav,omegav,conc_ret,ret_factor,par3v
      integer fm
      integer current_node, current_model, upcoming_node
      real*8 ps_print
      real*8 s_print
      integer izone, current_zone, icnl_subst
      integer lower_limit, upper_limit
      real*8 local_pos, local_spacing
      real*8 edtx, edty, edtz
      real*8 x2save, y2save, z2save, v_x, v_y, v_z, vnew, vcurrent
      real*8 ddxvsave, ddyvsave, ddzvsave
      real*8 x_ts, y_ts, z_ts
      integer ijkvsave, ijkvssave
      logical done, inside
      integer ntry_pos
      logical too_many_jumps
      character*4 summstring
      integer iwrite
      real*8  xwrite

c...... s kelkar 8/26/99
      integer ix,iy,iz
      real*8 xp,yp,zp
      real*8 vjx,vjy,vjz
      integer ijkvnp
      integer iskprt,isotropic_flag,leave_all
      integer iprint, iprcount
      real*8 tsubtract
      logical ex, totflag, btc_done, path_done, time_reset
      integer time_step_flag
c.........

      integer inp_adv

      real*8 vx_zheng,vy_zheng,vz_zheng,rw
c

      integer idebug,idum_count
      integer count_sptr2,inc_sptr2

      save count_sptr2,inc_sptr2

c......................

c	Definition of time variables

c	ddt - heat and mass transfer time step that the particle
c		tracking simulation is to proceed to (s)
c	dtmx - maximum time step for a particle to travel (s)
c	itm - maximum number of time steps to accomplish the
c		current time step
c	edt - current time step that the particles are taking (s)
c	tto - current time since entering the subroutine (s)
c		i.e. the time of the particle tracking simulation
c		where time 0 is the previous heat and mass transfer
c		time. The clock resets to 0 each time entering ptrac3
c	tt1_old - previous time of particle tracking solution inside
c		the time step loop (days)
c	tt1 - same as tto but in days instead of seconds
c       ttp1 - for each particle, the current time that the particle
c               has reached (days)
c	dt - for each particle, the current time step for computing
c		the next particle location (s)
c	dtx - for each particle, the time required to reach the next
c	cell face if it is in the x-direction (s)
c	dty - for each particle, the time required to reach the next
c	cell face if it is in the y-direction (s)
c	dtz - for each particle, the time required to reach the next
c	cell face if it is in the z-direction (s)

c....s kelkar  march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular node, not a source/sink
c     irray(i,0) = -i-2000 : regular node on a external boundary
c     irray(i,0) = -i : regular interior node that is a sink/source
c                        but not explicitly specified in sptr macro
c     irray(i,0) < -10000000 : regular interior node that is 
c                        specified as a sink/source in the sptr macro
c               in this case -(irray(i,0)+10000000) is the pointer
c                for the storage location in well_radius for this node
c                = -(i+2000000) : spring node
c                = -(i+1000000) : well-capture node on external bound
c                   similar to =-i case but with half space solution
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c.................................................................

      idum_count=0
      path_done = .false.
      count_sptr2=0
      inc_sptr2=0
      btc_done = .false.
      totflag = .false.
      if (.not. compute_flow .or. l .eq. 1) time_reset = .false.

c*******************************************************************
c s kelkar april 1 05 for debugging the random walk
c      omr_flag=.true.
c      open(unit=96,file='rand.dat')
c***************************************************************

c...s kelkar  nov 13 02..............................
c if freez_time is gt.0, then ptrac3 is called only at the end of the 
c flow calculations, and for a velocity frozen at that time,  
c particle tracks are calculated for freez_time days
      if(freez_time.eq.0.) then
         ddt=day*86400.
      else
         ddt=freez_time*86400.
      endif
c....................................................

      nelmm=neq
      n50=100000
      ep=1.e-7
      ep5=.5*ep
      x60=log(1./ep)
      time_step_flag=0
      if(dtmx.lt.0.) then
         time_step_flag=+1
         dtmx=abs(dtmx)
      endif

c............ skelkar 10/6/99
c note here that if iprt<0, the snapshot output is not done and
c the integrated concentration is output in the .trc file
c if iptr >0 then a snapshot output is made for every iprt iterations

      iskprt=abs(iprt)

      if(.not.unstruct_sptr .and. compute_flow) then

c
c call routine to load a_axy array into ggg array
c GAZ 5/12/98
c
         call load_omr_flux_array

      else

         call load_unstruct_flux_array

      end if
      
      if (omr_flag) then      
c initialize a bogus value in divd_omr and dpordx_omr to use as 
c flag if the node i has seen a particle or not
         do i=1,neq
            divd_omr(i,1)=-1.1e+20
            dpordx_omr(i,1)=-1.1e+20
         enddo
      endif

      idone=0
      if (time_btc .gt. 0.) then
         edt = min (dtmx, time_btc*86400.)
      else         
         edt = dtmn
      end if
      if(ddt.lt.edt) edt=ddt
      tto=0

      idebug=0

c......s kelkar 5/2/2001 snapshot output
         
      if(iskprt.eq.iprt.and.iprt.ne.0) then
         itime = 0
         if (idebug .ne. 0 .and. iptty .ne. 0) 
     &        write(iptty,*) 'Calling snapshot', tt1,itime
         call write_snapshot
      end if
c............................

      do itime=1,itm
c     Store old time
         tt1_old = tt1
         if (totflag .and. nzbtc .ne. 0) edt = min (edt*part_mult, dtmx)
         tto=tto+edt
c     *** see if this is last time step***
         if(tto.ge.ddt) then
            idone=1
            edt=ddt-tto+edt
         endif
         
         tt1=tt1+edt/86400.
         if (tt1_old .lt. time_btc .and. tt1 .gt. time_btc) then
            tt1 = time_btc
         end if
                  
c......s kelkar  7/23/01
c     Write to output file for plume
      if(nplum.gt.0) then
         call write_plum_values
      end if
c.............................
         
c     Set flag so that particles that are ahead of the current
c     time due to sorption or matrix diffusion are not advanced
c     or dispersed
         
         do np1 = 1, num_part
            if(tt1.ge.ttp1(np1)) then
               ioutt(np1) = 0
c If this is a wtsi problem, and the particle is in a partially 
c saturated or unsaturated node, move it down to the water table
c izone_free_nodes(i) = 0, non-wtsi node
c izone_free_nodes(i) = 1, fully saturated wtsi node
c izone_free_nodes(i) = 2, partially saturated wtsi node
c izone_free_nodes(i) = 3, fully unsaturated wtsi node
               if (ifree .ne. 0 .and. istop(np1) .eq. 0 .and.
     &              ijkv(np1) .ne. 0 ) then
                  if (izone_free_nodes(ijkv(np1)).gt.1) then
                     current_node = ijkv(np1)
                     call wtsi_ptrac3(np1)
c Output new location if the particle has been moved
                     if (ioutt(np1) .ne. -1 .and. iprto .ne. 0) 
     &                    call write_path_info
                  end if
               end if
               if (.not. pstart_out(np1) .and. ioutt(np1) .ne. -1) then
                  pstart_out(np1) = .true.
                  if (iprto .ne. 0) call write_path_info
               end if                                
            else
               ioutt(np1) = -1
            end if
         end do

         
c     *** do steps to accomplish time interval edt *****
         do istep=1,n50
            
            
c     Convert to loop over all particles BAR 1-13-99
            
            do np1 = 1, num_part
               
               if(ioutt(np1).ge.1) ioutt(np1) = ioutt(np1) + 1
               
c     Only do computations for particles that are not already
c     ahead of the time step clock due to sorption or matrix
c     diffusion or already out of the system
c zvd 10-Sep-07 added check for particles in dry cells that shouldn't be moved               
               if(tt1.gt.ttp1(np1) .and. istop(np1).eq.0
     &              .and. ioutt(np1).ne.-1) then

                  current_node = ijkv(np1)
                  current_model = itrc(current_node)
c zvd 13-May-09 Skip the following if dispersion is not invoked
                  if(tprpflag(current_model).eq.2.or.
     $                 tprpflag(current_model).eq.4.or.
     $                 tprpflag(current_model).eq.13.or.
     $                 tprpflag(current_model).eq.14) then
                     if (.not. omr_flag) then
c...s kelkar  april 7 04............................................
c the following are for dispersion/random-walk and not fixed for OMR
c hence skipped if omr flag is "TRUE"
cc...8/29/01 s kelkar ...................
cc find which quadrant of the control volume corresponding
cc to current_node does the particle reside in
cc note that x2,y2,z2=x1,y1,z1 at this stage
                        call find_quadrant(np1,current_node,ix,iy,iz,
     $                       xp,yp,zp, vjx,vjy,vjz)

cc       divd is calculated in dispersion_divergence 
cc       and grad(ps) in porosity_gradient
                        call dispersion_divergence(current_node,ix,iy,
     &                       iz,xp,yp,zp)
                        call porosity_gradient_log(current_node,ix,iy,
     &                       iz)

                     else
c for OMR grids, use least squares interpolation routine to estimate
c velocity derivatives. s kelkar march 30,05
c do this calculation only if the node has not seen a particle before
c use divd_omr(node,1)=-1.e+20 as a flag
                        if(divd_omr(current_node,1).le.-1.e+20) then
                           call dispersion_divergence_omr(current_node)
c 10/26/2006 zvd Remove call until correctly implemented
c                        call porosity_gradient_omr(current_node)
                        endif
                     end if
                  end if
c     Set velocities and time steps needed to perform calculations
                  
                  if(.not.unstruct_sptr) then
                     
                     call set_velocities_struct
                     call set_ts_struct
                     
                  else
                     
                     call set_velocities_unstruct
                     call set_ts_unstruct
                     
                  end if
                  
               end if
               
            end do 

            do np1 = 1, num_part
               if (axyzm(np1) .gt. 1.e-28) then
                  x61(np1)=x60/axyzm(np1)
                  if (x61(np1) .eq. 0.d0) then
c If our time step for this particle is 0 stop it
                     write (ierr,*) 'Particle ', np1, ' has 0 time'
                     istop(np1) = 10
                  end if
               else
                  x61(np1)=1.e28
               end if
            end do
                        
c     Convert to loop over all particles BAR 1-13-99
            
            do np1 = 1, num_part
               
c     Only do computations for particles that are not already
c     ahead of the time step clock due to sorption or matrix
c     diffusion
               
               if(tt1.gt.ttp1(np1).and.istop(np1).eq.0
     &              .and. ioutt(np1).ne.-1) then
                  
                  current_node = ijkv(np1)
                  current_model = itrc(current_node)
c zvd 13-May-09 Skip the following if dispersion is not invoked
                  if(tprpflag(current_model).eq.2.or.
     $                 tprpflag(current_model).eq.4.or.
     $                 tprpflag(current_model).eq.13.or.
     $                 tprpflag(current_model).eq.14) then
                     if (.not. omr_flag) then
c...s kelkar  april 7 04............................................
c the following are for dispersion/random-walk and not fixed for OMR
c hence skipped if omr flag is "TRUE"
c*** 4/16/02 these calls are being repeated, need to store this info
cc...8/29/01 s kelkar ...................
cc find which quadrant of the control volume corresponding
cc to current_node does the particle reside in
cc note that x2,y2,z2=x1,y1,z1 at this stage
                        call find_quadrant(np1,current_node,ix,iy,iz,
     $                       xp,yp,zp, vjx,vjy,vjz)

cc       divd is calculated in dispersion_divergence 
cc       and grad(ps) in porosity_gradient
                        call dispersion_divergence(current_node,ix,iy,
     &                    iz,xp,yp,zp)
                        call porosity_gradient_log(current_node,ix,iy,
     &                       iz)
                     else
c for OMR grids, use least squares interpolation routine to estimate
c velocity derivatives. s kelkar march 30,05
                        if(divd_omr(current_node,1).le.-1.e+20) then
                           call dispersion_divergence_omr(current_node)
                           call porosity_gradient_omr(current_node)
                        endif
                     end if
                  end if
cc.............................................
cc****
c 9/16/04 s kelkar time step calculations moved to compute_exit_new
c                    if(.not.unstruct_sptr) then
c                     call transtime_struct
c
c...  s kelkar  9/6/02
c                     if(time_step_flag.eq.+1) then
c                        dt(np1)=60.
c                     endif
c.........................................
c ..s kelkar  march 25, 04......3DOMR
                  call compute_exit_new
     $                 (np1,ep5,ep,small_slope,idum_count,
     $                 time_step_flag,inp_adv,
     $                 vx_zheng,vy_zheng,vz_zheng)
                  if(istop(np1).eq.1) then
                     tt(np1) = tt(np1) + dt(np1)
                     ttt(np1) = ttt(np1) + dt(np1)/86400.
                     ttp1(np1) = ttp1(np1) + dt(np1)/86400.
                     goto 99999
                  endif
c...s kelkar 9/17/04 next 4 lines not needed......................
c                  else
c                     call transtime_unstruct
c                     call compute_exit_unstruct
c                  end if
                  
c...  s kelkar  feb 7 2001
c skip rest of the calculations if captured by a pumping node
c in compute_exit_struct
                  if(istop(np1).ge.1) then
c..s kelkar march 1 02............................
                     call write_particle_exit(np1,isptr4,itime)
c.................................................
                     goto 99999
                  endif
c......................................
                  
c     Internal subroutine to determine if particle has previously
c     or just now reached a btc zone
                  
!                  call seek_btczone(np1)

c......... s kelkar  jul23 01
c     Internal subroutine to determine if particle has previously
c     or just now reached a plume  zone

                  call seek_plumzone(np1)
c................
                  
                  if(istop(np1).eq.1.or.ioutt(np1).ge.2) then
                     if(ijkv(np1).eq.0) then
                        current_node = ijkvs(np1)
                        upcoming_node = 0
                     else
                        ijkv(np1) = ijkvs(np1)
                        current_node = ijkvs(np1)
                        upcoming_node = ijkv(np1)
                     end if
                  else
                     x1(np1) = x2(np1)
                     y1(np1) = y2(np1)
                     z1(np1) = z2(np1)
                     tt(np1) = tt(np1) + dt(np1)
                     ttt(np1) = ttt(np1) + dt(np1)/86400.
                     
c     Compute particle delay during this time step and
c     add to current time of the particle (ttp1)
                     
                     tsubtract = tt(np1)-dt(np1)
                     call compute_part_delay(0)
                     
c perform random 
c (3/14/06 s kelkar) but only if using a random walk model 
                     current_model = itrc(current_node)
                     if(tprpflag(current_model).eq.2.or.
     2                    tprpflag(current_model).eq.4.or.
     2                    tprpflag(current_model).eq.13.or.
     $                    tprpflag(current_model).eq.14) then
                        
                        if(.not. omr_flag) then
                           if(.not.unstruct_sptr) then
                              call new_part_loc_struct
                           else
                              call new_part_loc_unstruct
                           end if
                        else
                           call new_loc_omr(itime,istep,np1,edt,
     $                          x_ts,y_ts,z_ts,time_step_flag,
     $                          vx_zheng,vy_zheng,vz_zheng,rw,inp_adv)
                           if(istop(np1).eq.1) goto 99999
                        end if
c     need to update x1,y1,z1 again bcs x2,y2,z2 might have changed in
c     random_walk
                        x1(np1) = x2(np1)
                        y1(np1) = y2(np1)
                        z1(np1) = z2(np1)
                     endif
c.........................................................

                     x3(np1) = x2(np1) + corn(ijkv(np1), 1)
                     y3(np1) = y2(np1) + corn(ijkv(np1), 2)
                     z3(np1) = z2(np1) + corn(ijkv(np1), 3)
                     
c...  s kelkar  feb 7 2001
c     skip rest of the calculations if captured by a pumping node
c     in new_part_loc
                     if(istop(np1).ge.1) then 
c..s kelkar march 1 02............................
                        call write_particle_exit(np1,isptr4,itime)
c.................................................
                        goto 99999
                     endif
                     
c     Do the time delay if a particle has just
c     bounced out of the cell, but not if the particle
c     just got to an interface above and has bounced 
c     into a new cell from the interface
                     
                     if(tt(np1).ne.0) then
                        tsubtract = tt(np1)
                        call compute_part_delay(1)
                     end if
                                       
c     As a final check, make sure that no particle is outside the cell
c     If it is, change the coordinate in error and place the particle in
c     the middle of the current cell
                     
                     if(.not.unstruct_sptr) then
                        call check_part_loc_struct
                     else
                        call check_part_loc_unstruct
                     end if
                     
                     
                     
c     ******* state change done *********
c     **** state is ijkv,x1,y1,z1 ***********
                     
                  end if
                  
c     ******* state change done *********
c     **** state is ijkv,x1,y1,z1 ***********
                  
c     **** set istop=1 if point out of domain****
                  
                  if(ijkv(np1).le.0) then
                     istop(np1) = 1

c..s kelkar march 1 02............................
                     call write_particle_exit(np1,isptr4,itime)
c.................................................

c                  elseif(irray(ijkv(np1),0).lt.0) then
c                     continue
c     istop(np1) = 1
c     ijkv(np1) = 0
                  elseif(ijkv(np1).gt.0.and.istop(np1).ne.3) then
                     istop(np1) = 0
                  end if
c     End of if block for determining if the particle is ahead
c     of the clock due to sorption or matrix diffusion
                  
c...  s kelkar  feb 7 2001
99999             continue
c.....................

c     Output of particle information changed: BAR 6-15-99
                  inc_sptr2=inc_sptr2+1
                  if(iprto.ne.0) then
c zvd 14-Aug-07
                     if (iprto .lt. 0) then
c For plumecalc output we always need to check if the particle is
c leaving its current cell
                        call write_path_info
                     else if(linelim.lt.0.and.inc_sptr2.eq.abs(linelim))
     &                        then
c Output only every inc_sptr2 th line
                        call write_path_info
                        inc_sptr2=0
                     elseif(linelim.eq.0) then
                        call write_path_info
                     elseif(linelim.gt.0) then
c Limit total number of lines in sptr2 output file
                        if(count_sptr2.gt.abs(linelim)) then
                           write(iptty,*)' exiting ptrac3 because  '
                           write(iptty,*)'count_sptr2 is gt linelim.'
                           write(ierr,*)' exiting ptrac3 because '
                           write(ierr,*)'count_sptr2 is gt linelim.'
                           goto 102
                        else
                           call write_path_info
                        endif
                     end if
                  endif
                  
                  
               end if
! Called from compute_exit_new, etc. now
!               call seek_btczone(np1)
c     End of loop over all particles               
            end do
            
c..debugging s kelkar 4/30/02
c            call write_snapshot
c            write(isptr1,*)itime,istep
c.....................................

c            write(*,*)itime,istep

c exit the outer time loop, output data na dexit ptrac3
c if the state of the particles is such that 
c every particle has either reached a boundary(ijkv(ip)=0) or a 
c pumping node(istop(ip)=2) or a stuck node (istop=3)
            leave_all=1
            do ip=1,num_part
!               if(istop(ip).ne.2.and.istop(ip).ne.3.and.ijkv(ip).ne.0) 
!     $              leave_all=0
               if (istop(ip) .eq. 0 .and. ijkv(ip) .ne. 0) leave_all=0
            enddo
            if(leave_all.eq.1) go to 102
            
c if some particles are still in the system (istop=0) and
c if all of these inner particles have gotten to edt(ioutt(ip).ne.0)
c then Exit the istep-loop for this time step edt and commence the next
c time step in the itime-loop 
            itloop=1
            do ip=1,num_part
               if(tt1.ge.ttp1(ip)) then
                  if((istop(ip).eq.0).and.(ioutt(ip).eq.0)) itloop=0
               end if
            enddo
            if(itloop.eq.1) go to 101

         enddo                  ! end step loop

 101     continue
         
c     compute concentrations
c not doing concentrations for now 5/10/04 s kelkar 3DOMR
         if (.not. omr_flag) call compute_concentrations
         
c......s kelkar 9/17 99 snapshot output
         
         if(iskprt.eq.iprt.and.iprt.ne.0) then
            iskprt = 0
            call write_snapshot
         end if
         
         iskprt=iskprt+1
c......................................
c     Write to output file of breakthrough curve
         
         if(.not. alt_btc .and. nzbtc.gt.0) then
            if (.not. btc_done .and. time_btc .le. tt1) then
               if (.not. time_reset) then
                  edt = dtmn
                  time_reset = .true.
               end if
               if (btc_flag) call write_btc_values
            end if
         end if

         if(idone.eq.1) go to 102

      enddo                     ! end time loop

 102  continue
            
c......s kelkar 9/17 99 snapshot output on exit
      if(iprt.gt.0) then
         call write_snapshot
      end if

c......................................
c     Write to output file of breakthrough curve one last time
         
         if (nzbtc.gt.0) then
            if (.not. btc_done) then
               if (.not. alt_btc .and. btc_flag) call write_btc_values
            end if
         end if

! Write final locations of particles still in the system for restart
      if (sptr_flag) call write_sptrs

! Modified by zvd 14-Sep-05, num_part is now ouput at the start of 
! the file when using the minimal output options
! Added by BAR for minimal output
!      if(iprto .eq. -1) then
!         summstring = 'summ'
!         write(isptr2,'(a4)') summstring
!         write(isptr2,*) num_part
!      elseif(iprto .eq. -2 .or. iprto .eq. -3) then
c     write 0's to denote end of particle tracking info
!         iwrite = 0
!         xwrite = 0.
!         write(isptr2) iwrite, xwrite, iwrite
!         write(isptr2) num_part
c     Only output final location for normal sptr2, not for minimal output
c      if (iprto .ne. 0) then
      if (iprto .gt. 0) then
         path_done = .true.
         do np1 = 1, num_part
            current_node = abs(ijkv(np1))
            call write_path_info
         end do
      end if

c     close(96)

      return

      contains

******************************************************************
******************************************************************

      subroutine massflux
      implicit none
      real*8  flow_out, flow_in
      integer i2
      integer add_fact, lbox, neq1, ntmp, ntmp1, f_m_box
        do i=1,n0
	  flow_in = 0.
          flow_out = 0.
          add_fact = 0
	  lbox=i
          neq1=neq+1
	  if(i.gt.neq)then
            lbox=i-neq
	    add_fact=nelm(neq1)-neq1
	  endif
          ntmp=nelm(lbox)-neq+add_fact 
          ntmp1=nelm(lbox+1)-neq1+add_fact 
          f_m_box=(nelm(neq1)-neq1)*2+lbox
          do i2=ntmp,ntmp1
            if(a_axy(i2).gt.0) then
              flow_out = flow_out+a_axy(i2)
            else
              flow_in = flow_in + a_axy(i2)
            end if
          enddo
C                                           take care of flow between 
C                                           fracture and matrix in dpdp model
          if(idpdp.ne.0) then
            if(add_fact.gt.0) then
              if(a_axy(f_m_box).gt.0) then
                flow_in=flow_in-a_axy(f_m_box)
              else
                flow_out=flow_out-a_axy(f_m_box)
              endif
            else
              if(a_axy(f_m_box).gt.0) then
                flow_out=flow_out+a_axy(f_m_box)
              else
                flow_in=flow_in+a_axy(f_m_box)
              endif
            endif
          endif

         
	  if( flow_out .ne. 0. ) then        
            if( flow_out .ge. -flow_in ) then
              flow_ot(i)= flow_out
            else
              flow_ot(i)= -flow_in
            end if
          else
            flow_ot(i)= 0.
          end if
       end do
      return

      end subroutine massflux

******************************************************************
******************************************************************
******************************************************************

      subroutine write_path_info
      implicit none

      real*8 dist, sptr_time, idum
      integer position_in_string, upcoming_node
      character*200 sptr_prop_values

      upcoming_node = abs(ijkv(np1))
      sptr_time = ttp1(np1)
      sptr_prop_values = ''

c     Since x1, y1, and z1 have been updated, corn has to
c     use the upcoming or new node to get the correct coordinate
      if(upcoming_node.ne.0) then
         xcoordw = x1(np1) + corn(upcoming_node,1)
         ycoordw = y1(np1) + corn(upcoming_node,2)
         if(icnl.eq.0) then
            zcoordw = z1(np1) + corn(upcoming_node,3)
         else
            zcoordw = 0.
         end if
      else
         xcoordw = x3(np1)
         ycoordw = y3(np1)
         zcoordw = z3(np1)
      end if
c     tmp for debug
      if (part_id(np1,1) .eq. 2064) then
         idum = 1
      end if

c
! If this time and location have already been output, return
      if (path_done .and. ttpo(np1) .eq. ttp1(np1) .and.
     &     xcoordw .eq. xo(np1) .and. ycoordw .eq. yo(np1)
     &     .and. zcoordw .eq. zo(np1) .and. (current_node .eq. 
     &     lastnode(np1) .or. current_node .eq. 0)) return
      
! If we are ahead of the current time don't output location unless
! we are at the end of the simulation and output is requested
      if (sptr_time .gt. tims .and. .not. output_end) return
      if (sptr_time .gt. days .and. days .lt. tims) return

      if (iprto .eq. -1) then
         if(upcoming_node.ne.current_node) then
            if (current_node .eq. lastnode(np1)) then
! We shouldn't be here because we should have moved somewhere else during the previous step
               write (ierr, *) 
     &              'Particle exiting same node as last output'
               write (ierr, 9002)part_id(np1,1),ttpo(np1),lastnode(np1),
     &              xo(np1), yo(np1), zo(np1)
               write (ierr, 9002) part_id(np1,1),sptr_time,current_node,
     &              xcoordw, ycoordw, zcoordw
            end if
            if (xyz_flag) then
               write(isptr2,9002) part_id(np1,1),sptr_time,current_node,
     &              xcoordw, ycoordw, zcoordw    
            else 
               write(isptr2,9001) part_id(np1,1),sptr_time,current_node
            end if
! Save values just written 
            count_sptr2=count_sptr2+1
            ttpo(np1) = ttp1(np1)
            xo(np1) = xcoordw
            yo(np1) = ycoordw
            zo(np1) = zcoordw
            lastnode(np1) = current_node
            call flush (isptr2)
         else if (path_done .and. xyz_flag) then
            write(isptr2,9002) part_id(np1,1),sptr_time,current_node,
     &              xcoordw, ycoordw, zcoordw    
         end if
      else if (iprto.eq.-2 .or. iprto.eq.-3) then
         if(upcoming_node.ne.current_node) then
            if (xyz_flag) then
               write(isptr2) part_id(np1,1),sptr_time,current_node,
     &              xcoordw, ycoordw, zcoordw 
            else   
               write(isptr2) part_id(np1,1),sptr_time,current_node
            end if
! Save values just written 
            count_sptr2=count_sptr2+1
            ttpo(np1) = ttp1(np1)
            xo(np1) = xcoordw
            yo(np1) = ycoordw
            zo(np1) = zcoordw
            lastnode(np1) = current_node
            call flush (isptr2)
         else if (path_done .and. xyz_flag) then
               write(isptr2) part_id(np1,1),sptr_time,current_node,
     &              xcoordw, ycoordw, zcoordw 
         end if
      else

         
! Check to see if particle has moved enough to be output
         dist = dsqrt ((xcoordw - xo(np1))**2 +  
     &        (ycoordw - yo(np1))**2 + (zcoordw - zo(np1))**2)
         if (path_done .or. istop(np1) .ne. 0) then
! If this is the last call or the particle is stopped 
! (but hasn't exited) we want to print the location
         else
! Check write control parameters
            if (tplim .ne. 0. .and. distlim .ne. 0.) then
! Use both time and distance limit to see if particle has moved far
! enough and long enough
               if ((ttp1(np1) - ttpo(np1)) .lt. tplim .and. 
     &              dist .lt. distlim) return
            else if (tplim .eq. 0. .and. distlim .ne. 0.) then
! Use distance limit to see if article has moved far enough
               if (dist .lt. distlim) return
            else if (tplim .ne. 0. .and. distlim .eq. 0.) then 
! Use time limit to see if particle has moved long enough
               if ((ttp1(np1) - ttpo(np1)) .lt. tplim) return
            end if  
         end if   

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

c     ps_print = ps_trac(current_node)
c     s_print = s(current_node)
         end if
         position_in_string = 1
         
         do i = 1, iprcount
            if (i .eq. write_prop(7) .or. i .eq. write_prop(8)) 
     &           then
               write (sptr_prop_values(position_in_string:
     &              position_in_string+9), '(i8,2x)') int(sptr_prop(i))
               position_in_string = position_in_string+10
            else
               write (sptr_prop_values(position_in_string:
     &              position_in_string+17), '(g16.9,2x)') 
     &              sptr_prop(i)
               position_in_string = position_in_string+18
            end if
         end do
         
         if (istop(np1) .eq. 2 .and. ijkv(np1) > 0) then
! Left system at spring node
            write (sptr_prop_values(position_in_string:
     &           position_in_string+18), '(i8,2x,sp,i8)')
     &           current_node, ijkv(np1)
         else if (istop(np1) .eq. 10) then
! Stopped in system at current_node 
            write (sptr_prop_values(position_in_string:
     &           position_in_string+18), '(i8,2x,sp,i9)')
     &           current_node, ijkv(np1)*(-10)
         else
! Left at flowing boundary or capture node
            write (sptr_prop_values(position_in_string:
     &           position_in_string+17), '(i8,2x,i8)')
     &           current_node, ijkv(np1)
         end if
         position_in_string = len_trim(sptr_prop_values)
!         if (sptr_time .gt. days) sptr_time = days
         write(isptr2,8001) part_id(np1,1), xcoordw, ycoordw, zcoordw,
     2        sptr_time, sptr_prop_values(1: position_in_string)
 
! Save values just written 
         count_sptr2=count_sptr2+1
         ttpo(np1) = ttp1(np1)
         xo(np1) = xcoordw
         yo(np1) = ycoordw
         zo(np1) = zcoordw
         lastnode(np1) = current_node
         call flush (isptr2)
      end if

 8001 format(1x, i8, 3(1x,g16.9), 1x, g21.14, 200a)
 9001 format(1x,i8,1x,g21.14,1x,i8)
 9002 format(1x,i8,1x,g21.14,1x,i8,3(1x,g16.9))

      return
      end subroutine write_path_info

******************************************************************
******************************************************************

      subroutine write_snapshot
      implicit none

      integer print_node

      write(isptr1,9795)'days=', tt1,itime
 9795 format(a6,g17.10,1x,i8)
      do np1=1,num_part
         upcoming_node = abs(ijkv(np1))
         if(upcoming_node.gt.0) then
            print_node=upcoming_node
            if(istop(np1).ge.2) print_node=-upcoming_node
            xcoordw = x1(np1) + corn(upcoming_node,1)
            ycoordw = y1(np1) + corn(upcoming_node,2)
            if(icnl.eq.0) then
               zcoordw = z1(np1) + corn(upcoming_node,3)
            else
               zcoordw = 0.
            end if
         elseif(upcoming_node.eq.0) then
            print_node=exit_node(np1)
            xcoordw = x3(np1)
            ycoordw = y3(np1)
            zcoordw = z3(np1)
         end if
c         write(isptr1,9796) np1, xcoordw, ycoordw, zcoordw,tt1
c     3        ,print_node,itime
         write(isptr1,9797) np1, xcoordw, ycoordw, zcoordw,
     3        print_node
 9796    format(1x, i8, 3(1x,g16.9), 1x, g21.14, 2(1x, i8))
 9797    format(1x, i8, 3(1x,g16.9), 1(1x, i8))
97977    format(1x, i8, 3(1x,g20.9), 1(1x, i8))
 
      enddo

      return
      end subroutine write_snapshot

******************************************************************
******************************************************************
      

      subroutine compute_concentrations
      implicit none

      anl = anv
      anv = 0.
      do np1 = 1, num_part
         current_node = ijkvs(np1)
         if(ijkv(np1).ne.0) then
            anv(current_node) = anv(current_node) + 1.
         end if
      end do

c     concentration in particles per fluid mass

      do i = 1, n0
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            denom = rolf(i)*sx1(i)*ps_trac(i)*s(i)
         else
            denom = rolf(i)*sx1(i)*ps_trac(i)
         end if
         if(denom.lt.1.e-20) then
            anv(i) = 0.
         else
            anv(i) = anv(i)/denom
         end if
         if(iprtr.lt.0) then
            an(i) = an(i) + 0.5*(tt1-tt1_old)*(anl(i)+anv(i))/num_part
         end if
      end do

      if(iprtr.ge.0) then
         an = anv
      end if
      days_save = days
      days = tt1
      npn = 0
      call plotc1(1,0)
      days = days_save

      return

      end subroutine compute_concentrations

******************************************************************
******************************************************************      

c      subroutine load_struct_flux_array
c made a seperate file for this sub. dec 13 01 s kelkar

******************************************************************      

      subroutine load_unstruct_flux_array
      implicit none

c     Not yet implemented

      return 
      end subroutine load_unstruct_flux_array

******************************************************************
******************************************************************      

      subroutine write_btc_values
      implicit none
      
      logical any_zone_ok, write_out, write_last
      integer delta_total, previous 
      real*8 temp_count

!      totalpart = 0
!cli updated the following two blocks so that when
!    there are no new particles that go out of the system, the
!    code will not output to sptr3. 08/21/05
      write_out=.FALSE.
      write_last=.FALSE.

cHari 01-Nov-2006 Colloid Diversity Model
c totalpart is now a real*8 to incorporate importance sampling
c which requires weighting the btc (ret_weight is normalized to num_part)
c zvd 24-Jul-07, changed totalpart back to an integer and added
c totalpart_ret to contain the weighted value (time scaling options
c were not working with real valued particle counts)
      do izone = 1, nzbtc
         temp_count=totalpart(izone)
         do np1 = 1, num_part
            if(izonebtc(izone,np1).eq.-1) then
               if(tt1.ge.ttbtc(izone,np1)) then
                  izonebtc(izone,np1) = 1
                  totalpart(izone) = totalpart(izone) + 1
                  if(div_flag) then
                     totalpart_ret(izone) = totalpart_ret(izone) +
     &                    ret_weight(np1,1)
                  endif
               end if
!            elseif(izonebtc(izone,np1).eq.1) then
!               totalpart(izone) = totalpart(izone) + 1
            end if
         end do
         if(temp_count .ne. totalpart(izone))then
            write_out=.TRUE.
            if (temp_count .eq. 0 .and. tt1 .gt. rsttime) 
     &           write_last = .TRUE.
         end if
      end do

      any_zone_ok = .FALSE.
      btc_done = .TRUE.
      totflag = .FALSE.
      do izone = 1, nzbtc
         if ( totalpart(izone) .gt. part_frac ) then
            delta_total = totalpart(izone) - totlast(izone, 1)
            if (delta_total .ge. delta_part) then
! reset counter
               totlast(izone, 2) = 1
               any_zone_ok = .TRUE.
            else
               if (totlast(izone, 2) .gt. part_steps) then
                  totlast(izone, 2) = 1
                  totflag = .TRUE.
               else
                  totlast(izone, 2) = totlast(izone, 2) + 1
               end if
            end if
! Check to see if all of the particles have left the btc zones
            if  ( totalpart(izone) .lt. num_part ) then
               btc_done = .FALSE.
            else if (totalpart(izone) .eq. num_part .and.
     &              totlast(izone, 2) .eq. 1) then
               btc_done = .FALSE.
            end if 
         else
! If enough particles haven't reached a zone yet we don't want to 
! increase time step (this could be a problem if you have a zone 
! particles never get to)
            any_zone_ok = .TRUE.
            btc_done = .FALSE.
! Check to see if time is greater than time_btc so we can reset edt to
! dtmn if it is large
            if (time_btc .ge. tt1) edt = dtmn
         end if
      end do

      if (any_zone_ok) totflag = .FALSE.
      if (btc_done) edt = dtmx

      if(write_out)then
! Reset btc_flag (only reset flag if we've written information)
         btc_flag = .false.
         if (div_flag) then
            if (write_last) write(isptr3,1500) max(0.d0, rsttime, 
     2           tt1-edt/86400.), (totlast_ret(izone), izone=1,nzbtc)
            write(isptr3,1500) tt1, (totalpart_ret(izone),
     2           izone=1,nzbtc)
         else
            if (write_last) write(isptr3,1510) max(0.d0, rsttime, 
     2           tt1-edt/86400.), (totlast(izone, 1), izone=1,nzbtc)
            write(isptr3,1510) tt1, (totalpart(izone),
     2           izone=1,nzbtc)
         end if
         if(iptty.ne.0) then
            if (div_flag) then
               if (write_last) write(iptty,1500) max(0.d0, rsttime, 
     2              tt1-edt/86400.), (totlast_ret(izone), izone=1,nzbtc)
               write(iptty,1500) tt1, (totalpart_ret(izone),
     2              izone=1,nzbtc)
            else
               if (write_last) write(iptty,1510) max(0.d0, rsttime, 
     2              tt1-edt/86400.), (totlast(izone, 1), izone=1,nzbtc)
               write(iptty,1510) tt1, (totalpart(izone),
     2              izone=1,nzbtc)
            end if
         end if
         do izone=1,nzbtc
            totlast(izone, 1) = totalpart(izone)
            if(div_flag) totlast_ret(izone) = totalpart_ret(izone)
         end do
         call flush(isptr3)
      end if
 1500 format(1x,g16.7,10(1x,g16.7))
 1510 format(1x,g16.7,10(1x,i7))

      return
      end subroutine write_btc_values

******************************************************************
******************************************************************      

      subroutine compute_part_delay(iflag1)
      implicit none
      integer iflag1
      integer numsptrspecies
c     ijkvs is used because we are delaying the particle in the
c     node the particle is in now, not the one it is about
c     to travel through 

              current_node = ijkvs(np1)
              upcoming_node = ijkv(np1)
c     Only do the time delay when a particle has just
c     reached an interface and is about to leave
              if(current_node.ne.upcoming_node) then
                 current_model = itrc(current_node)
                 ret_factor = rd_frac(current_model,1)
                 if(tprpflag(current_model).eq.1.or.
     2                tprpflag(current_model).eq.2) then
                    omegav = omega_partial(current_node)
                    sigmav = 0.

                 elseif(tprpflag(current_model).ge.11)then
cHari 01-Nov-2006 Colloid Diversity Model 
                    omegav = rcoll_div(np1,1)
                    sigmav = 0.             

                 
                 elseif(aperture(current_model).lt.0.) then
                    omegav = omega_partial(current_node)
                    sigmav = sigma_partial(current_node)
                 else
c     Total time since particle entered the cell is used for
c     computing the delay time parameters
                    sigmav = sigma_partial(current_node)
     2                   /sqrt(tt(np1))
                    omegav = omega_partial(current_node)
     2                   *sqrt(tt(np1))
                 end if
                 par3v = 1.
                 fm = 1
                 numsptrspecies = 1
                 call time_delay(tprpflag(current_model),
     2                numsptrspecies,current_node,sigmav,
     3                omegav,par3v,fm,ret_factor,rseed,tt(np1),
     4                conc_ret,dt_delay )
c     To get the correct particle time after applying the
c     delay, we need to substract off the time since the cell
c     was entered tsubtract from ttp1, and then we add
c     the total time with delay 
                 ttp1(np1) = ttp1(np1) + 
     2                (dt_delay-tsubtract)/86400.
c     Reset the time clock for a particle within the current cell
                 tt(np1) = 0.
              elseif(iflag1.eq.0) then
c     No delay yet because we are inside the cell. Just add 
c     the time step w/o delay to the running toal to keep it
c     current for output purposes for cases w/o delay
c     Only do this for the case where we are here not from bouncing
c     in from another cell
                 ttp1(np1) = ttp1(np1) + dt(np1)/86400.
              end if

              return
              end subroutine compute_part_delay


******************************************************************
******************************************************************      
              
      subroutine new_part_loc_struct

      use comsk

      implicit none

      integer inp1,icapture,ijkv_in
      real*8 dx,dy,dz,x,xx,y,yy,z,zz
      
      x2save = x2(np1)
      y2save = y2(np1)
      z2save = z2(np1)
      ijkvsave = ijkv(np1)
      ijkv_in=ijkv(np1)
      ijkvssave = ijkvs(np1)
      ddxvsave = ddxv(np1)
      ddyvsave = ddyv(np1)
      ddzvsave = ddzv(np1)
      
      done = .false.
c     Cell-based velocities
c     v_x = 0.5*(ggg(ijkvsave,+1)-ggg(ijkvsave,-1))
c     v_y = 0.5*(ggg(ijkvsave,+2)-ggg(ijkvsave,-2))
c     v_z = 0.5*(ggg(ijkvsave,+3)-ggg(ijkvsave,-3))
c     vcurrent = sqrt(v_x**2+v_y**2+v_z**2)
c     Local velocities
      v_x = (ggg(ijkvsave,+1)+ggg(ijkvsave,-1))/ddxv(np1)
      v_x=v_x*x2(np1)-ggg(ijkvsave,-1)
      v_y = (ggg(ijkvsave,+2)+ggg(ijkvsave,-2))/ddyv(np1)
      v_y=v_y*y2(np1)-ggg(ijkvsave,-2)
      v_z = (ggg(ijkvsave,+3)+ggg(ijkvsave,-3))/ddzv(np1)
      v_z=v_z*z2(np1)-ggg(ijkvsave,-3)
      vcurrent = sqrt(v_x**2+v_y**2+v_z**2)
      walk_counter = 0
      
      do while(.not.done)
         
         walk_counter = walk_counter + 1
         
         x2(np1) = x2save
         y2(np1) = y2save
         z2(np1) = z2save
         ijkv(np1) = ijkvsave
         ijkvs(np1) = ijkvssave
         ddxv(np1) = ddxvsave
         ddyv(np1) = ddyvsave
         ddzv(np1) = ddzvsave
         current_node = ijkvs(np1)
         
c     Perform random walk for this particle,
   
!         if(upcoming_node.ne.0.and.
!     2        irray(upcoming_node,0).gt.0) then
         if(upcoming_node.ne.0) then
            current_model = itrc(current_node)
c...s kelkar  april 7 04............................................
c the following are for dispersion/random-walk and not fixed for OMR
c hence skipped if omr_flag is "TRUE"
            if (.not. omr_flag) then
               call random_walk(1,np1,current_node,
     2              current_model,edt,x_ts,y_ts,z_ts,itime,istep
     3              ,time_step_flag)
            
cc for debugging 8/15/02.......................................
               call random_walk(0,np1,current_node,
     2              current_model,edt,x_ts,y_ts,z_ts,itime,istep
     3              ,time_step_flag)
            end if
cc............................................................
            
c...  check if near a capture node, 1/22/01 s kelkar 
            
            inp1=ijkv(np1)
            if(irray(inp1,0).eq.-inp1.or.irray(inp1,0).
     1           lt.-10000000) then
               dx=ddx(inp1)
               dy=ddy(inp1)
               
               x=x1(np1)-0.5*dx
               y=y1(np1)-0.5*dy
               xx=x2(np1)-0.5*dx
               yy=y2(np1)-0.5*dy
               
               call check_well_capture(inp1,x,y,xx,yy,icapture,
     1              vx_zheng,vy_zheng,vz_zheng,rw)               
               if(icapture.eq.1) then
c     particle has entered the extraction well, remove it from the system
                  exit_node(np1)=-ijkv(np1)
                  call capture_particle(np1)
                  goto 9999
                  
               else
c     particle has not entered the well, continue
                  
               endif
               
            endif

c.................
            
            
            inside = .false.
            ntry_pos = 0
            do while(.not.inside.and.ntry_pos.le.10)
               ntry_pos = ntry_pos + 1
c     Move particles to new cell in the x direction
               
               local_pos = x2(np1)
               local_spacing = ddxv(np1)
               call move_particle(np1, ddx, local_pos,
     2              local_spacing, too_many_jumps,1)
               if(too_many_jumps) goto 2000
               x2(np1) = local_pos
               ddxv(np1) = local_spacing
               
               
c     Move particles to new cell in the y direction
               
               local_pos = y2(np1)
               local_spacing = ddyv(np1)
               call move_particle(np1, ddy, local_pos,
     2              local_spacing, too_many_jumps,2)
               if(too_many_jumps) goto 2000
               y2(np1) = local_pos
               ddyv(np1) = local_spacing
               
               
c     Move particles to new cell in the z direction
               
               if(icnl.eq.0) then
                  
                  local_pos = z2(np1)
                  local_spacing = ddzv(np1)
                  call move_particle(np1, ddz, local_pos,
     2                 local_spacing, too_many_jumps,3)
                  if(too_many_jumps) goto 2000
                  z2(np1) = local_pos
                  ddzv(np1) = local_spacing
               end if
               
               x2nondim = x2(np1)/ddxv(np1)
               y2nondim = y2(np1)/ddyv(np1)
               z2nondim = z2(np1)/ddzv(np1)
               
               if(x2nondim.le.1..and.x2nondim.ge.0.) then
                  if(y2nondim.le.1..and.y2nondim.ge.0.) then
                     if(z2nondim.le.1..and.z2nondim.ge.0.) then
                        inside = .true.
                     end if
                  end if
               end if

            end do

c            if(.not.inside) then
c               if (iptty .ne. 0) write(iptty,*)
c     1              'Warning: Particle # ',np1,
c     2              'location determination required placing particle
c     3              in the middle of the cell #', ijkv(np1)
c            end if

c...  check if near a capture node, 1/22/01 s kelkar 
            
!            inp1=ijkv(np1)
!            if(irray(inp1,0).eq.-inp1.or.irray(inp1,0).
!     1           lt.-10000000) then
!               icapture=1
c particle has entered the extraction well, remove it from the system
!               exit_node(np1)=-ijkv(np1)
!               call capture_particle(np1)
!               goto 9999
!            endif

            
            
c     Now determine if the particle has artifically jumped into 
c     slow moving fluid
c     First, skip this section if the move_particle routine
c     returned too many jumps for the particle, because we want to
c     try over
 2000       continue
            if(.not.too_many_jumps) then
c     if(ijkv(np1).ne.ijkvsave) then
c     If new velocity is loo low, we're not done, go back and try again
c     v_x = 0.5*(ggg(ijkv(np1),+1)-ggg(ijkv(np1),-1))
c     v_y = 0.5*(ggg(ijkv(np1),+2)-ggg(ijkv(np1),-2))
c     v_z = 0.5*(ggg(ijkv(np1),+3)-ggg(ijkv(np1),-3))
               v_x = (ggg(ijkv(np1),+1)+ggg(ijkv(np1),-1))/
     2              ddxv(np1)
               v_x=v_x*x2(np1)-ggg(ijkv(np1),-1)
               v_y = (ggg(ijkv(np1),+2)+ggg(ijkv(np1),-2))/
     2              ddyv(np1)
               v_y=v_y*y2(np1)-ggg(ijkv(np1),-2)
               v_z = (ggg(ijkv(np1),+3)+ggg(ijkv(np1),-3))/
     2              ddzv(np1)
               v_z=v_z*z2(np1)-ggg(ijkv(np1),-3)
               vnew = sqrt(v_x**2+v_y**2+v_z**2)
               if(vnew/vcurrent.ge.
     2              abs(vratio(current_model))) then
c     Velocity in new cell is ok, we're done
                  done = .true.
               end if
            else
               done = .true.
            end if
         end if
         
         if(.not.done.and.walk_counter.ge.10) then
c     Can't get a valid jump, leave particle where it is and proceed
            x2(np1) = x2save
            y2(np1) = y2save
            z2(np1) = z2save
            ijkv(np1) = ijkvsave
            ijkvs(np1) = ijkvssave
            done = .true.
         end if
         
      end do
      
 9999 continue

c if  particle has made a random jump on the
c  previous call to new_part_loc_struct. In that case 
c  reset oldnode and oldnode2 to ijkv(np1) to avoid confusing
c this with a bouncing particle
      If(ijkv_in.ne.ijkv(np1)) then
         oldnode(np1)=ijkv(np1)
         oldnode2(np1)=ijkv(np1)
      endif

      return

      end subroutine new_part_loc_struct

******************************************************************
******************************************************************      

      subroutine new_part_loc_unstruct
      implicit none
      
c     Not yet implemented
      
      return
      end subroutine new_part_loc_unstruct
      
******************************************************************
******************************************************************      

      subroutine set_velocities_struct
      implicit none
      
      ijkvnp=ijkv(np1)
      vx2bv(np1)= ggg(ijkvnp, 1)
      vx1bv(np1)=-ggg(ijkvnp,-1)
      vy2bv(np1)= ggg(ijkvnp, 2)
      vy1bv(np1)=-ggg(ijkvnp,-2)
      vz2bv(np1)= ggg(ijkvnp, 3)
      vz1bv(np1)=-ggg(ijkvnp,-3)
      ddxv(np1)=ddx(ijkvnp)
      ddyv(np1)=ddy(ijkvnp)
      ddzv(np1)=ddz(ijkvnp)
c If dx, dy, or dz is zero we don't have enough neighbors
c so set stop for this particle
      if (ddxv(np1) .eq. 0.d0) then
         axv(np1)=0.d0
c         istop(np1) = 10
      else
         axv(np1)=(vx2bv(np1)-vx1bv(np1))/ddxv(np1)
      end if
      if (ddxv(np1) .eq. 0.d0) then
         ayv(np1)=0.d0
c         istop(np1) = 10         
      else
         ayv(np1)=(vy2bv(np1)-vy1bv(np1))/ddyv(np1)
      end if
      if (ddxv(np1) .eq. 0.d0) then
         ayv(np1)=0.d0
c         istop(np1) = 10
      else
         azv(np1)=(vz2bv(np1)-vz1bv(np1))/ddzv(np1)
      end if
      axyzm(np1)= abs(axv(np1))+1.e-28
      if(ayv(np1).gt.axyzm(np1)) axyzm(np1)=abs(ayv(np1))
      if(azv(np1).gt.axyzm(np1)) axyzm(np1)=abs(azv(np1))

      return
      end subroutine set_velocities_struct

******************************************************************
******************************************************************      

      subroutine set_velocities_unstruct
      implicit none


c     Not yet implemented
      
      
      return
      end subroutine set_velocities_unstruct

******************************************************************
******************************************************************      

      subroutine set_ts_struct
      implicit none


c     Subroutine to determine maximum time step to honor Courant
c     condition for random walk dispersion
      
      real*8 ddtx,ddty,ddtz
      
      current_node = ijkvnp
c.......................................
      
      if(current_node.ne.0) then
         current_model = itrc(current_node)
         if(tprpflag(current_model).eq.2.or.
     2        tprpflag(current_model).eq.4.or.
     2        tprpflag(current_model).eq.13.or.
     2        tprpflag(current_model).eq.14) then
            
            if(courant_factor.gt.0.) then
c     Velocity
               edtx = abs(0.5*(ggg(current_node,1)-
     2              ggg(current_node,-1)))
               edty = abs(0.5*(ggg(current_node,2)-
     2              ggg(current_node,-2)))
               edtz = abs(0.5*(ggg(current_node,3)-
     2              ggg(current_node,-3)))
c     Time step based on velocity only
               edtx = ddxv(np1)/max(1.d-30,edtx)
               edty = ddyv(np1)/max(1.d-30,edty)
               edtz = ddzv(np1)/max(1.d-30,edtz)
c find time-step size based on del.(D) terms
c               ddtx=1.e+30
c               ddty=1.e+30
c               ddtz=1.e+30
c               if (abs(divd_omr(current_node,1)).gt.0.)
c     $              ddtx=ddxv(np1)/abs(divd_omr(current_node,1))
c               if (abs(divd_omr(current_node,2)).gt.0.)
c     $              ddty=ddyv(np1)/abs(divd_omr(current_node,2))
c               if (abs(divd_omr(current_node,1)).gt.0.)
c     $              ddtz=ddzv(np1)/abs(divd_omr(current_node,3))
c
c     Find minimum ts, apply Courant condition
c               edtmax(np1) = courant_factor*min(edtx,edty,edtz,
c     $              ddtx,ddty,ddtz)
c zvd 13-May-09 If there is a zero velocity in x, y, or z ignore it when
c calculating time step
               if (edtx .ne. 0 .and. edty .ne. 0 .and. edtz .ne.0) then
                  edtmax(np1) = courant_factor*min(edtx,edty,edtz)
               else if (edtx .ne. 0 .and. edty .ne. 0) then
                  edtmax(np1) = courant_factor*min(edtx,edty)
               else if (edtx .ne. 0 .and. edtz .ne. 0) then
                  edtmax(np1) = courant_factor*min(edtx,edtz)
               else if (edty .ne. 0 .and. edtz .ne. 0) then
                  edtmax(np1) = courant_factor*min(edty,edtz)
               else
                  edtmax(np1) = courant_factor*max(edtx,edty,edtz)
               end if
               if (edtmax(np1) .eq. 0.d0) then
c There is a problem with this particle so stop it
                  write (ierr, *) 'Particle ', np1, ' time problem'
                  istop(np1) = 10
               end if

            else
c...s kelkar  april 7 04............................................
c the following are for dispersion/random-walk and not fixed for OMR
c hence skipped if omr_flag is "TRUE"
c
cc     time step based on estimate of magnitude of dispersion
cc     in x, y, and y directions

               call random_walk(0,np1,current_node,current_model,
     2              edt, x_ts, y_ts, z_ts,itime,istep,time_step_flag)
c................................................................

               edtmax(np1) = abs(courant_factor)*
     2              min(x_ts,y_ts,z_ts)
            end if
            
         else
            edtmax(np1) = 1.e30
         end if
      end if
      
      return
      end subroutine set_ts_struct

******************************************************************
******************************************************************      


      subroutine set_ts_unstruct
      implicit none

c     Not yet implemented
      
      return
      end subroutine set_ts_unstruct

******************************************************************
******************************************************************      

      subroutine transtime_struct
      implicit none

      real*8 epsilon, temp1,temp2,temp3,epomr

      epsilon=1.e-10
      epomr=1.e-4

c.....1/22/01 s kelkar particle capture
c irray(i,0)=-i if the particle in a capture cell, use Zheng formulation

      if(irray(ijkv(np1),0).eq.-ijkv(np1).or.irray(ijkv(np1),0).
     1           lt.-10000000) then
         call time_zheng(icnl,np1,vx_zheng,vy_zheng,vz_zheng,rw)
c......................................

      else         
         
c     Handle velocities in the x-direction
         
         vx1(np1)=axv(np1)*(x1(np1))+vx1bv(np1)
         
         
         if(vx1(np1).eq.0.) then
            dtx(np1)=1.e+30
         elseif(vx1(np1).gt.0.) then
            x_new(np1)=vx2bv(np1)/vx1(np1)
c s kelkar  april 13 04 3DOMR modified the next conditioner
c to exclude the case when particle starts on a face. this is a 
c quick fix. 
            if(abs(axv(np1)).le.small_slope.or.
     $           abs(x_new(np1)-1.).le.epomr) then
               dtx(np1)=(ddxv(np1)-x1(np1))/vx1(np1)
            else
               if(x_new(np1).lt.ep) then
                  dtx(np1)=1.e+30
               else
                  dtx(np1)=log(x_new(np1))/axv(np1)
               end if
            end if
         else
            x_new(np1)=vx1bv(np1)/vx1(np1)
            if(abs(axv(np1)).le.small_slope.or.
     $           abs(x_new(np1)-1.).le.epomr) then
               dtx(np1)=-x1(np1)/vx1(np1)
            else
               if(x_new(np1).lt.ep) then
                  dtx(np1)=1.e+30
               else
                  dtx(np1)=log(x_new(np1))/axv(np1)
               end if
            end if
         end if
         
c     Handle velocities in the y-direction
         
         vy1(np1)=ayv(np1)*(y1(np1))+vy1bv(np1)
         
         
         if(vy1(np1).eq.0.) then
            dty(np1)=1.e+30
         elseif(vy1(np1).gt.0.) then
            y_new(np1)=vy2bv(np1)/vy1(np1)
            if(abs(ayv(np1)).le.small_slope.or.
     $           abs(y_new(np1)-1.).le.epomr) then
               dty(np1)=(ddyv(np1)-y1(np1))/vy1(np1)
            else
               if(y_new(np1).lt.ep) then
                  dty(np1)=1.e+30
               else
                  dty(np1)=log(y_new(np1))/ayv(np1)
               end if
            end if
         else
            y_new(np1)=vy1bv(np1)/vy1(np1)
            if(abs(ayv(np1)).le.small_slope.or.
     $           abs(y_new(np1)-1.).le.epomr) then
               dty(np1)=-y1(np1)/vy1(np1)
            else
               if(y_new(np1).lt.ep) then
                  dty(np1)=1.e+30
               else
                  dty(np1)=log(y_new(np1))/ayv(np1)
               end if
            end if
         end if
         
c     Handle velocities in the z-direction
         
         vz1(np1)=azv(np1)*(z1(np1))+vz1bv(np1)
         
         
         if(vz1(np1).eq.0.) then
            dtz(np1)=1.e+30
         elseif(vz1(np1).gt.0.) then
            z_new(np1)=vz2bv(np1)/vz1(np1)
            if(abs(azv(np1)).le.small_slope.or.
     $           abs(z_new(np1)-1.).le.epomr) then
               dtz(np1)=(ddzv(np1)-z1(np1))/vz1(np1)
            else
               if(z_new(np1).lt.ep) then
                  dtz(np1)=1.e+30
               else
                  dtz(np1)=log(z_new(np1))/azv(np1)
               end if
            end if
         else
            z_new(np1)=vz1bv(np1)/vz1(np1)
            if(abs(azv(np1)).le.small_slope.or.
     $           abs(z_new(np1)-1.).le.epomr) then
               dtz(np1)=-z1(np1)/vz1(np1)
            else
               if(z_new(np1).lt.ep) then
                  dtz(np1)=1.e+30
               else
                  dtz(np1)=log(z_new(np1))/azv(np1)
               end if
            end if
         end if

c....s kelkar jan 29, 2001 ..............         
      endif
c..............
         
c     Set dt
c s kelkar  april 12 04 .. running into problems with dt being
c very small when the point lands on or near one of the three faces
c accounted in the above code. Hence find the min this way
      temp1=min(ddx(ijkv(np1)),ddy(ijkv(np1)),ddz(ijkv(np1)))
      temp2=max(abs(vx1(np1)),abs(vy1(np1)),abs(vz1(np1)))
      temp1=temp1/temp2
      temp1=temp1*epomr
      temp3=1.e+30
      if(dtx(np1).gt.temp1.and.dtx(np1).lt.temp3) temp3=dtx(np1)
      if(dty(np1).gt.temp1.and.dty(np1).lt.temp3) temp3=dty(np1)
      if(dtz(np1).gt.temp1.and.dtz(np1).lt.temp3) temp3=dtz(np1)

c.............................................
c      dt(np1) = min(x61(np1),dtx(np1),dty(np1),dtz(np1),
c     2     edtmax(np1))
      dt(np1) = min(x61(np1),temp3,
     2     edtmax(np1))
      
c     Reset dt to hit the end of the time step exactly if it 
c     will be in the current cell at that time. Set flag ioutt
c     so record this
      
      if(ioutt(np1).eq.0) then
         if(dt(np1)/86400..gt.tt1-ttp1(np1)) then
            dt(np1)=86400.*(tt1-ttp1(np1))
            ioutt(np1) = 1
         end if
      end if
      
c     Add check to make sure dt is positive
      
      dt(np1) = max(1.d-10,dt(np1))
      
      return

      end subroutine transtime_struct
         
******************************************************************
******************************************************************      

      subroutine transtime_unstruct
      implicit none

c     Not yet implemented

      return
      end subroutine transtime_unstruct

******************************************************************
******************************************************************      


******************************************************************
******************************************************************      

      subroutine compute_exit_unstruct
      implicit none

c     Not yet implemented

      return
      end subroutine compute_exit_unstruct


******************************************************************
******************************************************************      

      subroutine check_part_loc_struct
      implicit none

      integer flagout, ijkvnp1

c     As a final check, make sure that no particle is outside the cell
c     If it is, change the coordinate in error and place the particle in
c     the middle of the current cell

      flagout=0
      ijkvnp1=ijkv(np1)

      if(ijkvnp1.ne.0) then
         
         x2nondim = x2(np1)/ddx(ijkvnp1)
         if(x2nondim.gt.1.or.x2nondim.lt.0) then
            flagout=1
            x2(np1) = ddx(ijkvnp1)/2.
         end if
         
         y2nondim = y2(np1)/ddy(ijkvnp1)
         if(y2nondim.gt.1.or.y2nondim.lt.0) then
            flagout=1
            y2(np1) = ddy(ijkvnp1)/2.
         end if
         
         z2nondim = z2(np1)/ddz(ijkvnp1)
         if(z2nondim.gt.1.or.z2nondim.lt.0) then
            flagout=1
            z2(np1) = ddz(ijkvnp1)/2.
         end if
         
      end if
      
!      if(flagout.eq.1) then
!         write(ierr,*)'Warning: Particle # ',np1,
!     2        'location determination required placing particle ',
!     3        'in the middle of the cell #', ijkv(np1)
!         if (iptty .ne. 0) then
!            write(iptty,*)'Warning: Particle # ',np1,
!     2           'location determination required placing particle ',
!     3           'in the middle of the cell #', ijkv(np1)
!         endif
!      endif

      if(ijkvnp1.ne.0) then
         x1(np1) = x2(np1)
         y1(np1) = y2(np1)
         z1(np1) = z2(np1)
      else
         istop(np1) = 1
      end if
c      if(irray(ijkv(np1),0).lt.0) then
c         istop(np1) = 1
c         ijkv(np1) = 0
c      end if


      return
      end subroutine check_part_loc_struct


******************************************************************
******************************************************************      

      subroutine check_part_loc_unstruct
      implicit none

c     not yet implemented

      return
      end subroutine check_part_loc_unstruct


******************************************************************

**********************************************************


**********************************************************

******************************************************************
c......... s kelkar  jul23 01

      subroutine seek_plumzone(iparticle)
      implicit none
      integer iparticle


      upcoming_node = ijkv(iparticle)
      if(upcoming_node.ne.0) then
         current_zone = izonef(upcoming_node)
         zone_loop: do izone = 1, nplum
         if(plum(izone).eq.current_zone) then
            if(izoneplum(izone,iparticle).eq.0) then
c     Particle just got here, set flag to -1
c     Otherwise, the particle has already gotten here and
c     we want to leave the flag alone (-1 if particle time
c     has not yet been reached, or 1 if it already has)
               izoneplum(izone,iparticle) = -1
            end if
            exit zone_loop
         end if
      end do zone_loop
      end if
      return
      end subroutine seek_plumzone

******************************************************************      

      subroutine write_plum_values
      implicit none

      integer iskptr

      iskptr=isptr4

      call write_path_info

      plumpart = 0
      do izone = 1, nplum
         do np1 = 1, num_part
            if(izoneplum(izone,np1).eq.-1) then
               if(tt1.ge.ttp1(np1)) then
                  izoneplum(izone,np1) = 1
                  plumpart(izone) = plumpart(izone) + 1
               end if
            elseif(izoneplum(izone,np1).eq.1) then
               plumpart(izone) = plumpart(izone) + 1
            end if
         end do
      end do

c     Time is the time before taking the time step, since
c     the particle movement for this time step have not yet
c     been done
      
      write(isptr5,1500) tt1-edt/86400., (totalpart(izone),
     2     izone=1,nzbtc)
      if(iptty.ne.0) then
         write(iptty,1500) tt1-edt/86400., (totalpart(izone),
     2        izone=1,nzbtc)
      end if
 1500 format(1x,g16.7,10(1x,i7))
      return

      end subroutine write_plum_values

*********************************************************

      end subroutine ptrac3


******************************************************

      subroutine  time_zheng(icnl,np1,vx_zheng,vy_zheng,vz_zheng,rw)

      use comdi
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,icnl

      real*8 tfac,epsilon
      real*8 dx,dy,dz,x,y,z,xx,yy,zz,skz,rr
      real*8 qwp,qqw,a,fac,xyw,pi

      real*8 vx_zheng,vy_zheng,vz_zheng,rw

      pi=3.14159
      epsilon=1.e-20
      tfac = 0.1

c use equations from Zheng, 1994
c*****check with bruce about the correct use of ggg sk etc
      inp1=ijkv(np1)
      dx=ddx(inp1)
      dy=ddy(inp1)
      dz=ddz(inp1)
      x=x2(np1)-0.5*dx
      y=y2(np1)-0.5*dy
      z=z2(np1)-0.5*dz

c distance to the particle from the welbore, for timestep size
      rr=sqrt(x*x+y*y)

c convert sk from kg/sec to m^3/sec assuming rhow=1000kg/m^3
c To be consistant with Zheng's conventions, divide sk by dx*dy*dz and 
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
      else
         vx_zheng=fac*x/xyw+0.5*(+ggg(inp1,1)-ggg(inp1,-1))
      endif

      if(abs(vx_zheng).gt.epsilon) then
c using the distance to the wellbore, x, as the  sale
         dtx(np1)=abs(tfac*rr/vx_zheng)
      else
         dtx(np1)=1.e+30
      endif

      if(irevers.eq.1) then
         vy_zheng=fac*y/xyw+0.5*(-ggg(inp1,2)+ggg(inp1,-2))
      else
         vy_zheng=fac*y/xyw+0.5*(+ggg(inp1,2)-ggg(inp1,-2))
      endif

      if(abs(vy_zheng).gt.epsilon) then
c using the distance to the wellbore, x, as the  sale
         dty(np1)=abs(tfac*rr/vy_zheng)
      else
         dty(np1)=1.e+30
      endif
      
      if(icnl.eq.0) then

         if(irevers.eq.1) then
       vz_zheng=((-ggg(inp1,3)-ggg(inp1,-3))/dz)*(z-0.5*dz)+ggg(inp1,-3)
         else
       vz_zheng=((+ggg(inp1,3)+ggg(inp1,-3))/dz)*(z-0.5*dz)-ggg(inp1,-3)
         endif

         if(abs(vz_zheng).gt.epsilon) then
            dtz(np1)=abs(tfac*dz/vz_zheng)
         else
            dtz(np1)=1.e+30
         endif
      endif
      
      return
      
      end  subroutine  time_zheng

c.........................................................

      subroutine check_well_capture(inp1,x,y,xx,yy,icapture,
     1     vx_zheng,vy_zheng,vz_zheng,rw)

c     check if particle has entered the well. this is done
c     by solving the equation of a circle with radius rw
c     and the st line given by (x,y) and (xx,yy)

      use comai , only : neq
      use comdi
      use comsptr
      use comsk

      implicit none

      logical exit_flag

      integer icapture,inp1,indexrw

      real*8 x,y,z,xx,yy,zz,dis2
      real*8 ss,c,fac,rw,xx1,xx2,yy1,yy2

      real*8 vx_zheng,vy_zheng,vz_zheng

c....  rw as an input
      indexrw=-(irray(inp1,0)+10000000)
      if(indexrw.gt.0.and.indexrw.le.neq) then
         rw=well_radius(indexrw)
      else
         rw = 0.1
      endif
c first check if xx,yy are already inside the circle
      dis2=xx*xx+yy*yy
      if(dis2.le.rw*rw) then
         icapture=1
      else

         if(x.ne.xx) then
c     particle track not parallel to y axis, finite slope
            ss=(yy-y)/(xx-x)
            c=y-ss*x
            fac=rw*rw*(1+ss*ss)-c*c
            if(fac.gt.0.) then
c     equation has real roots, see if either is between x and xx
               fac=sqrt(fac)
               xx1=(-ss*c+fac)/(1.+ss*ss)
               xx2=(-ss*c-fac)/(1.+ss*ss)
               if((xx-xx1)*(x-xx1).lt.0.) then
c     particle has entered the extraction well, set flag
                  icapture=1
               elseif((xx-xx2)*(x-xx2).lt.0.) then
c     particle has entered the extraction well, set flag
                  icapture=1
               else
c     particle has not entered the well
                  icapture=0
               endif
            else
c     particle track does not intersect the well, set flag
               icapture=0
            endif
            
         else
            
c     particle track parallel to y axis
            fac=rw*rw-x*x
            if(fac.gt.0.) then
c     equation has real roots, see if either is between y and yy
               yy1=sqrt(fac)
               yy2=-yy1
               if((yy-yy1)*(y-yy1).lt.0.) then
c     particle has entered the extraction well, set flag
                  icapture=1
               elseif((yy-yy2)*(y-yy2).lt.0.) then
c     particle has entered the extraction well, set flag
                  icapture=1
               else
c     particle has not entered the well
                  icapture=0
               endif
            else
c     particle track does not intersect the well
               icapture=0
            endif
            
         endif
         
      endif

      return
      
      end subroutine check_well_capture
               
c..........................................................

      subroutine capture_particle(np1)

      use comsptr
      use comsk

      implicit none

      integer np1

      istop(np1)=2

      return

      end subroutine capture_particle

c.......................................................

      subroutine write_particle_exit(np1,isptr4,itime)
c..s kelkar march 1 02
c write the time and location for each exiting particle

      use comsptr

      implicit none

      integer print_node,np1,isptr4,itime

      real*8 xcoordw, ycoordw, zcoordw

      logical:: opnd = .FALSE.

      inquire (isptr4, opened=opnd)
      if (.not.opnd) return

      print_node=exit_node(np1)
      xcoordw = x3(np1)
      ycoordw = y3(np1)
      zcoordw = z3(np1)
      write(isptr4,9797)ttp1(np1), xcoordw, ycoordw, zcoordw,
     3     print_node,np1,itime
 9797 format(4(2x,g16.7), 3(2x, i8))
 
      return

      end subroutine write_particle_exit

c.....................................................................
c.............................................................

      subroutine find_quadrant(np1,i,ix,iy,iz,xp,yp,zp,
     $                       vjx,vjy,vjz)

c find which quadrant of the control volume the particle
c is in. Note that i here refers to the node # that the control
c volume corresponds to and x2,y2,z2, are the local coordinates
c ix,iy,iz are the nodes that define the other 3 corners of the cube

      use comai
      use comsptr
      use combi

      implicit none

      integer np1,i,ix,iy,iz

      real*8 xp,yp,zp,vjx,vjy,vjz

      xp= x1(np1) + corn(i,1)
      yp= y1(np1) + corn(i,2)

      if(xp.ge.cord(i,1))then
         ix=+1
         vjx=+ggg(i,+1)
      else
         ix=-1
         vjx=-ggg(i,-1)
      endif

      if(yp.ge.cord(i,2)) then
         iy=+1
         vjy=+ggg(i,+2)
      else
         iy=-1
         vjy=-ggg(i,-2)
      endif

      if(icnl.eq.0) then
         zp = z1(np1) + corn(i,3)
         if(zp.ge.cord(i,3)) then
            iz=+1
            vjz=+ggg(i,+3)
         else
            iz=-1
            vjz=-ggg(i,-3)
         endif
      else
         zp = 0.
         iz=i
      end if

      return

      end

c.......................................................

      subroutine  matinv2(aa,aainv)

      use comai , only : ierr
      implicit none

      real*8 aa(2,2),aainv(2,2),det,epsilon

      epsilon = 1.e-18

      det=aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
      
      if(abs(det).gt.epsilon) then
         aainv(1,1)=+aa(2,2)/det
         aainv(1,2)=-aa(1,2)/det
         aainv(2,1)=-aa(2,1)/det
         aainv(2,2)=+aa(1,1)/det
      else
         write(ierr,*)'ERROR in matinv2. STOP'
         write(ierr,*)aa
         call exit_ptrac_omr
      endif

      return
      
      end

c...................................................


      subroutine matinv3(a,c)
c calculate the inverse of a 3x3 matrix "a" and return
c it in "c"

      use comai , only : ierr
      implicit none

      integer i,j

      real*8 a(3,3),c(3,3)
      real*8 d(3,3),det,epsilon

      epsilon=1.e-22

      d(1,1)=+a(2,2)*a(3,3)-a(2,3)*a(3,2)
      d(1,2)=-a(2,1)*a(3,3)+a(2,3)*a(3,1)
      d(1,3)=+a(2,1)*a(3,2)-a(2,2)*A(3,1)
      d(2,1)=-a(1,2)*a(3,3)+a(1,3)*a(3,2)
      d(2,2)=+a(1,1)*a(3,3)-a(1,3)*a(3,1)
      d(2,3)=-a(1,1)*a(3,2)+a(1,2)*a(3,1)
      d(3,1)=+a(1,2)*a(2,3)-a(1,3)*a(2,2)
      d(3,2)=-a(1,1)*a(2,3)+a(1,3)*a(2,1)
      d(3,3)=+a(1,1)*a(2,2)-a(1,2)*a(2,1)

      det=a(1,1)*d(1,1)+a(1,2)*d(1,2)+a(1,3)*d(1,3)
      if(abs(det).gt.epsilon) then
         do i=1,3
            do j=1,3
               c(i,j)=d(j,i)/det
            enddo
         enddo
      else
         write(ierr,*)'error in matinv. det=0. stop'
         write(ierr,*)a
         call exit_ptrac_omr
      endif
         
      return

      end

c...................................


      subroutine seek_btczone(iparticle, upcoming_node)

      use combi, only : izonef
      use comsptr

      implicit none

      integer iparticle, upcoming_node, current_zone
      integer izone

      if(upcoming_node.ne.0) then
         current_zone = izonef(upcoming_node)
         zone_loop: do izone = 1, nzbtc
            if(zbtc(izone).eq.current_zone) then
               if(izonebtc(izone,iparticle).eq.0) then
c     Particle just got here, set flag to -1
c     Otherwise, the particle has already gotten here and
c     we want to leave the flag alone (-1 if particle time
c     has not yet been reached, or 1 if it already has)
                  izonebtc(izone,iparticle) = -1
                  ttbtc(izone,iparticle) = ttp1(iparticle) + 
     &                 dt(iparticle)/86400.
                  if (alt_btc) then
                     call write_alt_btc(izone, iparticle, upcoming_node)
                  else
                     btc_flag = .true.
                  end if
               end if
               exit zone_loop
            end if
         end do zone_loop
      end if
      return
      end subroutine seek_btczone

c...................................

      subroutine write_sptrs

      use comai, only : isptr8
      use compart, only : rseed
      use comsptr
      use comxi

      implicit none
      integer i, j
      logical particles_in_system, opnd

      inquire(unit = isptr8, opened = opnd)
      if (.not. opnd) then
         open(isptr8, file = nmfil(25), status = cstats(25))
         read(isptr8, *)
         read(isptr8, *)
      end if

      write (isptr8, *) rseed
      particles_in_system = .FALSE.
      do i = 1, num_part
         if (ijkv(i) .ne. 0 .and. istop(i) .eq. 0) then
! output only if particle is still in the system
            particles_in_system = .TRUE.
            write(isptr8, 100) part_id(i,2), x3(i), y3(i), z3(i),
     &           part_id(i,1), ijkv(i), ttp1(i)
         end if
      end do

      if (particles_in_system) then
         if (nzbtc .ne. 0) then
            write (isptr8, 101) nzbtc
            write (isptr8, *) (zbtc(i), i = 1, nzbtc)
            write (isptr8, *) (totalpart(i), i = 1, nzbtc)
            do i = 1, num_part
               if (ijkv(i) .ne. 0 .and. istop(i) .eq. 0) then
                  write (isptr8, *) (izonebtc(j, i), j = 1, nzbtc)
               end if
            end do
         end if
         write(isptr8,*)
      else
         write(isptr8,*) 'All particles have exited the system'
      end if
      close (isptr8)
 100  format (i8, 3(1x, g16.9), 2(1x, i8), 1x,  g21.14)
 101  format ('zbtc', i10)
      end subroutine write_sptrs

******************************************************************
******************************************************************
     
      subroutine write_alt_btc(izone, iparticle, upcoming_node)

      use comai, only : isptr3
      use combi, only : cord, izonef
      use comdi, only : pnx, pny, pnz
      use comsptr
     
      implicit none
      integer iparticle, izone, upcoming_node
      real*8 px, py, pz

c Instead of using the coordinate of the  upcoming_node, use the position where the particle enters the zone
      if (xyz_flag2) then
         if (write_prop(3) .ne. 0) then
            px = log10(1.d-6*pnx(upcoming_node))
c            py = log10(1.d-6*pny(upcoming_node))
c            pz = log10(1.d-6*pnz(upcoming_node))
c            write (isptr3, 445) cord(upcoming_node,1),
c     &           cord(upcoming_node,2), cord(upcoming_node,3),
            write (isptr3, 445) x3(iparticle), y3(iparticle),
     &           z3(iparticle),
     &           ttbtc(izone,iparticle), part_id(iparticle,1), 
     &           part_id(iparticle,2), izonef(upcoming_node), 
     &           upcoming_node, px
         else
c            write (isptr3, 445) cord(upcoming_node,1),
c     &           cord(upcoming_node,2), cord(upcoming_node,3),
            write (isptr3, 445) x3(iparticle), y3(iparticle),
     &           z3(iparticle),
     &           ttbtc(izone,iparticle), part_id(iparticle,1), 
     &           part_id(iparticle,2), izonef(upcoming_node), 
     &           upcoming_node
         end if
      else
         write (isptr3, 444) ttbtc(izone,iparticle),
     &        part_id(iparticle,1), part_id(iparticle,2), 
     &        izonef(upcoming_node), upcoming_node
      end if
      call flush(isptr3)

 444  format (1x, g21.14, 4(1x, i7))
 445  format (1x, 3(g16.9, 1x), g21.14, 4(1x, i7), 1(1x, g16.9))
      end subroutine write_alt_btc

******************************************************************
******************************************************************

      subroutine wtsi_ptrac3(np1)
c     s kelkar aug 29, 05
c     for water table nodes, if S(ijkv(np1))>Smin then set the initial
c     position below ddz*S, if S<Smin then move the particle vertically
c     downward until a node with S>Smin is encountered.

      use comdi, only : izone_free_nodes,s
      use comsptr
      use comsk

      implicit none

      integer inp1,np1,newnode
      real (kind=8) :: epsilontime = 1.e-10

      real*8 xp,yp,zp,dumm,zwt

      inp1 = ijkv(np1)
      zp=z1(np1)
      zwt=ddz(inp1)*s(inp1)
      if (deltawt .lt. 0.) then
         if (izone_free_nodes(ijkv(np1)) .eq. 3) then
            ioutt(np1) = -1
         else
            if(zp.gt.zwt) ioutt(np1) = -1
         end if
         return
      end if
      if(zp.ge.zwt) then
         xp=x1(np1)
         yp=y1(np1)
         newnode=inp1
         call wtsi_find_water(inp1,np1,xp,yp,zp,newnode)
         if (newnode .gt. 0) then
            call wtsi_displace_node(inp1,np1,xp,yp,zp,newnode)
            ijkv(np1)=newnode
            x1(np1)=xp
            y1(np1)=yp
            z1(np1)=zp
c Add epsilon time to particle time to indicate it has been moved
            ttp1(np1) =  ttp1(np1) + epsilontime
         else
            ioutt(np1) = -1
         end if
      endif

      end subroutine wtsi_ptrac3

******************************************************************
