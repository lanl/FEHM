      subroutine part_track(begin_time,end_time,hmon,lstep)
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
CD1
CD1  PURPOSE
CD1
CD1  This subroutine performs the particle tracking simulation.  
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              S. Henderson           Initial implementation.
CD2
CD2  $Log:   /pvcs.config/fehm90/src/part_track.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:36   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2 Made necessary changes and optimize the code for multiple species
CD2 simulations.				CLI, APR 18, 1997.
CD2
CD2    Rev 1.13   Wed Jun 12 09:49:06 1996   hend
CD2 Added abs check to pout option 3
CD2 
CD2    Rev 1.12   Wed May 08 14:16:02 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.11   Wed Mar 27 11:34:50 1996   hend
CD2 Updated Mixing Option (4)
CD2 Outputs now computed on entering a cell
CD2 
CD2    Rev 1.10   Tue Mar 05 09:34:22 1996   llt
CD2 Added new output option for radioactive Isotope Mixing (hend)
CD2 
CD2    Rev 1.9   Fri Mar 01 15:23:46 1996   hend
CD2 Pout=0 Option now Accounts for Cell Volume
CD2 
CD2    Rev 1.8   Tue Jan 16 14:06:56 1996   hend
CD2 Added output option for trace, pout=0
CD2 
CD2    Rev 1.7   Wed Jan 10 15:12:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.6   Mon Jan 08 10:43:50 1996   robinson
CD2 Algorithm no longer recomputes things once heat and mass is turned off
CD2 
CD2    Rev 1.5   09/06/95 10:01:46   robinson
CD2 New dispersion ellipsoid model and call to new dispersion routine
CD2
CD2    Rev 1.4   08/07/95 11:40:58   awolf
CD2 Fixed for DPDP (n0 instead of neq). Also changed for breakthrough output
CD2
CD2    Rev 1.3   04/03/95 08:47:42   robinson
CD2 Particles now can enter and leave system for nodes with ieos < 0
CD2 
CD2    Rev 1.2   03/16/95 16:10:48   llt
CD2 changed max function, so both arguments are double for IBM
CD2 
CD2    Rev 1.1   03/15/95 17:05:24   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.0   02/02/95 15:32:24   llt
CD2 made arguments for max function both double for IBM
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.5 Cell-based particle-tracking module
CD3  2.9   Interface with GoldSim
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
CPS BEGIN part_track
CPS
CPS IF (single_species)
CPS   IF(not_a_restart_run)
CPS     call set_ptrk to realease partciles
CPS     set restarting to true
CPS   ELSE
CPS     allocate arrays for re_starting run
CPS   ENDIF
CPS ELSE
CPS   IF(not_a_restart_run)
CPS     call inmptr to read in and set-up multiple species particle info.
CPS     set restarting to true
CPS   ENDIF
CPS ENDIF
CPS
CPS FOR each species, calculate the 10% value of number of particles
CPS     this value is used to calculate percentage of particles simulated.
CPS END FOR
CPS
CPS IF(array_sumdecayed not allocated)
CPS   allocate arrays for decay calculations
CPS ENDIF
CPS
CPS IF (GoldSim_coupled_run)
CPS   call getrip to get mass passed from GoldSim and allocated particles based 
CPS        information passed from GoldSim
CPS ENDIF
CPS
CPS FOR each species
CPS   call rearrange to re-arrange particle locations in memory to save space
CPS END FOR
CPS
CPS IF(heat_and_mass_run or first_time_step)
CPS   call massflow to calculate flux info. and retardation factor at each node 
CPS   set-up filtration flag, corresponding arrays, and pore size data
CPS ENDIF
CPS
CPS IF(decay-ingrowth calculation)
CPS   call ingrowth to perform decay and ingrowth calculation
CPS ENDIF
CPS
CPS IF(needs to output particle tracking info.)
CPS   write one line output indicating particle tracking will follow
CPS ENDIF
CPS
CPS
CPS FOR each species
CPS   call write_ptrk_info to print out particles tracking info
CPS   initialize pcount and pconcmax array for particle concentration calculation
CPS   FOR each segment of the particle array
CPS     write percentage of particles simulated
CPS     FOR each particle wihtin this segment
CPS       IF the particle is still in the system
CPS         IF the particle has already been injected into the system
CPS           call get_ptrk_params to derive particle related parameters
CPS           IF(transfer function approach for matrix diffusion)
CPS             call advect_move to move a particle into connected fracture/matrix node
CPS                  by advective flow only
CPS             call get_ptrk_params to derive particle related parameters at this new node
CPS             call f_m_mixing_trans_based_time to calculate diffusion related time and location
CPS           ELSE
CPS             IF (there is no outflow from the node)
CPS               record the current time and goes to next particle
CPS             ENDIF
CPS             IF(composite model run)
CPS               calculate the plug time for the fracture+matrix node
CPS             ELSE
CPS               calculate plug time for the particle at the current fractur or matrix node
CPS             ENDIF
CPS             call get_first_cell_time to derive the resident time of the particle
CPS           ENDIF
CPS           call set_ptrk_concens to calculate current concentration
CPS         
CPS           DO WHILE (current_time+timeleft < end_time)
CPS             update current_time so that current_time=current_time+timeleft
CPS             IF(composite model run)
CPS               call get_newloc_composite to determine the new location of a particle
CPS             ELSE
CPS               IF(transfer function approach for matrix diffusion)
CPS                 call f_m_mixing_trans_based_move to determine a new particle location
CPS               ELSE
CPS                 call get_newloc to determine a new particle location based on fraction of flow
CPS               ENDIF
CPS             ENDIF
CPS             reset frac_done
CPS             call check_filter to see whether a particle can be filtered or not
CPS             call check_partexit to check whether a particle exits the system or not
CPS             IF(this particle exits the system)goes to next particle
CPS             call get_ptrk_params to derive parameters related to this particle at the new location
CPS             call set_anl0 to update node concentration 
CPS             call function stuck_particle to check whether a particle is stucked or not
CPS             IF(the particle is not stuck)
CPS               IF(transfer function approach for matrix diffusion)
CPS                 call advect_move to move a particle into connected fracture/matrix node
CPS                                   by advective flow only
CPS                 call get_ptrk_params to derive particle related parameters at this new node
CPS                 call f_m_mixing_trans_based_time to calculate diffusion related time and location
CPS               ELSE
CPS                 IF(composite model run)
CPS                   calculate plug time based on fracture+matrix flux
CPS                 ELSE
CPS                   calculate plug time based on node flux data
CPS                 ENDIF
CPS                 call get_new_cell_time to calculate resident time at this new location
CPS               ENDIF
CPS             ENDIF
CPS             call set_anl to set concentrations in arrays
CPS           END DO WHIOLE LOOP
CPS           IF(a particle is not stuck and not moved out the system neither)
CPS             update frac_done to reflect time already spent in at a node
CPS           ENDIF
CPS           update current age of a particle and store it in timeleft
CPS         ENDIF a particle has not entered the system
CPS       ENDIF a particle is still in the system
CPS     END For each particle in this segment
CPS   END For each segment
CPS   IF(a particle has not entered the system)
CPS     set timeleft=end_time-start_time so that the age of a particle is negative
CPS   ENDIF
CPS   call decay_calc to perform single species decay calculation. Multiple species decay-ingrowth
CPS                   calculation is performed by ingrowth subroutine
CPS   call set_final_anl to calculate node concentration at the end of current time step for
CPS                       various output options
CPS   IF(ipzone >0)
CPS     call compute_ptrk_concs to calculate concentration within specified zones
CPS   ENDIF
CPS END FOR each species     
CPS END part_track           
CPS 
C***********************************************************************
c Subroutine for particle tracking, activated with macro ptrk or mptr. 
c Due to the choice of algorithm it is useful only for
c advection dominated problems. A residence time in a cell is 
c computed for a given particle as (mass of cell) / (total mass
c flow rate exiting cell). When this time is completed the particle
c changes cells probabilisticallly based on the exiting flow rates.
c 
c Date: Summer '94
c Author: Stephen Henderson, EES-5
c         Contact: Bruce Robinson, EES-4 
c                  George Zyvoloski, EES-5
c

      use compart
      use compfrac
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
      use comsk
      use comsptr,  only : tprpflag
      use comuserc, only : in

      implicit none

c ********** GLOBAL VARIABLES: **********
c timeleft(i) is the current residence time for particle i
c box(i) is the current node number of particle i
c ptrak is set to 1 if particle tracking is implemented
c trak_type is set to 1 for liquid tracking, 2 for vapor
c frac_done(i) records what fraction of the current node has been
c    completed -- the remaining residence time needs to be set 
c    accordingly
c start_time(i) records at what time in seconds particle i enters the
c    system
c pcnsk(i) records the number of particles to be injected at node i
c rseed stores the random number seed
c p is used in calculating new node to travel to (probability)
c flow_ot(i) stores the outgoing flow from node i
c Rf is used in a linear Sorption model
c pconc(i) is used for radioactive decay concentration

c ********** LOCAL VARIABLES: **********
c mass is used for the mass of any given node
c current_sum " "
c oldbox records the previous location of the particle -- this is
c    used in determining if the particle has exited the system
c i, i2, idpr, ifrearrange, and li are counters
c lbox is the box number used for nelm. This may differ from the 
c     actual box number by neq in dpdp.
c dpdpbox is the actual box number, equivalent to box(i)
c The variable dpdpbox is used to store the actual location (box) of the
c  particle. This is used to index material properties which are set for
c  all n0 elements in a dpdp problem. Lbox, however, always runs from 
c  zero to neq. If the problem is not dpdp, lbox will equal dpdpbox. 
c  Otherwise they will differ by neq, and lbox will be used to
c  reference the nelm array, and then indices will be adjusted by 
c  add_fact.
c current_time records how long the current particle has been in the
c    system
c hmon is 1 if heat and mass solution is active
c lstep, number of time steps
c tempneed, indicator used for node concentration calculation
c neq1, equals neq+1
c f_m_box, fracture_matrix interaction term location in the a_axy array
c nsizep, nsizep_sf hold the size of probability array p and p_sf, 
c    respectively
c ith, counter for species index
c jj, jj_frac: represents node layer and fracture layer index
c return_flag, indicator. 0 a particle still in the system. 
c   1 a particle is out 
c frac records the numerical value of 10% of injected particles
c plugtime_dbl, double precision of plug_time
c flow_fracture,
c prob_frac

      integer i,j,k,i2,idpr,ifrearrange,li,lbox,dpdpbox,add_fact,oldbox
      integer oldlbox,cur_part,hmon, lstep,tempneed,neq1, f_m_box
      integer nsizep, nsizep_sf,ith, jj, jj_frac, return_flag
      integer las,las1
      integer :: lbox_trans = 0
c bhl_12/8/05
      integer dum_psv
c bhl_12/8/05

      integer, allocatable:: frac(:)
	
      real*8 external pore_radius
      real*8 plugtime_dbl
      real*8 begin_time, end_time, thistime

      real plugtime,theta_disp,theta_diff,ran_sp
      real satmat,timecell, current_time, matpor

      real, parameter:: yearindays=365.25, dayinsecs=86400.

      logical filtered, filt_option, stuck
      logical sigma_flag, omega_flag, par3_flag
      real weight(4)
      integer points(4)

c....................................................................
c 3/15/07 s kelkar importace sampling for colloid diversity
      integer nsegs_temp,np_temp_max,np_temp, n_col_irv,itemp
      real*8 , allocatable :: ret_weight_temp(:)
      real*8 , allocatable :: rcoll_div_temp(:)
c....................................................................




      tempneed=0
      neq1 = neq + 1

 !pass the time step, lstep, to set up the particle tracking. cli

   	
      if(nspeci.eq.1) then
         if(.not.restarting)then
            cur_part = 1
            call set_ptrk(0,lstep,cur_part)
            theta = -99999.
            restarting=.true.
            tempneed=1
            if (ripfehm .eq. 1 .and. prnt_rst .ge. 40) 
     &           ibox_offset = 2 * neq + 1
         else
            if(idpdp.eq.0) then
               nsizep = nelm(neq+1) - neq - 1
            else
               nsizep = (nelm(neq+1)-neq-1)*2
            endif
            nsizep_sf = nelm(neq+1) - neq - 1
            if(.not. allocated(p)) then
               allocate(p(nsizep),p_sf(nsizep_sf))
               p = -5
               p_sf = -5
            end if
         end if
      else
         if(.not.restarting)then
c     zvd 02-May-08
c     Set rseed here for use in inmptr for particle size distribution
c     during a GoldSim simulation
            if (ripfehm .eq. 1) rseed = int(in(4)+0.1)
            call inmptr(in(3))
            restarting = .true.
            tempneed = 1
         end if
      end if

cbhl  5/2/08, 
c zvd 6/9/08 changed prnt_rst flag to >= 40
c zvd 7/01/08 Moved following from set_mptr
	if (.not. allocated (bconf_sav)) then    
          if (abs(prnt_rst) .ge. 40 .and. ripfehm .eq. 1) then            
               allocate (bconf_sav(max_particles,nspeci))
	    else
	         allocate (bconf_sav(1,nspeci))
	    end if
      endif
cbhl  5/2/08
      
       if(.not.allocated(frac))then
         allocate (frac(nspeci))
      end if

      do ith=1,nspeci
         if(nsegs(ith)>=astep(ith))then
	    write(ierr,*)''
	    write(ierr,*)'--------------------Warning------------------'
	    write(ierr,*)'For Species:',ith,'number of segments may be'
	    write(ierr,*)'bigger than set in mptr'
            write(ierr,*)'Current segments:',nsegs(ith)
	    write(ierr,*)'Set in mptr astep=',astep(ith)
	    write(ierr,*)'---------------------------------------------'
	    write(ierr,*)''
         end if
         frac(ith)=max_particles/10
         if (frac(ith).eq.0) frac(ith)=1
      enddo


 !added getrip to get source term information from RIP. cli 

      if(.not.allocated(sumdecayed))then
         allocate(sumdecayed(nspeci))
         sumdecayed=0
      endif
!      write(iptty,*) 'ripfehm =',ripfehm
      if(ripfehm.eq.1) then
         call getrip

      endif
c....................................................................
c 3/15/07 s kelkar importace sampling for irreversible colloid 
c diversity model
c we need to allocate ret_weight and rcoll_div as 3-d arrays
c first see if irreversible colloid species are present
c         n_col_irv=0
c         do ith=1,nspeci
c            if(flag_diversity(ith).and.flag_col_irrev(ith))then
c               n_col_irv=n_col_irv+1
c            endif
c         enddo
c         if(n_col_irv.gt.0) then
c            nsegs_temp=maxval(nsegs)
c             if(.not.allocated(ret_weight)) then
c                allocate(ret_weight(max_particles,n_col_irv))
c                allocate(rcoll_div(max_particles,n_col_irv))
c             endif
c            
c            do ith=1,nspeci
c              do li=1,nsegs(ith)
c                  las1=aidex(ith)+li
c                  las=las1-1
c
c	            np_temp=nsport(las1)-lsport(las)+1
c                  if( np_temp.gt. np_temp_max) np_temp_max=np_temp
c              enddo
c
c              allocate(rcoll_div_temp(np_temp))
c              allocate(ret_weight_temp(np_temp))
c              if(flag_diversity(ith).and.flag_col_irrev(ith))then
c                  do li=1,nsegs(ith)
c                     las1=aidex(ith)+li
c                     las=las1-1
c                     np_temp = nsport(las1)-lsport(las)+1
c                     rcoll_div_temp=0.
c                     ret_weight_temp=0.
c                     call impsample_ptrk(tprpflag(ith),ith,np_temp,
c     1                    rcoll_div_temp,ret_weight_temp, np_temp_max)
c     		     itemp=0
c                     do i=lsport(las),nsport(las1)
c                        itemp=itemp+1
c                        ret_weight(i,irrevs(divs(ith)))=
c     1                       ret_weight_temp(itemp)
c                        rcoll_div(i,irrevs(divs(ith)))=
c     1                       rcoll_div_temp(itemp)
c                     enddo
c                  enddo
c               endif
c              deallocate (ret_weight_temp)
c              deallocate (rcoll_div_temp)
c            enddo
c         endif
c end importance sampling changes...................................

      do ith=1,nspeci
         sumdecayed(ith)=0
         call rearrange
      enddo

      if(hmon.eq.1.or.lstep.eq.1) then
         call massflow

c     Need to reset itf_curve array for transfer functions if the
c     flow field has just been reset

         if(allocated(itf_curve)) then
            itf_curve = 0
            wt_curve = 0.
         end if

 !Determine if the maximum pore size filtration option is invoked
 !filt_option = .false.
 !do ith = 1, nspeci
 !  if(filter_flag(ith).gt.1) then
 !    filt_option = .true.
 !  end if
 !end do
 
 !Allocate (if not already in a previous time step) and call routine to 
 !compute maximum pore size that particles can pass through
          
 !if(filt_option) then
 !  if(.not.allocated(vg_poresize)) allocate(vg_poresize(n0))
 !  do i = 1, n0
 !    vg_poresize(i) = 2.*pore_radius(i)
 !  end do
 !end if
      endif

 !add the ingrowth calculation here. 7/9/97, cli.

      if(nreac.ne.0)then
         call ingrowth(begin_time,end_time)
      endif

      if (iptty.ne.0) write(iptty,303)
 303  format(1x,'**** particle tracking ****')

c ********** MAIN LOOP OVER ALL PARTICLES **********
      do 100 ith=1,nspeci

 !write information to screen and files
         call write_ptrk_info1

 !Initialize pcount
         do idpr=1,ipzone2
!            write(iptty,*) 'idpr',idpr
            pcount(idpr,ith) = 0
         enddo
		
 !Initialize pconcmax

         pconcmax=0
         do 110 li=1,nsegs(ith)
            ifrearrange=0
            las1=aidex(ith)+li
            las=las1-1
            do 200 i=lsport(las),nsport(las1)
 !modified the pcount calculation to account for variable confactors cli
 
               if(mod(i,frac(ith)).eq.0. .and. iptty .ne. 0 ) then
                  write(iptty,304) i*100/max_particles
 304              format(1x,'Ptrk is ',i3,'% complete')
               end if

 !If box(i,ith)<0, the particle has already exited the system, then 
 !move to the next particle. If the particle already in the system 
 !(start_time<end_time) then move the particle, otherwise go to next
 !one.  
 
               if(box(i,ith).gt.0) then
                  ifrearrange=ifrearrange+1
                  if(start_time(i,ith).le.end_time)then
 !derive particle related parameters
                     call get_ptrk_params
	      
 !The following if block is based on the transfer function algorithm
 !Bruce proposed for handling matrix diffusion.

                     if(diffflag(jj_frac,ith).le.-5)then
 !initialize lbox_trans for signaling the start of a new particle run
 !this is only used with the transfer function based run
                        lbox_trans=0
 !f-m interaction based on advective flux
                        if(frac_done(i,ith).eq.0.)call advect_move 
!get parameters for this particle after advection 
                        call get_ptrk_params
 !get resident time based on transfer function                
                        call f_m_mixing_trans_based_time 
!                        write(iptty,*) 'i2c=',i
                     else
                        if(diffflag(jj_frac,ith).ge.3)then
 !If there is no flow exiting move on to next particle -- this one 
 !isn't going anywhere.  Store the age of the particle first
                           if(flow_ottot(lbox).eq.0.)then
                              timeleft(i,ith) = end_time - 
     &                             start_time(i,ith)
                              goto 200
                           end if
 !for complete mixing fracture and matrix are one block
                           plugtime = (mass(lbox)+mass(lbox+neq))
     &                          /max(1.d-30,flow_ottot(lbox))
                        else
 !If there is no flow exiting move on to next particle -- this one 
 !isn't going anywhere.  Store the age of the particle first
                           if(flow_ot(dpdpbox).eq.0.)then
                              timeleft(i,ith) = end_time - 
     &                             start_time(i,ith)
                              goto 200
                           end if
                           plugtime = mass(dpdpbox)/
     &                          max(1.d-30,flow_ot(dpdpbox))
                        end if
              
 !compute the current residence time, timecell

                        call get_first_cell_time
                     end if
	        
                     current_time=max(dble(start_time(i,ith)),
     &                    begin_time)
 !first time -- need to set concentrations
                     call set_ptrk_concs(tempneed)
		
 !While the particles can move out of the current node within the 
 !current heat and mass time step, do so. 
 !Set the current residence time of the particle

                     particle_transit: do while 
     2                    (current_time+timeleft(i,ith).le.end_time)
 !Add the residence time to how long the particle has been 
 !in the system during this call to part_track

                     current_time=current_time+timeleft(i,ith)

 !Determine new node location for particle different routines 
 !depending on whether the composite f/m model is being used or not.

                     if(diffflag(jj_frac,ith).ge.3) then
                        call get_newloc_composite(i2)
                     else		
                        if(diffflag(jj_frac,ith).le.-5)then
                           call f_m_mixing_trans_based_move
                        else
                           call get_newloc(i2)
                        end if
                     end if

 !check to see whether this one can filtered out
                     call check_filter	

 !Only a filtered particle will have - box at this point.
 !We skip this exiting particle section for filtered particles

 !Check if the particle has exited the system. If so, oldbox 
 !will equal box.  Set box to negative value and record the 
 !the time it arrived at the node in timeleft.

                     call check_partexit
                     if(return_flag.eq.1) goto 200

 !set-up node concentration calculation
                     call set_anl0
              
 !if this particle if filtered then bypass it 
                     if(box(i,ith)<0) goto 200

 !reset the time fraction a particle stayed at a node
                     frac_done(i,ith)=0.
                     call get_ptrk_params

 !Check if particle is stuck in the new node
 !stuck = stuck_particle(i, ith)  !this line is removed to be consistent with V2.20

 !Now that we are in a new node, find the residence time and set
 !the particles time left here.                 
 !if(.not.stuck)then !this line is removed to be consistent with V2.20
                     if(diffflag(jj_frac,ith).le.-5)then
 !f-m interaction based on advective flux
                        call advect_move 
 !get parameters for this particle after advection 
                        call get_ptrk_params
 !get resident time based on transfer function   
                        call f_m_mixing_trans_based_time               
                     else                  
 !considering composite model and analytical f-m diffusion model
                        if(diffflag(jj_frac,ith).ge.3)then
 !If there is no flow exiting move on to next particle -- this one 
 !isn't going anywhere.  Store the age of the particle first
                           if(flow_ottot(lbox).eq.0.)then
                              timeleft(i,ith) = end_time - 
     &                             start_time(i,ith)
                              goto 200
                           end if  
 !for complete mixing fracture and matrix are one block
                           plugtime = (mass(lbox)+mass(lbox+neq))
     &                          /max(1.d-30,flow_ottot(lbox))
                        else
 !If there is no flow exiting move on to next particle -- this one 
 !isn't going anywhere.  Store the age of the particle first
                           if(flow_ot(dpdpbox).eq.0.)then
                              timeleft(i,ith) = end_time - 
     &                             start_time(i,ith)
                              goto 200
                           end if  
                           plugtime = mass(dpdpbox)/
     &                          max(1.d-30,flow_ot(dpdpbox))
                        end if
              
 !compute the current residence time, timecell
                        call get_new_cell_time
                     end if
 !end if  !corresponding to if(.not.stuck.) be consistent with V2.20
 !Update concentrations in arrays
                     call set_anl 
                  enddo particle_transit
             
 !The particle is unable to leave the current node within the current
 !heat and mass time step. Therefore, compute the fraction of the 
 !cell which it can complete, and store this so the residence time 
 !can be modified appropriately after the next heat and mass step.
 !Only do the calculation if the particle is not stuck. hen
	
 !if(.not.stuck)then   !this line is removed to be consistent with V2.20
                  frac_done(i,ith)=frac_done(i,ith)+(1-frac_done(i,ith))
     +                 *((end_time-current_time)/timeleft(i,ith))
 !end if   !corresponding to if(.not.stuck)

 !Record the current age of the particle, store in timeleft so that
 !it can be output
                  timeleft(i,ith)=end_time-start_time(i,ith) 
               endif            !corresponds to start_time<end_time
            endif               !corresponds to box(i,ith)>0
 200     continue
c bhl_5/15/08
c zvd 6/9/08 changed prnt_rst flag to < 40
      if (abs(prnt_rst) .lt. 40 .or. ripfehm.ne.1) then
         if(ifrearrange.eq.0)lsport(las)=nsport(las1)+1
      endif   
c bhl_5/15/08
 110  continue

 !If a particle has not yet entered the system, set timeleft
 !so that the age is recorded as a negative age
   
      do i = 1, num_particles(ith)
         if(start_time(i,ith).gt.end_time) then
	    timeleft(i,ith) = end_time - start_time(i,ith)
         end if
      end do

 !Calculate decay for single species or multiple
      call decay_calc

 !Determine final concentrations for the various pout options
      call set_final_anl

 !determine concentrations of fluid exiting specified zones
      if(ipzone.gt.0) then
         call compute_ptrk_concs
      end if
	
 !write out concentration data
      call write_ptrk_concs

 100  continue

      if(ipartout.eq.1 .and. iout .ne. 0) then
         write(iout,*) 'Particle tracking transfer function'
         write(iout,*) '           Curve Density'
         write(iout,*) 'i      j      k      density'
         do i = 1, nump1
            do j = 1, nump2
               do k = 1, nump3
                  write(iout,3500) i,j,k,param_density(i,j,k)
 3500             format(1x,i10,1x,i10,1x,i10,1x,i10)
               end do
            end do
         end do
      end if

c Check range of parameters used for transfer functions
      if (pfrac_read) then
         sigma_flag = .false.
         omega_flag = .false.
         par3_flag = .false.
         if ((sigma_low(1) .gt. sigma_low(2)) .or.
     &        sigma_high(1) .lt. sigma_high(2)) sigma_flag = .true.
         if ((omega_low(1) .gt. omega_low(2)) .or.
     &        omega_high(1) .lt. omega_high(2)) omega_flag = .true.
         if ((par3_low(1) .gt. par3_low(2)) .or.
     &        par3_high(1) .lt. par3_high(2)) par3_flag = .true.
         write (ierr, 10)
         write (ierr, 11) sigma_low(1), sigma_high(1), 
     2        omega_low(1), omega_high(1), par3_low(1), par3_high(1)
         write (ierr, 12)
         write (ierr, 11) sigma_low(2), sigma_high(2), 
     2        omega_low(2), omega_high(2), par3_low(2), par3_high(2)
         if (sigma_flag .or. omega_flag .or. par3_flag) then
            write (ierr, 14)
         end if
 10      format ('Diffusion model defined parameter space:')
 11      format ('min sigma = ', g16.9, ' max sigma = ', g16.9, /,
     2        'min omega = ', g16.9, ' max omega = ', g16.9, /,
     3        'min par3  = ', g16.9, ' max par3  = ', g16.9) 
 12      format ('Range needed (min/max parameter values used):')
 14      format ('***** WARNING ***** ',
     2        'Found values outside of defined parameter space')
      end if
      if (allocated(frac)) deallocate(frac)

      return

      contains
      
      subroutine getrip
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine is actived within particle tracking when RIP is linked.
CD1  It is called from subroutine part_trak to get the source terms
CD1  then call set_ptrk to initialize the position of the particles.
CD1
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2  FEHM Version 2.0, SC-194 (Fortran 90) 
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.3.5 Cell-based particle-tracking module
CD3 2.6   Provide Input/Output Data Files
CD3 3.0   INPUT AND OUTPUT REQUIREMENTS
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
c
c     Contact: Bruce Robinson, EES-5
c              George Zyvoloski, EES-5
cli   modified the code to handle variable conversion factors.  2/6/98
	
      implicit none

      integer:: idaughter,ithp,istp,ip,j,k,memleft,lns,lns1,tmp_ptindex
      integer:: M_fine, N_large, nodes_left,selected_node,newly_failed
      integer:: trans_node, in_buffer,M_add_i,M_add_N,index_N_large
      integer:: index_mass,ioverlap,region,nodes_in_region,N_failed
      integer:: gspeci,index_ith_mass, index_fracflow, number_of_species
      integer:: index_in_species, nflow_frac, ibin, icount, ispecies
      integer:: lans,lans1,last_sign,lsn,lsn1
c bhl_5/15/08
      integer:: ithpn,ithp_old
c bhl_5/15/08
      integer,save::out_unit
      
c     the in array index for M_fine groups value

c     index_M_fine was 4 prior to V 2.23, but we have added 4
c     values before it in in() for the random number seeds and
c     the flag to determine if fracture fraction data are passed
c     BAR 2-9-2005

cHari water table rise modification 3/1/07
      integer, parameter:: index_M_fine= 8 

      !set allocatable arrays for storing failed nodes and regions data
      ! Moved to compart
C     integer, allocatable:: M_N_region(:), M_N_failed_nodes(:)
C     integer, allocatable:: selected_in_region(:)

      !real variables

      real:: yearsleft,scaleconpart,sumdecayout,sumdecayin
      real:: decayin,decayout,rescale,conpart,adjconpart,conpartd4
      real:: x_fine, y_fine, x_dis, y_dis, xy_dis, xy_least_dis 
      real:: ran_sp

      real*8:: pinmass,sumpinmass,tmp_pcnsk

      real*8, allocatable:: region_mass(:)
      
      !set up parameters. 

      real, parameter:: xhalf=0.5, xquarter=0.25
      real, parameter:: large_dis=1.E9,long_time=1.E20
      
      !logical parameter for inquire statement
      
      !logical:: unit_not_available		

      !get injection time in days

      t1sk(1) = begin_time/dayinsecs
      t2sk(1) = end_time/dayinsecs

c     For Version 2.23 and later, we now are retrieving the random
c     number seed from GoldSim for cases in which FEHM and GoldSim are 
c     coupled. In inptrk and inmptr, the random number seeds have been
c     reset to 0 for this case (values in the fehm input file are
c     overwritten). Here, we are checking to see if the values are
c     0 and using the GoldSim supplied values instead. This means
c     that the values only get get to in(4) and in(5) once, and
c     after that the sequence of random numbers are controlled by
c     the FEHM simulation. BAR 2-9-2005

      if(rseed.eq.0) then
c     force real*8 value to integer, accounting for roundoff
         rseed = int(in(4)+0.1)
      end if
      if(rseed_release.eq.0) then
         rseed_release = int(in(5)+0.1)
      end if

      !for random ebs release, get M_fine_group and N_large_group #
      
      M_fine=in(index_M_fine)
      index_N_large=index_M_fine+M_fine+M_fine+1
      N_large=in(index_N_large)
c     As of V 2.23, we now use a flag to decide whether there are
c     nflow_frac inputs of fracture fractional flow data to
c     handle, or if the array skips that input and proceeds
c     directly to number_of_species. This flag is used in the 
c     if block below. BAR 2-9-2005

      if(int(in(6)).eq.1) then
c     flow fraction data exists
         index_fracflow = index_N_large+N_large+1
         nflow_frac = int(in(index_fracflow))
         index_in_species=index_N_large+int(in(index_N_large))+
     2        nflow_frac+2
         number_of_species = int(in(index_in_species))
         index_mass=index_N_large+N_large+N_large*number_of_species+3
	 if(.not.allocated(fracflow)) then
            allocate(fracflow(number_of_species, N_large))
         endif
c     
c     Assign fracflow array values from in()
c     
         icount = 1
         do ibin = 1, N_large
	    do ispecies = 1, number_of_species
               fracflow(ispecies,ibin)=in(index_fracflow+icount)
               icount = icount + 1
	    enddo
         enddo

      else

c     no flow fraction data
c     Simply get number_of_species and assign default fracflow
c     values of 1

         index_in_species=index_N_large+int(in(index_N_large))+1
         number_of_species = int(in(index_in_species))
	 if(.not.allocated(fracflow))then
            allocate(fracflow(number_of_species, N_large))
	 endif
         index_mass=index_N_large+N_large+2
         fracflow = 1
      end if

      !get # of buffer infromation from RIP input array in()
c     index_in_species=index_N_large+int(in(index_N_large))+1
c     number_of_species = int(in(index_in_species))

c     if(number_of_species.gt.45)then 
c     nflow_frac = number_of_species
c     index_in_species=index_N_large+int(in(index_N_large))+
c     2      nflow_frac+2
c     number_of_species = int(in(index_in_species))
c     index_mass=index_N_large+N_large+N_large*number_of_species+3
c     index_fracflow = index_N_large+N_large+1
c     if(.not.allocated(fracflow)) then
c     allocate(fracflow(number_of_species, N_large))
c     endif
c     icount = 1

c     do ibin = 1, N_large
c     do ispecies = 1, number_of_species
c     fracflow(ispecies,ibin)=in(index_fracflow+icount)
C     Hari write(6,*)ispecies,ibin, 
C     Hari     2	icount,in(index_fracflow+icount)
c     icount = icount + 1
c     enddo
c     enddo
C     Hari	   write(6,*)'assigned fracflow'
c     else
c     if(.not.allocated(fracflow))then
c     write(6,*)'allocating fracflow'
c     write(6,*)number_of_species,N_large
c     allocate(fracflow(number_of_species, N_large))
c     endif
c     index_mass=index_N_large+N_large+2
c     fracflow = 1
c     endif
      
      
      
      if(in(index_mass) == 0.0.and.lstep /= 1)return
      gspeci=in(index_mass-1)
      index_mass=index_mass+3
C     Hari


      !Moved calculation outside loop to handle iterations other than
      !the first one
      !allocate arrays for EBS random release 

      M_add_N=M_fine+N_large
      if(.not. allocated(M_N_region))then
         allocate (M_N_region(M_add_N), M_N_failed_nodes(M_add_N))
         allocate (selected_in_region(N_large))
         M_N_region=0
         M_N_failed_nodes=0
         selected_in_region=0
      endif
      
      if(.not.allocated(region_mass))then
         allocate(region_mass(N_large))
         region_mass=0.
      endif

      !Select release nodes for random EBS release

      !Get M_fine group (x,y) and locate the node index and region
      
      if(lstep == 1) then
         in_buffer=index_M_fine
         do i=1,M_fine
	    ioverlap=0
	    in_buffer=in_buffer+1
	    x_fine=in(in_buffer)
	    in_buffer=in_buffer+1
	    y_fine=in(in_buffer)
	    xy_least_dis=large_dis
	    do j=1,N_large
               do k=insnode(j)+1,insnode(j+1)
                  x_dis=cord(ptindex(k),1)-x_fine
                  y_dis=cord(ptindex(k),2)-y_fine
                  xy_dis=sqrt(x_dis*x_dis+y_dis*y_dis)
                  if(xy_dis < xy_least_dis)then
                     M_N_failed_nodes(i)=-k
                     M_N_region(i)=j
                     xy_least_dis=xy_dis
                  endif
               enddo
	    enddo
            
      !check for nodes that may be selected twice due to close locations
            
	    do j=1,i-1
               if(M_N_failed_nodes(i) == M_N_failed_nodes(j))then
                  ioverlap=ioverlap+1
               endif
	    enddo

 !if found a node is selected twice, write error message, then, stop
            
	    if(ioverlap /= 0)then
               write(ierr,*)''
               write(ierr,*)
     2              'Error Message: For EBS random release in FEHM'
               write(ierr,*)
     2              'Code stoped in SUBROUTINE getrip in part_track'
               write(ierr,*)''
               write(ierr,*)'Found',ioverlap, '   Overlapped Nodes'
               write(ierr,*)'Overlap at (',x_fine,',',y_fine,' )'
               write(ierr,*)
     2              'Overlaped node:',ptindex(-M_N_failed_nodes(i))
               write(ierr,*)
     2              'Program terminated, reassign M_fine coordinates'
               write(iptty,*)''
               write(iptty,*)
     2              'Error Message: For EBS random release in FEHM'
               write(iptty,*)
     2              'Code stoped in SUBROUTINE getrip in part_track'
               write(iptty,*)''
               write(iptty,*)'Found',ioverlap, '   Overlapped Nodes'
               write(iptty,*)'Overlap at (',x_fine,',',y_fine,' )'
               write(iptty,*)
     2              'Overlaped node:',ptindex(-M_N_failed_nodes(i))
               write(iptty,*)
     2              'Program terminated,reassign M_fine coordinates'
               write(iptty,*)
     2              'See Error Message File: fehmn.err for detail'
	    endif
            
 !reshuffle the ptindex and pcnsk array to put selected M_fine 
 !nodes at the end of the arrays.
            
	    region=M_N_region(i)
	    selected_in_region(region)=0
	    do j=1,i-1
               if(M_N_region(j) == region)then
                  selected_in_region(region)=
     2                 selected_in_region(region)+1
               endif
	    enddo
	    selected_node=insnode(region+1)-selected_in_region(region)
	    trans_node=-M_N_failed_nodes(i)
	    tmp_ptindex=ptindex(selected_node)
	    ptindex(selected_node)=ptindex(trans_node)
	    ptindex(trans_node)=tmp_ptindex
	    tmp_pcnsk=pcnsk(selected_node)
	    pcnsk(selected_node)=pcnsk(trans_node)
	    pcnsk(trans_node)=tmp_pcnsk
	    M_N_failed_nodes(i)=-selected_node	   
         enddo
         
         !detect number of nodes selected in M_fine groups.
         
         do i=1,N_large
	    selected_in_region(i)=0
	    do j=1,M_fine
               if(M_N_region(j) == i)then
                  selected_in_region(i)=selected_in_region(i)+1
               endif
	    enddo
         enddo
      endif
      
      !For N_large groups: randomly select N_failed package nodes
      
      in_buffer=index_N_large
      do i=1, N_large
         in_buffer=in_buffer+1
         N_failed=in(in_buffer)
         M_add_i=M_fine+i
         nodes_in_region=insnode(i+1)-insnode(i)-selected_in_region(i)
         newly_failed=N_failed-M_N_failed_nodes(M_add_i)
         nodes_left=nodes_in_region-M_N_failed_nodes(M_add_i)
         if(N_failed < nodes_in_region)then 
	    M_N_region(M_add_i)=i
	    do j=1,newly_failed
               if(nodes_left > 1)then
                  do 
                     k=int(ran_sp(rseed_release)*nodes_left+xhalf)
                     if(k /= 0)exit
                  enddo
               else
                  k=1
               endif
               selected_node=insnode(i)+M_N_failed_nodes(M_add_i)+k
               trans_node=insnode(i)+M_N_failed_nodes(M_add_i)+1
               tmp_ptindex=ptindex(trans_node)
               ptindex(trans_node)=ptindex(selected_node)
               ptindex(selected_node)=tmp_ptindex
               tmp_pcnsk=pcnsk(trans_node)
               pcnsk(trans_node)=pcnsk(selected_node)
               pcnsk(selected_node)=tmp_pcnsk
               nodes_left=nodes_left-1
               M_N_failed_nodes(M_add_i)=M_N_failed_nodes(M_add_i)+1
	    enddo
         else
	    M_N_failed_nodes(M_add_i)=nodes_in_region
	    M_N_region(M_add_i)=i
         endif
      enddo	
      
      !determine number of particles to be injected into the system

      yearsleft=(enday-t1sk(1))/yearindays+xhalf
      conpart=(t2sk(1)-t1sk(1))/yearindays
      conpartd4=conpart*xhalf
      
      if(nspeci > 1)then
         do ith=1, nspeci
	    istp=1
	    idaughter=0
	    decayin=0.
	    decayout=0.
	    sumdecayin=0.
	    sumdecayout=0.
	    ithp=num_particles(ith)
c     write(6,*)'species ',ith
c     write(6,*)'number of particles',num_particles(ith)  
      !determine mass to be injected in a region for each species
c     Subtract 1 to get to the correct spot
	    index_ith_mass=index_mass+ith-1
 	    do i=1,N_large
               region_mass(i)=0
               do j=1,M_add_N
                  if(M_N_region(j) == i)then
		     region_mass(i)=region_mass(i)+
     &                    in(index_ith_mass+(j-1)*gspeci)
                  endif
               enddo
	    enddo
            
            !search for daughter species & rescale particles injected
            
	    do j=1,nreac
               if(ith.eq.ori(j))then
                  idaughter=obj(j)
               endif
	    enddo
	    if(t1sk(1) >= conftime(ith))then
               if(newconfactor(ith) /= 0.)then
                  lns1=nsegs(ith)+1
                  lans=aidex(ith)+nsegs(ith)
                  !cli added a statement for array checking
                  if(lns1.gt.astep(ith))then
                     write(ierr,*) 'Array size may be out of bound ',
     &                    'for :', ith
                     write(ierr,*) 'Begin time (day):',
     &                    begin_time/dayinsecs
                     write(ierr,*) 'Current segment size:', nsegs(ith)
                     write(ierr,*) 'Set for:',astep(ith)
                     write(ierr,*) 'Begin time:',
     &                    begin_time/(yearindays*dayinsecs)
                     write(ierr,*) 'End   time:',
     &                    end_time/(yearindays*dayinsecs)
                     write(ierr,*) 'Time step: ',lstep
                  end if
                  if(newconfactor(ith) <0.)then
                     confactor(lans)=-newconfactor(ith)
                     bconfactor(lans)=1.D0/confactor(lans)
                     lans1=lans
                     do j=lns1+1,astep(ith)-1
                        lans1=lans1+1
                        confactor(lans1)=confactor(lans)
                        bconfactor(lans1)=bconfactor(lans)
                     enddo
                 !cli made the change so that user could have the option
                 !of using updated confactors or the original one
                     lans1=lans1+1
                     confactor(lans1)=newconfactor(ith)
                     bconfactor(lans1)=bconfactor(lans)
                     ioconfactor(ith)=0
                  else
                     lans=aidex(ith)+nsegs(ith)
                     bconfactor(lans)=1.D0/newconfactor(ith)
                     lans1=lans
                     do j=lns1,astep(ith)
                        confactor(lans1)=newconfactor(ith)
                        bconfactor(lans1)=bconfactor(lans)
                        lans1=lans1+1
                     enddo
                     ioconfactor(ith)=0
                  end if
               else
                  ioconfactor(ith)=1
               endif
               conftime(ith)=long_time
            endif
c     bhl_12/8/05
c     if(ioconfactor(ith) == 1)then
	    if(ioconfactor(ith) == 1.and.cfraction(ith)<=1.0e-8)then
c     bhl_12/8/05
c     write(6,*)'using memleft'
 !cli. Start comment out the following lines to reinstate the use of p_fraction
 !cli memleft=max_particles-ithp
 !cli do j=1,nreac
 !cli if(idaughter > 0)then
 !cli if(idaughter == obj(j))then
 !cli if(num_particles(idaughter) > ithp)then
 !cli memleft=max_particles-num_particles(idaughter)
 !cli endif
 !cli endif
 !cli endif
 !cli if(ith == obj(j))then
 !cli sumdecayin=sumdecayin+sumdecayed(ori(j))
 !cli endif
 !cli enddo
 !cli sumdecayout=sumdecayed(ith)
 !cli if(sumdecayout >= conpartd4)then
 !cli decayout=xquarter
 !cli else
 !cli decayout=1-sumdecayout/conpart
 !cli endif
 !cli if(sumdecayin.ge.conpartd4)then
 !cli decayin=xhalf*conpartd4/sumdecayin
 !cli else
 !cli decayin=1-sumdecayin/conpart
 !cli endif
 !cli rescale=decayin
 !cli if(decayout < rescale)rescale=decayout
 !cli scaleconpart=rescale*memleft/yearsleft
 !cli adjconpart=scaleconpart*conpart/N_large+xhalf
 !cli if(adjconpart < 1.)adjconpart=1.
 !cli c   Subtract 1 to get to the correct spot
 !cli index_ith_mass=index_mass+ith-1
 !cli do j=1,M_add_N
 !cli region=M_N_region(j)	!cli_replaced itmp++ with M_N_region(j)	
 !cli pinmass=in(index_ith_mass+(j-1)*gspeci)
 !cli if(pinmass /= 0.D0)then
 !cli dum_p(istp)=adjconpart*pinmass/
 !cli &			        !cli region_mass(region)
 !cli region=region+1
 !cli if(dum_p(istp) == 0)dum_p(istp)=1
 !cli call set_mptr(region,M_N_failed_nodes(j),ith,
 !cli &		                istp,ithp,idaughter)	!cli_added M_N_failed_nodes(j) for random ebs
 !cli lns=nsegs(ith)+aidex(ith)-1
 !cli confactor(lns)=dum_p(istp)*gmol(ith)/pinmass
 !cli bconfactor(lns)=1.D0/confactor(lns)
 !cli endif
 !cli enddo
 !cli end of commented out lines	      
               
 !CHari pfraction part of the code crashes the full Las Vegas Model 4/20/01
 !CHari revert back to the old conservative way
               
 !cli reinstated the use of p_fraction in the code. cli, 02/14/2005
 !cli per Barry Lester's request I put a check to make sure that before 
 !the last time step, p_fraction is reset to 0.25 to prevent an ill-set 
 !p_fraction cause the code to relase too many particles at the last step
               
               if(p_fraction(ith)>xquarter)then
                  if((enday-t2sk(1))<=(t2sk(1)-t1sk(1)))then
                     p_fraction(ith)=xquarter
                     write(ierr,*)''
                     write(ierr,*)'Close to last time step, ',
     &                    'reset p_fraction to: ', xquarter
                  end if
               end if 
               
                                !check for memory left    
               memleft=max_particles-ithp
               do j=1,nreac
                  if(idaughter > 0)then
                     if(idaughter == obj(j))then
                        if(num_particles(idaughter)>ithp)then
                           memleft=max_particles - 
     &                          num_particles(idaughter)
                        end if
                     end if
                  end if
                  if(ith == obj(j))then
                     sumdecayin=sumdecayin+sumdecayed(ori(j))
                  endif
               enddo

 !check to make sure we do have positive number of particles left
               if(memleft<0)then
                  write(ierr,*)''
                  write(ierr,*)'For Species: ',ith
                  write(ierr,*)'Inside part_track getrip subroutine'
                  write(ierr,*)'memory left may not be enough for ', 
     &                 'left time'
                  write(ierr,*)'at time step:', lstep
                  write(ierr,*)'Begin time:',
     &                 begin_time/(yearindays*dayinsecs)
                  write(ierr,*)'End   time:',
     &                 end_time/(yearindays*dayinsecs)
                  write(ierr,*)'Consider increasing parent p_fraction'
                  write(ierr,*)'or reduce daughter p_fraction'
                  write(ierr,*)'reset memleft to positive value and ',
     &                 'coninue'
                  memleft=-memleft
               end if        
               sumdecayout=sumdecayed(ith)
               if(sumdecayout >= conpartd4)then
                  decayout=p_fraction(ith)
               else
                  decayout=1-sumdecayout/conpart
               endif
               if(sumdecayin >= conpartd4)then
                  decayin=p_fraction(ith)
               else
                  decayin=1-sumdecayin/conpart
               endif
               rescale=decayin
               if(decayout < rescale)rescale=decayout
               if(rescale > p_fraction(ith))rescale=p_fraction(ith)
               scaleconpart=rescale*memleft/yearsleft
               adjconpart=scaleconpart*conpart/N_large+xhalf
               if(adjconpart < 1.)adjconpart=1.
                                !Subtract 1 to get to the correct spot
               index_ith_mass=index_mass+ith-1
               do j=1,M_add_N
                  region=M_N_region(j) 
 !cli_replaced itmp++ with M_N_region(j)		
                  pinmass=in(index_ith_mass+(j-1)*gspeci)
                  if(pinmass /= 0.D0)then
                     dum_p(istp)=adjconpart*pinmass/region_mass(region)
                     region=region+1
c bhl_11/3/06
c                    if(dum_p(istp) == 0)dum_p(istp)=1
                     if(dum_p(istp) < min_part) dum_p(istp)=min_part
c bhl_11/3/06
 !check to make sure that nsegs is less than astep
                     if((nsegs(ith)+1)>= astep(ith))then
                        write(ierr,*)'Array size may be out of bound ',
     &                       'for',ith
                        write(ierr,*)'Begin time (day):',
     &                       begin_time/dayinsecs
                        write(ierr,*)'Current segment size:', nsegs(ith)
                        write(ierr,*)'Set size:',astep(ith)
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)
                        write(ierr,*)'Time step: ',lstep
                     end if
 !check if number of particles is going to be bigger than max_particles
                     if((ithp+dum_p(istp)) >= max_particles)then
                        write(ierr,*)'Max num. of particles could be ',
     &                       'reached'
                        write(ierr,*)'Current number of particles:',ithp
                        write(ierr,*)'Particles to be injected:',
     &                       dum_p(istp)
                        write(ierr,*)'in region:',M_N_region(j)
                        write(ierr,*)'Max. num particles set:',
     &                       max_particles
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs) 
                        write(ierr,*)'at time step:',lstep           
                     end if
                                !set up partciles
c bhl_5/15/08
                     ithp_old=ithp
c bhl_5/15/08
                     call set_mptr(region,M_N_failed_nodes(j),ith,
     &                    istp,ithp,idaughter)
                 if(flag_diversity(ith).and.flag_col_irrev(ith))then
                   lns1=aidex(ith)+nsegs(ith)
		   lns=lns1-1
      	           np_temp=nsport(lns1)-lsport(lns)+1
                    if( np_temp.gt. np_temp_max) np_temp_max=np_temp
                      allocate(rcoll_div_temp(np_temp))
                      allocate(ret_weight_temp(np_temp))
                      rcoll_div_temp=0.
                      ret_weight_temp=0.
                      call impsample_ptrk(tprpflag(ith),ith,np_temp,
     1                rcoll_div_temp,ret_weight_temp, np_temp_max)
           	      itemp=0
                      do i=lsport(lns),nsport(lns1)
                         itemp=itemp+1
                         ret_weight(i,irrevs(divs(ith)))=
     1                   ret_weight_temp(itemp)
                         rcoll_div(i,irrevs(divs(ith)))=
     1                   rcoll_div_temp(itemp)
                         enddo
                    deallocate (ret_weight_temp)
                    deallocate (rcoll_div_temp)
                 endif
     
 !cli_added M_N_failed_nodes(j) for random ebs
                     lns=nsegs(ith)+aidex(ith)-1
                     confactor(lns)=dum_p(istp)*gmol(ith)/pinmass
                     bconfactor(lns)=1.D0/confactor(lns)
                  endif
               enddo	      
c     bhl_12/8/05
c     else
	    elseif (cfraction(ith).le.1.0e-8) then
c     bhl_12/8/05
c     write(6,*)'using confactor'
                                !Subtract 1 to get to the correct spot
               index_ith_mass=index_mass+ith-1
               do j=1,M_add_N
                  region=M_N_region(j)
                  pinmass=in(index_ith_mass+(j-1)*gspeci)
                  if(pinmass /= 0.D0)then
 !check to make sure nsegs is less than astep
                     if((nsegs(ith)+1)>= astep(ith))then
                        write(ierr,*)'Array size may be out of bound ',
     &                       'for',ith
                        write(ierr,*)'Begin time (day):',
     &                       begin_time/dayinsecs
                        write(ierr,*)'Current segment size:',nsegs(ith)
                        write(ierr,*)'Set size:',astep(ith)
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)
                        write(ierr,*)'Time step: ',lstep
                     end if
 !cli made the change so that user have the choice of either using 
 !updated confactor or the original one
                     last_sign=aidex(ith)+astep(ith)-1
                     lans=aidex(ith)+nsegs(ith)
                     if(confactor(last_sign) <0.)then
                        confactor(lans)=-confactor(last_sign)
                     end if
                     dum_p(istp)=pinmass*confactor(lans)/gmol(ith)
                     region=region+1
c bhl_11/3/06
c                    if(dum_p(istp) == 0)dum_p(istp)=1
                     if(dum_p(istp) < min_part) dum_p(istp)=min_part
c bhl_11/3/06
 !check if number of particles is going to be bigger than max_particles
                     if((ithp+dum_p(istp)) >= max_particles)then
                        write(ierr,*)'Max num. of particles could be ',
     &                       'reached'
                        write(ierr,*)'Current number of particles:',ithp
                        write(ierr,*)'Particles to be injected:',
     &                       dum_p(istp)
                        write(ierr,*)'in region:',M_N_region(j)
                        write(ierr,*)'Max. num particles set:',
     &                       max_particles
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)	        
                        write(ierr,*)'at time step:',lstep            
                     end if
                     call set_mptr(region,M_N_failed_nodes(j),ith,
     &                    istp,ithp,idaughter) 

                 if(flag_diversity(ith).and.flag_col_irrev(ith))then
                   lns1=aidex(ith)+nsegs(ith)
		   lns=lns1-1
      	           np_temp=nsport(lns1)-lsport(lns)+1
                    if( np_temp.gt. np_temp_max) np_temp_max=np_temp
                      allocate(rcoll_div_temp(np_temp))
                      allocate(ret_weight_temp(np_temp))
                      rcoll_div_temp=0.
                      ret_weight_temp=0.
                      call impsample_ptrk(tprpflag(ith),ith,np_temp,
     1                rcoll_div_temp,ret_weight_temp, np_temp_max)
           	      itemp=0
                      do i=lsport(lns),nsport(lns1)
                         itemp=itemp+1
                         ret_weight(i,irrevs(divs(ith)))=
     1                   ret_weight_temp(itemp)
                         rcoll_div(i,irrevs(divs(ith)))=
     1                   rcoll_div_temp(itemp)
                         enddo
                    deallocate (ret_weight_temp)
                    deallocate (rcoll_div_temp)
                 endif


 !cli_added M_N_failed_nodes(j) for random ebs
                     confactor(lans)=dum_p(istp)*gmol(ith)/pinmass
                     bconfactor(lans)=1.D0/confactor(lans)
c bhl_5/15/08
c zvd 6/9/08 changed prnt_rst flag to >= 40
                     if (abs(prnt_rst) .ge. 40) then
                        do ithpn=ithp_old+1,ithp
                           bconf_sav(ithpn,ith)=bconfactor(lans)   
                        enddo
                     endif
c bhl_5/15/08
                  endif
               enddo
c     bhl_12/8/05 outer
c     bhl_12/8/05
	    else
c     write(6,*)'using both memleft and confactor'
c     bhl_12/8/05
               
               if(p_fraction(ith)>xquarter)then
                  if((enday-t2sk(1))<=(t2sk(1)-t1sk(1)))then
                     p_fraction(ith)=xquarter
                     write(ierr,*)''
                     write(ierr,*)'Close to last time step, ',
     &                    'reset p_fraction to: ',xquarter
                  end if
               end if 
               
c     bhl_12/8/05
 !check for non mass-based memory left    
               memleft=max_particles*(1.0-cfraction(ith))-
     &              (ithp-num_part_mem(ith))
c     bhl_12/8/05
               do j=1,nreac
                  if(idaughter > 0)then
                     if(idaughter == obj(j))then
                        if(num_particles(idaughter)>ithp)then
c     bhl_12/8/05
                           memleft=max_particles*(1.0-cfraction(ith))-
     &                          num_particles(idaughter)
c     bhl_12/8/05
                        end if
                     end if
                  end if
                  if(ith == obj(j))then
                     sumdecayin=sumdecayin+sumdecayed(ori(j))
                  endif
               enddo

 !check to make sure we do have positive number of particles left
               if(memleft<0)then
                  write(ierr,*)''
                  write(ierr,*)'For Species: ',ith
                  write(ierr,*)'Inside part_track getrip subroutine'
                  write(ierr,*)'memory left may not be enough for ',
     &                 'left time'
                  write(ierr,*)'at time step:', lstep
                  write(ierr,*)'Begin time:',
     &                 begin_time/(yearindays*dayinsecs)
                  write(ierr,*)'End   time:',
     &                 end_time/(yearindays*dayinsecs)
                  write(ierr,*)'Consider increase parent p_fraction or'
                  write(ierr,*)'reduce daughter p_fraction'
                  write(ierr,*)'reset memleft to positive value and ',
     &                 'coninue'
                  memleft=-memleft
               end if        
               sumdecayout=sumdecayed(ith)
               if(sumdecayout >= conpartd4)then
                  decayout=p_fraction(ith)
               else
                  decayout=1-sumdecayout/conpart
               endif
               if(sumdecayin >= conpartd4)then
                  decayin=p_fraction(ith)
               else
                  decayin=1-sumdecayin/conpart
               endif
               rescale=decayin
               if(decayout < rescale)rescale=decayout
               if(rescale > p_fraction(ith))rescale=p_fraction(ith)
               scaleconpart=rescale*memleft/yearsleft
               adjconpart=scaleconpart*conpart/N_large+xhalf
               if(adjconpart < 1.)adjconpart=1.
                                !Subtract 1 to get to the correct spot
               index_ith_mass=index_mass+ith-1
               do j=1,M_add_N
                  region=M_N_region(j) 
 !cli_replaced itmp++ with M_N_region(j)		
                  pinmass=in(index_ith_mass+(j-1)*gspeci)
                  if(pinmass /= 0.D0)then
                     dum_p(istp)=adjconpart*pinmass/region_mass(region)
                     region=region+1
c bhl_11/3/06
c                    if(dum_p(istp) == 0)dum_p(istp)=1
                     if(dum_p(istp) < min_part) dum_p(istp)=min_part
c bhl_11/3/06
 !check to make sure that nsegs is less than astep
                     if((nsegs(ith)+1)>= astep(ith))then
                        write(ierr,*)'Array size may be out of bound ',
     &                       'for', ith
                        write(ierr,*)'Begin time (day):',
     &                       begin_time/dayinsecs
                        write(ierr,*)'Current segment size:',nsegs(ith)
                        write(ierr,*)'Set size:',astep(ith)
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)
                        write(ierr,*)'Time step: ',lstep
                     end if
 !check if number of particles is going to be bigger than max_particles
                     if((ithp+dum_p(istp)) >= max_particles)then
                        write(ierr,*)'Max num. of particles could be ',
     &                       'reached'
                        write(ierr,*)'Current number of particles:',ithp
                        write(ierr,*)'Particles to be injected:',
     &                       dum_p(istp)
                        write(ierr,*)'in region:',M_N_region(j)
                        write(ierr,*)'Max. num particles set:',
     &                       max_particles
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs) 
                        write(ierr,*)'at time step:',lstep	        
                     end if
 !set up particles
c     bhl_12/8/05
c     call set_mptr(region,M_N_failed_nodes(j),ith,
c     &               istp,ithp,idaughter)	!cli_added M_N_failed_nodes(j) for random ebs
c     lns=nsegs(ith)+aidex(ith)-1
c     confactor(lns)=dum_p(istp)*gmol(ith)/pinmass
c     bconfactor(lns)=1.D0/confactor(lns)
c     bhl_12/8/05
                  endif
c     bhl_12/8/05
c     write(6,*)'using confactor to supplement particles'
 !Subtract 1 to get to the correct spot
                  index_ith_mass=index_mass+ith-1
c     bhl_12/8/05
                  region=M_N_region(j)
                  pinmass=in(index_ith_mass+(j-1)*gspeci)
                  if(pinmass /= 0.D0)then
 !check to make sure nsegs is less than astep
                     if((nsegs(ith)+1)>= astep(ith))then
                        write(ierr,*)'Array size may be out of bound ',
     &                       'for', ith
                        write(ierr,*)'Begin time (day):',
     &                       begin_time/dayinsecs
                        write(ierr,*)'Current segment size:',nsegs(ith)
                        write(ierr,*)'Set size:',astep(ith)
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)
                        write(ierr,*)'Time step: ',lstep
                     end if
 !cli made the change so that user have the choice of either using 
 !updated confactor or the original one
                     last_sign=aidex(ith)+astep(ith)-1
                     lans=aidex(ith)+nsegs(ith)
                     if(confactor(last_sign) <0.)then
                        confactor(lans)=-confactor(last_sign)
                     end if
c     bhl_12/8/05
                     dum_psv=dum_p(istp)
                     dum_p(istp)=dum_p(istp)+pinmass*cfraction(ith)*
     &                    confactor(lans)/gmol(ith)
c     bhl_12/8/05
                     region=region+1
                     if(dum_p(istp) == 0)dum_p(istp)=1
 !check if number of particles is going to be bigger than max_particles
                     if((ithp+dum_p(istp)) >= max_particles)then
                        write(ierr,*)'Max num. of particles could be ',
     &                       'reached'
                        write(ierr,*)'Current number of particles:',ithp
                        write(ierr,*)'Particles to be injected:',
     &                       dum_p(istp)
                        write(ierr,*)'in region:',M_N_region(j)
                        write(ierr,*)'Max. num particles set:',
     &                       max_particles
                        write(ierr,*)'Begin time:',
     &                       begin_time/(yearindays*dayinsecs)
                        write(ierr,*)'End   time:',
     &                       end_time/(yearindays*dayinsecs)           
                        write(ierr,*)'at time step:',lstep	 
                     end if
c bhl_5/15/08
                     ithp_old=ithp
c bhl_5/15/08
                     call set_mptr(region,M_N_failed_nodes(j),ith,
     &                    istp,ithp,idaughter) 

                 if(flag_diversity(ith).and.flag_col_irrev(ith))then
                   lns1=aidex(ith)+nsegs(ith)
		   lns=lns1-1
      	           np_temp=nsport(lns1)-lsport(lns)+1
                    if( np_temp.gt. np_temp_max) np_temp_max=np_temp
                      allocate(rcoll_div_temp(np_temp))
                      allocate(ret_weight_temp(np_temp))
                      rcoll_div_temp=0.
                      ret_weight_temp=0.
                      call impsample_ptrk(tprpflag(ith),ith,np_temp,
     1                rcoll_div_temp,ret_weight_temp, np_temp_max)
           	      itemp=0
                      do i=lsport(lns),nsport(lns1)
                         itemp=itemp+1
                         ret_weight(i,irrevs(divs(ith)))=
     1                   ret_weight_temp(itemp)
                         rcoll_div(i,irrevs(divs(ith)))=
     1                   rcoll_div_temp(itemp)
                         enddo
                    deallocate (ret_weight_temp)
                    deallocate (rcoll_div_temp)
                 endif


 !cli_added M_N_failed_nodes(j) for random ebs
                     confactor(lans)=dum_p(istp)*gmol(ith)/pinmass
                     bconfactor(lans)=1.D0/confactor(lans)
c bhl_5/15/08
c zvd 6/9/08 changed prnt_rst flag to >= 40
                     if (abs(prnt_rst) .ge. 40) then
                        do ithpn=ithp_old+1,ithp
                           bconf_sav(ithpn,ith)=bconfactor(lans)   
                        enddo
                     endif
c bhl_5/15/08
c     bhl_12/14/05
                     num_part_mem(ith)=num_part_mem(ith)+dum_psv
c     bhl_12/14/05
                  endif
               enddo
            endif
c     bhl_12/8/05 outer
            num_particles(ith)=ithp
         enddo
      else
         do ith=1,nspeci
c     Subtract 1 to get to the correct spot
	    index_ith_mass=index_mass+ith-1
	    do j=1,M_add_N
               cur_part=num_particles(1)
               num_particles(1)=in(index_ith_mass+j-1)*
     &              confactor(1)/gmol(ith)
               call set_ptrk(0,2,cur_part)
	    enddo
         enddo
      endif
      
      !deallocate the memory used by region_mass

      deallocate(region_mass)

      return
      end subroutine getrip

      !function pore_radius is used to calculate matrix pore radius in dual continua model
      real*8 function pore_radius(i)
        implicit none
        integer i
        real*8 conversion_factor
        real*8 cutoff_radius,cutoff_pcp
      
        parameter(cutoff_radius = 1.e4)      
        parameter(conversion_factor = 144.5)     
        parameter(cutoff_pcp = conversion_factor/cutoff_radius)

        !Routine uses capillary pressure and pore radius equation of
        !Marshall et al. (1996) to compute pore radius subject to
        !the Young Laplace equation
        !Equation: r = 2*gamma*cos(theta)/Pc
        !where r is the radius in meters
        !gamma is surface tension  = 0.07225 J/m2
        !cos(theta) is assumed to be 1 (contact angle theta = 0)
        !Pc is capillary pressure (fehm variable pcp is in Mpa: Pa = 10^6*MPa)

        !pore radius units converted to nm by multiplying by 1e9

        !pore radius(microns) = 1.e9*2*(0.07225)/(pcp*1e6) = 144.5/pcp
        !Assume maximum radius is 1e4 nm, cutoff pcp = 1.445d-2

        if(pcp(i).ge.cutoff_pcp) then
          pore_radius = conversion_factor/pcp(i)
        else
          !for low capillary pressures, use the cutoff value of 1.e4 nm
          pore_radius = cutoff_radius
        end if

        return
      end function pore_radius

 	!subroutine for rearrange the ***(nstep,nspeci) array so that
 	!the size of the corresponding array will not be exceeded.
	subroutine rearrange
	implicit none
	integer j,j1,j2,k,l_rearr,kk,kk1,kk2,lnsr,lnsr1

        integer np_temp,inp_temp

	k=0
	l_rearr=0
	do j=1,nsegs(ith)
	  lnsr1=aidex(ith)+j
	  lnsr=lnsr1-1
	  if(lsport(lnsr).gt.nsport(lnsr1))then
	    k=k+1
	  endif
	  l_rearr=l_rearr+1
	  if(k.ne.l_rearr)goto 999
	enddo

999	if(k.gt.0)then
	  kk1=0
	  do j=k+1,nsegs(ith)
	    j1=j+aidex(ith)-1
	    j2=j1+j1
	    kk=kk1+aidex(ith)
	    kk2=kk+kk
	    kk1=kk1+1
	    
          ivdt(kk)=ivdt(j1)
	    lsport(kk)=lsport(j1)
		  nsport(kk)=nsport(j1)
	    nsport(kk+1)=nsport(j1+1)
	    nprevd(kk)=nprevd(j1)
	    tmsport(kk2)=tmsport(j2)
	    tmsport(kk2-1)=tmsport(j2-1)
	    confactor(kk)=confactor(j1)
	    bconfactor(kk)=bconfactor(j1)

	  enddo
	  nsegs(ith)=kk1
	endif

	return
	end subroutine rearrange


      subroutine write_ptrk_info1
        implicit none

	  integer lans

	  if(iptty.ne.0)write(iptty,*)'Processing Species: ',ith
	  if(pout.le.0.and.pout.ne.-7)then
	    write(istrc,*)'Species: ',ith
	  endif
	  
	  if(nsegs(ith).ge.astep(ith))then
	    write(ierr,*)''
	    write(ierr,*)'Error: Size of the NSEGS array exceeded'
	    write(ierr,*)'Species: ',ith,' nsegs=',nsegs(ith),' Set for',
     &				   astep(ith)
	    write(ierr,*)'Begin time:',
     &                  begin_time/(yearindays*dayinsecs)
	    write(ierr,*)'End   time:',
     &                  end_time/(yearindays*dayinsecs)
		  write(ierr,*)'Time step: ',lstep
	    write(ierr,*)''
	  endif
	  
	  lans=aidex(ith)+nsegs(ith)
	  if(nsport(lans) > max_particles) then
	    write(ierr,*)''
	    write(ierr,*)'Warning: Size of particle [] may be exceeded'
	    write(ierr,*)'Species: ',ith,':   Particles injected ',
     &                  nsport(lans)
	  endif
      return
      end subroutine write_ptrk_info1

	!this function is used to test whether a particle is stuck or not
      logical function stuck_particle(ipart, ipspec)
      implicit none
      integer ipart, ipspec

      stuck_particle = .false.
      if(filter_flag(ipspec).gt.1) then
        if(vg_poresize(box(ipart,ipspec)).le.
     2	   partsize(ipart,sizes(ipspec))) then
          stuck_particle = .true.
          timeleft(ipart,ipspec) = 1.e30
        end if      
      end if
      return
      end function stuck_particle
	

      !this subroutine is used to ge the parmeters valuse related to a particle
      subroutine get_ptrk_params
      implicit none

      !If node is a matrix node, determine corresponding fracture
      !node and set factor that get you to the right position in the flux
      !array (add_fact)

      dpdpbox=box(i,ith)
      if (dpdpbox.gt.neq) then
         lbox=dpdpbox-neq
         add_fact=nelm(neq+1)-neq-1
      else
         lbox=dpdpbox
         add_fact=0
      endif

      !set the model index number, fracture and matrix

      jj = itrc(dpdpbox)
      jj_frac = itrc(lbox)


      !For no matrix diffusion, porosity of the appropriate
      !medium is set in matpor

c     For matrix diffusion, a continuum model requires the
c     matrix porosity to be set from matrix_por, and the 
c     saturation is assumed to be unity
c     For dpdp, we make sure we are using the matrix node
c     for setting the properties

      if( diffflag(jj,ith) .eq. 0 ) then
        !No diffusion, need to store porosity in matpor
        matpor = ps(dpdpbox)
      else
	  if(idpdp.eq.0) then
          !porosity of matrix set in matrix_por, assumed saturated
          matpor = matrix_por(jj)
          satmat = 1.
        else
          !composite model and overlapping cells model
          if(diffflag(jj_frac,ith).ge.3) then
            matpor = ps(lbox+neq)
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               satmat = s(lbox+neq)
            else
               satmat = 1.
            end if
          else
            !set parameters based on matrix values of dual perm model
            matpor = ps(lbox+neq)
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               satmat = s(lbox+neq)
            else
               satmat = 1.
            end if
          end if
        end if
      end if

      return
      end subroutine get_ptrk_params

      !subroutine set_ptrk_concs is used to record the node concentration
	!values for Cl36, C14 simulations
      subroutine set_ptrk_concs(tempneed)
        implicit none
        integer tempneed
        real getconc14, getconc36

        if ((tempneed.eq.1)) then
          if(abs(pout).eq.4) then
            anl(dpdpbox)=anl(dpdpbox)+1
            anv(dpdpbox)=anv(dpdpbox)+
     +           (getconc36(current_time+timeleft(i,ith))-
     +           getconc36(current_time))/
     +           (timeleft(i,ith))
          elseif(abs(pout).eq.6) then
            anl(dpdpbox)=anl(dpdpbox)+1
            anv(dpdpbox)=anv(dpdpbox)+
     +           (getconc14(current_time+timeleft(i,ith))-
     +           getconc14(current_time))/
     +           (timeleft(i,ith))
          elseif(pout.eq.5) then
            anl(dpdpbox)=anl(dpdpbox)+1
            anv(dpdpbox)=anv(dpdpbox)+partconc(i)
          end if
        endif
        return
      end subroutine set_ptrk_concs


      !subroutine used to get the cell resident time at the begining of the loop	
	subroutine get_first_cell_time
	implicit none
      
	real*8  rtemp

	if(volume(dpdpbox) .eq. 0. ) then
	  !zero volume cell, set resident time=0
        timeleft(i,ith) = 0.
      else
        !if(start_time(i,ith).le.end_time)then 
	    !Check if particle is stuck via filtration mechanism within a unit
	    !stuck = stuck_particle(i, ith)  !this line removed to be consistent with V2.20
	    !if(.not.stuck) then
	      if( theta(i,ith) .eq. -99999. ) then    
	        !compute diffusion RTTF factor
	        plugtime_dbl = plugtime
		      if(diffflag(jj_frac,ith).ge.3) then
	          !composite model get the delay factor
		        theta_disp = 1.
		        call get_comptime(theta(i,ith))
		      else
	          !dual permeability model get the delay factor
		        call get_difftime(plugtime_dbl,theta(i,ith))
		      end if
	      end if
	    !end if  !corresponds to (.not.stuck) this line removed to be consistent with V2.20
	    timecell = plugtime*theta(i,ith)
	    !calculate the amount time the particle needs to spent.
	    timeleft(i,ith)=timecell*(1.-frac_done(i,ith))
	  !end if
	end if
	return
	end subroutine get_first_cell_time

      
	!this subroutine calculate the resident time of a particle entering 
	!into a new node.
	subroutine get_new_cell_time
	  plugtime_dbl=plugtime	
        if(diffflag(jj_frac,ith).ge.3) then
	    !Composite model, dispersion calculation
          call get_disptime(plugtime,theta_disp,thistime)
	    !compute delay factor due to adsorption
          call get_comptime(theta_diff)
        else
          !for dual permeability model 
          if((dpdpbox.eq.oldbox+neq).or.
     2        (dpdpbox+neq.eq.oldbox)) then
            thistime = plugtime
            theta_disp = 1.
          else
            !Compute dispersion RTTF factor
            call get_disptime(plugtime,theta_disp, thistime)
          end if
	    !Get diffusion delay factor using current travel time of
	    !dispersed particle (thistime)       
	    call get_difftime(thistime,theta_diff)
        end if

        theta(i,ith) = theta_disp*theta_diff
        timecell = plugtime*theta(i,ith)
						
        timeleft(i,ith) = timecell
	  return
	end subroutine get_new_cell_time

      !subroutine get_comptime is used to calculate the particle overall retardation factor
      subroutine get_comptime(theta_comp)
        implicit none
	  integer jj
        real theta_comp

        !Compute theta based only on matrix properties
        jj = itrc(lbox+neq)
        if(flag_diversity(ith))then
           if(flag_col_irrev(ith))then
              theta_comp = (Rf(lbox+neq,ith)+
     2             kcoll(jj,ith)*rcoll_div(i,irrevs(divs(ith))))/
     3             (1.+kcoll(jj,ith))
           else
              theta_comp = (Rf(lbox+neq,ith)+
     2             kcoll(jj,ith)*mean_rcoll_reve(reves(divs(ith))))/
     3             (1.+kcoll(jj,ith))
           endif
        else
           theta_comp = (Rf(lbox+neq,ith)+
     2          kcoll(jj,ith)*rcoll(jj,ith))/(1.+kcoll(jj,ith))
        endif
        return
      end subroutine get_comptime


      !subroutine get_difftime is used to calculate diffusion delay for S-F approach for
	!matrix diffusion
      subroutine get_difftime(thistime,theta_diff)
        implicit none

        integer fm
        real theta_diff, fact_term, apwid, fracrd, spacing, snorm
        real*8 thistime, fracrd_dbl, conc_ret, theta_dbl
        real*8 sigmav, omegav, par3v

        if( diffflag(jj,ith) .eq. 1 ) then
          if( trak_type(ith) .eq. 1 ) then
            fact_term = (matpor*satmat)**2*
     2                  diffmfl(ith,jj)*Rf(dpdpbox,ith)
          else
            fact_term = (matpor*(1.-satmat))**2*
     2                  diffmfl(ith,jj)*Rf(dpdpbox,ith)
          end if
          apwid = aperture(jj)*(1.+kcoll(jj,ith))

          if(flag_diversity(ith))then
          
             if(flag_col_irrev(ith))then
                fracrd =(rd_frac(jj,ith)+
     2               kcoll(jj,ith)*rcoll_div(i,irrevs(divs(ith))))/
     3               (1.+kcoll(jj,ith))
             else
                fracrd =(rd_frac(jj,ith)+
     2               kcoll(jj,ith)*mean_rcoll_reve(reves(divs(ith))))/
     3               (1.+kcoll(jj,ith))
             endif
          else
             fracrd =(rd_frac(jj,ith)+
     2            kcoll(jj,ith)*rcoll(jj,ith))/(1.+kcoll(jj,ith))

          endif
          call time_diff(real(thistime),rseed,apwid,
     2                  fracrd,fact_term,theta_diff)

        !Sudicky and Frind model
        elseif( diffflag(jj,ith) .le. -1 ) then 
          apwid = aperture(jj)*(1.+kcoll(jj,ith))
          if(idpdp.eq.0) then
             if (irdof .ne. 13 .or .ifree .ne. 0) then
                spacing = abs(aperture(jj))/(ps(lbox)*s(lbox))
             else
                spacing = abs(aperture(jj))/(ps(lbox))
             end if
          else
            if(diffflag(jj,ith) .lt. -1) then
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  snorm = (s(lbox)-sresidual(jj))/(1.-sresidual(jj))
               else
                  snorm = (1.-sresidual(jj))/(1.-sresidual(jj))
               end if
               snorm = max(1.e-10,snorm,sresidual(jj))
               spacing = 2.*apuv1(lbox)*snorm**(-gamma_afm(jj))
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  apwid = apwid*s(lbox)/snorm**(1.+gamma_afm(jj))
               else
                  apwid = apwid/snorm**(1.+gamma_afm(jj))
               end if  
            else
               spacing = 2.*apuv1(lbox)
               if (irdof .ne. 13 .or. ifree .ne. 0) 
     &              apwid = apwid*s(lbox)
            end if
          end if
          if(flag_diversity(ith))then
             if(flag_col_irrev(ith))then
                fracrd_dbl=(rd_frac(jj,ith)+kcoll(jj,ith)*
     2             rcoll_div(i,irrevs(divs(ith))))/(1.+kcoll(jj,ith))
             else
                fracrd_dbl=(rd_frac(jj,ith)+kcoll(jj,ith)*
     2               mean_rcoll_reve(reves(divs(ith))))/
     3               (1.+kcoll(jj,ith))
             endif
          else
             fracrd_dbl=(rd_frac(jj,ith)+kcoll(jj,ith)*rcoll(jj,ith))/
     2            (1.+kcoll(jj,ith))
          endif

          omegav =matpor*satmat*sqrt(plugtime*Rf(dpdpbox,ith)*
     2            diffmfl(ith,jj))/(0.5*apwid)
          sigmav =sqrt(Rf(dpdpbox,ith)/(plugtime*diffmfl(ith,jj)))*0.5*
     2            (spacing-apwid)
          par3v = 1.
	    fm=1
          call time_delay(3,ith,lbox,sigmav,omegav,par3v,fm,
     2                    fracrd_dbl,rseed,plugtime_dbl,
     3                    conc_ret,theta_dbl)

cHari 01-Nov-06 for the colloid diversity model set itf_curve to 0
cHari since the closest curve needs to be found for
cHari every particle
          if(flag_col_irrev(ith))then
             itf_curve(ith,i,1) = 0
          endif
          theta_diff = theta_dbl/plugtime  
        else
        


          !theta_diff is the delay factor for pure sorption w/o matrix diffusion
           if(flag_diversity(ith))then
              if(flag_col_irrev(ith))then
                 theta_diff =(Rf(dpdpbox,ith)+kcoll(jj,ith)*
     2             rcoll_div(i,irrevs(divs(ith))))/(1.+kcoll(jj,ith))
              else
                 theta_diff =(Rf(dpdpbox,ith)+kcoll(jj,ith)*
     2                mean_rcoll_reve(reves(divs(ith))))/
     3                (1.+kcoll(jj,ith))
              endif
           else
              theta_diff =(Rf(dpdpbox,ith)+kcoll(jj,ith)*
     2             rcoll(jj,ith))/(1.+kcoll(jj,ith))
           endif
        end if
      
        return
      end subroutine get_difftime


	!subroutine get_disptime is used for dispersion calculation
      subroutine get_disptime(plugtime, theta_disp, thistime)
        implicit none
        integer tempoldbox
        real*8 xcoord_old,xcoord_new,ycoord_old,ycoord_new
        real*8 zcoord_old,zcoord_new,deltax, deltay, deltaz
        real pe, dimc,  plugtime, theta_disp, ran_sp
        real*8 thistime

        if( dispflag(jj,ith) .ne. 0 ) then
! If oldbox hasn't been set yet, start it at 1. This is fixed outside this rounte
!          if (oldbox.eq.0) oldbox = 1
          if (oldbox.gt.neq) then
            tempoldbox=oldbox-neq
          else
            tempoldbox=oldbox
          endif
          xcoord_old = cord(tempoldbox,1)
          xcoord_new = cord(lbox,1)
          ycoord_old = cord(tempoldbox,2)
          ycoord_new = cord(lbox,2)
          zcoord_old = cord(tempoldbox,3)
          zcoord_new = cord(lbox,3)
          deltax = xcoord_old-xcoord_new
          deltay = ycoord_old-ycoord_new
          deltaz = zcoord_old-zcoord_new
          if ( icnl .eq. 0 ) then
            pe=(deltax*deltax)/(tclx(1,jj)*tclx(1,jj))
            pe=pe+(deltay*deltay)/(tcly(1,jj)*tcly(1,jj))
            pe=pe+(deltaz*deltaz)/(tclz(1,jj)*tclz(1,jj))
            pe=sqrt(pe)
          else
            pe=(deltax*deltax)/(tclx(1,jj)*tclx(1,jj))
            pe=pe+(deltay*deltay)/(tcly(1,jj)*tcly(1,jj))
            pe=sqrt(pe)
          end if
          dimc=ran_sp(rseed)
          call time_disp2(pe,dimc,theta_disp)
          thistime = theta_disp*plugtime
        else
          theta_disp = 1.
          thistime = plugtime
        end if
        return
      end subroutine get_disptime


      !subroutine get_newloc_composite is used to select the node location for a particle to
	!move to in a composite model
      subroutine get_newloc_composite(i2)
        implicit none
        integer ii1, ii2, i2
        real r, probability, ran_sp

        !Select random number for use in determining new location
        r=ran_sp(rseed)
                  
        !Select new location
        ii1 = (nelm(lbox)-neq)
        ii2 = (nelm(lbox+1)-neq-1)

        particle_exit: do i2=ii1,ii2
          if (r.lt.p_sf(i2)) exit particle_exit
        enddo particle_exit
                  
        !Determine if particle is entering the new fracture or matrix node
        oldbox=box(i,ith)

        !Store old value of fracture node
        if(oldbox.le.neq) then
          oldlbox = oldbox
        else
          oldlbox = oldbox - neq
        end if
        
        box(i,ith)=nelm(i2+neq+1)

        !If particle left at a sink, set oldbox = box
        if(oldlbox.eq.nelm(i2+neq+1)) then
          oldbox = box(i,ith)
        end if
      return
      end subroutine get_newloc_composite

      !subroutine get_newloc is used to get the new node location for a particle
      subroutine get_newloc(i2)
        implicit none
        integer ii1, ii2, i2
        real r, ran_sp
        real*8 prob_fm,fact_excl

        !Select random number for use in determining new location
        r=ran_sp(rseed)
                  
        !Select new location
        ii1 = (nelm(lbox)-neq)+add_fact
        ii2 = (nelm(lbox+1)-neq-1)+add_fact
        filtered = .false.
        
        !Handle size exclusion option for colloidal particles
        if(fcoll(jj,ith).lt.1.) then
          if(add_fact.eq.0) then
            prob_fm = 1.-p(ii2)
            prob_fm = prob_fm*fcoll(jj,ith)
            if(p(ii2).gt.0.) then
              fact_excl = (1.-prob_fm)/p(ii2)
            else
              fact_excl = 1.
              if(r.gt.prob_fm) then
                filtered = .true.
              end if
            end if
          else
            fact_excl = 1.
          end if
        else
          fact_excl = 1.
        end if

        do i2=ii1,ii2
          if (r.lt.fact_excl*p(i2)) goto 10
        enddo
                  
        !moving between fracture and matrix
        oldbox=box(i,ith)
        if (dpdpbox.gt.lbox) then
          box(i,ith)=lbox
        else
          box(i,ith)=box(i,ith)+neq
        endif

        !Particle is filtered at f/m interface
        if(filtered) then
          box(i,ith)=-oldbox
          timeleft(i,ith) = -current_time
        end if
        goto 27

        !moving between nodes staying within either fracture or matrix
 10     continue
        oldbox=box(i,ith)
        if (dpdpbox.gt.lbox) then
          box(i,ith)=nelm(i2-add_fact+neq+1)+neq
        else
          box(i,ith)=nelm(i2-add_fact+neq+1)
        endif

 27     continue
      return
      end subroutine get_newloc


      !this subroutine is used to determine whether a colloid can be filter
	!If a particle is filtered, set box to negative value and record the 
      !simulation time that it arrived at the node in timeleft.            
      subroutine check_filter
      implicit none
      integer kindex1, kindex2, jbin
      real r, pore_size, ran_sp
      real*8  fact_excl

      !Check to see if particle can get filtered
      if(filtered)goto 27
      if(filter_flag(ith).eq.1) then
	  !if ftn_factor is negative, do the pore size option instead of a simple factor
        kindex1=istrw_cold(i2)
        if(ftn_factor(kindex1).lt.0.) then
          kindex2=itfcsize(kindex1)
          jbin=1
          r=ran_sp(rseed)
          do while(itfcprobsize(jbin,kindex2).lt.r)
            jbin=jbin+1
          end do

	    !Interpolate to get exact size          
	    pore_size=itfcporsize(jbin-1,kindex2)+
     2               (itfcporsize(jbin,kindex2)-
     2               itfcporsize(jbin-1,kindex2))*
     3               (r-itfcprobsize(jbin-1,kindex2))/
     4               (itfcprobsize(jbin,kindex2)-
     5               itfcprobsize(jbin-1,kindex2))
     
	    !Filtered if particle is as big or bigger than pore              
	    if(pore_size.le.partsize(i,sizes(ith))) then
            filtered = .true.
            box(i,ith)=-oldbox
            timeleft(i,ith) = -current_time
            goto 27
          end if
        else
          fact_excl=ftn_factor(istrw_cold(i2))
          r=ran_sp(rseed)
          if(r.gt.fact_excl) then
            filtered = .true.
            box(i,ith)=-oldbox
            timeleft(i,ith) = -current_time
            goto 27
          end if
        end if
      end if

      !Check if want special output to .trc file for tracing the particles 
      !through a certain node
      !now set for entering -- if want leaving exchange box(i) and oldbox
      
      if(pout.lt.0) then
        if(pout.eq.-7) then
c          write(istrc,*) i,oldbox,izonef(oldbox),timeleft(i,ith)
          write(istrc,*) i,oldbox,izonef(oldbox),timeleft(i,ith),
     &          current_time
        else
          if((box(i,ith).eq.nskw(1)).or.(box(i,ith)-neq.eq.nskw(1))) 
     +      write(istrc,*) box(i,ith),current_time-start_time(i,ith)
        end if
      endif
    
27    continue
      return
      end subroutine check_filter


      !this subroutine is used to check whether a particle has exited the system
      subroutine check_partexit
      implicit none
      integer inmatrix, infractu, tempbox, idpr, idpr2, itmpbox
      real*8 ret_weight_fac
      real*8 :: zmax = -9e30

      return_flag = 0

      !if a particle has exited the system, oldbox equal box.  
      if(box(i,ith).gt.0) then
c....................................................................
c 3/15/07 s kelkar importace sampling for colloid diversity
c multiplying bconfactor by ret_weight
         ret_weight_fac = 1.
         if(flag_diversity(ith).and.flag_col_irrev(ith))then
            ret_weight_fac = ret_weight(i,irrevs(divs(ith)))
         endif
         if(flag_col_daughter(ith) .ne. 0)then
            ret_weight_fac = ret_weight_daughter(i,divs_d(ith))
         endif
c....................................................................
         if(box(i,ith).eq.oldbox .or. 
     &        volume(box(i,ith)) .eq. 0.) then
            if(box(i,ith).gt.neq) then
               !update node mass data for matrix ndoe
               inmatrix=1
               tempbox = box(i,ith) -neq
               loop1: do idpr=1,ipzone
               if(izonef(tempbox).eq.idzone(idpr))then
                  idpr2=ipzone1+idpr
                  pcount(idpr2,ith)=pcount(idpr2,ith)+
     2                 ret_weight_fac*bconfactor(las)
                  pconcmax(box(i,ith)) = pconcmax(box(i,ith)) +
     2                 ret_weight_fac*bconfactor(las)
                  inmatrix=0
                  exit loop1
               endif
            enddo loop1
            if(inmatrix.eq.1)pcount(ipzone2,ith)= pcount(ipzone2,ith)
     &           +ret_weight_fac*bconfactor(las)
         else
            !update node mass data for fracture data
            infractu=1
            tempbox=box(i,ith)
            loop2: do idpr=1,ipzone
            if(izonef(tempbox).eq.idzone(idpr))then
               pcount(idpr,ith)=pcount(idpr,ith)+
     2              ret_weight_fac*bconfactor(las)
               pconcmax(tempbox)= pconcmax(box(i,ith)) + 
     2              ret_weight_fac*bconfactor(las)
               infractu=0
               exit loop2
            endif
         enddo loop2
         if(infractu.eq.1)pcount(ipzone1,ith)=pcount(ipzone1,ith)
     &        +ret_weight_fac*bconfactor(las)
      end if
              !this particle left the system, set box to negative value and 
              !record the simulation time that it arrived at the node in timeleft.
              if(box(i,ith).gt.neq) then
                itmpbox = box(i,ith) - neq
              else
                itmpbox = box(i,ith)
              end if
              zmax = max(zmax,cord(itmpbox,3))
! Increment count of # particles that have left the system at this node
              num_exit(box(i,ith),ith) = num_exit(box(i,ith),ith) + 1
              box(i,ith)=-box(i,ith)
              timeleft(i,ith) = current_time
              return_flag = 1
            endif
          end if
        return
      end subroutine check_partexit


      !subroutine advect_move is used to move particles between 
	!fracture or matrix nodes at the start of each node simulation
	!based on B.A. Robinson's transfer function algorithm.
	subroutine advect_move
	  implicit none
	  real r,pfm_advect,pmf_advect,ran_sp

	  !flow_ot is used to calculate the fraction of f-m exchange which
	  !is used to determine the probability of a particle moves into
	  !the corresponding fracture or matrix node by advection only.
        dpdpbox=box(i,ith)
	  if(dpdpbox.gt.neq)then
          lbox=dpdpbox-neq
	    f_m_box=(nelm(neq1)-neq1)*2+lbox
          if(a_axy(f_m_box).lt.0)then
	      pmf_advect=-a_axy(f_m_box)/flow_ot(dpdpbox)
	      r=ran_sp(rseed)
	      if(r.le.pmf_advect)then
		      oldbox=box(i,ith)
	        box(i,ith)=lbox
	      end if
	    end if
	  else
	    lbox=dpdpbox
	    f_m_box=(nelm(neq1)-neq1)*2+lbox
          if(a_axy(f_m_box).gt.0)then
	      pfm_advect=a_axy(f_m_box)/flow_ot(lbox)
	      !considering colloid size exclusing
	      if(fcoll(jj,ith).lt.1)then
		      pfm_advect=pfm_advect*fcoll(jj,ith)
            end if
	      r=ran_sp(rseed)
	      if(r.le.pfm_advect)then
	        oldbox=box(i,ith)
	        box(i,ith)=lbox+neq
	      end if
	    end if
	  end if
	end subroutine advect_move

	!subroutine f_m_mixing_trans_based_time
	!this subroutine is used to move particles between fracture
	!and matrix based on transfer function algorithm developed
	!by Bruce Robinson. This step only considers matrix diffusion
	!mixing caused by advection is taken care of in advec_move
	subroutine f_m_mixing_trans_based_time        
	implicit none

	integer fm, lbox_neq,tempoldbox
	real*8 par1v,par2v,par3v,concv,prob_con,ret_factor,tau_zero
	real*8 dummy_a, plugtime_f,plugtime_m,half_spacing,half_aperture
	real*8 prob_fm, timeleft_dbl, flowfrac, composite_time
        real r, ran_sp, snorm
        integer box_subst

        ret_factor=0.
        tau_zero=1.

	!calculate fracture aperture and spacing
        half_aperture=0.5*aperture(jj_frac)*(1.+kcoll(jj_frac,ith))
        if(diffflag(jj_frac,ith) .eq. -5) then
           if (irdof .ne. 13 .or. ifree .ne. 0) then
              snorm = (s(lbox)-sresidual(jj_frac))/
     &             (1.-sresidual(jj_frac))
           else
              snorm = (1.-sresidual(jj_frac))/
     &             (1.-sresidual(jj_frac))
           end if
           snorm = max(1.e-10,snorm,sresidual(jj_frac))
           half_spacing = apuv1(lbox)*snorm**(-gamma_afm(jj_frac))
           half_aperture = half_aperture/snorm
        else
           half_spacing = apuv1(lbox)
        end if
	  
        !calculate resident time of fracture
        if(flow_ot(lbox).eq.0.)then
           plugtime_f = 1.d20
        else
          plugtime_f=mass(lbox)/max(1.d-30,flow_ot(lbox))
        end if
	  !calculate resident time of matrix
	  lbox_neq=lbox+neq
	  if(flow_ot(lbox_neq).eq.0.)then
            plugtime_m = 1.d20
	  else
	    plugtime_m=mass(lbox_neq)/max(1.d-30,flow_ot(lbox_neq))
	  end if
	  !calculate the 3 dimensionless parameters for the type curve
c
c     Make sure that when diffusion is very low, we are computing
c     par1v and par2v based on the limiting values we chose for
c     the transfer function curves. This forces the search into
c     the correct part of the parameter space
c     Parameters used in current transfer functions are:
c     Dmin = 1.e-20 m2/s
c     B = 0.5 m
c     b = 5.e-4 m
c     mat porosity = 0.2
c     tauf = 1.5e6 s
c     Rf = 1
c     par1v min = (1.e-20)*(1.5e6)*1/(0.5^2*Rm)
c      = 6.e-14/Rm
c     par2v = (1.e-20)*(1.5e6)*0.2/(5.e-4*0.5*frac por*frac sat)
c      = 1.2e-11/(frac por*frac sat)
c
c     Another correction made to these calculations is to
c     handle the end-member case in which a large matrix retardation
c     factor is input. For the purposes of finding the correct
c     transfer function curve, the code pegs the matrix retardation
c     factor at 1000. This way, the code finds an appropriate 
c     transfer function curve, rather than compensating for the
c     high Rm with a low diffusion coefficient. The composite time
c     is computed with the original Rm, so that the correct travel
c     time is computed in the end.
c
c     Maximum Rm used in the transfer function curves was 1000, so
c     the dimensionless parameters use a max of 1000 in the statements
c     below.

          if(diffmfl(ith,jj_frac).lt.1.e-18) then
             par1v = 6.e-14/min(1000.d0,Rf(lbox_neq,ith))
             par2v = 1.2e-11/(ps(lbox)*s(lbox))
          else
             dummy_a=diffmfl(ith,jj_frac)*plugtime_f
     &            /(half_spacing*half_spacing)
             par1v=dummy_a*Rf(lbox,ith)/min(1000.d0,Rf(lbox_neq,ith))
             par2v=dummy_a*ps(lbox_neq)*s(lbox_neq)*half_spacing
     &            /(ps(lbox)*s(lbox)*half_aperture)
          end if

c     BAR 8-26-2005
c     Change to 99.99% to accomodate the new transfer functions
c     generated up to 0.9999
c     ZVD 12-07-2005 Implemented with input parameter ffmax 

c     Make sure we are computing par3v based on 99% fracture flow
c     even if the actual percentage is higher. This ensures
c     that transfer functions at 99% fracture flow will be
c     used as the end member

          flowfrac = flow_ot(lbox)/(flow_ot(lbox)+flow_ot(lbox_neq))
          if(flowfrac.gt.ffmax) then
c     Peg the value of par3v for purposes of searching for the
c     correct transfer function curve
             par3v = par1v*(1.-ffmax)/(ffmax*par2v)
          else
             par3v=plugtime_f*Rf(lbox,ith)/
     2            (plugtime_m*min(1000.d0,Rf(lbox_neq,ith)))
          end if

            !generate random number for deciding whether a particle can diffuse 
            !into a different media and calculate resident time
          r=ran_sp(rseed)
          if(box(i,ith).le.neq)then 
             box_subst = box(i,ith)
             if(frac_done(i,ith).ne.0.)then
            !this is a returning particle skip the f-m moving process
                fm=1
             else
                prob_fm=fconc(ith,box_subst,par1v,par2v,par3v,2)
                                !handling colloid size exclusion
                if(fcoll(jj,ith).lt.1.)then
                   prob_fm=prob_fm*fcoll(jj,ith)
                end if
                if(r.le.prob_fm)then
                                !this particle will move into the matrix
                   if (nump3 .ne. 1) then
                      fm=2
                   else
                      fm=1
                   end if
                   oldbox=box(i,ith) 
                   box(i,ith)=box(i,ith)+neq
                   add_fact=nelm(neq+1)-neq-1
                   call get_ptrk_params
                else
                   fm=1
                end if
             end if 
          else
             box_subst = box(i,ith) - neq
             if(frac_done(i,ith).ne.0.)then
                                !this is a returning particle skip the m-f moving process
                fm=4
                add_fact=nelm(neq1)-neq-1
             else
                
                prob_fm=fconc(ith,box_subst,par1v,par2v,par3v,3)
                if(r.le.prob_fm)then   
                                !this particle will move from matrix to fracture
                   fm=3
                   oldbox=box(i,ith) 
                   box(i,ith)=box(i,ith)-neq
                   call get_ptrk_params
                else
                   fm=4
                   add_fact=nelm(neq1)-neq-1
                end if
             end if
          end if
c
c     Composite residence time is used in the normalization
c

          composite_time = (mass(lbox)*Rf(lbox,ith)+
     2         mass(lbox_neq)*Rf(lbox_neq,ith))
     3         /(flow_ot(lbox)+ flow_ot(lbox_neq))
c     If flow is virtually all in the matrix, used composite time,
c     which reverts to the matrix time
          if(flowfrac.lt.1.e-4) then
             timeleft(i,ith)=composite_time
          else
             timeleft_dbl = timeleft(i,ith)
             call time_delay(3, ith, box_subst, 
     2            par1v, par2v, par3v, fm, 
     3            ret_factor,rseed,tau_zero, concv, timeleft_dbl)
             
cHari 01-Nov-06 for the colloid diversity model set itf_curve to 0
cHari since the closest curve needs to be found for
cHari every particle
           
             if(flag_col_irrev(ith))then
                itf_curve(ith,i,1) = 0
             endif
         

c     Correct for the slight possibility that interpolation errors
c     give rise to negative travel times through the cell - BAR 6-21-2005

             timeleft_dbl = max(0.d0, timeleft_dbl)
             timeleft(i,ith)=(composite_time-plugtime_f*Rf(lbox,ith))*
     2            timeleft_dbl+plugtime_f*Rf(lbox,ith)
          end if
	    !decide whether disperson should be calculated
	    if(lbox_trans .ne.0)then
            if(lbox.ne.lbox_trans)then
	        tempoldbox=oldbox
	        oldbox=lbox_trans
              call get_disptime(timeleft(i,ith),theta_disp,thistime)
	        oldbox=tempoldbox
	        lbox_trans=lbox
	      else
	        theta_disp=1.  
	      end if
	    else
	      theta_disp=1.
	    end if
          timeleft(i,ith)=timeleft(i,ith)*theta_disp*
     2         (1-frac_done(i,ith))
          return    	       	 
	end subroutine f_m_mixing_trans_based_time

! subroutine fconc is used to interpolate the type curve final 
! concentrations based on the dimensionless parameters related to a 
! particle in the field.
        real*8 function fconc(ith,cur_node,parv1,parv2,parv3,fm)
c        use compfrac, only : curve_structure, numparams, itf_curve
        implicit none
	integer ilow,ihigh,jlow,jhigh,klow,khigh,fm
	real*8 parv1,parv2,parv3
        
	real*8 a1, a2, av, b1, b2, bv, c1, c2, cv
        real*8 t1, t2, t3, t4, t5, t6, t7, t8
        real*8 d1, d2, d3, d4, d5, d6, d7, d8, totald
        integer ith, cur_node

        
c     Transfer functions in free format requires finding the closest
c     one (if it hasn't been found already) and setting the indices
c     accordingly
        
        if(curve_structure.eq.1) then
           if(itf_curve(ith,cur_node,1).eq.0) then
c     Nearest neighbor search
              call find_closest_curve(0,parv1,parv2,parv3,
     2             ilow,jlow,klow,points,weight)
              itf_curve(ith,cur_node,1) = ilow
              t1=conc(ilow,jlow,klow,fm,nump(ilow,jlow,klow,fm))
              fconc = t1
              return
           else
              ilow = itf_curve(ith,cur_node,1)
              jlow = 1
              klow = 1
              t1=conc(ilow,jlow,klow,fm,nump(ilow,jlow,klow,fm))
              fconc = t1
              return
           end if
        elseif(curve_structure.gt.1) then
           if(itf_curve(ith,cur_node,1).eq.0) then
c     SVD scheme to interpolate between points
c     Routine finds the points and the weights
              call find_closest_curve(1,parv1,parv2,parv3,
     2             ilow,jlow,klow,points,weight)

              itf_curve(ith,cur_node,1) = points(1)
              itf_curve(ith,cur_node,2) = points(2)
              itf_curve(ith,cur_node,3) = points(3)
              wt_curve(ith,cur_node,1) = weight(1)
              wt_curve(ith,cur_node,2) = weight(2)
              wt_curve(ith,cur_node,3) = weight(3)
              if(numparams.gt.2) then
                 itf_curve(ith,cur_node,4) = points(4)
                 wt_curve(ith,cur_node,4) = weight(4)
              end if
           
c     Points are already found, just do the interpolated mean

           end if
           t1 = 0.
           t1=t1+wt_curve(ith,cur_node,1)*
     2          conc(itf_curve(ith,cur_node,1),1,1,
     3          fm,nump(itf_curve(ith,cur_node,1),1,1,fm))
           t1=t1+wt_curve(ith,cur_node,2)*
     2          conc(itf_curve(ith,cur_node,2),1,1,
     3          fm,nump(itf_curve(ith,cur_node,2),1,1,fm))
           t1=t1+wt_curve(ith,cur_node,3)*
     2          conc(itf_curve(ith,cur_node,3),1,1,
     3          fm,nump(itf_curve(ith,cur_node,3),1,1,fm))
           if(numparams.gt.2) then
              t1=t1+wt_curve(ith,cur_node,4)*
     2             conc(itf_curve(ith,cur_node,4),1,1,
     3             fm,nump(itf_curve(ith,cur_node,3),1,1,fm))
           end if
           if(t1.lt.0.) then
              fconc = 0.
           elseif(t1.gt.1.) then
              fconc = 1.
           else
              fconc = t1
           end if
           return
         else
c     regular structure for parameters, do linear interpolation
c     with simple search algorithm
           call indices(ilow,ihigh,parv1,param1,nump1)
           call indices(jlow,jhigh,parv2,param2,nump2)
           if(nump3.gt.1)then
              call indices(klow,khigh,parv3,param3,nump3)
           else
              klow=1
              khigh=1
              fm=1
           end if
        
           a1=param1(ilow)
           a2=param1(ihigh)
           av=parv1
           
           b1=param2(jlow)
           b2=param2(jhigh)
           bv=parv2
           
           c1=param3(klow)
           c2=param3(khigh)
           cv=parv3
           
           t1=conc(ilow,jlow,klow,fm,nump(ilow,jlow,klow,fm))
           t2=conc(ilow,jlow,khigh,fm,nump(ilow,jlow,khigh,fm))
           t3=conc(ilow,jhigh,klow,fm,nump(ilow,jhigh,klow,fm))
           t4=conc(ilow,jhigh,khigh,fm,nump(ilow,jhigh,khigh,fm))
           t5=conc(ihigh,jlow,klow,fm,nump(ihigh,jlow,klow,fm))
           t6=conc(ihigh,jlow,khigh,fm,nump(ihigh,jlow,khigh,fm))
           t7=conc(ihigh,jhigh,klow,fm,nump(ihigh,jhigh,klow,fm))
           t8=conc(ihigh,jhigh,khigh,fm,nump(ihigh,jhigh,khigh,fm))

c           if (nump3 .eq. 1) fm = fmtmp

           d1 = dsqrt( (a1 - av)**2 +(b1 - bv)**2 + (c1 - cv)**2 )
           
           d2 = dsqrt( (a1 - av)**2 +(b1 - bv)**2 + (c2 - cv)**2)
           
           d3 = dsqrt( (a1 - av)**2 +(b2 - bv)**2 + (c1 - cv)**2)
           
           d4 = dsqrt( (a1 - av)**2 +(b2 - bv)**2 + (c2 - cv)**2)
           
           d5 = dsqrt( (a2 - av)**2 +(b1 - bv)**2 + (c1 - cv)**2)
           
           d6 = dsqrt( (a2 - av)**2 +(b1 - bv)**2 + (c2 - cv)**2)
           
           d7 = dsqrt( (a2 - av)**2 +(b2 - bv)**2 + (c1 - cv)**2)
           
           d8 = dsqrt( (a2 - av)**2 +(b2 - bv)**2 + (c2 - cv)**2)
           
           if (d1 .eq. 0) then
              fconc = t1
              return
           end if
           
           if (d2 .eq. 0) then
              fconc = t2
              return
           end if
           
           if (d3 .eq. 0) then
              fconc = t3
              return
           end if
           
           if (d4 .eq. 0) then
              fconc = t4
              return
           end if
           
           if (d5 .eq. 0) then
              fconc = t5
              return
           end if
           
           if (d6 .eq. 0) then
              fconc = t6
              return
           end if
           
           if (d7 .eq. 0) then
              fconc = t7
              return
           end if
           
           if (d8 .eq. 0) then
              fconc = t8
              return
           end if
           
           d1 = 1./d1
           d2 = 1./d2
           d3 = 1./d3
           d4 = 1./d4
           d5 = 1./d5
           d6 = 1./d6
           d7 = 1./d7
           d8 = 1./d8
           
           totald = 1. /(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8)

           fconc = totald * (d1*t1 + d2*t2 + d3*t3 + d4*t4
     .          + d5*t5 + d6*t6 + d7*t7 + d8*t8)
        end if
        return
        end function fconc


      !subroutine f_m_mixing_trans_based_move moves a particle based on flow
	!distribution into a connected fracture or matrix node
	subroutine f_m_mixing_trans_based_move
c     Changing to box(i,ith)
c        if(dpdpbox.gt.neq)then
        if(box(i,ith).gt.neq)then
	    call get_newloc_mm
        else
	    call get_newloc_ff
        end if
        return
	end subroutine f_m_mixing_trans_based_move

      !code used to redistribut particles from fracture to fracture node
      subroutine get_newloc_ff
      implicit none
      real  r, ran_sp 
      integer ii1, ii2

      !select random number for use in determining new location      
      r=ran_sp(rseed)
                  
      !Select new location      
      ii1 = (nelm(lbox)-neq)
      ii2 = (nelm(lbox+1)-neq-1)
	            
	r=r*p(ii2)
      do i2=ii1,ii2
        if(r.lt.p(i2))then
          oldbox=box(i,ith)
		  lbox_trans=lbox
          box(i,ith)=nelm(i2+neq+1)
	    exit
	  end if
      enddo

      return
      end subroutine get_newloc_ff


      !code used to redistribute particles from matrix to matrix node
      subroutine get_newloc_mm
      implicit none
      real  r, ran_sp
      integer ii1, ii2, i2

      !select random number for use in determining new location      
      r=ran_sp(rseed)
                  
      !select new location      
      ii1 = (nelm(lbox)-neq)+add_fact
      ii2 = (nelm(lbox+1)-neq-1)+add_fact
	                 
	r=r*p(ii2)
      do i2=ii1,ii2
        if(r.le.p(i2))then
          oldbox=box(i,ith)
		  lbox_trans=lbox
          box(i,ith)=nelm(i2-add_fact+neq+1)+neq
	    exit
	  end if
      enddo
      
	return
      end subroutine get_newloc_mm

      !subroutine set_anl0 is used to calculate total solute mass based on particle info.
      subroutine set_anl0
        implicit none
        integer itmp
        real*8 denom


        if (ripfehm .eq. 1 .and. prnt_rst .ge. 40 .and. 
     &       abs(box(i,ith)) .ge. ibox_offset) then
           dpdpbox=abs(box(i,ith)) - ibox_offset
        else
           dpdpbox=abs(box(i,ith))
        end if
                
        if (dpdpbox.gt.neq) then
          lbox=dpdpbox-neq
          add_fact=nelm(neq+1)-neq-1
        else
          lbox=dpdpbox
          add_fact=0
        endif	  

        if (pout.eq.0) then
          itmp=dpdpbox+npt(ith)
          if (trak_type(ith) .eq. 1) then
             if (irdof .ne. 13 .or. ifree .ne. 0) then
                denom = (rolf(dpdpbox)*
     2               ps(dpdpbox)*s(dpdpbox)*sx1(dpdpbox))
             else
                denom = (rolf(dpdpbox)*
     2               ps(dpdpbox)*sx1(dpdpbox))
             end if
            if(denom .ne. 0.) then
              an(itmp)=an(itmp)+1/denom
            else
              an(itmp) = 0.
            end if
            anv(itmp) = an(itmp)
          else
            denom = (rovf(dpdpbox)*
     2              ps(dpdpbox)*(1-s(dpdpbox))*sx1(dpdpbox))
            if(denom .ne. 0.) then
              an(itmp)=an(itmp)+1/denom
            else
              an(itmp) = 0.
            end if
            anv(itmp) = an(itmp)
          endif
        endif
        return
      end subroutine set_anl0

      !subroutine set_anl is used to calculate node concentration for Cl36 C14 simulations
      subroutine set_anl
        implicit none
        real getconc14, getconc36

        if (abs(pout).eq.4) then
          anl(dpdpbox)=anl(dpdpbox)+1
          anv(dpdpbox)=anv(dpdpbox)+
     +        (getconc36(current_time+timeleft(i,ith))-
     +        getconc36(current_time))/(timeleft(i,ith))
        elseif (abs(pout).eq.6) then
          anl(dpdpbox)=anl(dpdpbox)+1
          anv(dpdpbox)=anv(dpdpbox)+
     +        (getconc14(current_time+timeleft(i,ith))-
     +        getconc14(current_time))/(timeleft(i,ith))
        elseif(pout.eq.5) then
          anl(dpdpbox)=anl(dpdpbox)+1
          anv(dpdpbox)=anv(dpdpbox)+partconc(i)
        end if
        return
      end subroutine set_anl

      
      subroutine decay_calc
      implicit none
      integer inode, ipart,lns,lns1,ibindex

      !Store positions for avs output in concentration array (_con).
      !Currently the concentration is the number of particles in any 
      !given node divided by the total volume of the node.
      !Before calculating concentrations, allow particles to decay.

      if(abs(pout).ne.4.and.abs(pout).ne.6.and.pout.ne.5) then
        do inode=1,n0
          pconc(inode)=0.
        enddo
        
	if(nspeci.eq.1) then
	  !ptrk, one species - use kfact(1) for decay
           do li=1,nsegs(ith)
	      lns1=aidex(ith)+li
	      lns=lns1-1
              do ipart=lsport(lns),nsport(lns1)
                 if(box(ipart,ith).gt.0 .and.
     2		      start_time(ipart,ith).le. end_time) then
                    pconc(box(ipart,ith))=pconc(box(ipart,ith))+
     2                   exp(-kfact(1)*(end_time-start_time(ipart,ith)))
                 elseif(box(ipart,ith).lt.0 .and.
     2                   timeleft(ipart,ith).lt.0.) then
                    if (ripfehm .eq. 1 .and. prnt_rst .ge. 40 .and. 
     &                   abs(box(i,ith)) .ge. ibox_offset) then
                       ibindex = abs(box(ipart,ith)) - ibox_offset
                    else
                       ibindex = abs(box(ipart,ith))
                    end if
                    pconc(-box(ipart,ith))=pconc(-box(ipart,ith))+
     2                   exp(-kfact(1)*(end_time-start_time(ipart,ith)))
                 end if
              enddo
           enddo       
        else
           do li=1,nsegs(ith)
	      lns1=aidex(ith)+li
	      lns=lns1-1
c Why isn't the second variable nsport?
c              do ipart=lsport(lns),lsport(lns1)
              do ipart=lsport(lns),nsport(lns1)
                 if(box(ipart,ith).gt.0 .and.
     2		      start_time(ipart,ith).le. end_time) then
                    pconc(box(ipart,ith))=pconc(box(ipart,ith))+
     2                   gmol(ith)*bconfactor(lns)
                 elseif(box(ipart,ith).lt.0 .and.
     2                   timeleft(ipart,ith).lt.0.) then
                    if (ripfehm .eq. 1 .and. prnt_rst .ge. 40 .and. 
     &                   abs(box(i,ith)) .ge. ibox_offset) then
                       ibindex = abs(box(ipart,ith)) - ibox_offset
                    else
                       ibindex = abs(box(ipart,ith))
                    end if
                    pconc(-box(ipart,ith))=pconc(-box(ipart,ith))+
     2                   gmol(ith)*bconfactor(lns)
                 end if
              enddo
           enddo
        end if
      endif
           
      return
      end subroutine decay_calc


 	!this subroutine calculate node concentration at the end of each time step
      subroutine set_final_anl         
      implicit none
      integer inode, itmp

      !Output #/ Total Volume      
      if(abs(pout).eq.1) then
        do inode=1,n0
          itmp=inode+npt(ith)
          an(itmp)=pconc(inode)/sx1(inode)
          anv(itmp) = an(itmp)
        enddo
     
      !Output #/ kg fluid (liquid or vapor, depending on trak_type)
      else if (abs(pout).eq.2) then
        if( trak_type(ith) .eq. 1 ) then
          do inode=1,n0
            itmp=inode+npt(ith)
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               an(itmp)=pconc(inode)/max(1.d-20,rolf(inode)*ps(inode)*
     2              s(inode)*sx1(inode))
            else
               an(itmp)=pconc(inode)/max(1.d-20,rolf(inode)*ps(inode)*
     2              sx1(inode))
            end if
            anv(itmp) = an(itmp)
          enddo
        else
          do inode=1,n0
            itmp=npt(ith)+inode
            an(itmp)=pconc(inode)/max(1.d-20,rovf(inode)*ps(inode)*
     2		     (1-s(inode))*sx1(inode))
            anv(itmp) = an(itmp)
          enddo
        end if
     
      !Output Mixed Conc for Bomb Pulse Studies      
	else if (abs(pout).eq.4.or.abs(pout).eq.6.or.pout.eq.5) then
        do inode=1,n0
          if (anl(inode).eq.0) then
            an(inode)=1e-20
          else
            an(inode)=anv(inode)/anl(inode)
          endif
        enddo
     
	!Output total #                    
	else if (abs(pout).eq.3) then 
        do inode=1,n0
          itmp=npt(ith)+inode
          an(itmp)=pconc(inode)
          anv(itmp) = an(itmp)
        enddo

!Output total # currently in a cell or have left from given cell     
	else if (abs(pout).eq.10) then
! Zero out an and anv for new count
           do inode = 1, n0
              itmp=inode+npt(ith)
              an(itmp) = 0.
              anv(itmp) = 0.
           end do
! Count number currently in a cell or that have exited system from cell
           do inode = 1, num_particles(ith)
              if (ripfehm .eq. 1 .and. prnt_rst .ge. 40 .and. 
     &             abs(box(i,ith)) .ge. ibox_offset) then
                 itmp=abs(box(inode,ith)) - ibox_offset + npt(ith)
              else
                 itmp=abs(box(inode,ith)) + npt(ith)
              end if
              if (box(inode,ith) .lt. 0) then
                 if (timeleft(inode,ith) .ge. 0.) then
                    an(itmp)=an(itmp) + 1
                 end if
              else if (box(inode,ith) .gt. 0) then
                 anv(itmp) = anv(itmp) + 1
              end if
           enddo
      endif
           
      do inode = 1, n0
        if( volume(inode) .eq. 0. ) then
          itmp=inode+npt(ith)
          an(itmp) = 0.
          anv(itmp) = 0.
        end if
      end do
      return
      end subroutine set_final_anl

      !subroutine compute_ptrk_concs is used to calculate node concentration at the 
	!end of the current time step
      subroutine compute_ptrk_concs
        implicit none
        integer idpr, add_fact, tempbox
        real fluxnode
!        real*8 end_time, begin_time

        do i = 1, ipzone2
          idcavg(i,ith) = 0.
          idcmax(i,ith) = 0.
        end do

        !For first species, compute total water flux out of each zone
        !These values then are used for all species
        if(ith.eq.1) then
          idflux = 0.
          do i = 1, neq
            loop4: do idpr = 1,ipzone 
              if(izonef(i).eq.idzone(idpr)) then
                fluxnode = a_axy(nelmdg(i)-neq-1)
                if(fluxnode.gt.0.and.fluxnode.lt.1.e8) then
                  idflux(idpr) = idflux(idpr)+fluxnode
                end if
                exit loop4
              end if
            end do loop4
          end do

          !calculate matrix node
          do i = neq+1, n0
	      !loop through each zone
            loop5: do idpr = 1,ipzone 
              if(izonef(i-neq).eq.idzone(idpr)) then
                add_fact = nelm(neq+1)-neq-1
                fluxnode = a_axy(add_fact+nelmdg(i-neq)-neq-1)
                if(fluxnode.gt.0.and.fluxnode.lt.1.e8) then
                  idflux(ipzone1+idpr)=idflux(ipzone1+idpr)+
     2                                 fluxnode
                end if
                exit loop5
              end if
            end do loop5
          end do
        end if
      
	  !calculate con. for fracture node
        do i = 1, n0
          if(i.gt.neq) then
            tempbox = i-neq
            add_fact = nelm(neq+1)-neq-1
            fluxnode = a_axy(add_fact+nelmdg(tempbox)-neq-1)
            add_fact = ipzone1
          else
            tempbox = i
            fluxnode = a_axy(nelmdg(tempbox)-neq-1)
            add_fact = 0
          end if
          if(fluxnode.eq.0) then
            fluxnode = 1.e30
          end if
	    !loop through each zone
          loop3: do idpr = 1,ipzone
            if(izonef(tempbox).eq.idzone(idpr)) then
              idcmax(idpr+add_fact,ith) = 
     2            max(idcmax(idpr+add_fact,ith),gmol(ith)*pconcmax(i)/
     3            ((end_time-begin_time)*fluxnode))
              exit loop3
            end if
          end do loop3
        end do
      
        return
      end subroutine compute_ptrk_concs

      !subroutine write_ptrk_concs is used to output node concentration data
      subroutine write_ptrk_concs
      implicit none
      integer idpr, idpr2
      real, parameter:: yearindays=365.25, dayinsecs=86400.
 	
      if (iout .ne. 0) then
         write(iout, '("Species: ",i4)') ith
         write(iout,*)'at time: ',end_time/(yearindays*dayinsecs)
         if(idpdp.ne.0)then
	    !for dual permeability model output both fracture and matrix node con.
            write(iout,*) 'Mass leaving each zone'
            do idpr=1,ipzone
               idpr2=ipzone1+idpr
               write(iout,*) idzone(idpr),pcount(idpr,ith)*gmol(ith),
     &              pcount(idpr2,ith)*gmol(ith)
            enddo
          
            write(iout,*)'rest: ',pcount(ipzone1,ith)*gmol(ith),
     &           pcount(ipzone2,ith)*gmol(ith)
            write(iout,*) 'Average concentration leaving each zone'
            do idpr=1,ipzone
               idpr2=ipzone1+idpr
               if(idflux(idpr).le.0.) idflux(idpr) = 1.e30
               if(idflux(idpr2).le.0.) idflux(idpr2) = 1.e30
               idcavg(idpr,ith) =pcount(idpr,ith)*gmol(ith)/
     2              ((end_time-begin_time)*idflux(idpr))
               idcavg(idpr2,ith)=pcount(idpr2,ith)*gmol(ith)/
     2              ((end_time-begin_time)*idflux(idpr2))
               write(iout,*)idzone(idpr),idcavg(idpr,ith),
     2              idcavg(idpr2,ith)
            enddo
            
            write(iout,*) 'Maximum concentration leaving each zone'
            do idpr=1,ipzone
               idpr2=ipzone1+idpr
               write(iout,*)idzone(idpr),idcmax(idpr,ith),
     2              idcmax(idpr2,ith)
            enddo
         else
	    !output node concentration for single continum model 
            write(iout,*) 'Mass leaving each zone'
            do idpr=1,ipzone
               write(iout,*)idzone(idpr),pcount(idpr,ith)*gmol(ith)
            enddo
            write(iout,*)'rest: ',pcount(ipzone1,ith)*gmol(ith)
            write(iout,*) 'Average concentration leaving each zone'
            do idpr=1,ipzone
               if(idflux(idpr).le.0.) idflux(idpr) = 1.e30
               idcavg(idpr,ith)= pcount(idpr,ith)*gmol(ith)/
     2              ((end_time-begin_time)*idflux(idpr))
               write(iout,*)idzone(idpr),idcavg(idpr,ith)
            enddo
            write(iout,*) 'Maximum concentration leaving each zone'
          
            do idpr=1,ipzone
               write(iout,*)idzone(idpr),idcmax(idpr,ith)
            enddo
         endif
      end if
      return
      end subroutine write_ptrk_concs

      end subroutine part_track
