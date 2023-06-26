      subroutine userr(iz)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To allow for user specified changes in distcoeff by zone as 
CD1 a function of time.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 3.3
CD2
CD2 Initial implementation: 15-SEP-2022, Programmer: J. Ortiz 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/userc.f_a  $
!D2 
!D2    Rev 2.0   Fri May 07 14:47:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C***********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 N/A
CD3 
C***********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 This user subroutine has one option:
CD4    1 -- Use a given set of concentration v. time data to control
CD4         the injection concentration with time. To use this option
CD4         the time data must range from 0 seconds ago to x seconds 
CD4         ago, and the negative inlet concnetration must be the time
CD4         ago in years to start the simulation. 
CD4
C***********************************************************************
CD5
CD5 Local Subprograms
CD5
CD5 userr_model_setup 
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      use comci
      use compart
      use comrxni, only : scl, dsccl
      use comuserr  !jpo
      implicit  none

      integer jl,ju,jm,iz,i2
c      real*8 slope,getconc,timein,getflux
c      integer, optional :: i               !JPO
c      real(8), optional :: rc_ss, drc_ss   !JPO
      integer ispecies, j, k, k1, curcolumn
      integer, allocatable :: nsindex(:)
      integer nsnodes, iread1, nszones, mi, mim
c     SPC
      integer flag_kd2, irip, jjj, sflag
      real*8  deltat, in3_old
      
      integer flag_kd, i1,j1, ii, jj,  tf1
      
      logical used
      integer icount, open_file
      character*100 filename
      integer imsg(16)
      real*8 xmsg(16)
      integer msg(16)
      integer nwds
      character*32 cmsg(16)
      character*80 input_msg
      character*80 input_msg1
      character key_word*80  
      integer :: num_models, num_times





      if (iz.eq.0) then
            if (iout .ne. 0) write (iout,*)
     &         'rxn userr subroutine is invoked'
            if (iatty.gt.0) write(iatty,*)
     &         'rxn userr subroutine is invoked'
            filename = ''
c           Allow user to enter name of file for data
            read (inpt,'(a100)') filename
            if (filename(1:4) .eq. 'file') then
               read (inpt,'(a100)') filename
            else
               backspace inpt
               filename = ''
               filename(1:14) = 'userr_data.dat' !default 
            end if

       
            iread1 = open_file(trim(filename),'old')
            if (iout .ne. 0) write (iout,*)
     &         'rxn userr subroutine read from optional input file: ',
     &         filename
            if (iatty.gt.0) write(iatty,*)
     &         'rxn userr subroutine read from optional input file: ',
     &         filename

c           Read 1st line of userr_data.dat (NUM_TIMES, NUM_MODELS)
            read(iread1,*)num_times,num_models

            call userr_model_setup(iread1,iptty,iout,ierr, 
     &                             num_times, num_models)


            close(iread1)
        else 
        endif



c        return

        end subroutine userr


!************************************************************ 
!   MODIFID VERSION OF model_setup.f SUBROUTINE 
!************************************************************ 
contains
      subroutine userr_model_setup(iread1, iptty, iout, ierr, 
     &      num_times, num_models) 

      use comuserr
      implicit none


      integer i,j,iread1,iptty,iout,ierr
      integer num_times, num_models
      character*80 key
c      ! Hard-code in a "time()" array for testing
      ! real*8 time(5,2)   !jpo hard-code for now
c      ! Hard-code in a "node_ch_model()" array for testing
c      ! ! real*8 node_ch_model(4,*) 
c      ! real*8 node_ch_model(4) 
c      ! Hard-code in a "time_cycle)" array for testing
c      ! real*8 time_cycle(4,*) 
c      real*8 time_cycle(4) 
c      logical boun_out
      integer irmod,irsubmod,mrmodel
      integer maxrmodel
c      integer, allocatable :: time_interpolate(:) !jpo
c      ! allocatable :: time_interpolate(:) !jpo
      integer time_type(num_models)  !jpo
c      ! ---- For the zone readin part ------
      integer iarray,ipoint,inumber,ja,jb,jc
      logical null1,readflag
      integer inode,max_arrays
      integer ii,izunit,nin, monum
      integer zarray(3)  ! array for individual ja, jb, jc specs


cc     jpo - DEBUG (iatty = 6)
c      write(6,*) 'userr_model_setup'
c      write(6,*) 'num_times  = ', num_times  
c      write(6,*) 'num_models = ', num_models 
      if(.not.allocated(distcoeffs)) then
          allocate(distcoeffs(num_times,num_models), 
     &             userrtime(num_times,num_models),
     &             distcoeff_zones(3,num_models))
      endif

      irmod=0
      irsubmod=0
      maxrmodel=1e4  ! max number of rxn models allowed
      mrmodel=num_models ! jpo - not sure if this is right..
     
      !------------------------------------------------------- 
      !----- MODEL NUMBERS -----------------------------------
      !------------------------------------------------------- 
 10   continue
      read(iread1,'(a80)') key
c      print *, key !jpo
      if(key(1:3).eq.'mod') then
         irmod=irmod+1  !increment the model number
         irsubmod=1
c         ntimes=1  !jpo - try turning this one off
         if(irmod.gt.maxrmodel) go to 20
         read(iread1,'(a9)') key
c         write(iout,*) key     !jpo - debugging
      else
         irsubmod=irsubmod+1
      endif

c      !------------------------------------------------------- 
c      !----- TIME STEPS FOR BC CHANGES -----------------------
c      !   NOTE: may want to move this outside of subroutine since
c      !         I only want to allow one set of times for all
c      !         of the models (1 model per zone)
c      !------------------------------------------------------- 
      if(key(1:2).eq.'ti') then
         read(iread1,*) (userrtime(i,irmod),i=1,num_times)
c         time_cycle(imod)=time(ntimes,imod)
         time_type(irmod)=-num_times
c       jpo - Eventually may need to add this stuff back
c         if(key(1:9).eq.'ti_linear') then
c            time_interpolate(irmod) = 1
c            if (iout .ne. 0 .and. boun_out) then
c               write(iout,*)
c     &              'Linear interpolation between time intervals'
c               write(iout,*) 'Set time_interpolate = 1'
c               write(iout,*) irmod, time_interpolate(irmod)
c            end if
cc            ! if (iptty .ne. 0 .and. boun_out) then
cc               ! write(iptty,*) &
cc     ! &              'Linear interpolation between time intervals'
cc               ! write(iptty,*) 'Set time_interpolate = 1'
cc               ! write(iptty,*) imod, time_interpolate(imod)
cc            ! end if   
c         else          
c            time_interpolate(imod) = 0
c            if (iout .ne. 0 .and. boun_out) then
c               write(iout,*) 'Constant between time intervals'
c               write(iout,*) 'Set time_interpolate = 0'
c               write(iout,*) imod, time_interpolate(imod)
c            end if
cc            ! if (iptty .ne. 0 .and. boun_out) then
cc               ! write(iptty,*) 'Constant between time intervals'
cc               ! write(iptty,*) 'Set time_interpolate = 0'
cc               ! write(iptty,*) imod, time_interpolate(imod)
cc            ! end if
c         endif

c      !------------------------------------------------------- 
c      !----- BC TYPE  OR BC CHANGES --------------------------
c      !------------------------------------------------------- 
c
c     !-------------------------------------------------- 
c     ! NEW BC CHANGE FOR KEQ VALUES    jpo
c     !    - fills array distcoeffs(ntimes, num_models)
c     !-------------------------------------------------- 
c     ! else if(key(1:9).eq.'distcoeff' .or. key(1:9).eq.'DISTCOEFF') then 
      else if(key(1:9).eq.'distcoeff') then
         read(iread1,*) (distcoeffs(i,irmod),i=1,num_times)
c         ! ! jpo - below do loop may not be necessary
c         ! do i=1,ntimes
c             ! if(distcoeffs(i,imod).eq.0.0) then
c                ! distcoeffs(i,imod)=-9999999.0 
c             ! endif
c         ! enddo
c         distcoeffs_type(irmod)=1
c         ! if(isubmod.le.1) go to 30  ! jpo - this may not be necessary




      else if(key(1:3).eq.'end'.or.key(1:3).eq.'   ') then
         mrmodel=irmod
         go to 40 
      else 
         if(iptty.ne.0) write(iptty, 100) trim(key)
         if(iout.ne.0) write(iout, 100) trim(key) 
         write(ierr, 100) trim(key)
         stop
      endif
 100  format ('illegal keyword in macro boun, stopping:', /, a)

      go to 10
 20   continue
      if(iptty.ne.0) write(iptty, 200) 
      if(iout.ne.0) write(iout, 200)
      write(ierr, 200)
      stop
 200  format ('exceeded storage for number of models, stopping')
 30   continue
      if(iptty.ne.0) write(iptty, 300)
      if(iout.ne.0) write(iout, 300 )
      write(ierr, 300)
      stop
 300  format ('time change was not first keyword, stopping')
 40   continue
c     adjust unit of time
c      do imod = 1, mmodel
c	 ! if(tunit_type(imod).eq.1) then
c! ! c     sec to days
c            ! 
c            ! do j = 1, iabs(time_type(imod))
c               ! time(j,imod) = fac_sec_days*time(j,imod)
c               ! if(timestep_type(imod).ne.0) then
c                  ! timestep(j,imod) =fac_sec_days*timestep(j,imod)
c               ! endif 
c            ! enddo
c	 ! else if(tunit_type(imod).eq.2) then
c! ! c     min to days
c            ! do j = 1, iabs(time_type(imod))
c               ! time(j,imod) = fac_min_days*time(j,imod)
c               ! if(timestep_type(imod).ne.0) then
c                  ! timestep(j,imod) =fac_min_days*timestep(j,imod)
c               ! endif 
c            ! enddo
c	 ! else if(tunit_type(imod).eq.4) then
c! ! c     years to days
c            ! do j = 1, iabs(time_type(imod))
c               ! time(j,imod) = fac_year_days*time(j,imod)
c               ! if(timestep_type(imod).ne.0) then
c                  ! timestep(j,imod) =fac_year_days*timestep(j,imod)
c               ! endif 
c            ! enddo
c      enddo

c      !------------------------------------------------------- 
c      !----- ZONES FOR DISTCOEFF CHANGES ---------------------
c      !------------------------------------------------------- 
c      jpo - eventually make this more general so zone numbers
c            can be read in (not just nodes as JA JB JC)

c      !--------------------------
c      !  JPO - from initdata2.f
c      !--------------------------
c      ! Should be at the first blank line before zones
c      ! Read num_models number of lines after that first blank line
c      print *
c      print *, 'AFTER DISTCOEFFS HAVE BEEN READ'
c      print *, 'Looping through ',num_models, 'models...'
c          print *, 'NODES (JA JB JC) FOR EACH MODEL:'
      do i = 1,num_models 
          read(iread1, *) ja, jb, jc, monum
c          print *, '(model',monum,')', ja, jb, jc
          zarray = (/ja,jb,jc/)
          do j = 1,3
              distcoeff_zones(j,i) = zarray(j)
          enddo
      enddo

c     Loop over all lines of data to read
c! 1000 continue
c         ! read(in_number, '(a80)') strtot
c         ! if( null1(strtot) .eqv. .TRUE. ) goto 2000
c         ! backspace in_number
c         ! read(in_number, *, ERR = 3000) ja, jb, jc, 
c     ! 2        (values(iarray), iarray = 1, narrays )
c         ! if (ja .eq. 0) goto 2000
c         ! goto 4000
c! 3000    continue
c         ! inumber = ireturn + 1
c         ! ireturn = -2
c         ! goto 9000
c ! 4000    continue
c         ! ireturn = ireturn + 1
c! c     Input by zones is first, or else input is by node
c         ! if( ja .lt. 0 ) then
    

      return  !return from userr_model_setup
      end subroutine userr_model_setup

