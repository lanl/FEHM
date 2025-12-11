       subroutine model_setup(inpt, iptty, iout, ierr)
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
CD1 PURPOSE
CD1 
CD1 To read in parameters for time-dependent boundary models.  
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2
CD2 Initial implementation: NOV-96, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/model_setup.f_a  $
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:08   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:10   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C***********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 2.3.7 Sources and sinks
CD3 2.6   Provide Input/Output Data Files
CD3 3.0   INPUT AND OUTPUT REQUIREMENTS
CD3 
C***********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4 
C***********************************************************************
c       
c
c read in model parameters
c
c
      use comai, only : boun_out
      use comdti
      use comdi

      implicit none

      integer i,j,ntimes,inpt,iptty,iout,ierr
      real*8 tol_boun
      parameter(tol_boun=1.d-30)
      character*80 key
c     
c     imod=0
      imod=mmodel 
      isubmod=0
c     
 10   continue
      read(inpt,'(a80)') key
      if(key(1:3).eq.'mod') then
         imod=imod+1
         isubmod=1
         ntimes=1
         if(imod.gt.maxmodel) go to 20
         read(inpt,'(a9)') key
      else
         isubmod=isubmod+1
      endif
      if(key(1:2).eq.'ti') then
         read(inpt,*) ntimes,(time(i,imod),i=1,ntimes)
         time_cycle(imod)=time(ntimes,imod)
         time_type(imod)=-ntimes
         if(key(1:9).eq.'ti_linear') then
            time_interpolate(imod) = 1
            if (iout .ne. 0 .and. boun_out) then
               write(iout,*)
     &              'Linear interpolation between time intervals'
               write(iout,*) 'Set time_interpolate = 1'
               write(iout,*) imod, time_interpolate(imod)
            end if
            if (iptty .ne. 0 .and. boun_out) then
               write(iptty,*)
     &              'Linear interpolation between time intervals'
               write(iptty,*) 'Set time_interpolate = 1'
               write(iptty,*) imod, time_interpolate(imod)
            end if
         else
            time_interpolate(imod) = 0
            if (iout .ne. 0 .and. boun_out) then
               write(iout,*) 'Constant between time intervals'
               write(iout,*) 'Set time_interpolate = 0'
               write(iout,*) imod, time_interpolate(imod)
            end if
            if (iptty .ne. 0 .and. boun_out) then
               write(iptty,*) 'Constant between time intervals'
               write(iptty,*) 'Set time_interpolate = 0'
               write(iptty,*) imod, time_interpolate(imod)
            end if
         endif
      else if(key(1:4).eq.'tran') then
c gaz 062620 turn off tran but leave in
         if(key(6:8).ne.'off') then
          steady_type(imod)=-1
         endif
      else if(key(1:3).eq.'sec') then
         tunit_type(imod)=1
      else if(key(1:3).eq.'min') then
         tunit_type(imod)=2
      else if(key(1:3).eq.'day') then
         tunit_type(imod)=3
      else if(key(1:4).eq.'year') then
         tunit_type(imod)=4
      else if(key(1:2).eq.'cy') then
         read(inpt,*) ntimes,(time(i,imod),i=1,ntimes)
         time_cycle(imod)=time(ntimes,imod)
         time_type(imod)=ntimes
         if(key(1:9).eq.'cy_linear') then
            time_interpolate(imod) = 1
            if (iout .ne. 0 .and. boun_out) then
               write(iout,*) 
     &              'Linear interpolation between time intervals'
               write(iout,*) 'Set time_interpolate = 1'
               write(iout,*) imod, time_interpolate(imod)
            end if
            if (iptty .ne. 0 .and. boun_out) then
               write(iptty,*) 
     &              'Linear interpolation between time intervals'
               write(iptty,*) 'Set time_interpolate = 1'
               write(iptty,*) imod, time_interpolate(imod)
            end if
         else
            time_interpolate(imod) = 0
            if (iout .ne. 0 .and. boun_out) then
               write(iout,*) 'Constant between time intervals'
               write(iout,*) 'Set time_interpolate = 0'
               write(iout,*) imod, time_interpolate(imod)
            end if
            if (iptty .ne. 0 .and. boun_out) then
               write(iptty,*) 'Constant between time intervals'
               write(iptty,*) 'Set time_interpolate = 0'
               write(iptty,*) imod, time_interpolate(imod)
            end if
         endif
      else if(key(1:4).eq.'chmo') then
         read(inpt,*) (node_ch_model(i,imod),i=1,ntimes)
         node_ch_model_type(imod)=1
         if(isubmod.le.1) go to 30
c     new initial value input here 12/3/98 GAZ
      else if(key(1:3).eq.'tmi') then
         read(inpt,*) (temperature_ini(i,imod),i=1,ntimes)
         temperature_ini_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'si') then
         read(inpt,*) (saturation_ini(i,imod),i=1,ntimes)
         saturation_ini_type(imod)=1
         if(isubmod.le.1) go to 30 
      else if(key(1:3).eq.'pai') then
         read(inpt,*) (pressurea_ini(i,imod),i=1,ntimes)
         pressurea_ini_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'pwi') then
         read(inpt,*) (pressurew_ini(i,imod),i=1,ntimes)
         pressurew_ini_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'sfh') then
         read(inpt,*) (seepfac(i,imod),i=1,ntimes)
         seepfac_type(imod)=1
         do i=1,ntimes
            seepfac(i,imod)= -(seepfac(i,imod)+tol_boun)
         enddo
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'sf') then
         read(inpt,*) (seepfac(i,imod),i=1,ntimes)
         seepfac_type(imod)=1
         do i=1,ntimes
            seepfac(i,imod)= seepfac(i,imod)+tol_boun
         enddo
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'fd') then
         read(inpt,*) (drainar(i,imod),i=1,ntimes)
         drainar_type(imod)=1
         do i=1,ntimes
            drainar(i,imod)= drainar(i,imod)+tol_boun
         enddo 
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'dfd') then
         read(inpt,*) (drainar(i,imod),i=1,ntimes)
         drainar_type(imod)=-1
         do i=1,ntimes
            drainar(i,imod)= drainar(i,imod)+tol_boun
         enddo
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'fxa') then
c fxa is air mass fraction of inflow
c gaz I think I can use array sourcea because air flowrate
c and air fraction cannot both be set for the same node
c gaz 111418 can exist in different model        
         read(inpt,*) (sourcea(i,imod),i=1,ntimes)
         sourcea_type(imod)=-2
         do i=1,ntimes
            sourcea(i,imod)= sourcea(i,imod)+tol_boun
         enddo
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'sa') then
         read(inpt,*) (sourcea(i,imod),i=1,ntimes)
         sourcea_type(imod)=1
         do i=1,ntimes
            sourcea(i,imod)= sourcea(i,imod)+tol_boun
         enddo
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'swf') then
         do i=1,n0
            qw0(i) = sk(i)
         enddo
         read(inpt,*) (sourcef(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcef(i,imod)= sourcef(i,imod)+tol_boun
         enddo
         sourcef_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'sw') then
         read(inpt,*) (sourcew(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcew(i,imod)= sourcew(i,imod)+tol_boun
         enddo
         sourcew_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'se') then
         read(inpt,*) (sourcee(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcee(i,imod)= sourcee(i,imod)+tol_boun
         enddo
         sourcee_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:4).eq.'sco2') then
         read(inpt,*) (sourceco2(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourceco2(i,imod)= sourceco2(i,imod)+tol_boun
         enddo
         sourceco2_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'dsa') then
         read(inpt,*) (sourcea(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcea(i,imod)= sourcea(i,imod)+tol_boun
         enddo
         sourcea_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'dsw') then
         read(inpt,*) (sourcew(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcew(i,imod)= sourcew(i,imod)+tol_boun
         enddo
         sourcew_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'dse') then
         read(inpt,*) (sourcee(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourcee(i,imod)= sourcee(i,imod)+tol_boun
         enddo
         sourcee_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:5).eq.'dsco2') then
         read(inpt,*) (sourceco2(i,imod),i=1,ntimes)
         do i=1,ntimes
            sourceco2(i,imod)= sourceco2(i,imod)+tol_boun
         enddo
         sourceco2_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:4).eq.'tco2') then
         read(inpt,*) (esourceco2(i,imod),i=1,ntimes)
         do i=1,ntimes
            esourceco2(i,imod)= esourceco2(i,imod)+tol_boun
         enddo
         esourceco2_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'s ') then
         read(inpt,*) (saturation(i,imod),i=1,ntimes)
         saturation_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'pwo') then
         read(inpt,*) (pressurew(i,imod),i=1,ntimes)
         do i=1,ntimes
            if(pressurew(i,imod).eq.0.0) then
               pressurew(i,imod)=-9999999.0
            endif
         enddo
         pressurew_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'pao') then
         read(inpt,*) (pressurea(i,imod),i=1,ntimes)
         do i=1,ntimes
            if(pressurea(i,imod).eq.0.0) then
               pressurea(i,imod)=-9999999.0
            endif
         enddo
         pressurea_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'pw') then
         read(inpt,*) (pressurew(i,imod),i=1,ntimes)
         do i=1,ntimes
            if(pressurew(i,imod).eq.0.0) then
               pressurew(i,imod)=-9999999.0
            endif
         enddo
         pressurew_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'hdo') then
         read(inpt,*) (temperature(i,imod),i=1,ntimes)
         do i=1,ntimes
            if(temperature(i,imod).eq.0.0) then
               temperature(i,imod)=-9999999.0
            endif
         enddo
         temperature_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'hd') then
         read(inpt,*) (temperature(i,imod),i=1,ntimes)
         do i=1,ntimes
            if(temperature(i,imod).eq.0.0) then
               temperature(i,imod)=-9999999.0
            endif
         enddo
         temperature_type(imod)=1
         if(isubmod.le.1) go to 30
c gaz new 061816
      else if(key(1:3).eq.'huf') then
         read(inpt,*) (humid(i,imod),i=1,ntimes)
         humid_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'hu') then
         read(inpt,*) (humid(i,imod),i=1,ntimes)
         humid_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'ph') then
         read(inpt,*) (phumid(i,imod),i=1,ntimes)
         phumid_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'th') then
         read(inpt,*) (thumid(i,imod),i=1,ntimes)
         thumid_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:1).eq.'h') then
         read(inpt,*) (saturation(i,imod),i=1,ntimes)
         do i=1,ntimes
            saturation(i,imod)= saturation(i,imod)+tol_boun
         enddo
         saturation_type(imod)=-1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'pa') then
         read(inpt,*) (pressurea(i,imod),i=1,ntimes)
         do i=1,ntimes
            pressurea(i,imod)= pressurea(i,imod)+tol_boun
         enddo
         pressurea_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'if') then
         read(inpt,*) (impedance(i,imod),i=1,ntimes)
         impedance_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:3).eq.'fen') then
         read(inpt,*) (enthalpy(i,imod),i=1,ntimes)
         enthalpy_type(imod)=1
         if(isubmod.le.1) go to 30         
      else if(key(1:2).eq.'en' .and. key(1:3).ne.'end') then
         read(inpt,*) (enthalpy(i,imod),i=1,ntimes)
         enthalpy_type(imod)=1
         if(isubmod.le.1) go to 30
c gaz 012122 added equation tolerance    
      else if(key(1:4).eq.'eqto') then
         read(inpt,*) (eqtolerance(i,imod),i=1,ntimes)
         eqtol_type(imod)=1
         if(isubmod.le.1) go to 30
c gaz 012122 added max timestep      
      else if(key(1:4).eq.'tsma') then
         read(inpt,*) (timestepmax(i,imod),i=1,ntimes)
         tsmax_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'ts') then
         read(inpt,*) (timestep(i,imod),i=1,ntimes)
         timestep_type(imod)=1
         if(isubmod.le.1) go to 30         
      else if(key(1:1).eq.'t') then
         read(inpt,*) (temperature(i,imod),i=1,ntimes)
         do i=1,ntimes
            temperature(i,imod) = temperature(i,imod) + tol_boun
         enddo
         temperature_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2).eq.'ft') then
         read(inpt,*) (enthalpy(i,imod),i=1,ntimes)
         do i=1,ntimes
            enthalpy(i,imod)=-enthalpy(i,imod)
         enddo
         enthalpy_type(imod)=1
         if(isubmod.le.1) go to 30
! Add permeability change model (from FEHM V2.26)
      else if(key(1:2) .eq. 'kx') then
         read(inpt,*) (permx(i,imod),i=1,ntimes)
         permx_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2) .eq. 'ky') then
         read(inpt,*) (permy(i,imod),i=1,ntimes)
         permy_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:2) .eq. 'kz') then
         read(inpt,*) (permz(i,imod),i=1,ntimes)
         permz_type(imod)=1
         if(isubmod.le.1) go to 30
      else if(key(1:len_trim(key)).eq.'wgtp') then
         weight_type(imod)=2
      else if(key(1:len_trim(key)).eq.'wgtr') then
         weight_type(imod)=4
      else if (key(1:len_trim(key)).eq.'wgww') then
         weight_type(imod)=5
      else if (key(1:len_trim(key)).eq.'wgtx') then
         weight_type(imod)=6
      else if (key(1:len_trim(key)).eq.'wgty') then
         weight_type(imod)=7
      else if (key(1:len_trim(key)).eq.'wgtz') then
         weight_type(imod)=8
      else if (key(1:len_trim(key)).eq.'wgtpx') then
         weight_type(imod)=9
      else if (key(1:len_trim(key)).eq.'wgtpy') then
         weight_type(imod)=10
      else if (key(1:len_trim(key)).eq.'wgtpz') then
         weight_type(imod)=11
      else if(key(1:len_trim(key)).eq.'wgtu') then
         weight_type(imod)=3
      else if(key(1:len_trim(key)).eq.'wgt') then
         weight_type(imod)=1
      else if(key(1:3).eq.'end'.or.key(1:3).eq.'   ') then
         mmodel=imod
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
      do imod = 1, mmodel
	 if(tunit_type(imod).eq.1) then
c     sec to days
            
            do j = 1, iabs(time_type(imod))
               time(j,imod) = fac_sec_days*time(j,imod)
               if(timestep_type(imod).ne.0) then
                  timestep(j,imod) =fac_sec_days*timestep(j,imod)
               endif 
            enddo
	 else if(tunit_type(imod).eq.2) then
c     min to days
            do j = 1, iabs(time_type(imod))
               time(j,imod) = fac_min_days*time(j,imod)
               if(timestep_type(imod).ne.0) then
                  timestep(j,imod) =fac_min_days*timestep(j,imod)
               endif 
            enddo
	 else if(tunit_type(imod).eq.4) then
c     years to days
            do j = 1, iabs(time_type(imod))
               time(j,imod) = fac_year_days*time(j,imod)
               if(timestep_type(imod).ne.0) then
                  timestep(j,imod) =fac_year_days*timestep(j,imod)
               endif 
            enddo
         endif
      enddo
      return
      end                          
