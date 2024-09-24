       subroutine flow_boun(iz,n,ico2,idof
     * ,days0,days,day,daynew,vol_nd
     * ,inpt,iptty,iout,ierr,l,igrav,ihead)
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
CD1 To read and apply boundary and source conditions for FEHM.   
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 8-21-94, Programmer: G. Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/flow_boun.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:04:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2                                     
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 iz           int     I       Control parameter
CD3 n            int     I       Total number of unknowns
CD3 ico2         int     I       physics model identifier             
CD3 idof         int     I       degree of freedom identifier           
CD3 node_model   int     I       array of boundary model identifers
CD3 modmin       int     O       boundary model that caused
CD3                              time step change
CD3 days0        real*8  I       time(days) at previous time step
CD3 days         real*8  I       time(days) at current time step
CD3 day          real*8  I       current timestep size          
CD3 vol_nd       real*8  I       array of nodal volumes
CD3 qa           real*8  O       array of nodal air source terms
CD3 qw           real*8  O       array of nodal water source terms
CD3 qh           real*8  O       array of nodal enthalpy source terms
CD3 sat          real*8  O       array of nodal saturation boundary terms
CD3 enth         real*8  O       array of nodal saturation boundary terms
CD3 temp         real*8  O       array of nodal temperture boundary terms
CD3 presa        real*8  O       array of nodal air pressure boundary terms
CD3 presw        real*8  O       array of nodal water pressure boundary terms
CD3 imped        real*8  O       array of nodal impedance terms
CD3 inpt         int     I       unit number for input file
CD3 iptty       int     I       unit number for ouput messages
CD3 iout        int     I       unit number for ouput messages
CD3 l           int     I       current time step number        
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 inpt                   I    File contains input data
CD3 iptty                  I    File contains output data
CD3 iout                   I    File contains output data
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4
CD4 NONE
CD4 
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier     Type        Description
CD5 
CD5 dmaxmodel       int         maximum number of boundary models
CD5 dmaxtimes       int         maximum number of boundary time changes
CD5 
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier      Type        Description
CD5 
CD5 sourcea_type    int         identifer of air flux term being used  
CD5 sourcew_type    int         identifer of water flux term being used 
CD5 sourcef_type    int         identifer of water flux factor used 
CD5 sourcee_type    int         identifer of enthalpy flux term being
CD5                             used    
CD5 saturation_type int         identifer of saturation boundary being
CD5                             used    
CD5 pressurea_type  int         identifer of air pressure boundary      
CD5                             being used    
CD5 pressurew_type  int         identifer of water pressure boundary  
CD5                             being used    
CD5 enthalpy_type   int         identifer of enthalpy boundary      
CD5                             being used    
CD5 temperture_type int         identifer of temperture boundary      
CD5                             being used    
CD5 timestep_type   int         identifer of time step sze change     
CD5                             being used    
CD5 time_type       int         position in time change array
CD5 time_interpolate int        Interpolation method:
CD5                             = 0 constant value within time interval
CD5                             = 1 linear interpolateion between time intervals
CD5 impedance_type  int         identifer of impedance factor           
CD5                             being used    
CD5 
CD5 Local Subprograms
CD5
CD5 model_setup
CD5 time_adjust
CD5
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.3.7 Sources and sinks
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN flow_boun
CPS 
CPS IF data are read in on this call
CPS   
CPS   read keyword
CPS   if keyword is not null 
CPS    CALL model_setup
CPS
CPS ENDIF
CPS   
CPS   ENDLOOP
CPS   
CPS   Set parameters for reading region numbers
CPS     
CPS   initdata - read region numbers and set at each node
CPS   
CPS   Set flag to denote that the rlp macro has been called
CPS       
CPS 
CPS   FOR each node
CPS   
CPS END flow_boundary_conditions
CPS
C**********************************************************************
c
c
c  routine to set boundary conditions for all modules
c
c
c key words glossary
c
c arrays for source and variable specification
c
c sourcea - source term for air flow (kg/s)
c sourcew - source term for water flow (kg/s)
c sourcee - source term for enthalpy flow (kg/s)
c saturation -specified saturation for source node (mass fraction)
c impedance - specified impedance          
c pressurea - specified air pressure (Mpa)
c pressurew - specified water pressure (Mpa)
c enthalpy - specified enthalpy (Mj/kg)
c timestep - specified timestep size(days )
c node_ch_model - specified model number change 
c
c arrays for identifying model attribrutes
c
c sourcea_type - source term for air flow (kg/s)
c sourcew_type - source term for water flow (kg/s)
c sourcee_type - source term for enthapy flow (kg/s)
c saturation_type -specified saturation for source node (mass fraction)
c impedance_type - specified impedance
c pressurea_type - specified air pressure (Mpa)
c pressurew_type - specified water pressure (Mpa)
c enthalpy_type - specified enthalpy (Mj/kg)
c timestep_type - specified timestep size(days )
c temperature_type - specified temperature(deg C)
c time_type - change any of above with time
c node_ch_model_type - specified model number change 

c
c node_model - model number for each node
c
      use comai, only : boun_out, neq_primary, daymaxboun, daymax,
     &      tmch,tmchboun
      use combi
      use comci
      use comdi
      use comdti, only : zero_t
      use comsplitts, only : dtot_next
      use comwt

      implicit none

      integer ico2,iz,i,j,n,idof,ntimes,imodel      
      integer inpt,iptty,iout,ierr,l,ihead
      integer iimodel,istea,pump_node(1000),zmax(1000)
      integer igrav
c gaz 090820
      integer mod1,mod2
      
      real*8 day,days0,days,daynew,daym1,vfac,tmdum
      real*8 time_factor1, time_factor2, if_time_interpolate
c     
      real*8 vol_nd(*)
      real*8 perm_tot,dd,rldum,sf, v11norm
c     
      logical checkpa,checkpw,checke,checks,checkt,checkqa,checkqw
      logical time_interpolate_used
      character key_word*4,timchar*9

      zmax=-1000
      sf = 1.0d0
      v11norm = 0.0d0
c     
      if(iz.eq.0) then
c     
c     read model        
c     
         read(inpt,'(a4)') key_word
         if(key_word.ne.'    ') then         
            backspace inpt
c     
c     increment model number
c     
            call model_setup(inpt,iptty,iout,ierr)
         else
            return
         endif
c_______________________________________________________________________
      else if(iz.eq.1) then           
c     
c     adjust times so ending time is very long
c     
         do imod=1,mmodel
            ntimes=abs(time_type(imod))
c     gaz 9-18-00
            ntimes=ntimes+1
            if(time(1,imod).ne.0.0) then
               do i=ntimes,2,-1
                  time(i,imod)=time(i-1,imod)
                  if(iqa.ne.0) sourcea(i,imod)=sourcea(i-1,imod)
                  if(iqw.ne.0) sourcew(i,imod)=sourcew(i-1,imod)
                  if(iqco2.ne.0) sourceco2(i,imod)=sourceco2(i-1,imod)
                  if(iqf.ne.0) sourcef(i,imod)=sourcef(i-1,imod)
                  if(isf.ne.0) seepfac(i,imod)=seepfac(i-1,imod)
                  if(ifd.ne.0) drainar(i,imod)=drainar(i-1,imod)
                  if(iqenth.ne.0) sourcee(i,imod)=sourcee(i-1,imod)
                  if(isatb.ne.0)
     &                 saturation(i,imod)=saturation(i-1,imod)
                  if(isatb_ini.ne.0) 
     &                 saturation_ini(i,imod)=saturation(i-1,imod)
                  if(ipresw.ne.0) pressurew(i,imod)=pressurew(i-1,imod)
                  if(ipresw_ini.ne.0) 
     &                 pressurew_ini(i,imod)=pressurew(i-1,imod)
                  if(ipresa.ne.0) pressurea(i,imod)=pressurea(i-1,imod)
                  if(ipresa_ini.ne.0) 
     &                 pressurea_ini(i,imod)=pressurea(i-1,imod)
                  if(imped.ne.0) impedance(i,imod)=impedance(i-1,imod)
                  if(ienth.ne.0) enthalpy(i,imod)=enthalpy(i-1,imod)
                  timestep(i,imod)=timestep(i-1,imod)
                  if(itempb.ne.0) 
     &                 temperature(i,imod)=temperature(i-1,imod)
                  if(itempb_ini.ne.0) 
     &                 temperature_ini(i,imod)=temperature(i-1,imod)
                  if(ixperm.ne.0) permx(i,imod)= permx(i-1,imod)
                  if(iyperm.ne.0) permy(i,imod)= permy(i-1,imod)
                  if(izperm.ne.0) permz(i,imod)= permz(i-1,imod)
                  if(icm.ne.0) 
     &                 node_ch_model(i,imod)=node_ch_model(i-1,imod)
                  if(iha.ne.0) humid(i,imod)=humid(i-1,imod)
                  if(ipha.ne.0) phumid(i,imod)=phumid(i-1,imod)
                  if(itha.ne.0) thumid(i,imod)=thumid(i-1,imod)
               enddo
               ntimes=ntimes+1
               if(iqa.ne.0) sourcea(1,imod)=0.0
               if(iqw.ne.0) sourcew(1,imod)=0.0
               if(iqco2.ne.0) sourceco2(1,imod)=0.0
               if(isf.ne.0) seepfac(1,imod)=0.0
               if(ifd.ne.0) drainar(1,imod)=0.0
               if(iqf.ne.0) sourcef(1,imod)=0.0
               if(iqenth.ne.0) sourcee(1,imod)=0.0
               if(isatb.ne.0) saturation(1,imod)=0.0
               if(ipresw.ne.0) pressurew(1,imod)=0.0
               if(ipresa.ne.0) pressurea(1,imod)=0.0
               if(imped.ne.0) impedance(1,imod)=0.0
               if(ienth.ne.0) enthalpy(1,imod)=0.0
               timestep(1,imod)=day
               if(itempb.ne.0) temperature(1,imod)=0.0
               if(ixperm.ne.0) permx(1,imod)= 0.
               if(iyperm.ne.0) permy(1,imod)= 0.
               if(izperm.ne.0) permz(1,imod)= 0.
               if(icm.ne.0) node_ch_model(i,imod)=0   
                  if(iha.ne.0) humid(i,imod)= 0.0
                  if(ipha.ne.0) phumid(i,imod)= 0.0
                  if(itha.ne.0) thumid(i,imod)= 0.0
            endif
            if(time_type(imod).gt.0) then
               if(iqa.ne.0) sourcea(ntimes,imod)=sourcea(1,imod)
               if(iqw.ne.0) sourcew(ntimes,imod)=sourcew(1,imod)
               if(iqco2.ne.0) sourceco2(ntimes,imod)=sourceco2(1,imod)
               if(isf.ne.0) seepfac(ntimes,imod)=seepfac(1,imod)
               if(ifd.ne.0) drainar(ntimes,imod)=drainar(1,imod)
               if(iqf.ne.0) sourcef(ntimes,imod)=sourcef(1,imod)
               if(iqenth.ne.0) sourcee(ntimes,imod)=sourcee(1,imod)
               if(isatb.ne.0) 
     &              saturation(ntimes,imod)=saturation(1,imod)
               if(isatb_ini.ne.0) 
     &              saturation_ini(ntimes,imod)=saturation(1,imod)
               if(ipresw.ne.0) pressurew(ntimes,imod)=pressurew(1,imod)
               if(ipresw_ini.ne.0) 
     &              pressurew_ini(ntimes,imod)=pressurew(1,imod)
               if(ipresa.ne.0) pressurea(ntimes,imod)=pressurea(1,imod)
               if(ipresa_ini.ne.0) 
     &              pressurea_ini(ntimes,imod)=pressurea(1,imod)
               if(imped.ne.0) impedance(ntimes,imod)=impedance(1,imod)
               if(ienth.ne.0) enthalpy(ntimes,imod)=enthalpy(1,imod)
               timestep(ntimes,imod)=timestep(1,imod)
               if(itempb.ne.0) 
     &              temperature(ntimes,imod)=temperature(1,imod)
               if(itempb_ini.ne.0) 
     &              temperature_ini(ntimes,imod)=temperature(1,imod)
               if(ixperm.ne.0) permx(ntimes,imod)=permx(1,imod) 
               if(iyperm.ne.0) permy(ntimes,imod)=permy(1,imod)
               if(izperm.ne.0) permz(ntimes,imod)=permz(1,imod) 
               if(icm.ne.0) 
     &              node_ch_model(ntimes,imod)=node_ch_model(1,imod)
                  if(iha.ne.0) humid(ntimes,imod)=humid(1,imod)
                  if(ipha.ne.0) phumid(ntimes,imod)=phumid(1,imod)
                  if(itha.ne.0) thumid(ntimes,imod)=thumid(1,imod)
            endif
            if(timestep(1,imod).le.0.0) then
               timestep(1,imod)=day
            endif
            if(time_type(imod).lt.0) then
c     time(ntimes+1,imod)=1.e30
c     GAZ 2-23-00
               time(ntimes,imod)=1.e30
               time_cycle(imod)=1.e30 
            endif
            if(time_type(imod).lt.0) then
               time_type(imod)=-1
c     time(1,imod)=0.0
c     GAZ 2-23-00
c     time(1,imod)=days
               time(1,imod)=0.0
            else
               time_type(imod)=1
               time(1,imod)=time_cycle(imod)
            endif
         enddo
c     
c     initialize lchange
c     
         lchange=0
c     
c     first set check parameters false for compatibility check
c     
         checke=.false.
         checks=.false.
         checkt=.false.
         checkqa=.false.
         checkqw=.false.
         checkpa=.false.
         checkpw=.false.
c     
c     check for consistency
c     
         if(idof.eq.1) then
c     
c     heat only conduction solution
c     
c     
c     check for consistency
            
            if(itempb.ne.0) then
               checkt=.true.      
            endif
            if(checkt) then
            else
               if (iout .ne. 0) write(iout, 100)
               if (iptty .ne. 0) write(iptty, 100)
 100           format(1x,'>> Warning: boundary conditions incompatible'
     &              ,' with heat conduction')
c     stop        
            endif      
         else if(ico2.lt.0 ) then
c     
c     air water isothermal solution                     
c     
c     check for consistency
c     
            if(isatb.ne.0) then
               checks=.true.
            endif
            if(ipresa.ne.0) then
               checkpa=.true.
            endif
            if(ipresw.ne.0) then
               checkpw=.true.
            endif
c     head info is contained in temperature arrays for this physics
c     module
            if(itempb.ne.0) then
               checkt=.true.      
            endif
            if(checkpa.or.checks) then
            else if(checkpw.or.checks) then
            else if(checkt.or.checks) then
            else
               if (iout .ne. 0) write(iout, 101)
               if (iptty .ne. 0) write(iptty, 101)
 101           format(1x, '>> Warning: boundary conditions ',
     &              'incompatible with iso air/water <<')
c     stop        
            endif
         else if(ico2.eq.0 ) then
c     
c     water(liquid and vapor) and heat
c     
c     check for consistency
c     
            if(isatb.ne.0) then
               checks=.true.
            endif
            if(ipresw.ne.0) then
               checkpw=.true.
            endif
            if(ienth.ne.0) then
               checke=.true.
            endif
            if(itempb.ne.0) then
               checkt=.true.      
            endif
            if(checkpw.and.checks) then
            else if(checkpw.and.checke) then
            else if(checkpw.and.checkt) then
            else if(checkpa.and.checke) then
            else if(checkpw.and.checkt) then
            else if(checkpa.and.checks) then
            else
               if (iout .ne. 0) write(iout ,102)
               if (iptty .ne. 0) write(iptty, 102)
 102           format(1x, '>> Warning: boundary conditions ',
     &              'incompatible with vapor/water <<')
            endif
         else if(ico2.gt.0) then
c     
c     air and water(liquid and vapor) and heat
c     non condensible gas
c     
c     check for consistency
c     
            if(isatb.ne.0) then
               checks=.true.
            endif
            if(ipresa.ne.0) then
               checkpa=.true.
            endif
            if(ipresw.ne.0) then
               checkpw=.true.
            endif
            if(ienth.ne.0) then
               checke=.true.
            endif
            if(itempb.ne.0) then
               checkt=.true.      
            endif
            if(checkpa.and.checks.and.checke) then
            else if(checkpa.and.checks.and.checkt) then
            else if(checkpa.and.checkpw.and.checke) then
            else if(checkpa.and.checkpw.and.checkt) then
            else
               if (iout .ne. 0) write(iout ,103)
               if (iptty .ne. 0) write(iptty, 103)
 103           format(1x, '>> Warning: boundary conditions ',
     &              'incompatible with noniso air/water <<')
c     stop        
            endif
         endif
c_______________________________________________________________________
      else if (iz.eq.3) then
c     
c     apply model
c     
c     find minumum time
c     
c     before the first timestep is a special case
         if(l.gt.0.and.lchange.eq.0) then
            call time_adjust(mmodel,modmin,time_type,time_cycle,
     *           time,days0,day,days,daynew,maxtimes,iptty,iout,l,   
     *           lchange,node_model,steady_type,isty,n)
            do i=1,mmodel
C     
C     Trigger model updates in the case of
C     models with time interpolation.
C     
               if (time_interpolate(i) .ne. 0) then
                  lchange = l
               end if
            enddo
         else if(l.eq.0) then
            lchange=0
            do i=1,mmodel
               if(time_type(i).lt.0) then
c     GAZ 9-20-00
c     do j=1,maxtimes
                  do j=2,maxtimes
                     if(days.le.time(j,i)) then
                        go to 810
                     endif
                  enddo
 810              continue
                  time_type(i)=-(j-1)
               else
                  time_type(i)=1
                  days0=0.0d00
c     GAZ 9-20-00
c     do j=1,maxtimes
                  do j=2,maxtimes
                     if(days.le.time(j,i)) then
                        go to 820
                     endif
                  enddo
 820              continue
c     time_type(i)=j
                  time_type(i)=j-1
               endif
            enddo
         endif
c     
c     if change has occurred(modmin ne 0) then update model
c     
         if(lchange.eq.l) then
            lchange=0
c     
c     update all models that change at the beginning of this time step
c     
            min_model = 0
            do i=1,mmodel
               if(l.eq.0) then
c     daym1=0.0d00
c     GAZ 2-15-00
c     daym1=days    
                  daym1=0.0d00
               else if(time_type(i).lt.0) then
                  daym1=time(abs(time_type(i)-1),i)
               else
                  daym1=time(abs(time_type(i)+1),i)
               endif
               if((daym1-days0)/max(1.d-6,days0).le.1.d-6) then
                  if(time_type(i).lt.0.and.l.ne.0) then
                     time_type(i)=time_type(i)-1
                     min_model(i) = i
                  else if(time_type(i).gt.0.and.l.ne.0) then
                     time_type(i)=time_type(i)+1
                     min_model(i) = i
                  endif
                  if(abs(time_cycle(i)-days0)/max(1.d-6,days0).le.
     &                 1.d-6) then
                     do j=2,time_type(i)
                        time(j,i)=time(j,i)+time(1,i)
                     enddo
                     time_cycle(i)=time_cycle(i)+time(1,i)
                     time_type(i)=1
                  endif
               endif
            enddo
            do i=1,mmodel
               if(isty.ne.0) then
                  istea = steady_type(i)
               else
                  istea = 1
               endif
               if(istea.gt.0) then
                  if(min_model(i).ne.0) then
                     if(tunit_type(i).eq.1) then
                        tmdum = days0/fac_sec_days
                        timchar = ' seconds '
                     else if(tunit_type(i).eq.2) then
                        tmdum = days0/fac_min_days
                        timchar = ' minutes '
                     else if(tunit_type(i).eq.4) then
                        tmdum = days0/fac_year_days
                        timchar = '   years '
                     endif
	             if(tunit_type(i).ne.3) then
                        if (iout .ne. 0 .and. boun_out) then
                           write(iout ,*) ' '
                           write(iout ,10) i,days0,timchar,tmdum
                           write(iout ,*) ' '
                        end if
                        if(iptty.ne.0 .and. boun_out) then
                           write(iptty ,*) ' '
                           write(iptty,10) i,days0,timchar,tmdum
                           write(iptty ,*) ' '
                        endif
                   else
c gaz 090820 added output for boun keyword chmo(correction 092920)
                    if(icm.ne.0) then
                     mod2 = 0
                     mod1 = 0                         
                     if(node_ch_model_type(min_model(i)).ne.0) then
                      mod2=node_ch_model(abs(time_type(i)),min_model(i))
                      mod1=node_ch_model(abs(time_type(i))-1,
     &                  min_model(i))
                     endif
                    else
                     mod2 = 0
                     mod1 = 0 
                    endif
                        if (iout .ne. 0 .and. boun_out) then
                           write(iout ,*) ' '
                           write(iout ,11) i,days0
                           if(mod2.ne.0) then
                            write(iout ,12) i,mod2, mod1  
                           endif
                           write(iout ,*) ' '
                        end if
                        if(iptty.ne.0 .and. boun_out) then
                           write(iptty ,*) ' '
                           write(iptty,11) i,days0
                           if(mod2.ne.0) then
                            write(iptty,12) i,mod2, mod1  
                           endif
                           write(iptty ,*) ' '
                        endif
                     endif
 10                  format(1x,'BOUNDARY CONDITION CHANGE : model # ',
     *                    i3, ' time(days) ',g20.13, a9, g12.5)
 11                  format(1x,'BOUNDARY CONDITION CHANGE : model # ',
     *                    i3, ' time(days) ',g20.13)
 12                  format(1x,'MODEL ',i5,' : BOUN MODEL CHANGE :'
     &                      , ' model ',i5,' replaces model ',i5)
                  endif
               end if
            enddo
            if(l.ne.0) then
               do i=1,n
                  if(abs(node_model(i)).eq.modmin) then
                     node_model(i) = abs(node_model(i))
                  endif
               enddo
            endif

            do i=1,mmodel
c     gaz steady state management 10-30-04
               if(isty.ne.0) then
                  istea = steady_type(i)
               else
                  istea = 1
               endif
               if(istea.gt.0) then
                  if(iqw.ne.0) then
                     if(sourcew_type(i).lt.0) then
                        vtotw(i)=0.0d00
                     endif
                  endif
c     zvd 16-Jul08, all weighting will use vtotw for now
                  if(iqenth.ne.0) then
                     if(sourcee_type(i).lt.0) then
                        vtotw(i)=0.0d00
c     vtote(i)=0.0d00
                     endif
                  endif
c gaz 012819 multiple places distinguish between dsa and fxa by -1 and -2 respectively                   
                  if(iqa.ne.0) then
                     if(sourcea_type(i).eq.-1) then
                        vtotw(i)=0.0d00
c     vtota(i)=0.0d00
                     endif
                  endif            
                  if(iqco2.ne.0) then
                     if(sourceco2_type(i).lt.0) then
                        vtotw(i)=0.0d00
c     vtotco2(i)=0.0d00
                     endif
                  endif
                  if(ifd.ne.0) then
                     if(drainar_type(i).lt.0) then
                        vtotw(i)=0.0d00
c     atotd(i)=0.0d00
                     endif
                  endif
c     
c     check for timestep size change
c     
                  if(timestep_type(i).ne.0) then
                     day=timestep(abs(time_type(i)),i)
c gaz 111720 modification for first 'ts'  at time = 0
                     daynew=day
                     if(l.ne.0) then
                      days=days0+day
                      dtot_next = 0.0
                     endif
                  endif
c gaz 012322
c     
c     check for equation tolerance change
c
                  if(ieqtol.ne.0) then
                   if(eqtol_type(i).ne.0) then
                     tmchboun=eqtolerance(abs(time_type(i)),i)
                     tmch = tmchboun
                   endif
                  endif                   
c gaz 012122
c     
c     check for timestep maximum size change
c
                  if(itsmax.ne.0) then
                   if(tsmax_type(i).ne.0) then
                     daymaxboun=timestepmax(abs(time_type(i)),i)
                     daymax = daymaxboun
                     daynew=min(day,daymaxboun)
                     day = daynew
                     if(l.ne.0) then
                      days=days0+day
                      dtot_next = 0.0
                     endif
                   endif
                  endif                   
               endif
            enddo

            do i=1,n
               iimodel=node_model(i)
               if(iimodel.gt.0) then
                  iwght = weight_type(iimodel)
                  if(iwght.eq.5) then
                     if(cord(i,igrav).gt.zmax(iimodel)) then
                        zmax(iimodel)=cord(i,igrav)
                        pump_node(iimodel)=i
                     endif
                  endif
               endif
	    end do

            do i=1,n
               iimodel=node_model(i)
               if(iimodel.gt.0) then
                  if(icm.ne.0) then
                     if(node_ch_model_type(iimodel).gt.0) then
                        imodel=node_ch_model(abs(time_type(iimodel)),
     &                       iimodel)
                     else
                        imodel = iimodel
                     endif
                  else
                     imodel = iimodel
                  endif
               else
                  imodel = iimodel
               endif
c     gaz steady state management 10-30-04
               if(isty.ne.0.and.imodel.gt.0) then
                  istea = steady_type(imodel)
               else
                  istea = 1
               endif
               if(imodel.gt.0.and.istea.gt.0) then
c     zvd 16-Oct-08 Variable perms from V2.26 (update perms before 
c     weights in case using permeability weighting)
                  if(ixperm.ne.0) then
                     if(permx_type(imodel).ne.0) then
c     If permx value is zero, don't update
                        if(permx(abs(time_type(imodel)),imodel) .ne. 0.)
     &                       then
                           pnx(i)=permx(abs(time_type(imodel)),imodel)
                           if(pnx(i).ge.0.0d00) then
                              pnx(i) = max (zero_t, pnx(i))
                           else
                              pnx(i) = 1.0d00/10.0d00**(abs(pnx(i)))
                           endif
                           pnx(i) = pnx(i) * 1.d+06
                        end if
                     end if
                  end if
                  if(iyperm.ne.0) then
                     if(permy_type(imodel).ne.0) then
c     If permy value is zero, don't update
                        if(permy(abs(time_type(imodel)),imodel) .ne. 0.)
     &                       then
                           pny(i)=permy(abs(time_type(imodel)),imodel)
                           if(pny(i).ge.0.0d00) then
                              pny(i) = max (zero_t, pny(i))
                           else
                              pny(i) = 1.0d00/10.0d00**(abs(pny(i)))
                           endif
                           pny(i) = pny(i) * 1.d+06
                        end if
                     end if
                  end if
                  if(izperm.ne.0) then
                     if(permz_type(imodel).ne.0) then
c     If permz value is zero, don't update
                        if(permz(abs(time_type(imodel)),imodel) .ne. 0.)
     &                       then     
                           pnz(i)=permz(abs(time_type(imodel)),imodel)
                           if(pnz(i).ge.0.0d00) then
                              pnz(i) = max (zero_t, pnz(i))
                           else
                              pnz(i) = 1.0d00/10.0d00**(abs(pnz(i)))
                           endif
                           pnz(i) = pnz(i) * 1.d+06
                        end if
                     end if
                  end if
c     gaz new weighting capability 10-26-04
                  iwght = weight_type(imodel)
                  if(iwght.le.1) then
                     vtotw(imodel)=vtotw(imodel)+sf*vol_nd(i)
                  elseif(iwght.eq.2) then
                     perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                     vtotw(imodel)=vtotw(imodel)+perm_tot*sf*vol_nd(i)
                  elseif(iwght.eq.4) then
                     if(ihead.ne.0) then
	                rldum = max(0.d0,rlxyf(i)-rlptol)
	                perm_tot = sqrt(rldum**2*
     &                       (pnx(i)**2+pny(i)**2+pnz(i)**2))
                     else
                        perm_tot = sqrt(rlf(i)**2*
     &                       ((pnx(i)**2+pny(i)**2+pnz(i)**2)))
                     endif
                     vtotw(imodel)=vtotw(imodel)+perm_tot*sf*vol_nd(i)
                  elseif(iwght.eq.5) then
c     this takes into account distance from pump
                     rldum = max(0.d0,rlxyf(i)-rlptol)
                     perm_tot = sqrt(rldum**2*
     &                    (pnx(i)**2+pny(i)**2+pnz(i)**2))
                     dd=zmax(imodel)-cord(i,igrav)
                     vtotw(imodel)=vtotw(imodel)+perm_tot*sf*vol_nd(i)*
     &                    exp(-2.45-.01*dd)
C     
C     09/15/06 CWG
C     case 6 - 11
C     Volume/distance weighting, where distance is the [xyz]_d[xyz] = [xyz]_max - [xyz]_min of
C     the bounding box of the midpoint of all edges incident upon node i.
C     The extra if statement inside each case is to avoid a divide by
C     zero if [xyz]_d[xyz] is very small.
C     09/15/06 CWG
C     
                  elseif(iwght.eq.6) then
                     if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        vtotw(imodel)=vtotw(imodel)+
     &                       (sf*vol_nd(i)/dx_bbox(i))
                     endif
                  elseif(iwght.eq.7) then
                     if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        vtotw(imodel)=vtotw(imodel)+
     &                       (sf*vol_nd(i)/dy_bbox(i))
                     endif
                  elseif(iwght.eq.8) then
                     if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        vtotw(imodel)=vtotw(imodel)+
     &                       (sf*vol_nd(i)/dz_bbox(i))
                     endif
                  elseif(iwght.eq.9) then
                     if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                        vtotw(imodel) = vtotw(imodel)+
     &                       perm_tot*(sf*vol_nd(i)/dx_bbox(i))
                     endif
                  elseif(iwght.eq.10) then
                     if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                        vtotw(imodel)=vtotw(imodel)+
     &                       perm_tot*(sf*vol_nd(i)/dy_bbox(i))
                     endif
                  elseif(iwght.eq.11) then
                     if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)then
                        perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                        vtotw(imodel)=vtotw(imodel)+
     &                       perm_tot*(sf*vol_nd(i)/dz_bbox(i))
                     endif
                  endif
                  if(isf.ne.0) then
                     if(seepfac_type(imodel).gt.0) then
                        sp(i)=seepfac(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(ifd.ne.0) then
                     if(drainar_type(imodel).gt.0) then
                        drain(i)=drainar(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(iqw.ne.0) then
                     if(sourcew_type(imodel).gt.0) then
                        qw(i)=sourcew(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(iqco2.ne.0) then
                     if(sourceco2_type(imodel).gt.0) then
                        qco2b(i)=sourceco2(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(itempb2.ne.0) then
                     if(esourceco2_type(imodel).gt.0) then
                        eco2b(i)=esourceco2(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(iqf.ne.0) then
                     if(sourcef_type(imodel).gt.0) then
                        qw(i)=qw0(i)*sourcef(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(iqa.ne.0) then
                     if(sourcea_type(imodel).gt.0) then
                        qa(i)=sourcea(abs(time_type(imodel)),imodel)
                     endif
                  endif
c gaz 111418 need separate arrays for "air fraction of sw" and "sa"  : uses lt.0              
                  if(ixa.ne.0) then
                     if(sourcea_type(imodel).eq.-2) then
                        qaxf(i)=sourcea(abs(time_type(imodel)),imodel)
                     endif
                  endif                  
                  if(iha.ne.0) then
                     if(humid_type(imodel).lt.0) then
                        huma(i)=-humid(abs(time_type(imodel)),imodel)
                        phuma(i) = pressure_std
                        thuma(i) = temperature_std
                     endif
                  endif
                  if(iha.ne.0) then
                     if(humid_type(imodel).gt.0) then
                        huma(i)=humid(abs(time_type(imodel)),imodel)
                        phuma(i) = pressure_std
                        thuma(i) = temperature_std
                     endif
                  endif
                  if(ipha.ne.0) then
                     if(phumid_type(imodel).gt.0) then
                        phuma(i)=phumid(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(itha.ne.0) then
                     if(thumid_type(imodel).gt.0) then
                        thuma(i)=thumid(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(iqenth.ne.0) then
                     if(sourcee_type(imodel).gt.0) then
                        qenth(i)=sourcee(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(isatb.ne.0) then
                     if(saturation_type(imodel).lt.0) then
                        satb(i)=-saturation(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(isatb.ne.0) then
                     if(saturation_type(imodel).gt.0) then
                        satb(i)=saturation(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(ipresw.ne.0) then
                     if(pressurew_type(imodel).lt.0) then
                        presw(i)=-pressurew(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(ipresw.ne.0) then
                     if(pressurew_type(imodel).gt.0) then
                        presw(i)=pressurew(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(ipresa.ne.0) then
                     if(pressurea_type(imodel).lt.0) then
                        presa(i)=-pressurea(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(ipresa.ne.0) then
                     if(pressurea_type(imodel).gt.0) then
                        presa(i)=pressurea(abs(time_type(imodel)),
     &                       imodel)
                     endif
                  endif
                  if(imped.ne.0) then
                     if(impedance_type(imodel).gt.0) then
                        wellim(i)=impedance(abs(time_type(imodel)),
     &                       imodel)*1.d06
                     endif
                  endif
                  if(ienth.ne.0) then
                     if(enthalpy_type(imodel).ne.0) then
                        enth(i)=enthalpy(abs(time_type(imodel)),imodel)
                     endif
                  endif
                  if(itempb.ne.0) then
                     if(temperature_type(imodel).ne.0) then
                        tempb(i)=temperature(abs(time_type(imodel)),
     &                       imodel)
                     endif         
                  endif         
c     added initial value capability 12/3/98 GAZ
                  if(itempb_ini.ne.0.and.node_model(i).gt.0) then
                     if(temperature_ini_type(imodel).ne.0) then
                        tempb_ini(i)=temperature_ini(abs(time_type
     &                       (imodel)),imodel)
                     endif         
                  endif         
                  if(isatb_ini.ne.0.and.node_model(i).gt.0) then
                     if(saturation_ini_type(imodel).ne.0) then
                        satb_ini(i)=saturation_ini(abs(time_type
     &                       (imodel)),imodel)
                     endif
                  endif
                  if(ipresw_ini.ne.0.and.node_model(i).gt.0) then
                     if(pressurew_ini_type(imodel).ne.0) then
                        presw_ini(i)=pressurew_ini(abs(time_type
     &                       (imodel)),imodel)
                     endif
                  endif
                  if(ipresa_ini.ne.0.and.node_model(i).gt.0) then
                     if(pressurea_ini_type(imodel).ne.0) then
                        presa_ini(i)=pressurea_ini(abs(time_type
     &                       (imodel)),imodel)
                     endif
                  endif
               endif
            enddo


            do i=1,n
               imodel=node_model(i)
               if(imodel.gt.0) then
c     gaz new weighting capability 10-26-04
c     gaz steady state management 10-30-04
                  if(isty.ne.0.and.imodel.gt.0) then
                     istea = steady_type(imodel)
                  else
                     istea = 1
                  endif
               endif
               if(imodel.gt.0.and.istea.gt.0) then
                  iwght = weight_type(imodel)
                  if(iwght.le.1) then
                     vfac = sf*vol_nd(i)/vtotw(imodel)
                  elseif(iwght.eq.2) then
                     perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                     vfac = sf*vol_nd(i)*perm_tot/vtotw(imodel)
                  elseif(iwght.eq.3) then
                     vfac = wgt_area(i)
                  elseif(iwght.eq.4) then
                     rldum = max(0.d0,rlxyf(i)-rlptol)
                     perm_tot = sqrt(rldum**2*
     &                    (pnx(i)**2+pny(i)**2+pnz(i)**2))
                     vfac = sf*vol_nd(i)*perm_tot/vtotw(imodel)
                  elseif(iwght.eq.5) then
                     rldum = max(0.d0,rlxyf(i)-rlptol)
                     perm_tot = sqrt(rldum**2*
     &                    (pnx(i)**2+pny(i)**2+pnz(i)**2))
                     dd=zmax(imodel)-cord(i,igrav)
                     vfac = sf*vol_nd(i)*perm_tot*
     &                    exp(-2.45-.01*dd)/vtotw(imodel)
                  elseif((iwght.ge.6).and.(iwght.le.11))then
C     Begin Volume/distance weighted cases
                     if(iwght.eq.6) then
                        if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)
     &                       then
                           vfac = (sf*vol_nd(i)/dx_bbox(i))/
     &                          vtotw(imodel)
                        else
                           vfac = 0.0d0
                        endif
                     elseif(iwght.eq.7) then
                        if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)
     &                       then
                           vfac = (sf*vol_nd(i)/dy_bbox(i))/
     &                          vtotw(imodel)
                        else
                           vfac = 0.0d0
                        endif
                     elseif(iwght.eq.8) then
                        if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*1.e-7)
     &                       then
                           vfac = (sf*vol_nd(i)/dz_bbox(i))/
     &                          vtotw(imodel)
                        else
                           vfac = 0.0d0
                        endif
                     elseif((iwght.ge.9).and.(iwght.le.11)) then
C     Begin perm weighted cases
                        perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                        if(iwght.eq.9) then
                           if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                             pnz(i)**2)
                              vfac = 
     &                             perm_tot*(sf*vol_nd(i)/dx_bbox(i))/
     &                             vtotw(imodel)
                           else
                              vfac = 0.0d0
                           endif
                        elseif(iwght.eq.10) then
                           if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                             pnz(i)**2)
                              vfac = 
     &                             perm_tot*(sf*vol_nd(i)/dy_bbox(i))/
     &                             vtotw(imodel)
                           else
                              vfac = 0.0d0
                           endif
                        elseif(iwght.eq.11) then
                           if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                             pnz(i)**2)
                              vfac = 
     &                             perm_tot*(sf*vol_nd(i)/dz_bbox(i))/
     &                             vtotw(imodel)
                              if(v11norm .le. 0) v11norm = vfac
                           else
                              vfac = 0.0d0
                           endif
                        endif
C     End perm weighted cases
                     endif
C     End Volume/distance weighted cases
                  endif
                  if(ifd.ne.0) then
                     if(drainar_type(imodel).lt.0) then
                        drain(i)=drainar(abs(time_type(imodel)),imodel)
     &                       *vfac
                     endif
                  endif
                  if(iqw.ne.0) then
                     if(sourcew_type(imodel).lt.0) then
                        qw(i)=sourcew(abs(time_type(imodel)),imodel)*
     &                       vfac
                     endif
                  endif
                  if(iqco2.ne.0) then
                     if(sourceco2_type(imodel).lt.0) then
                        qco2b(i)=sourceco2(abs(time_type(imodel)),
     &                       imodel)*vfac
                     endif
                  endif
                  if(itempb2.ne.0) then
                     if(esourceco2_type(imodel).lt.0) then
                        eco2b(i)=esourceco2(abs(time_type(imodel)),
     &                       imodel)*vfac
                     endif
                  endif
                  if(iqa.ne.0) then
                     if(sourcea_type(imodel).eq.-1) then
c     vfac = sf*vol_nd(i)/vtota(imodel)
                        qa(i)=sourcea(abs(time_type(imodel)),imodel)*
     &                       vfac
                     endif
                  endif
                  if(iqenth.ne.0) then
                     if(sourcee_type(imodel).lt.0) then
                        qenth(i)=sourcee(abs(time_type(imodel)),imodel)
     &                       *vfac
                     endif
                  endif
C     
               endif
            enddo
C     
C     Adjust time step if necessary.
C     
            call time_adjust(mmodel,modmin,time_type,time_cycle,
     *           time,days0,day,days,daynew,maxtimes,iptty,iout,l,
     *           lchange,node_model,steady_type,isty,n)
C     
C     Enter a second loop over all the nodes to determine
C     source/sink terms in cases where there is time interpolation
C     of the source/sink terms.
C     
            if_time_interpolate = 0
            do i=1,mmodel
               if(time_interpolate(i) .ne. 0)then
                  if_time_interpolate = if_time_interpolate + 1
               endif
            enddo

            if(if_time_interpolate .gt. 0)then
               do i=1,n
                  imodel=node_model(i)
                  if(imodel.gt.0) then
c     gaz new weighting capability 10-26-04
c     gaz steady state management 10-30-04
                     if(isty.ne.0.and.imodel.gt.0) then
                        istea = steady_type(imodel)
                     else
                        istea = 1
                     endif
                  endif
                  if(imodel.gt.0.and.istea.gt.0) then
                     iwght = weight_type(imodel)
                     if(iwght.le.1) then
                        vfac = sf*vol_nd(i)/vtotw(imodel)
                     elseif(iwght.eq.2) then
                        perm_tot = sqrt(pnx(i)**2+pny(i)**2+pnz(i)**2)
                        vfac = sf*vol_nd(i)*perm_tot/vtotw(imodel)
                     elseif(iwght.eq.3) then
                        vfac = wgt_area(i)
                     elseif(iwght.eq.4) then
                        rldum = max(0.d0,rlxyf(i)-rlptol)
                        perm_tot = sqrt(rldum**2*
     &                       (pnx(i)**2+pny(i)**2+pnz(i)**2))
                        vfac = sf*vol_nd(i)*perm_tot/vtotw(imodel)
                     elseif(iwght.eq.5) then
                        rldum = max(0.d0,rlxyf(i)-rlptol)
                        perm_tot = sqrt(rldum**2*
     &                       (pnx(i)**2+pny(i)**2+pnz(i)**2))
                        dd=zmax(imodel)-cord(i,igrav)
                        vfac = sf*vol_nd(i)*perm_tot*
     &                       exp(-2.45-.01*dd)/vtotw(imodel)
                     elseif((iwght.ge.6).and.(iwght.le.11))then
C     Begin Volume/distance weighted cases
                        if(iwght.eq.6) then
                           if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              vfac = (sf*vol_nd(i)/dx_bbox(i))/
     &                             vtotw(imodel)
                           else
                              vfac = 0.0d0
                           endif
                        elseif(iwght.eq.7) then
                           if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              vfac = (sf*vol_nd(i)/dy_bbox(i))/
     &                             vtotw(imodel)
                           else
                              vfac = 0.0d0
                           endif
                        elseif(iwght.eq.8) then
                           if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                          1.e-7) then
                              vfac = (sf*vol_nd(i)/dz_bbox(i))/
     &                             vtotw(imodel)
                           else
                              vfac = 0.0d0
                           endif
                        elseif((iwght.ge.9).and.(iwght.le.11)) then
C     Begin perm weighted cases
                           perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                          pnz(i)**2)
                           if(iwght.eq.9) then
                              if(dx_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                             1.e-7) then
                                 perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                                pnz(i)**2)
                                 vfac = 
     &                                perm_tot*(sf*vol_nd(i)/dx_bbox(i))
     &                                /vtotw(imodel)
                              else
                                 vfac = 0.0d0
                              endif
                           elseif(iwght.eq.10) then
                              if(dy_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                             1.e-7) then
                                 perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                                pnz(i)**2)
                                 vfac = 
     &                                perm_tot*(sf*vol_nd(i)/dy_bbox(i))
     &                                /vtotw(imodel)
                              else
                                 vfac = 0.0d0
                              endif
                           elseif(iwght.eq.11) then
                              if(dz_bbox(i) .gt. (sf*vol_nd(i)**(1/3))*
     &                             1.e-7) then
                                 perm_tot = sqrt(pnx(i)**2+pny(i)**2+
     &                                pnz(i)**2)
                                 vfac = 
     &                                perm_tot*(sf*vol_nd(i)/dz_bbox(i))
     &                                /vtotw(imodel)
                              else
                                 vfac = 0.0d0
                              endif
                           endif
C     End perm weighted cases
                        endif
C     End Volume/distance weighted cases
                     endif
C     
C     Linear interpolation between time and flux values.
C     This is turned on at input by use of ti_linear or cy_linear
C     instead of ti or cy keywords.
C     
                     if(time_interpolate(imodel) .eq. 1)then
C     
C     Compute the proportion of the time cycle that has
C     completed.
C     
C     days = time_start, time_factor1=1, time_factor2=0
C     days = time_end,   time_factor1=0, time_factor2=1
C     
C     time_start = time(abs(time_type(imodel))  ,imodel)
C     time_end   = time(abs(time_type(imodel))+1,imodel)
C     time_factor2 = 
                         
                        time_interpolate_used = .false.
                        if(days .ge. time(abs(time_type(imodel)),
     &                       imodel)) then
                           time_factor2 = (days - 
     &                          time(abs(time_type(imodel)),imodel))/
     &                          (time(abs(time_type(imodel))+1,imodel)
     &                          - time(abs(time_type(imodel)),imodel))
                           time_factor1 = 1.0d0 - time_factor2
                        else
C     
C     No action if the start time of the model has not been reached.
C     
                           time_factor1 = 0.0d0
                           time_factor2 = 0.0d0
                        endif

                        if(ifd.ne.0) then
                           if(drainar_type(imodel).lt.0) then
                              drain(i)= time_factor1*drain(i) + 
     &                             time_factor2*drainar(abs
     &                             (time_type(imodel))+1,imodel)*vfac
                              time_interpolate_used = .true.
                           endif
                        endif
                        if(iqw.ne.0) then
                           if(sourcew_type(imodel).lt.0) then
                              qw(i)   = time_factor1*qw(i) + 
     &                             time_factor2*sourcew(abs
     &                             (time_type(imodel))+1,imodel)*vfac
                              time_interpolate_used = .true.            
c gaz 012819 added  sw ("gt.0" change) (removed vfac term)                                
                           else if(sourcew_type(imodel).gt.0) then
                              qw(i)   = time_factor1*qw(i) + 
     &                             time_factor2*sourcew(abs
     &                             (time_type(imodel))+1,imodel) 
                              time_interpolate_used = .true. 
                           endif
                        endif
                        if(iqco2.ne.0) then
                           if(sourceco2_type(imodel).lt.0) then
                              qco2b(i)   = time_factor1*qw(i) + 
     &                             time_factor2*sourceco2(abs
     &                             (time_type(imodel))+1,imodel)*vfac
                              time_interpolate_used = .true. 
                           endif
                        endif
                        if(iqa.ne.0) then
                           if(sourcea_type(imodel).eq.-1) then
                              vfac = sf*vol_nd(i)/vtota(imodel)
                              qa(i)   = time_factor1*qa(i) + 
     &                             time_factor2*sourcea(abs
     &                             (time_type(imodel))+1,imodel)*vfac
                              time_interpolate_used = .true. 
c gaz 012819 added  sw ("gt.0" change) (removed vfac term)                                
                           elseif(sourcea_type(imodel).gt.0) then
                              qa(i)   = time_factor1*qa(i) + 
     &                             time_factor2*sourcea(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                           endif
                        endif
                        if(iqenth.ne.0) then
                           if(sourcee_type(imodel).lt.0) then
                              qenth(i)= time_factor1*qenth(i) + 
     &                             time_factor2*sourcee(abs
     &                             (time_type(imodel))+1,imodel)*vfac
                              time_interpolate_used = .true. 
c gaz 012819 added  se("gt.0" change) (removed vfac term)                                
                           elseif(sourcee_type(imodel).gt.0) then
                              qenth(i) = time_factor1*qenth(i) + 
     &                             time_factor2*sourcee(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true.                              
                           endif
                        endif
c gaz 012819 added  air pressure (removed vfac term)   
                        if(ipresa.ne.0) then
                            if(pressurea_type(imodel).ne.0) then
                              presa(i)= time_factor1*presa(i) + 
     &                             time_factor2*pressurea(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                            continue  
                           endif
                        endif                         
c gaz 012819 added  pressure (removed vfac term)   
                        if(ipresw.ne.0) then
                            if(pressurew_type(imodel).ne.0) then
                              presw(i)= time_factor1*presw(i) + 
     &                             time_factor2*pressurew(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                            continue  
                           endif
                        endif                          
c gaz 012819 added  temperature (removed vfac term) 
c gaz 021619 added  head (uses temperature temperature arrays(can't use both)
c should be valid for both 'hd and hdo'                        
                        if(itempb.ne.0) then
                            if(temperature_type(imodel).ne.0) then
                              tempb(i)= time_factor1*tempb(i) + 
     &                             time_factor2*temperature(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                            continue  
                           endif
                        endif 
c gaz 021619 enthalpy and flowing enthalpy                        
                        if(ienth.ne.0) then
                            if(enthalpy_type(imodel).ne.0) then
                              enth(i)= time_factor1*enth(i) + 
     &                             time_factor2*enthalpy(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                            continue  
                           endif
                        endif  
c gaz 021619 saturation                        
                        if(isatb.ne.0) then
                            if(saturation_type(imodel).ne.0) then
                              satb(i)= time_factor1*satb(i) + 
     &                             time_factor2*saturation(abs
     &                             (time_type(imodel))+1,imodel)
                              time_interpolate_used = .true. 
                            continue  
                           endif
                        endif                           
                     endif
                  endif
               enddo
               if(iout.ne.0 .and. boun_out) write(iout,*) ' '
               if(iptty.ne.0 .and. boun_out) write(iptty,*) ' '              
               do i=1,mmodel
c gaz 012819  write details only if used                   
                  if(time_interpolate(i) .ne. 0 .and.
     &                 time_interpolate_used) then
                     if(days .ge. time(abs(time_type(i)),i))then
                        time_factor2 = 
     &                       (days - time(abs(time_type(i)),  i))/
     &                       (time(abs(time_type(i))+1,i) - 
     &                       time(abs(time_type(i))  ,i))
                        time_factor1 = 1.0d0 - time_factor2
                     else
                        time_factor1 = 0.0d0
                        time_factor2 = 0.0d0
                     endif
                     if (iout .ne. 0 .and. boun_out) then
                        write(iout,20)' ti_linear t=',days,
     *                       ' model #=',i,
     *                       ' t1=',time(abs(time_type(i))  ,i),
     *                       ' t2=',time(abs(time_type(i))+1,i),
     *                       ' f1=',time_factor1,
     *                       ' f2=',time_factor2
                     endif
                     if(iptty .ne. 0 .and. boun_out) then
                        write(iptty,20)' ti_linear t=',days,
     *                       ' model #=',i,
     *                       ' t1=',time(abs(time_type(i))  ,i),
     *                       ' t2=',time(abs(time_type(i))+1,i),
     *                       ' f1=',time_factor1,
     *                       ' f2=',time_factor2
                     endif
c gaz 012819  write warning if ti_linear enabled with know appropriate variables                    
                  elseif(time_interpolate(i) .ne. 0 .and.
     &                 .not.time_interpolate_used) then 
                   if (iout .ne. 0 .and. boun_out) then
                     write(iout,21) ' ti_linear t=',days,' model #=',i,
     &                ' no appropriate keywords found (disabled)',
     &                ' only keywords allowed with ti_linear:',
     &                ' fd, sw, dsw, dsco2, sa, dsa, se, dse, t, ft,',
     &                ' hd, hdo, s'      
                   endif
                   if (iptty .ne. 0 .and. boun_out) then
                     write(iptty,21) ' ti_linear t=',days,' model #=',i,
     &                ' no appropriate keywords found (disabled)',
     &                ' only keywords allowed with ti_linear:',
     &                ' fd, sw, dsw, dsco2, sa, dsa, se, dse, t, ft,',
     &                ' hd, hdo, s'   
                   endif
                  endif
 20               format(a,g15.10,a,i4,a,g13.6,a,g13.6,a,g13.6,a,g13.6)
 21               format(a,g15.10,a,i4,a,/,a,a,a)                 
               enddo
            endif
         endif
      endif
      return
      end
