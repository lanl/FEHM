      subroutine  inhist
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
!D1 Read control information for writing history files.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 6-FEB-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inhist.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai
      use combi, only : prnt_flxzvar
      use comco2, only : icarb
      use comdi
c gaz 041620 need n0 and pci      
      use comfi, only : pci
      use comdti, only : n0
      use comrlp, only : ishisrlp, delta_sat, num_sat, sat_out
      use comsi, only : ihms
      use comwt, only : wt_flag
      use comxi
      use comzone
      use davidi
      implicit none

      character*5 flxzone
      character*8 hissfx, trcsfx, chtmp(3)
      character*32 cmsg(4), cdum
      character*80 chdum
      character*120 fname, root
      logical null1, opend
      integer i, iroot, msg(4), imsg(4), nwds, ishisfzz, ishiscfzz
      real*8 xmsg(4)

      hist_flag = .TRUE.

! Default time in days
      time_flag = 2
! Default output format
      form_flag = 0
      hissfx = '.his'
      trcsfx = '.trc'
! Default print for each time step
      nhist = 1
! Default print every 1 million days
      histime = 1.e6
! Other flags initialized to 0 -> no output
      ishisg = 0
      ishisp = 0
      ishisp2 = 0
      ishisp3 = 0
      ishist = 0
      ishishd = 0
      ishiss = 0
      ishisf = 0
      ishisfa = 0
      ishise = 0
      ishisef = 0
      ishisd = 0
      ishisv = 0
      ishishm = 0
      ishisfz = 0
      ishiswt = 0
      ishisc = 0
      ishiswc = 0
      ishiscm = 0
      ishiscmf = 0
      ishiscmd = 0
      ishiscfz = 0
      ishiscsl = 0
      ishiscsg = 0
      ishisdisx = 0
      ishisdisy = 0
      ishisdisz = 0
      ishisstr = 0
      ishisstrx = 0
      ishisstry = 0
      ishisstrz = 0
      wt_flag = 0

      if (null1(root_name)) then
! Find root of history file name or if not present, input file name
         if (nmfil(8) .ne. nmfily(3) .and. nmfil(8) .ne. ' ') then
            call file_prefix(nmfil(8), iroot)
            if (iroot .gt. 100) iroot = 100
            root(1:iroot) = nmfil(8)(1:iroot)
         else
            if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
               call file_prefix(nmfil(5), iroot)
               if (iroot .gt. 100) iroot = 100
               root(1:iroot) = nmfil(5)(1:iroot)
            else
               if (nmfil(2)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil2 file: ', nmfil(2),
     .                 ' unable to determine history file prefix'
                  stop
               else
                  call file_prefix(nmfil(2), iroot)
                  if (iroot .gt. 100) iroot = 100
                  root(1:iroot) = nmfil(2)(1:iroot)
               end if
            end if
         endif
      else
         iroot = len_trim (root_name)
         if (iroot .gt. 100) iroot = 100
         root(1:iroot) = root_name(1:iroot)
      end if

      ishis = nufilb(8)
      fname =  root(1:iroot) // trim(hissfx)
      inquire (unit = ishis, opened = opend)
      if (.not. opend) then
         open (unit=ishis, file=fname, form='formatted')
         write(ishis, 6000) verno, jdate, jtime, trim(wdd)
      end if

! Backspace and reread macro line to see if an alternate output 
! format is specified
      backspace inpt
      read (inpt, '(a80)') chdum
      call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
      if (nwds .ge. 2 .and. msg(2) .eq. 3) then
         if (cmsg(2)(1:3) .eq. 'tec' .or. cmsg(2)(1:3) .eq. 'TEC') then
            form_flag = 1
            hissfx = '_his.dat'
         else if (cmsg(2)(1:3) .eq. 'sur' .or. cmsg(2)(1:3) .eq. 'SUR'
     &           .or. cmsg(2)(1:3) .eq. 'csv' .or. cmsg(2)(1:3) .eq. 
     &           'CSV') then
            form_flag = 2
            hissfx = '_his.csv'
         end if
      end if
! While loop to read history output selections
      do
         read (inpt, '(a80)') chdum
         if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or. 
     &        chdum(1:3) .eq. 'END') exit
         fname = ''
         if ((chdum(1:3) .eq. 'dis' .or. chdum(1:3) .eq. 'DIS' .or.
     &        chdum(1:3) .eq. 'str' .or. chdum(1:3) .eq. 'STR') 
     &        .and. ihms .ne. -3) then
            if (iout .ne. 0) write(iout, 6080) istrs_coupl, 
     &           trim(chdum)
            if (iptty .ne. 0) write(iptty, 6080) istrs_coupl, 
     &           trim(chdum)
            chdum = ''
         end if

         flags : select case (chdum(1:3))
         case ('tec', 'TEC')
! Output using tecplot format
               form_flag = 1
               hissfx = '_his.dat'
         case ('sur', 'SUR', 'csv', 'CSV')
! Output using surfer (csv) format
               form_flag = 2
               hissfx = '_his.csv'
         case ('yea','day','sec','hrs','hou','YEA','DAY','SEC','HRS',
     &              'HOU')
            if (chdum(1:3) .eq. 'yea' .or. chdum(1:3) .eq. 'YEA') then
! Output time in years
               time_flag = 1
            else if (chdum(1:3) .eq. 'day' .or. chdum(1:3) .eq. 'DAY') 
     &              then
! Output time in days
               time_flag = 2
            else if (chdum(1:3) .eq. 'sec' .or. chdum(1:3) .eq. 'SEC')  
     &              then
! Output time in seconds
               time_flag = 3
            else if (chdum(1:3) .eq. 'hrs' .or. chdum(1:3) .eq. 'HRS' 
     &              .or. chdum(1:3) .eq. 'hou' .or. chdum(1:3) .eq. 
     &              'HOU') then 
! Output time in hours
               time_flag = 4
            end if
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .eq. 3) then
! Set printout timestep interval
               nhist = imsg(2)
! Set printout time interval
               if (msg(3).eq.1) then
                  histime=imsg(3)
               else
                  histime = xmsg(3)
               endif
            end if
         case ('mpa', 'MPa', 'MPA','pre','PRE')
! Output pressures in MPa
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .gt. 1) then
! Determine which pressures to output
               nwd : select case (nwds)
               case (2)
                  if (cmsg(2)(1:3) .eq. 'tot' .or. cmsg(2)(1:3) .eq. 
     &                 'TOT' .or. cmsg(2)(1:3) .eq. 'wat' .or.
     &                 cmsg(2)(1:3) .eq. 'WAT' ) then
! Output total/water pressure
                     pres_flag = 1
                  else if ((cmsg(2)(1:3) .eq. 'air' .or. cmsg(2)(1:3) 
     &                    .eq. 'AIR') .and. ihead .ne. 1) then
! Output air/vapor pressure
                     pres_flag = 4
                  else if ((cmsg(2)(1:3) .eq. 'cap' .or. cmsg(2)(1:3) 
     &                    .eq. 'CAP') .and. (irdof .ne. 13 .or.
     &                    ifree .ne. 0)) then
! Output capillary pressure
                     pres_flag = 5
                  else if ((cmsg(2)(1:3) .eq. 'co2' .or. cmsg(2)(1:3) 
     &                    .eq. 'CO2') .and. icarb .ne. 0) then
! Output CO2 pressure
                     pres_flag = 8
                  else
                     pres_flag = 1
                     if (iout .ne. 0) write(iout, 6030)
                     if (iptty .ne. 0) write(iptty, 6030)
                  end if
               case (3)
                  if (ihead .eq. 1) then
                     pres_flag = 1
                     if (iout .ne. 0) write(iout, 6030)
                     if (iptty .ne. 0) write(iptty, 6030)      
                  else if ((cmsg(2)(1:3) .eq. 'tot' .or. cmsg(3)(1:3) 
     &                    .eq. 'tot' .or. cmsg(2)(1:3) .eq. 'TOT' 
     &                    .or. cmsg(3)(1:3) .eq. 'TOT') .and. 
     &                    (cmsg(2)(1:3) .eq. 'air' .or. cmsg(3)(1:3) 
     &                    .eq. 'air' .or. cmsg(2)(1:3) .eq. 'AIR'
     &                    .or. cmsg(3)(1:3) .eq. 'AIR')) then
                     pres_flag = 2
                  else if ((cmsg(2)(1:3) .eq. 'tot' .or. cmsg(3)(1:3) 
     &                    .eq. 'tot'.or. cmsg(2)(1:3) .eq. 'TOT' 
     &                    .or. cmsg(3)(1:3) .eq. 'TOT') .and.  
     &                    (cmsg(2)(1:3) .eq. 'cap' .or. cmsg(3)(1:3) 
     &                    .eq. 'cap' .or. cmsg(2)(1:3) .eq. 'CAP'
     &                    .or. cmsg(3)(1:3) .eq. 'CAP') .and. 
     &                    (irdof .ne. 13 .or. ifree .ne. 0)) then
                     pres_flag = 6
                  else if ((cmsg(2)(1:3) .eq. 'air' .or. cmsg(3)(1:3) 
     &                    .eq. 'air' .or. cmsg(2)(1:3) .eq. 'AIR'
     &                    .or. cmsg(3)(1:3) .eq. 'AIR') .and. 
     &                    (cmsg(2)(1:3) .eq. 'cap'  .or. cmsg(3)(1:3) 
     &                    .eq. 'cap' .or. cmsg(2)(1:3) .eq. 'CAP'
     &                    .or. cmsg(3)(1:3) .eq. 'CAP') .and. 
     &                    (irdof .ne. 13 .or. ifree .ne. 0))  then
                     pres_flag = 7
                  else
                     pres_flag = 2
                     if (iout .ne. 0) write(iout, 6040)
                     if (iptty .ne. 0) write(iptty, 6040)
                  end if
               case (4)
                  if (ico2 .gt. 0 .or. irdof .ne. 13 .or. ifree .ne. 0)
     &                 then
! Output total, air/vapor and capillary pressure
                     pres_flag = 3
                  else
                     if (ihead .eq. 1) then
! Output total/water pressure
                        pres_flag = 1
                        if (iout .ne. 0) write(iout, 6030)
                        if (iptty .ne. 0) write(iptty, 6030)
                     else
! Output Total/water and  air/vapor pressure
                        pres_flag = 2
                        if (iout .ne. 0) write(iout, 6040)
                        if (iptty .ne. 0) write(iptty, 6040)
                     end if
                  end if
               end select nwd
            else
! Defaults determined from run parameters
               if (icarb .ne. 0) then
! Output CO2 and water pressure
                  pres_flag = 9
               else if (ico2 .gt. 0 .or. irdof .ne. 13) then
! Output total, air/vapor and capillary pressure
                  pres_flag = 3
               else
                  if (ihead .eq. 1) then
! Output total/water pressure
                     pres_flag = 1
                  else
! Output Total/water and  air/vapor pressure
                     pres_flag = 2
                  end if
               end if
               end if
c gaz 041520 make sure pci() iz allocated(n) for pres_flag 1,2,3 (pres_flag 4 OK)
            if(allocated(pci).and.ico2.le.0) then
             deallocate(pci)
             allocate(pci(n0))
             pci= 0.0d0
            else if(.not.allocated(pci)) then
             allocate(pci(n0))
             pci= 0.0d0                
            endif
            ishisp = ishis + 100
            chtmp = ''
            select case (pres_flag)
            case (1)
! One output variable, total/water pressure
               chtmp(1) = '_presWAT'
            case (2)
! Two output variables, total/water and air/vapor pressure
               chtmp(1) = '_presWAT'
               chtmp(2) = '_presAIR'
            case (3)
! Three output variables (total, air/vapor and capillary pressure)
               chtmp(1) = '_presWAT'
               chtmp(2) = '_presVAP'
               chtmp(3) = '_presCAP'
            case (4)
! One output variable, air/vapor pressure
               chtmp(1) = '_presAIR'
            case (5)
! One output variable, capillary pressure
               chtmp(1) = '_presCAP'
            case (6)
! Two output variables, total/water and capillary pressure
               chtmp(1) = '_presTOT'
               chtmp(2) = '_presCAP'
            case (7)
! Two output variables, air/vapor and capillary pressure
               chtmp(1) = '_presVAP'
               chtmp(2) = '_presCAP'
            case (8)
! One output variable, CO2 pressure
               chtmp(1) = '_presCO2'
            case (9)
! Two output variables, CO2 and water pressure
               chtmp(1) = '_presCO2'
               chtmp(2) = '_presWAT'
            end select
            fname =  root(1:iroot) // trim(chtmp(1)) // hissfx
            open (unit=ishisp, file=fname, form='formatted')
            if (pres_flag .ne. 1 .and. pres_flag .ne. 4 .and.
     &           pres_flag .ne. 5 .and. pres_flag .ne. 8) then
               ishisp2 = ishisp + 1
               fname = ''
               fname =  root(1:iroot) // trim(chtmp(2)) // hissfx
               open (unit=ishisp2, file=fname, form='formatted')
               if (pres_flag .eq. 3) then
                  ishisp3 = ishisp2 + 1
                  fname = ''
                  fname =  root(1:iroot) // trim(chtmp(3)) // hissfx
                  open (unit=ishisp3, file=fname, form='formatted')
               end if
            end if
            select case (form_flag)
            case (0)
               write(ishisp, 6000) verno, jdate, jtime, trim(wdd)
               if (ishisp2 .ne. 0) write(ishisp2, 6000) verno, jdate, 
     &              jtime, trim(wdd)
               if (ishisp3 .ne. 0) write(ishisp3, 6000) verno, jdate, 
     &              jtime, trim(wdd)
            case (1)
               write(ishisp, 6005) verno, jdate, jtime
               if (ishisp2 .ne. 0) write(ishisp2, 6005) verno, jdate, 
     &              jtime
               if (ishisp3 .ne. 0) write(ishisp3, 6005) verno, jdate, 
     &              jtime
            end select
            if (out_zones)  then
               ozflag = ozflag + 1
               pflag = ozflag
            end if

         case ('deg', 'DEG', 'tem', 'TEM')
! Output temperature in degrees C
            fname =  root(1:iroot) // '_temp' // hissfx
            ishist = ishis + 110
            open (unit=ishist, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishist, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishist, 6005) verno, jdate, jtime
            end select
            if (out_zones) then
               ozflag = ozflag + 1
               tflag = ozflag
            end if
         case ('met', 'MET', 'hea', 'HEA')
! Output heads in m
            if (ihead .ne. 0 .or. ichead .ne. 0) then
               if (ishishd .ne. 0) then
                  if (iout .ne. 0) write(iout, 6025) 'feet'
                  if (iptty .ne. 0) write(iptty, 6025) 'feet'
               else
                  fname =  root(1:iroot) // '_head' // hissfx
                  ishishd = ishis + 120
                  call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
                  if (nwds .ge. 2) then
                     if (msg(2) .eq. 3) then
                        if (cmsg(2) .eq. 'nodes') ishishd = ishishd + 1
                        if (nwds .ge. 3) then
                           if (msg(3) .eq. 1) wt_flag = imsg(3)
                        end if
                     else if (msg(2) .eq. 1) then
                        wt_flag = imsg(2)
                     end if
                  end if
                  if (out_zones) then
                     ozflag = ozflag + 1
                     hflag = ozflag
                  end if
                  open (unit=ishishd, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishishd, 6000) verno, jdate, jtime, trim(wdd)
                  case (1)
                     write(ishishd, 6005) verno, jdate, jtime
                  end select
               end if
            else
               if (iout .ne. 0) write(iout, 6020)
               if (iptty .ne. 0) write(iptty, 6020)
            end if
         case ('fee', 'FEE')
! Output heads in ft
            if (ihead .ne. 0 .or. ichead .ne. 0) then
               if (ishishd .ne. 0) then
                  if (iout .ne. 0) write(iout, 6025) 'meters'
                  if (iptty .ne. 0) write(iptty, 6025) 'meters'
               else
                  fname =  root(1:iroot) // '_head' // hissfx
                  ishishd = ishis + 125
                  open (unit=ishishd, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishishd, 6000) verno, jdate, jtime, trim(wdd)
                  case (1)
                     write(ishishd, 6005) verno, jdate, jtime
                  end select
                  if (out_zones) then
                     ozflag = ozflag + 1
                     hflag = ozflag
                  end if
               end if
            else
               if (iout .ne. 0) write(iout, 6020)
               if (iptty .ne. 0) write(iptty, 6020)
            end if
         case ('sat', 'SAT')
! Output saturations
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               fname =  root(1:iroot) // '_satr' // hissfx
               ishiss = ishis + 130
               open (unit=ishiss, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishiss, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishiss, 6005) verno, jdate, jtime
               end select
            else 
               if (iout .ne. 0) write(iout, 6010)
               if (iptty .ne. 0) write(iptty, 6010)
            end if 
         case ('wco', 'WCO')
! Output water content
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               fname =  root(1:iroot) // '_wcon.his'
               ishiswc = ishis + 135
               open (unit=ishiswc, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishiswc, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishiswc, 6005) verno, jdate, jtime
               end select
            else 
               if (iout .ne. 0) write(iout, 6010)
               if (iptty .ne. 0) write(iptty, 6010)
            end if 
         case ('fla', 'FLA')
! Output flow in kg/s
! Not available for water and water vapor (ico2=0)
          if(ico2.ne.0) then
            fname =  root(1:iroot) // '_floa' // hissfx
            ishisfa = ishis + 137
            open (unit=ishisfa, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishisfa, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishisfa, 6005) verno, jdate, jtime
            end select
          else
           if(iout.ne.0) write (iout,*) 
     &      '>>>> not output for gas src/sink <<<<'
           if(iptty.ne.0) write (iptty,*) 
     &      '>>>> not output for gas src/sink <<<<'         
          endif 
         case ('kgs ', 'KGS ', 'flo', 'FLO')
! Output flow in kg/s
            fname =  root(1:iroot) // '_flow' // hissfx
            ishisf = ishis + 140
            open (unit=ishisf, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishisf, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishisf, 6005) verno, jdate, jtime
            end select
         case ('ent', 'ENT')
! Output enthalpy in MJ/kg
            fname =  root(1:iroot) // '_enth' // hissfx
            ishise = ishis + 150
            open (unit=ishise, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishise, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishise, 6005) verno, jdate, jtime
            end select
c            if (out_zones) then
c               ozflag = ozflag + 1
c               eflag = ozflag
c            end if
         case ('efl', 'EFL', 'mjs', 'MJS')
! Output enthalpy flow in MJ/s
            fname =  root(1:iroot) // '_eflow' // hissfx
c           ishise = ishis + 150
            ishise = ishis + 155
            open (unit=ishisef, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishisef, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishisef, 6005) verno, jdate, jtime
            end select
            if (out_zones) then
               ozflag = ozflag + 1
               eflag = ozflag
            end if

         case ('den', 'DEN')
! Output density
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
! Determine which densities to output
            select case (nwds)
            case (2)
               if (cmsg(2)(1:3) .eq. 'wat' .or. 
     &              cmsg(2)(1:3) .eq. 'WAT') then
! Output water density
                  den_flag = 2
               else if (cmsg(2)(1:3) .eq. 'air' .or. 
     &                 cmsg(2)(1:3) .eq. 'AIR') then
! Output air/vapor density
                  den_flag = 3
               else if (cmsg(2)(1:3) .eq. 'co2' .or. 
     &                 cmsg(2)(1:3) .eq. 'CO2') then
                  if (icarb .ne. 0) then
                     if (cmsg(2)(1:4) .eq. 'co2l' .or. 
     &                    cmsg(2)(1:3) .eq. 'CO2l') then
! Output CO2 liquid density
                        den_flag = 5
                     else if (cmsg(2)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(2)(1:3) .eq. 'CO2g') then
! Output CO2 gas density
                        den_flag = 6
                     else
! Output CO2 liquid and gas density
                        den_flag = 9
                     end if
                  else
                     den_flag = 0
                     if (iout .ne. 0) write(iout, 6067) 'density'
                     if (iptty .ne. 0) write(iptty, 6067) 'density'
                  end if
               end if
            case (3)
               if ((cmsg(2)(1:3) .eq. 'wat' .or. cmsg(3)(1:3) 
     &              .eq. 'wat' .or. cmsg(2)(1:3) .eq. 'WAT' 
     &              .or. cmsg(3)(1:3) .eq. 'WAT') .and. 
     &              (cmsg(2)(1:3) .eq. 'air' .or. cmsg(3)(1:3) 
     &              .eq. 'air' .or. cmsg(2)(1:3) .eq. 'AIR'
     &              .or. cmsg(3)(1:3) .eq. 'AIR')) then
! Output water and air/vapor density
                  den_flag = 1
               else if ((cmsg(2)(1:3) .eq. 'wat' .or. cmsg(3)(1:3) 
     &                 .eq. 'wat' .or. cmsg(2)(1:3) .eq. 'WAT' 
     &                 .or. cmsg(3)(1:3) .eq. 'WAT') .and.
     &                 (cmsg(2)(1:3) .eq. 'co2' .or. cmsg(3)(1:3) 
     &                 .eq. 'co2' .or. cmsg(2)(1:3) .eq. 'CO2'
     &                 .or. cmsg(3)(1:3) .eq. 'CO2')) then
                  if (icarb .ne. 0) then
                     if (cmsg(2)(1:4) .eq. 'co2l' .or. cmsg(3)(1:4) 
     &                    .eq. 'co2l' .or. cmsg(2)(1:4) .eq. 'CO2l'
     &                    .or. cmsg(3)(1:4) .eq. 'CO2l') then
! Output water and CO2 liquid density
                        den_flag = 7
                     else if (cmsg(2)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(3)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(2)(1:4) .eq. 'CO2g' .or. 
     &                       cmsg(3)(1:4) .eq. 'CO2g') then
! Output water and CO2 gas density
                        den_flag = 8
                     else
! Output water and CO2 densities
                        den_flag = 4
                     end if
                  else
                     den_flag = 2
                     if (iout .ne. 0) write(iout, 6067) 
     &                    'co2 density'
                     if (iptty .ne. 0) write(iptty, 6067) 
     &                    'co2 density'
                  end if
               else if ((cmsg(2)(1:4) .eq. 'co2l' .or. cmsg(3)(1:4) 
     &                 .eq. 'co2l' .or. cmsg(2)(1:4) .eq. 'CO2l' 
     &                 .or. cmsg(3)(1:4) .eq. 'CO2l') .and.
     &                 (cmsg(2)(1:4) .eq. 'co2g' .or. cmsg(3)(1:4) 
     &                 .eq. 'co2g' .or. cmsg(2)(1:4) .eq. 'CO2g'
     &                 .or. cmsg(3)(1:4) .eq. 'CO2g')) then
                  if (icarb .ne. 0) then
                     den_flag = 9
                  else
                     if (iout .ne. 0) write(iout, 6067) 
     &                    'co2 density'
                     if (iptty .ne. 0) write(iptty, 6067) 
     &                    'co2 density'
                  end if
               end if
            case default
! Defaults determined from run parameters
               if (icarb .ne. 0) then
! Output CO2 and water densities
                  den_flag = 4
               else if (ico2 .gt. 0 .or. irdof .ne. 13 .or. 
     &                 irdof .ne. 11) then
! Output water and air/vapor density
                  den_flag = 1
               else if (irdof .eq. 13) then
! Output water density
                  den_flag = 2
               else if (irdof .eq. 11) then
! Output air density
                  den_flag = 3
               end if
            end select

            if (den_flag .ne. 0) then
               ishisd = ishis + 170
               chtmp(1) = ''
               select case (den_flag)
               case (1)
                  chtmp(1) = '_denWAT'
                  chtmp(2) = '_denAIR'
               case(2)
                  chtmp(1) = '_denWAT'
               case(3)
                  chtmp(1) = '_denAIR'
               case(4)
                  chtmp(1) = '_denWAT'
                  chtmp(2) = '_denCO2l'
                  chtmp(3) = '_denCO2g'
               case(5)
                  chtmp(1) = '_denCO2l'
               case(6)
                  chtmp(1) = '_denCO2g'
               case(7)
                  chtmp(1) = '_denWAT'
                  chtmp(2) = '_denCO2l'
               case(8)
                  chtmp(1) = '_denWAT'
                  chtmp(2) = '_denCO2g'
                case(9)
                  chtmp(1) = '_denCO2l'
                  chtmp(2) = '_denCO2g'
               end select
               fname =  root(1:iroot) // trim(chtmp(1)) // hissfx
               open (unit=ishisd, file=fname, form='formatted')
               if (den_flag .eq. 1 .or. den_flag .eq. 4 .or. 
     &              den_flag .eq. 7 .or. den_flag .eq. 8 .or.
     &              den_flag .eq. 9) then
                  ishisd2 = ishisd + 1
                  fname = ''
                  fname =  root(1:iroot) // trim(chtmp(2)) // hissfx
                  open (unit=ishisd2, file=fname, form='formatted')
               end if
               if (den_flag .eq. 4) then
                  ishisd3 = ishisd2 + 1
                  fname = ''
                  fname =  root(1:iroot) // trim(chtmp(3)) // hissfx
                  open (unit=ishisd3, file=fname, form='formatted')
               end if

               select case (form_flag)
               case (0)
                  write(ishisd, 6000) verno, jdate, jtime, trim(wdd)
                  if (ishisd2 .ne. 0) write(ishisd2, 6000) verno, 
     &                 jdate, jtime, trim(wdd)
                  if (ishisd3 .ne. 0) write(ishisd3, 6000) verno, 
     &                 jdate, jtime, trim(wdd)
               case (1)
                  write(ishisd, 6005) verno, jdate, jtime
                  if (ishisd2 .ne. 0) write(ishisd2, 6005) verno, 
     &                 jdate, jtime
                  if (ishisd3 .ne. 0) write(ishisd3, 6005) verno, 
     &                 jdate, jtime
               end select
            end if

         case ('vis', 'VIS')
! Output viscosity
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
! Determine which viscosities to output
            select case (nwds)
            case (2)
               if (cmsg(2)(1:3) .eq. 'wat' .or. 
     &              cmsg(2)(1:3) .eq. 'WAT') then
! Output water viscosity
                  vis_flag = 2
               else if (cmsg(2)(1:3) .eq. 'air' .or. 
     &                 cmsg(2)(1:3) .eq. 'AIR') then
! Output air/vapor viscosity
                  vis_flag = 3
               else if (cmsg(2)(1:3) .eq. 'co2' .or. 
     &                 cmsg(2)(1:3) .eq. 'CO2') then
! CO2 viscosity 
                  if (icarb .ne. 0) then
                     if (cmsg(2)(1:4) .eq. 'co2l' .or. 
     &                    cmsg(2)(1:3) .eq. 'CO2l') then
                        vis_flag = 5
                     else if (cmsg(2)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(2)(1:3) .eq. 'CO2g') then
                        vis_flag = 6
                     else
                        vis_flag = 4
                     end if
                  else
                     vis_flag = 0
                     if (iout .ne. 0) write(iout, 6067) 'viscosity'
                     if (iptty .ne. 0) write(iptty, 6067) 'viscosity'
                  end if
               end if
            case (3)
               if ((cmsg(2)(1:3) .eq. 'wat' .or. cmsg(3)(1:3) 
     &              .eq. 'watt' .or. cmsg(2)(1:3) .eq. 'WAT' 
     &              .or. cmsg(3)(1:3) .eq. 'WAT') .and. 
     &              (cmsg(2)(1:3) .eq. 'air' .or. cmsg(3)(1:3) 
     &              .eq. 'air' .or. cmsg(2)(1:3) .eq. 'AIR'
     &              .or. cmsg(3)(1:3) .eq. 'AIR')) then
! Output water and air/vapor viscosity
                  vis_flag = 1
               else if ((cmsg(2)(1:3) .eq. 'wat' .or. cmsg(3)(1:3) 
     &                 .eq. 'watt' .or. cmsg(2)(1:3) .eq. 'WAT' 
     &                 .or. cmsg(3)(1:3) .eq. 'WAT') .and.
     &                 (cmsg(2)(1:3) .eq. 'co2' .or. cmsg(3)(1:3) 
     &                 .eq. 'co2' .or. cmsg(2)(1:3) .eq. 'CO2'
     &                 .or. cmsg(3)(1:3) .eq. 'CO2')) then
                  if (icarb .ne. 0) then
                     if (cmsg(2)(1:4) .eq. 'co2l' .or. cmsg(3)(1:4) 
     &                    .eq. 'co2l' .or. cmsg(2)(1:4) .eq. 'CO2l'
     &                    .or. cmsg(3)(1:4) .eq. 'CO2l') then
! Output water and CO2 liquid viscosity
                        vis_flag = 7
                     else if (cmsg(2)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(3)(1:4) .eq. 'co2g' .or. 
     &                       cmsg(2)(1:4) .eq. 'CO2g' .or. 
     &                       cmsg(3)(1:4) .eq. 'CO2g') then
! Output water and CO2 gas viscosity
                        vis_flag = 8
                     else
! Output water and CO2 viscosities
                        vis_flag = 4
                     end if
                  else
                     vis_flag = 2
                     if (iout .ne. 0) write(iout, 6067) 
     &                    'co2 viscosity'
                     if (iptty .ne. 0) write(iptty, 6067) 
     &                    'co2 viscosity'
                  end if
               else if ((cmsg(2)(1:4) .eq. 'co2l' .or. cmsg(3)(1:4) 
     &                 .eq. 'co2l' .or. cmsg(2)(1:4) .eq. 'CO2l' 
     &                 .or. cmsg(3)(1:4) .eq. 'CO2l') .and.
     &                 (cmsg(2)(1:4) .eq. 'co2g' .or. cmsg(3)(1:4) 
     &                 .eq. 'co2g' .or. cmsg(2)(1:4) .eq. 'CO2g'
     &                 .or. cmsg(3)(1:4) .eq. 'CO2g')) then
                  if (icarb .ne. 0) then
                     vis_flag = 9
                  else
                     if (iout .ne. 0) write(iout, 6067) 
     &                    'co2 viscosity'
                     if (iptty .ne. 0) write(iptty, 6067) 
     &                    'co2 viscosity'
                  end if
               end if
            case default
! Defaults determined from run parameters
               if (icarb .ne. 0) then
! Output CO2 and water viscosities
                  vis_flag = 4
               else if (ico2 .gt. 0 .or. irdof .ne. 13 .or. 
     &                 irdof .ne. 11) then
! Output water and air/vapor viscosity
                  vis_flag = 1
               else if (irdof .eq. 13) then
! Output water viscosity
                  vis_flag = 2
               else if (irdof .eq. 11) then
! Output air viscosity
                  vis_flag = 3
               end if
            end select

            if (vis_flag .ne. 0) then
               ishisv = ishis + 175
               chtmp(1) = ''
               select case (vis_flag)
               case (1)
                  chtmp(1) = '_visWAT'
                  chtmp(2) = '_visAIR'
               case(2)
                  chtmp(1) = '_visWAT'
               case(3)
                  chtmp(1) = '_visAIR'
               case(4)
                  chtmp(1) = '_visWAT'
                  chtmp(2) = '_visCO2l'
                  chtmp(3) = '_visCO2g'
               case(5)
                  chtmp(1) = '_visCO2l'
               case(6)
                  chtmp(1) = '_visCO2g'
               case(7)
                  chtmp(1) = '_visWAT'
                  chtmp(2) = '_visCO2l'
               case(8)
                  chtmp(1) = '_visWAT'
                  chtmp(2) = '_visCO2g'
                case(9)
                  chtmp(1) = '_visCO2l'
                  chtmp(2) = '_visCO2g'
               end select
               fname =  root(1:iroot) // trim(chtmp(1)) // hissfx
               open (unit=ishisv, file=fname, form='formatted')
               if (vis_flag .eq. 1 .or. vis_flag .eq. 4 .or. 
     &              vis_flag .eq. 7 .or. vis_flag .eq. 8 .or.
     &              vis_flag .eq. 9) then
                  ishisv2 = ishisv + 1
                  fname = ''
                  fname =  root(1:iroot) // trim(chtmp(2)) // hissfx
                  open (unit=ishisv2, file=fname, form='formatted')
               end if
               if (vis_flag .eq. 4) then
                  ishisv3 = ishisv2 + 1
                  fname = ''
                  fname =  root(1:iroot) // trim(chtmp(3)) // hissfx
                  open (unit=ishisv3, file=fname, form='formatted')
               end if

               select case (form_flag)
               case (0)
                  write(ishisv, 6000) verno, jdate, jtime, trim(wdd)
                  if (ishisv2 .ne. 0) write(ishisd2, 6000) verno, 
     &                 jdate, jtime, trim(wdd)
                  if (ishisv3 .ne. 0) write(ishisd3, 6000) verno, 
     &                 jdate, jtime, trim(wdd)
               case (1)
                  write(ishisv, 6005) verno, jdate, jtime
                  if (ishisv2 .ne. 0) write(ishisd2, 6005) verno, 
     &                 jdate, jtime
                  if (ishisv3 .ne. 0) write(ishisd3, 6005) verno, 
     &                 jdate, jtime
               end select
            end if

         case ('hum', 'HUM')
! Output humidity
            fname =  root(1:iroot) // '_humd' // hissfx
            ishishm = ishis + 160
            open (unit=ishishm, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishishm, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishishm, 6005) verno, jdate, jtime
            end select
         case ('zfl', 'ZFL')
! Output zone fluxes
            if (nflxz .ne. 0) then
               call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
               if (nwds .gt. 1) then
! Select fluxes to output
                  prnt_flxzvar = .FALSE.
                  do i = 2, nwds
                     if (msg(i) .eq. 3) then
                        cdum = cmsg(i)
                        select case (cdum(1:3))
                        case ('sou', 'SOU')   ! Source
                           prnt_flxzvar(1) = .TRUE.
                        case ('sin', 'SIN')   ! Sink
                           prnt_flxzvar(2) = .TRUE.
                        case ('net', 'NET')   ! Net
                           prnt_flxzvar(3) = .TRUE.
                        case ('bou', 'BOU')   ! Boundary
                           prnt_flxzvar(4) = .TRUE.
                        case ('hco', 'HCO')   ! Heat Conduction
                           prnt_flxzvar(5) = .TRUE.                
                        end select
                     end if
                  end do
               end if
               ishisfz = ishis + 500
               do i = 1, nflxz
                  if (wflux_flag) then
                     fname = ''
                     flxzone = ''
                     write (flxzone, '(i5.5)') iflxz(i)
                     fname =  root(1:iroot) // '_wflxz' // flxzone // 
     &                    hissfx
                     ishisfzz = ishisfz +i
                     open (unit=ishisfzz, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishisfzz, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishisfzz, 6005) verno, jdate, jtime
                     end select
                  end if
                  if (vflux_flag) then
                     fname = ''
                     flxzone = ''
                     write (flxzone, '(i5.5)') iflxz(i)
                     fname =  root(1:iroot) // '_vflxz' // flxzone // 
     &                    hissfx
                     ishisfzz = ishisfz + i + 400
                     open (unit=ishisfzz, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishisfzz, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishisfzz, 6005) verno, jdate, jtime
                     end select
                  end if
                  if (eflux_flag) then
                     fname = ''
                     flxzone = ''
                     write (flxzone, '(i5.5)') iflxz(i)
                     fname =  root(1:iroot) // '_eflxz' // flxzone // 
     &                    hissfx
                     ishisfzz = ishisfz + i + 800
                     open (unit=ishisfzz, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishisfzz, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishisfzz, 6005) verno, jdate, jtime
                     end select
                  end if                  
               end do
            else
               if (iout .ne. 0) write(iout, 6050)
               if (iptty .ne. 0) write(iptty, 6050)
            end if
         case ('con ', 'CON')
! Output concentrations
            if (iccen .ne. 0) then
               if (istrc .eq. 0) then
                  istrc = nufilb(9)
                  fname =  root(1:iroot) // trim(trcsfx)
                  inquire (unit = istrc, opened = opend)
                  if (.not. opend) then
                     open (unit=istrc, file=fname, form='formatted')
                     write(istrc, 6000) verno, jdate, jtime, trim(wdd)
                  end if
               end if
               ishisc = istrc + 200
! Files for output will be opened on first call to plotc1 once it is 
! determined which species will be output
            else
               if (iout .ne. 0) write(iout, 6060)
               if (iptty .ne. 0) write(iptty, 6060)
            end if               
         case ('wt ', 'WT ')
! Output water table in m
            fname =  root(1:iroot) // '_wt' // hissfx
            ishiswt = ishis + 190
            open (unit=ishiswt, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishiswt, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishiswt, 6005) verno, jdate, jtime
            end select
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .ge. 2) then
               if (msg(2) .eq. 1) wt_flag = imsg(2)
            end if
            if (out_zones) then
               ozflag = ozflag + 1
            end if
! RJP 08/09/07 added below for outputing CO2 related states
         case ('co2', 'CO2', 'sco', 'SCO')
            if (chdum(1:3) .eq. 'sco' .or. chdum(1:3) .eq. 'SCO') then
! recognize old keyword and replace with new for saturation output
               chdum = ''
               chdum = 'co2s'
            end if               
            if (icarb .ne. 0 ) then
               select case (chdum(4:4))
               case ('m', 'M')
                  select case (chdum(5:5))
                  case ('t','T')
! Output total CO2 mass in kg
                     ishiscm = ishis + 210
                  case ('f','F')
! Output free CO2 mass fraction
                     ishiscmf = ishis + 211
                  case ('d','D')
! Output dissolved CO2 mass fraction
                     ishiscmd = ishis + 212
                  case default
                     ishiscm = ishis + 210
                     ishiscmf = ishis + 211
                     ishiscmd = ishis + 212
                  end select
                  if (ishiscm .ne. 0) then
                     fname =  root(1:iroot) // '_co2mt' // hissfx
                     open (unit=ishiscm, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishiscm, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishiscm, 6005) verno, jdate, jtime
                     end select
                     if (out_zones) then
                        ozflag = ozflag + 1
                        carbflag = ozflag
                     end if
                  end if
                  if (ishiscmf .ne. 0) then
                     fname =  root(1:iroot) // '_co2mf' // hissfx
                     open (unit=ishiscmf, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishiscmf, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishiscmf, 6005) verno, jdate, jtime
                     end select
                  end if
                  if (ishiscmd .ne. 0) then
                     fname =  root(1:iroot) // '_co2md' // hissfx
                     open (unit=ishiscmd, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishiscmd, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishiscmd, 6005) verno, jdate, jtime
                     end select
                  end if
               case ('s', 'S')
! Output CO2 saturations
                  select case (chdum(5:5))
                  case ('l', 'L')
! Liquid CO2
                     ishiscsl = ishis + 230
                  case ('g', 'G')
! Gaseous CO2
                     ishiscsg = ishis + 231
                  case default
! Both liquid and gas
                     ishiscsl = ishis + 230
                     ishiscsg = ishis + 231
                  end select
                  if (ishiscsl .ne. 0) then
                     fname =  root(1:iroot) // '_co2sl' // hissfx
                     open (unit=ishiscsl, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishiscsl, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishiscsl, 6005) verno, jdate, jtime
                     end select
                  end if
                  if (ishiscsg .ne. 0) then
                     fname =  root(1:iroot) // '_co2sg' // hissfx
                     open (unit=ishiscsg, file=fname, form='formatted')
                     select case (form_flag)
                     case (0)
                        write(ishiscsg, 6000) verno, jdate, jtime, 
     &                       trim(wdd)
                     case (1)
                        write(ishiscsg, 6005) verno, jdate, jtime
                     end select
                  end if
               case default
                  if (iout .ne. 0) write(iout, 6066) trim(chdum)
                  if (iptty .ne. 0) write(iptty, 6066) trim(chdum)
               end select
            else
               if (iout .ne. 0) write(iout, 6065)
               if (iptty .ne. 0) write(iptty, 6065)               
            end if
! RJP 07/05/07 added below for outputing CO2 flux
         case ('cfl', 'CFL')
! Output CO2 zone fluxes in kg/s
            if (nflxz .ne. 0 .and. icarb .ne. 0) then
               ishiscfz = ishis + 700
               do i = 1, nflxz
                  fname = ''
                  flxzone = ''
                  write (flxzone, '(i5.5)') iflxz(i)
                  fname =  root(1:iroot) // '_co2flx' // flxzone // 
     &                 hissfx
                  ishiscfzz = ishiscfz + i
                  open (unit=ishiscfzz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishiscfzz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishiscfzz, 6005) verno, jdate, jtime
                  end select
               end do
            else
               if (icarb .eq. 0) then
                  if (iout .ne. 0) write(iout, 6065)
                  if (iptty .ne. 0) write(iptty, 6065)
               end if              
               if (nflxz .eq. 0) then
                  if (iout .ne. 0) write(iout, 6050)
                  if (iptty .ne. 0) write(iptty, 6050)
               end if
            end if
         case ('dis', 'DIS')
! Output displacements in m
            select case (chdum(4:4))
            case ('x', 'X')
               ishisdisx = ishis + 240
               fname =  root(1:iroot) // '_disx' // hissfx
               open (unit=ishisdisx, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisdisx, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisdisx, 6005) verno, jdate, jtime
               end select
            case ('y', 'Y')
               ishisdisy = ishis + 241
               fname =  root(1:iroot) // '_disy' // hissfx
               open (unit=ishisdisy, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisdisy, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisdisy, 6005) verno, jdate, jtime
               end select
             case ('z', 'Z')
               if (icnl .eq. 0 ) then
                  ishisdisz = ishis + 242
                  fname =  root(1:iroot) // '_disz' // hissfx
                  open (unit=ishisdisz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishisdisz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishisdisz, 6005) verno, jdate, jtime
                  end select
               else
                  if (iout .ne. 0) write(iout, 6070) 'Z displacement'
                  if (iptty .ne. 0) write(iptty, 6070) 'Z displacement'
               end if
            case default
! All displacements will be output
               ishisdisx = ishis + 240
               fname =  root(1:iroot) // '_disx' // hissfx
               open (unit=ishisdisx, file=fname, form='formatted')
               ishisdisy = ishis + 241
               fname =  root(1:iroot) // '_disy' // hissfx
               open (unit=ishisdisy, file=fname, form='formatted')
               if (icnl .eq. 0) then
                  ishisdisz = ishis + 242
                  fname =  root(1:iroot) // '_disz' // hissfx
                  open (unit=ishisdisz, file=fname, form='formatted')
               end if
               select case (form_flag)
               case (0)
                  write(ishisdisx, 6000) verno, jdate, jtime, trim(wdd)
                  write(ishisdisy, 6000) verno, jdate, jtime, trim(wdd)
                  if (ishisdisz .ne. 0) write(ishisdisz, 6000) 
     &                    verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisdisx, 6005) verno, jdate, jtime
                  write(ishisdisy, 6005) verno, jdate, jtime
                  if (ishisdisz .ne. 0) write(ishisdisz, 6005) 
     &                 verno, jdate, jtime
               end select
            end select
         case ('str', 'STR')
! Output stress / strain
            select case (chdum(4:5))
            case ('ai', 'AI')
               ishisstr = ishis + 250
               fname =  root(1:iroot) // '_strain' // hissfx
               open (unit=ishisstr, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisstr, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisstr, 6005) verno, jdate, jtime
               end select
            case ('x ', 'X ')
               ishisstrx = ishis + 251
               fname =  root(1:iroot) // '_strx' // hissfx
               open (unit=ishisstrx, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisstrx, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisstrx, 6005) verno, jdate, jtime
               end select
            case ('y ', 'Y ')
               ishisstry = ishis + 252
               fname =  root(1:iroot) // '_stry' // hissfx
               open (unit=ishisstry, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisstry, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisstry, 6005) verno, jdate, jtime
               end select
            case ('xy', 'XY')
               ishisstrxy = ishis + 253
               fname =  root(1:iroot) // '_strxy' // hissfx
               open (unit=ishisstrxy, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishisstrxy, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishisstrxy, 6005) verno, jdate, jtime
               end select
            case ('z ', 'Z ')
               if (icnl .eq. 0) then
                  ishisstrz = ishis + 254
                  fname =  root(1:iroot) // '_strz' // hissfx
                  open (unit=ishisstrz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishisstrz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishisstrz, 6005) verno, jdate, jtime
                  end select
               else
                  if (iout .ne. 0) write(iout, 6070) 'Z stress'
                  if (iptty .ne. 0) write(iptty, 6070) 'Z stress'
               end if
            case ('xz', 'XZ')
               if (icnl .eq. 0) then
                  ishisstrxz = ishis + 255
                  fname =  root(1:iroot) // '_strxz' // hissfx
                  open (unit=ishisstrxz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishisstrxz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishisstrxz, 6005) verno, jdate, jtime
                  end select
               else
                  if (iout .ne. 0) write(iout, 6070) 'XZ stress'
                  if (iptty .ne. 0) write(iptty, 6070) 'XZ stress'
               end if
            case ('yz', 'YZ')
               if (icnl .eq. 0) then
                  ishisstryz = ishis + 256
                  fname =  root(1:iroot) // '_stryz' // hissfx
                  open (unit=ishisstryz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishisstryz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishisstryz, 6005) verno, jdate, jtime
                  end select
               else
                  if (iout .ne. 0) write(iout, 6070) 'YZ stress'
                  if (iptty .ne. 0) write(iptty, 6070) 'YZ stress'
               end if
            case default
! All stresses will be output
               ishisstrx = ishis + 251
               fname =  root(1:iroot) // '_strx' // hissfx
               open (unit=ishisstrx, file=fname, form='formatted')
               ishisstry = ishis + 252
               fname =  root(1:iroot) // '_stry' // hissfx
               open (unit=ishisstry, file=fname, form='formatted')
               ishisstrxy = ishis + 253
               fname =  root(1:iroot) // '_strxy' // hissfx
               open (unit=ishisstrxy, file=fname, form='formatted')
               if (icnl .eq. 0) then
                  ishisstrz = ishis + 254
                  fname =  root(1:iroot) // '_strz' // hissfx
                  open (unit=ishisstrz, file=fname, form='formatted')
                  ishisstrxz = ishis + 255
                  fname =  root(1:iroot) // '_strxz' // hissfx
                  open (unit=ishisstrxz, file=fname, form='formatted')
                  ishisstryz = ishis + 256
                  fname =  root(1:iroot) // '_stryz' // hissfx
                  open (unit=ishisstryz, file=fname, form='formatted')
               end if
               select case (form_flag)
               case (0)
                  write(ishisstrx, 6000) verno, jdate, jtime, trim(wdd)
                  write(ishisstry, 6000) verno, jdate, jtime, trim(wdd)
                  write(ishisstrxy, 6000) verno, jdate, jtime, trim(wdd)
                  if (icnl .eq. 0) then
                     write(ishisstrz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                     write(ishisstrxz, 6000) verno, jdate, jtime,  
     &                    trim(wdd)
                     write(ishisstryz, 6000) verno, jdate, jtime,  
     &                    trim(wdd)
                  end if
               case (1)
                  write(ishisstrx, 6005) verno, jdate, jtime
                  write(ishisstry, 6005) verno, jdate, jtime
                  write(ishisstrxy, 6000) verno, jdate, jtime
                  if (icnl .eq. 0) then
                     write(ishisstrz, 6005) verno, jdate, jtime
                     write(ishisstrxz, 6005) verno, jdate, jtime
                     write(ishisstryz, 6005) verno, jdate, jtime
                  end if
               end select
            end select
         case ('rel', 'REL')
! Output table of relative permeability and cappillary pressure values
            fname = root(1:iroot) // '_rlp'
c            select case (form_flag)
c            case (0)
c                fname = trim(fname) // '.tbl'
c            case (1)
c               fname = trim(fname) // '_tbl.dat'
c            case (2)
c               fname = trim(fname) // '_tbl.csv'
c            end select 
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
! Set defaults
            num_sat = 0
            delta_sat = 0.05
            if (nwds .gt. 1) then
               if ( msg(2) .eq. 2) then
! Read spacing for saturation output
                  delta_sat = xmsg(2)
               else if ( msg(2) .eq. 3 .and. cmsg(2) .eq. 'list') then
! Read a list of saturation values for output
                 read (inpt, *) num_sat
                 allocate (sat_out(num_sat))
                 backspace inpt
                 read (inpt, *) num_sat, (sat_out(i), i = 1, num_sat)
               end if
            end if
            ishisrlp = ishis + 500
            open (unit=ishisrlp, file=fname, form='formatted')
c gaz 071623 first line now starts in check_rlp
c            select case (form_flag)
c            case (0)
c               write(ishisrlp, 6000) verno, jdate, jtime, trim(wdd)
c            case (1)
c               write(ishisrlp, 6005) verno, jdate, jtime
c            end select            
         case ('glo', 'GLO')
! Output global variables
            fname =  root(1:iroot) // '_glob' // hissfx
            ishisg = ishis + 180
            open (unit=ishisg, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishisg, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishisg, 6005) verno, jdate, jtime
            end select
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .eq. 1) then
! Set defaults
               if (ico2.ge.0 .or. ice.ne.0) then
! Global Mass & Energy Balances
                  glob_flag = 1
               else if (ico2.lt.0 .and. ice.eq.0) then
! Global Water & Air Balances
                  glob_flag = 2
               else if (ico2.lt.0 .and. ice.ne.0) then
                  glob_flag = 7
               end if
            else if (nwds .eq. 2) then
               globs: select case (cmsg(2)(1:3))
               case ('mas', 'wat', 'MAS', 'WAT')
! Output mass or water balances (excluding steam)
                  glob_flag = 3
               case ('ste', 'STE')
! Output mass or water balances (including steam)
                  glob_flag = 4
               case ('air', 'AIR')
! Output air or vapor balances
                  glob_flag = 5
               case ('ene', 'ENE')
! Output energy balances
                  glob_flag = 6
               end select globs
            else
            end if
! Total water in system kg
! Total mass of steam in system kg
! Water discharge kg (kg/s)
! Water input kg (kg/s)
! Net kg water discharge
! Total water discharge kg (kg/s)
! Total water input kg (kg/s)
! Total air in system kg
! Air discharge kg (kg/s)
! Air input kg (kg/s)
! Net kg air discharge
! Total air discharge kg (kg/s)
! Total air input kg (kg/s)
! Total enthalpy in system MJ
! Enthalpy discharge MJ (MJ/s)
! Enthalpy input MJ (MJ/s)
! Net MJ enthalpy discharge
! Total enthalpy discharge MJ (MJ/s)
! Total enthalpy input MJ (MJ/s)

         end select flags
      enddo

 6000 format(a30, 2x, a10, 2x, a8, /, a)
 6005 format('TITLE = "', a30, 2x, a10, 2x, a8,'"')
 6010 format('Liquid only problem, saturation will not be output')
 6020 format('Head not specified for the problem, will not be output')
 6025 format('Head already selected, output will be in ', a6)
 6030 format('Head problem, only total/water pressure will be output')
 6040 format('Isothermal air-water problem, only water and air ', 
     &     'pressure will be output')
 6050 format('Zone fluxes not specified for the problem, ', 
     &     'will not be output')
 6060 format('Transport not specified for the problem, ',      
     &     'will not be output')
 6065 format('CO2 not specified for the problem, ',      
     &     'will not be output')
 6066 format('Invalid CO2 output option, ', a,      
     &     'no output will be generated')
 6067 format('CO2 not specified for the problem, ', a,      
     &     ' will not be output')
 6070 format('2-D problem, ', a, ' will not be output')
 6080 format('ihms =', i3, ' Time history for ', a, 
     & ' will not be output')

      end
