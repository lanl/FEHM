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
      use comdi
      use comsi, only : ihms
      use comxi
      use comzone
      use comwt, only : wt_flag
      use davidi
      implicit none

      character*5 flxzone
      character*8 hissfx, trcsfx
      character*32 cmsg(4), cdum
      character*80 chdum
      character*120 fname, root
      logical null1, opend
      integer i, iroot, msg(4), imsg(4), nwds, ishisfzz, ishiscfzz
      integer ishiscs
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
      ishisp = 0
      ishist = 0
      ishishd = 0
      ishiss = 0
      ishisf = 0
      ishise = 0
      ishishm = 0
      ishisfz = 0
      ishiswt = 0
      ishisc = 0
      ishiswc = 0
      ishiscm = 0
      ishiscfz = 0
      ishiscs = 0
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
         case ('yea','day','sec','hrs','YEA','DAY','SEC','HRS')
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
            else if (chdum(1:3) .eq. 'hrs' .or. chdum(1:3) .eq. 'HRS')  
     &              then
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
            fname =  root(1:iroot) // '_pres' // hissfx
            ishisp = ishis + 100
            open (unit=ishisp, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishisp, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishisp, 6005) verno, jdate, jtime
            end select
            if (out_zones)  then
               ozflag = ozflag + 1
               pflag = ozflag
            end if
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .gt. 1) then
! Determine which pressures to output
               nwd : select case (nwds)
               case (2)
                  if (cmsg(2)(1:3) .eq. 'tot' .or. 
     &                 cmsg(2)(1:3) .eq. 'TOT' ) then
! Output total/water pressure
                     pres_flag = 1
                  else if ((cmsg(2)(1:3) .eq. 'air' .or. cmsg(2)(1:3) 
     &                .eq. 'AIR') .and. ihead .ne. 1) then
! Output air/vapor pressure
                     pres_flag = 4
                  else if ((cmsg(2)(1:3) .eq. 'cap' .or. cmsg(2)(1:3) 
     &                    .eq. 'CAP') .and. ico2 .gt. 0) then
! Output capillary pressure
                     pres_flag = 5
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
     &                    ico2 .gt. 0) then
                     pres_flag = 6
                  else if ((cmsg(2)(1:3) .eq. 'air' .or. cmsg(3)(1:3) 
     &                    .eq. 'air' .or. cmsg(2)(1:3) .eq. 'AIR'
     &                    .or. cmsg(3)(1:3) .eq. 'AIR') .and. 
     &                    (cmsg(2)(1:3) .eq. 'cap'  .or. cmsg(3)(1:3) 
     &                    .eq. 'cap' .or. cmsg(2)(1:3) .eq. 'CAP'
     &                    .or. cmsg(3)(1:3) .eq. 'CAP') .and. 
     &                    ico2 .gt. 0)  then
                     pres_flag = 7
                  else
                     pres_flag = 2
                     if (iout .ne. 0) write(iout, 6040)
                     if (iptty .ne. 0) write(iptty, 6040)
                  end if
               case (4)
                  if (ico2 .gt. 0) then
                     pres_flag = 3
                  else
                     if (ihead .eq. 1) then
                        pres_flag = 1
                        if (iout .ne. 0) write(iout, 6030)
                        if (iptty .ne. 0) write(iptty, 6030)
                     else
                        pres_flag = 2
                        if (iout .ne. 0) write(iout, 6040)
                        if (iptty .ne. 0) write(iptty, 6040)
                     end if
                  end if
               end select nwd
            else
! Defaults determined from run parameters
               if (ico2 .gt. 0) then
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
            if (ihead .ne. 0 ) then
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
            if (ihead .ne. 0 ) then
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
            fname =  root(1:iroot) // '_flow' // hissfx
            ishise = ishis + 150
            open (unit=ishise, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishise, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishise, 6005) verno, jdate, jtime
            end select
            if (out_zones) then
               ozflag = ozflag + 1
               eflag = ozflag
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
                        case ('vap', 'VAP')   ! Vapor
                           prnt_flxzvar(5) = .TRUE.
                        end select
                     end if
                  end do
               end if
               ishisfz = ishis + 500
               do i = 1, nflxz
                  fname = ''
                  flxzone = ''
                  write (flxzone, '(i5.5)') iflxz(i)
                  fname =  root(1:iroot) // '_flxz' // flxzone // hissfx
                  ishisfzz = ishisfz +i
                  open (unit=ishisfzz, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishisfzz, 6000) verno, jdate, jtime, 
     &                    trim(wdd)
                  case (1)
                     write(ishisfzz, 6005) verno, jdate, jtime
                  end select
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
! RJP added below for outputing CO2 mass
         case ('co2', 'CO2')
! Output CO2 mass in kg
            fname =  root(1:iroot) // '_co2m' // hissfx
            ishiscm = ishis + 210
            open (unit=ishiscm, file=fname, form='formatted')
            select case (form_flag)
            case (0)
               write(ishiscm, 6000) verno, jdate, jtime, trim(wdd)
            case (1)
               write(ishiscm, 6005) verno, jdate, jtime
            end select
            if (out_zones) then
               ozflag = ozflag + 1
               carbflag = ozflag
            end if
! RJP 07/05/07 added below for outputing CO2 flux
         case ('cfl', 'CFL')
! Output CO2 zone fluxes in kg/s
            if (nflxz .ne. 0) then
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
               if (iout .ne. 0) write(iout, 6050)
               if (iptty .ne. 0) write(iptty, 6050)
            end if
         case ('sco', 'SCO')
! Output CO2 saturations
            select case (chdum(4:4))
            case ('l', 'L')
! Liquid CO2
               ishiscsl = ishis + 230
               fname =  root(1:iroot) // '_co2sl' // hissfx
               open (unit=ishiscsl, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishiscsl, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishiscsl, 6005) verno, jdate, jtime
               end select
c               if (out_zones) then
c                  ozflag = ozflag + 1
c                  carbflag2 = ozflag
c               end if
            case ('g', 'G')
! Gaseous CO2
               ishiscs = ishis + 231
               fname =  root(1:iroot) // '_co2sg' // hissfx
               open (unit=ishiscsg, file=fname, form='formatted')
               select case (form_flag)
               case (0)
                  write(ishiscsg, 6000) verno, jdate, jtime, trim(wdd)
               case (1)
                  write(ishiscsg, 6005) verno, jdate, jtime
               end select
c               if (out_zones) then
c                  ozflag = ozflag + 1
c                  carbflag2 = ozflag
c               end if
            case default
! Both liquid and gas
               do i = 0, 1
                  select case (i)
                  case (0)
                     fname =  root(1:iroot) // '_co2sl' // hissfx
                     ishiscsl = ishis + 230
                     ishiscs = ishiscsl
                  case (1)
                     fname =  root(1:iroot) // '_co2sg' // hissfx
                     ishiscsg = ishis + 231
                     ishiscs = ishiscsg
                  end select
                  open (unit=ishiscs, file=fname, form='formatted')
                  select case (form_flag)
                  case (0)
                     write(ishiscs, 6000) verno, jdate, jtime, trim(wdd)
                  case (1)
                     write(ishiscs, 6005) verno, jdate, jtime
                  end select
c               if (out_zones) then
c                  ozflag = ozflag + 1
c                  carbflag2 = ozflag
c               end if
               end do
            end select
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
 6060 format('Transport not specified for the problem,',      
     &     'will not be output')
 6070 format('2-D problem, ', a, ' will not be output')
 6080 format('ihms =', i3, ' Time history for ', a, 
     & ' will not be output')

      end
