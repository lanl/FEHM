      subroutine plot_new(igf, totalflin, totalein, curinflow,
     &     cureinflow)
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
!D1 Write time history files.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
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
      use combi
      use comci, only : enlf, envf, rolf, dil, rlf, rovf, div, rvf
      use comco2, only : icarb, denco2h, fl, fg, fw, phico2, xc, yc,
     &     co2_prop, wat_prop
      use comdi
      use comdti
      use comfi
      use comii, only : crl
      use comrlp, only : rlpnew
      use comsi, only : du, dv , dw, str_x, str_y, str_z, str_xy,
     &  str_xz, str_yz, vol_strain, du_ini, dv_ini, dw_ini, idisp_rel
      use comwt
      use comzone
      use davidi
      implicit none

      integer i, ic1, ic2, igf, mi, length, last_step, j, k, ndummy
      integer count, start, end, var_num, ishisfzz, ishiscfzz, ishisnode
      integer, allocatable :: dumlist(:)
      real*8 pwatersat, psat, dpdummy, headdum, ptime, last_time
      real*8 totalflin, totalein, curinflow, cureinflow
      real*8 dae, dmass, start_ae, start_mass
      real*8 pdum, tdum, rolconv, dumconv, dumconv1
      real*8, allocatable :: dumv(:)
      character*400 glob_string
      character*80 title_string, formh_string 
      character*80 form1_string, form2_string, formp_string
!      character*110 formz_string
      character*120 formz_string
      character*500 info_string
      character*90 :: formf_string
      character*28 vt_string
      character*80 formcs_string
      character*14 time_string, file_format, dumv_string
      character*18, allocatable :: var_string(:), var_tmp(:)
      character*10 zone_string
      character*3 dls
      logical time2print

      save last_step, last_time, dumv, start_ae, start_mass
      save form1_string, form2_string, formh_string
      save formp_string, formcs_string

      ishisnode = ishisp + ishist + ishishd + ishiss + ishisf + ishisfa 
     &     + ishise+ ishisef + ishisd + ishisv + ishishm + ishiswt 
     &     + ishisc + ishiswc + ishiscm + ishiscmf + ishiscmd + ishiscsl 
     &     + ishiscsg + ishisdisx + ishisdisy + ishisdisz
     &     + ishisstr + ishisstrx + ishisstry + ishisstrz
      if (m .eq. 0 .and. node_azones .eq. 0) then
         if (ishisfz .eq. 0) then
            write (ierr, 100)
            if (iout .ne. 0) write (iout, 100)
            if (iptty .ne. 0 ) write (iptty, 100)
            stop
         else if (ishisnode .ne. 0) then
            if (l .eq. 1) write (ierr, 120)
            if (iout .ne. 0) write (iout, 120)
            if (iptty .ne. 0 ) write (iptty, 120)
         end if
      end if
      if (out_zones .and. ozflag .eq. 0) then
         write (ierr, 110)
         if (iout .ne. 0) write (iout, 110)
         if (iptty .ne. 0 ) write (iptty, 110)
         stop
      endif
 100  format ('***** STOPPING *****', /, 
     &     'nodes and/or zones must be defined for history plot output')
 110  format ('***** STOPPING *****', /, 
     &     'node macro for zone averaging must precede hist macro')
 120  format ('***** WARNING *****', /, 
     &     'nodes and/or zones must be defined for history plot output',
     &     ' other than zone fluxes')

            
      if (igf .eq. 0)  then

         if(ico2.gt.0) then
            write(ishis, '(a4)')  'ngas'
         elseif(ico2.lt.0) then
            write(ishis, '(a4)')  'airw'
         else
c            write(ishis, '(a4)')  '    '
         end if

         if (iccen .ne. 0)  then
            write(ishis, '(a4)')  'trac'
         else
c            write(ishis, '(a4)')  '    '
         end if

         if ( istrs .ne. 0 )  then
            write(ishis, '(a4)')  'strs'
         else
c            write(ishis, '(a4)')  '    '
         end if

         if (time_flag .eq. 1) then
            time_string = 'Time (years)'
         else if (time_flag .eq. 2) then
            time_string = ' Time (days)'
         else if (time_flag .eq. 3) then
            time_string = 'Time(seconds)'
         else if (time_flag .eq. 4) then
            time_string = 'Time (hours)'
         end if

         info_string = ''
         form1_string = ''
         form2_string = ''
         formh_string = ''
         formp_string = ''
         formz_string = ''
         title_string = ''
         length = 0
         ic1 = 1
         ic2 = 1
         if (out_zones) then
            if (.not. allocated(dumlist)) 
     &           allocate (dumlist(node_azones))
            zone_ptr => zone_head
            do i = 1, node_azones
               dumlist(i) = zone_ptr%node_number
               zone_ptr =>zone_ptr%nnp
            end do
         end if

         var_num = max (m + node_azones, 1)
         allocate (var_string(var_num), var_tmp(var_num))
c gaz 081023 adjust spacing for easier reading
         do i = 1, m
            dumv_string = ''
            write (dumv_string(1:8), '(i8)') nskw(i)
            var_string(i) = '    Node' // dumv_string(1:8)
            continue
         end do
         if (out_zones) then
            do i = m+1, var_num
               dumv_string = ''
               write (dumv_string, '(i8)') dumlist(i-m)
               var_string(i) = 'Zone' // trim(dumv_string)
            end do
         end if
         var_tmp = var_string

         if (form_flag .eq. 1) then
            file_format = 'Tecplot'
            write (form1_string, 200) max (m, 1)
            write (form2_string, 200) var_num
         else if (form_flag .eq. 2) then
            file_format = 'Surfer (csv)'
            write (form1_string, 210) max (m, 1)
            write (form2_string, 210) var_num
         else
            file_format = 'Standard text'
            write (form1_string, 220) max(m, 1)
            write (form2_string, 220) var_num
         end if
 200     format ("('variables = ',", "'", '"', "', a, '", '" ', "',", 
     &        i5, "('", '"', "', a, '", '" ', "'))")
 210     format ("( a,", i5, '(", ", a))')
 220     format ("( a,", i5, "(1x, a))")
 230     format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')

         if (ishisp .ne. 0 ) then
! Output pressures in MPa
            select case (pres_flag)
            case (1)
               title_string = 'Water Pressure (MPa)'
               call plot_header(ishisp,var_num,form2_string)
               info_string = info_string(ic1:ic2) // 'Water pressure'
            case (2)
               title_string = 'Water Pressure (MPa)'
               call plot_header(ishisp,var_num,form2_string)
               title_string = 'Air pressure (MPa)'
               call plot_header(ishisp2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and air pressure'
           case (3)
               title_string = 'Water Pressure (MPa)'
               call plot_header(ishisp,var_num,form2_string)
               title_string = 'Vapor pressure (MPa)'
               call plot_header(ishisp2,m,form1_string)
               title_string = 'Capillary pressure (MPa)'
               call plot_header(ishisp3,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water, vapor, and capillary pressure'
            case (4)
               title_string = 'Air Pressure (MPa)'
               call plot_header(ishisp,m,form1_string)
                info_string = info_string(ic1:ic2) // 'Air pressure'
           case (5)
               title_string = 'Capillary Pressure (MPa)'
               call plot_header(ishisp,m,form1_string)
                info_string = info_string(ic1:ic2) // 
     &              'Capillary pressure'
            case (6)
               title_string = 'Total Pressure (MPa)'
               call plot_header(ishisp,var_num,form2_string)
               title_string = 'Capillary pressure (MPa)'
               call plot_header(ishisp2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Total and capillary pressure'
            case (7)
               title_string = 'Vapor Pressure (MPa)'
               call plot_header(ishisp,m,form1_string)
               title_string = 'Capillary pressure (MPa)'
               call plot_header(ishisp2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Vapor and capillary pressure'
            case (8)
               title_string = 'CO2 Pressure (MPa)'
               call plot_header(ishisp,m,form1_string)
               info_string = info_string(ic1:ic2) // 'CO2 pressure'
            case (9)
               title_string = 'CO2 Pressure (MPa)'
               call plot_header(ishisp,m,form1_string)
               title_string = 'Water pressure (MPa)'
               call plot_header(ishisp2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 pressure'
            end select
            ic2 = len_trim(info_string) + 1
         end if
         if (ishist .ne. 0 ) then
! Output temperature in degrees C
            info_string = info_string(ic1:ic2) // 'Temperature '
            ic2 = len_trim(info_string) + 1
            title_string = 'Temperature (C)'
            call plot_header(ishist,var_num,form2_string)
         end if
         if (ishishd .ne. 0 ) then
            info_string = info_string(ic1:ic2) // 'Head '
            ic2 = len_trim(info_string) + 1
            if (ishishd .eq. ishis + 120) then
! Output heads in m
               title_string = 'Head (m)' 
               call plot_header(ishishd,var_num,form2_string)
            else if (ishishd .eq. ishis + 125) then
! Output heads in ft
               title_string = 'Head (ft)'
               call plot_header(ishishd,var_num,form2_string)
            else if (ishishd .eq. ishis + 121) then
               title_string = 'Head (m)'
               if (out_zones) then
                  dumv_string = 'Node/Zone'
               else
                  dumv_string = 'Node'
               end if
               if (form_flag .eq. 1) then
                  formh_string = "variables = " // '"' // 
     &                 trim(time_string) // '" "' // trim(dumv_string) 
     &                 // '" "Head(m)"'
                  write (ishishd, '(a)') trim(formh_string)
                  write(ishishd,230) 50., 95., trim(wdd)
                  formh_string = "(g16.9, 1x, i7, 1x, g16.9)"
               else if (form_flag .eq. 2) then
                  formh_string = trim(time_string) // ", " // 
     &                 trim(dumv_string) // ", Head(m)"
                  write (ishishd, '(a)') trim(formh_string)
                  formh_string = '(g16.9, ", ", i7, ", ", g16.9)'
               else
                  formh_string = trim(time_string) // " " //
     &                trim(dumv_string) // " Head(m)" 
                  write (ishishd, '(a)') trim(formh_string)
                  formh_string = "(g16.9, 1x, i7, 1x, g16.9)"
               end if
            end if
            if (.not. allocated(dumv)) then
               if (m .gt. 0) then
                  allocate (dumv(m))
               else
                  allocate (dumv(1))
               end if
            end if
         end if
         if (ishiss .ne. 0 ) then
! Ouput saturations
            info_string = info_string(ic1:ic2) // 'Saturation '
            ic2 = len_trim(info_string) + 1
            title_string = 'Saturation'
            call plot_header(ishiss, m, form1_string)
            continue
         end if
	   if (ishiswc .ne. 0 ) then
! Ouput water content
            info_string = info_string(ic1:ic2) // 'Water_content '
            ic2 = len_trim(info_string) + 1
            title_string = 'Water_content'
            call plot_header(ishiswc, m, form1_string)
         end if
         if (ishisfa .ne. 0 ) then
! Ouput flow in kg/s
            info_string = info_string(ic1:ic2) // 'floa '
            ic2 = len_trim(info_string) + 1
            title_string = 'Floa (kg/s)'
            call plot_header(ishisfa, m, form1_string)
            if (.not. allocated(dumv)) then
               if (m .gt. 0) then
                  allocate (dumv(m))
               else
                  allocate (dumv(1))
               end if
            end if
         end if
         if (ishisf .ne. 0 ) then
! Ouput flow in kg/s
            info_string = info_string(ic1:ic2) // 'flow '
            ic2 = len_trim(info_string) + 1
            title_string = 'Flow (kg/s)'
            call plot_header(ishisf, m, form1_string)
            if (.not. allocated(dumv)) then
               if (m .gt. 0) then
                  allocate (dumv(m))
               else
                  allocate (dumv(1))
               end if
            end if
         end if
         if (ishise .ne. 0 ) then
! Ouput enthalpy in MJ/kg
            info_string = info_string(ic1:ic2) // 'Enthalpy '
            ic2 = len_trim(info_string) + 1
            title_string = 'Enthalpy (MJ/kg)'
            call plot_header(ishise,var_num,form2_string)
         end if
         if (ishisef .ne. 0 ) then
! Ouput enthalpy in MJ/s
            info_string = info_string(ic1:ic2) // 'Enthalpy Flow'
            ic2 = len_trim(info_string) + 1
            title_string = 'Enthalpy Flow (MJ/s)'
            call plot_header(ishise,var_num,form2_string)
         end if
         if (ishisd .ne. 0) then
! Output density in kg/m^3
            select case (den_flag)
            case (1)
               title_string = 'Water Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               title_string = 'Air Density (kg/m^3)'
               call plot_header(ishisd2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and air density'
            case (2)
               title_string = 'Water Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water density'
            case (3)
               title_string = 'Air Density (kg/m^3)'
               call plot_header(ishisd,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Air density'
            case (4)
               title_string = 'Water Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               title_string = 'CO2 liquid Density (kg/m^3)'
               call plot_header(ishisd2,m,form1_string)
               title_string = 'CO2 gas Density (kg/m^3)'
               call plot_header(ishisd3,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 liquid and gas density'
            case (5)
               title_string = 'CO2 liquid Density (kg/m^3)'
               call plot_header(ishisd,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 liquid density'
            case (6)
               title_string = 'CO2 gas Density (kg/m^3)'
               call plot_header(ishisd,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 gas density'
            case (7)
               title_string = 'Water Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               title_string = 'CO2 liquid Density (kg/m^3)'
               call plot_header(ishisd2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 liquid density'
             case (8)
               title_string = 'Water Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               title_string = 'CO2 gas Density (kg/m^3)'
               call plot_header(ishisd2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 gas density'
              case (9)
               title_string = 'CO2 liquid Density (kg/m^3)'
               call plot_header(ishisd,var_num,form2_string)
               title_string = 'CO2 gas Density (kg/m^3)'
               call plot_header(ishisd2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 liquid and gas density'
           end select
            ic2 = len_trim(info_string) + 1
         end if
         if (ishisv .ne. 0) then
! Output viscosity in Pa s
            select case (vis_flag)
            case (1)
               title_string = 'Water Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               title_string = 'Air Viscosity (Pa s)'
               call plot_header(ishisv2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and air viscosity'
            case (2)
               title_string = 'Water Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water viscosity'
            case (3)
               title_string = 'Air Viscosity (Pa s)'
               call plot_header(ishisv,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Air viscosity'
            case (4)
               title_string = 'Water Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               title_string = 'CO2 liquid Viscosity (Pa s)'
               call plot_header(ishisv2,m,form1_string)
               title_string = 'CO2 gas Viscosity (Pa s)'
               call plot_header(ishisv3,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 liquid and gas viscosity'
            case (5)
               title_string = 'CO2 liquid Viscosity (Pa s)'
               call plot_header(ishisv,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 liquid viscosity'
            case (6)
               title_string = 'CO2 gas Viscosity (Pa s)'
               call plot_header(ishisv,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 gas viscosity'
            case (7)
               title_string = 'Water Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               title_string = 'CO2 liquid viscosity (Pa s)'
               call plot_header(ishisv2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 liquid viscosity'
             case (8)
               title_string = 'Water Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               title_string = 'CO2 gas Viscosity (Pa s)'
               call plot_header(ishisv2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'Water and CO2 gas viscosity'
              case (9)
               title_string = 'CO2 liquid Viscosity (Pa s)'
               call plot_header(ishisv,var_num,form2_string)
               title_string = 'CO2 gas Viscosity (Pa s)'
               call plot_header(ishisv2,m,form1_string)
               info_string = info_string(ic1:ic2) // 
     &              'CO2 liquid and gas viscosity'
           end select
           ic2 = len_trim(info_string) + 1
         end if
         if (ishishm .ne. 0 ) then
! Ouput humidity
            info_string = info_string(ic1:ic2) // 'Humidity '
            ic2 = len_trim(info_string) + 1
            title_string = 'Humidity'
            call plot_header(ishishm, m, form1_string)
            if (.not. allocated(dumv)) then
               if (m .gt. 0) then
                  allocate (dumv(m))
               else
                  allocate (dumv(1))
               end if
            end if
         end if
c RJP 04/30/07 added following for outputting time-dependent CO2 mass
         if (ishiscm .ne. 0 ) then
! Output CO2 mass kg
            info_string = info_string(ic1:ic2) // 'CO2 mass '
            ic2 = len_trim(info_string) + 1
            title_string = 'CO2 mass (Kg)'
            call plot_header(ishiscm, var_num, form2_string)
         end if
         if (ishiscmf .ne. 0 ) then
! Output free CO2 mass fraction
            info_string = info_string(ic1:ic2) // 
     &           'Free CO2 mass fraction'
            ic2 = len_trim(info_string) + 1
            title_string = 'Free CO2 mass fraction'
            call plot_header(ishiscmf, var_num, form2_string)
         end if
         if (ishiscmd .ne. 0 ) then
! Output dissolved CO2 mass fraction
            info_string = info_string(ic1:ic2) // 
     &           'Dissolved CO2 mass fraction'
            ic2 = len_trim(info_string) + 1
            title_string = 'Dissolved CO2 mass fraction'
            call plot_header(ishiscmd, var_num, form2_string)
         end if
         if(icarb.ne.0) then
            if (ishiscsl .ne. 0 ) then
! Output CO2 Liquid Saturation
               info_string = info_string(ic1:ic2) // 
     &              'CO2 liquid saturation '
               ic2 = len_trim(info_string) + 1
               title_string = 'CO2 Liquid Saturation'
               call plot_header(ishiscsl, m, form1_string)
            end if
            if (ishiscsg .ne. 0 ) then
! Output CO2 Gas Saturation
               info_string = info_string(ic1:ic2) // 
     &              'CO2 gaseous saturation '
               ic2 = len_trim(info_string) + 1
               title_string = 'CO2 Gaseous Saturation'
               call plot_header(ishiscsg, m, form1_string)
            end if
         endif
         if (ishisfz .ne. 0) then
! Output zone fluxes
 240        format ('(a, ', i3, '(a))')
            info_string = info_string(ic1:ic2) // 'Zone flux '
            ic2 = len_trim(info_string) + 1
            formz_string = ''
            if (form_flag .eq. 1) then
               dls = '" "'
               vt_string ='variables = "' // trim(time_string)
            else if (form_flag .eq. 2) then
               dls = ", "
               vt_string = trim(time_string)
            else
               dls = " "
               vt_string = trim(time_string)
            end if
            formf_string = ''
            if (prnt_flxzvar(1)) formf_string = 
     &           trim(formf_string) // dls //  'Source'
            if (prnt_flxzvar(2)) formf_string = 
     &           trim(formf_string) // dls // 'Sink'
            if (prnt_flxzvar(3)) formf_string = 
     &           trim(formf_string) // dls // 'Net'
            if (prnt_flxzvar(4)) formf_string = 
     &           trim(formf_string) // dls // 'Boundary'
            if (form_flag .eq. 1) 
     &           formf_string = trim(formf_string) // '"'
            write(formz_string, 240) nflxz
            do i = 1, nflxz
               if (wflux_flag) then
                  ishisfzz = ishisfz + i
                  write(zone_string, '("Zone ", i5.5, 1x)') iflxz(i)
                  title_string = trim (zone_string) // 
     &                 ' Water Flux (kg/s)'
                  if (form_flag .eq. 0)
     &                 write(ishisfzz, '(a)') trim(title_string)
                  if (form_flag .eq. 2) then
                     write(ishisfzz, '(a, ": ",  a, a)') 
     &                    trim(zone_string), 
     &                    trim(vt_string), trim(formf_string)
                  else
                     write(ishisfzz, formz_string) trim(vt_string), 
     &                    trim(formf_string)
                  end if
                  if (form_flag .eq. 1) then
                     write(ishisfzz, 230) 50., 95., trim(wdd)
                     write(ishisfzz, 230) 50., 90., trim(title_string)
                  endif
               end if
               if (vflux_flag) then
                  ishisfzz = ishisfz + i + 400
                  write(zone_string, '("Zone ", i5.5, 1x)') iflxz(i)
                  title_string = trim (zone_string) // 
     &                 ' Vapor Flux (kg/s)'
                  if (form_flag .eq. 0)
     &                 write(ishisfzz, '(a)') trim(title_string)
                  if (form_flag .eq. 2) then
                     write(ishisfzz, '(a, ": ", a, a)') 
     &                    trim(zone_string), 
     &                    trim(vt_string), trim(formf_string)
                  else
                     write(ishisfzz, formz_string) trim(vt_string), 
     &                    trim(formf_string)
                  end if
                  if (form_flag .eq. 1) then
                     write(ishisfzz, 230) 50., 95., trim(wdd)
                     write(ishisfzz, 230) 50., 90., trim(title_string)
                  endif
               end if
            end do
         end if
         if (ishiscfz .ne. 0) then
! Output CO2 zone fluxes
            info_string = info_string(ic1:ic2) // 'Zone CO2 flux '
            ic2 = len_trim(info_string) + 1
            formz_string = ''
            do i = 1, nflxz
               ishiscfzz = ishiscfz + i
               zone_string = ''
               write(zone_string, '("Zone ", i5.5, 1x)') iflxz(i)
               title_string = trim(zone_string) // ' CO2 Flux (kg/s)'
               if (form_flag .eq. 1) then
                  formz_string = 'variables = "' // trim(time_string) //
     &                 '"' // ' "Source" "Sink" "Net In" "Net Out" ' //
     &                 '"Boundary" "SourceG" "SinkG" "NetGIn" "NetGOut"'
                  write(ishiscfzz, '(a)') trim(formz_string)
                  write(ishiscfzz, 230) 50., 95., trim(wdd)
                  write(ishiscfzz, 230) 50., 90., trim(title_string)
               else if (form_flag .eq. 2) then
                  formz_string = trim(zone_string) // ", " //
     &                 trim(time_string) // 
     &                 ", Source, Sink, Net In, Net Out, Boundary, " //
     &                 ", SourceG, SinkG, NetGIn NetGOut"
                  write(ishiscfzz, '(a)') trim(formz_string)
               else
                  formz_string = trim(time_string) //
     &                 " Source Sink NetIn NetOut Boundary SourceG" //
     &				 " SinkG NetGIn NetGOut"
                  write(ishiscfzz, '(a)') trim(title_string)
                  write(ishiscfzz, '(a)') trim(formz_string)
               end if
               
            end do
         end if
         if (ishisc .ne. 0) then
! Output concentrations
            info_string = info_string(ic1:ic2) // 'Concentration '
            ic2 = len_trim(info_string) + 1
         end if
         if (ishisdisx .ne. 0) then
! Output x displacement
            info_string = info_string(ic1:ic2) // 'X displacement '
            ic2 = len_trim(info_string) + 1
            title_string = 'X Displacement (m)'
            call plot_header(ishisdisx,var_num,form2_string)
         end if
         if (ishisdisy .ne. 0) then
! Output y displacement
            info_string = info_string(ic1:ic2) // 'Y displacement '
            ic2 = len_trim(info_string) + 1
            title_string = 'Y Displacement (m)'
            call plot_header(ishisdisy,var_num,form2_string)
         end if
         if (ishisdisz .ne. 0) then
! Output z displacement
            info_string = info_string(ic1:ic2) // 'Z displacement '
            ic2 = len_trim(info_string) + 1
            title_string = 'Z Displacement (m)'
            call plot_header(ishisdisz,var_num,form2_string)
         end if
         if (ishisstr .ne. 0) then
! Output strain
            info_string = info_string(ic1:ic2) // 'Volume strain  '
            ic2 = len_trim(info_string) + 1
            title_string = 'Volume Strain'
            call plot_header(ishisstr,var_num,form2_string)
         end if
         if (ishisstrx .ne. 0) then
! Output x stress
            info_string = info_string(ic1:ic2) // 'X stress '
            ic2 = len_trim(info_string) + 1
            title_string = 'X Stress (MPa)'
            call plot_header(ishisstrx,var_num,form2_string)
         end if
         if (ishisstry .ne. 0) then
! Output y stress
            info_string = info_string(ic1:ic2) // 'Y stress '
            ic2 = len_trim(info_string) + 1
            title_string = 'Y Stress (MPa)'
            call plot_header(ishisstry,var_num,form2_string)
         end if
         if (ishisstrxy .ne. 0) then
! Output xy stress
            info_string = info_string(ic1:ic2) // 'XY stress '
            ic2 = len_trim(info_string) + 1
            title_string = 'XY Stress (MPa)'
            call plot_header(ishisstrxy,var_num,form2_string)
         end if
         if (ishisstrz .ne. 0) then
! Output z stress
               info_string = info_string(ic1:ic2) // 'Z stress '
               title_string = 'Z Stress (MPa)'
               call plot_header(ishisstrz,var_num,form2_string)
            ic2 = len_trim(info_string) + 1
         end if
         if (ishisstrxz .ne. 0) then
! Output xz stress
            info_string = info_string(ic1:ic2) // 'XY stress '
            ic2 = len_trim(info_string) + 1
            title_string = 'XZ Stress (MPa)'
            call plot_header(ishisstrxz,var_num,form2_string)
         end if
         if (ishisstryz .ne. 0) then
! Output yz stress
            info_string = info_string(ic1:ic2) // 'Y stress '
            ic2 = len_trim(info_string) + 1
            title_string = 'YZ Stress (MPa)'
            call plot_header(ishisstryz,var_num,form2_string)
         end if
         if (ishiswt .ne. 0) then
! output water table elevations
            formz_string = ''
            info_string = info_string(ic1:ic2) // 'Water table '
            ic2 = len_trim(info_string) + 1
            title_string = "Water table"
            if (form_flag .eq. 1) then
               if (ifree .eq. 0) then
                  formz_string = 'variables = "' // trim(time_string) //
     &                 '"' // ' "X (m)" "Y (m)" "Z (m)" "WT elev (m)" '
     &                 // '"Porosity" "Node" "Out node"'
               else
                  formz_string = 'variables = "' // trim(time_string) //
     &                 '"' // ' "X (m)" "Y (m)" "Z (m)" "WT elev (m)" '
     &                 // '"WT elev2 (m)" "Porosity" "Node" "Out node"'
               end if               
               write(ishiswt, '(a)') trim(formz_string)
               write(ishiswt, 230) 50., 95., trim(wdd)
               write(ishiswt, 230) 50., 90., trim(title_string)
            else if (form_flag .eq. 2) then
               if (ifree .eq. 0) then
                  formz_string = trim(time_string) //
     &                 ", X (m), Y (m), Z (m), WT elev (m), "
     &                 // "Porosity, Node, Out node"
               else
                  formz_string = trim(time_string) //
     &                 ", X (m), Y (m), Z (m), WT elev (m), "
     &                 // "WT elev2 (m), Porosity, Node, Out node"
               end if
               write(ishiswt, '(a)') trim(formz_string)
            else
               write(ishiswt, '(a)') trim(title_string)
               if (ifree .eq. 0) then  
                  write(ishiswt, 1000)
               else
                  write(ishiswt, 1005)
               end if
            end if
         end if

 1000    format('Time (days)', 6x, 'X (m)', 10x, 'Y (m)', 10x, 'Z (m)',
     &        11x, 'WT elev (m)', 5x, 'Porosity', 
     &        10x, 'Node', 10x, 'Out node')
 1005    format('Time (days)', 6x, 'X (m)', 10x, 'Y (m)', 10x, 'Z (m)',
     &        11x, 'WT elev (m)', 5x, 'WT elev2 (m)', 5x, 'Porosity', 
     &        10x, 'Node', 10x, 'Out node')
         if (ishisg .ne. 0 ) then
! Output global variables
            glob_string = ''
            info_string = info_string(ic1:ic2) // 'global '
            if (form_flag .eq. 1) then
               dls = '" "'
               glob_string = 'VARIABLES = "'
            else if (form_flag .eq. 2) then
               dls = ", "
            else
               dls = " "
            end if
            select case (glob_flag)
            case (1)
! Global Mass & Energy Balances
               glob_string = trim(glob_string) //
     &              trim(time_string) // dls //
     &              'Total mass in system (kg)' // dls //
     &              'Total mass of steam in system (kg)' // dls //
     &              'Water discharge (kg)' // dls //
     &              'Water input (kg)' // dls //
     &              'Total water discharge (kg)' // dls //
     &              'Total water input (kg)' // dls //
     &              'Net water discharge (kg)' // dls //
     &              'Delta mass (kg)' // dls //
     &              'Total enthalpy in system (MJ)' // dls //
     &              'Enthalpy discharge (MJ)' // dls //
     &              'Enthalpy input (MJ)' //
     &              'Total enthalpy discharge (MJ)' // dls //
     &              'Total enthalpy input (MJ)' // dls //
     &              'Net enthalpy discharge (MJ)' // dls //
     &              'Delta enthalpy (MJ)'

            case (2)
! Global Water & Air Balances
               glob_string = trim(glob_string) //
     &              trim(time_string) // dls // 
     &              'Total water in system (kg)' // dls //
     &              'Total mass of steam in system (kg)' // dls //
     &              'Water discharge (kg)' // dls //
     &              'Water input (kg)' // dls //
     &              'Total water discharge (kg)' // dls //
     &              'Total water input (kg)' // dls // 
     &              'Net water discharge (kg)' // dls //
     &              'Delta mass (kg)' // dls //
     &              'Total air in system (kg)' // dls //
     &              'Air discharge (kg)' // dls // 
     &              'Air input (kg)' // dls // 
     &              'Total air discharge (kg)' // dls // 
     &              'Total air input kg (kg/s)' // dls //
     &              'Net air discharge (kg)' // dls //
     &              'Delta air mass (kg)'
           case (3)
! Output mass or water balances (excluding steam)
               glob_string = trim(glob_string) //
     &              trim(time_string) // dls //
     &              'Total water in system (kg)' // dls //
     &              'Water discharge (kg) ' // dls //
     &              'Water input (kg)' // dls //
     &              'Total water discharge (kg)' // dls //
     &              'Total water input (kg)' // dls //
     &              'Net (kg) water discharge' // dls //
     &              'Delta mass (kg)'
            case (4)
! Output mass or water balances (including steam)
               glob_string = trim(glob_string) //
     &              trim(time_string) // dls // 
     &              'Total water in system (kg)' // dls //
     &              'Total mass of steam in system (kg) ' // dls //
     &              'Water discharge (kg)' // dls //
     &              'Water input (kg)' // dls //
     &              'Total water discharge (kg)' // dls //
     &              'Total water input (kg)' // dls //
     &              'Net water discharge (kg)' // dls //
     &              'Delta mass (kg)'
            case (5)
! Output air or vapor balances
               glob_string = trim(glob_string) //
     &              trim(time_string) // dls // 
     &              'Total air in system (kg)' // dls //
     &              'Air discharge (kg)' // dls //
     &              'Air input (kg)' //  dls //
     &              'Total air discharge (kg)' // dls // 
     &              'Total air input kg (kg/s)' // dls //
     &              'Net air discharge (kg)' // dls //
     &              'Delta air mass (kg)'
            case (6)
! Output energy balances
               glob_string = trim(glob_string) //
     &              trim(time_string) //  dls //
     &              'Total enthalpy in system (MJ)' // dls //
     &              'Enthalpy discharge (MJ)' // dls //
     &              'Enthalpy input (MJ)' // dls //
     &              'Total enthalpy discharge MJ' // dls //
     &              'Total enthalpy input (MJ)' // dls //
     &              'Net enthalpy discharge (MJ)' // dls //
     &              'Delta enthalpy (MJ)'
            end select
            if (form_flag .eq. 1) 
     &           glob_string = trim(glob_string) // '"'
            length = len_trim (glob_string)
            write(ishisg, '(a)') glob_string(1:length)

         end if
         write (ishis, 6000) trim(file_format)
         write (ishis, 6005)
         length = len_trim (info_string)
         write (ishis,'(a)') info_string(1:length)
         if (out_zones) then
            write (ishis, 6015)
         else
            write (ishis, 6010)
         end if
 6000    format ('Output file format: ', a)
 6005    format ('Parameters written to individual history files:')
 6010    format ('for the following nodes:')
 6015    format ('for the following nodes and zones:')
 6020    format (i8, 3(1x,g16.9))
 6025    format ('Number of output ', a,  ': ', i4)
 6030    format (i7, ' Nodes averaged in Zone ', i4)

c**** write number of plot nodes, node numbers and coordinates ****

         write (ishis, 6025)  'nodes', m
         do i = 1, m
            mi = nskw(i)
            if (mi .gt. neq*2) then
               mi = mi - neq*2
            else if (mi .gt. neq) then 
               mi = mi - neq
            end if
            write(ishis, 6020) nskw(i), corz(mi,1), corz(mi,2),
     *           corz(mi,3)
         end do

! List after nodes and coordinates number of nodes in averaged zones
         if (out_zones) then
            write (ishis, 6025) 'zones', node_azones
            node_ptr_num => node_head_num
            node_ptr => node_head
            start = 0
            end = 0
            do j = 1, node_azones
               i = node_ptr_num%node_number
               start = end + 1
               end = end + i
               write(ishis, 6030) i, dumlist(j)
               ic1 = 1

               ic2 = 8
               count = 0
               do k = start, end
                  count = count + 1
                  if (count .lt. 10) then
                     write (info_string(ic1:ic2), '(i7.7, x)') 
     &                    node_ptr%node_number
                     ic1 = ic2 + 1
                     ic2 = ic2 + 8
                  else
                     count = 0
                     write (info_string(ic1:ic2), '(i7.7)')
     &                    node_ptr%node_number
                     ic1 = 1
                     ic2 = 8
                     write (ishis, '(a80)') info_string
                  end if
                  node_ptr =>node_ptr%nnp
               end do
               if (count .ne. 0)  then
                  write (ishis, '(80a)') info_string(1:ic2-8)
               end if
               node_ptr_num =>node_ptr_num%nnp
               
            end do
         end if
         write(ishis, *)
         close (ishis)
         deallocate (var_string, var_tmp)
         if (allocated(dumlist)) deallocate (dumlist)
! define form_string for numerical output
         if (form_flag .le. 1) then
! Standard or Tecplot
            write (form1_string, 300) max (m, 1)
            write (form2_string, 300) var_num
         else if (form_flag .eq. 2) then
! Surfer
            write (form1_string, 310) max (m, 1)
            write (form2_string, 310) var_num
         end if
 300     format ("(g16.9, ", i5, "(1x, g16.9))")
 310     format ("(g16.9, ", i5, '(", ", g16.9))')
      else if (igf .eq. 2) then  
         if (allocated(dumlist)) deallocate (dumlist)
         if (allocated(dumv)) deallocate (dumv)
      end if

      if (time_flag .eq. 1) then
         ptime = abs(days / 365.25d00)
      else if (time_flag .eq. 2) then
         ptime = abs(days)
      else if (time_flag .eq. 3) then
         ptime = abs(days * 86400.d00)
      else if (time_flag .eq. 4) then
         ptime = abs(days * 24.d00)
      end if

      if (l .eq. 0) then
         time2print = .TRUE.
      else if (ifinsh .ne. 0) then
         if (ptime .ne. last_time) then
            time2print = .TRUE.
         else
            time2print = .FALSE.
         end if
      else if (l .eq. 1 .and. last_step .ne. 0) then
! Need to print if starting a transient solution following
! a steady state (last_step & last_time will be reset at end
         time2print = .TRUE.
      else 
         if (ifinsh .ne. 2 .and. l .lt. (last_step + nhist) .and. 
     &        ptime .lt. (last_time + histime)) then
            time2print = .FALSE.
         else
            time2print = .TRUE.
         end if
      end if

      if (.not. time2print) return

! If zone average values are requested
      if (out_zones) call hstz

      if (.not. allocated(dumv)) then
         if (m .gt. 0) then
            allocate (dumv(m))
         else
            allocate (dumv(1))
         end if
      end if
      if (ishisp .ne. 0 ) then
! Output pressures in MPa
         if (pres_flag .eq. 3 .or. pres_flag .eq. 5 .or. 
     &        pres_flag .eq. 6 .or. pres_flag .eq. 7) then
            do i = 1, m
               if(abs(pcp(nskw(i))).lt.1.d-98) then
                  dumv(i) = 0.0
               else
                  dumv(i) = pcp(nskw(i))
               endif
            end do
         end if
         if (pres_flag .eq. 1 .or. pres_flag .eq. 2 .or. 
     &        pres_flag .eq. 3) then
! Output total/water pressure
            if (out_zones) then
               write(ishisp, form2_string) ptime,  
     &              (max(phi(nskw(i))-phi_inc,1.d-98), i=1,m),
     &              (avg_values(j,pflag), j=1,node_azones)
            else
               write(ishisp, form2_string) ptime,  
     &              (max(phi(nskw(i))-phi_inc,1.d-98), i=1,m)
            end if
            if (pres_flag .eq. 2 .or. pres_flag .eq. 3) then
c gaz 080923 added pa output for isothermal 2 phase
             if(ico2.gt.0) then
              write(ishisp, form1_string) ptime,  
     &           (max(pci(nskw(i)),1.d-98), i=1,m)
             else
              write(ishisp, form1_string) ptime,  
     &           (max(phi(nskw(i))-pcp(nskw(i)),1.d-98), i=1,m)
             endif
            if (pres_flag .eq. 3)
     &             write(ishisp3, form1_string) ptime, (dumv(i), i=1,m)
            end if
         else if (pres_flag .eq. 4) then
c gaz 080923 added pa output for isothermal 2 phase
            if(ico2.gt.0) then
             write(ishisp, form1_string) ptime,  
     &           (max(pci(nskw(i)),1.d-98), i=1,m)
            else
             write(ishisp, form1_string) ptime,  
     &           (max(phi(nskw(i))-pcp(nskw(i)),1.d-98), i=1,m)
            endif
         else if (pres_flag .eq. 5) then
            write(ishisp, form1_string) ptime, (dumv(i), i=1,m)
         else if (pres_flag .eq. 6 .or. pres_flag .eq. 7) then
            if (pres_flag .eq. 6) then
               write(ishisp, form1_string) ptime, (max(phi(nskw(i)),
     &              1.d-98), i=1,m)
            else
c gaz 080923 added pa output for isothermal 2 phase
            if(ico2.gt.0) then
             write(ishisp, form1_string) ptime,  
     &           (max(pci(nskw(i)),1.d-98), i=1,m)
            else
             write(ishisp, form1_string) ptime,  
     &           (max(phi(nskw(i))-pcp(nskw(i)),1.d-98), i=1,m)
            endif
            end if
            write(ishisp2, form1_string) ptime, (dumv(i), i=1,m)
         else if (pres_flag .eq. 8 .or. pres_flag .eq. 9) then
            write(ishisp, form1_string) ptime, (max(phico2(nskw(i)),
     &              1.d-98), i = 1, m)
            if (pres_flag .eq. 9) write(ishisp2, form1_string) ptime, 
     &           (max(phi(nskw(i)), 1.d-98), i = 1, m)
         end if
         call flush(ishisp)
         if (ishisp2 .ne. 0) call flush(ishisp2)
         if (ishisp3 .ne. 0) call flush(ishisp3)
      end if
      if (ishist .ne. 0 ) then
! Output temperature in degrees C
         if (out_zones) then
            write(ishist, form2_string) ptime, 
     &           (max(t(nskw(i)), 1.d-98), i= 1, m),
     &           (avg_values(j,tflag), j=1,node_azones)
         else
            write(ishist, form2_string) ptime, 
     &           (max(t(nskw(i)), 1.d-98), i= 1, m)
         end if
         call flush(ishist)
      end if
      if(ishiswt.ne.0) then
! Output water table elevations
         if (l .gt. 0) then
            if (ifree .eq. 0) then
               if (.not. allocated(col)) then
                  if (.not. allocated(izone_free_nodes)) then
                      allocate(izone_free_nodes(n0))
                      izone_free_nodes=0
                   end if
                   call wtsi_column
                end if
                call airctr(12,0)
            else
               call wtsictr(11)
            end if
            call flush(ishiswt)
         end if
      endif
      if (ishishd .ne. 0 ) then
         if (ihead .ne. 0) then
            do i = 1, m
               call headctr(4,nskw(i),phi(nskw(i)),headdum)
               dumv(i)=max(headdum,0.0d00)
c gaz 7-22-05
               if (ifree .ne. 0) then 
                  if(rlxyf(nskw(i)).lt.rlptol+sattol) dumv(i) = 0.0d0
               endif
            end do
         else if (ichead .ne. 0) then
            ihead=1
            dumconv = crl(1,1)
            dumconv1 = crl(4,1)
            pdum = pres0+rol0*head0*(-grav)
            tdum = temp0        
            call water_density(tdum,pdum,rolconv)
            crl(1,1)=rolconv
            crl(4,1)=pres0
            rho1grav = rolconv*9.81d-6
            do i = 1, m
               call headctr(4,nskw(i),phi(nskw(i)),headdum)
               dumv(i)=max(headdum,0.0d00)
            end do
            crl(1,1)= dumconv
            crl(4,1)= dumconv1
            ihead=0
            if(ico2.lt.0) then
               rho1grav = crl(1,1)*(9.81d-6)
            else
               rho1grav = rol0*9.81d-6
            endif
         end if
         if (ishishd .eq. ishis + 120) then
! Output in columns (default)
            if (out_zones) then
               write(ishishd, form2_string) ptime, (dumv(i), i= 1, m) ,
     &              (avg_values(j,hflag), j=1,node_azones)
            else
               write(ishishd, form2_string) ptime, (dumv(i), i= 1, m)
            end if    
         else if (ishishd .eq. ishis + 121) then
c first, allocate
            if (.not. allocated(dumlist)) 
     &           allocate (dumlist(node_azones))
            zone_ptr => zone_head
            do i = 1, node_azones
               dumlist(i) = zone_ptr%node_number
               zone_ptr =>zone_ptr%nnp
            end do
c
! Output with node/zone number in lines
            do i=1,m
               write(ishishd, formh_string) ptime, nskw(i), dumv(i)
            end do
            if (out_zones) then
               do j=1,node_azones
                  write(ishishd, formh_string) ptime, dumlist(j),
     &                 avg_values(j,hflag)
               end do
            end if
         else if (ishishd .eq. ishis + 125) then
! Output heads in ft
            if (out_zones) then
               write(ishishd, form2_string) ptime, (dumv(i)/0.3048, 
     &              i= 1, m),(avg_values(j,hflag)/0.3048, 
     &              j=1,node_azones)
            else
               write(ishishd, form2_string) ptime, (dumv(i)/0.3048, 
     &              i= 1, m)
            end if    
         end if
         call flush(ishishd)
      end if
      if (ishiss .ne. 0 ) then
c gaz 081223 s zero threshold set to 1.d-98
! Output saturations
         do i = 1, m
            if (ps(nskw(i)) .le. 0.) then
c zero porosity node
               dumv(i) = 0.d0
            else if (icarb .ne. 0) then
c water saturation for CO2 problem
               dumv(i) = max(fw(nskw(i)),1.d-98)
            else if (irdof .ne. 13 .or. ifree .ne. 0) then
c               dumv(i) = max(s(nskw(i)), rlptol)
c               if (dumv(i) .le. rlptol) dumv(i) = 0.d0
c saturations are never zeroed out, report what is in array
               dumv(i) = max(s(nskw(i)),1.d-98)
            else
c water only problem
               dumv(i) = 1.0d0
            end if
         end do
         write(ishiss, form1_string) ptime, (dumv(i), i= 1, m)
         call flush(ishiss)
      end if
      if (ishiswc .ne. 0 ) then
! Output water contents
         do i = 1, m
            if (ifree .ne. 0) then
               dumv(i) = max(s(nskw(i)), rlptol)
	         dumv(i) = dumv(i)*ps(nskw(i))
               if (dumv(i) .le. rlptol) dumv(i) = 0.d0
            else
               dumv(i) = max(s(nskw(i)), 1.d-98)
	         dumv(i) = dumv(i)*ps(nskw(i))
            end if
         end do
         write(ishiswc, form1_string) ptime, (dumv(i), i= 1, m)
         call flush(ishiswc)
      end if
      if (ishisfa .ne. 0 ) then
! Output gas flow in kg/s
         do i = 1, m
           if(ico2.lt.0) then
c isothermal case
            if (abs(qh(nskw(i))).lt.1.d-98) then
               dumv(i) = 0.0
            else
               dumv(i)=qh(nskw(i))
            endif
           else if(ico2.gt.0) then
c non isothermal case
            if (abs(qc(nskw(i))).lt.1.d-98) then
               dumv(i) = 0.0
            else
               dumv(i)=qc(nskw(i))
            endif
           endif
         end do
         write(ishisfa, form1_string) ptime, (dumv(i), i= 1, m)
         call flush(ishisfa)
      end if
      if (ishisf .ne. 0 ) then
! Output flow in kg/s
         do i = 1, m
            if (abs(sk(nskw(i))).lt.1.d-98) then
               dumv(i) = 0.0
            else
               dumv(i)=sk(nskw(i))
            endif
         end do
         write(ishisf, form1_string) ptime, (dumv(i), i= 1, m)
         call flush(ishisf)
      end if
      if (ishise .ne. 0 ) then
! Output enthalpy in MJ/kg
         if (out_zones) then
c average currently computed for qh
c            write(ishise, form2_string) ptime,
c     &           (enlf(nskw(i)), i= 1, m),
c     &           (avg_values(j,eflag), j=1,node_azones)
         else
            write(ishise, form2_string) ptime,
     &           (enlf(nskw(i)), i= 1, m)
         end if
         call flush(ishise)
      end if
      if (ishisef .ne. 0 ) then
! Output enthalpy flow in MJ/s
         if (out_zones) then
            write(ishisef, form2_string) ptime,
     &           (qh(nskw(i)), i= 1, m),
     &           (avg_values(j,eflag), j=1,node_azones)
         else
            write(ishisef, form2_string) ptime,
     &           (qh(nskw(i)), i= 1, m)
         end if
         call flush(ishisef)
      end if
      if (ishisd .ne. 0) then
! Output density
         if (den_flag .eq. 1 .or. den_flag .eq. 2) then
            write(ishisd, form1_string) ptime, (rolf(nskw(i)), i = 1, m)
            if (den_flag .eq. 1)
     &           write(ishisd2, form1_string) ptime, 
     &           (rovf(nskw(i)), i = 1, m)
         end if
         if (den_flag .eq. 3) then
            write(ishisd, form1_string) ptime, (rovf(nskw(i)), i = 1, m)
         end if
         if (den_flag .eq. 4 .or. den_flag .eq. 7 .or. 
     &        den_flag .eq. 8) then
            write(ishisd, form1_string) ptime, 
     &           (wat_prop(nskw(i)), i = 1, m)
            if (den_flag .eq. 4 .or. den_flag .eq. 7) 
     &           write(ishisd2, form1_string) ptime, 
     &           (co2_prop(nskw(i)), i = 1, m)
            if (den_flag .eq. 8) write(ishisd2, form1_string) ptime, 
     &           (co2_prop(9*neq + nskw(i)), i = 1, m)
            if (den_flag .eq. 4) write(ishisd3, form1_string) ptime, 
     &           (co2_prop(9*neq + nskw(i)), i = 1, m)
         end if
         if (den_flag .eq. 5 .or. den_flag .eq. 9) then
            write(ishisd, form1_string) ptime, 
     &           (co2_prop(nskw(i)), i = 1, m)
            if (den_flag .eq. 9) write(ishisd2, form1_string) ptime, 
     &           (co2_prop(9*neq + nskw(i)), i = 1, m)
         end if
         if (den_flag .eq. 6) then
             write(ishisd, form1_string) ptime, 
     &           (co2_prop(9*neq + nskw(i)), i = 1, m)
         end if
         call flush(ishisd)
         if (ishisd2 .ne. 0) call flush(ishisd2) 
         if (ishisd3 .ne. 0) call flush(ishisd3) 
      end if
      if (ishisv .ne. 0) then
! Output viscosity
         if (vis_flag .le. 3) then
            ndummy = 0
            if (rlpnew) then
               call rlp_cap(ndummy)
            else
               call rlperm(ndummy,1)
            end if
         end if
         if (vis_flag .eq. 1 .or. vis_flag .eq. 2) then
            do i = 1, m
               dumv(i) = rlf(nskw(i)) * rolf(nskw(i)) / dil(nskw(i))
            end do
            write(ishisv, form1_string) ptime, (dumv(i), i = 1, m)
         end if
         if (vis_flag .eq. 1 .or. vis_flag .eq. 3) then
            do i = 1, m
               dumv(i) = rvf(nskw(i)) * rovf(nskw(i)) / div(nskw(i))
            end do
            if (vis_flag .eq. 1) then
               write(ishisv2, form1_string) ptime, (dumv(i), i = 1, m)
            else
               write(ishisv, form1_string) ptime, (dumv(i), i = 1, m)
            end if
         end if
         if (vis_flag .eq. 4 .or. vis_flag .eq. 7 .or. 
     &        vis_flag .eq. 8) then
            write(ishisv, form1_string) ptime, 
     &           (wat_prop(8*neq + nskw(i)), i = 1, m)
            if (vis_flag .eq. 4 .or. vis_flag .eq. 7) 
     &           write(ishisv2, form1_string) ptime, 
     &           (co2_prop(6*neq + nskw(i)), i = 1, m)
            if (vis_flag .eq. 8) write(ishisv2, form1_string) ptime, 
     &           (co2_prop(15*neq + nskw(i)), i = 1, m)
            if (vis_flag .eq. 4) write(ishisv3, form1_string) ptime, 
     &           (co2_prop(15*neq + nskw(i)), i = 1, m)
         end if
         if (vis_flag .eq. 5 .or. vis_flag .eq. 9) then
            write(ishisv, form1_string) ptime, 
     &           (co2_prop(6*neq + nskw(i)), i = 1, m)
            if (vis_flag .eq. 9) write(ishisv2, form1_string) ptime, 
     &           (co2_prop(15*neq + nskw(i)), i = 1, m)
         end if
         if (vis_flag .eq. 6) then
             write(ishisv, form1_string) ptime, 
     &           (co2_prop(15*neq + nskw(i)), i = 1, m)
         end if
         call flush(ishisv)
         if (ishisv2 .ne. 0) call flush(ishisv2) 
         if (ishisv3 .ne. 0) call flush(ishisv3) 
      end if
      if (ishishm .ne. 0 ) then
! Output humidity
c gaz 021421  humidity available only for ngas)          
       if(ico2.le.0.or..not.allocated(humida)) then
          write(ishishm,*) 
     &     'humidity output not available-need "ngas" or eqiv)'
       endif
         do i = 1, m          
c            pwatersat =  psat(t(nskw(i)),dpdummy,0)            
c            dumv(i) = phi(nskw(i))-pci(nskw(i))/pwatersat
         dumv(i) = humida(nskw(i))   
         end do
         write(ishishm, form1_string) ptime, (max(dumv(i), 1.d-98), 
     &        i= 1, m) 
         call flush(ishishm)
      end if
C RJP 04/30/07
      if (ishiscm .ne. 0 ) then
! Output CO2 mass in Kg
         if (out_zones) then
            write(ishiscm, form2_string) ptime,
     &           (denco2h(nskw(i))*volume(nskw(i)),i=1,m),
     &           (avg_values(j,carbflag), j=1,node_azones)
         else
            write(ishiscm, form2_string) ptime,
     &           (denco2h(nskw(i))*volume(nskw(i)),i=1,m)
         end if
         call flush(ishiscm)
      end if
      if (ishiscmf .ne. 0 ) then
! Output free CO2 mass fraction
            write(ishiscmf, form1_string) ptime,
     &           (max(xc(nskw(i)),0.d0), i=1,m)
         call flush(ishiscmf)
      end if
      if (ishiscmd .ne. 0 ) then
! Output dissolved CO2 mass fraction
            write(ishiscmd, form1_string) ptime,
     &           (max(yc(nskw(i)),0.d0), i=1,m)
         call flush(ishiscmd)
      end if
C RJP 08/09/07
      if (ishiscsl .ne. 0 ) then
! Output CO2 liquid sats only for nodes if present
         write(ishiscsl, form1_string) ptime, 
     &        (max(fl(nskw(i)),0.d0), i=1,m)
      end if
      if (ishiscsg .ne. 0 ) then
! Output CO2 gaseous sats only for nodes if present
         write(ishiscsg, form1_string) ptime,
     &        (max(fg(nskw(i)),0.d0),i=1,m)
      end if
      if (ishisfz .ne. 0) then
! Output zone fluxes
         call flxz (2, ptime)
         call eflxz(2, ptime)
      end if
! RJP 7/5/07 added following for co2 zone fluxes
      if (ishiscfz .ne. 0) then
! Output zone fluxes
         call flxz (3, ptime)
      end if
      if (ishisdisx .ne. 0 ) then
! Output x displacement
       if(idisp_rel .eq. 0) then
          write(ishisdisx, form1_string) ptime, (du(nskw(i)), i= 1, m) 
       else
          write(ishisdisx, form1_string) ptime, 
     &         (du(nskw(i))-du_ini(nskw(i)), i= 1, m)        
       endif  
         call flush(ishisdisx)
      end if
      if (ishisdisy .ne. 0 ) then
! Output y displacement
         if(idisp_rel .eq. 0) then
            write(ishisdisy, form1_string) ptime, (dv(nskw(i)), i= 1, m)
         else
            write(ishisdisy, form1_string) ptime, 
     &           (dv(nskw(i))-dv_ini(nskw(i)), i= 1, m)        
         endif  
         call flush(ishisdisy)
      end if
      if (ishisdisz .ne. 0 ) then
! Output z displacement
         if(idisp_rel .eq. 0) then
            write(ishisdisz, form1_string) ptime, (dw(nskw(i)), i= 1, m)
         else
            write(ishisdisz, form1_string) ptime, 
     &           (dw(nskw(i))-dw_ini(nskw(i)), i= 1, m)        
         endif  
         call flush(ishisdisz)
      end if
      if (ishisstr .ne. 0 ) then
! Output strain
         write(ishisstr, form1_string) ptime, (vol_strain(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstr)
      end if
      if (ishisstrx .ne. 0 ) then
! Output x stress
         write(ishisstrx, form1_string) ptime, (str_x(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstrx)
      end if
      if (ishisstry .ne. 0 ) then
! Output y stress
         write(ishisstry, form1_string) ptime, (str_y(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstry)
      end if
      if (ishisstrxy .ne. 0 ) then
! Output xy stress
         write(ishisstrxy, form1_string) ptime, (str_xy(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstrxy)
      end if
      if (ishisstrz .ne. 0 ) then
! Output z stress
         write(ishisstrz, form1_string) ptime, (str_z(nskw(i)), 
     &        i= 1, m)
         call flush(ishisstrz)
      end if
      if (ishisstrxz .ne. 0 ) then
! Output xz stress
         write(ishisstrxz, form1_string) ptime, (str_xz(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstrxz)
      end if
      if (ishisstryz .ne. 0 ) then
! Output yz stress
         write(ishisstryz, form1_string) ptime, (str_yz(nskw(i)), 
     &        i= 1, m) 
         call flush(ishisstryz)
      end if

      if (ishisg .ne. 0 ) then
! Output global variables
         if (l .eq. 0) then
            start_mass = amass
            start_ae = aener
         end if
         if (days .ge. 0) then
            dmass = amass - start_mass
            dae = aener - start_ae
            select case (glob_flag)
            case (1) ! Mass & Energy
               write(ishisg, '(16(g16.9, 1x))') ptime, amass, asteam, 
     &              qtoti, curinflow, toutfl, totalflin, qt, dmass,
     &              aener, qtotei, cureinflow, teoutf, totalein, qte,
     &              dae
            case (2)  ! Water & Air
               write(ishisg, '(16(g16.9, 1x))') ptime, amass, asteam,  
     &              qtoti, curinflow, toutfl, totalflin, qt, dmass,
     &              aener, qtotei, cureinflow, teoutf, totalein, qte,
     &              dae
            case (3)  ! Water (no steam)
               write(ishisg, '(8(g16.9, 1x))') ptime, amass, qtoti, 
     &              curinflow, toutfl, totalflin, qt, dmass
            case (4)  ! Water & Steam
               write(ishisg, '(9(g16.9, 1x))') ptime, amass, asteam, 
     &              qtoti, curinflow, toutfl, totalflin, qt, dmass
            case (5)  ! Air
               write(ishisg, '(8(g16.9, 1x))') ptime, aener, qtotei, 
     &              cureinflow, teoutf, totalein, qte, dae
            case (6)  ! Energy
               write(ishisg, '(8(g16.9, 1x))') ptime, aener, qtotei,
     &               cureinflow, teoutf, totalein, qte, dae
            end select
         end if
         call flush(ishisg)
      end if

      last_step = l
      last_time = ptime

      contains
      subroutine plot_header(iunit, varnum, form_string) 

      implicit none
      integer iunit, varnum
      character(*) form_string

      select case (form_flag)
      case (1)
! Tecplot
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(i)),i = 1, varnum)
         write(iunit,230) 50., 95., trim(wdd)
         write(iunit,230) 50., 90., trim(title_string)     
      case (2)
! Surfer               
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(i)), i = 1, varnum)
      case default
         write(iunit, '(a)') trim(title_string)
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(i)), i = 1, varnum)
      end select

 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')

      end subroutine plot_header

      end subroutine plot_new
