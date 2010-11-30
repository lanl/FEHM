      subroutine inrlp
!***********************************************************************
! Copyright 2009 Los Alamos National Security, LLC  All rights reserved
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
!D1
!D1 PURPOSE
!D1
!D1 Read parameters for relative permeability and capillary pressure 
!D1 models.
!D1
!***********************************************************************

      use comai
      use comco2, only : icarb
      use comcomp, only : ioil
      use comdi, only : irlp
      use comdti, only : n0
      use comki
      use comrlp

      implicit none

      integer :: lastmodel = 0
      integer :: lastrun = 0
      integer :: lasttbl = 0
      integer i, it, j, j2, jt, k, k2, maxphase, maxcpl, ndx, nparams
      integer table_unit, open_file, cn, tblnum, ip, ip2
      integer nwds, imsg(20), msg(20)
      integer, allocatable :: irlptmp(:)
      real*8 xmsg(20)
      real*8 scutm, smcut
      parameter(scutm = 1.d-03, maxphase = 10, maxcpl = 10)
      character*32 cmsg(20), dum_string
      character*80 chdum
      character*200 table_file
      logical null1, null_new

      save lastmodel, lastrun, lasttbl

      maxrp = 3
      maxcp = 2
      msg = 0
      imsg = 0
      xmsg = 0.
      cmsg = ''
      macro = 'rlpm'

      if (lastrun .ne. irun) then
         lastrun = irun
         lastmodel = 0
      end if

      allocate (irlptmp(n0))
      if (.not. allocated (rlp_group)) then
         allocate (rlp_group(nrlp), rlp_phase(nrlp, maxrp))
         allocate (rlp_type(nrlp, maxrp), rlp_pos(nrlp, maxphase))
         allocate (rlp_param(nrlp, maxrp * max_rp))
         allocate (rlp_fparam(nrlp, maxrp * max_rpf))
         rlp_pos = 0
         rlp_phase = 0
         rlp_type = 0
         rlp_param = 0.
         rlp_fparam = 0.
         allocate (cap_coupling(nrlp, maxcp), cap_type(nrlp, maxcp))
         allocate (cap_param(nrlp, max_cp*maxcp)) 
         allocate (cap_fparam(nrlp, max_cp*maxcp))
         allocate (cap_pos(nrlp, maxcpl))
         cap_pos = 0
         cap_coupling = 0
         cap_type = 0
         cap_param = 0.
         if (ntable .ne. 0) then
            allocate (tblindx(ntable, 2), rlp_table(ntblines, 4))
            tblindx = 0
            rlp_table = 0.
         end if
      end if
      
      i = lastmodel
      j = 0
      it = lasttbl
      do
         read (inpt, '(a80)') chdum
         if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or. 
     &        chdum(1:3) .eq. 'END') exit
         call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
         select case (cmsg(1))
         case ('group', 'GROUP')
c i is the current working group number
            i = i + 1
            if (j .eq. 1) then
               write (ierr, 50) rlp_group(lastmodel)
               stop
            else
c j is the current phase number
               j = 1
            end if
            k2 = 1
            lastmodel = i
            rlp_group(i) = imsg(2)
            goto 100
         case ('table', 'TABLE')
            if (nwds .lt. 5) write(ierr, 60) 5, cmsg(2), rlp_group(i)
            rlp_type(i, j) = 8
            cap_type(i, k2) = 8
            it = it + 1
            if (it .eq. 1) then
               tblindx(it,1) = 1
               ndx = 0
            else
               tblindx(it, 1) = tblindx(it - 1, 2) + 1
               ndx = tblindx(it - 1, 2)
            end if
            tblnum = imsg(2)
            nparams = imsg(3)
            if (nparams .eq. 3) then
               ip = 1
            else if (nparams .eq. 4) then
               ip = 2
            end if
            do ip2 = 1, ip
               if (ip2 .eq. 2) then
                  j = j + 1
                  rlp_type(i, j) = 8
               end if
               select case (cmsg(3 + ip2))
               case ('h2o_liquid','H2O_LIQUID','water','WATER')
                  if (icarb .eq. 0 .and. ioil.eq.0) then
                     rlp_phase(i, j) = 1
                  else
                     rlp_phase(i, j) = -1
                  end if
               case ('air','AIR')
                  rlp_phase(i, j) = 2
               case ('h2o_gas','H2O_GAS','vapor','VAPOR')
                  rlp_phase(i, j) = 5
               case ('co2_gas','CO2_GAS')
                  rlp_phase(i, j) = 4
               case ('co2_liquid','CO2_LIQUID')
                  rlp_phase(i, j) = 3
               case ('co2_sc','CO2_SC')
                  rlp_phase(i, j) = 3
               case ('methane_hydrate','METHANE_HYDRATE')
                  rlp_phase(i, j) = 6
               case ('oil', 'OIL')
                  rlp_phase(i, j) = 7
               case ('gas', 'GAS')
                  rlp_phase(i, j) = 8
               case ('oil_water', 'OIL_WATER')
                  rlp_phase(i, j) = 9
               case ('oil_gas', 'OIL_GAS')
                  rlp_phase(i, j) = 10
               end select
               rlp_pos(i, abs(rlp_phase(i, j))) = j
               k = (j - 1) * max_rp
               rlp_param(i, k + 1) = it
               rlp_param(i, k + 2) = ip2 + 1
               if (ip2 .eq. 1) then
                  rlp_param(i, k + 3) = rlp_phase(i, j)
               else
                  rlp_param(i, k + 3) = rlp_phase(i, j - 1)
               end if
            end do
            select case (cmsg(3 + ip2))
            case ('air/water', 'water/air')
               cap_coupling(i, k2) = 1
            case ('water/co2_liquid', 'co2_liquid/water')
               cap_coupling(i, k2) = 2
            case ('water/co2_gas', 'co2_gas/water')
               cap_coupling(i, k2) = 3
            case ('co2_liquid/co2_gas', 'co2_gas/co2_liquid')
               cap_coupling(i, k2) = 4
            case ('water/vapor', 'vapor/water')
               cap_coupling(i, k2) = 5
            case ('air/vapor', 'vapor/air')
               cap_coupling(i, k2) = 6
            case ('water/oil', 'oil/water')
               cap_coupling(i, k2) = 7
            case ('gas/water', 'water/gas')
               cap_coupling(i, k2) = 8
            case ('gas/oil', 'oil/gas')
               cap_coupling(i, k2) = 9
            case ('liquid/gas', 'gas/liquid')
               cap_coupling(i, k2) = 10
            end select
            cap_pos(i, cap_coupling(i, k2)) = k2
            k = (k2 - 1) * max_cp
            cap_param(i, k + 1) = it
            k2 = k2 + 1
            read (inpt, '(a80)') chdum
            if (chdum(1:4) .eq. 'file') then
c Data is read from file
               read (inpt, '(a200)') table_file
               table_unit = open_file(table_file,'old')
c Read past any header lines in the table (header lines should start with a character)
               do 
                  read (table_unit, '(a80)') chdum
                  call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
                  if (msg(1) .ne. 3) then
                     backspace (table_unit)
                     exit
                  end if
               end do
            else
c Data is found on the following lines
               table_unit = inpt
               backspace (inpt)
            end if
            do
               read (table_unit, '(a80)', end = 5) chdum
c Input is terminated with a blank line or 'end' or end-of-file)
               if (null_new(chdum) .or. chdum(1:3) .eq. 'end') exit
               ndx = ndx + 1
               read (chdum, *) (rlp_table(ndx,cn), cn = 1, nparams)
c Capillary pressure is put into position 4
               if (nparams .eq. 3) then
                  rlp_table(ndx,4) = rlp_table(ndx,3)
                  rlp_table(ndx,3) = 0.
               end if
            end do
 5          lasttbl = it
            tblindx(it, 2) = ndx
            if (table_unit .ne. inpt) close (table_unit)
            cmsg(2) = 'tabular'
c End of case 'tabular'
         case ('air', 'AIR')
            rlp_phase(i, j) = 2
         case ('h2o_gas','H2O_GAS','vapor','VAPOR')
            rlp_phase(i, j) = 5
         case ('h2o_liquid','H2O_LIQUID','water','WATER')
            if (icarb .eq. 0 .and. ioil .eq. 0) then
               rlp_phase(i, j) = 1
            else
               rlp_phase(i, j) = -1
            end if
         case ('co2_gas','CO2_GAS')
            rlp_phase(i, j) = 4
         case ('co2_liquid','CO2_LIQUID')
            rlp_phase(i, j) = 3
         case ('co2_sc','CO2_SC')
            rlp_phase(i, j) = 3
         case ('methane_hydrate','METHANE_HYDRATE')
            rlp_phase(i, j) = 6
         case ('oil', 'OIL')
            rlp_phase(i, j) = 7
         case ('gas', 'GAS')
            rlp_phase(i, j) = 8
         case ('oil_water', 'OIL_WATER')
            rlp_phase(i, j) = 9
         case ('oil_gas', 'OIL_GAS')
            rlp_phase(i, j) = 10
         case ('cap','CAP')
            select case (cmsg(2))
            case ('air/water', 'water/air')
               cap_coupling(i, k2) = 1
            case ('water/co2_liquid', 'co2_liquid/water')
               cap_coupling(i, k2) = 2
            case ('water/co2_gas', 'co2_gas/water')
               cap_coupling(i, k2) = 3
            case ('co2_liquid/co2_gas', 'co2_gas/co2_liquid')
               cap_coupling(i, k2) = 4
            case ('water/vapor', 'vapor/water')
               cap_coupling(i, k2) = 5
            case ('methane_hydrate/water', 'water/methane_hydrate')
               cap_coupling(i, k2) = 6
            case ('oil/water', 'water/oil')
               cap_coupling(i, k2) = 7
            case ('gas/water', 'water/gas')
               cap_coupling(i, k2) = 8
            case ('oil/gas', 'gas/oil')
               cap_coupling(i, k2) = 9
            case ('liquid/gas', 'gas/liquid')
               cap_coupling(i, k2) = 10
            case default
               write (ierr, 40) trim(cmsg(2)), rlp_group(i)
               stop
            end select
         case default
            write (ierr, 10) trim(cmsg(1)), rlp_group(i)
            stop
         end select

         if (cmsg(1) .eq. 'cap' .or. cmsg(1) .eq. 'CAP') then
            cap_pos(i, cap_coupling(i, k2)) = k2
c Read parameters associated with capillary model
            k = (k2 - 1) * max_cp
            select case (cmsg(3))
            case ('linear')
               if (nwds .lt. 4) write(ierr, 60) 4, cmsg(3), rlp_group(i)
               cap_type(i, k2) = 1
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
            case ('linear_for')
               if (nwds .lt. 4) write(ierr, 60) 4, cmsg(3), rlp_group(i)
               cap_type(i, k2) = 2
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
               if (cap_param(i, k + 2) .gt. 0.) then
                  cap_param(i, k + 1) = cap_param(i, k + 1) / 
     &                 cap_param(i, k + 2)
               else
                  cap_param(i, k + 1) = 0.
               end if
            case ('exponential')
               if (nwds .lt. 5) write(ierr, 60) 5, cmsg(3), rlp_group(i)
               cap_type(i, k2) = 3
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
               cap_param(i, k + 3) = xmsg(6) + imsg(6)
            case ('brooks-corey')
               if (nwds .lt. 6) write(ierr, 60) 6, cmsg(3), rlp_group(i)
               cap_type(i, j) = 4
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
               cap_param(i, k + 3) = xmsg(6) + imsg(6)
               cap_param(i, k + 4) = xmsg(7) + imsg(7)
               cap_param(i, k + 5) = xmsg(8) + imsg(8)
               cap_param(i, k + 6) = xmsg(9) + imsg(9)
            case ('vg', 'vg_cp', 'vg_cap', 'vg_capeq', 'vg_cpeq',
     &              'vg_cpneq', 'vg_capneq', 'vg_ek')
               if (nwds .lt. 6) write(ierr, 60) 6, cmsg(3), rlp_group(i)
               if (cmsg(3) .eq. 'vg') then
                  cap_type(i, k2) = 5
               else if (cmsg(3) .eq. 'vg_cpneq' .or. cmsg(3) .eq. 
     &                 'vg_capneq') then
                  cap_type(i, k2) = 7
               else if (cmsg(3) .eq. 'vg_ek') then
                  cap_type(i, k2) = 9
               else
                  cap_type(i, k2) = 6
               end if
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
               cap_param(i, k + 3) = xmsg(6) + imsg(6)
               cap_param(i, k + 4) = xmsg(7) + imsg(7)
               cap_param(i, k + 5) = xmsg(8) + imsg(8)
               cap_param(i, k + 6) = xmsg(9) + imsg(9)
               smcut = (cap_param(i, k + 6) - cap_param(i, k + 1)) /
     &              (cap_param(i, k + 2) - cap_param(i, k + 1))
               cap_param(i, k + 6) = max (smcut, scutm)
               call vgcap_fit2 (0, i, k2)
               dum_string = trim(cmsg(3)) // ' fracture'
               read (inpt, '(a80)') chdum
               if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or. 
     &           chdum(1:3) .eq. 'END') exit
               call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
               if (cmsg(1) .eq. 'fracture') then
                  if (nwds .lt. 7) write(ierr, 60) 6, dum_string, 
     &                 rlp_group(i)
                  cap_fparam(i, k + 1) = xmsg(2) + imsg(2)
                  cap_fparam(i, k + 2) = xmsg(3) + imsg(3)
                  cap_fparam(i, k + 3) = xmsg(4) + imsg(4)
                  cap_fparam(i, k + 4) = xmsg(5) + imsg(5)
                  cap_fparam(i, k + 5) = xmsg(6) + imsg(6)
                  cap_fparam(i, k + 6) = xmsg(7) + imsg(7)
                  smcut = (cap_fparam(i, k + 6) - cap_fparam(i, k + 1))
     &                  / (cap_fparam(i, k + 2) - cap_fparam(i, k + 1))
                  cap_fparam(i, k + 6) = max (smcut, scutm)
                  call vgcap_fit2 (1, i, k2)
               else
                  backspace inpt
               end if
            case default
               write (ierr, 20) trim(cmsg(3)), rlp_group(i)
               stop
            end select
            k2 = k2 + 1
         else if (cmsg(1) .ne. 'table' .and. cmsg(1) .ne. 'TABLE') then
            rlp_pos(i, abs(rlp_phase(i, j))) = j
c     Read parameters associated with rlp phase
            k = (j - 1) * max_rp
            select case (cmsg(2))
            case ('constant')
               if (nwds .lt. 3) write(ierr, 60) 1, cmsg(2), rlp_group(i)
               rlp_type(i, j) = 1
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
            case ('linear')
               if (nwds .lt. 4) write(ierr, 60) 2, cmsg(2), rlp_group(i)
               rlp_type(i, j) = 2
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
               rlp_param(i, k + 2) = xmsg(4) + imsg(4)
            case ('exponential')
               if (nwds .lt. 5) write(ierr, 60) 3, cmsg(2), rlp_group(i)
               rlp_type(i, j) = 3
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
               rlp_param(i, k + 2) = xmsg(4) + imsg(4)
               rlp_param(i, k + 3) = xmsg(5) + imsg(5)
            case ('corey')
c     Restrict to two-phase models
               if (nwds .lt. 4) write(ierr, 60) 2, cmsg(2), rlp_group(i)
               if (j .eq. 3) then
                  write (ierr, 30) trim(cmsg(2)), rlp_group(i)
                  stop
               end if
               rlp_type(i, j) = 4
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
               rlp_param(i, k + 2) = xmsg(4) + imsg(4)
            case ('brooks-corey')
               if (nwds .lt. 5) write(ierr, 60) 3, cmsg(2), rlp_group(i)
               rlp_type(i, j) = 5
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
               rlp_param(i, k + 2) = xmsg(4) + imsg(4)
               rlp_param(i, k + 3) = xmsg(5) + imsg(5)
            case ('vg', 'vg_cp', 'vg_cap')
c     Restrict to two-phase models
               if (nwds .lt. 6) write(ierr, 60) 4, cmsg(2), rlp_group(i)
               if (j .eq. 3) then
                  write (ierr, 30) trim(cmsg(3)), lastmodel
                  stop
               end if
               if (cmsg(2) .eq. 'vg') then
                  rlp_type(i, j) = 6
               else
                  rlp_type(i, j) = 7
               end if
               rlp_param(i, k + 1) = xmsg(3) + imsg(3)
               rlp_param(i, k + 2) = xmsg(4) + imsg(4)
               rlp_param(i, k + 3) = xmsg(5) + imsg(5)
               rlp_param(i, k + 4) = xmsg(6) + imsg(6)
               dum_string = trim(cmsg(2)) // ' fracture'
               read (inpt, '(a80)') chdum
               if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or. 
     &           chdum(1:3) .eq. 'END') exit
               call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
               if (cmsg(1) .eq. 'fracture') then
                  if (nwds .lt. 8) write(ierr, 65) 7, dum_string, 
     &                 rlp_group(i)
                  rlp_fparam(i, k + 1) = xmsg(2) + imsg(2)
                  rlp_fparam(i, k + 2) = xmsg(3) + imsg(3)
                  rlp_fparam(i, k + 3) = xmsg(4) + imsg(4)
                  rlp_fparam(i, k + 4) = xmsg(5) + imsg(5)
                  rlp_fparam(i, k + 5) = xmsg(6) + imsg(6)
                  rlp_fparam(i, k + 6) = xmsg(7) + imsg(7)
                  rlp_fparam(i, k + 7) = xmsg(8) + imsg(8)
                  if (nwds .lt. 9) then
                     rlp_fparam(i, k + 8) = 0.
                  else
                     rlp_fparam(i, k + 8) = xmsg(9) + imsg(9)
                  end if
                  if (rlp_fparam(i, k + 8) .ne. 0.) 
     &                 call rlp_frac2(0, 0, 0.0d00, 0.0d00, 0.0d00,
     &                 0.0d00, 0.0d00, 0.0d00, rlp_fparam(i, k + 8))
               else
                  backspace inpt
               end if
            case ('same')
               rlp_type(i, j) = rlp_type(i, j - 1)
               do j2 = 1, max_rp
                  rlp_param(i, k + j2) = rlp_param(i, j2)
               end do
               k = (j - 1) * max_rpf
               do j2 = 1, max_rpf
                  rlp_fparam(i, k + j2) = rlp_fparam(i, j2)
               end do
            case ('stone')
c For 3-phase water/oil/gas, oil rel perm in terms of krow and krog
c Table corresponding to oil_water, oil_gas phase and
c connate water saturation
               rlp_type(i, j) = 9
               rlp_param(i, k+1) = it - 1
               rlp_param(i, k+2) = it
               rlp_param(i, k+3) = xmsg(3) + imsg(3)
            case ('tabular')
c If type is table store table number in 1st rlp parameter position
               rlp_param(i, k+1) = it
            case default
               write (ierr, 20) trim(cmsg(3)), rlp_group(i)
               stop
            end select
            j = j + 1
         end if
 100  end do

      narrays = 1
      itype(1) = 4
      default(1) = 0
      macro = "rlpm"
      igroup = 2
      call initdata2( inpt, ischk, n0, narrays,
     2     itype, default, macroread(7), macro, igroup, ireturn,
     3     i4_1=irlptmp(1:n0) )

      do i = 1, n0
         do j = 1, nrlp
            if (rlp_group(j) .eq. irlptmp(i)) then
               irlp(i) = j
               exit
            end if
         end do
      end do
      
      macroread(7) = .TRUE.

      deallocate (irlptmp)

 10   format ('Unrecognized phase ', a, ' in rlpm, check model ', i4) 
 20   format ('Unrecognized type ', a, ' in rlpm, check model ', i4)
 30   format ('Invalid model type ', a, ' in rlpm, check model ', i4)
 40   format ('Invalid coupling ', a, ' in rlpm, check model ', i4)
 50   format ('At least two phases must be specified in rlpm, ',
     &     'check model ', i4)
 60   format (i1, 'parameters are required for model type', a,
     & ' in rlpm, check model ', i4)
 65   format ('At least ', i1, 'parameters are required for model type',
     &     a, ' in rlpm, check model ', i4)

      end subroutine inrlp
