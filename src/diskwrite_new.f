      subroutine  diskwrite_new
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
!D1 To write information from disk file for restart purposes.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/diskwrite.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3 2.7 Provide Restart Capability
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
!D6
!D6 GLOBAL OBJECTS
!D6
!D6 Global Constants
!D6
!D6   None
!D6
!D6 Global Types
!D6
!D6   None
!D6
!D6 Global Variables
!D6 
!D6 isave, jdate, jtime, nmfil, verno, wdd, days, ico2, iccen, istrs,
!D6 idpdp, idualp, to, s, pho, pcio, iread, wdd1, neq, t, phi, pci,
!D6 ieos, so
!D6
!D6 Global Subprograms
!D6
!D6 Identifier    Type     Description
!D6
!D6 dated         N/A      Determines current date and time
!D6 concen        N/A      Called to read or write concentration data
!D6 stress        N/A      Called to read or write stress data
!D6
!***********************************************************************
!D7
!D7 LOCAL IDENTIFIERS
!D7
!D7 Local Constants
!D7
!D7   None
!D7
!D7 Local Types
!D7
!D7 Local variables
!D7
!D7   Identifier      Type     Description
!D7
!D7   mi              INT      Index for writing array values
!D7   ncount          INT      Index for writing array values
!D7   nc              INT      Number of node sets in current run
!D7   nc1             INT      Index used for setting initial values
!D7   nc2             INT      Index used for setting initial values
!D7   nc3             INT      Index used for setting initial values
!D7   nd1             INT      Index used for setting initial values
!D7   nd2             INT      Index used for setting initial values
!D7   nd3             INT      Index used for setting initial values
!D7   j               INT      Do loop index
!D7   nx              INT      Index used for setting initial values
!D7   ny              INT      Index used for setting initial values
!D7   ny              INT      Index used for setting initial values
!D7   ii              INT      Do loop index
!D7   i               INT      Index used for setting initial values
!D7
!D7 Local Subprograms
!D7
!D7   None
!D7
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS BEGIN diskwrite
!PS
!PS   IF the run specifies a restart file to be written
!PS     Rewind tape
!PS     Write FEHM version number and current date and time
!PS     Write time of simulation
!PS     Write flags to identify options used in this run
!PS
!PS     IF this is a noncondensible gas run with heat
!PS       Write temperature, saturation, total pressure, and ...
!PS       ... gas pressure values
!PS     ELSEIF this is a water/water vapor run with heat
!PS       Write temperature, saturation, total pressure, pressure values
!PS     ELSEIF this is an isothermal air/water run
!PS       Write saturation and pressure values
!PS     ENDIF
!PS
!PS     concen - write concentration values to the file
!PS     stress - write stress information to the file
!PS
!PS   ENDIF
!PS
!PS END diskwrite
!PS
!***********************************************************************

      use comai
      use combi
      use comci, only : rolf
      use comdi
      use comdti
      use comfi
      use comflow
      use comii
      use commeth
      use compart
      use comco2
      use comriv
      use comsi
      use comxi
      use davidi
c gaz 031124      
      use com_prop_data, only : ihenryiso
      implicit none
c
      real*8 dummyreal,tolw,sat_dum
      character*4 dum_type
      character*11 dummy_string
      parameter(tolw=1.d-99,sat_dum= 1.d00)
      integer jx,mi,ncount,nc,nc1,nc2,nc3,nd1,nd2,nd3,j,nx,ny,ii,i
      logical :: write_trac = .FALSE.
      logical :: write_ptrk = .FALSE.
      logical :: write_temp = .FALSE.
      logical :: write_pres = .FALSE.
      logical :: write_gasp = .FALSE.
      logical :: write_sat  = .FALSE.
      logical :: write_mass = .FALSE.
      logical :: write_flux = .FALSE.
      logical :: write_co2  = .FALSE.
      logical :: write_pini = .FALSE.
      logical :: write_por = .FALSE.
c gaz 031124 added restart with isothermal 2 phase w solubility      
      logical :: write_mfrac_iso = .FALSE.
      logical :: write_disp(3) = .FALSE.
      logical :: write_strs(6) = .FALSE.

      call dated (jdate,jtime)
      inquire (isave, opened = ex)
      if (ex) close (isave)
      open (isave,file=nmfil(7),status=cstats(7),form=cform(7))
           
      rewind  isave
      if(ice.ne.0) then
         write(isave, 6000)  verno, jdate, jtime, wdd
         write(isave ,*)  days
         call icectr (-20,0)
         return
      endif

c For compatibility with previous versions 
c Liquid fluxes should be written to restart file
      if (abs(prnt_rst) .eq. 2) flux_flag = 'liquid flux'
           
c Determine size of flux arrays
      if (numflux .ne. 0. .and. .not. compute_flow) then
c Use value of numflux that was read in
      else if (idpdp .ne. 0) then
         numflux = 2*ldna + neq
      else
         numflux = ldna
      end if

      if (.not. allocated(rstw)) then
         allocate(rstw(1)) 
         rstw(1) = 'all'
         rstw_num = 1
      end if

      do i = 1, rstw_num
         select case (rstw(i))
         case ('all')
            if (icarb .eq. 1) then
               write_co2 = .TRUE.
            else
               if (ico2 .ge. 0) write_temp = .TRUE.
               write_pres = .TRUE.
               if (irdof .ne. 13 .and. ihead .eq. 0) write_sat  = .TRUE.
               if (ico2 .gt. 0) write_gasp = .TRUE.
               if (iporos .ne. 0) write_por = .TRUE.
c gaz 031124 added mass frac for isothermal   
               if (ihenryiso .ne. 0) write_mfrac_iso = .TRUE.
            end if
            if (iccen .ne. 0) write_trac = .TRUE.
            if (ptrak) write_ptrk = .TRUE.
            if (istrs .ne. 0) then
               if (icnl .eq. 0) then
                  write_strs = .TRUE.
                  write_disp = .TRUE.
               else
                  write_strs(1) = .TRUE.
                  write_strs(2) = .TRUE.
                  write_strs(4) = .TRUE.
                  write_disp(1) = .TRUE.
                  write_disp(2) = .TRUE.
               end if
            end if
         case ('temp')
            write_temp = .TRUE.
         case ('pres')
            write_pres = .TRUE.
         case ('satu')
            write_sat  = .TRUE.
         case ('trac')
            write_trac = .TRUE.
         case ('ptrk')
            write_ptrk = .TRUE.
         case ('gasp')
            write_gasp = .TRUE.
         case ('co2')
            if (icarb .eq. 1) then
               write_co2 = .TRUE.
            else
               if (iout .ne. 0) then
                  write(iout, 100) 
     &                 'Non-co2 problem, co2 parameters'
               end if  
            end if
         case ('mass') 
            write_mass = .TRUE.
         case ('pini')
            write_pini = .TRUE.
         case ('poro')
            write_por = .TRUE.
         case ('disp')
            write_disp = .TRUE.
            if (icnl .ne. 0) write_disp(3) = .FALSE.
         case ('disx')
            write_disp(1) = .TRUE.
         case ('disy')
            write_disp(2) = .TRUE.
         case ('disz')
            if (icnl .eq. 0) write_disp(3) = .TRUE.
         case ('strs', 'stre')
            if (icnl .eq. 0) then
               write_strs = .TRUE.
            else
               write_strs(1) = .TRUE.
               write_strs(2) = .TRUE.
               write_strs(4) = .TRUE.
            end if
         case ('strx')
            write_strs(1) = .TRUE.
         case ('stry')
            write_strs(2) = .TRUE.
         case ('strz')
            if (icnl .eq. 0 ) write_strs(3) = .TRUE.
         case ('stxy')
            write_strs(4) = .TRUE.
         case ('stxz')
            if (icnl .eq. 0 ) write_strs(5) = .TRUE.
         case ('styz')
            if (icnl .eq. 0 ) write_strs(6) = .TRUE.
         case ('none')
c Values from the restart file will not be used
         end select
      end do
      
      if (istrs .eq. 0) then
         do i = 1, 3
            if (write_disp(i)) then
               write_disp = .FALSE.
               if (iout .ne. 0) then
                  write(iout, 100) 
     &                 'Non-stress problem, displacements'
               end if 
               exit
            end if
         end do
         do i = 1, 6
            if (write_strs(i)) then
               write_strs = .FALSE.
               if (iout .ne. 0) then
                  write(iout, 100) 
     &                 'Non-stress problem, stresses'
               end if 
               exit
            end if
         end do
      end if

      if (idualp .ne. 0) then
         dum_type = 'dual'
      else if (idpdp .ne. 0) then
         dum_type = 'dpdp'
      else
         dum_type = 'nddp'
      end if

      if (cform(7) .eq. 'formatted') then
c Formatted output
         write(isave, 6000)  verno, jdate, jtime, wdd
 6000    format(a30, 3x, a11, 3x, a8, /, a80)
         write(isave, *)  days
         write(isave, '(i9, 1x, a4)') n, dum_type
         if (write_temp) then
            write(isave, '(a11)') 'temperature'
c   gaz 012111 removed tolw for printing neg temperatures (ice models)            
c            write(isave, 6002)  (max(to(mi),tolw),   mi=1,n )
           write(isave, 6002)  (to(mi),   mi=1,n )
         end if
         if (write_sat .and. irdof .ne. 13 .and. ihead .eq. 0) then
            write(isave, '(a11)') 'saturation '
            if(idoff.eq.-1) then
             write(isave, 6002)  (1.0d0,    mi=1,n )
            else
             write(isave, 6002)  (max(s(mi),tolw),    mi=1,n )
            endif
         else
         end if
         if (write_pres) then
            write(isave, '(a11)') 'pressure   '
            if (ico2 .lt. 0) then
               write(isave, 6002) (pho(mi)-phi_inc,  mi=1,n )
            else
               write(isave, 6002)  (max(pho(mi),tolw),  mi=1,n )
            end if
         end if
         if (write_gasp .and. irdof .ne. 13) then
            write(isave, '(a11)') 'gaspressure'
            write(isave, 6002)  (max(pcio(mi),tolw), mi=1,n )
         else
         end if
c gaz 031124 added restart with mass fraction         
         if (write_mfrac_iso .and. irdof .ne. 13) then
           write(isave, '(a10)') 'mfrac_iso '  
c gaz 040424, fo now, just simple isothermal 
c           write(isave, 6002)  (max(frac_gas_iso(mi),tolw), mi=1,n )
           write(isave, 6002)  (max(cnlf(mi),tolw), mi=1,n )
         endif    
         if (write_co2) then
            write(isave, '(a11)') 'co2temperat'
            write(isave, 6002)  (max(toco2(mi),tolw),   mi=1,n )
            write(isave, '(a11)') 'co2pressure'
            write(isave, 6002)  (max(phoco2(mi),0.d0),  mi=1,n )
            write(isave, '(a11)') 'wsaturation'
            write(isave, 6002)  (max(fow(mi),0.d0),    mi=1,n )
            write(isave, '(a11)') 'lco2saturat'
            write(isave, 6002)  (max(fol(mi),0.d0),    mi=1,n )
            write(isave, '(a11)') 'dissolvdco2'
            write(isave, 6002)  (max(yc(mi),0.d0),    mi=1,n )
            write(isave, '(a11)') 'eoswater   '
            write(isave, 6003)  (ieoso(mi),    mi=1,n )
            write(isave, '(a11)') 'eosco2     '
            write(isave, 6003)  (iceso(mi),    mi=1,n )
            write(isave, '(a11)') 'eosdc      '
            write(isave, 6003)  (ico2diso(mi),    mi=1,n )
         end if
         if (write_mass) then
            if (.not. allocated(mass_var)) allocate(mass_var(n))
            write(isave, '(a11)') 'mass       '
            do mi = 1, n
               mass_var(mi) = sx1(mi) * ps(mi) * rolf(mi)
            end do
            write(isave, 6002) ( mass_var(mi), mi=1,n ) 
         end if
         if (write_pini) then
            write(isave, '(a11)') 'initemperat'
            write(isave, 6002)  ( tini (mi) , mi=1,n )
            write(isave, '(a11)') 'inipressure'            
            write(isave, 6002)  ( phini (mi) , mi=1,n )
         end if
         if (write_por) then
            write(isave, '(a11)') 'porosity   '
            write(isave, 6002)  ( ps (mi) , mi=1,n )
         end if
         if (write_disp(1)) then
            write(isave, '(a11)') 'xdisplacmnt'
            write(isave, 6002)  ( du (mi) , mi=1,n )
         end if
         if (write_disp(2)) then
            write(isave, '(a11)') 'ydisplacmnt'
            write(isave, 6002)  ( dv (mi) , mi=1,n )
         end if
         if (write_disp(3)) then
            write(isave, '(a11)') 'zdisplacmnt'
            write(isave, 6002)  ( dw (mi) , mi=1,n )
         end if
         if (write_strs(1)) then
            write(isave, '(a11)') 'xstress    '
            write(isave, 6002)  ( str_x (mi) , mi=1,n )
         end if
         if (write_strs(2)) then
            write(isave, '(a11)') 'ystress    '
            write(isave, 6002)  ( str_y (mi) , mi=1,n )
         end if
         if (write_strs(4)) then
            write(isave, '(a11)') 'xystress   '
            write(isave, 6002)  ( str_xy (mi) , mi=1,n )
         end if
         if (icnl .eq. 0) then
            if (write_strs(3)) then
               write(isave, '(a11)') 'zstress    '
               write(isave, 6002)  ( str_z (mi) , mi=1,n )
            end if
            if (write_strs(5)) then
               write(isave, '(a11)') 'xzstress   '
               write(isave, 6002)  ( str_xz (mi) , mi=1,n )
            end if
            if (write_strs(6)) then
               write(isave, '(a11)') 'yzstress   '
               write(isave, 6002)  ( str_yz (mi) , mi=1,n )
            end if
         end if  
         if (flux_flag(1:2) .eq. 'no') then
            write (isave, '(a11)') flux_flag
         else
            if (idoff .lt. 0) then
               flux_flag = 'no fluxes  '
               write (isave, '(a11)') 'no fluxes  '
               if (iout .ne. 0) then
                  write (iout, *) 'WARNING: Heat conduction only ', 
     &                 'problem'
                  write (iout, *) 'fluxes can not be written to ',
     &                 'restart file'
               end if
            else
               if (flux_flag(1:3) .eq. 'all') then
                  if (irdof .ne. 13) then
                     write (isave, '(a11)') flux_flag
                     write (isave, *) numflux
                     write(isave, 6002) (a_axy(i), i = 1,numflux)
                     if(jswitch.eq.0) then
                      write(isave, 6002) (a_vxy(i), i = 1,numflux)
                     else
                      write(isave, 6002) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'liquid flux'
                     write (isave, '(a11)') flux_flag
                     write (isave, *) numflux
                     write (isave, 6002) (a_axy(i), i = 1,numflux)
                     if (iout .ne. 0) then
                        write(iout, *) 'WARNING: Liquid only problem,',
     &                    ' vapor flux'
                        write(iout, *) 'can not be written to restart',
     &                       ' file'
                     end if
                  end if
               else if (flux_flag(1:3) .eq. 'liq') then
                  write (isave, '(a11)') flux_flag
                  write (isave, *) numflux
                  write(isave, 6002) (a_axy(i), i = 1,numflux)
               else if (flux_flag(1:3) .eq. 'vap') then
                  if (irdof .ne. 13) then
                     write (isave, '(a11)') flux_flag
                     write (isave, *) numflux
                     if(jswitch.eq.0) then
                      write(isave, 6002) (a_vxy(i), i = 1,numflux)
                     else
                      write(isave, 6002) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'no fluxes  '
                     write (isave, '(a11)') flux_flag
                     if (iout .ne. 0) then
                        write(iout, *) 'WARNING: Liquid only problem,',
     &                    ' vapor flux'
                        write(iout, *) 'can not be written to restart',
     &                       ' file'
                     end if
                  end if
               end if
            end if
         end if

 6002    format(4g25.16)
 6003    format(30i3)
        
      else
c Unformatted output
         write(isave)  verno, jdate, jtime, wdd
         write(isave)  days
         write(isave) n, dum_type
         if (write_temp) then
            dummy_string = 'temperature'
            write(isave) dummy_string
c   gaz 012111 removed tolw for printing neg temperatures (ice models)              
c            write(isave)  (max(to(mi),tolw),   mi=1,n )
          write(isave)  (to(mi),   mi=1,n )
         end if
         if (write_sat .and. irdof .ne. 13 .and. ihead .eq. 0) then
            dummy_string = 'saturation '
            write(isave) dummy_string
            write(isave)  (max(s(mi),tolw),    mi=1,n )
         else
         end if
         if (write_pres) then
            dummy_string = 'pressure   '
            write(isave) dummy_string
            if (ico2 .lt. 0) then
               write(isave) (pho(mi)-phi_inc,  mi=1,n )
            else
               write(isave)  (max(pho(mi),tolw),  mi=1,n )
            end if
         end if
         if (write_gasp .and. irdof .ne. 13) then
            dummy_string = 'gaspressure'
            write (isave) dummy_string
            write(isave)  (max(pcio(mi),tolw), mi=1,n )
         else
         end if
         if (write_co2) then
            dummy_string = 'co2temperat'
            write (isave) dummy_string
            write(isave)  (max(toco2(mi),tolw),   mi=1,n )
            dummy_string = 'co2pressure'
            write (isave) dummy_string
            write(isave)  (max(phoco2(mi),tolw),  mi=1,n )
            dummy_string = 'wsaturation'
            write (isave) dummy_string
            write(isave)  (max(fow(mi),tolw),    mi=1,n )
            dummy_string = 'lco2saturat'
            write (isave) dummy_string
            write(isave)  (max(fol(mi),tolw),    mi=1,n )
            dummy_string = 'dissolvdco2'
            write (isave) dummy_string
            write(isave)  (max(yc(mi),tolw),    mi=1,n )
            dummy_string = 'eoswater   ' 
            write (isave) dummy_string
            write(isave)  (ieoso(mi),    mi=1,n )
            dummy_string = 'eosco2      '
            write(isave)  (iceso(mi),    mi=1,n )
         end if
         if (write_mass) then
            dummy_string = 'mass       '           
            write(isave) dummy_string
            if (.not. allocated(mass_var)) allocate(mass_var(n))
            do mi = 1, n
               mass_var(mi) = sx1(mi) * ps(mi) * rolf(mi)
            end do
            write(isave) ( mass_var(mi), mi=1,n ) 
         end if
         if (write_pini) then
            dummy_string = 'initemperat'
            write (isave) dummy_string
            write(isave)  ( tini (mi) , mi=1,n )
            dummy_string = 'inipressure'            
            write (isave) dummy_string
            write(isave)  ( phini (mi) , mi=1,n )
         end if
         if (write_por) then
            dummy_string = 'porosity   '
            write (isave) dummy_string
            write(isave)  ( ps (mi) , mi=1,n )
         end if
          if (write_disp(1)) then
            dummy_string = 'xdisplacmnt'
            write (isave) dummy_string
            write(isave)  ( du (mi) , mi=1,n )
         end if
         if (write_disp(2)) then
            dummy_string = 'ydisplacmnt'
            write (isave) dummy_string
           write(isave)  ( dv (mi) , mi=1,n )
         end if
         if (write_disp(3)) then
            dummy_string = 'zdisplacmnt'
            write (isave) dummy_string
            write(isave)  ( dw (mi) , mi=1,n )
         end if
         if (write_strs(1)) then
            dummy_string = 'xstress    '
            write (isave) dummy_string
            write(isave)  ( str_x (mi) , mi=1,n )
         end if
         if (write_strs(2)) then
            dummy_string = 'ystress    '
            write (isave) dummy_string
            write(isave)  ( str_y (mi) , mi=1,n )
         end if
         if (write_strs(4)) then
            dummy_string = 'xystress   '
            write (isave) dummy_string
            write(isave)  ( str_xy (mi) , mi=1,n )
         end if
         if (icnl .eq. 0) then
            if (write_strs(3)) then
               dummy_string = 'zstress    '
               write (isave) dummy_string
               write(isave)  ( str_z (mi) , mi=1,n )
            end if
            if (write_strs(5)) then
               dummy_string = 'xzstress   '
               write (isave) dummy_string
               write(isave)  ( str_xz (mi) , mi=1,n )
            end if
            if (write_strs(6)) then
               dummy_string = 'yzstress   '
               write (isave) dummy_string
               write(isave)  ( str_yz (mi) , mi=1,n )
            end if
         end if  

         if (flux_flag(1:2) .eq. 'no') then
            write (isave) flux_flag
         else
            if (idoff .lt. 0) then
               flux_flag = 'no fluxes  '
               write (isave) 'no fluxes  '
               if (iout .ne. 0) then
                  write (iout, *) 
     &                 'WARNING: Heat conduction only problem'
                  write (iout, *) 
     &                 'fluxes can not be written to restart file'
               end if
            else
               if (flux_flag(1:3) .eq. 'all') then
                  if (irdof .ne. 13) then
                     write (isave) flux_flag
                     write (isave) numflux
                     write (isave) (a_axy(i), i = 1,numflux)
                     if(jswitch.eq.0) then
                      write (isave) (a_vxy(i), i = 1,numflux)
                     else
                      write (isave) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'liquid flux'
                     write (isave) flux_flag
                     write (isave) numflux
                     write (isave) (a_axy(i), i = 1,numflux)
                     if (iout .ne. 0) then
                        write(iout, 100) 
     &                       'Liquid only problem, vapor flux'
                     end if
                  end if
               else if (flux_flag(1:3) .eq. 'liq') then
                  write (isave) flux_flag
                  write (isave) numflux
                  write (isave) (a_axy(i), i = 1,numflux)
               else if (flux_flag(1:3) .eq. 'vap') then
                  if (irdof .ne. 13) then
                     write (isave) flux_flag
                     write (isave) numflux
                     if(jswitch.eq.0) then
                      write (isave) (a_vxy(i), i = 1,numflux)
                     else
                      write (isave) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'no fluxes  '
                     write (isave, '(a11)') flux_flag
                     if (iout .ne. 0) then
                        write(iout, 100) 
     &                       'Liquid only problem, vapor flux'
                     end if
                  end if
               end if
            end if
         end if
      endif

      nsave  =  1

      if (iccen.eq.1 .and. write_trac) call diskc(write_trac)
      if (ptrak .and. write_ptrk) call diskp(write_ptrk)
      if (istrs.ne.0) call  stress  ( 4 )

      close  ( isave )
      
 100  format ('WARNING: ', a, ' can not be written to restart file')

      return
      end
      
