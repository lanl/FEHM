      subroutine  diskwrite
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
      use comco2
      use comdi
      use comdti
      use comfi
      use comflow
      use comii
      use commeth
      use compart
      use comriv
      use comxi
      use davidi
      implicit none

      real*8 dummyreal,tolw,sat_dum
      character*11 file_format
      parameter(tolw=1.d-99,sat_dum= 1.d00)
      integer jx,mi,ncount,nc,nc1,nc2,nc3,nd1,nd2,nd3,j,nx,ny,ii,i
      logical log_flag

      if (header_flag .eq. 'new' .and. iriver .eq. 0) then
c Force old style output if wellbore macro is invoked for now
         call diskwrite_new
         return
      else
c Look at rstw array to determine if initial pressures and temperatures
c should be output. All original variables will be output for 'old'
         if (iriver .ne. 0) then
            header_flag = 'old'
            write (ierr, 10)
            if (iout .ne. 0) write(iout, 10)
            if (iptty .ne. 0) write(iptty, 10)
         end if
         if (.not. allocated(rstw)) rstw_num = 0
         do i = 1, rstw_num
            if (rstw(i) .eq. 'pini') then
               ipini = 1
               exit
            end if
         end do
      end if
 10   format (/, 'WARNING: New format restart output not supported for ',
     &     'models using well / river macro', /)

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

cc
c if necessary, correct for head increment before saving
c only converted here for 'write' density not denfined on
c read portion 
c correction now defined only once in input 16-Aug-06
c
c      if(ihead.ne.0.and.head0.gt.0) then
c         phi_inc = head0*crl(1,1)*(-grav)
c      endif

! For compatibility with previous versions 
! Liquid fluxes should be written to restart file
      if (abs(prnt_rst) .eq. 2) flux_flag = 'liquid flux'
           
! Determine size of flux arrays
      if (numflux .ne. 0. .and. .not. compute_flow) then
! Use value of numflux that was read in
      else if (idpdp .ne. 0) then
         numflux = 2*ldna + neq
      else
         numflux = ldna
      end if

      if (cform(7) .eq. 'formatted') then
c Formatted output
         write(isave, 6000)  verno, jdate, jtime, wdd
 6000    format(a30, 3x, a11, 3x, a8, /, a80)
         write(isave ,*)  days
c     write descriptors
         if (iriver .ne. 0) then
            if (icarb .eq. 1) then
                write(isave,'(a4)') 'wcar'
            else
               if (ico2.gt.0) write(isave,'(a4)') 'wnga'
               if (ico2.lt.0) write(isave,'(a4)') 'wair'
               if (ico2.eq.0) write(isave,'(a4)') 'wh20'
            end if
         else
            if (icarb .eq. 1) then
               write(isave,'(a4)') 'carb'
            else
               if (ico2.gt.0) write(isave,'(a4)') 'ngas'
               if (ico2.lt.0) write(isave,'(a4)') 'air '
               if (ico2.eq.0) write(isave,'(a4)') 'h20 '
            end if
         end if
         if ((iccen.eq.0).and..not.ptrak) write(isave,'(a4)') 'ntra'
         if (iccen.ne.0) write(isave,'(a4)') 'trac'
         if (ptrak) write(isave,'(a4)') 'ptrk'
         if (istrs.eq.0) write(isave,'(a4)') 'nstr'
         if (istrs.ne.0) write(isave,'(a4)') 'strs'
         if (idpdp.eq.0) write(isave,'(a4)') 'ndpd'
         if (idpdp.ne.0) write(isave,'(a4)') 'dpdp'
         if (idualp.eq.0) write(isave,'(a4)') 'ndua'
         if (idualp.ne.0) write(isave,'(a4)') 'dual'
         if (ico2.gt.0) then                      
            write(isave ,6002)  (max(to(mi),tolw),   mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave ,6002)  (max(s(mi),tolw),    mi=1,n )
            else
               write(isave ,6002)  (sat_dum,    mi=1,n )
            endif
            write(isave ,6002)  (max(pho(mi),tolw),  mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave ,6002)  (max(pcio(mi),tolw), mi=1,n )
            else
               write(isave ,6002)  (max(pho(mi),tolw),  mi=1,n )
            endif
         else if(ico2.eq.0 .or. icarb .eq. 1) then
            write(isave ,6002)  (max(to(mi),tolw),   mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave ,6002)  (max(s(mi),tolw),    mi=1,n )
            else
               write(isave ,6002)  (sat_dum,    mi=1,n )
            endif
            write(isave ,6002)  (max(pho(mi),tolw),  mi=1,n )
            if (icarb .eq. 1) then
               write(isave ,6002)  (max(toco2(mi),tolw),   mi=1,n )
               write(isave ,6002)  (max(phoco2(mi),tolw),  mi=1,n )
               write(isave ,6002)  (max(fow(mi),tolw),    mi=1,n )
               write(isave ,6002)  (max(fog(mi),tolw),    mi=1,n )
               write(isave ,6002)  (max(fol(mi),tolw),    mi=1,n )
               write(isave ,6003)  (ieoso(mi),    mi=1,n )
               write(isave ,6003)  (iceso(mi),    mi=1,n )
            end if
         else if(ico2.lt.0) then
            if(ihead.ne.0 .or. (irdof .eq. 13 .and. 
     &           abs(ifree) .ne. 1)) then
               write(isave ,6002)  (sat_dum,    mi=1,n )
            else
               write(isave ,6002)  (max(s(mi),tolw),    mi=1,n )
            endif
            write(isave ,6002) (pho(mi)-phi_inc,  mi=1,n )
c     write(isave ,6002) (max(pho(mi)-phi_inc,tolw),  mi=1,n )
         endif
         
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
c     write descriptors
         if (iriver .ne. 0) then
            if (icarb .eq. 1) then
               write(isave) 'wcar'
            else
               if (ico2.gt.0) write(isave) 'wnga'
               if (ico2.lt.0) write(isave) 'wair'
               if (ico2.eq.0) write(isave) 'wh20'
            end if
         else
            if (icarb .eq. 1) then
               write(isave) 'carb'
            else
               if (ico2.gt.0) write(isave) 'ngas'
               if (ico2.lt.0) write(isave) 'air '
               if (ico2.eq.0) write(isave) 'h20 '
            end if
         end if
         if ((iccen.eq.0).and..not.ptrak) write(isave) 'ntra'
         if (iccen.ne.0) write(isave) 'trac'
         if (ptrak) write(isave) 'ptrk'
         if (istrs.eq.0) write(isave) 'nstr'
         if (istrs.ne.0) write(isave) 'strs'
         if (idpdp.eq.0) write(isave) 'ndpd'
         if (idpdp.ne.0) write(isave) 'dpdp'
         if (idualp.eq.0) write(isave) 'ndua'
         if (idualp.ne.0) write(isave) 'dual'
         if (ico2.gt.0) then                      
            write(isave)  (max(to(mi),tolw),   mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave)  (max(s(mi),tolw),    mi=1,n )
            else
               write(isave)  (sat_dum,    mi=1,n )
            endif
            write(isave)  (max(pho(mi),tolw),  mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave)  (max(pcio(mi),tolw), mi=1,n )
            else
               write(isave)  (max(pho(mi),tolw),  mi=1,n )
            endif
         else if(ico2.eq.0 .or. icarb .eq. 1) then
            write(isave)  (max(to(mi),tolw),   mi=1,n )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               write(isave)  (max(s(mi),tolw),    mi=1,n )
            else
               write(isave)  (sat_dum,    mi=1,n )
            endif
            write(isave)  (max(pho(mi),tolw),  mi=1,n )
            if (icarb .eq. 1) then
               write(isave)  (max(toco2(mi),tolw),   mi=1,n )
               write(isave)  (max(phoco2(mi),tolw),  mi=1,n )
               write(isave)  (max(fow(mi),tolw),    mi=1,n )
               write(isave)  (max(fog(mi),tolw),    mi=1,n )
               write(isave)  (max(fol(mi),tolw),    mi=1,n )
               write(isave)  (ieoso(mi),    mi=1,n )
               write(isave)  (iceso(mi),    mi=1,n )
            end if
         else if(ico2.lt.0) then
            if(ihead.ne.0 .or. (irdof .eq. 13 .and. 
     &           abs(ifree) .ne. 1)) then
               write(isave)  (sat_dum,    mi=1,n )
            else
               write(isave)  (max(s(mi),tolw),    mi=1,n )
            endif
            write(isave) (pho(mi)-phi_inc,  mi=1,n )
c     write(isave) (max(pho(mi)-phi_inc,tolw),  mi=1,n )
         endif
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
                      write(isave, 6002) (a_vxy(i), i = 1,numflux)
                     else
                      write(isave, 6002) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'liquid flux'
                     write (isave) flux_flag
                     write (isave) numflux
                     write (isave) (a_axy(i), i = 1,numflux)
                     if (iout .ne. 0) then
                        write(iout,*) 'WARNING: Liquid only problem, ',
     &                       'vapor flux'
                        write(iout,*) 'can not be written to restart ',
     &                       'file'
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
                      write(isave) (a_vxy(i), i = 1,numflux)
                     else
                      write(isave) (0.0d0, i = 1,numflux)
                     endif
                  else
                     flux_flag = 'no fluxes  '
                     write (isave, '(a11)') flux_flag
                     if (iout .ne. 0) then
                        write(iout, *) 'WARNING: Liquid only problem,',
     &                       ' vapor flux'
                        write(iout, *) 'can not be written to restart',
     &                       ' file'
                     end if
                  end if
               end if
            end if
         end if
      endif

      nsave  =  1

      log_flag = .TRUE.
      if (iccen.eq.1) call diskc(log_flag)
      if (ptrak) call diskp(log_flag)

c     if (istrs.ne.0) call  stressctr  (16,0 )
c write initial pressures and temperatures if enabled via restart macro
      if(ipini .ne. 0) then
         if (cform(7) .eq. 'formatted') then
            write(isave,*)  ( tini (mi) , mi=1,n )
            write(isave,*)  ( phini (mi) , mi=1,n )
         else
            write(isave)  ( tini (mi) , mi=1,n )
            write(isave)  ( phini (mi) , mi=1,n )
         endif
      endif
      close  ( isave )
      
      return
      end
      
