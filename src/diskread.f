      subroutine  diskread
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
!D1 To read information from disk file for restart purposes.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/diskread.f_a  $
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
!PS BEGIN diskread
!PS
!PS   IF a restart file is specified
!PS     Read the starting time
!PS     Read descriptors for the options used in the run that ...
!PS     ... produced this restart file
!PS     Set number of values to be read for each array
!PS
!PS     IF the previous run was noncondensible gas with heat
!PS       Read values for temperature, saturation, pressure, and ...
!PS       ... air pressure
!PS     ELSEIF the previous run was water/vapor with heat
!PS       Read values for temperature, saturation, and pressure
!PS     ELSEIF the previous run was isothermal air/water
!PS       Read values for saturation and pressure
!PS     ENDIF
!PS
!PS     IF this restart is a dual porosity/dual permeability run
!PS
!PS       IF the previous run was also dpdp
!PS         Set indexes for initializing arrays
!PS       ELSEIF the previous run was dual porosity
!PS         Set indexes for initializing arrays
!PS       ELSE the previous run was single porosity
!PS         Set indexes for initializing arrays
!PS       ENDIF
!PS     
!PS     ELSEIF this restart is a dual porosity run
!PS
!PS       IF the previous run was also dpdp
!PS         Set indexes for initializing arrays
!PS       ELSEIF the previous run was dual porosity
!PS         Set indexes for initializing arrays
!PS        ELSE the previous run was single porosity
!PS         Set indexes for initializing arrays
!PS       ENDIF
!PS     
!PS     ELSE this restart is a single porosity run
!PS
!PS       IF the previous run was also dpdp
!PS         Set indexes for initializing arrays
!PS       ELSEIF the previous run was dual porosity
!PS         Set indexes for initializing arrays
!PS       ELSE the previous run was single porosity
!PS         Set indexes for initializing arrays
!PS       ENDIF
!PS     
!PS     ENDIF
!PS
!PS     FOR each node set in current run
!PS       Set indexes used to initialize variable values
!PS
!PS       FOR each node
!PS         IF this node is to be initialized at the previous value
!PS           Initialize values at this node
!PS
!PS           IF there is no liquid
!PS             Set flag to indicate vapor only, set saturation to 0
!PS           ELSEIF there is no vapor
!PS             Set flag to indicate liquid only, set saturation to 1
!PS           ELSE there is liquid and vapor
!PS             Set flag to indicate two phase
!PS           ENDIF
!PS
!PS         ENDIF
!PS       ENDFOR each node
!PS
!PS     ENDFOR each node set
!PS
!PS     concen - read concentration values for restart
!PS     stress - read stress information for restart
!PS
!PS   ENDIF
!PS
!PS ENDIF
!PS 
!PS END diskread
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
      use comsi
      use comxi
      use davidi
c gaz 031124      
      use com_prop_data, only : ihenryiso
      implicit none

c      real*8 pl, tl, dum1, dumb, dumc(9)
      real*8 pl, tl, dum1, dumb
      real*8 dummyreal,tolw,sat_dum,satr
      parameter(tolw=1.d-99,sat_dum= 1.d00)
      integer jx,mi,ncount,nc,nc1,nc2,nc3,nd1,nd2,nd3,j,nx,ny,ii,i
      integer modneq, modneq_primary, rerr, dummyint, duma
      character*80 dumtitle
      character*30 dumver
      character*11 dumdate, dumflag
      character*8 dumtime
      character*4 :: geom_type = ''
      logical gdpm_read
      logical, dimension (2) :: read_trac = .FALSE.
      logical, dimension (2) :: read_ptrk = .FALSE.
      logical, dimension (2) :: read_temp = .FALSE.
      logical, dimension (2) :: read_pres = .FALSE.
      logical, dimension (2) :: read_gasp = .FALSE.
      logical, dimension (2) :: read_sat  = .FALSE.
      logical, dimension (2) :: read_flux = .FALSE.
      logical, dimension (2) :: read_mass = .FALSE.
      logical, dimension (2) :: read_co2 = .FALSE.
      logical, dimension (2) :: read_pini = .FALSE.
      logical, dimension (2) :: read_por = .FALSE.
      logical, dimension (2) :: read_disx = .FALSE.
      logical, dimension (2) :: read_disy = .FALSE.
      logical, dimension (2) :: read_disz = .FALSE.
      logical, dimension (2) :: read_strx = .FALSE.
      logical, dimension (2) :: read_stry = .FALSE.
      logical, dimension (2) :: read_strz = .FALSE.
      logical, dimension (2) :: read_strxy = .FALSE.
      logical, dimension (2) :: read_strxz = .FALSE.
      logical, dimension (2) :: read_stryz = .FALSE.
c gaz 040924
      if(.not.allocated(read_mfrac_iso)) then
       allocate(read_mfrac_iso(2))
       read_mfrac_iso(1) = .FALSE.
       read_mfrac_iso(2) = .FALSE.
      endif
      residual_stress = .FALSE.

      mass_read = .FALSE.
      co2_read = .FALSE.
      gdpm_read = .FALSE.
      inquire (iread, opened = ex)
      if (ex) close (iread)
      open (iread,file=nmfil(6),status=cstats(6),form=cform(6))

      if (iout .ne. 0) write (iout, 10) trim(nmfil(6))
      if (iptty .ne. 0) write (iptty, 10) trim(nmfil(6))
 10   format ('Reading initial condition data from file ', a)
c the first block of code defines what variables are eligible for reading
c the second block of code does the read
c now we handle the case where there was not a 'read' keyword followed by a list of variables to read
      if (.not. allocated(rstr)) then
         allocate(rstr(1))
         rstr(1) = 'all'
         rstr_num = 1
      end if

      do i = 1, rstr_num
         select case (rstr(i))
         case ('all')
c handle the default case here
c always read pressure and porosity
               read_pres(1) = .TRUE.
               pres_read = .TRUE.
               read_por(1) = .TRUE.
c read other variables depending on the physics
            if (icarb .eq. 1) then 
               read_co2(1) = .TRUE.
            else
               if (ico2 .ge. 0) read_temp(1) = .TRUE.
            endif
c gaz 031924 
            if(ihenryiso.ne.0) then
                read_mfrac_iso(1) = .true.   
            endif    
c              read_pres(1) = .TRUE.
c              pres_read = .TRUE.
               if (irdof .ne. 13) read_sat(1)  = .TRUE.
               if (ico2 .gt. 0) read_gasp(1) = .TRUE.
c              if (iporos .ne. 0) read_por(1) = .TRUE.
c           end if
            if (iccen .ne. 0) read_trac(1) = .TRUE.
            if (ptrak) read_ptrk(1) = .TRUE.
            if (istrs .eq. 1) then
               residual_stress = .TRUE.
               read_disx(1) = .TRUE. 
               read_disy(1) = .TRUE.
               read_disz(1) = .TRUE. 
               read_strx(1) = .TRUE. 
               read_stry(1) = .TRUE. 
               read_strxy(1) = .TRUE.
               if (icnl .eq. 0) then
                  read_strz(1) = .TRUE.
                  read_strxz(1) = .TRUE.
                  read_stryz(1) = .TRUE.
               end if
            end if
c handle the case where there was a 'read' keyword followed by specific variables
         case ('temp')
            read_temp(1) = .TRUE.
         case ('pres')
            read_pres(1) = .TRUE.
            pres_read = .TRUE.
         case ('satu')
            read_sat(1)  = .TRUE.
         case ('trac')
            read_trac(1) = .TRUE.
         case ('ptrk')
            read_ptrk(1) = .TRUE.
         case ('gasp')
            read_gasp(1) = .TRUE.
         case ('co2')
            read_co2(1) = .TRUE.
         case ('mass')
            read_mass(1) = .TRUE.
            mass_read = .TRUE.
         case ('pini')
            read_pini(1) = .TRUE.
            ipini = 1
         case ('poro')
            read_por(1) = .TRUE.
         case ('disp')
            read_disx(1) = .TRUE. 
            read_disy(1) = .TRUE.
            if (icnl .eq. 0 ) read_disz(1) = .TRUE.             
         case ('disx')
            read_disx(1) = .TRUE.
         case ('disy')
            read_disy(1) = .TRUE.
         case ('disz')
            if (icnl .eq. 0 ) read_disz(1) = .TRUE.
         case ('strs', 'stre')
               residual_stress = .TRUE.
               read_strx(1) = .TRUE. 
               read_stry(1) = .TRUE. 
               read_strxy(1) = .TRUE.
               if (icnl .eq. 0) then
                  read_strz(1) = .TRUE.
                  read_strxz(1) = .TRUE.
                  read_stryz(1) = .TRUE.
               end if
         case ('strx')
            read_strx(1) = .TRUE.
         case ('stry')
            read_stry(1) = .TRUE.
         case ('strz')
            if (icnl .eq. 0 ) read_strz(1) = .TRUE.
         case ('stxy')
            read_strxy(1) = .TRUE.
         case ('stxz')
            if (icnl .eq. 0 ) read_strxz(1) = .TRUE.
         case ('styz')
            if (icnl .eq. 0 ) read_stryz(1) = .TRUE.
         case ('none')
c Values from the restart file will not be used
         end select
      end do

      if (read_strx(1)) then
         allocate (str_x0(n0))
         allocate (str_y0(n0))
         allocate (str_xy0(n0))
         if (icnl .eq. 0) then
            allocate (str_z0(n0))
            allocate (str_xz0(n0))
            allocate (str_yz0(n0))
         end if
      end if

      if (.not. compute_flow) read_flux(1) = .TRUE.

      if (.not. read_pres(1)) then
         phi = pho
         pres_read = .FALSE.
      end if

      rewind  iread
      if (cform(6) .eq. 'formatted') then
c Reading formatted file
         read (iread ,   *) dumver, dumdate, dumtime
         read (iread ,   *) dumtitle
         read (iread ,   *)  days

         if(ice.ne.0) then
            call icectr (20,0)
            return
         endif

c     
c     read descriptors
c     
         read(iread,'(a4)') wdd1(1:4)
         if (wdd1(1:4) .eq. 'h20 ' .or. wdd1(1:4) .eq. 'wh20' .or.
     &        wdd1(1:4) .eq. 'air ' .or. wdd1(1:4) .eq. 'wair' .or.
     &        wdd1(1:4) .eq. 'ngas' .or. wdd1(1:4) .eq. 'wnga' .or.
     &        wdd1(1:4) .eq. 'carb' .or. wdd1(1:4) .eq. 'wcar') then 
c this is old style input
            read(iread,'(a4)') wdd1(5:8)
            read(iread,'(a4)') wdd1(9:12)
            read(iread,'(a4)') wdd1(13:16)
            read(iread,'(a4)') wdd1(17:20)
            ncount=neq
            if(wdd1(13:16).eq.'dpdp') ncount=neq+neq
c            if(wdd1(13:16).eq.'dpdp') gdpm_read = .TRUE.
            if(wdd1(17:20).eq.'dual') ncount=neq+neq+neq
            if (wdd1(1:4).eq.'ngas' .or. wdd1(1:4).eq.'wnga') then
               if(iriver.ne.0 .and. wdd1 .eq. 'ngas') 
     &              ncount = ncount - npoint_riv
               if (read_temp(1)) then
                  read(iread ,*)  ( t   (mi) , mi=1,ncount )
                  read_temp(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'temperatures'
                  if (iptty .ne. 0) write (iptty, 400) 'temperatures'
               end if
               if (read_sat(1)) then
                  read(iread ,*)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
               if (read_pres(1)) then
                  read(iread ,*)  ( phi (mi) , mi=1,ncount )
                  read_pres(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (read_gasp(1)) then
                  read(iread ,*)  ( pci (mi) , mi=1,ncount )
                  read_gasp(2) = .TRUE.
c     
c     check for negative partial pressures
c     
                  do i=1,ncount
                     pci(i)=max(0.0d+0,pci(i))
                  enddo
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'gas pressures'
                  if (iptty .ne. 0)write (iptty, 400) 'gas pressures'
               end if
               if (iriver .ne. 0 .and. wdd1(1:4) .eq. 'ngas') then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_temp(1)) 
     &                    t(i) = t(iriver_con_node(i-ncount))
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
                     if (read_gasp(1)) 
     &                    pci(i) = pci(iriver_con_node(i-ncount))
		  enddo
               end if

            else if(wdd1(1:4).eq.'h20 ' .or. wdd1(1:4).eq.'wh20' .or.
     &              wdd1(1:4).eq.'carb' .or. wdd1(1:4).eq.'wcar') then
               if (iriver .ne. 0 .and. (wdd1(1:4) .eq. 'h20 ' .or. 
     &              wdd1(1:4).eq.'carb')) ncount = ncount - npoint_riv
               if (read_temp(1) .or. read_co2(1)) then
                  read(iread ,*)  ( t   (mi) , mi=1,ncount )
                  read_temp(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'temperatures'
                  if (iptty .ne. 0) write (iptty, 400) 'temperatures'
               end if
               if (read_sat(1) .and. irdof .ne. 13) then
                  read(iread ,*)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.                           
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
c gaz 031924                  
              if (read_mfrac_iso(1) .and. irdof .ne. 13) then
                  read(iread ,*)  (frac_gas_iso(mi), mi=1,ncount )
                  read_mfrac_iso(2) = .TRUE.  
               endif   
               if (read_pres(1) .or. read_co2(1)) then
                  read(iread ,*)  ( phi (mi) , mi=1,ncount )
                  pres_read = .TRUE.
                  read_pres(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (wdd1(1:4).eq.'carb' .or. wdd1(1:4) .eq. 'wcar') then
                  read(iread ,*) ( tco2 (mi), mi = 1,ncount )
                  read(iread ,*) ( phico2 (mi), mi = 1,ncount )
                  read(iread ,*) ( fow (mi), mi = 1,ncount )
                  read(iread ,*) ( fog (mi), mi = 1,ncount )
                  read(iread ,*) ( fol (mi), mi = 1,ncount )
                  read(iread ,*) ( ieoso (mi),  mi = 1,ncount )
                  read(iread ,*) ( iceso (mi),  mi = 1,ncount )
                  read(iread ,*) ( ico2diso (mi),  mi = 1,ncount )
                  read_co2(2) = .TRUE.
                  do i = 1, ncount
                     if(inico2flg(i).eq.1) then
c Use values from input
                        fow(i) = fw(i)
                        fol(i) = fl(i)
                        fog(i) = fg(i)
                        yoc(i) = yc(i)
                     endif
                  enddo
                  fw = fow
                  fl = fol
                  fg = fog
                  yc = yoc
                  ieos = ieoso
                  ices = iceso
                  ico2dis = ico2diso
               else
                  if (icarb .eq. 1) then
                     do i = 1, ncount + npoint_riv
                        if(inico2flg(i).ne.1) then
                           fw(i) = 1.d0
                           fg(i) = 0.d0
                           fl(i) = 0.d0
                        endif
                     enddo
                  end if
               end if
               if (iriver .ne. 0 .and. (wdd1(1:4) .eq. 'h20 ' .or. 
     &              wdd1(1:4).eq.'carb')) then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_temp(1) .or. read_co2(1)) 
     &                    t(i) = t(iriver_con_node(i-ncount))
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1) .or. read_co2(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
                     if (wdd1(1:4).eq.'carb') then
                        tco2(i) = tco2(iriver_con_node(i-ncount))
                        phico2(i) = phico2(iriver_con_node(i-ncount))
                        fw(i) = fw(iriver_con_node(i-ncount))
                        fg(i) = fg(iriver_con_node(i-ncount))
                        fl(i) = fl(iriver_con_node(i-ncount))
                        ieoso(i) = ieoso(iriver_con_node(i-ncount))
                        iceso(i) = iceso(iriver_con_node(i-ncount))
                       ico2diso(i) = ico2diso(iriver_con_node(i-ncount))
                     end if
		  enddo
               end if
               if (icarb .eq. 1) then
                  if (wdd1(1:4).eq.'h20 ' .or. wdd1(1:4).eq.'nh20') then
                     tco2 = t
                     phico2 = phi
                     do i = 1, ncount + npoint_riv
                        pl = phico2(i)
                        tl = tco2(i)
                        call co2_properties(1,duma,pl,tl,dum1,
     &                       ices(i),dumb,dumc,duma)
                     enddo
                     iceso = ices
                     ieos = 1
                     ieoso = ieos
				   ico2diso = ico2dis
                  else
                     ieos = ieoso
                     ices = iceso
				   ico2dis = ico2diso
                  end if
                  toco2 = tco2
                  phoco2 = phico2
                  fow = fw
                  fog = fg
                  fol = fl                
               end if

            else if(wdd1(1:4).eq.'air ' .or. wdd1(1:4).eq.'wair') then
               if(iriver.ne.0 .and. wdd1(1:4) .eq. 'air ') 
     &              ncount = ncount - npoint_riv
               if (read_sat(1) .and. irdof .ne. 13) then
                  read(iread ,*)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
               if (read_pres(1)) then
                  read(iread ,*)  ( phi (mi) , mi=1,ncount )
                  read_pres(2) = .TRUE.
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (iriver .ne. 0 .and. wdd1(1:4) .eq. 'air ') then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
		  enddo
               end if
            endif
            
         else
            rerr = 0
            backspace (iread)
            read (iread, *) ncount, geom_type
c gaz 063020 mod is not correct here changing to difference       
c            modneq = mod(ncount,neq)
            modneq = ncount - neq
            if(geom_type.eq.'gdpm') gdpm_read = .TRUE.
            if (gdpm_read) then
               modneq_primary = 0
            else
               modneq_primary = mod(ncount,neq_primary)
            end if
c gaz 063020 diskwrite_new has no gdpm or gdkm            
c            if (modneq .ne. 0 .or. modneq_primary .ne. 0) goto 2000
            if (modneq .ne. 0) goto 2000
            do
               read (iread ,'(a11)', end = 100) dumflag
               select case (dumflag(1:4))
               case ('temp')
                  if (read_temp(1) .or. read_co2(1)) then
                     read (iread, *) ( t(mi), mi = 1,ncount )
                     read_temp(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'temperatures'
                     if (iptty .ne. 0) write (iptty, 400) 'temperatures'
                  end if
               case ('pres')
                  if (read_pres(1) .or. read_co2(1)) then
                     read (iread, *) ( phi(mi), mi = 1,ncount )
                     pres_read = .TRUE.
                     read_pres(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'pressures'
                     if (iptty .ne. 0) write (iptty, 400) 'pressures'
                  end if
               case ('satu', 'sat ')
                  if (read_sat(1) .and. irdof .ne. 13) then
                     read (iread, *) ( s(mi), mi = 1,ncount )
                     read_sat(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'saturations'
                     if (iptty .ne. 0) write (iptty, 400) 'saturations'
                  end if
               case ('mfra', 'mfr ')
                  if (read_mfrac_iso(1) .and. irdof .ne. 13) then
                     read (iread, *) (frac_gas_iso(mi), mi = 1,ncount)
                     read_mfrac_iso(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'saturations'
                     if (iptty .ne. 0) write (iptty, 400) 'saturations'
                  end if                  
               case ('gasp', 'gas ')
                  if (read_gasp(1)) then
                     read (iread, *) ( pci(mi), mi = 1,ncount )
                     read_gasp(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'gas pressures'
                     if (iptty .ne. 0)write (iptty, 400) 'gas pressures'
                  end if
               case ('poro')
                  if (read_por(1)) then
c gaz debug 110314 set psini to input file values
c                    psini = ps
                     read (iread, *) ( ps(mi), mi = 1,ncount )
                     read_por(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'porosities'
                     if (iptty .ne. 0) write (iptty, 400) 'porosities'
                  end if
c CO2 variables
               case ('co2t')
c CO2 temperature
                  if (read_co2(1)) then
                     read (iread, *) (tco2(mi), mi = 1,ncount )
                     read_co2(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)'co2 temperatures'
                     if (iptty .ne. 0) write (iptty, 400)
     &                    'co2 temperatures'
                     read_co2(2) = .FALSE.
                  end if
                  
               case ('co2p')
c CO2 pressure
                  if (read_co2(1)) then
                     read (iread, *) (phico2(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'co2 pressures'
                     if (iptty .ne. 0)write (iptty, 400) 'co2 pressures'
                     read_co2(2) = .FALSE.
                  end if
               case ('wsat')
c Water saturation
                  if (read_co2(1)) then
                     read (iread, *) (fow(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'co2 water sats'
                     if (iptty .ne. 0)write (iptty, 400)'co2 water sats'
                     read_co2(2) = .FALSE.
                  end if
               case ('lco2')
c Liquid CO2 saturation
                  if (read_co2(1)) then
                     read (iread, *) (fol(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400) 'liquid co2 sats'
                     if (iptty .ne. 0)write(iptty, 400)'liquid co2 sats'
                     read_co2(2) = .FALSE.
                  end if
               case ('diss')
c Dissolved CO2
                  if (read_co2(1)) then
                     read (iread, *) (yc(mi), mi = 1,ncount )
                      do i = 1, ncount
                        if(inico2flg(i).eq.1) then
c Use values from input
                         yc(i) = yc_tmp(i)
                         yoc(i) = yc(i)
                        endif
                      enddo
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400) 'dissolvled co2'
                     if (iptty .ne. 0)write(iptty, 400) 'dissolved co2'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosw')
c Phase-state of water
                  if (read_co2(1)) then
                     read (iread, *) (ieoso(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'phase-state of water'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'phase-state of water'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosc')
c Phase-state of CO2
                  if (read_co2(1)) then
                     read (iread, *) (iceso(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'phase-state of co2'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'phase-state of co2'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosd')
c Array of dissolved CO2 flag
                  if (read_co2(1)) then
                     read (iread, *) (ico2diso(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'check of equilibrium co2 solubility'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'check of equilibrium co2 solubility'
                     read_co2(2) = .FALSE.
                  end if
               case ('mass')
                  if (dumflag(1:9) .ne. 'mass flux') then
                     if (read_mass(1)) then
                        if (.not. allocated(mass_var)) 
     &                       allocate(mass_var(n))
                        read (iread, *) (mass_var(mi), mi = 1,ncount )
                        read_mass(2) = .TRUE.
                     else
                        read (iread, *) ( dummyreal, mi = 1,ncount )
                        if (iout .ne. 0) write (iout, 400) 'masses'
                        if (iptty .ne. 0)write (iptty, 400) 'masses'
                     endif
                  else
                     backspace (iread)
                     exit
                  end if
               case ('pini')
c Initial pressures and temps 
                  if (read_pini(1)) then
                     read(iread, *)  ( tini (mi) , mi=1,ncount )
                     read(iread, *)  ( phini (mi) , mi=1,ncount )
                     read_pini(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'initial pressure and temperature'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'initial pressure and temperature'
                  end if
               case ('xdis')
                  if (read_disx(1)) then
                     read(iread, *)  ( du (mi) , mi=1,ncount )
                     read_disx(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'x displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'x displacement'
                  end if
               case ('ydis')
                  if (read_disy(1)) then
                     read(iread, *)  ( dv (mi) , mi=1,ncount )
                     read_disy(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'y displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'y displacement'
                  end if
               case ('zdis')
                  if (read_disz(1)) then
                     read(iread, *)  ( dw (mi) , mi=1,ncount )
                     read_disz(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'z displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'z displacement'
                  end if
               case ('xstr')
                  if (read_strx(1)) then
                     if(.not.allocated(str_x0)) allocate (str_x0(n0))
                     read(iread, *)  ( str_x0 (mi) , mi=1,ncount )
                     read_strx(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'x stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'x stress'
                  end if
               case ('ystr')
                  if (read_stry(1)) then
                  if(.not.allocated(str_y0)) allocate (str_y0(n0))
                     read(iread, *)  ( str_y0 (mi) , mi=1,ncount )               
                     read_stry(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'y stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'y stress'
                  end if
               case ('xyst')
                  if (read_strxy(1)) then
                  if(.not.allocated(str_xy0))allocate (str_xy0(n0))
                     read(iread, *)  ( str_xy0 (mi) , mi=1,ncount )
                     read_strxy(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'xy stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'xy stress'
                  end if
               case ('zstr')
                  if (read_strz(1)) then
                  if(.not.allocated(str_z0))allocate (str_z0(n0))
                     read(iread, *)  ( str_z0 (mi) , mi=1,ncount )
                     read_strz(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'z stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'z stress'
                  end if
               case ('xzst')
                  if (read_strxz(1)) then
                  if(.not.allocated(str_xz0))allocate (str_xz0(n0))
                     read(iread, *)  ( str_xz0 (mi) , mi=1,ncount )
                     read_strxz(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'xz stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'xz stress'
                  end if
               case ('yzst')
                  if (read_stryz(1)) then
                  if(.not.allocated(str_yz0)) allocate (str_yz0(n0))
                     read(iread, *)  ( str_yz0 (mi) , mi=1,ncount )
                     read_stryz(2) = .TRUE.
                  else
                     read (iread, *) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'yz stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'yz stress'
                  end if
              case default
                  backspace (iread)
                  exit
               end select
            end do
         end if
         read(iread ,'(a11)', end = 100) dumflag
         if (dumflag .ne. 'no fluxes  ') then
            select case (dumflag(1:3))
            case ('all')
               read(iread, *, end = 100) numflux
               read(iread, 6002) (a_axy(i), i = 1,numflux)
               if(jswitch.eq.0) then
                read(iread, 6002) (a_vxy(i), i = 1,numflux)
               else
                read(iread, 6002) (dum1, i = 1,numflux)
               endif
               read_flux(2) = .TRUE.
            case ('liq')
               read(iread, *, end = 100) numflux
               read(iread, 6002) (a_axy(i), i = 1,numflux)
               read_flux(2) = .TRUE.
            case ('mas')
               read(iread, *, end = 100) numflux
               if (size(a_axy) .lt. numflux) then
                  deallocate (a_axy)
                  allocate (a_axy(numflux))
               end if
               read(iread, 6001) (a_axy(i), i = 1,numflux)
               read_flux(2) = .TRUE.
            case ('vap')
               read(iread, *, end = 100) numflux
               if(jswitch.eq.0) then
                read(iread, 6002) (a_vxy(i), i = 1,numflux)
               else
                read(iread, 6002) (dum1, i = 1,numflux)
               endif
               read_flux(2) = .TRUE.
            case default
               backspace iread
            end select
c***  add water table rise subroutine
            if (wtrise_flag) call wtrise
c***  add water table rise subroutine           
            if (reverse_flow) then
               a_axy = -1. * a_axy
               a_vxy = -1. * a_vxy
            end if
         else
            if (.not. compute_flow) then
               write (ierr, *) 'Fluxes not found in restart file',
     &              'while using rflo option'
               if (iptty .ne. 0) then
                  write (iptty, *) 'Fluxes not found in restart file',
     &                 'while using rflo option: STOPPING'
               end if
               stop
            end if
         end if
 6001    format(5g15.8)
 6002    format(4g25.16)
      else
c Reading unformatted file
         read (iread) dumver, dumdate, dumtime, dumtitle
         read (iread)  days
c     
c     read descriptors
c     
         read(iread) wdd1(1:4)
         if (wdd1(1:4) .eq. 'h20 ' .or. wdd1(1:4) .eq. 'wh20' .or.
     &        wdd1(1:4) .eq. 'air ' .or. wdd1(1:4) .eq. 'wair' .or.
     &        wdd1(1:4) .eq. 'ngas' .or. wdd1(1:4) .eq. 'wnga' .or.
     &        wdd1(1:4) .eq. 'carb' .or. wdd1(1:4) .eq. 'wcar') then 
            read(iread) wdd1(5:8)
            read(iread) wdd1(9:12)
            read(iread) wdd1(13:16)
            read(iread) wdd1(17:20)
            ncount=neq
            if(wdd1(13:16).eq.'dpdp') ncount=neq+neq
            if(wdd1(17:20).eq.'dual') ncount=neq+neq+neq
            if (wdd1(1:4).eq.'ngas' .or. wdd1(1:4).eq.'wnga') then
               if(iriver.ne.0 .and. wdd1 .eq. 'ngas') 
     &              ncount = ncount - npoint_riv
               if (read_temp(1)) then
                  read(iread)  ( t   (mi) , mi=1,ncount )
                  read_temp(2) = .TRUE.
               else
                  read(iread)  ( dummyreal , mi=1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'temperatures'
                  if (iptty .ne. 0) write (iptty, 400) 'temperatures'
               end if
               if (read_sat(1)) then
                  read(iread)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.
               else
                  read (iread) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
               if (read_pres(1)) then
                  read(iread)  ( phi (mi) , mi=1,ncount )
                  read_pres(2) = .TRUE.
               else
                  read (iread) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (read_gasp(1)) then
                  read(iread)  ( pci (mi) , mi=1,ncount )
                  read_gasp(2) = .TRUE.
c     
c     check for negative partial pressures
c     
                  do i=1,ncount
                     pci(i)=max(0.0d+0,pci(i))
                  enddo
               else
                  read (iread, *) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'gas pressures'
                  if (iptty .ne. 0)write (iptty, 400) 'gas pressures'
               end if
               if (iriver .ne. 0 .and. wdd1(1:4) .eq. 'ngas') then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_temp(1)) 
     &                    t(i) = t(iriver_con_node(i-ncount))
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
                     if (read_gasp(1)) 
     &                    pci(i) = pci(iriver_con_node(i-ncount))
		  enddo
               end if

            else if(wdd1(1:4).eq.'h20 ' .or. wdd1(1:4).eq.'nh20') then
               if(iriver.ne.0 .and. wdd1(1:4) .eq. 'h20 ') 
     &              ncount = ncount - npoint_riv
               if (read_temp(1) .or. read_co2(1)) then
                  read(iread)  ( t   (mi) , mi=1,ncount )
                  read_temp(2) = .TRUE.
               else
                  read (iread) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'temperatures'
                  if (iptty .ne. 0) write (iptty, 400) 'temperatures'
               end if
               if (read_sat(1) .and. irdof .ne. 13) then
                  read(iread)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.
               else
                  read(iread)  ( dummyreal , mi=1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
               if (read_pres(1) .or. read_co2(1)) then
                  read(iread)  ( phi (mi) , mi=1,ncount )
                  pres_read = .TRUE.
                  read_pres(2) = .TRUE.
               else
                  read (iread) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (wdd1(1:4).eq.'carb' .or. wdd1(1:4) .eq. 'wcar') then
                  read(iread) ( tco2 (mi), mi = 1,ncount )
                  read(iread) ( phico2 (mi), mi = 1,ncount )
                  read(iread) ( fow (mi), mi = 1,ncount )
                  read(iread) ( fog (mi), mi = 1,ncount )
                  read(iread) ( fol (mi), mi = 1,ncount )
                  read(iread) ( ieoso (mi),  mi = 1,ncount )
                  read(iread) ( iceso (mi),  mi = 1,ncount )
                  read(iread) ( ico2diso (mi),  mi = 1,ncount )
                  read_co2(2) = .TRUE.
                  do i = 1, ncount
                     if(inico2flg(i).eq.1) then
c Use values from input
                        fow(i) = fw(i)
                        fol(i) = fl(i)
                        fog(i) = fg(i)
                     endif
                  enddo
                  fw = fow
                  fl = fol
                  fg = fog
                  ieos = ieoso
                  ices = iceso
                  ico2dis = ico2diso
               else
                  if (icarb .eq. 1) then
                     do i = 1, ncount + npoint_riv
                        if(inico2flg(i).ne.1) then
                           fw(i) = 1.d0
                           fg(i) = 0.d0
                           fl(i) = 0.d0
                        endif
                     enddo
                  end if
               end if
               if (iriver .ne. 0 .and. (wdd1(1:4) .eq. 'h20 ' .or. 
     &              wdd1(1:4).eq.'carb')) then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_temp(1) .or. read_co2(1)) 
     &                    t(i) = t(iriver_con_node(i-ncount))
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1) .or. read_co2(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
                     if (wdd1(1:4).eq.'carb') then
                        tco2(i) = tco2(iriver_con_node(i-ncount))
                        phico2(i) = phico2(iriver_con_node(i-ncount))
                        fw(i) = fw(iriver_con_node(i-ncount))
                        fg(i) = fg(iriver_con_node(i-ncount))
                        fl(i) = fl(iriver_con_node(i-ncount))
                        ieoso(i) = ieoso(iriver_con_node(i-ncount))
                        iceso(i) = iceso(iriver_con_node(i-ncount))
                       ico2diso(i) = ico2diso(iriver_con_node(i-ncount))
                     end if
		  enddo
               end if
               if (icarb .eq. 1) then
                  if (wdd1(1:4).eq.'h20 ' .or. wdd1(1:4).eq.'nh20') then
                     tco2 = t
                     phico2 = phi
                     do i = 1, ncount + npoint_riv
                        pl = phico2(i)
                        tl = tco2(i)
                        call co2_properties(1,duma,pl,tl,dum1,
     &                       ices(i),dumb,dumc,duma)
                     enddo
                     iceso = ices
                     ieos = 1
                     ieoso = ieos
				   ico2diso = ico2dis
                  else
                     ieos = ieoso
                     ices = iceso
				   ico2dis = ico2diso
                   end if
                  toco2 = tco2
                  phoco2 = phico2
                  fow = fw
                  fog = fg
                  fol = fl
               end if

            else if(wdd1(1:4).eq.'air ' .or. wdd1(1:4).eq.'nair') then
               if(iriver.ne.0 .and. wdd1(1:4) .eq. 'air ') 
     &              ncount = ncount - npoint_riv
               if (read_sat(1) .and. irdof .ne. 13) then
                  read(iread)  ( s   (mi) , mi=1,ncount )
                  read_sat(2) = .TRUE.
               else
                  read(iread)  ( dummyreal , mi=1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'saturations'
                  if (iptty .ne. 0) write (iptty, 400) 'saturations'
               end if
               if (read_pres(1)) then
                  read(iread)  ( phi (mi) , mi=1,ncount )
                  read_pres(2) = .TRUE.
               else
                  read (iread) ( dummyreal, mi = 1,ncount )
                  if (iout .ne. 0) write (iout, 400) 'pressures'
                  if (iptty .ne. 0) write (iptty, 400) 'pressures'
               end if
               if (iriver .ne. 0 .and. wdd1(1:4) .eq. 'air ') then
                  do i = ncount+1,ncount+npoint_riv
                     j = iriver_con_node(i-ncount)
                     if (read_sat(1)) 
     &                    s(i) = s(iriver_con_node(i-ncount))
                     if (read_pres(1)) 
     &                    phi(i) = phi(iriver_con_node(i-ncount))
		  enddo
               end if
            endif

            dumflag = '           '
            read(iread , end = 100) dumflag

         else
            rerr = 0
            rewind (iread)
            read (iread) dumver, dumdate, dumtime, dumtitle
            read (iread)  days
            read (iread) ncount, geom_type
            modneq = mod(ncount,neq)
            if (modneq .ne. 0) goto 2000
            do
               dumflag = '           '
               read (iread , end = 100) dumflag
               select case (dumflag(1:4))
               case ('temp')
                  if (read_temp(1) .or. read_co2(1)) then
                     read (iread) ( t(mi), mi = 1,ncount )
                     read_temp(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'temperatures'
                     if (iptty .ne. 0) write (iptty, 400) 'temperatures'
                  end if
               case ('pres')
                  if (read_pres(1) .or. read_co2(1)) then
                     read (iread) ( phi(mi), mi = 1,ncount )
                     pres_read = .TRUE.
                     read_pres(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'pressures'
                     if (iptty .ne. 0) write (iptty, 400) 'pressures'
                  end if
               case ('satu', 'sat ')
                  if (read_sat(1) .and. irdof .ne. 13) then
                     read (iread) ( s(mi), mi = 1,ncount )
                     read_sat(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'saturations'
                     if (iptty .ne. 0) write (iptty, 400) 'saturations'
                  end if
               case ('gasp', 'gas ')
                  if (read_gasp(1)) then
                     read (iread) ( pci(mi), mi = 1,ncount )
                     read_gasp(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'gas pressures'
                     if (iptty .ne. 0)write (iptty, 400) 'gas pressures'
                  end if
               case ('poro')
                  if (read_por(1)) then
c gaz debug 110314 set initial porosity to input file values
c                    psini = ps
                     read (iread) ( ps(mi), mi = 1,ncount )
                     read_por(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'porosities'
                     if (iptty .ne. 0) write (iptty, 400) 'porosities'
                  end if
c CO2 variables
               case ('co2t')
c CO2 temperature
                  if (read_co2(1)) then
                     read (iread) (tco2(mi), mi = 1,ncount )
                     read_co2(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)'co2 temperatures'
                     if (iptty .ne. 0) write (iptty, 400)
     &                    'co2 temperatures'
                     read_co2(2) = .FALSE.
                  end if
                  
               case ('co2p')
c CO2 pressure
                  if (read_co2(1)) then
                     read (iread) (phico2(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'co2 pressures'
                     if (iptty .ne. 0)write (iptty, 400) 'co2 pressures'
                     read_co2(2) = .FALSE.
                  end if
               case ('wsat')
c Water saturation
                  if (read_co2(1)) then
                     read (iread) (fow(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'co2 water sats'
                     if (iptty .ne. 0)write (iptty, 400)'co2 water sats'
                     read_co2(2) = .FALSE.
                  end if
               case ('lco2')
c Liquid CO2 saturation
                  if (read_co2(1)) then
                     read (iread) (fol(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400) 'liquid co2 sats'
                     if (iptty .ne. 0)write(iptty, 400)'liquid co2 sats'
                     read_co2(2) = .FALSE.
                  end if
               case ('diss')
c Dissolved CO2
                  if (read_co2(1)) then
                     read (iread) (yc(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400) 'dissolvled co2'
                     if (iptty .ne. 0)write(iptty, 400) 'dissolved co2'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosw')
c Phase-state of water
                  if (read_co2(1)) then
                     read (iread, *) (ieoso(mi), mi = 1,ncount )
                  else
                     read (iread, *) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'phase-state of water'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'phase-state of water'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosc')
c Phase-state of CO2
                  if (read_co2(1)) then
                     read (iread) (iceso(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'phase-state of co2'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'phase-state of co2'
                     read_co2(2) = .FALSE.
                  end if
               case ('eosd')
c Phase-state of CO2
                  if (read_co2(1)) then
                     read (iread) (ico2diso(mi), mi = 1,ncount )
                  else
                     read (iread) ( dummyint, mi = 1,ncount )
                     if (iout .ne. 0) write(iout, 400)
     &                    'check of equilibrium co2 solubility'
                     if (iptty .ne. 0)write(iptty, 400)
     &                    'check of equilibrium co2 solubility'
                     read_co2(2) = .FALSE.
                  end if

               case ('mass')
                  if (read_mass(1)) then
                     read (iread) (mass_var(mi), mi = 1,ncount )
                     read_mass(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 'masses'
                     if (iptty .ne. 0)write (iptty, 400) 'masses'
                  endif
               case ('pini')
c Initial pressures and temps 
                  if (read_pini(1)) then
                     read(iread)  ( tini (mi) , mi=1,ncount )
                     read(iread)  ( phini (mi) , mi=1,ncount )
                     read_pini(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400)  
     &                    'initial pressure and temperature'
                     if (iptty .ne. 0)write (iptty, 400)  
     &                    'initial pressure and temperature'
                  end if
               case ('xdis')
                  if (read_disx(1)) then
                     read(iread)  ( du (mi) , mi=1,ncount )
                     read_disx(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'x displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'x displacement'
                  end if
               case ('ydis')
                  if (read_disy(1)) then
                     read(iread)  ( dv (mi) , mi=1,ncount )
                     read_disy(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'y displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'y displacement'
                  end if
               case ('zdis')
                  if (read_disz(1)) then
                     read(iread)  ( dw (mi) , mi=1,ncount )
                     read_disz(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'z displacement'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'z displacement'
                  end if
               case ('xstr')
                  if (read_strx(1)) then
                     read(iread)  ( str_x0 (mi) , mi=1,ncount )
                     read_strx(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'x stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'x stress'
                  end if
               case ('ystr')
                  if (read_stry(1)) then
                     read(iread)  ( str_y0 (mi) , mi=1,ncount )
                     read_stry(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'y stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'y stress'
                  end if
               case ('xyst')
                  if (read_strxy(1)) then
                     read(iread)  ( str_xy0 (mi) , mi=1,ncount )
                     read_strxy(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'xy stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'xy stress'
                  end if
               case ('zstr')
                  if (read_strz(1)) then
                     read(iread)  ( str_z0 (mi) , mi=1,ncount )
                     read_strz(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'z stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'z stress'
                  end if
               case ('xzst')
                  if (read_strxz(1)) then
                     read(iread)  ( str_xz0 (mi) , mi=1,ncount )
                     read_strxz(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'xz stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'xz stress'
                  end if
               case ('yzst')
                  if (read_stryz(1)) then
                     read(iread)  ( str_yz0 (mi) , mi=1,ncount )
                     read_stryz(2) = .TRUE.
                  else
                     read (iread) ( dummyreal, mi = 1,ncount )
                     if (iout .ne. 0) write (iout, 400) 
     &                    'yz stress'
                     if (iptty .ne. 0)write (iptty, 400) 
     &                    'yz stress'
                  end if
               case default
                  exit
               end select
            end do
         end if
         if (dumflag .ne. 'no fluxes  ') then
            read(iread) numflux
            if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &           'liq') read(iread) (a_axy(i), i = 1,numflux)
            if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &           'vap') then
               if(jswitch.eq.0) then
                read(iread) (a_vxy(i), i = 1,numflux)
               else
                read(iread) (dum1, i = 1,numflux)
               endif
            endif
            if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &           'liq' .or. flux_flag(1:3) .eq. 'vap')
     &           read_flux(2) = .TRUE.
            if (reverse_flow) then
               a_axy = -1. * a_axy
               a_vxy = -1. * a_vxy
            end if
         else
            if (.not. compute_flow) then
               if (dumflag .eq. 'no fluxes  ' .or. dumflag .eq. '') then
                  write (ierr, *) 'Fluxes not found in restart file',
     &                 'while using rflo option'
                  if (iptty .ne. 0) then
                     write (iptty, *) 'Fluxes not found in restart ',
     &                    'file while using rflo option: STOPPING'
                  end if
                  stop
               end if
            end if
         end if
      end if
 400  format ('Not using ', a, ' found in restart file')
 500  format (a, ' not found in restart file')
 600  format (a, ' found in restart file will be used')
 100  continue
      if (icarb .eq. 1 .and. .not. read_co2(2)) then
c Set CO2 parameters using t and phi
         tco2 = t
         phico2 = phi
         do i = 1, ncount + npoint_riv
            pl = phico2(i)
            tl = tco2(i)
            call co2_properties(1,duma,pl,tl,dum1,
     &           ices(i),dumb,dumc,duma)
         enddo
         iceso = ices
         ieos = 1
         ieoso = ieos
         ico2diso = ico2dis
         do i = 1, ncount + npoint_riv
            if(inico2flg(i).ne.1) then
               fw(i) = 1.d0
               fg(i) = 0.d0
               fl(i) = 0.d0
            endif
         enddo
         fow = fw
         fog = fg
         fol = fl
         toco2 = tco2
         phoco2 = phico2
      end if
      if (read_co2(2)) then
         co2_read = .TRUE.
         do i = 1, n0
            if(inico2flg(i).eq.1) then
               fow(i) = fw(i)
               fol(i) = fl(i)
               fog(i) = fg(i)
            else
               fog(i) = max(1.d0-fow(i)-fol(i), 0.d0)
            end if
         end do
         fw = fow
         fl = fol
         fg = fog
         toco2 = tco2
         phoco2 = phico2
         ieos = ieoso
         ices = iceso
         ico2dis = ico2diso
         if (.not. read_temp(2)) then
            t = tco2
            to = tco2
         end if
         if (.not. read_pres(2)) then
            phi = phico2
            pho = phico2
         end if
      end if
      if (read_temp(1) .and. .not. read_temp(2)) then
         if (iout .ne. 0) write (iout, 500) 'Temperatures'
         if (iptty .ne. 0) write (iptty, 500) 'Temperatures'
      end if
      if (read_pres(1) .and. .not. read_pres(2)) then
         if (iout .ne. 0) write (iout, 500) 'Pressures'
         if (iptty .ne. 0) write (iptty, 500) 'Pressures'
      end if
      if (read_pres(1) .and. read_mass(1)) then
         if (read_mass(2)) then
c Mass option currently not supported, swapped formats 400/600
c   and set mass_read flag to false
            mass_read = .FALSE.
            if (iout .ne. 0) write (iout, 400) 'Masses'
            if (iptty .ne. 0) write (iptty, 400) 'Masses'
            if (read_pres(2)) then
               if (iout .ne. 0) write (iout, 600) 'Pressures'
               if (iptty .ne. 0) write (iptty, 600) 'Pressures'
            end if
         else if (read_pres(2)) then
            if (iout .ne. 0) write (iout, 600) 'Pressures'
            if (iptty .ne. 0) write (iptty, 600) 'Pressures'
         end if
      end if
      if (read_mass(1) .and. .not. read_mass(2)) then
         if (iout .ne. 0) write (iout, 500) 'Masses'
         if (iptty .ne. 0) write (iptty, 500) 'Masses'
         mass_read = .FALSE.
      end if
      if (read_gasp(1) .and. .not. read_gasp(2)) then
         if (iout .ne. 0) write (iout, 500) 'Gas pressures'
         if (iptty .ne. 0) write (iptty, 500) 'Gas pressures'
      end if
      if (read_sat(1) .and. .not. read_sat(2)) then
         if (iout .ne. 0) write (iout, 500) 'Saturations'
         if (iptty .ne. 0) write (iptty, 500) 'Saturations'
      end if
      if (read_co2(1) .and. .not. read_co2(2)) then
         if (iout .ne. 0) write (iout, 500) 'CO2 parameters'
         if (iptty .ne. 0) write (iptty, 500) 'CO2 parameters'         
      end if
      if (read_flux(1) .and. .not. read_flux(2)) then
         if (iout .ne. 0) write (iout, 500) 'Fluxes'
         if (iptty .ne. 0) write (iptty, 500) 'Fluxes'
      end if

      nc=1         
c     
c     sort out multiple grids
c 
      if(idpdp.ne.0) then
         if(wdd1(13:16).eq.'dpdp') then
            nc=2
            nc1=0
            nc2=neq
            nd1=0
            nd2=neq
         else if(wdd1(17:20).eq.'dual') then
            nc=2
            nc1=0
            nc2=neq
            nd1=0
            nd2=neq+neq
         else
            nc=2
            nc1=0
            nc2=neq
            nd1=0
            nd2=0
         endif
      else if(idualp.ne.0) then
         if(wdd1(13:16).eq.'dpdp') then
            nc=3
            nc1=0
            nc2=neq
            nc3=neq+neq
            nd1=0
            nd2=neq
            nd3=neq
         else if(wdd1(17:20).eq.'dual') then
            nc=3
            nc1=0
            nc2=neq
            nc3=neq+neq
            nd1=0
            nd2=neq
            nd3=neq+neq
         else
            nc=3
            nc1=0
            nc2=neq
            nc3=neq+neq
            nd1=0
            nd2=0
            nd3=0
         endif
      else
         if(wdd1(13:16).eq.'dpdp') then
            nc=1
            nc1=0
            nd1=neq
         else if(wdd1(17:20).eq.'dual') then
            nc=1
            nc1=0
            nd1=neq+neq
         else
            nc=1
            nc1=0
            nd1=0
         endif
      endif
      do j=1,nc
         if(j.eq.1) nx=nc1
         if(j.eq.2) nx=nc2
         if(j.eq.3) nx=nc3
         if(j.eq.1) ny=nd1
         if(j.eq.2) ny=nd2
         if(j.eq.3) ny=nd3

         do ii=1,neq
            i=(j-1)*neq +ii
            if(ieos(i).ge.0) then
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  satr = s(ii+ny)
               else
                  satr = sat_dum
               end if
               if (idof .eq. 1) then
c                  if (s(ii+ny) .lt. 1.0) then
! heat conduction only
                  if(ico2.ge.0) to(ii+nx)=t(ii+ny) 
c                  end if
               else
                  if(ico2.ge.0) to(ii+nx)=t(ii+ny)
                  pho(ii+nx)=phi(ii+ny)
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     pcio(ii+nx)=pci(ii+ny)
                     if(s(ii+ny).le.tolw) then
                        ieos(i)=3
                        so(ii+nx)=0.0 
                     else if(s(ii+ny).ge.1.0) then
                        ieos(i)=1
                        so(ii+nx)=1.0 
                     else
                        ieos(i)=2
                        so(ii+nx)=s(ii+ny)
                     endif
                  endif
               endif
            endif
         enddo
      enddo

      nsave  =  0
      
      if(wdd1(5:8).eq.'trac' .or. read_trac(1)) then
         if (iccen .eq. 1) then
            call diskc(read_trac(1))
         else if (read_trac(1)) then
            write (ierr, 1002)
            if (iout .ne. 0) write (iout, 1002)
            if (iptty .ne. 0) write (iptty, 1002)
         else
            write (ierr, 1000)
            if (iout .ne. 0) write (iout, 1000)
            if (iptty .ne. 0) write (iptty, 1000)
        end if
      end if
    
      if(wdd1(5:8).eq.'ptrk' .or. read_ptrk(1)) then
         if (ptrak) then
           call diskp(read_ptrk(1))
         else if (read_ptrk(1)) then
            write (ierr, 1003)
            if (iout .ne. 0) write (iout, 1003)
            if (iptty .ne. 0) write (iptty, 1003)
         else
            write (ierr, 1001)
            if (iout .ne. 0) write (iout, 1001)
            if (iptty .ne. 0) write (iptty, 1001)
        end if
      end if     
c      if(wdd1(9:12).eq.'strs') call  stress  ( 4 )
c**** save initial pressures and temps to restart file ****
      if (read_pini(1) .and. .not. read_pini(2) ) then
         if (cform(6) .eq. 'formatted') then
            read(iread, *, END=3000)  ( tini (mi) , mi=1,ncount )
            read(iread, *, END=3000)  ( phini (mi) , mi=1,ncount )
         else
            read(iread, END=3000)  ( tini (mi) , mi=1,ncount )
            read(iread, END=3000)  ( phini (mi) , mi=1,ncount )
         endif
      endif     
c*** EK overwrite initial dissolved CO2 concentrations with tracer output
        if(icarb.eq.1) then
        do i=1,ncount
        if(inico2flg(i).eq.2) then
          if(.not.read_trac(1)) then
           write(ierr,*) 'no tracer data in restart file'
           write(ierr,*) 'incompatable with co2frac inico2flg=2'
                                stop
          endif
          if(carbon_tracer.eq.0) then
           write(ierr,*) 'component H2CO3 required for carbon coupling'
                                stop
          endif

 666            format(3i10,e10.3)
                yc(i)=anl((carbon_tracer-1)*ncount+i)*44E-3
        endif
        end do
        endif
 1000 format (1x, 'Tracer data found in restart file for non-trac',
     &     ' problem, data will not be used')
 1001 format (1x, 'Particle tracking data found in restart file for',
     &     ' non-ptrk problem, data will not be used')
 1002 format (1x, 'Tracer data requested for non-trac problem')
 1003 format (1x, 'Particle tracking data requested for non-ptrk',
     &     ' problem')
      close (iread)
      return

 2000 continue
      if (geom_type .ne. 'dual' .or. geom_type .ne. 'dpdp' .or.
     &     geom_type .ne. 'gdpm') geom_type = ''
      write (ierr, 1004) geom_type, ncount, neq
      if (iout .ne. 0) write (iout, 1004) geom_type, ncount, neq
      if (iptty .ne. 0) write (iptty, 1004) geom_type, ncount, neq
 1004 format(1x,a, 'Problem has ', i7, ' data elements for grid with ', 
     &     i10, ' nodes', /, 'Stopping')
      stop

 3000 continue
      if (iout .ne. 0) write (iout, 500) 
     &     'Initial pressures and temperatures'
      if (iptty .ne. 0) write (iptty, 500) 
     &     'Initial pressures and temperatures'
      stop

      end
      
