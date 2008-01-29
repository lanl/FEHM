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
      use comdi
      use comdti
      use comfi
      use comflow
      use comii
      use compart
      use comxi
      use davidi
      implicit none

      real*8 dummyreal,tolw,sat_dum,satr
      parameter(tolw=1.d-99,sat_dum= 1.d00)
      integer jx,mi,ncount,nc,nc1,nc2,nc3,nd1,nd2,nd3,j,nx,ny,ii,i
      character*80 dumtitle
      character*30 dumver
      character*11 dumdate, dumflag
      character*8 dumtime

      inquire (iread, opened = ex)
      if (ex) close (iread)
      open (iread,file=nmfil(6),status=cstats(6),form=cform(6))
           
      rewind  iread
      if (cform(6) .eq. 'formatted') then
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
         read(iread,'(a4)') wdd1(5:8)
         read(iread,'(a4)') wdd1(9:12)
         read(iread,'(a4)') wdd1(13:16)
         read(iread,'(a4)') wdd1(17:20)
         ncount=neq
         if(wdd1(13:16).eq.'dpdp') ncount=neq+neq
         if(wdd1(17:20).eq.'dual') ncount=neq+neq+neq
         if (wdd1(1:4).eq.'ngas') then
            read(iread ,*)  ( t   (mi) , mi=1,ncount )
            read(iread ,*)  ( s   (mi) , mi=1,ncount )
            read(iread ,*)  ( phi (mi) , mi=1,ncount )
            read(iread ,*)  ( pci (mi) , mi=1,ncount )
c     
c     check for negative partial pressures
c     
            do i=1,ncount
               pci(i)=max(0.0d+0,pci(i))
            enddo

         else if(wdd1(1:4).eq.'h20 ') then
            read(iread ,*)  ( t   (mi) , mi=1,ncount )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               read(iread ,*)  ( s   (mi) , mi=1,ncount )
            else
               read(iread ,*)  ( dummyreal , mi=1,ncount )
            end if
            read(iread ,*)  ( phi (mi) , mi=1,ncount )
         else if(wdd1(1:4).eq.'air ') then
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               read(iread ,*)  ( s   (mi) , mi=1,ncount )
            else
               read(iread ,*)  ( dummyreal , mi=1,ncount )
            end if
            read(iread ,*)  ( phi (mi) , mi=1,ncount )
         endif
         read(iread ,'(a11)', end = 100) dumflag
         if (dumflag .ne. 'no fluxes  ') then
            select case (dumflag(1:3))
            case ('all')
               read(iread, *, end = 100) numflux
               read(iread, 6002) (a_axy(i), i = 1,numflux)
               read(iread, 6002) (a_vxy(i), i = 1,numflux)
            case ('liq')
               read(iread, *, end = 100) numflux
               read(iread, 6002) (a_axy(i), i = 1,numflux)
            case ('mas')
               read(iread, *, end = 100) numflux
               if (size(a_axy) .lt. numflux) then
                  deallocate (a_axy)
                  allocate (a_axy(numflux))
               end if
               read(iread, 6001) (a_axy(i), i = 1,numflux)
            case ('vap')
               read(iread, *, end = 100) numflux
               read(iread, 6002) (a_vxy(i), i = 1,numflux)
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
         read (iread) dumver, dumdate, dumtime, dumtitle
         read (iread)  days
c     
c     read descriptors
c     
         read(iread) wdd1(1:4)
         read(iread) wdd1(5:8)
         read(iread) wdd1(9:12)
         read(iread) wdd1(13:16)
         read(iread) wdd1(17:20)
         ncount=neq
         if(wdd1(13:16).eq.'dpdp') ncount=neq+neq
         if(wdd1(17:20).eq.'dual') ncount=neq+neq+neq
         if (wdd1(1:4).eq.'ngas') then
            read(iread)  ( t   (mi) , mi=1,ncount )
            read(iread)  ( s   (mi) , mi=1,ncount )
            read(iread)  ( phi (mi) , mi=1,ncount )
            read(iread)  ( pci (mi) , mi=1,ncount )
c     
c     check for negative partial pressures
c     
            do i=1,ncount
               pci(i)=max(0.0d+0,pci(i))
            enddo

         else if(wdd1(1:4).eq.'h20 ') then
            read(iread)  ( t   (mi) , mi=1,ncount )
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               read(iread)  ( s   (mi) , mi=1,ncount )
            else
               read(iread)  ( dummyreal , mi=1,ncount )
            end if
            read(iread)  ( phi (mi) , mi=1,ncount )
         else if(wdd1(1:4).eq.'air ') then
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               read(iread)  ( s   (mi) , mi=1,ncount )
            else
               read(iread)  ( dummyreal , mi=1,ncount )
            end if
            read(iread)  ( phi (mi) , mi=1,ncount )
         endif
         dumflag = '           '
         read(iread , end = 100) dumflag
         if (dumflag .ne. 'no fluxes  ') then
            read(iread) numflux
            if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &           'liq') read(iread) (a_axy(i), i = 1,numflux)
            if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &           'vap') read(iread) (a_vxy(i), i = 1,numflux)
            if (reverse_flow) then
               a_axy = -1. * a_axy
               a_vxy = -1. * a_vxy
            end if
         end if
      end if
c     
c     sort out multiple grids
c 
 100  continue
      if (.not. compute_flow) then
         if (dumflag .eq. 'no fluxes  ' .or. dumflag .eq.  '') then
            write (ierr, *) 'Fluxes not found in restart file',
     &           'while using rflo option'
            if (iptty .ne. 0) then
               write (iptty, *) 'Fluxes not found in restart file',
     &              'while using rflo option: STOPPING'
            end if
            stop
         end if
      end if
      nc=1         
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
                  if (s(ii+ny) .lt. 1.0) then
! heat conduction only
                     if(ico2.ge.0) to(ii+nx)=t(ii+ny) 
                  end if
               else
                  if(ico2.ge.0) to(ii+nx)=t(ii+ny)
                  pho(ii+nx)=phi(ii+ny)
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     pcio(ii+nx)=pci(ii+ny)
                     if(s(ii+ny).le.0.0) then
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
      
      if(wdd1(5:8).eq.'trac') then
         if (iccen .eq. 1) then
            call diskc
         else
            write (ierr, 1000)
            if (iout .ne. 0) write (iout, 1000)
            if (iptty .ne. 0) write (iptty, 1000)
        end if
      end if     
      if(wdd1(5:8).eq.'ptrak') then
         if (ptrak) then
           call diskp 
         else
            write (ierr, 1001)
            if (iout .ne. 0) write (iout, 1001)
            if (iptty .ne. 0) write (iptty, 1001)
        end if
      end if     
      if(wdd1(9:12).eq.'strs') call  stress  ( 4 )

 1000 format (1x, 'Tracer data found in restart file for non-trac',
     &     ' problem, data will not be used')
 1001 format (1x, 'Particle tracking  data found in restart file for ',
     &     ' non-ptrk problem, data will not be used')
      close (iread)
      return
      end
      
