      subroutine steady(iflg,inflow,inflowe)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
CD1 Determime if steady state solution achieved                        
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 18-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/steady.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Feb 05 11:30:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:55:52   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:18   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   isave                    I/O  Unit number for restart file (to write)
CD3   
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   am0             REAL*8   faar   Initial mass in problem
CD4   ame             REAL*8   faar   ?
CD4   astmo           REAL*8   faar   Initial steam mass in problem
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   deneh           REAL*8   fdd    Last time step energy accumulation term
CD4                                     at each node
CD4   denei           REAL*8   fcc    Energy accumulation term
CD4   denej           REAL*8   fdd    Last time step energy accumulation time
CD4                                     derivative at each node
CD4   denh            REAL*8   fdd    Last time step mass accumulation term at
CD4                                     each node
CD4   deni            REAL*8   fcc    Mass accumulation term
CD4   denj            REAL*8   fdd    Last time step mass accumulation time
CD4                                     derivative at each node
CD4   dstm            REAL*8   fcc    Steam mass
CD4   dtot            REAL*8   faar   Current time step size in seconds
CD4   eflow           REAL*8   fdd    Energy flow at each source node
CD4   grav            REAL*8   faar   Value of gravity
CD4   icnl            INT      faai   Problem dimension
CD4   isave           INT      faai   Unit number for restart file (to write)
CD4   ka              INT      fbb    Contains boundary type information for
CD4                                     each node
CD4   n               INT      faai   Total number of nodes
CD4   pflow           REAL*8   fdd    Flowing pressure at each source node
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   phini           REAL*8   fdd2   Initial pressure at each node
CD4   pho             REAL*8   fdd    Last time step pressure at each node
CD4   pnx             REAL*8   fdd    ?Permeability in the x-direction, liquid
CD4                                     velocity in the x-direction, vapor
CD4                                     velocity in the x-direction
CD4   pny             REAL*8   fdd    ?Permeability in the y-direction liquid
CD4                                     velocity in the y direction, vapor
CD4                                     velocity in the y-direction
CD4   pnz             REAL*8   fdd    ?Permeability in the z-direction, liquid
CD4                                     velocity in the z-direction, vapor
CD4                                     velocity in the z-direction
CD4   ps              REAL*8   fdd    Porosity at each node
CD4   psini           REAL*8   fdd2   Initial porosity at each node
CD4   qflux           REAL*8   fdd    Heat flux at each node
CD4   qflxm           REAL*8   fdd    Heat flux impedance at each node
CD4   qh              REAL*8   fdd    Energy source term at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   so              REAL*8   fdd    Last time step saturation at each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4   thx             REAL*8   fdd    Thermal conductivity x-direction
CD4   thy             REAL*8   fdd    Thermal conductivity y-direction
CD4   thz             REAL*8   fdd    Thermal conductivity z-direction
CD4   tini            REAL*8   fdd2   Initial temperature at each node
CD4   to              REAL*8   fdd    Last time step temperature at each node
CD4   volume          REAL*8   fdd    Volume associated at each node
CD4   wellim          REAL*8   fdd    Well impedance at each source node
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   bnswer                   Initiate thermodynamics and solution routines
CD4   enthp           REAL*8   Calculate enthalpy as a function of pressure
CD4                              and temperature   
CD4                              
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   cztop           REAL*8   Maximum coordinate 
CD5   etoll           REAL*8   ?Coordinate tolerance
CD5   i               INT      Loop index
CD5   ir              INT      Problem dimension
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This routine will be used in the future but is not currently 
CD8 used by the code.
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 Not Applicable.  See Special Comments.
CD9
C***********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comfi
      use comsplitts
      use comsteady
      use davidi
      implicit none

      integer i, iflg, itt
      real*8 divisor, divisore, inflow, inflowe, btol_save, tolde_save
      real*8 flow_max, ratio2, balfac1, balfac2, balfac3, balfac12
      real*8 balfac13, tol_divisor, tmch1_min, tmch1_max
      character*20 dummy
      logical null1
      parameter (tol_divisor= 1.d-8)

      if(isteady.eq.0) return 
      if(iflg .eq. 0) then
c     Reading steady input parameters
         toldh = tolerance
         toldp = tolerance
         tolds = tolerance
         toldt = tolerance
         toldc = tolerance
         tolde = tolerance
         tacc = tolerance
         tol_str = tolerance
         spercent = .false.
         divisor = 1.
         divisore = 1.
         smult = 1.
         sday = 0.
         shtl = 0.0
         sdmx = 0.0
         stmch = 0.0
         balance_tol = tolerance
         tmch_old = 1.d20
         svar_flag = .false.
         sflux_flag = .false.
         info_string = 'absolute difference'
         time_ss = 1.d20
         sminstep = 2
         if (ifree .ne. 0)  then
            snstep = 150
            smult = 2.0
         else if (ico2 .lt. 0 ) then
            snstep = 500
            smult = 2.0
         else if (ico2 .eq. 0) then
            snstep = 200
            smult = 2.0
         else if (ico2 .gt. 0) then
            snstep = 200
            smult = 1.5
         end if

         do
            read (inpt, '(a80)') wdd1
            if ( null1(wdd1) .or. wdd1(1:3) .eq. 'end' ) exit
            select case (wdd1(1:4))
            case ('shea')
               read (wdd1, *) dummy, toldh
               toldp = toldh * 997. * 9.8d-6
               svar_flag = .true.
            case('spre')
               read (wdd1, *) dummy, toldp
               svar_flag = .true.
            case ('stem')
               read (wdd1, *) dummy, toldt
               svar_flag = .true.
            case('ssat')
               read (wdd1, *) dummy, tolds
               svar_flag = .true.
            case('sair')
               read (wdd1, *) dummy, toldc
               svar_flag = .true.
            case('sflu')
               read (wdd1, *) dummy, balance_tol
               sflux_flag = .true.
            case('sent')
               read (wdd1, *) dummy, tolde
               sflux_flag = .true.             
            case ('stim')
               read (wdd1, *) dummy, time_ss
            case ('smul')
               read (wdd1, *) dummy, smult
            case('sday')
               read (wdd1, *) dummy, sday
            case('sdmx')
               read (wdd1, *) dummy, sdmx               
            case('snst')
               read (wdd1, *) dummy, snstep
            case('smst')
               read (wdd1, *) dummy, sminstep
            case('shtl')
               read (wdd1, *) dummy, shtl
            case('sacc')
               read (wdd1, *) dummy, tacc
            case('stmc')
               read (wdd1, *) dummy, stmch
            case('smas')
               read (wdd1, *) dummy, smass  
            case('tstr')
               read (wdd1, *) dummy, tol_str                 
            case('stss')
               read (wdd1, *) dummy, stsstr                
            case('sper')
               spercent = .true.
               info_string = 'fractional change'
            end select
         end do
      
      else if(iflg.eq.-1) then
c 'steady macro overrides usual stopping time tims (gaz 032405)
       tims_trans = tims
       daymax_save = daymax
       tims = time_ss
       if(sdmx.gt.0.0) then
        daymax = sdmx
       else 
        daymax = tims
       endif
	 daycs_save=daycs
	 daycf_save=daycf
	 dayhf_save=dayhf
	 dayhs_save=dayhs 
       dayhf = 1.d50
	 dayhs = 1.d50
c gaz 033020       
       daycs = 1.d50
       daycf = 1.d50
	 day_save_ss = day
         if (sday .ne. 0) day = sday
         if (shtl .ne. 0.0) shtl0 = head_tol
         if (stmch .ne. 0.0) stmch0 = tmch
         nstep_save = nstep
         nstep = snstep
         aiaa_save = aiaa
         aiaa = smult
c gaz 010222 added  tmch_save
         tmch_save = tmch

      else if(iflg.eq.1) then
        divisor = 1.
        divisore = 1.
        i_pdiff = 0
        pmax_i = 0.
         if (svar_flag) then
            isteady=-1
         else
            isteady = 1
            return
         endif
         if(irdof.eq.13) then
            pdifmax = 0.0
            accmax = 0.0
            accdif_i = 0
            do i = 1, n
               if (spercent) then
                  if (pho(i) .ne. 0.) then
                     divisor = pho(i)
                  else if (phi(i) .ne. 0.) then
                     divisor = phi(i)
                  else
                     divisor = 1.
                  end if
               end if
               pdifmax = max(abs((phi(i)-pho(i))/divisor),pdifmax)
c 082319 gaz                
               accdif_i = abs((deni(i)*sx1(i))*dtot)
               if(accdif_i.gt.accmax) then
                i_accdif = i
                accmax = accdif_i
               endif  
            end do
            if(pdifmax.gt.toldp) go to 100 
            if(accmax.gt.tacc) go to 100     
         else if(ico2.lt.0) then
            pdifmax = 0.0
            sdifmax = 0.0
            accmax = 0.0
            amass = 0.0
            i_pdiff = 0
            i_sdiff = 0
            i_accdif = 0
            do i = 1, n
               if (spercent) then
                  if (pho(i) .ne. 0) then
                     divisor = pho(i)
                  else if (phi(i) .ne. 0) then
                     divisor = phi(i)
                  else
                     divisor = 1.
                  end if
               end if
               pdiff_i = abs((phi(i)-pho(i))/divisor)
               if(pdiff_i.gt.pdifmax) then
                i_pdiff = i
                pmax_i = phi(i)
                pmax_io = pho(i)
                pdifmax = pdiff_i
               endif
               if (spercent) then
                     divisor = 1.
               end if
               sdiff_i = abs((s(i)-so(i))/1.)
               if(sdiff_i.gt.sdifmax) then
                i_sdiff = i
                smax_i = s(i)
                smax_io = so(i)
                sdifmax = sdiff_i
               endif 
c 082319 gaz                
               accdif_i = abs((deni(i)*sx1(i))*dtot)
               if(accdif_i.gt.accmax) then
                i_accdif = i
                accmax = accdif_i
               endif                
            amass = amass + denh(i) * volume(i)
            end do
            if(l.eq.1) then
             amass0 = am0
             amass_ch =(amass-amass0)
             amass0 = amass
            else
             amass_ch =(amass-amass0)
             amass0 = amass
            endif
            if(pdifmax.gt.toldp) go to 100
            if(sdifmax.gt.tolds) go to 100 
            if(accmax.gt.tacc) go to 100   
         else if(ico2.eq.0) then   
            pdifmax = 0.0
            sdifmax = 0.0
            tdifmax = 0.0
            i_accdif = 0
            accmax = 0.0
            amass = 0
            do i = 1, n
               if (spercent) then
                  if (pho(i) .ne. 0.) then
                     divisor = pho(i)
                  else if (phi(i) .ne. 0.) then
                     divisor = phi(i)
                  else
                     divisor = 1.
                  end if
               end if
               pdifmax = max(abs((phi(i)-pho(i))/divisor),pdifmax) 
               if (spercent) then
                     divisor = 1.
               end if
               sdifmax = max(abs((s(i)-so(i))/1.),sdifmax)   
               if (spercent) then
                  if (to(i) .ne. 0.) then
                     divisor = to(i)
                  else if (t(i) .ne. 0.) then
                     divisor = t(i)
                  else
                     divisor = 1.
                  end if
               end if
               tdifmax = max(abs((t(i)-to(i))/divisor),tdifmax)
c 082319 gaz                
               accdif_i = abs((deni(i)*sx1(i))*dtot)
               if(accdif_i.gt.accmax) then
                i_accdif = i
                accmax = accdif_i
               endif                
              amass = amass + denh(i) * volume(i)
            end do
            if(l.eq.1) then
             amass0 = am0
             amass_ch =(amass-amass0)
             amass0 = amass
            else
             amass_ch =(amass-amass0)
             amass0 = amass
            endif
            if(pdifmax.gt.toldp) go to 100   
            if(sdifmax.gt.tolds) go to 100   
            if(tdifmax.gt.toldt) go to 100  
            if(accmax.gt.tacc) go to 100     
         else if(ico2.gt.0) then
            pdifmax = 0.0
            sdifmax = 0.0
            tdifmax = 0.0
            pcidifmax = 0.0
            accmax = 0.0
            do i = 1, n
               if (spercent) then
                  if (pho(i) .ne. 0.) then
                     divisor = pho(i)
                  else if (phi(i) .ne. 0.) then
                     divisor = phi(i)
                  else
                     divisor = 1.
                  end if
               end if
               pdifmax = max(abs((phi(i)-pho(i))/divisor),pdifmax) 
               if (spercent) then
                     divisor = 1.
               end if
               sdifmax = max(abs((s(i)-so(i))/1.),sdifmax)   
               if (spercent) then
                  if (to(i) .ne. 0.) then
                     divisor = to(i)
                  else if (t(i) .ne. 0.) then
                     divisor = t(i)
                  else
                     divisor = 1.
                  end if
               end if
               tdifmax = max(abs((t(i)-to(i))/divisor),tdifmax)
               if (spercent) then
                  if (pcio(i) .ne. 0.) then
                     divisor = pcio(i)
                  else if (pci(i) .ne. 0.) then
                     divisor = pci(i)
                  else
                     divisor = 1.
                  end if
               end if
               pcidifmax = max(abs((pci(i)-pcio(i))/divisor),pcidifmax)
               accmax = max(abs(deni(i)*sx1(i)),accmax) 
            end do
            if(pdifmax.gt.toldp) go to 100   
            if(sdifmax.gt.tolds) go to 100   
            if(tdifmax.gt.toldt) go to 100   
            if(pcidifmax.gt.toldc) go to 100  
            if(accmax.gt.tacc) go to 100   
         endif
         isteady=1
 100     continue
      else if(iflg.eq.2) then
         if(time_ss.gt.0.0) then
	    itt = 1
         else
	    itt = 0
         endif
         if(days.ge.time_ss.and.itt.eq.1) then
            isteady = 2
            btol_save = balance_tol
            balance_tol = tolerance
            tolde_save = tolde
            tolde = tolerance
            go to 20
         endif
         if(.not. sflux_flag) then
            if (isteady .eq. 1) then
               isteady = 2
               itt = 0
            end if
            go to 10
         endif

 20      continue
         if(ico2.lt.0) then
            if (spercent) then
               if (qtoti .ne. 0) then
                  divisor = abs(qtoti)
               else if (inflow .ne. 0.) then
                  divisor = abs(inflow)
               else
                  divisor = 1.
               end if
            else
                divisor = day*8.64e4
            end if
c gaz 082819 always go through flowrate adjustment of tmch            
            flow_rate = (abs(inflow) - abs(qtoti))/divisor
            if(abs(flow_rate) .le. balance_tol) then
               if(isteady.eq.-1)then 
                  isteady = 1
               else if(isteady.eq.1) then
                  isteady = 2
                  itt = 0
               endif
            endif   

               if (shtl.ne.0) then
c gaz 082619                   
c                  ratio = max(1.0d0,abs(flow_rate/balance_tol))
c                  head_tol = min(ratio*shtl0,shtl)
                   ratio = 1. 
                   head_tol = min(ratio*shtl0,shtl)
               endif
               if (stmch.ne.0) then
                  flow_max = abs(flow_rate)
c                  ratio = flow_max/balance_tol
                  balfac13 = 1.5*balance_tol
                  balfac12 = 2.*balance_tol
                  balfac1 = 4*balance_tol
                  balfac2 = 100.*balance_tol
                  balfac3 = 1000.*balance_tol
                  if(flow_max.gt.balfac3.or
     &                   .days.lt.tol_str) then
                   tmch1 = 100.*stmch0
          write(ierr,*) 'balfac3 flow_max ', flow_max,' tmch1 ',tmch1
                  else if(flow_max.gt.balfac2) then
                   tmch1_min = 50.*stmch0
                   tmch1_max = 100.*stmch0
                   tmch1 = ((tmch1_max-tmch1_min)/(balfac3 -balfac2))*
     &               flow_max + tmch1_min
          write(ierr,*) 'balfac2 flow_max ', flow_max,' tmch1 ',tmch1     
                  else if(flow_max.gt.balfac1) then 
                   tmch1_min = 10.*stmch0
                   tmch1_max = 50.*stmch0
                   tmch1 = ((tmch1_max-tmch1_min)/(balfac2 -balfac1))*
     &               flow_max + tmch1_min 
          write(ierr,*) 'balfac2 flow_max ', flow_max,' tmch1 ',tmch1      
                  else if(flow_max.gt.balfac12) then
                   tmch1_min = 5.*stmch0
                   tmch1_max = 10.*stmch0
                   tmch1 = ((tmch1_max-tmch1_min)/(balfac1 -balfac12))*
     &               flow_max + tmch1_min  
          write(ierr,*) 'balfac12 flow_max ', flow_max,' tmch1 ',tmch1      
                  else if(flow_max.gt.balfac13) then
                   tmch1_min = 1.5*stmch0
                   tmch1_max = 5.*stmch0
                   tmch1 = ((tmch1_max-tmch1_min)/(balfac12 -balfac13))*
     &               flow_max + tmch1_min 
            write(ierr,*) 'balfac13 flow_max ', flow_max,' tmch1 ',tmch1       
                  else if(flow_max.gt.balance_tol) then
                   tmch1_min = stmch0
                   tmch1_max = 1.5*stmch0
                   tmch1=((tmch1_max-tmch1_min)/(balfac13-balance_tol))*
     &               flow_max + tmch1_min   
            write(ierr,*) 'bal_tol flow_max ', flow_max,' tmch1 ',tmch1       
                  else                   
                   tmch1 = stmch0                    
                  endif 
c gaz 080519 tolerance should decrease monotonically                     
                  tmch = min(tmch1,tmch_old)   
                  write(ierr,*) ' tmch ', tmch,' tmch_old ',tmch_old   
                  tmch_old = tmch  
               endif			 		 
c            endif
         else if(ico2.ge.0) then
            if (spercent) then
               if (qtoti .ne. 0) then
                  divisor = abs(qtoti)
               else if (inflow .ne. 0.) then
                  divisor = abs(inflow)
               else
                  divisor = 1.
               end if
               if (qtotei .ne. 0) then
                  divisore = abs(qtotei)
               else if (inflowe .ne. 0.) then
                  divisore = abs(inflowe)
               else
                  divisore = 1.
               end if
            else
               divisor = day*8.64e4
               divisore = day*8.64e4
            end if
            flow_rate = (abs(inflow) - abs(qtoti))/divisor
            enth_rate = (abs(inflowe) - abs(qtotei))/divisore
            if(abs(flow_rate) .le. balance_tol .and. 
     &           abs(enth_rate) .le. tolde) then
               if(isteady.eq.-1)then 
                  isteady = 1
               else if(isteady.eq.1) then
                  isteady = 2
                  itt = 0
               endif
            else
               if (stmch.ne.0) then
                  ratio = max(1.0d0,abs(flow_rate/balance_tol),
     &                 abs(enth_rate/tolde))
                  tmch1 = min(ratio*stmch0,stmch)
                  ratio = tims/days
                  tmch2 = min(ratio*stmch0,stmch)
                  tmch = min(tmch1,tmch2)
               end if

            endif
         end if

 10      continue
         if(isteady .eq. 2) then
            if (l .lt. sminstep) then
               isteady = 1
c               return
            end if
c     If we have reached steady state or reached our maximum
c     time limit, make sure we output information
c gaz 033020             
            if(itt .eq. 0.and.isteady.ge.1) then
               if (iout .ne. 0) write(iout,50) days, trim(info_string)
               if(iptty.ne.0) write(iptty,50) days, trim(info_string) 
               if(isteady.eq.1) write(iout,57)
               if(isteady.eq.1) write(iptty,57)
            else
               if (iout .ne. 0) then
                  write (iout, 52) snstep, time_ss, days
                  write (iout, 53) trim(info_string)

               end if
               if(iptty.ne.0) then
                  write (iptty, 52) snstep, time_ss, days
                  write (iptty, 53) trim(info_string)
                  if(isteady.eq.1) write(iptty,57)
               end if
               balance_tol = btol_save
               tolde = tolde_save
            endif
         end if
           
         call write_steady_state

         if (isteady .eq. 2) then
            ntty = 2
            
            if (isty .eq. 0) then
               ifinsh = 1
               isteady = 0
               istea_pest = 1
            else 
               ifinsh = 2
               isteady = 0
               tims = tims_trans               
               daycs=daycs_save
               daycf=daycf_save
               dayhf=dayhf_save
               dayhs=dayhs_save 
               day = day_save_ss
               daymax = daymax_save
               nstep = nstep_save               
               aiaa = aiaa_save
               tmch = tmch_save
            end if
         end if
      else if(iflg.eq.3) then
         if(isteady .ne. 2) then
            if (iout .ne. 0) write(iout,53) trim(info_string)
            if (iptty .ne. 0) write(iptty,53) trim(info_string)

            call write_steady_state

            if (iout .ne. 0) write(iout,54) snstep, time_ss, days
            if (iptty.ne.0) write(iptty,54) snstep, time_ss, days
         endif
      endif

      return

 50   format(/,'>>>>> Steady State achieved at ', g16.9, 
     &     ' days (', a, ') <<<<<')
 52   format(/,'>>>>> Steady state time exceeded <<<<<',/,1p,
     &     'Max timesteps: ', i6, ' Specified Time: ',
     &     g15.6,' Simulated Time: ',g15.6)
 53   format(/,'>>>>> Steady State NOT achieved with (', a, ') <<<<<')
 54   format(/,'>>>>> Time control information <<<<<',/,1p,
     &     'Max timesteps: ', i6, ' Specified Time: ',
     &     g15.6,' Simulated Time: ',g15.6, /)
 56   format('>>>>> Total Mass ',1x,1p,g15.6
     &     'Max timesteps: ', i6, ' Specified Time: ',
     &     g15.6,' Simulated Time: ',g15.6, /)   
 57   format(/,'>>>>> Min steady state timesteps not met <<<<<')
      contains 

      subroutine write_steady_state

      if (toldh .ne. tolerance)
     &     hdifmax = pdifmax /( 997. * 9.8d-6)               
      if (iout .ne. 0) then
         if (isteady  .ne. 2 .and. iflg .ne. 3) 
     &        write (iout, 41) trim(info_string)
         if (toldp .ne. tolerance)
     &        write(iout, 40) 'pressure', pdifmax, toldp
         if (toldp .ne. tolerance. and. i_pdiff. ne. 0) 
     &        write(iout, 43) i_pdiff, pmax_i, pmax_io,izonef(i_pdiff)
         if (toldh .ne. tolerance)
     &        write(iout, 40) 'head', hdifmax, toldh
         if (tolds .ne. tolerance)
     &        write(iout, 40) 'saturation', sdifmax, tolds
         if (tolds .ne. tolerance. and. i_sdiff. ne. 0) 
     &        write(iout, 44) i_sdiff, smax_i, smax_io, izonef(i_sdiff)         
         if (toldt .ne. tolerance)
     &        write(iout, 40) 'temperature', tdifmax, toldt
         if (tacc .ne. tolerance)
     &        write(iout, 40) 'nodal mass', accmax, tacc
         if (tacc .ne. tolerance. and. i_accdif. ne. 0) 
     &    write(iout,45)i_accdif,accmax/(sx1(i_accdif)*denh(i_accdif)),
     &            amass_ch,izonef(i_accdif)         
         if (toldc .ne. tolerance)         
     &        write(iout, 40) 'air pressure', pcidifmax, toldc         
         if (balance_tol .ne. tolerance) write(iout, 40) 
     &        'flow balance (IN-OUT)', flow_rate, balance_tol
         if (tolde .ne. tolerance) write(iout, 40) 
     &        'enthalpy balance (IN-OUT)', enth_rate, tolde
            write(iout,42) tmch
c gaz 033020            
         if (isty .ne. 0.and.isteady.eq.2) write(iout,51)
      end if
      if(iptty.ne.0) then
         if (isteady .ne. 2 .and. iflg .ne. 3) 
     &        write (iptty, 41) trim(info_string)
         if (toldp .ne. tolerance)
     &        write(iptty, 40) 'pressure', pdifmax, toldp
         if (toldp .ne. tolerance. and. i_pdiff. ne. 0) 
     &        write(iptty, 43) i_pdiff, pmax_i, pmax_io, izonef(i_pdiff)
         if (toldh .ne. tolerance)
     &        write(iptty, 40) 'head', hdifmax, toldh
         if (tolds .ne. tolerance)
     &        write(iptty, 40) 'saturation', sdifmax, tolds
         if (tolds .ne. tolerance. and. i_sdiff. ne. 0) 
     &        write(iptty, 44) i_sdiff, smax_i, smax_io, izonef(i_sdiff) 
         if (toldt .ne. tolerance)         
     &        write(iptty, 40) 'temperature', tdifmax, toldt
         if (tacc .ne. tolerance)
     &        write(iptty, 40) 'nodal mass', accmax, tacc
         if (tacc .ne. tolerance. and. i_accdif. ne. 0) 
     &    write(iptty,45)i_accdif,accmax/(sx1(i_accdif)*denh(i_accdif)),
     &            amass_ch,izonef(i_accdif)         
         if (toldc .ne. tolerance)
     &        write(iptty, 40) 'air pressure', pcidifmax, toldc
         if (balance_tol .ne. tolerance) write(iptty, 40) 
     &        'flow balance (IN-OUT)', flow_rate, balance_tol
         if (tolde .ne. tolerance) write(iptty, 40) 
     &        'enthalpy balance (IN-OUT)', enth_rate, tolde
         write(iptty,42) tmch
         if (isty .ne. 0.and.isteady.eq.2) write(iptty,51)
      endif 
      
 40   format ('>>>>> Max change in ', a, 1pg12.4, ' tolerance',
     &     1pg12.4, ' <<<<<')
 43   format ('>>>>> Gridblock Max Pres change ',i8,
     &     ' Pres(new) ',1p,g12.4,' Pres(old) ',g12.4,
     &      ' Zone ',i6,' <<<<<')
 44   format ('>>>>> Gridblock Max Sat change ',i8,
     &     ' Sat(new) ',1p,g12.4,' Sat(old) ',g12.4,
     &     ' Zone ',i6,' <<<<<')
  45  format ('>>>>> Gridblock Max Mass change ',i8,' Frac change ',1p, 
     &       g12.4,' Total Mass change ',g14.6,' Zone ',i6,' <<<<<')     
 41   format (/, '>>>>> Steady state check (', a, ') <<<<<')
 42   format ('>>>>> Equation tolerance (',1p,g12.4, ') <<<<<')    
 51   format(/,'>>>>> Starting transient simulation',
     &           ' after steady solution <<<<<',/)
     
      end subroutine write_steady_state

      end subroutine steady
