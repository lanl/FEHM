      subroutine timcrl
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
CD1 Procedure to adjust timestep.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 25-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/timcrl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:18   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Fri Jun 28 09:04:36 1996   zvd
CD2 Added timestep error termination write to ierr
CD2 
CD2    Rev 1.7   Tue Jun 04 09:22:50 1996   gaz
CD2 modified time step cut in TS failure so it would always decrease
CD2 
CD2    Rev 1.6   Fri May 10 08:48:54 1996   gaz
CD2 changed 'ditnd =  1.0d+20' to 'ditnd =  1.0d+80'
CD2 
CD2    Rev 1.5   Fri Feb 02 12:53:40 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   08/16/95 16:28:24   robinson
CD2 Changed name of variable to set print out interval
CD2 
CD2    Rev 1.3   06/02/95 10:30:46   llt
CD2 read in absolute value of MAXIT in ctrl macro (gaz)
CD2 
CD2    Rev 1.2   01/28/95 13:56:20   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 15:56:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:56   pvcs
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
CD3   iatty                    O    File used for check information
CD4   iout                     O    File used for general code output
CD3   iptty                    O    Output to the tty
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
CD4   aiaa            REAL*8   faar   Time step multiplication factor
CD4   aiar            REAL*8   faar   Time step reduction factor
CD4   aw              REAL*8   faar   Time step weighting parameter for heat and
CD4                                     mass solution (ctrl)
CD4   awt             REAL*8   faar   Value of implicitness factor
CD4   ay              REAL*8   faar   Time step weighting parameter for tracer
CD4   day             REAL*8   faar   Current time step size in days
CD4   daycf           REAL*8   faar   Time at which tracer solution stops
CD4   daycs           REAL*8   faar   Time at which tracer solution starts
CD4   dayhf           REAL*8   faar   Time at which flow solution stops
CD4   dayhs           REAL*8   faar   Time at which flow solution starts
CD4   daymax          REAL*8   faar   Maximum time step allowed
CD4   daymin          REAL*8   faar   Minimum time step allowed
CD4   daynew          REAL*8   faar   Parameter used in time step control
CD4   days            REAL*8   faar   Current simulation time
CD4   daysi           REAL*8   faar   Simulation time at last time step
CD4   daysp           REAL*8   faar   Time at next time step
CD4   dit             REAL*8   fdd1   Array containing time step changes
CD4   ditnd           REAL*8   faar   Next time step change time
CD4   iac             INT      faai   Counter for print-out interval
CD4   iad             INT      faai   Current iteration number in flow solution
CD4   iamm            INT      faai   Maximum iterations allowed for time step
CD4                                     increase (heat and mass solution)
CD4   iatty           INT      faai   Unit number for check file
CD4   iccen           INT      faai   ?
CD4   icgts           INT      faai   Parameter controlling the time of solution
CD4                                     parameter changes
CD4   icontr          INT      faai   Parameter used in contour plot management
CD4   ics             INT      faai   Parameter indicating status of tracer
CD4                                     solution
CD4   icf             INT      faai   Parameter indicating status of tracer
CD4                                     solution
CD4   ifinsh          INT      faai   Indicates if the finishing criteria for
CD4                                     the simulation is achieved
CD4   ihf             INT      faai   Parameter indicating status of flow
CD4                                     solution
CD4   ihs             INT      faai   Parameter indicating status of flow
CD4                                     solution
CD4   iprtout         INT      faai   Print-out interval, number of time steps
CD4   ilt             INT      faai   Parameter used in time step control
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for tty
CD4   l               INT      faai   Current time step number
CD4   maxit           INT      faai   Maximum number of iterations allowed
CD4                                     before time step is halved
CD4   nicg            INT      faai   Parameter used in time step control
CD4   nsave           INT      faai   Indicates if a restart file will be
CD4                                     created in the current problem
CD4   tims            REAL*8   faar   Ending simulation time
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   co2ctr                   Control co2 simulation
CD4   contr                    Write contour plot file output
CD4   daycrl          REAL*8   Compute new time step size using adjusted
CD4                              timestep multipier if necessary
CD4   disk                     Write information to save file  
CD4   dual                     Control dual porosity simulation
CD4   plot                     Write history plot file output
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
CD5   dafac           REAL*8   ?
CD5   pravg           REAL*8   ?   
CD5   tfacd           REAL*8   ?
CD5   tmavg           REAL*8   ?
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
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.5.1 Implement time-step mechanism
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************

      use comdi
      use compart
      use comdti
      use comsplitts
      use comai
      use combi, only : ka
      use avsio, only : dit_flag
      implicit none

      real*8 dafac, daycrl, pravg, tfacd, tmavg
c gaz 032121 added integer i for do loop
      integer i 
      
      tfacd = max(5.0d00,aiaa*aiar)
      if (l .eq. 1 .and. iad .lt. abs(maxit))  then
         daynew =  daynew / aiaa
         iac    =  iprtout
         if (abs(daycs) .lt. zero_t .and. iccen .ne. 0)  then
            ics    =  1
            if (iout .ne. 0) write(iout  ,6000)
            if (iatty .gt. 0)  write(iatty ,6000)
 6000       format(/,1x,'tracer started at days = 0.0')
         end if
      end if
      
      if (daynew .lt. daymin)  day    =  daymin
      daysi  =  days
      if (iad    .gt. abs(maxit) )  then
         daynew =  daynew/tfacd
         dtot_next = 0.0
      endif
      if (iad    .le. iamm  )  then
         daynew =  daycrl(0)
      endif
      if (daynew .gt. daymax)  daynew =  daymax
      day    =  daynew
      
      if (day .lt. daymin)  then
         days   =  daysi
         if (iout .ne. 0) write(iout  ,6010)  l , day , days
         write(ierr  ,6010)  l , day , days
         if (iptty .gt. 0)  write(iptty ,6010)  l , day , days
 6010    format(/,1x,'timestep less than daymin',i5,2g20.6)
c**** reset variables to old ts values ****
         nsave  =  1
         call  dual    (1)
c gaz change nicg to nicg-1 below
c zvd 13-Jun-07 checked to make sure nicg > 1
         if (allocated(itc) .and. nicg .gt. 1) then
            if (itc(nicg-1).gt.0) then
               if (isave .ne. 0) call diskwrite
            endif
         endif
c**** set co2 pressure ****
         call  co2ctr  (2)
         days   = -1.0
         if (hist_flag) then
            call plot_new (1, 0.0d00, 0.0d00, 0.0d00, 0.0d00)
         else
            call plot (1, tmavg, pravg)
         end if
         call  contr   (-1)
         call  contr   ( 1)
         stop
      end if
      
      icontr =  0
      ifinsh =  0
c**** adjust time step ****
      daysi  =  days
      ilt    =  0
      daysp  =  days + day

c     only use dits if transient)
      if(isteady.eq.0) then
         if (daysp .ge. ditnd)  then
            day    =  ditnd - days
            ilt    =  1
            daysp  =  ditnd
c If dit_flag is true, contour data will not be written for each dit
c New cont keyword "nodit"
            if (.not. dit_flag) icontr =  1
            if (icgts .ne. 0)  then
               iprtout    =  iabs(itc(nicg))
               iac    =  iprtout
               dafac  =  dit(nicg+icgts)
               if (dafac  .lt. 0.0)  aiaa   = -dafac
               if (aiaa.lt.1.0) aiaa = 1.0
               if (dafac  .ge. 0.0)  daynew = dafac
               if (daynew .gt. daymax)  daynew = daymax
               aw = dit(icgts + icgts + nicg)
c     Add change to include changing max time step
               if(dit(3*icgts + nicg) .gt. 0.) then
                  daymax = dit(3*icgts + nicg)
               end if
               if (aw .lt. 1.0)  aw =  1.0
               if (aw .gt. 1.0)  aw =  1.5
               if (l .eq. 1)  then
                  awt    =  aw
                  aw     =  1.0
               end if
               ay     =  1.0 - aw
               nicg   =  nicg+1
               ditnd =  dit(nicg)
               if (nicg .gt. icgts)  ditnd =  1.0d+80
            end if
         end if
      endif

      if ((iccen .ne. 0).or.ptrak)  then
         if (daysp .ge. dayhf .and. ihf .eq. 0)  then
            if (iout .ne. 0) write(iout  ,6020)  l
            if (iatty .gt. 0)  write(iatty ,6020)  l
 6020       format(/,1x,'note>>h.t. solution stopped after time step',
     *           i7)
            ihf    =  1
            iac    =  iprtout
            if (l .ne. 1)  then
               day    =  dayhf - days
               daysp  =  dayhf
            end if
         end if

         if (daysp .ge. dayhs .and. ihs .eq. 0)  then
            if (iout .ne. 0) write(iout  ,6021)  l
            if (iatty .gt. 0)  write(iatty ,6021)  l
 6021       format(/,1x,'note>>h.t. solution started on time step',i7)
            ihs    =  1
            iac    =  iprtout
            if (l .ne. 1)  then
               day    =  dayhs - days
               daysp  =  dayhs
            end if
         end if

         if (daysp .ge. daycs .and. ics .eq. 0)  then
            if (iout .ne. 0) write(iout  ,6022)  l
            if (iatty .gt. 0)  write(iatty ,6022)  l
 6022       format(/,1x,'note>>tracer solution started on time step',i7)
            ics    =  1
            iac    =  iprtout
            if (l .ne. 1)  then
               if (abs(daycs) .gt. zero_t)  then
                  day    =  daycs - days
                  daysp  =  daycs
               end if
            end if
         end if

         if (daysp .ge. daycf .and. icf .eq. 0.)  then
            if (iout .ne. 0) write(iout  ,6023)  l
            if (iatty .gt. 0)  write(iatty ,6023)  l
 6023       format(/,1x,'note>>tracer solution stopped after time step',
     *           i7)
            icf    =  1
            iac    =  iprtout
            if (l .ne. 1)  then
               day    =  daycf - days
               daysp  =  daycf
            end if
 
         end if
      end if

c     
c     adjust for time varying boundary conditions
c     
      days=daysp
c gaz 112920 add coding to reset water source to original
c subroutine thrmwc (AWH) divides water source into water and air 
c gaz 032121 (only use sk0 for AWH (ico2 gt 0) and specified source)
      if(iboun.ne.0) then
       call flow_boundary_conditions(3)
      else if(ico2.gt.0) then
        do i = 1, n0  
         if(ka(i).eq.1) sk(i) = sk0(i) 
        enddo
      endif
      daysp=days
      days=daysp-day
      if (nicg .gt.1) then
         if (daysp .lt. dit(nicg - 1)) then
            nicg = nicg - 1
            ditnd =  dit(nicg)
         end if
      end if
         

      if (daysp .ge. tims)  then
c     Took out check for rip changes - need to determine why it was there
c     if (l .ne. 1)  then
         day    =  tims - days
         daysp  =  tims
         iac    =  iprtout
c     end if
         ifinsh =  1
      end if
      days   =  daysp
      if (day .lt. daymin)  then
         days   =  daysi + daymin
         day    =  daymin
      end if

      return
      end
