      subroutine  concen  ( iz, lstep )
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To control the status of the tracer simulation, calling routines
CD1 that either perform calculations or read or write information.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/concen.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Mon Mar 04 16:01:30 1996   hend
CD2 Removed uneccessary calculations from coneq1 and added trac input option
CD2 
CD2    Rev 1.8   Thu Feb 15 09:25:14 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.7   Mon Jan 29 13:53:48 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   Mon Jan 08 10:44:02 1996   robinson
CD2 Algorithm no longer recomputes things once heat and mass is turned off
CD2 
CD2    Rev 1.5   08/16/95 16:26:24   robinson
CD2 Changed call to wrtcon for writing solute info to .out file
CD2 
CD2    Rev 1.4   05/08/95 12:48:14   robinson
CD2 Corrected write to .trc for ptrack simulation
CD2 
CD2    Rev 1.3   03/15/95 17:04:24   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.2   01/28/95 14:03:34   llt
CD2 modified for new particle tracking module
CD2 
CD2    Rev 1.1   03/18/94 16:15:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:24   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name       Type        Description
CD3 
CD3 iz          I          Control parameter determining the
CD3                        subroutine to call and processing to perform
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 iccen, days, daycs, daycf, ics, icf
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
CD4 rdcon   N/A      Read in solute input
CD4 csolve  N/A      Perform solute computations until next stopping
CD4                  point is found
CD4 wrtcon  N/A      Write solute information to output file
CD4 diskc   N/A      Write solute information to restart file
CD4 contrc  N/A      Write solute information to contour files
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following functions are carried out in this routine:
CD6 
CD6   The code checks to see if there is a tracer simulation, and if
CD6   there is:
CD6   
CD6     The code uses an IF-ELSEIF construct to decide among the
CD6     following actions:
CD6     
CD6       Reading of solute input information.
CD6       
CD6       Compute concentrations at current time (only if tracer
CD6       solution is enabled at this time).
CD6       
CD6       Writing of solute information to output file.
CD6       
CD6       Reading or writing of solute information to (from) restart
CD6       file.
CD6       
CD6       Writing of solute information to contour file.
CD6       
CD6   The code then returns to the calling routine.
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.3.4 Solute-transport equations
CD9 2.7   Provide Restart Capability
CD9 2.7.3 Resume the calculation
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN concen
CPS 
CPS IF a solute simulation is being performed
CPS 
CPS   IF we a reading input information at this point in the run
CPS     rdcon - read solute input information
CPS   ELSEIF we a potentially computing solute transport at the...
CPS   ... current time
CPS   
CPS     IF the tracer solution is enabled
CPS       csolve - compute solute concentrations at the current time
CPS     ENDIF
CPS     
CPS   ELSEIF we are writing solute information
CPS     wrtcon - write solute information to output file
CPS   ELSEIF we are reading or writing to the restart file
CPS     diskc - read or write solute information to restart file
CPS   ELSEIF we are writing to the contour file
CPS     contrc - write solute information to the contour file
CPS   ENDIF
CPS   
CPS ENDIF
CPS   
CPS END concen
CPS 
C**********************************************************************
c version FEHM5.1J changes
c 24-june-91
c modify routine CONTRC to also for free or unformatted output
c 19-july-91
c modify routine thermc so iscnd ge 0 means liquid phase tracer
c 22-july-91
c modified max,min in rdcon so both arguments are double precision
c 2-august-91
c changed coding in diskc there were errors
c 6-august-91
c changed coding for new gravity term
c 7-august-91
c put in dispersion term=diffm+alpha*velocity
c changes in rdcon thermc and coneq1
c 19-august-91
c changed equiv to reflect smaller size of common /fcc/
c 4-oct-91
c checked for min value of diffmf(=1.e-20)
c 8-oct-91
c major input changes for adsoprtion changes to be spacially varying
c changed subscripts for tcx,tcy,tcz
c 11-oct-91
c added coding in cnswer so dpdp solution is available
c removed call in cnswer to thrmc for dual solution
c 15-nov-91
c changed gencon,coneq1 to accommodate symmetry programming
c changed usage of awc,ayc to be like heat-mass
c set dencj=denci
c 16-nov-91
c if(icns(nsp)... out of place in coneq1,corrected it
c added iau,ial def to coneq1
c added call to dualp in gencoN
c 25-nov-91
c kb>kz corrections in coneq1
c ps(i)>ps(iz)
c a(jmia)>a(jmia+nmat(1))
c 3-dec-91
c changed to generic intrisic functions
c 19-mar-92
c changed internal reads to conform to strict f77
c defined swi in coneq1
c IMPORTANT-changed def of tracer by +1(ie conservative tracer=iadsor of 1)
c 26-mar-92
c delete backspace around line 164 in rdcon
c 26-may-92
c made t9=fid always
c 3-june-92
c took out backspace isave in diskc
c 7-aug-92
c took out subroutines thermc and rdcon,they are separate in trad_br.f
c 11-aug-92
c read model data from single group
c in routine rdcon
c 1-sept-92
c set read not read to an0 in diskc
c 21-sept-92
c changed source designation for wrtout
c called cbon in CSOLVE to remove boundary nodes from balance calcs
c 27-jan-93
c added h,s,c,to solve input line
c 3-may-93
c changed sk(ja-npn) in loop 122 to rc(ja)
c substituted srmim for sk(mdd) in routine wrtcon
c made irdof for tracers enabled for irdof=9(ie never)
c 11-may-93
c set iter=3*north in gencon
c 18-may-93
c added printout headings for different matrix levels
c 7-may-93 llt
c changed equivalences to pointers
c****--------------------------------------------------------------****c
c**** trad.f                                                       ****c
c****--------------------------------------------------------------****c
c
      use combi
      use comci
      use comdi
      use comgi
      use compart
      use comdti
      use comai
	use comrxni, only : cden_flag
      implicit none
      
      integer iz, hmon, lstep, i , i2
      real(8) begin_time, end_time
c hmon is 1 if heat and mass solution is active

! zvd 22-Jan-04 diskc and diskp are no longer called from this routine
!               check for use of both trac and ptrk moved to scanin     
      if (ihf.ne.ihs .or. .not. compute_flow) then
         hmon=0
      else
         hmon=1
      endif
      if (ptrak) then
         if (iz.eq.1) then
            if ((days.le.daycs).or.(days.gt.daycf)) return
            if(lstep.lt.0) then
               hmon = 1
               lstep = -lstep
            end if
            begin_time=8.64d4*days-dtot
            end_time=8.64d4*days
            call part_track(begin_time,end_time,hmon,lstep)

         else if (iz.eq.2) then
            if ((days.le.daycs).or.(days.gt.daycf)) return
            call wrtptrk
            npn = 0
            call plotc1(1,0)
         endif
      else if (iccen.ne.0) then
         if ( iz .le. 0 )  then
            call  rdcon   ( iz )
         else if ( iz .eq. 1 ) then
            if (iz.eq.1) then
               if ( days .le. daycs .or. days .gt. daycf  ) return
               call  csolve(hmon)
            endif
         else if ( iz .eq. 2 )  then
c gaz 033020 check for steady solution if transient
c -----------------------------------------------------              
c gaz 041220 need last printout          
            if ( ics .ne. 0 .or. icf.ne. 0) then   
              if(daysp .le. daycf) then
               call  wrtcon(1)
              endif
            end if
         else if ( iz .eq. 5 )  then
            call  contrc
         else if ( iz .eq. 6 )  then
          if(cden.and.cden_flag.eq.0) then
           do i = 1,neq
            i2 = i + (ispcden - 1)*n0
            anl(i2) = anlo(i2)
           enddo
          endif
         endif
      end if

      return
      end
