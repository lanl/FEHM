      function daycrl(iflg) 
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
CD1 Compute new time step size using adjusted timestep multipier if necessary. 
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/daycrl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Feb 15 09:59:30 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.2   Mon Jan 29 14:48:16 1996   hend
CD2 Updated Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:53:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:56   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   day             REAL*8   I    Current time step size in days
CD3   aiaa            REAL*8   I    Time step multiplication factor
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 None
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
CD5   aiaab           REAL*8   Computed time step multiplication factor
CD5   daylog          REAL*8   Log base 10 of current timestep size
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C*******************************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
CD6
C*******************************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C*******************************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This subroutine controls the time step size only.
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
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN  daycrl
CPS 
CPS   IF the timestep multiplier is greater than zero
CPS      compute new time step size
CPS   ELSE
CPS      get log of time step size and use it to compute new multiplier . . .
CPS      . . . then compute new timestep size
CPS   ENDIF
CPS      
CPS END  daycrl
CPS
C***********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comgi
      use compart
      use davidi 
      use comsplitts
      implicit none


      integer i,iflg,icd,iz,neqp1
      integer ifac,izfac,jj,ncyc
      real*8  aiaab, daycrl, daylog
      real*8  stest,s_up,stol
      real*8  sfac,sfac1,sfac2,sfac3,sfac_old
      real*8  dtot_sfac,stest_ch
      logical dayset
      parameter(stol=0.5d00,s_up=0.999d00,ncyc=1)

      dayset = .false.
      if (icgts .ne. 0) then
         if (days .eq. dit(icgts) .and. dit(icgts+1) .gt. 0) then
! use value set from dit2
            daycrl = daynew
            dayset = .true.
         end if
      end if
      if (dayset) then
! Set daycrl above from dit
      else if (impf .ne. 0) then
         if(fimp.gt.1.0d00) then
           daycrl = 
     &      max(daymin,daynew,daynew*(aiaa/fimp),dtot_next/86400.)
         else
           daycrl =  max(daynew * aiaa, dtot_next/86400.)
         endif
 
      elseif (aiaa .ge. 1.0) then

c        daycrl =  daynew * aiaa
         daycrl =  max(daynew * aiaa, dtot_next/86400.)

      elseif (aiaa .eq. 0.0) then

         daylog =  log10 (daynew)

         aiaab  =  1.0 + 10.0**(-daylog**3 / 60.0 - daylog**2 / 12.0
     *        - daylog * 4.0 / 15.0 - 1.0) * 0.3333

         daycrl =  daynew * aiaab

      else

c don't change timestep
         write (iout, *) '*** daycrl not changed = ', daycrl, ' ***'

      end if
 
      daycrl = min (daycrl,daymax)
      dtot_next = daycrl*86400.

      if(isplitts.eq.+1) then
c check for timestep limitation from possible phase change
c get fluxes


      do jj=1,ncyc 
      dtot = daycrl*86400.
       call airctr(3,0)
         bp = 0.0d00
         do i=1,neq
           call flux_net(i)
           bp(i+nrhs(1)) = bp(i+nrhs(1)) - sx1(i)*deni(i)
     &                     - sk(i)
         enddo

      sfac1=1.0d00
      sfac2=1.0d00
      sfac3=1.0d00
      sfac=1.0d00
      sfac_old = 1.0d00
c
c form constants for i>neq
c
      neqp1=neq+1

      do i=1,n

      if(i.gt.neq.and.idualp.eq.0.and.idpdp.eq.0) then
         icd=neq
         iz=i - neq
      else
         icd=0
         iz=i
      endif

       stest = (dtot*bp(iz+nrhs(1))/sx1(i)+ denh(i))
     &       /(rolf(i)*ps(i))
c
c compute saturation solution
c
c      s(i) = min(1.0d00,max(stest,0.0d00))
       ieos(i) = 2
       if(stest.le.0.0.and.so(i).gt.0.0) then
        sfac1 = (so(i)-stest)/so(i)
        ieos(i) = 3
       endif
       if(stest.ge.1.0.and.so(i).lt.s_up) then
        sfac2 = (so(i)-1.0)/(so(i)-stest)
        ieos(i) = 1
       endif
       if(abs(stest-so(i)).gt.stol) then
        sfac3 = stol/abs(stest-so(i))
       endif
         sfac = min(sfac1,sfac2,sfac3)                       
         if(sfac.lt.sfac_old) then
           sfac_old = sfac
           ifac = i
           izfac = i
         endif
      enddo
     
      if(sfac.lt.1.0d00) then
       dtot_sfac = sfac*dtot

       if(explicit_factor.gt.0.0d00) then
         isplitts = -2
       else
         isplitts = +2
       endif

      else
       go to 999
      endif
      enddo
999   daycrl = sfac*daycrl         

      else if(isplitts.eq.-2.and.mlz.eq.0) then
       isplitts = +2
c     daycrl = explicit_factor
      else if(isplitts.eq.+2.and.mlz.eq.0) then
       isplitts = +1
      endif

      end
