      subroutine bnswer_part
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To call routines to assemble finite element equations and solve for
!D1 the Newton-Raphson equations. 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/bnswer_part.f_a  $
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.5.2 Solve nonlinear equation set at each time step  
!D3 2.3.1 Heat-conduction equations
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
!D3
!**********************************************************************
!D4
!D4 SPECIAL COMMENTS
!D4
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
!D4    FEHM Application Version 2.20
!D4
!**********************************************************************

      use combi
      use comdi
      use comgi
      use comei
      use comdti
      use comai
      use comsplitts
      use com_part
      use comco2
      implicit none

      integer  iad_ck
      integer zone
      real*8 fdum_tot,strd_part
      real*4 caz(2)   
      real*8 tyming,tstart
      parameter (strd_part=1.0)

      itert=0
      minkt=0
      strd=1.
c gaz 110919 moved iad_min to comai (global variable now)     
      if(g1.lt.0.0.and.jswitch.eq.0) then
	 iad_min = 0.0
      else if(maxit.ge.0) then
         iad_min=1
      else
         iad_min=2
      endif
      if(iflux_ts.ne.0.and.ico2.lt.0) then
         call cascade_sat(0)
      endif

      iad_part = 0
      itert_part = 0
      timing_part = 0.0
c     
c     check iterations against maximum
c     
 1000 continue

      fdum_tot = 0.0
      fdum_part = 0.0
      f0_part = 0.0
      iad_ck = 0
      
c     loop on every zone, must pass the zone number into everything
c     so that the correct arrays and counts are used
c     must locate each routine that uses nelm because that is different
c     if the routine only uses neq along with other arrays it is OK
c     don't forget to look into routines that are called within a routine
c     
      do zone = 1, nzone_para
         
c     get relative CPU times for each zone
         tstart = tyming(caz)

         if(iad_part(zone).gt.abs(maxit)) then
            mlz=-1
            goto 2000
         endif

         if(iad_part(zone).gt.abs(maxit)/2) then
            strd=strd_part
         endif
c     
c     update variable states and call EOS packages
c     

c     uses neq but not nelm so don't change
c     now passes zone to airctr_part
         call varchk_part(0,0,zone)

c     call appropriate sub to generate equations
         if(idof.le.1) then
            call gensl3
         else
            if(ico2.gt.0.and.ice.eq.0) then
               if(idpdp.eq.0) then
                  call gensl4
               else
                  call dpdp(3)
               endif
            endif

            if(ico2.lt.0.and.ice.eq.0) then
c     single porosity
               if(idpdp.eq.0) then
                  if(isplitts.ge.-1) then
                     call airctr_part(4,0,zone)

                  else if(isplitts.eq.-2) then
                     call airctr_part(-4,0,zone)
c     explicit iterations finished in airctr so return
                     go to 2000
                  endif
               else
c     dpdp enabled
                  call dpdp(3)
               endif
            endif

c     gaz 10-15-2001 (multi-componet system)
            if(ice.ne.0) then
c     single porosity
               if(idpdp.eq.0) then
                  call icectr(4,0)
               else
c     dpdp not enabled
c     
               endif
            endif

c     RJP 04/10/07 added following make sure works.
            if(icarb.eq.1) then
               call icectrco2(4,0)
            endif

            if(ico2.eq.0.and.ice.eq.0) then
               if(idpdp.eq.0) then
                  call gensl1
               else
                  call dpdp(3)
               endif
            endif

         endif

         if (fdum_part(zone) .ne. -1.0) then
            fdum_tot = fdum_tot + fdum_part(zone)
         endif

         iad_part(zone) = iad_part(zone) +1
         if(iad_part(zone).ge.iad_min) iad_ck=iad_min
         itotal_part(zone) = itotal_part(zone) +1
         timing_part(zone) = timing_part(zone) +
     &        tyming(caz)-tstart              

c     end loop on zones
      enddo

      if(fdum_tot.le.0.0.and.iad_ck.ge.iad_min) goto 2000
c     turn off initial parent BC if appropriate
c     
c     apply nr corrections (ie update variables)
c     
      do zone = 1, nzone_para
         iad=iad_part(zone)
         call varchk_part(1,0,zone)
      enddo
      call paractr(6)
      if(epe.le.0.) goto 2000
c     
c     check if varibles are out of bounds
c     

      call outbnd
      if(mlz.ne.0) goto 2000
      goto 1000
 2000 continue
c     
c     for bookkeeping  and output put bp_part back into bp 
c     
      iad=0
      do zone = 1, nzone_para
         iad=iad+iad_part(zone)
         itotal=itotal+iad_part(zone)
      enddo

      iad = iad/nzone_para
      do zone = 1, nzone_para
         itotals = itotals + itotals_part(zone)
      enddo
      call paractr(2)

      end
