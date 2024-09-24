      subroutine varchk(ifl,ndummy)
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
CD1 To determine variable set and make n-r corrections. 
CD1 gaz 042323 cleaned up comments
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY                                                 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                     However, previous non-YMP
CD2                                     versions of FEHM exist, and
CD2                                     the current version differs
CD2                                     from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/varchk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:10   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Thu Dec 18 15:23:50 1997   gaz
CD2 correction of ieosd in varchk 
CD2 
CD2    Rev 1.8   Fri Sep 26 15:14:06 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.7   Mon Mar 31 08:45:00 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.6   Fri Apr 26 16:39:44 1996   gaz
CD2 small changes in N-R step length
CD2 
CD2    Rev 1.5   Fri Feb 16 13:19:18 1996   zvd
CD2 Added/modified requirements.
CD2 
CD2    Rev 1.4   Wed Feb 07 10:05:06 1996   gaz
CD2 changed storage of NR steplength strd
CD2 
CD2    Rev 1.3   06/02/95 10:33:52   llt
CD2 gaz changes
CD2 
CD2    Rev 1.2   01/28/95 13:56:28   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 15:45:20   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:10   pvcs
CD2 original version in process of being certified
c 1/23/95 gaz for ngas and change from 2-phase to gas , remove eosml*pci
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier     Type    Use      Description
CD3 
CD3 ifl            int      I       Flag to determine the purpose of
CD3                                     this call to the routine
CD3 ndummy         int      I       Index directing code to correct
CD3                                     node number for dual porosity
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
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 ico2, neq, phi, t, s, ieos, ps, idof, pcp, dpcef, pnx, pci, ice,
CD4 ices, sii, iad, nr1, nrhs, bp, tmelt, strd
CD4 
CD4 Global Subprograms
CD4 
CD4 None
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5 
CD5 psatmn       real*8      minimum saturation
CD5 eosmg        real*8      eos parameter used in calculation
CD5 eosml        real*8      eos parameter used in calculation
CD5 eostol       real*8      eos parameter used in calculation
CD5 stepl        real*8      SOR underelaxation parameter used when
CD5                              phase changes occur
CD5 pcimin       real*8      Minimum gas pressure
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop parameter over all nodes
CD5 pl           real*8      Current pressure
CD5 ij           int         Current node number
CD5 x            real*8      Parameter used in calculation
CD5 tl           real*8      Current temperature
CD5 sl           real*8      Current saturation
CD5 ieosd        int         Current equation of state flag
CD5 ieosdc       int         Current equation of state flag
CD5 tboil        real*8      Boiling temperature
CD5 dtsatp       real*8      Parameter returned from call to psatl
CD5                             not used)
CD5 dpsats       real*8      Parameter returned from call to psatl
CD5                             not used)
CD5 dpsatt       real*8      Parameter returned from call to psatl
CD5                             not used)
CD5 pcl          real*8      Gas pressure
CD5 pboil        real*8      Vapor pressure at given temperature
CD5 pvapor       real*8      Vapor pressure
CD5 iced         int         Current flag for ice saturation
CD5 icedc        int         Current flag for ice saturation
CD5 siid         real*8      Current ice saturation
CD5 nr1          int         Index for solution arrays
CD5 nr2          int         Index for solution arrays
CD5 nr3          int         Index for solution arrays
CD5 i1           int         Index used in computing new variable values
CD5 i2           int         Index used in computing new variable values
CD5 i3           int         Index used in computing new variable values
CD5 
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 N/A
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
CD9 2.3.1 Heat-conduction equations
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9 2.4.3 Equation-of-state models
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
CPS BEGIN varchk
CPS 
CPS IF this call is to determine the eos status of each node
CPS 
CPS   IF this is a pure water simulation
CPS 
CPS     FOR each node
CPS       IF a phase change check is unnecessary
CPS         Set eos flag parameter
CPS       ELSEIF the fluid was previously liquid only
CPS         psatl - compute the boiling temperature
CPS         IF flashing of steam occurred
CPS           Set eos flag parameter, temperature, and saturation
CPS         ENDIF
CPS       ELSEIF the fluid was previously two-phase
CPS         IF the saturation is computed to be greater than unity
CPS           Set eos flag to pure liquid value
CPS           psatl - compute temperature
CPS           Set saturation to unity
CPS         ELSEIF the saturation is computed to be less than zero
CPS           Set eos flag to pure steam value
CPS           psatl - compute temperature
CPS           Set saturation to zero
CPS         ENDIF
CPS       ELSEIF the fluid was previously pure steam
CPS         psatl - compute boiling temperature
CPS         IF the calculation shows that some steam liquifies
CPS           Set eos flag to two-phase value
CPS           Set saturation and temperature
CPS         ENDIF
CPS       ENDIF
CPS       IF a change occurred in the state of the fluid
CPS         Set parameter
CPS       ELSE
CPS         Set parameter
CPS       ENDIF
CPS       Set eos flag in array
CPS     ENDFOR
CPS     
CPS   ELSEIF this is a noncondensible gas solution with heat
CPS   
CPS     FOR each node
CPS       IF a phase change check is unnecessary
CPS         Set eos flag parameter
CPS       ELSEIF the fluid was previously liquid only
CPS         psatl - compute the saturation vapor pressure
CPS         IF flashing of steam occurred
CPS           Set eos flag parameter, pressure, and saturation
CPS         ENDIF
CPS       ELSEIF the fluid was previously two-phase
CPS         IF the saturation is computed to be greater than unity
CPS           Set eos flag to pure liquid value
CPS           psatl - compute pressure
CPS           Set saturation to unity
CPS         ELSEIF the saturation is computed to be less than zero
CPS           Set eos flag to pure steam value
CPS           psatl - compute pressure
CPS           Set saturation to zero
CPS         ENDIF
CPS       ELSEIF the fluid was previously pure steam
CPS         psatl - compute vapor pressure
CPS         IF the calculation shows that some steam liquifies
CPS           Set eos flag to two-phase value
CPS           Set saturation, temperature, and pressure
CPS         ENDIF
CPS       ENDIF
CPS       IF a change occurred in the state of the fluid
CPS         Set parameter
CPS       ELSE
CPS         Set parameter
CPS       ENDIF
CPS       Set eos flag in array
CPS     ENDFOR
CPS     
CPS   ELSE this is an isothermal air-water simulation
CPS   
CPS     airctr - compute pressures and saturations
CPS     
CPS   ENDIF
CPS   
CPS   IF this is an ice solution
CPS     FOR each node
CPS     
CPS       Set flags
CPS       
CPS       IF no ice is currently present
CPS         IF ice has just formed
CPS           Set ice saturation and temperature
CPS         ENDIF
CPS       ELSEIF a partial melt solution exists
CPS         IF the calculation shows that all ice has melted
CPS           Set temperature and ice saturation
CPS         ELSEIF the calculation shows that the entire solution is ice
CPS           Set temperature and ice saturation
CPS         ENDIF
CPS       ELSEIF a solid ice solution exists
CPS         IF partial melt has just occurred
CPS           Set ice saturation and temperature
CPS         ENDIF
CPS       ENDIF
CPS       
CPS       Set flags
CPS       
CPS     ENDFOR
CPS   ENDIF
CPS   
CPS   vcon - explicitly update thermal conductivities on the first...
CPS   ... iteration
CPS   
CPS   IF this is a noncondensible gas simulation
CPS     thrmwc - compute thermodynamic parameters and derivatives
CPS   ELSEIF this is an isothermal air-water simulation
CPS     airctr - compute thermodynamic parameters and derivatives
CPS   ELSE this is a pure water simulation
CPS     thermw - compute thermodynamic parameters and derivatives
CPS   ENDIF
CPS   
CPS ELSE this call is to make the corrections to the variables
CPS 
CPS   IF this is a pure water simulation
CPS   
CPS     FOR each node
CPS       IF the fluid is pure water
CPS         Compute new temperature and pressure
CPS       ELSEIF the fluid is two-phase
CPS         Compute new pressure and saturation
CPS       ELSEIF the fluid is pure vapor
CPS         Compute new pressure and temperature
CPS       ELSEIF this is a temperature only solution
CPS         Compute new temperature
CPS     ENDFOR
CPS     
CPS   ELSEIF this is a noncondensible gas solution with heat
CPS   
CPS     FOR each node
CPS       IF the fluid is pure water
CPS         Compute new temperature, pressure, and gas pressure
CPS       ELSEIF the fluid is two-phase
CPS         Compute new pressure, saturation, and temperature
CPS       ELSEIF the fluid is pure vapor
CPS         Compute new pressure, temperature, and gas pressure
CPS       ELSEIF this is a temperature only solution
CPS         Compute new temperature
CPS     ENDFOR
CPS     
CPS   ELSEIF this is an isothermal air-water simulation
CPS     airctr - compute new pressures and saturations
CPS   ENDIF
CPS   
CPS   IF this is an ice solution
CPS   
CPS     FOR each node
CPS       IF this node is at partial melt conditions
CPS         IF this is a two-phase liquid solution
CPS           Correct temperature and compute ice saturation
CPS         ELSE
CPS           Correct temperature and compute ice saturation
CPS         ENDIF
CPS       ENDIF
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS 
CPS ENDIF
CPS 
CPS END varchk
CPS 
C**********************************************************************

      use comai
      use combi
      use comci
      use comco2
      use comdi
      use comdti
      use comei
      use comfi
      use comgi
      use comii
      use com_prop_data, only : pl_last, tl_last, pcl_last
      use comtable, only : tableFLAG
      use davidi
c gaz 011021 added com_exphase global variables    
      use com_exphase, only : iphase_chk
c gaz 080622       
      use commass_AWH, only : ivar_mass
  
      implicit none

      real*8 psatmn
      real*8 eosmg
      real*8 eosml
      real*8 eostol 
      real*8 stepl
c gaz 121918 added pci0
      real*8 pci0
      real*8 pcimin
      real*8 phase_mult
      real*8 phase_sat 
      real*8 satml 
c gaz 101419       
      real*8 xdiff_tol
c gaz 041823
      real*8 fdum_0
      save fdum_0
      parameter(psatmn=0.0001)
      parameter(eosmg=1.0001)
c      parameter(eosml=0.95)
c      parameter(eosml=0.99)
      parameter(eostol=0.0001)
c      parameter(stepl=0.95)
      parameter(pcimin=0.0)
c      parameter(phase_mult=1.00)
      parameter(phase_sat=1.0d-9)

c gaz 101419      
      parameter(xdiff_tol=1.0d-4)
c gaz 121718   0923219 h2o_crit   moved to comai
c      real*8 pcrit_h2o, tcrit_h2o
c      parameter(pcrit_h2o=22.00d0, tcrit_h2o=373.95)
c ich_max should be odd or even but don't know which
c      parameter(ich_max = 2)
      real*8 psatl
      integer i
      real*8 pl
      real*8 tl
      real*8 sl
      real*8 x
      real*8 dum
      integer ij
      integer ieosd
      integer ieosdc
      real*8 tboil
      real*8 pboil
      real*8 dtsatp
      real*8 dpsats
      real*8 dpsatt
      real*8 pcl
      real*8 pvapor
      integer iced
      integer icedc
      real*8 siid
      integer nr1
      integer nr2
      integer nr3
      integer i1
      integer i2
      integer i3
      integer ifl
      integer ndummy,k
c gaz 120821      
      integer eval_test_h2o_0
      real*8 stepl_hm, phase_mult_hm, satml_hm, eosml_hm
      real*8 stepl_hma, phase_mult_hma, satml_hma, eosml_hma
c gaz 112021 (water EOS table may have different pcrit_h2o and tcrit_h2o)       
      real*8 xnl,dxnlp,dxnlpc,dxnlt,air_mass1, air_mass_l, pcl_gas
      integer ipr_phase_detail
      parameter(ipr_phase_detail = 1)
      real*8  p0,t0,pc0,s0
      real*8  wgtchng1, wgtchng2, tNR_min
      character chariad*12
c gaz 111621 stepl_hm = 0.95 to stepl_hm = 0.85     
c gaz 111721 phase_mult_hm = 1.0001 to 1.0      
c      parameter(stepl_hm = 0.85, phase_mult_hm = 1.0001)
c gaz debug 112721   stepl_hm = 0.96 to 0.9    
       parameter(stepl_hm = 0.95, phase_mult_hm = 1.000)
c gaz debug 051516 (optimized for geothermal;satml = 1.0d-2 may be too large)      
      parameter(satml_hm = 1.0d-6, eosml_hm = 0.9999)
c gaz 101419  102119 ...mult_hma = 1.00001 is best
      parameter(stepl_hma = 0.95, phase_mult_hma = 1.00001)
      parameter(satml_hma = 1.0d-7, eosml_hma = 0.99)
      parameter(wgtchng1 = 0.9d0, wgtchng2 = 0.9d0,TNR_min=5.0d0)
c
      
      if(ico2.eq.0) then
c pure water and heat
        satml = satml_hm
        phase_mult = phase_mult_hm
        eosml = eosml_hm
      else
c ngas (air,methane,co2),water, heat (AWH)
        satml = satml_hma
        phase_mult = phase_mult_hma
        eosml = eosml_hma
      endif
      ich_m1 = 10
      ich_m2 = 10
      ich_max = 6   
c
c determine nr correction multiplier
c
      if(nr_stop.eq.0) then
         if(ico2.eq.0) then
          stepl = stepl_hm
         else
          stepl = stepl_hma
         endif
      else
         stepl = strd_iter
      endif
      if(var_stop.and.iad.ge.1.and.ifl.eq.1)then
         call nr_stop_ctr(1)
      endif
c gaz debug 090515  
c new code to save mass and energy in each block for ngas 
       if(ico2.gt.0) then    
         denei_ch = 0.0
         deni_ch = 0.0
         denpci_ch = 0.0
         ieos_bal = 0.0
       endif
      if(ifl.eq.0) then
c
c     evaluate eos status of each node. change status as prescribed
c     by phase change instructions
            
c set newton raphson step length to 1.0 at beginning of timestep
         if(iad.eq.0) then
            strd =1.0
            ieos_ch = 0
         else if(iad.ne.0.and.iter_intrvl.ne.0) then   
             if(mod(iad,iter_intrvl).eq.0) then
              strd =1.0
              ieos_ch = 0
             endif
         endif
cc
c
c new loop to check on degree of freedom 
c if idof=1, the problem is either heat conduction or saturated only
c 
c
         if(idof.ne.1) then
            if(ico2.eq.0) then
c
c     determine phase state for pure water
c
c gaz 111015 added super crtical phase
               if(icarb.eq.0) then
                  do i=1,neq
                     ij=i+ndummy
                     pl=phi(ij)
                     x=pl
                     tl=t(ij)
                     sl=s(ij)
                     ieosd=ieos(ij)
                     if(ieosd.eq.4) then
                      if(tableFLAG.eq.0.and.iwater_table.eq.0) then
c gaz 112020 stop if SC phase and no table write error msg and stop
                      if(iptty.ne.0) write(iptty,30)
                      if(iout.ne.0) write(iout,30) 
                      write(ierr,30)
30    format('>>>>>>>>>> super critical water phase but no EOS table',
     &       ' , stopping <<<<<<<<<<<<')           
                      endif   
                     endif
                     if(ieosd.eq.0) ieosd=1
                     ieosdc=ieosd
c
c     check if eos change is necessary
c
                     if(ps(ij).eq.0.0.or.idof.le.1) then
c                 ieosdc=-1 gaz 10-17-2001 do nothing
                     elseif(ieosd.eq.1) then                    
c
c     liquid only state
c
                     tboil = 2000.
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                      else if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3 
                        s(ij) = 0.0
                      else if(pl.lt.pcrit_h2o) then
                         tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                      endif
         
C*****
C***** AF 11/15/10
C****
                       if(tableFLAG.EQ.1) tboil = 1201. !  phs 4/27/99
C*****
c     change to 2-phase
                        if(tl.ge.tboil*phase_mult.
     &                       and.days.ge.time_ieos(ij)) then
                           ieosdc=2
c gaz debug 121817                           
c eosml = 0.9999 (042323) 
                           s(ij) = eosml
                           t(ij)=tboil
                           time_ieos(ij) = days + time_ch
c 
                        endif
c
                     elseif(ieosd.eq.2) then
c
c     2-phase conditions
c
c  gaz 102917 update temperature earlier in iteration    
                         
c  gaz 112121 check for SC conditions (also check for ieos = 1
                     tboil = 2000.
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                        s(ij) = 1.0
                      else if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3 
                        s(ij) = 0.0
                      else if(pl.ge.pcrit_h2o.and.tl.le.tcrit_h2o) then
                        ieosdc = 1 
                        s(ij) = 1.0                        
                      else if(pl.lt.pcrit_h2o) then
                         tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                         t(ij) = tboil
                      endif                         
                         
                         
                         
c gaz  112121    see tboil above                      
c                           t(ij)=psatl(pl,pcp(ij),dpcef(ij),
c     &                          dtsatp,dpsats,1,an(ij))  
c
                         
                        if(sl.gt.1..and.days.ge.time_ieos(ij)) then
c     change to liquid only conditions
                           ieosdc=1
c gaz 102917 calculated above                            
c                           t(ij)=psatl(pl,pcp(ij),dpcef(ij),
c     2                          dtsatp,dpsats,1,an(ij))     
                           s(ij)=1.0
                           time_ieos(ij) = days + time_ch
                        elseif(s(ij).le.0.0.and.days.ge.time_ieos(ij))
     &                          then
c     change to gas only
                           ieosdc=3
c                           t(ij)=psatl(pl,pcp(ij),dpcef(ij),
c     2                          dtsatp,dpsats,1,an(ij))*eosmg
c gaz debug 121817                           
c                           t(ij) = t(ij)*eosmg
                           t(ij) = t(ij)*eosmg
                           s(ij)=0.0
                           time_ieos(ij) = days + time_ch
                        endif
c     
                     elseif(ieosd.eq.3) then
c
c     gas conditions
c 
c gaz 072120 added ieos_sc() ,102521 added  ieosdc = 4
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieos_sc(ij) = 1
                        ieosdc = 4
c gaz 112121 should be s(ij) = 1 for pure water                     
c                        s(ij) = 0
                        s(ij) = 1.0
                      else if(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        ieos_sc(ij) = 3  
                        ieosdc = 1
                        s(ij) = 1.0
                      else 
c gaz 101719 modify if block (leave for now) 
c gaz 072120 pl now in if block                          
                       tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                        if(tl.le.tboil/phase_mult.
     &                       and.days.ge.time_ieos(ij)) then
c     change to 2-phase
                           s(ij)=satml         
                           t(ij)=tboil
                           ieosdc=2
                           time_ieos(ij) = days + time_ch
                           endif
                      endif     
                     elseif(ieosd.eq.4) then
c gaz 072120 ieosd.eq.4 equivalent to ieos_sc(ij).eq.1 (both SC)
c     sc conditions 
c 
c sufficient  to go to gas or liquid
                      if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieos_sc(ij) = 2  
                        ieosdc = 3
                        s(ij) = 0.0
                      else if(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        ieos_sc(ij) = 3   
                        ieosdc = 1
                        s(ij) = 1.0
                      else if(pl.lt.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                       tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                       if(tl.ge.tboil) then
                        ieosdc = 3
                        s(ij) = 0.0
                        ieos_sc(ij) = 0
                       else
                        ieosdc = 1
                        s(ij) = 1.0
                        ieos_sc(ij) = 0
                       endif 
                     else
                      ieosdc = 4
                      ieos_sc(ij) = 1
                     endif
                   endif
c
c
c     remember if danl eos change occured
c     tally eos numbers
c                                    
                     if(ieosd.ne.ieosdc) then
                        strd    =stepl
                        ieos_ch(ij) = ieos_ch(ij) +1
c  gaz connented out 111621                      if (ieos_ch(ij).gt.3) 
c     &                       write (ierr,233) 
c     &                       ij, ieos_ch(ij), (cord(ij,k),k=1,3)
                      endif
                     ieos(ij)=ieosdc
                  enddo
c 233              format('$$$$$$$  >>>>>> ', 2i8,1p,3g14.4)
          else if(icarb.ne.0) then
            call icectrco2(-1,0)
c            ieos = 1
          endif 
c     
         else if(ico2.gt.0.and.ivar_mass.gt.0) then
c gaz 082622 if nr_completed gt.0 then call  
c gaz 082922 nr_completed not working right             
c          if(nr_completed.gt.0) then
         if(l.gt.0) then
           call varchk_AWH(1)
c            nr_completed = 0
         endif
         else if(ico2.gt.0) then
c
c     determine phase state for water and noncondensible
c
c gaz 072020 modified code to use partial pressures or pci
c gaz 111520   next 2 lines only for debuging   
                ij = fdum+l+iad
c gaz 122920 allocate memory and initialize variables in explicit phase change
                  call explicit_phase_update(0,0)
                  if(.not.allocated(descrip)) then
                   allocate(descrip(neq,4))
                   descrip = '              '
                  endif
                  nr1=nrhs(1)
                  nr2=nrhs(2)
                  nr3=nrhs(3)
                  do i=1,neq
                     ij=i+ndummy
                     pl=phi(ij)
                     pci0 = pcio(ij)
                     pcl=pci(ij)
                     pci(ij)=pcl
                     x=pl-pcl
                     tl=t(ij)
                     sl=s(ij)
c gaz 040423 
                  ieosd=ieos(ij)      
                  s0 = sl
                  pc0 = pcl
                  if(ieosd.ne.2) then
                   p0 = bp(ij+nr1)
                   t0 = bp(ij+nr2)
                   pc0 = bp(ij+nr3)               
                  else
                   p0 = bp(ij+nr1)
                   t0 = bp(ij+nr3)
                   s0 = bp(ij+nr2)                    
                  endif  
                     if(ieosd.eq.4) then
                      if(tableFLAG.eq.0.and.iwater_table.eq.0) then
c gaz 112020 stop if SC phase and no table write error msg and stop
                      if(iptty.ne.0) write(iptty,30)
                      if(iout.ne.0) write(iout,30) 
                      write(ierr,30)       
                      endif   
                     endif                     
                     if(ieosd.eq.0) ieosd=1
                     ieosdc=ieosd
c     
c     check if eos change is necessary
c
                     if(ps(ij).eq.0.0.or.idof.le.1) then
c                 ieosdc=-1 gaz 10-17-2001 do nothing
c gax 120818 need a goto statement here???                         
                     endif
c 
                     if(ieosd.eq.1) then
c
c     liquid only state
c
c gaz 120818 add transition to supercritical state
c SC state occurs when total pressure (pl here)  reaches Pcrit 
c pure water uses total pressure and same here because it is liquid phase
c or leaving liquid phase
c gaz 072020 pl changed to  pl-pcl whwn comparing to pl 

c                          
                     tboil = tcrit_h2o
c gaz 072020 pl to pl-pcl, finished 072020 with ioes_sc()                
                      if(pl-pcl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o
     &                  .and.pcl.lt.pl) then
                        ieosdc = 4
                        ieos_sc(ij) = 3
                        descrip(ij,1) = '1 to 4 pcrit'
                        go to 101
c gaz 110220 gaz added ifblock for physical 2 phase                        
                      else if(pl-pcl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o
     &                  .and.pcl.lt.pl) then
                        ieosdc = 3 
                        ieos_sc(ij) = 3
                        s(ij) = 0.0
                        descrip(ij,1) = '1 to 3 pcrit'
                        go to 101
                       else if(pcl.ge.pl) then
c gaz 040823 
                         sl =1.0d0                              
                         s(ij) = sl
c gaz 042423   (varchk_12d)
c                         phi(ij) = pcl
c                         ieosdc = 2
                         ieosdc = 2
                         ieos_sc(ij) = 1
                         pci(ij) = pcl*wgtchng1 + (1.d0-wgtchng1)*pc0
                         t(ij)  = tl*wgtchng1 + (1.d0-wgtchng1)*t0
                         descrip(ij,1) = '1 to 2 pcl>=pl'   
                        go to 101
                   else if(pl.lt.pcrit_h2o) then
c gaz 121918 pci0 instead of pcl   
c use total pressure because air is dissolved 
c gaz 091920 is boiling point changed with dissolved air
                         tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))        
                   endif                       
                       pboil=psatl(tl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             0,an(ij))

c     change to 2-phasec                        
c                           if((x.le.pboil/phase_mult.or.
c     &                        tl.gt.tboil*phase_mult).
                        
c                          if(tl.gt.tboil*phase_mult. 
                           if(tl.gt.tboil.
     &                     and.days.ge.time_ieos(ij)) then                              
                           ieosdc=2
                           sl = 1.0
                           ieos_sc(ij) = 1
c gaz 100420 change tboil based on pcl                        
                         tboil=psatl(pl-pcl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))                        
                           t(ij) = tboil
                           s(ij) = sl
                           pci(ij) = pcl

                           descrip(ij,1) = '1 to 2 mass chk'
                           time_ieos(ij) = days + time_ch
                        endif
101                  continue  
                     endif
c 
                     if(ieosd.eq.2) then
c
c  
c gaz 070720 pl changed pl-pcl 
                      ieos_sc(ij) = 2   
                      if(pl-pcl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                        ieos_sc(ij) = 4
                        s(ij) = 1.
                        go to 201
                      else if(pl-pcl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o)
     &                   then
c gaz_101220 (iflg = 2 is for phase 2 to 3)
                        call phase_chk(2,ij,pcl,sl,ieosd,'pcl')   
                        ieosdc = 3 
                        pci(ij) = pcl
                        ieos_sc(ij) = 3
                        s(ij) = 0.0
                        descrip(ij,2) = '2 to 3 mass chk'
                        go to 201                      
c                     else if(pcl.lt.-0.0001) then 
                       else if(pcl.lt.0.0) then 
c gaz 122018 another phase change (physical ?)                         
                        ieosdc =3
                        s(ij) = 0.0 
                        ieos_sc(ij) = 3
                        descrip(ij,2) = '2 to 3 pcl<0'
                        go to 201
                      else if(pl-pcl.gt.pcrit_h2o.and.tl.lt.tcrit_h2o)
     &                   then      
                        ieosdc = 2
                        ieos_sc(ij) = 1
                       else if(pl-pcl.lt.pcrit_h2o) then                         
                         tboil=psatl(pl-pcl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))
                      endif                       
c     2-phase conditions
                        if(sl.ge.1.and.days.ge.time_ieos(ij)) then 
c     change to liquid only conditions
                           pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
c gaz 011116 testing
                         if(sl.ge.1.000) then
                           ieosdc=1
                           ieos_sc(ij) = 1
                           descrip(ij,2) = '2 to 1 sl>=1'
                           if(pl-pcl.lt.pboil) then
c gaz040823
c                            pcl = pl-pboil
                            ieosdc = 2
                            ieos_sc(ij) = 2
                            pci(ij)=pcl 
c gaz 111120
c gaz 040823
                             s(ij) = 1.0d0
                             descrip(ij,2) = '2 to 2 pl-pcl.lt.pboil'
                            go to 201
                           else
c gaz 040823
c                            pcl = pl-pboil
                            pci(ij)=pcl 
                            ieosdc=1
                            ieos_sc(ij) = 1
                            s(ij) = 1.0
                            descrip(ij,2) = '2 to 1 pl-pcl.ge.pboil'
c                            strd = stepl
                           endif                        

                           time_ieos(ij) = days + time_ch
                         endif
                        endif 
c end new addition on 092020
c change to gas only
                        if(sl.le.0.00.and.days.ge.time_ieos(ij))then
c gaz debug 090515
                           pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
 
                          if(sl.le.0.0) then                                  

                           ieosdc=3 
                           ieos_sc(ij) = 3
                           s(ij)=0.0
                           descrip(ij,2) = '2 to 3 sl<=0.0'
                           strd = stepl


                          endif
                        endif
c gaz debug 120714
                        if(x.lt.-xdiff_tol.and.sl.le.0.0.and.ieosdc.eq.2
     &                     .and.days.ge.time_ieos(ij)) then
                           ieosdc=3
                           ieos_sc(ij) = 3
                           s(ij)=0.0
                           phi(ij) = pcl
                           descrip(ij,2) = '2 to 3 x.. p = pc'
                           strd = stepl
                           time_ieos(ij) = days + time_ch
                        endif
                        pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
                        if(x.lt.pboil.and.sl.le.1.e-8.and.ieosdc.
     &                     eq.2.and.days.ge.time_ieos(ij)) then
                           ieosdc=3
                           ieos_sc(ij) = 3
                           s(ij)=0.0
                           t(ij) = 1.0001*tl 
                           descrip(ij,2) = '2 to 3 x.. p = pc'
c                           phi(ij) = pci(ij) + pboil*0.999
c gaz 012819                            
c                            pci(ij) = phi(ij) - pboil*0.999
                            phi(ij) = pci(ij) + pboil
                           strd = stepl
                           time_ieos(ij) = days + time_ch
                        endif
                     endif
c 
201                  continue
c                     if(ieosd.eq.3.or.ieosdc.eq.3) then
c gaz debug 120814
                     if(ieosd.eq.3) then

c gaz 120818 add SC transition logic (use pv (x here) instead of total pressure

                        ieos_sc(ij) = 3
                        if(pcl.ge.pl) then
 
                         pci(ij) = pcl
                         continue
                        endif
                        if(pl-pcl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                         ieosdc = 4
                         ieos_sc(ij) = 3
                         s(ij) = 0.0
                         go to 500
                        elseif(pl-pcl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o)
     &                      then
                         ieosdc = 2
                         ieos_sc(ij) = 1
                         s(ij) = 1.0
                         descrip(ij,3) = '3 to 2 pl-pcl>=pcr'
                         go to 500
                        elseif(pl-pcl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o)
     &                      then
                          ieosdc = 3
                          ieos_sc(ij) = 3
                          pvapor = pcrit_h2o
                          go to 500
                        else
                         pvapor=psatl(tl,pcp(ij),dpcef(ij),dpsatt,
     &                           dpsats,0,an(ij))                     
                        endif  
               
c     check vapor pressure against saturated vapor pressure
c     change if lower
                       if(x.ge.pvapor*phase_mult.and.days.
     &                    ge.time_ieos(ij)) then 
                           s(ij)=satml  
                           tboil=psatl(x,pcp(ij),dpcef(ij),dtsatp,
     &                            dpsats,1,an(ij))
c                        
c gaz 101819 
c  
c gaz debug 102420                           
c 
                           t(ij) = tboil
                           s(ij) = sl
                           pci(ij) = pl-pvapor
                           time_ieos(ij) = days + time_ch
                           ieosdc = 2
                           ieos_sc(ij) = 2
                           descrip(ij,3) = '3 to 2 phas_chk'
                      endif
500                  continue
                     endif  
c
c     sc conditions 
c                   
                  if(ieosd.eq.4) then 
c gaz 101920 check for pcl.ge.pl and estimate new pcl
                      ieos_sc(ij) = 3
                      if(pcl.ge.pl) then
                       call phase_chk(5,ij,pcl,sl,ieosd,'pcl') 
                       pci(ij) = pcl
                       continue
                      endif
                      if(pl-pcl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieos_sc(ij) = 3
                        ieosdc = 3
                        s(ij) = 0.0
                      else if(pl-pcl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o)
     &                     then
c gaz 102520 phase_chk returns a saturation (no air solubility)   
c check if dissolved in liquid (if so,ieosdc =1, if not, then ieosdc =2)
                         pvapor=psatl(tl,pcp(ij),dpcef(ij),dpsatt,
     &                           dpsats,0,an(ij))   
c gaz following call saves 7 % iterations on 400 day problem               
                        ieosdc = 2
                        ieos_sc(ij) = 1                       
                        s(ij) = max(sl,0.0)+0.001
c                        pci(ij) = pl-pvapor
c gaz 111120 added phi(ij) = ...                        
c                        phi(ij) = pcl+pvapor (did not work)
c                        strd = 0.98 is best
                        strd = 0.98
                      else if(pl-pcl.lt.pcrit_h2o.and.tl.lt.tcrit_h2o)
     &                  then
                        tboil=psatl(pl-pcl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))
                       if(tl.ge.tboil) then
                        ieosdc = 3
                        ieos_sc(ij) = 3
                        s(ij) = 0.0
                       else
c gaz 111220 trying to this 4-3 or 4-2   
c phase_chk(6,ij,pcl,sl,ieosd,'sl ')  estimates sl assuming no liquid solubility                           
c                        call phase_chk(6,ij,pcl,sl,ieosd,'sl ')  
                        ieos_sc(ij) = 2
                        ieosdc = 2
c gaz 111520   commented out next line                     
c                        pci(ij) = pl-pvapor
c                        s(ij) = 0.
c                        s(ij) = sl
                        s(ij) = max(sl,0.0)+0.01
                       endif 
                      else
                       ieosdc = 4
                       ieos_sc(ij) = 3
                      endif                     
                  endif  
c
c     remember if danl eos change occured
c     tally eos numbers
c 
c gaz 012919 **************** very new limiting change to pci*********************************   
                     pci(ij)=max(0.0d00,pci(ij))
c gaz 040622                     
c                     pci(ij)=min(phi(ij),pci(ij))
                      if(phi(ij).lt.0.0d0) phi(ij) = 0.0d0                    
                     if(ieosd.ne.ieosdc) then
                        iphase_chng = iphase_chng +1
                        ieos_ch(ij) = ieos_ch(ij) +1
c gaz 122920 save phase change node info (countchanges                    
                      call explicit_phase_update(1,ij)
                     endif
                                  
c gaz 040323 print out detailed phase change information
                if(ipr_phase_detail.ne.0) then
                 if(iad.eq.0) then
                  chariad(1:6) = '      '
                  fdum_0 = -1
                 else
                  if(fdum_0.gt.fdum) then
                   chariad(1:6) = 'lower '
                   fdum_0 = fdum
                  else
                   chariad(1:6) = 'higher'
                   fdum_0 = fdum
                  endif
                 endif
                  if(ij.eq.1) then
                   if(iad.gt.maxit) then
                    chariad(7:12) = ' MAXIT'
                   else
                    chariad(7:12) = '      '
                   endif
                   write(ierr,555) l, day, days, iad, iphase_chng, fdum,
     &              chariad
                   if(iphase_chng.ge.0)  write(ierr,556)
                  endif
                 if(descrip(ij,ieosd)(1:13).ne.'            ') then
                  write(ierr,557) ij,ieosd,ieosdc,p0,phi(ij),t0,t(ij), 
     &            pc0,pci(ij),s0,s(ij),descrip(ij,ieosd)
                 endif
c                  if(ij.eq.neq) write(ierr,558) iphase_chng, strd
555     format('time step ',i6,' day ',1p,g14.6,' days ',1p,g14.6,
     &   ' iteration ',i5,' phase changes: ',i5,1x,'fdum ',g14.6,1x,a12)
556     format('node',t9,'ieos0',t15,'ieos',t21,'press0',t35,'press',
     &         t49,'t0',t63,'t',t77,'pc0',t91,'pci',t105,'s0',t119,'s') 
557     format(i7,t9,i6,t15,i6,t21,1p,g14.5,t35,g14.5,
     &         t49,g14.5,t63,g14.5,t77,g14.5,t91,g14.5,
     &         t105,g14.5,t119,g14.5,t133,1x,a22)
558    format('phase changes: ',i5,1x,' strd ',0p,f9.4)
                  endif  
                     ieos(ij)=ieosdc                       
                  enddo  
c gaz 050323 print out nodal masses and energy, 050523
                call  awh_accumulation_calc(1,1,neq,0)
                call  awh_accumulation_calc(4,0,0,0)
c gaz 122920 added explicit update
c this solves for new gridblock mass and energy conservation 
                call explicit_phase_update(41,0)
                call explicit_phase_update(2,0)
                call explicit_phase_update(42,0)  
102       continue                
               else if(ico2.lt.0.and.ice.eq.0)then
c
c     determine phase state for isothermal air-water flow
c
                  call airctr(1,ndummy)
               else if(ico2.lt.0.and.ice.ne.0)then
c
c     determine phase state for low-temperature solid-liquid-gas system
c
                  call icectr(1,ndummy)
c RJP 04/10/07 modified for CO2
               else if(icarb.eq.1) then
                  call icectrco2(1,ndummy)
                  call icectrco2(-34,ndummy)
               endif
c end block for idof.ne.1
            endif
c
c     call eos routines
c
c
c     one call to vcon to explicity update thermal conductivities
c
            if (iad.eq.0) call vcon(1,ndummy)
c
            if(ico2.gt.0) then
c gaz 120821 make sure all properties are evaluated for all nodes initially                
              if(l.eq.0) then  
               eval_test_h2o_0 = eval_test_h2o
               eval_test_h2o = 1                
                call thrmwc(ndummy)
               eval_test_h2o = eval_test_h2o_0 
              else
                call thrmwc(ndummy)
              endif
            else if(ico2.eq.0) then
c RJP 04/10/07 
             if(icarb.eq.1) then
               call icectrco2(-34,ndummy)
c              call icectrco2(3,ndummy)
c              call icectrco2(-3,ndummy)
               call icectrco2(33,ndummy)
               call icectrco2(-35,ndummy)
c gaz 120821 make sure all properties are evaluated for all nodes initially                
              else if(l.eq.0) then  
               eval_test_h2o_0 = eval_test_h2o
               eval_test_h2o = 1                
                call thermw(ndummy)
               eval_test_h2o = eval_test_h2o_0 
              else
                call thermw(ndummy)
              endif

            else if(ico2.lt.0.and.ice.eq.0)then
               call airctr(3,ndummy)
            else if(ico2.lt.0.and.ice.ne.0)then
               call icectr(-34,ndummy)
               call icectr(3,ndummy)
               call icectr(-3,ndummy)
               call icectr(33,ndummy)
               call icectr(-35,ndummy)
            end if
c     
c
c      NR corrections
c
         else
            if(ico2.eq.0) then
c
c     pure water
c
c           strd is passed through common
               nr1=nrhs(1)
               nr2=nrhs(2)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  ieosd=ieos(i)
c gaz 112115
                  if(ieosd.eq.4) ieosd = 1
                  if(ps(i).eq.0.0.or.ieosd.eq.0) then
c gaz 10-18-2001
                     t(i)=t(i)-bp(i2)*strd
                  elseif(ieosd.eq.1) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                  elseif(ieosd.eq.2) then
                     phi(i)=phi(i)-bp(i1)*strd
                     s(i)=s(i)-bp(i2)*strd
                  elseif(ieosd.eq.3) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                  endif
               enddo    
               else if(ico2.gt.0) then
c
c     NR corrections for water and noncondensible
c
c gaz 032722 need to save true variable change (with strd etc)            
c gaz 082622 added variable nr_completed                   
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
c  gaz   041323 
               descrip = '                      '
               iphase_chng = 0   
c gaz 043023 calculate  mass water, energy, mass ngas
             if(imass_chk.eq.0.and.l.eq.1)then
c check  if mass check file exists
                call awh_accumulation_calc(3,0,0,0)
                if(imass_chk.ne.0) call awh_accumulation_calc(0,1,neq,0)
             endif
             if(imass_chk.ne.0.and.iad.gt.0)then
c gaz 051023 
c               call awh_accumulation_calc(2,1,neq,0)
               call awh_accumulation_calc(0,1,neq,0)
c gaz 050523 move after phase chng
c               call awh_accumulation_calc(1,1,neq,0)
             else if(imass_chk.eq.0) then
               
             endif
          
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  ieosd=ieos(i)
c gaz 032722                  
                  pl_last = phi(i )
                  tl_last = t(i)
                  pcl_last = pci(i) 
                  sl = s(i)
                  if(ps(i).eq.0.0) then
c     gaz 10-18-2001
                     t(i)=t(i)-bp(i2)*strd
c gaz 121718 added ieosd 4 (SC)                  
                  elseif(ieosd.eq.1.or.ieosd.eq.4) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                     pci(i)=pci(i)-bp(i3)*strd
                  elseif(ieosd.eq.2) then
                     phi(i)=phi(i)-bp(i1)*strd
                     s(i)=s(i)-bp(i2)*strd
c  GAZ 5/1/98          
                     t(i)=t(i)-bp(i3)*strd
c gaz 122418                      
                         pvapor=psatl(t(i),pcp(i),dpcef(i),dpsatt,
     &                           dpsats,0,an(i))  
c gaz debug 093019(stays in! )                      
                     if(pvapor.ge.phi(i)) then
                          t(i)=psatl(phi(i),pcp(i),dpcef(i),
     &                          dtsatp,dpsats,1,an(i)) 
                        pvapor = phi(i)
                     endif     
c gaz 041023                    pci(i) = phi(i)-pvapor
                  elseif(ieosd.eq.3) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                     pci(i)=pci(i)-bp(i3)*strd
                  endif
c gaz 041323 
                   if(phi(i).lt.0.0d00) then
                    if(ieosd.ne.1) iphase_chng = iphase_chng +1
                    if(ieosd.ne.1) descrip(i,ieosd)
     &                     = '2 to 1 pl <= 0 NR '   
                    if(ieosd.eq.1) descrip(i,ieosd)
     &                     = '1 to 1 pl <= 0 NR ' 
                    ieos(i) = 1
                    ieosd = 1
                    t(i) = tl_last
                    phi(i) = pl_last  
                   endif 
c gaz 040923 changed min-max to include phase state
c gaz 041323   ' lt instead of le ',  ' gt instead of ge '    OK     
c gaz 041923 changing for pci < 0      
                  if(pci(i).lt.0.0d00) then
                   iphase_chng = iphase_chng +1
                   if(ieosd.eq.1)  then
                    pci(i) = 0.0d0
                    descrip(i,ieosd) =  '1 to 1 pcl <= 0 NR '
                   else if(ieosd.eq.2)  then
                    pci(i) = 0.0d0
                    ieos(i) = 1
                    descrip(i,ieosd) =  '2 to 1 pcl <= 0 NR '          
                   else
                    pci(i) = 0.0d0
                    ieos(i) = 3
                    descrip(i,ieosd) =  '3 to 3 pcl <= 0 NR '
                    s(i) =0.0d0
                   endif
                  endif
                  if(pci(i).gt.phi(i)) then
c                   if(ieosd.ne.1) iphase_chng = iphase_chng +1
c gaz 041923
c                   ieos(i) = 1                                     
c                   ieosd = 1
                   pci(i) = phi(i)
                   if(ieosd.eq.2) descrip(i,ieosd)
     &                     = '2 to 2 pcl >= pl NR '   
                  if(ieosd.eq.1) descrip(i,ieosd)
     &                     = '1 to 1 pcl >= pl NR ' 
                  if(ieosd.eq.3) descrip(i,ieosd)
     &                     = '3 to 3 pcl >= pl NR ' 
                  endif
c gaz 041323   ' lt instead of le ',  ' gt instead of ge '     OK
                  if(s(i).lt.0.0d00) then
                   if(ieosd.ne.3) iphase_chng = iphase_chng +1
                   ieos(i) = 3
                   ieosd = 3
c gaz debug 042623
                   t(i) = (tl_last+t(i))*0.5d0
                   descrip(i,ieosd) = '2 to 3 sl <= 0 NR '
                   s(i) = 0.0d00
                  elseif(s(i).gt.1.0d00) then
                   if(ieosd.ne.1) iphase_chng = iphase_chng +1
                   ieos(i) = 1
                   ieosd = 1
                   descrip(i,ieosd) = '2 to 1 sl >= 1 NR '
                   s(i) = 1.0d00
                  endif
                  if(t(i).le.tNR_min) then
c                   t(i) = tl_last
                  endif
c gaz 032722  load variable changes into bp
                  if(ieosd.ne.2) then
                   bp(i1) = pl_last
                   bp(i2) = tl_last
                   bp(i3) = pcl_last                  
                  else
                   bp(i1) = pl_last
                   bp(i3) = tl_last
                   bp(i2) = sl                     
                  endif                  
               enddo 
c gaz 051023
                  call awh_accumulation_calc(2,1,neq,0)
c gaz 082622 added nr_completed = 3              
              nr_completed = 3
            else if(ico2.lt.0.and.ice.eq.0) then
c     
c     make corrections for isothermal air-water mixture
c     
               call airctr(2,ndummy)
            else if(ico2.lt.0.and.ice.ne.0) then
c     
c     make corrections for solid-liquid-gas mixture
c     
               call icectr(2,ndummy)
c RJP 04/10/07 added CO2 part
c this is now in bnswer 111410 (gaz)
c         else if(icarb.eq.1) then
c            call icectrco2(2,ndummy)
            endif
c
         endif


      return
      end
c gaz auxilliary subroutines 1. prop_phase_change 2. awh_accumulation_calc
      subroutine prop_phase_change(iflg,i,ieosd1,ieosd2,
     & tl_new,phi_new,pcl_new)
c
c calculate phase properties
c compare props for possible phase change   
c output new variables      
c
c gaz 122118
c
      use comai
      use comci
      use comdi
      use comfi
      implicit none
      integer iflg,i,ieosd1,ieosd2
      real*8 tl,pl,pcl,denpci0
      real*8 roc,drocpc,dtin
      real*8 tl_new,phi_new,pcl_new,pcl0
      if(iflg.eq.1) then
c check props for change from state 2 to state 3      
       if(ieosd1.eq.2.and.ieosd2.eq.3) then
         if(ps(i).gt.0.0) then
          tl = t(i)
          pl = phi(i)
          pcl = pci(i)
          dtin=1.0/dtot
          denpci0 = (denpci(i)/dtin+denpch(i))/ps(i)
          roc0=1.292864
          pcl0=0.101325 
          drocpc=roc0*(273./(tl+273.))/pcl0    
c roc is mass of air/vol (enough to compare)          
          roc=drocpc*pcl
          pcl_new = denpci0/drocpc
          continue         
         else
c error condition : no pore space             
         endif  
       endif
      else
      endif
      return
      end

      subroutine awh_accumulation_calc(iflg,i1,i2,ndummy)
c
c compar differences in accumulation terms
c
      use comai
      use combi,only : sx1,nelm,nelmdg
      use comci
      use comdi
      use comdti,only : n0
      use comei,only : a
      use comfi
      use comgi,only : bp
      use davidi,only : nmat,nrhs
      use com_exphase,only : i_ex_update, ieq_ex
      use com_prop_data
      implicit none
      real*8 dtin, ztol, energy_norm, frac_ngas, sx1d
      real*8 delmax(3), tol_eq1, tol_eq2, tol_eq3
      real*8 fdum_phase_calc, fdum_phase_prev,
     &        tol_phase_calc, fdum_raw(3)
      real*8 dtot_orig, ts_fac
      real*8 strd_smpl
      real*8 phi_low,phi_high,t_low,t_high,pci_low,pci_high
      real*8 s_low,s_high
      logical test_091322, exists
      integer ndex(3)
      integer inr, max_inr,ic_tol_count, inr_fail, ifail_phase_chk
c    
      integer i_ex_update_save, ieq_ex_save, ieos_save
      integer  ieos_in
      integer i, id, i1, i2, ij, iflg, neqp1, ndummy 
      real*8  weightavg, weight_21, dmpf_new, dmpf_old, dmpf_avg
      real*8 rolw_awh, s_over, drolw_awh_dp 
c gaz 050123
      character*80 accum_AWH_title
      integer open_file
      integer ieosd, ieosdc
      parameter (weight_21 = 1.0)
      parameter(ztol = 1.d-14, max_inr =10, tol_phase_calc = 1.d-6)
      parameter(ts_fac = 0.5d0)
      parameter(phi_low =-0.01,phi_high = 100., t_low = -1.)
      parameter(t_high = 1000., pci_low = -1, pci_high = 100.)
      parameter(s_low = -10, s_high = 10.)

       if(iflg.eq.0) then
c allocate memory  
        if(.not.allocated(zntotf)) allocate(zntotf(n0)) 
        if(.not.allocated(mass_h2o)) then 
         allocate(mass_h2o(n0))
         allocate(mass_ngas(n0))
         allocate(energy_tot(n0))
         allocate(mass0_h2o(n0))
         allocate(mass0_ngas(n0))
         allocate(energy0_tot(n0))
         allocate(dum_calc_eos(n0))
         allocate(dtot_test(n0))
         allocate(ieos0(n0))
c calculate initial nodal mass and energy totals (thrmwc already called )
           do i = 1, nawh
            id = nawh_nodes(i)
            sx1d = sx1(id)
            mass_h2o(id) = (deni(id)*dtot + denh(id))*sx1d
            mass_ngas(id) = (denpci(id)*dtot + denpch(id))*sx1d
            energy_tot(id) = (denei(id)*dtot + deneh(id))*sx1d
c
           enddo
        endif
      else if(iflg.eq.1) then
c
c gaz for phase change 
c
        neqp1 = neq+1
        do i = 1, nawh
         id = nawh_nodes(i)
         ieq_ex = id
         sx1d = sx1(id) 
          call thrmwc(0)          
         mass_h2o(id) = (deni(id)*dtot + denh(id))*sx1d
         mass_ngas(id) = (denpci(id)*dtot + denpch(id))*sx1d
         energy_tot(id) = (denei(id)*dtot + deneh(id))*sx1d
         frac_ngas = mass_ngas(id)/(mass_ngas(id)+mass_h2o(id)+ztol)
         zntotf(id) = frac_ngas
c
        enddo
      else if(iflg.eq.2) then
        do i = 1, nawh
         id = nawh_nodes(i)
         mass0_h2o(id)   = mass_h2o(id)
         mass0_ngas(id)  = mass_ngas(id)
         energy0_tot(id) = energy_tot(id)
         ieos0(id) = ieos(id)
        enddo
      else if(iflg.eq.3) then
c  read test nodes if file exists        
          imass_unit = 0
          inquire(file = 'acc_AWH_nodes.txt', exist = exists)
          if(exists)  then
           imass_chk = 1
           imass_unit = open_file('acc_AWH_nodes.txt', 'old')
           read(imass_unit,'(a80)')  accum_AWH_title
           read(imass_unit,*) nawh
           allocate(nawh_nodes(nawh))
           backspace imass_unit
           read(imass_unit,*) nawh,(nawh_nodes(i),i = 1, nawh)
           imout = open_file('acc_AWH_write.txt', 'unknown')
          else
           imass_chk = 0
          endif
          
      else if(iflg.eq.4) then
c      printout mass change
            do i = 1,nawh
               ij = nawh_nodes(i)
                 ieosd = ieos0(ij)
                 ieosdc = ieos(ij)
                  if(i.eq.1) then
                   write(imout,555) l, day, days, iad, iphase_chng,
     &              fdum, '            '
                   write(imout,556)
                  endif
c gaz 050423 debug
c                  
                  write(imout,557) ij,ieosd,ieosdc,mass0_h2o(ij),
     &            mass_h2o(ij),energy0_tot(ij),energy_tot(ij), 
     &            mass0_ngas(ij),mass_ngas(ij),descrip(ij,ieosd)(1:22)
c                endif
            enddo
555     format('time step ',i6,' day ',1p,g14.6,' days ',1p,g14.6,
     &   ' iteration ',i5,' phase changes: ',i5,1x,'fdum ',g14.6,1x,a12)
556     format('node',t10,'ieos0',t17,'ieos',t26,'mass0',t41,'mass',
     &         t52,'energy0',t67,'energy',t85,'mngas0',t100,'mngas',
     &         t114,'discription')
557     format(i7,t9,i6,t15,i6,t21,1p,g14.5,t35,g14.5,
     &         t49,g14.5,t63,g14.5,t77,g14.5,t91,g14.5,t114,a22)
c     &         t105,g14.5,t119,g14.5,t135,' description ',a22)
558    format('phase changes: ',i5,1x,' strd ',0p,f9.4)
       endif
      return 
      end