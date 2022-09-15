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
      use comtable, only : tableFLAG
      use davidi
  
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
      real*8 xdiff_tol
      parameter(psatmn=0.0001)
      parameter(eosmg=1.0001)
c      parameter(eosml=0.95)
c      parameter(eosml=0.99)
      parameter(eostol=0.0001)
c      parameter(stepl=0.95)
      parameter(pcimin=0.0)
c      parameter(phase_mult=1.00)
      parameter(phase_sat=1.0d-9)
c gaz debug 092115
c      parameter(satml=1.0d-4)
c      parameter(satml=1.0d-6)
c      parameter(phase_mult=1.01)
      parameter(xdiff_tol=1.0d-4)
      real*8 pcrit_h2o, tcrit_h2o
      parameter(pcrit_h2o=22.00d0, tcrit_h2o=373.95)
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
      real*8 stepl_hm, phase_mult_hm, satml_hm, eosml_hm
      real*8 stepl_hma, phase_mult_hma, satml_hma, eosml_hma
      parameter(stepl_hm = 0.95, phase_mult_hm = 1.00001)
c gaz debug 051516 (optimized for geothermal;satml = 1.0d-2 may be too large)      
      parameter(satml_hm = 1.0d-2, eosml_hm = 0.9)
      parameter(stepl_hma = 0.95, phase_mult_hma = 1.05)
      parameter(satml_hma = 1.0d-7, eosml_hma = 0.99)
c
c
c ich_m1 and ich_m2 are in comai but are adjusted here
c these are maybe best for geothermal systems
c
      if(ico2.eq.0) then
        satml = satml_hm
        phase_mult = phase_mult_hm
        eosml = eosml_hm
      else
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
      if(nr_stop.ne.0.and.iad.ge.1.and.ifl.eq.1)then
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
c                         tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
c     &                             1,an(ij))
                          tboil = 3000.   !JPO  test (Mars)
                          ieosdc = 1      !JPO  test (Mars)
c***** JPO 09/12/22
c      (Mars) Force liquid phase if tsa0 (from sther.f) = 1000
c                   ! else if(pl.lt.pcrit_h2o.and.tsa0.eq.1000) then
c                           do nothing (keep tboil=2000.) 
c                      ! tboil = 3000.
c                      ! ieosdc = 1
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
c     s(ij)=eosml
c gaz debug 121817                           
c                            s(ij)=0.999    
                           s(ij) = eosml
                           t(ij)=tboil
                           time_ieos(ij) = days + time_ch
c     call phase_balance_ctr(1,ij,phi(ij),t(ij),s(ij))
                        endif
c
                     elseif(ieosd.eq.2) then
c
c     2-phase conditions
c
c  gaz 102917 update temperature earlier in iteration                         
                           t(ij)=psatl(pl,pcp(ij),dpcef(ij),
     &                          dtsatp,dpsats,1,an(ij))  
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
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                        s(ij) = 1
                      else if(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        ieosdc = 1
                        s(ij) = 1.0
                      endif  
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
                     elseif(ieosd.eq.4) then
c
c     sc conditions 
c 
c sufficient  to go to gas or liquid
                      if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3
                        s(ij) = 0.0
                      else if(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        ieosdc = 1
                        s(ij) = 1.0
                      else if(pl.lt.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                       tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                       if(tl.ge.tboil) then
                        ieosdc = 3
                        s(ij) = 0.0
                       else
                        ieosdc = 1
                        s(ij) = 1.0
                       endif 
                     else
                      ieosdc = 4
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
                        if (ieos_ch(ij).gt.3) 
     &                       write (ierr,233) 
     &                       ij, ieos_ch(ij), (cord(ij,k),k=1,3)
                     endif
                     ieos(ij)=ieosdc
                  enddo
 233              format('$$$$$$$  >>>>>> ', 2i8,1p,3g14.4)
          else if(icarb.ne.0) then
            call icectrco2(-1,0)
c            ieos = 1
          endif 
c     
               else if(ico2.gt.0) then
c
c     determine phase state for water and noncondensible
c
c gaz debug 082714
c gaz debug added saltctr calls 091414
                  continue
                  i = l
                  i = fdum
                  do i=1,neq
                     ij=i+ndummy
                     pl=phi(ij)
                     pci0 = pcio(ij)
                     pcl=pci(ij)
c                     pcl=max(pcl,pcimin)
                     pci(ij)=pcl
                     x=pl-pcl
                     tl=t(ij)
                     sl=s(ij)
                     ieosd=ieos(ij)
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
                         
                     tboil = 2000.
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                        go to 101
                      else if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3 
                        s(ij) = 0.0
                        go to 101
                      else if(pcl.ge.pl) then
c gaz 121918  add physical phase change 
c gaz 01                          
                        ieosdc = 2
                        s(ij) = 0.9
                        strd = stepl
                        pboil=psatl(tl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             0,an(ij))
                        phi(ij)= pcl
                        pci(ij) = phi(ij)-pboil
                        go to 101
                   else if(pl.lt.pcrit_h2o) then
c gaz 121918 pci0 instead of pcl   
c use total pressure
                         tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))
                      endif                       
                       pboil=psatl(tl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             0,an(ij))

c     change to 2-phase
c gaz 121118 and 121818(changed to "t"
c                        if((x.le.pboil/phase_mult.or.
c     &                        tl.gt.tboil*phase_mult).
                        
                          if(tl.gt.tboil*phase_mult.  
     &                     and.days.ge.time_ieos(ij)) then
                           ieosdc=2
c                           strd=stepl
                           strd = 1.
c                           pboil=psatl(tl,pcp(ij),dpcef(ij),dtsatp,
c     &                             dpsats,0,an(ij))
c gaz 121218 and 121818                          
c                           pci(ij)=(max(pl-pboil,0.0d00)+pcl)/2.
                           s(ij)=0.95
                           time_ieos(ij) = days + time_ch
                        endif
101                  continue  
                     endif
c 
                     if(ieosd.eq.2) then
c
c                                  
                      if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 4
                        s(ij) = 1.
                        go to 201
                      else if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3 
                        s(ij) = 0.0
                        go to 201
                      else if(pcl.lt.-0.0001) then 
c gaz 122018 another phase change                          
                        ieosdc =3
                        s(ij) = 0.0                               
                        call prop_phase_change(1,ij,2,3,
     &                  tl,pl,pcl)
                        go to 201
                      else if(pl.lt.pcrit_h2o) then
c gaz 121918 pcl0 instead of pcl                          
                         tboil=psatl(pl-pcl,pcp(ij),dpcef(ij),dtsatp,
     &                             dpsats,1,an(ij))
                      endif                       
c     2-phase conditions
                        if(sl.ge.1.and.days.ge.time_ieos(ij)) then 
c     change to liquid only conditions
                           pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
c gaz 011116 testing
c                         if(sl.ge.1.000.and.so(ij).gt.0.95) then
                         if(sl.ge.1.000) then
                           ieosdc=1
                           pci(ij)=pl-pboil
                           s(ij)=1.0
                           strd = stepl
                           time_ieos(ij) = days + time_ch
c                         else
c                           s(ij) = 0.999
                         endif
                        endif 
c change to gas only
                        if(sl.le.0.00.and.days.ge.time_ieos(ij))then
c gaz debug 090515
                           pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
c                          if(sl.le.0.0.and.so(ij).lt.0.05) then   
                          if(sl.le.0.0) then         
                           denei_ch(ij) = denei(ij)
                           deni_ch(ij) = deni(ij)
                           denpci_ch(ij) = denpci(ij)
                           ieos_bal(ij) = ieos(ij)
                           dum = so(ij)
                           ieosdc=3 
                           s(ij)=0.0
                           strd = stepl
                           time_ieos(ij) = days + time_ch
c                          else
c                           s(ij) = 0.0001
                          endif
                        endif
c gaz debug 120714
                        if(x.lt.-xdiff_tol.and.sl.le.0.0.and.ieosdc.eq.2
     &                     .and.days.ge.time_ieos(ij)) then
                           ieosdc=3
                           s(ij)=0.0
                           phi(ij) = pcl
                           strd = stepl
                           time_ieos(ij) = days + time_ch
                        endif
                        pboil = psatl(tl,pcp(ij),dpcef(ij),
     &                          dpsatt,dpsats,0,an(ij))
                        if(x.lt.pboil.and.sl.le.1.e-8.and.ieosdc.
     &                     eq.2.and.days.ge.time_ieos(ij)) then
                           ieosdc=3
                           s(ij)=0.0
                           t(ij) = 1.0001*tl 
c                           phi(ij) = pci(ij) + pboil*0.999
                            pci(ij) = phi(ij) - pboil*0.999
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
c gaz 121118 use pl for liquid phase(s)
c                       if(x.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        if(pl.ge.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                         ieosdc = 4
                         s(ij) = 1
                         go to 500
c                      else if(x.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        elseif(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                         ieosdc = 1
                         s(ij) = 1.0
                         go to 500
                        elseif(tl.ge.tcrit_h2o) then
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
c                           t(ij) = 0.5*(tboil+tl)
c                           if(ieosd.eq.3) then
c                             s(ij)= satml
c                           endif
                           time_ieos(ij) = days + time_ch
                           strd = stepl
                           ieosdc = 2
                      endif
500                  continue
                     endif
c gaz 121118 add SC   
c
c     sc conditions 
c 
c sufficient  to go to gas or liquid
                  if(ieosd.eq.4) then 
                      if(pl.lt.pcrit_h2o.and.tl.ge.tcrit_h2o) then
                        ieosdc = 3
                        s(ij) = 0.0
                      else if(pl.ge.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                        ieosdc = 1
                        s(ij) = 1.0
                      else if(pl.lt.pcrit_h2o.and.tl.lt.tcrit_h2o) then
                       tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                             1,an(ij))
                       if(tl.ge.tboil) then
                        ieosdc = 3
                        s(ij) = 0.0
                       else
                        ieosdc = 1
                        s(ij) = 1.0
                       endif 
                      else
                       ieosdc = 4
                      endif                     
                   endif  
c
c     remember if danl eos change occured
c     tally eos numbers
c 
                     if(ieosd.ne.ieosdc) then
c gaz debug 120814
c                        strd    =stepl
c                        if(ieosdc.eq.2) strd = 1.
                        ieos_ch(ij) = ieos_ch(ij) +1
                        
c                        if (ieos_ch(ij).gt.3) 
c     &                       write (iout,233) 
c     &                       ij, ieos_ch(ij), (cord(ij,k),k=1,3)
                     endif
                     ieos(ij)=ieosdc
                  enddo    
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
               call thrmwc(ndummy)
            else if(ico2.eq.0) then
c RJP 04/10/07 
            if(icarb.eq.1) then
               call icectrco2(-34,ndummy)
c              call icectrco2(3,ndummy)
c              call icectrco2(-3,ndummy)
               call icectrco2(33,ndummy)
               call icectrco2(-35,ndummy)
               else
                  call thermw(ndummy)
               end if
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
c     strd is passed through common
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  ieosd=ieos(i)
                  if(ps(i).eq.0.0) then
c     gaz 10-18-2001
                     t(i)=t(i)-bp(i2)*strd
c gaz 121718 added ieosd 4                     
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
c                     tboil=psatl(phi(i),pcp(i),dpcef(i),dtsatp,
c     &                            dpsats,1,an(i))
                         pvapor=psatl(t(i),pcp(i),dpcef(i),dpsatt,
     &                           dpsats,0,an(i))                      
                     pci(i) = phi(i)-pvapor
                  elseif(ieosd.eq.3) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                     pci(i)=pci(i)-bp(i3)*strd
                  endif
c big change gaz 11/26/96
c gaz debug 120714 (phi correction)
c                  if(ieosd.ne.2.and.phi(i).lt.pci(i)) then
c                   if(phi(i).lt.pci(i)) phi(i) = pci(i)
c gaz debug 011415 try different approach pci correction
c gaz debug 011415 this worked well (so far)
c
c                   pci(i) = phi(i)
c                   strd = stepl
c                  endif
                  pci(i)=max(0.0d00,pci(i))
                  pci(i)=min(phi(i),pci(i))
                  s(i)=min(1.d00,max(0.0d00,s(i)))
               enddo      
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
