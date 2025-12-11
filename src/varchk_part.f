      subroutine varchk_part(ifl,ndummy,zone)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************

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

      use comdti
      use comai
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comii
      use davidi
      implicit none

      real*8 psatmn
      real*8 eosmg
      real*8 eosml
      real*8 eostol
      real*8 stepl
      real*8 pcimin
      real*8 phase_mult
      real*8 phase_sat 
      real*8 satml 
      parameter(psatmn=0.0001)
      parameter(eosmg=1.0001)
      parameter(eosml=0.999)
      parameter(eostol=0.0001)
      parameter(stepl=0.95)
      parameter(pcimin=0.0)
      parameter(phase_mult=1.00)
      parameter(phase_sat=1.0d-9)
      parameter(satml=1.0d-5)
      real*8 psatl
      integer i
      integer zone
      real*8 pl
      real*8 tl
      real*8 sl
      real*8 x
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
      integer ndummy
c
c determine variable set
c
      if(ifl.eq.0) then
c
c     evaluate eos status of each node. change status as prescribed
c     by phase change instructions
         
         
         
c set newton raphson step length to 1.0 at beginning of timestep
        if(iad.eq.0) strd    =1.0
cc
c
c new loop to check on degree of freedom 
c if idof=1, the problem is either heat conduction or saturated only
c
        if(idof.ne.1) then
         if(ico2.eq.0.and.ice.ge.0) then
c
c     determine phase state for pure water
c
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
                  tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                        1,an(ij))
c     change to 2-phase
                  if(tl.ge.tboil*phase_mult) then
                     ieosdc=2
                     s(ij)=eosml
                     t(ij)=tboil
                  endif
c
               elseif(ieosd.eq.2) then
c
c     2-phase conditions
c
                  if(sl.ge.1.) then
c     change to liquid only conditions
                     ieosdc=1
                     t(ij)=psatl(pl,pcp(ij),dpcef(ij),
     2                    dtsatp,dpsats,1,an(ij))*eosml
                     s(ij)=1.0
                  elseif(s(ij).le.0.0) then
c     change to gas only
                     ieosdc=3
                     t(ij)=psatl(pl,pcp(ij),dpcef(ij),
     2                    dtsatp,dpsats,1,an(ij))*eosmg
                     s(ij)=0.0
                  endif
c     
               elseif(ieosd.eq.3) then
c
c     gas conditions
c     
                  tboil=psatl(pl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                        1,an(ij))
                  if(tl.le.tboil/phase_mult) then
c     change to 2-phase
                     s(ij)=satml         
                     t(ij)=tboil
                     ieosdc=2
                  endif
               endif
c
c
c     remember if danl eos change occured
c     tally eos numbers
c     
               if(ieosd.ne.ieosdc) then
                  strd    =stepl
               endif
               ieos(ij)=ieosdc
            enddo     
c     
         else if(ico2.gt.0.and.ice.eq.0) then
c
c     determine phase state for water and noncondensible
c
            do i=1,neq
               ij=i+ndummy
               pl=phi(ij)
               pcl=pci(ij)
               pcl=max(pcl,pcimin)
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
               endif
c 
               if(ieosd.eq.1) then
c
c     liquid only state
c
                  tl=t(ij)
                  pboil=psatl(tl,pcp(ij),dpcef(ij),dtsatp,dpsats,
     &                        0,an(ij))
c     change to 2-phase
                  if(x.le.pboil/phase_mult) then
                     ieosdc=2
                     pci(ij)=max(pl-pboil,0.0d00)
                     s(ij)=eosml
c
                  endif
               endif
c 
               if(ieosd.eq.2) then
c
c     2-phase conditions
                  if(sl.ge.1.) then
c     change to liquid only conditions
                     ieosdc=1
c                    tl=t(ij)*eosml
c                    t(ij)=tl
                     pci(ij)=pl-psatl(tl,pcp(ij),dpcef(ij),
     2                    dpsatt,dpsats,0,an(ij))
c                    pci(ij)=max(pl-pboil,0.0d00)
                     s(ij)=1.0
                  endif
c     change to gas only
                  if(sl.le.0.0) then
                   ieosdc=3
                   s(ij)=0.0
                   strd    =stepl
                  endif
               endif
c 
               if(ieosd.eq.3.or.ieosdc.eq.3) then
c
c     gas conditions
c
c                 pvapor=psatl(to(ij),pcp(ij),dpcef(ij),dpsatt,dpsats,0)
                  pvapor=psatl(tl,pcp(ij),dpcef(ij),dpsatt,dpsats,
     &                         0,an(ij))
c     check vapor pressure against saturated vapor pressure
c     change if lower
                  if(x.ge.pvapor) then
                     s(ij)=0.0   
                     if(ieosd.eq.3) then
                      s(ij)= satml
                     endif
                     ieosdc=2
                     pci(ij)=max(pl-pvapor,pcimin)
                     strd    =stepl
                  endif
               endif
c     
c
c     remember if danl eos change occured
c     tally eos numbers
c
               if(ieosd.ne.ieosdc) then
c                 pnx(1+n)=stepl
                  strd    =stepl
               endif
               ieos(ij)=ieosdc
            enddo    
         else if(ico2.lt.0.and.ice.eq.0)then
c
c     determine phase state for isothermal air-water flow
c
            call airctr_part(1,ndummy,zone)
         endif
         if(ice.ne.0)then
c
c     determine phase state for low-temperature solid-liquid-gas system
c
c     water state
            call icectr(-1,ndummy)
c     clathrate state
            call icectr(1,ndummy)
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
         if(ico2.gt.0.and.ice.ge.0) then
            call thrmwc(ndummy)
         else if(ico2.eq.0.and.ice.eq.0) then
            call thermw(ndummy)
         else if(ico2.lt.0.and.ice.eq.0)then
            call airctr_part(3,ndummy,zone)
         end if
         if(ice.ne.0)then
            call icectr(3,ndummy)
         endif
c     
c
c     Newton-Raphson corrections
c
      else
         if(ico2.eq.0.and.ice.ge.0) then
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
               if(ps(i).eq.0.0) then
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
         else if(ico2.gt.0.and.ice.eq.0) then
c
c     water and noncondensible
c
c           strd is passed through common
            nr1=nrhs(1)
            nr2=nrhs(2)
            nr3=nrhs(3)
            do i=1,neq
               i1=i+nr1
               i2=i+nr2
               i3=i+nr3
               ieosd=ieos(i)
               if(ps(i).eq.0.0) then
c gaz 10-18-2001
                  t(i)=t(i)-bp(i2)*strd
               elseif(ieosd.eq.1) then
                  phi(i)=phi(i)-bp(i1)*strd
                  t(i)=t(i)-bp(i2)*strd
                  pci(i)=pci(i)-bp(i3)*strd
               elseif(ieosd.eq.2) then
                  phi(i)=phi(i)-bp(i1)*strd
                  s(i)=s(i)-bp(i2)*strd
c  GAZ 5/1/98          
                  t(i)=t(i)-bp(i3)*strd
               elseif(ieosd.eq.3) then
                  phi(i)=phi(i)-bp(i1)*strd
                  t(i)=t(i)-bp(i2)*strd
                  pci(i)=pci(i)-bp(i3)*strd
               endif
c big change gaz 11/26/96
            pci(i)=max(0.0d00,pci(i))
            s(i)=min(1.d00,max(0.0d00,s(i)))
            enddo      
         else if(ico2.lt.0.and.ice.eq.0) then
c     
c     make corrections for isothermal air-water mixture
c     
            call airctr_part(2,ndummy,zone)
         endif
         if(ice.ne.0) then
c     
c     make corrections for solid-liquid-gas mixture
c     
            call icectr(2,ndummy)
         endif
c
      endif
      return
      end

