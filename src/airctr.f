      subroutine airctr(iflg,ndummy)
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
CD1 To manage the isothermal air-water calculations. 
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
CD2 $Log:   /pvcs.config/fehm90/src/airctr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.21   Mon Mar 31 08:28:24 1997   gaz
CD2 new iteration parameters
CD2 
CD2    Rev 1.20   Wed Jun 12 16:44:12 1996   zvd
CD2 Added missing comma to format statement
CD2 
CD2    Rev 1.19   Wed Jun 12 15:20:58 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.18   Mon Jun 03 11:17:42 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.17   Fri May 31 10:33:08 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.16   Wed May 08 13:51:52 1996   hend
CD2 Rearranged and Added Output
CD2 
CD2    Rev 1.15   Fri Apr 26 14:45:10 1996   gaz
CD2 minor changes to phase change criteria
CD2
CD2    Rev 1.14   Wed Feb 14 10:19:00 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.13   Wed Feb 07 10:07:44 1996   gaz
CD2 step length changes and phase change stradegy
CD2 
CD2    Rev 1.12   Mon Jan 29 13:10:44 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.11   Mon Jan 29 10:02:32 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.10   11/29/95 14:02:18   gaz
CD2 format change
CD2 
CD2    Rev 1.9   11/15/95 10:13:32   gaz
CD2 complimentary changes from air_rdof.f
CD2 
CD2    Rev 1.8   08/18/95 09:50:04   llt
CD2 irlp and irlpt were already defined, removed for cray
CD2 
CD2    Rev 1.7   06/02/95 10:16:58   llt
CD2 removed upwinding commons plus gaz changes
CD2 
CD2    Rev 1.6   05/01/95 15:16:14   gaz
CD2  phase change control modified (gaz)
CD2 
CD2    Rev 1.4   03/24/95 13:54:56   gaz
CD2 gaz added phi(I) lt crl(4,1) to phase change criteria
CD2 
CD2    Rev 1.3   03/23/95 18:14:52   gaz
CD2 gaz now determine and save phase state
CD2 also enable ico2 to determine irdof behavior
CD2 
CD2    Rev 1.1   03/18/94 15:40:54   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:12   pvcs
CD2 original version in process of being certified
CD2
c 17-mar-94
c got rid of b array access(done in wrtout.f
c 2/9/95 gaz initialized liquid pressure
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 iflg                    int       I   Parameter used to control
CD3                                           the execution of the
CD3                                           routine
CD3 ndummy                  int       I   Parameter used to obtain
CD3                                           correct node number for
CD3                                           dual porosity calculations
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name           Description
CD3
CD3 file number inpt
CD3                Contains all input data from FEHMN macros in ASCII
CD3                form
CD3 file number iout
CD3                Contains output data
CD3 file number iatty
CD3                Contains output data
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
CD4 ipsx, n0, ico2, qtc, qtotc, amc, inpt, pmin, pmax, tmin, tmax,
CD4 neq, dtot, ps, ieos, iieos, iwelb, s, phi, t, crl, dmpf, rolf, 
CD4 dmef, dq, ieos, s, dil, thx, thy, thz, nrhs, bp, ntty, iout, idualp,
CD4 m, idpdp, nskw, sk, qh, pcp, dte, dife, qtote, qtotei, n, qc, to,
CD4 tini, so, b, strd, irlp, irlpt
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 ico2d        int         Flag denoting the type of two-phase
CD5                              simulation
CD5 tref         real*8      Reference temperature for isothermal air-
CD5                              water calculation
CD5 pref         real*8      Reference pressure for isothermal air-
CD5                              water calculation
CD5 ndum         int         Temporary storage for variable neq
CD5 pssv         real*8      porosity
CD5 ieoss        int         Equation of state parameter
CD5 iieoss       int         Equation of state parameter
CD5 iwelbs       int         Denotes if wellbore solution is enabled
CD5 ssv          real*8      Saturation
CD5 phisv        real*8      Pressure
CD5 dmpfd        real*8      Derivative of mass accumulation term
CD5                              with respect to pressure
CD5 dmefd        real*8      Derivative of mass accumulation term
CD5                              with respect to energy variable
CD5 dqd          real*8      Derivative of mass source term with
CD5                              respect to pressure
CD5 i            int         Do loop parameter
CD5 mid          int         Do loop parameter
CD5 mi           int         Current node number
CD5 nr1          int         Flag for right hand side of equation
CD5 nr2          int         Flag for right hand side of equation
CD5 i1           int         Index for unknown number
CD5 i2           int         Index for unknown number
CD5 ilev         int         Number of sets of nodes
CD5 mlev         int         Number of nodes at which information is
CD5                              written
CD5 il           int         Do loop index over all node levels
CD5 md           int         Current node number being written
CD5 rqd          real*8      Term used in output calculation
CD5 qcd          real*8      Mass balance term
CD5 irlpsv       int         Temporary storage for irlp value
CD5 irlptsv      int         Temporary storage for irlpt value
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 
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
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations 
CD9 2.5.1 Implement time-step mechanism
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
CPS BEGIN airctr
CPS 
CPS IF there is air present
CPS 
CPS   IF this call is for initialization
CPS   
CPS     Read model flag
CPS     IF isothermal air-water simulation is not specified
CPS       Set flag parameter
CPS     ELSE isothermal air-water simulation is not specified
CPS       Set flag parameter
CPS       Read reference temperature and pressure
CPS       Set minimum temperature and pressure values
CPS       Initialize parameter values
CPS       thermw - calculate density, compressibility, and viscosity...
CPS       ... of water
CPS       Set parameters for reference values
CPS       FOR each node
CPS         Set low values for thermal conductivities
CPS       ENDFOR
CPS       
CPS     ENDIF
CPS   
CPS   ELSEIF this call is to determine the phase state
CPS   
CPS     FOR each node
CPS       Determine flag for phase state
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to update the solution
CPS   
CPS     FOR each node
CPS       Compute new pressure, saturation
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to compute air water thermodynamics
CPS   
CPS     thrair - compute isothermal air water thermodynamic parameters
CPS     
CPS   ELSEIF this call is for output
CPS   
CPS     IF air output is being written
CPS     
CPS       IF this is a dual porosity simulation
CPS         Set parameter values for writing
CPS       ELSEIF this is a dpdp simulation
CPS         Set parameter values for writing
CPS       ELSE it is a single porosity simulation
CPS         Set parameter values for writing
CPS       ENDIF
CPS       
CPS       FOR each set of nodes
CPS         Write header denoting which set of nodes
CPS         FOR each node being written
CPS           Determine node number
CPS           Determine sink term written during output
CPS           Write information for this node
CPS         ENDFOR
CPS       ENDFOR
CPS       Write overall mass balance information
CPS       IF output is also going to a second file
CPS         Write overall mass balance information
CPS       ENDIF
CPS       
CPS     ENDIF
CPS     
CPS   ELSEIF this call is for initialization
CPS     
CPS     FOR all nodes
CPS       Set temperatures, source flow rates, and EOS flag
CPS       IF this is a compressed liquid
CPS         Set saturation to unity
CPS       ENDIF
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to zero out enthalpy in equations for...
CPS   ... isotherms problems
CPS   
CPS     FOR each node
CPS       Zero out residual terms for enthalpy equations
CPS     ENDFOR
CPS     
CPS   ENDIF ends all choices for control parameter
CPS   
CPS ENDIF there is air present
CPS 
CPS END airctr
CPS 
C**********************************************************************
c
c subroutine to manage isothermal air calculations
c
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comxi
      use comwt
      use davidi
      use comflow 
      use comsplitts 
      use com_prop_data, only : xnl_ngas, ihenryiso, xnl_max, xnl_ini,
     &   xnl_chng_low,xnl_chng_high, den_h2o, visc_h2o
      use com_nondarcy 
      implicit none

      integer iflg,ndummy,ico2d,ndum,nndum,ieoss,iieoss,iwelbs,i,mid
c gaz 121821 added n0dum for iso water properties using fluid_control_prop.f
      integer n0dum
      integer mi,i1,i2,ilev,mlev,il,md,irlpsv,irlptsv,icapsv,nr1
      integer nr2,irdofsv 
      integer  ii,ij,im,inode,iwm,j,ja,k,kb
c gaz 110819 removed tref, pref (now global)       
      real*8 pssv,ssv,phisv,dmpfd,dmefd,dqd,rqd,qcd
      real*8 inflow_thstime,inen_thstime,denht,deneht
      real*8 dels,delp,dfdum11,dfdum12,dfdum21,dfdum22
      real*8 dfdum11i,dfdum12i,dfdum21i,dfdum22i,detdf
      real*8 fdum01,fdum02,sx1d,phidum,phi_dif,phi_1,phi_2
      real*8 hmax, hmin, hmid
      real*8 cden_correction
      character*80 form_string
      character*80 dum_air
      real*8 pref_1_2,pref_2_1,s_1_2
      real*8 t_low
c gaz 103019 added strd_satneg for under relaxation when neg saturations 
c gaz 081623   
      real*8 strd_satneg, strd_old, phi_unsat_to_sat, p_uzmin
      real*8 strd1, strd2, phi_old
c gaz 042224 
      real*8 xnl_tol, s_xnl_tol
c gaz 081823 ieosd
      integer ieosd,nr_test
      integer i_t_bad 
c gaz debug 112119      
      integer ieos_c
c gaz 071223
      integer ico2_sv
c gaz 092723
      integer isotherup    
      parameter (isotherup = 0)   
      parameter(t_low = 5.,strd_satneg = 0.85, nr_test = 2)
      parameter (p_uzmin = 0.099, phi_unsat_to_sat = 0.001)
c      parameter (pchng = 0.005,schng = 0.005)
c      gaz pchng and schng in comai 103009
c gaz 110819 tref,pref now global      
c      save tref,pref
      if (jswitch.ne.0) strd_iter = strd_rich
      
      k = l +iad 
      ico2d = 0
c     return if no air present 
      if(ico2.lt.0) then
         if(iflg.eq.0) then
c      read data and initialization
            qtc=0.0
            qtotc=0.0
            amc=0.0
c     
c     read in reference pressure and temperature
c
c air is default gas
            itype_air = 1
            itype_meth = 0
            itype_co2 = 0
            itype_h2 = 0
c gaz 030924 ihenryiso set scannin           
c            ihenryiso = 0
            read(inpt,'(a80)')  dum_air    
            read(dum_air,*) tref, pref
c gaz 121323 added read Henry's constant
            do i = 1, 77
             if(dum_air(i:i+3).eq.'meth'.or.
     &        dum_air(i:i+3).eq.'METH') then
              itype_meth = 1
              itype_air = 0
              do j = i+2,80
               if(dum_air(j:j+4).eq.'henry'.or.
     &          dum_air(j:j+4).eq.'HENRY') then
                read(dum_air(j+5:80),*) alpha_meth
                go to 1000
               endif
              enddo
              go to 1000
             else if(dum_air(i:i+2).eq.'co2'.or.
     &        dum_air(i:i+2).eq.'CO2') then
              itype_co2 = 1
              itype_air = 0
              do j = i+2,80
               if(dum_air(j:j+4).eq.'henry'.or.
     &          dum_air(j:j+4).eq.'HENRY') then
                read(dum_air(j+5:80),*) alpha_co2
                go to 1000
               endif
              enddo
              go to 1000
             else if(dum_air(i:i+1).eq.'h2'.or.
     &        dum_air(i:i+1).eq.'H2') then
              itype_h2 = 1
              itype_air = 0
              do j = i+2,80
               if(dum_air(j:j+4).eq.'henry'.or.
     &          dum_air(j:j+4).eq.'HENRY') then
                read(dum_air(j+5:80),*) alpha_h2
                go to 1000
               endif
              enddo
              go to 1000
             else if(dum_air(i:i+2).eq.'air'.or.
     &        dum_air(i:i+2).eq.'AIR') then
              itype_air = 1
c gaz 051224  "i+2,80" to "i+2,76"  works for Release version          
              do j = i+2,76
               if(dum_air(j:j+4).eq.'henry'.or.
     &          dum_air(j:j+4).eq.'HENRY') then
                read(dum_air(j+5:80),*) alpha_air
                go to 1000
               endif
              enddo
              go to 1000
             endif
            enddo
1000    continue
c    density(kg/m**3) and viscosity(Pa-sec) are referenced for NIST 
c     at 1.0 Mpa and 0 C for the perfect gas equation
        if(itype_meth.eq.1) then
c methane props at 0.1 Mpa and 20 C
         roc0 = 0.7081d0
         visc_gas = 1.037e-5
         if(iout.ne.0) write(iout,*) 
     #   '>>>> gas phase fluid changed to methane (perfect gas) <<<'   
         if(iptty.ne.0) write(iptty,*) 
     #   '>>>> gas phase fluid changed to methane (perfect gas) <<<'   
        else if(itype_co2.eq.1) then
c co2 props at 0.1 Mpa and 20 C
         roc0 = 1.951d0
         visc_gas = 1.371e-05
         if(iout.ne.0) write(iout,*) 
     #   '>>>> gas phase fluid changed to co2 (perfect gas) <<<'   
         if(iptty.ne.0) write(iptty,*) 
     #   '>>>> gas phase fluid changed to co2 (perfect gas) <<<'   
        else if(itype_h2.eq.1) then
c h2 props at 0.1 Mpa and 20 C
c values from QM Bui 2018
         roc0 = 0.083D0
         visc_gas = 0.9e-05
         if(iout.ne.0) write(iout,*) 
     #   '>>>> gas phase fluid changed to H2 (perfect gas) <<<'   
         if(iptty.ne.0) write(iptty,*) 
     #   '>>>> gas phase fluid changed to H2 (perfect gas) <<<'   
        else if(itype_air.eq.1) then
c air props at 0.1 Mpa and 0 C
         roc0 = 1.292864d0
         visc_gas = 1.758e-05
        else
c default is air
         roc0 = 1.292864d0
         visc_gas = 1.758e-05
        endif
c
c gaz 121923 initializing cnlf need to program input   <<<<<<<<<<< 
c
c gaz 031024 initial mass fracion done in pres
      go to 450        
       read(inpt,'(a80)')  dum_air  
        if(dum_air(1:8).eq.'massfrac') then
         read(dum_air(9:80),*) xnl_ini
         xnl_ini = max(xnl_ini,1.d-20)
        else
         backspace inpt
         xnl_ini = 1.d-20
        endif
c gaz 021324
        if(ihenryiso.ne.0) then
         do i =1,n0                   
          cnlf(i) = xnl_ini
          dclf(i) = 0.0d0
          dclef(i) = 0.0d0  
         enddo
        endif
450   continue
c     
c     set max and min values of p and t
c     
            pmin(1)=-1.e15
            pmax(1)=1.e15
            tmin(1)=-1.e05
            tmax(1)=1.e05

c     Read in NAPL properties for SZ NAPL case

            if(ico2.eq.-3) then
               read(inpt,*) dennapl, viscnapl
            end if
         elseif(iflg.eq.-1) then

c     
c     calculate density,compressibility,and viscosity of water
c     
            ndum=neq
            nndum = n
c gaz 121821  added n0dum for fluid_control_prop.f     
            n0dum = n0
            n0 = 1
            irdofsv=irdof
            irdof=0
            neq=1
            n = 1
            dtot=1.0
            pssv=ps(1)
            ps(1)=1.0
            ieoss=ieos(1)
            iwelbs=iwelb
            iwelb=0
            ssv=s(1)
            irlpsv=irlp(1)
            irlptsv=irlpt(1)
            icapsv = icapp
            icapp = 0
            s(1)=1.0
c gaz  050320  keep alternate eos          
c            iieoss=iieos(1)
c            iieos(1)=1
            ieos(1)=1
            phisv=phi(1)
            phi(1)=pref
            t(1)=tref
            irlp(1)=1
            irlpt(1)=1

            call thermw(0)

            irdof=irdofsv
            neq=ndum
            n = nndum
            n0 = n0dum
            ieos(1)=ieoss
c gaz 050320 keep altermate eos            
c            iieos(1)=iieoss
            phi(1)=phisv
            ps(1)=pssv
            s(1)=ssv
            iwelb=iwelbs
            irlp(1)=irlpsv
            irlpt(1)=irlptsv
            icapp = icapsv
c     reference density is in crl(1,1)
c     reference liquid viscosity is in crl(2,1)
c     reference compressibility is in crl(3,1)
c     reference pressure in crl(4,1)
c     reference air viscosity is in crl(5,1)
c     subtract initial density calculation (added back in thrair.f)
            if(cden) then
             crl(1,1)=rolf(1)  - cden_correction(1)
            else
             crl(1,1)=rolf(1)
            endif
            crl(2,1)=1.0/(dil(1)/rolf(1))
            crl(3,1)=dmpf(1)/rolf(1)
c gaz 110819 pref, tref (global) read in scanin  (and above)            
c            crl(4,1)=pref
            crl(5,1)=visc_gas
c            crl(6,1)=tref
            if(ico2d.eq.-1) then
               crl(7,1)=1.0
            else
               crl(7,1)=0.0
            endif
c gaz 110424 check for non-darcy flow (nd_flow)
      if(nd_flow) then
       if(allocated(den_h2o)) then
        deallocate(den_h2o)
        allocate(den_h2o(n0,6))
        deallocate(visc_h2o)
        allocate(visc_h2o(n0,6))        
       endif
      endif

c     initialize 2-phase regions             
c     
            do i=1,n0
               ieos(i)=2
            enddo
         elseif(iflg.eq.-2) then
c     
c     for rich eq make "a" smaller
c     eliminate a_vxy
c     
            if(jswitch.ne.0) then
c               if (idpdp .eq. 0 .and. allocated(a_vxy)) 
c     &              deallocate(a_vxy)
               if (allocated(a_vxy)) deallocate(a_vxy)
               ndum = nelm(neq+1) - (neq+1)
               if(allocated(a)) deallocate(a)
               if (idpdp .eq. 0) then
                  allocate(a(2*ndum))
               else
                  allocate(a(8*ndum))
               end if
            end if
         elseif(iflg.eq.1) then
c     
c     determine phase state
c     
            if(ifree.ne.0) then
c gaz 090819 mass conservative phase change
c             call phase_change_mass_conv(0,0)
               do  mid=1,neq
                  mi=mid+ndummy
                  ieos(mi) = 2
               enddo
               continue
            else if(irdof.ne.13) then   
c gaz 103019 might need to set strd = 1 for every iad 
c               if (iad.eq.0) strd = 1.0 
               strd = 1
c gaz 110819 pref, tref (global) read in scanin                 
c               pref = crl(4,1) 
c gaz 082719               
               pref_1_2 =  pref - pchng
               pref_2_1 =  pref + pchng  
               s_1_2 = 1.0 - schng
c gaz debug 112119 added ieos_c
	       if (jswitch.ne.0) then
c determine phase state Rich Eq                 
                  ieos_c = 0
                  do  mid=1,neq
                     mi=mid+ndummy
c     
c gaz 103019 testing pref here  (pref_1_2 to pref) if(phi(mi).lt.pref_1_2.and.ieos(mi) no change
c   s_1_2 to 1.               best =       1. - 0.000001
c gaz 103019 set strd_old to strd                       
                     strd_old = strd
                     if(phi(mi).lt.pref.and.ieos(mi)
     &                 .eq.1.and.days.ge.time_ieos(mi)) then
                        strd = strd_iter
                        s(mi) = s_1_2
                        phi(mi) = pref
                        ieos(mi)=2 
                        time_ieos(mi) = days + time_ch
                        ieos_c = ieos_c + 1
                     else if(s(mi).lt.1.0-tol_phase.and.ieos(mi) 
     &                   .eq.1.and.days.ge.time_ieos(mi)) then 
                        strd = strd_iter
                        s(mi) = s_1_2
                        ieos(mi)=2 
                        time_ieos(mi) = days + time_ch
                        ieos_c = ieos_c + 1
                     else if(s(mi).gt.tol_phase.and.ieos(mi)
     &                .eq.3.and.days.ge.time_ieos(mi)) then 
c   changed eq.2 to eq.3 gaz 103009     
                        strd = strd_iter
                        ieos(mi)=2   
                        ieos_c = ieos_c + 1                  
                     else if(s(mi).ge.1.0+tol_phase.and.ieos(mi)
     &                   .eq.2.and.days.ge.time_ieos(mi)) then     
                        strd = strd_iter
c gaz debug 110119       phi_unsat_to_sat (0.001) instead of     phi(mi) = pref_2_1
                        phi(mi) = pref + phi_unsat_to_sat
                        s(mi) = 1.
                        ieos(mi)=1
                        time_ieos(mi) = days + time_ch
                        ieos_c = ieos_c + 1
                     else if(s(mi).le.-tol_phase.and.ieos(mi)
     &                   .eq.2.and.days.ge.time_ieos(mi)) then 
c gaz 071919 uncommented  s(mi)=0.0     
c gaz 072019 commented out again
c gaz 103019 and again   s(mi)=0.0
c                        s(mi)=0.0
                        strd = strd_satneg
c     ieos(mi)=3
                     endif
c big deal?    
c gaz 103019 does not look like min strd is used
                    strd = min(strd_old,strd)
                  enddo
c gaz debug 112119
c                write(ierr,878) 'rich',l,iad,ieos_c,strd,fdum
878             format(a4,' l ',i5,' iad ',i5,' ieos_c ',i3,
     &             ' strd ',f13.3,' fdum ', f15.8)           
               else
c     
c     determine phase state uzsz
c     
                  if(irdof.ne.13) then
c     pnx used to be set in the following loops -- currently removed
c gaz 081723 try skipping section
                    if(nr_test.eq.0) then
c  do none
                     go to 778
                    else if(nr_test.eq.1) then
c  do 1st part only 
                    else if(nr_test.eq.2) then
c  do both
                    endif
                     ieos_c = 0
c gaz 100723  use default schng = 1.e-3 (sub data.f)
c                     schng = 1.d-2
                     do  mid=1,neq
                        mi=mid+ndummy
c     
                        strd_old = strd
                        if(s(mi).lt.1.0-schng.and.ieos(mi).eq.1) then
                           s(mi) = 1.0-schng                           
                           ieos(mi)=2
                           ieos_c = ieos_c +1  
                        else if(s(mi).gt.0.0.and.ieos(mi).eq.3) then
                           strd = strd_iter  
                           ieos(mi)=2
                           ieos_c = ieos_c +1 
                        else if(s(mi).gt.1.0.and.ieos(mi).eq.2) then
                          if(ihenryiso.eq.0) then
                           ieos(mi)=1
c gaz 120919                            
                           s(mi) = 1.0d0
                           phi(mi) = phi(mi) + pchng
c gaz 010623
                           pcp(mi) = 0.0d0
                           strd = min(strd_iter,strd)
                           ieos_c = ieos_c +1
                          endif
c  gaz 081923 le to lt  
                        else if(s(mi).lt.0.0d0.and.ieos(mi).eq.2) then
                           s(mi)=0.0d0
                           strd = strd_iter
                           ieos(mi)=3
                           ieos_c = ieos_c +1  
                        endif
                       if(ihenryiso.ne.0) then
                        xnl_tol =1.d-09
                        s_xnl_tol = 1.d-6
                        if(s(mi).ge.1.d0.and.ieos(mi).eq.2) then
c gaz 042124 
                         call solubility_isothermal(1,mi)
                         xnl_max = xnl_ngas(mi,1)
                         if(abs(cnlf(mi)-xnl_max).gt.xnl_tol) then 
                          cnlf(mi)= xnl_max 
                          s(mi) = 1.d0 
                          ieos(mi) = 1
                          strd = strd_mass 
                          strd = 0.995
                         else
                          s(mi) = 1.0d0-s_xnl_tol
                          strd = 0.995d0
                         endif                    
                        else if(s(mi).eq.1.d0.and.ieos(mi).eq.1) then
                         call solubility_isothermal(1,mi)
                          xnl_max = xnl_ngas(mi,1)
c gaz  042824 test
c                         if(cnlf(mi).gt.xnl_max+xnl_tol) then
                          if(abs(cnlf(mi)-xnl_max).gt.xnl_tol.and.
     &                       cnlf(mi).gt.xnl_max) then
                            cnlf(mi)= xnl_max
                            ieos(mi)=2
c gaz 041724 calc s(mi) in  phase_change_mass_conv
                            s(mi) = 1.d0- s_xnl_tol
c
                            strd = 0.995d0
                          else
                            ieos(mi)=1
                          endif
                        else if(s(mi).lt.1.0d0.and.ieos(mi).eq.2) then
                          call solubility_isothermal(1,mi) 
                          xnl_max = xnl_ngas(mi,1)
                          if(cnlf(mi).lt.xnl_max-xnl_tol) then
c                         call phase_change_mass_conv(3,mi,mi)
c                            strd = strd_mass
c                            strd = 0.95
c                            cnlf(mi)= xnl_max
c                            ieos(mi) = 2
                          endif                         
                        endif
                       endif
                             strd = min(strd_old,strd)
                     enddo
c gaz 110319 could use some cleanup  
c gaz 121923 good place to check phase gas solubility
777               continue
                  if(nr_test.eq.1) goto 778                   
                     do  mid=1,neq
                        mi=mid+ndummy
                        strd_old = strd
c gaz 120919 debug  
c                        if(phi(mi).lt.p_uzmin) then
                        if(phi(mi).lt.p_uzmin.and.ieos(mi).eq.1) then
                         strd = strd_iter
c gaz 112419 
                         phi(mi) = p_uzmin
c                         s(mi) =0.999
                         s(mi) = 1.0-schng
                         ieos(mi) = 2
c  gaz 120919 
c                           if (so(mi).ge.1.0) then
                         if (so(mi).ge.100.0) then
                              s(mi)= 0.999
                              phi(mi) = p_uzmin
                              strd = strd_iter
                           endif
                        endif
c     
                        if(s(mi).lt.-schng) then                            
                          strd = strd_iter
c gaz 120919 13:50                          
                          ieos(mi) = 3
                          s(mi) = 0.0
c     strd = strd_iter
                        endif
c                       
                        if(s(mi).gt.1.0+schng) then
                         strd = strd_iter
                         ieos(mi) = 1
                         s(mi) = 1.
                        endif
c     
                        strd = min(strd_old,strd)
                     enddo
c                   write(ierr,878) 'uzsz',l,iad,ieos_c,strd,fdum
c gaz 081723 (end of nr_test eq 0)
778               continue
                  endif
               endif
            endif
         elseif(iflg.eq.2) then
c gaz 121923 enable variable switching
           if(ihenryiso.ne.0) then
            call phase_change_mass_conv(2,1,n0)
            return
           endif
            if(abs(irdof).ne.14) then
c     
c     update solution
c     
               nr1=nrhs(1)
               nr2=nrhs(2)
c gaz 081723 testing 
                  strd  = 1.0d0
                  strd1 = 1.0d0
                  strd2 = 1.d0-strd1
c     strd is passed through common
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2  
                  phi_old = phi(i)
                  phi(i)=phi_old-bp(i1)*strd
                  if (irdof .ne. 13 ) 
     &                 s(i)=s(i)-bp(i2)*strd
c gaz 081623  make first correction and phase change here   
                if(isotherup.ne.0) then
                  ieosd = ieos(i)
                 if(nr_test.lt.1) then
                  if(s(i).lt.0.0d0.and.ieosd.eq.2) then
c less pressure change 2 to 3
                   strd1 = 0.95d0
                   strd2 = 1.d0-strd1
                   ieos(i) = 3                                            
                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = 0.0d0
                  else if(s(i).gt.1.0d0.and.ieosd.eq.2) then
c more pressure change 2 to 1
                   strd1 = 1.0d0
                   strd2 = 1.d0-strd1
                   ieos(i) = 1
                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = 1.0d0
                  else if(s(i).le.1.0d0.and.ieosd.eq.1) then
c less pressure change 1 to 2
                   strd1 = 0.95d0
                   strd2 = 1.d0-strd1
                   ieos(i) = 2
                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = max(0.0d0,s(i))
                  else if(s(i).ge.0.0d0.and.ieosd.eq.3) then
c more pressure change 3 to 2
                   strd1 = 0.95d0
                   strd2 = 1.d0-strd1
                   ieos(i) = 2
                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = min(1.0d0,s(i))
                  else if(s(i).lt.0.0d0.and.ieosd.eq.3) then
c check for min saturation
c                   strd1 = 1.0d0
c                   strd2 = 1.d0-strd1
c                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = min(0.0d0,s(i))
                  else if(s(i).gt.1.0d0.and.ieosd.eq.1) then
c check for max saturation
c                   strd1 = 1.00d0
c                   strd2 = 1.d0-strd1
c                   phi(i)=(phi_old*strd2+phi(i)*strd1)
                   s(i) = max(1.0d0,s(i))
                  endif  
                endif
c gaz 092723 end isoco
                endif                                                       
               enddo
c     
c     call cascade redistribution if requested
c     
               if(iflux_ts.ne.0) then
                  call cascade_sat(1)
               endif
c     
            else
c     
c     update solution(with exact mass balance)
c     
               nr1=nrhs(1)
               nr2=nrhs(2)
c     strd is passed through common
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  phi(i)=phi(i)-bp(i1)*strd
               enddo
c     
c     call wellrate(0,0) 
c     call thrair(0)
               call interblock_iso(0) 
               
            endif
         elseif(iflg.eq.3) then
            if(ico2.eq.-3) then
c     call isothermal SZ NAPL thermodynamics
               call thrsznapl(ndummy)
            else
c     call isothermal air water thermodynamics
               if(ifree.ne.0) then
                  call wtsictr(2)
                  call thrair(ndummy)
c     call wtsictr(4)
               else
                  call thrair(ndummy)
               endif
            end if
         elseif(iflg.eq.4) then
c     call equation generation and load a array
c     new:zero a_axy (anisotropy requires accumulation)
            a_axy=0.0d00
            if(jswitch.eq.1) then
               call gensl2_switch
            else 
               call gensl2
            endif
         elseif(iflg.eq.-4) then
c     solve explicit equations
            call gensl2_explicit       
c     
         elseif(iflg.eq.5 .and. irdof .ne. 13) then
            if(ntty.eq.2) then
c     
c     output for air
c    
c gaz 111823 write out gas name
            if(itype_air.ne.0) then
             dum_air(1:9) = 'Air      '
            else if(itype_meth.ne.0) then
             dum_air(1:9) = 'Methane  '
            else if(itype_co2.ne.0) then
             dum_air(1:9) = 'CO2      '
            else if(itype_h2.ne.0) then
             dum_air(1:11) = 'Hydrogen '
            endif

 
               if (iout .ne. 0) write(iout,803)  dum_air(1:9)
               if (iatty .ne. 0) write(iatty,803) dum_air(1:9)
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
c     
               if(idualp.ne.0) then
                  ilev=3
                  mlev=m/3
               else if(idpdp.ne.0) then
                  ilev=2
                  mlev=m/2
               else
                  ilev=1
                  mlev=m
               endif
               
               do il=1,ilev
                  if(il.ne.1) then
                     if (iout .ne. 0) write(iout,702) il
                     if (iatty .gt. 0) write(iatty,702) il
 702                 format(2x,'Matrix Level = ',i1)
                  endif
                  do i=1,mlev
                     md=  nskw(i+(il-1)*mlev)
                     rqd= sk(md)
                     qcd=qh(md)
                     if(ihead.ne.0.and.ifree.ne.0) then
c     wtsi solution
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref                          
                        if(s(md).le.0.0) then
                           phidum = pref - phi_inc
                        else if(s(md).lt.1.0) then
                           phi_1 = head12(md,1)
                           phi_2 = head12(md,2)
                           phi_dif = phi_2-phi_1
                           phidum = pref + s(md)*(phi_dif)
                        else
                           phi_1 = head12(md,1)
                           phi_2 = head12(md,2)
                           phidum = phi(md) - phi_inc + 
     &                          0.5*(phi_2-phi_1)
                        endif
                        if (iout .ne. 0) write(iout,804) 
     &                       md,phidum,0.0,phidum,qcd
                        if (iatty .ne. 0) 
     &                       write(iout,804) md,phidum,0.0,phidum,qcd
                     else if (ihead.ne.0 .or. 
     &                       (irdof .eq. 13 .and. abs(ifree) .ne. 1)) 
     &                       then
c     head solution
                        if (iout .ne. 0) write(iout,804) 
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                        if (iatty .ne. 0) write(iatty,804)
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                     else
c     two-phase solution
c gaz 121923 output total mass fraction, liquid mass fraction
                        call phase_change_mass_conv(1,md,md)
                        if (iout .ne. 0) write(iout,804) 
     &                     md,phi(md),pcp(md),
     &                     phi(md)-pcp(md),qcd,1.-s(md),
     &                     frac_gas_iso(md),cnlf(md)
                        if (iatty .ne. 0)
     &                     write(iatty,804) md,phi(md),pcp(md),
     &                     phi(md)-pcp(md),qcd,1.-s(md),
     &                     frac_gas_iso(md),cnlf(md)
                     endif
                  enddo
               enddo
 803           format(/,20x,'Nodal Information (Vapor)',/, 9x, 
     &              a9,26x, 'source/sink', /, 3x,'Node',1x,
     &              '  P (MPa)',3x,'P Cap (MPa)',1x,'P Liq (MPa)'
     &              ,1x,'Vapor (kg/s)',5x,'S gas',8x,
     &              'Tot (mfrac)',4x,'Liq (mfrac)')
 804           format(i7,1p,1x,g11.4,1x,g11.4,1x,g11.4,1x,g11.4,3x,
     &                g13.5,3x,g12.5,3x,g12.5)
c     calculate global mass and energy flows
               if (iout .ne. 0) then
                  write(iout,703) 
                  write(iout,704) qtotei
                  write(iout,705) qtote
                  write(iout,706) qte
                  write(iout,707) dife
               end if
c     calculate global mass and energy flows
               if(iatty.ne.0) then
                  write(iatty,703) 
                  write(iatty,704) qtotei
                  write(iatty,705) qtote
                  write(iatty,706) qte
                  write(iatty,707) dife
               endif
 703           format(/,20x,'Global Mass Balances (Vapor)')
 704           format(1x,'Vapor discharge this time step: ',e14.6,' kg')
 705           format(1x,'Total vapor discharge: ',9x,e14.6,' kg')
 706           format(/,1x,'Net kg vapor discharge (total out-total ',
     &              'in): ',e14.6)
 707           format(1x,'Conservation Error: ',25x,e14.6)
            endif
         elseif(iflg.eq.6) then
c     
c     store sk in qc
c     initialize t,to,tini,iieos
c  
c gaz 110919 tref now global   
c            tref=crl(6,1)
            do i=1,n
               qc(i)=sk(i)
               if(to(i).eq.0.0d00) then
                  t(i)=tref
                  to(i)=tref
                  tini(i)=tref
               endif
               iieos(i)=1
               if(ieos(i).eq.1) then
c gaz 042224 
c                  ieos(i)=2
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     s(i)=1.0d0
                     so(i)=1.0d0
                  end if
               endif
               if(abs(irdof).eq.14) then
                  denj(i)=1.0d0-s(i)
               endif
            enddo
c     
c     check also for free surface calcs
c     
            if(ifree.ne.0) then
               call wtsictr(1)	
c               call wtsictr(12)
c     call wtsictr(9)
            endif
         elseif(iflg.eq.7) then
c     
c     this only applies if we have a moving water table, etc.
c     
            call headctr(3,0,0.0,0.0)

         elseif(iflg.eq.8) then
c     
c     convert from pressure to head
c     
            call headctr(2,0,0.0,0.0)
         elseif(iflg.eq.9) then
c     
c     convert from head to pressure boundary value
c     
            call headctr(6,0,0.0,0.0)

         elseif(iflg.eq.10) then
c     
c     convert from head to pressure boundary and initial 
c     conditions (based on pair = pref at max height)
c     correct for negative pressures
c     
            call head_2phase(0)           
            
         elseif(iflg.eq.11) then
c     
c     calculate the gridblock length in the gravity direction
c     
            if(.not.allocated(dzrg))then
               allocate(dzrg(neq))
            else
               deallocate(dzrg)
               allocate(dzrg(neq))               
            endif 
            do i = 1,neq_primary
               i1=nelm(i)+1
               i2=nelm(i+1)
               hmid=cord(i,igrav)
               hmin=0.
               hmax=0.
               do ii =i1,i2
                  kb=nelm(ii)
                  hmax=max(cord(kb,igrav)-hmid,hmax)
                  hmin=min(cord(kb,igrav)-hmid,hmin)
               enddo
c     distinguish between block and edge centered
               if(ivf.eq.-1) then
                  dzrg(i) = max(hmax,abs(hmin))
               else
                  dzrg(i) = abs(hmax-hmin)/2.	            
               endif      
            enddo
            
         elseif(iflg.eq.12) then
c     write wt outpt
c     do i=1,n_wt_cols
            iwm=0
            do ij=1,m
               i=wcol(nskw(ij)) 
               do k=1,iwm
                  if(i.eq.col_out(k)) goto 566
               end do
               iwm=iwm+1
               col_out(iwm)=i
               do im=n_col(i),1,-1
                  inode=col(i,im)
                  if(s(inode).lt.1.or.im.eq.1) then
                     wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                    cord(inode,igrav)
                     if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                        inode=col(i,im+1)
                        wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                       cord(inode,igrav)
                     endif
                     if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                        if (iptty .ne. 0) write(iptty, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                        if (iout .ne. 0) write(iout, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                     endif
                     goto 4009
                  endif
               end do
 4009          continue
               if (form_flag .eq. 2) then
                  write(ishiswt,4004) days, cord(inode,1),
     &                 cord(inode,2), cord(inode,3), wt_elev,
     &                 ps(inode), inode, nskw(ij)
               else
                  write(ishiswt,4005) days, cord(inode,1),
     &                 cord(inode,2), cord(inode,3), wt_elev,
     &                 ps(inode), inode, nskw(ij)
               end if
 566        end do
 4004       format(1x,6(g16.9,', '),i8,', ',i8)
 4005       format(1x,6(g16.9,1x),2(i8,1x))
 4006       format(1x,'all nodes are dry ', i8, 1x, 2(g16.9, 1x))
         else if (iflg.eq.13) then
c check for bad temperatures
           i_t_bad = 0
           do i = 1,n
            if(to(i).lt.t_low) then
             i_t_bad = i_t_bad + 1
             write(ierr,*) 'node ', i, ' T changed from ',to(i),' to ',
     &                     tref
             to(i) = tref
             t(i) = tref
            endif
           enddo
           if(i_t_bad.gt.0) then
            if(iptty.ne.0) write(iptty,*)'************************'
            if(iout.ne.0) write(iout,*)'************************'
            if(iptty.ne.0) write(iptty,*) i_t_bad,' temps are too low ',
     &                    '(set to Tref) - see error file'
            if(iout.ne.0) write(iout,*) i_t_bad,' temps are too low',
     &                    '(set to Tref) - see error file'
            if(iptty.ne.0) write(iptty,*)'************************'
            if(iout.ne.0) write(iout,*)'************************'
           endif
         else if (iflg.eq.14) then
c     write wt output for contours
            if (altc(1:4) .eq. 'avsx') then
               write (form_string, 4015) ' : ', ' : '
               write (isconwt, 4019)
            else if (altc(1:3) .eq. 'sur') then
               write (form_string, 4015) ', ', ', '
               write (isconwt, 4020)
            else if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'tec') 
     &              then
               write (form_string, 4015) ' ', ' '
               if (altc(1:3) .eq. 'tec') then
                  write (isconwt, 4018) days
               else
                  write (isconwt, '("04 1 1 1 1")')
                  write (isconwt, '(a)') 'X coordinate (m), (m)'
                  write (isconwt, '(a)') 'Y coordinate (m), (m)'
                  write (isconwt, '(a)') 'Z coordinate (m), (m)'
                  write (isconwt, '(a)') 
     &                 'Water table elevation (m), (m)'
               end if
            end if
            do i=1,n_wt_cols
               do im=n_col(i),1,-1
                  inode=col(i,im)
                  if(s(inode).lt.1.or.im.eq.1) then
                     wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                    cord(inode,igrav)
                     if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                        inode=col(i,im+1)
                        wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                       cord(inode,igrav)
                     endif
                     if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                        if (iptty .ne. 0) write(iptty, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                        if (iout .ne. 0) write(iout, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                     endif
                     goto 4010
                  endif
               end do
 4010          continue
               write(isconwt,form_string) cord(inode,1), cord(inode,2),
     &              cord(inode,3), izonef(inode), wt_elev, 0.0
            end do

 4015       format("(1x, 3(g16.9, '", a, "'), i4, 2('", a, "', g16.9))")
 4020       format(1x, 'X (m), Y (m), Z (m), Zone, WT elev (m), ', 
     &           'WT elev2 (m)')
 4019       format(1x, 'X (m) : Y (m) : Z (m) : Zone : WT elev (m) : ',
     &           'WT elev2 (m)')           
 4018       format('variables = "X (m)" "Y (m)" "Z (m)" " Zone" ', 
     &           '"WT elev (m) "', '"WT elev2 (m)"', / 
     &           'zone t = "Simulation time ', g16.9, ' days"') 
         endif
      endif       

      return
      end
      subroutine phase_timcrl(days,day,dtot,dtotdm)
c     
c     adjusts timestep when phase change occurs
c     
      implicit none
      real*8 days,dtot,day,dtotdm
c     
      return
      if(dtot.gt.dtotdm) then
         days=days-day
         dtot=dtotdm  
         day = dtot/86400.0d00
         days=days+day
      endif
      return
      end
