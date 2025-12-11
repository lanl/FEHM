          real*8 function psatl(tl,pcaps,dpcaps,dpsatt,dpsats,
     &                          isatf,salt_con)
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
CD1 To calculate the saturation temperature or pressure.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-30-93     G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/psatl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Feb 05 11:21:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:10:12   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:34   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 tl           real*8  I       Temperature (or pressure if
CD3                                  saturation temperature is computed)
CD3 pcaps        real*8  O       Capillary pressure
CD3 dpcaps       real*8  O       Derivative of cap. pressure with
CD3                                  saturation
CD3 dpsatt       real*8  O       Derivative of saturation pressure
CD3                                  with temperature
CD3 dpsats       real*8  O       Derivative of saturation pressure
CD3                                  with saturation
CD3 isatf        int     I       Control parameter to decide if
CD3                                 saturation temperature or pressure
CD3                                 is computed
CD3                                  
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 File with number iout  O    File used to write warning and error
CD3                             messages (screen output if specified
CD3                             as 6)
CD3 File with number iptty O    File used to write warning and error
CD3                             messages (screen output if specified
CD3                             as 6)
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 ice, ipsat, psta1, psta2, psta3, psta4, psa0, pstb1, pstb2, pstb3,
CD4 pstb4, psb0, ivapl, tspa1,tspa2,tspa3,tspa4, tsa0, tspb1,tspb2,
CD4 tspb3,tspb4, ico2, iout, iptty
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 vaporl          N/A      Computes vapor pressure lowing terms
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
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 x            real*8      Temperature in pressure calculation
CD5 x2           real*8      Square of x
CD5 x3           real*8      Cube of x
CD5 x4           real*8      x to fourth power
CD5 pfun         real*8      Intermediate term in calculation
CD5 pfunn        real*8      Intermediate term in calculation
CD5 pfund        real*8      Intermediate term in calculation
CD5 dpst         real*8      Derivative of saturation pressure with
CD5                              temperature
CD5 dptsn        real*8      Intermediate term in calculation
CD5 dpstd        real*8      Intermediate term in calculation
CD5 psatl0       real*8      Vapor pressure before vapor pressure
CD5                              lowering correction
CD5 delp         real*8      Vapor pressure lowering
CD5 ddelt        real*8      Derivative of vapor pressure with
CD5                              temperature
CD5 ddels        real*8      Derivative of vapor pressure with
CD5                              saturation
CD5 tfun         real*8      Intermediate term in calculation
CD5 tfunn        real*8      Intermediate term in calculation
CD5 tfund        real*8      Intermediate term in calculation
CD5 maxitp       int         Maximum number of iterations allowed
CD5 k            int         Current iteration number in do loop
CD5 pfun0        real*8      Intermediate term in calculation
CD5 resid        real*8      Residual in iterative calculation
CD5 drlp         real*8      Intermediate term in calculation
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
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
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
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN psatl
CPS 
CPS Initialize parameters
CPS 
CPS IF this is a ice solution
CPS   ERROREXIT
CPS ENDIF
CPS 
CPS IF saturation pressure and derivative are being calculated
CPS 
CPS   IF temperature is low
CPS     Set saturation pressure to limiting value
CPS   ELSEIF temperature is high
CPS     Set saturation pressure to limiting value
CPS   ELSEIF saturation pressure is not automatically set to 0
CPS     Compute saturation pressure and derivatives before vapor...
CPS     ... pressure lowering
CPS     
CPS     IF vapor pressure lowing is accounted for
CPS       vaporl - compute vapor pressure lowering contribution
CPS       Compute saturation pressure and derivatives
CPS     ENDIF
CPS     
CPS   ENDIF
CPS 
CPS ELSE saturation temperature is being calculated
CPS 
CPS   IF pressure is low
CPS     Set temperature to limiting value
CPS   ELSEIF pressure is high
CPS     Set temperature to limiting value
CPS   ELSE we calculate temperature
CPS     Compute temperature as though no air is present
CPS     IF air phase is present
CPS     
CPS       FOR all iterations
CPS         Compute temperature without vapor pressure lowering...
CPS         ... contribution
CPS         
CPS         IF vapor pressure lowering is included
CPS           Compute vapor pressure lowering contribution
CPS         ENDIF
CPS         Compute error in iterative solution for temperature
CPS       EXITIF convergence is achieved
CPS       ENDFOR
CPS       IF convergence was not achieved
CPS         Write error message
CPS       ENDIF
CPS     ENDIF
CPS     Set saturation temperature
CPS   ENDIF
CPS     
CPS   ENDIF
CPS ENDIF
CPS 
CPS END psatl
CPS
C**********************************************************************
c version FEHM5.1J changes
c 2-apr-92
c broke out subroutines to work on vapor pressure lowering
c initial guess for temp may be important
c 07-apr-92
c changed tl to tl+273 in sub vaporl
c 26-mar-93
c put in thermc for latest ther_gz.f module
c made sane changes to hflx term as in ther_gz.f
c 30-mar-93
c zero out flow terms and derivatives
c 31-mar-93
c zreo out xnl ,roc as in ther_gz if necessary
c 7-may-93 llt
c changed equivalences to pointers
C***********************************************************************

      use comii
      use comdti
      use comai
      use property_interpolate_1
      implicit none

      integer isatf,k,maxitp,ifail
      real*8 tl,pcaps,dpcaps,dpsatt,dpsats
      real*8 x,x2,x3,x4,pfun,pfunn,pfund,dpst,dptsn,dpstd,psatl0,delp
      real*8 ddelt,ddels,tfun,tfunn,tfund,pfun0,resid,drlp
      real*8 salt_con,pv_sc,dsct,dscc
      real*8 dtps,dtpsn,dtpsd
      psatl=0.0d0
      dpsatt=0.0d0
      if(ice.ne.0) then
         goto 9000
      end if
c ev3 is the reference value for vapor density,initialized(1.) in main.s
c sat pressure as function of sat temp
      if(isatf.le.0.and.ipsat.eq.0) then
c check for limiting values
        if(iwater_table.ne.1) then
c gaz debug 112717        
c         if(tl.lt.10.) then
c            psatl=0.00123
         if(tl.lt.5.0) then
            psatl=0.000752
         else if(tl.gt.340.0) then
            psatl=14.5941
         else
            x=tl
            x2=x*x
            x3=x2*x
            x4=x3*x
            pfunn=psa0+psta1*x+psta2*x2+psta3*x3+psta4*x4
            pfund=psb0+pstb1*x+pstb2*x2+pstb3*x3+pstb4*x4
            pfun=pfunn/pfund
            dptsn=((psta1+2.*psta2*x+3.*psta3*x2+4.*psta4*x3)*pfund)-
     &           (pfunn*(pstb1+2.*pstb2*x+3.*pstb3*x2+4.*pstb4*x3))
            dpstd=pfund**2
            dpst=dptsn/dpstd
            psatl=pfun
            dpsatt=dpst
            dpsats=0.0
         endif
        else
c gaz 102420 make sure sat p and t defined above h2o crit point 
          if(tl.lt.tcrit_h2o) then
           call get_h2o_sat_pressure(ifail,tl,psatl,dpst) 
           dpsatt=dpst
           dpsats=0.0 
          else 
           psatl = pcrit_h2o
            dpsatt=0.0
            dpsats=0.0
          endif
c          call get_h2o_sat_pressure(ifail,tl,psatl,dpst)
        endif
c
c get vapor pressure lowering (salt concentration)
c 
             if(isalt.ne.0.and.ivaprsalt.gt.1) then 
              call vaporl_salt(tl,salt_con,pv_sc,dsct,dscc)
c              psatl = psatl + pv_sc
c              dpsatt= dpsatt + dsct 
c   gaz debug 060316 (embedded sparrow)
c
              psatl = pv_sc
              dpsatt= dsct             
             endif
c
c get vapor pressure lowering (capillary pressure)
c
            if(ivapl.gt.0) then
               psatl0=psatl
               call vaporl(tl,pcaps,dpcaps,ivapl,delp,ddelt,ddels)
               psatl=psatl0*delp
               dpsatt=dpsatt*delp+psatl0*ddelt
               dpsats=psatl0*ddels
            endif

c
c sat temp as function of sat pres
c
      else
c here tl=is the pressure
       if(iwater_table.ne.1) then
         x=tl
         if(x.lt.0.00123) then
            psatl=10.0
         else if (x.gt.14.5941) then
            psatl=340.0
         else
            x2=x*x
            x3=x2*x
            x4=x3*x
            tfunn=tsa0+tspa1*x+tspa2*x2+tspa3*x3+tspa4*x4
            tfund=tsb0+tspb1*x+tspb2*x2+tspb3*x3+tspb4*x4
            tfun=tfunn/tfund
c calculate derivative wrt p          
            dtpsn=((tspa1+2.*tspa2*x+3.*tspa3*x2+4.*tspa4*x3)*tfund)-
     &           (tfunn*(tspb1+2.*tspb2*x+3.*tspb3*x2+4.*tspb4*x3))
            dtpsd=tfund**2
            dtps=dtpsn/dtpsd  
            dpsatt = dtps
c psatl is the saturation temp
            x=tfun  
c do this for co2 systems
            if(ico2.gt.0) then
c this is the initial guess
               maxitp=10
               eps=1.e-6
               do 20 k=1,maxitp
                  x2=x*x
                  x3=x2*x
                  x4=x3*x
                  pfunn=psa0+psta1*x+psta2*x2+psta3*x3+psta4*x4
                  pfund=psb0+pstb1*x+pstb2*x2+pstb3*x3+pstb4*x4
                  pfun=pfunn/pfund
                  dptsn=((psta1+2.*psta2*x+3.*psta3*x2+
     2              4.*psta4*x3)*pfund)-
     3              (pfunn*(pstb1+2.*pstb2*x+3.*pstb3*x2+4.*pstb4*x3))
                  dpstd=pfund**2
                  dpst=dptsn/dpstd
                  if(ivapl.gt.0) then
                     call vaporl(x,pcaps,dpcaps,ivapl,delp,ddelt,ddels)
                     pfun0=pfun
                     pfun=pfun0*delp
                     dpst=dpst*delp+pfun0*ddelt
                  endif
                  resid=pfun-tl
                  drlp=dpst
                  delp=-resid/drlp
                  x=x+delp
                  if (abs(resid).le.eps) go to 1000
 20            continue
c failed to converge
 1000          continue
               dpsatt=1./dpst
               if(k .gt. maxitp) then
                  if(iout.ne.0)
     &                 write(iout,*)'failed to converge in psatl'
                  if(iptty.ne.0)
     &                 write(iptty,*) 'failed to converge in psatl'
               end if
            endif
            psatl=x
         end if
         else
c tl is the pressure,psatl = temp
c gaz debug 092119   (removed 042521          
c          if(l.eq.1) then
c            write(ierr,*) ' psatl chk ', ifail,tl,psatl,dpsatt
c          endif
c here tl=is the pressure   
c gaz 102420          
          if(tl.lt.pcrit_h2o) then
           call get_h2o_sat_temperature(ifail,tl,psatl,dpsatt)
            dpsats=0.0  
          else
           psatl= tcrit_h2o
           dpsats = 0.0 
           dpsatt = 0.0
         endif
        endif
      endif
 9000 continue
      return
      end

