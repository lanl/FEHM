          real*8 function psat(tl,dpsatt,isatf)
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
CD2 ?            G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/psat.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Feb 05 11:19:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:11:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:32   pvcs
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
CD3 dpsatt       real*8  O       Derivative of saturation pressure
CD3                                  with temperature
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
CD4 pstb4, psb0,  tspa1,tspa2,tspa3,tspa4, tsa0, tspb1,tspb2,
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
CD5 delp         real*8      Term used in iterative solution
CD5 tfun         real*8      Intermediate term in calculation
CD5 tfunn        real*8      Intermediate term in calculation
CD5 tfund        real*8      Intermediate term in calculation
CD5 maxitp       int         Maximum number of iterations allowed
CD5 k            int         Current iteration number in do loop
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
CPS BEGIN psat
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
CPS     Compute saturation pressure and derivatives
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
CPS         Compute temperature
CPS         
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
CPS END psat
CPS
C**********************************************************************
c version FEHM5.1J changes
c 20-june-91 changed Gangi model so that if iporos=-2
c then the initial porosity and permeability are calculated
c so that when the gangi parameters are used the input
c porosity and permeabilities are used
c 25-june-91
c corrected some derivatives for gangi model in POROSI
c 25-june-91
c initialized gangi parameters so if no defined on input
c the gangi model will behave like constant porosity
c and permeability
c 29-june-91
c in VFCAL checked for wellbore solution first
c 2-july-91
c if VFCAL deleted references to wellbore solution
c 2-july-91
c this will be handled in sub WELBOR
c put call in to welbor in routine rlperm
c 2-july-91
c put call to welbor for call to rock in porosi
c 2-july-91
c put call to welbor in thermw for calc of volume factors
c 12-july-91
c corrected the corey curves (sl-rp1)
c 15-july-91
c in thermw took out if blocks that checked for phase in
c calculating transmissibilities
c 6-august-91
c changed thrmwc to lump rlperms in dil,etc
c 13-august-91
c changed psatf(t) to tsatf(p) or viceversa
c 14-august-91
c put return in psatf for pure watre systems
c the co2 stuff uses psa coeff.,pure water uses tsa coef.
c 19-august-91
c equivalenced rlperms to arrays in commom/fcc/ in COMCI
c preparatory to removing the storage places
c 20-august-91
c moved evaluation on rlperms to beginning of 100 loop
c in ther routines so the would not be overwritten
c when vfcal called need new dil ,dilp,dile
c 19-sept-91
c changed the input sections so blanks can be read
c 8-oct-91
c changed iinp to inpt and got rid of iinp=inpt
c 24-nov-91
c in zone input changed neq to n
c 29-nov-91
c changed rl s so applicable to DP and DPDP
c 3-dec-91
c looked(didnt find) for intrinic functions
c 20-jan-92
c in psat put in limits after else (if x.lt.0.00123,etc)
c 11-mar-92
c added subroutine vcon for nonlinear thermal conductivities
c 27-mar-92
c took away double precision function (let implicit work)
c 3-apr-92
c added call to thrmwcl for vapor pressure lowering
c 14-apr-92
c made rlperm go through rl etc if van Genucten model used
c 2-may-92
c moved specified noncondensible bc out of if block
c 3-may-92
c added capability in thrmwc and thermw to hlod sl constant
c enabled  in macro hflx with qflxm(mi)<0
c put temp bc in thrmwc if ieosd=2 and pbounc<0
c 21-may-92
c sub rlperm: added var iupk>0 to upwind intrisic perm
c 17-aug-92
c calculated smcut as normalized saturations for use in spline fit
c changed back sl to ss in spline fits(not noted before!)
c 27-aug-92
c parameterizes smcutm ans supm and iupk in rlperm
c 18-nov-92
c put cutoff in linear cap model in routine cappr
c 19-jan-93
c added derivative to air-viscosity
c 25-jan-93
c added subroutine dvacalc to end and put call at end of thrmwc
c 26-jan-93
c added enva in thrmwc for air enthalpy and derivatives
c 27-jan-93
c added enva-ens in thrmwc for air enthalpy and derivatives
c 2-feb-93 
c put return in dvacalc for iadif=0
c 12-feb-93
c added houseworth fix ro rlperm
c changed definition of equivalent permeability
c did not implement change for it=0 ,all
c node must have a relative perm model
c need to check for this in datchk
c 16-feb-93
c changed if(kq.eq.0) deqh(mi)=0.0 to if(..and.qflux(mi).eq.0.0 in thermw
c 26-mar-93
c changed derivative of hflux in thrmwc
c 30-mar-93
c c zero out flow terms thermo routines(and derivatives)
c 31-mar-93
c put min concentrations of air =0 in both phases and zwreod derivatives.
c this may slow convergence
c 5-april-93
c commeneted out cutoff for roc,xnl,and psat
c 7-may-93 llt
c changed equivalences to pointers
c 7-may-93 llt
c replaced gotos embedded in code with a goto at end of subroutine
C***********************************************************************

      use comii
      use comdti
      use comai
      implicit none

      integer isatf,k,maxitp
      real*8 tl,dpsatt,x,x2,x3,x4,pfun,pfunn,pfund,dpst,dptsn,dpstd
      real*8 delp,tfun,tfunn,tfund,resid,drlp

      psat=0.0
      dpsatt=0.0
      if(ice.ne.0) then
         goto 9000
      end if
c ev3 is the reference value for vapor density,initialized(1.) in main.s
c sat pressure as function of sat temp
      if(isatf.le.0) then
c check for limiting values
         if(tl.lt.10.) then
            psat=0.00123
         else if(tl.gt.340.0) then
            psat=14.5941
         elseif(ipsat.eq.0) then
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
            psat=pfun
            dpsatt=dpst
         end if
c
c sat temp as function of sat pres
c
      else
c here tl=is the pressure
         x=tl
         if(x.lt.0.00123) then
            psat=10.0
         else if (x.gt.14.5941) then
            psat=340.0
         else
            x2=x*x
            x3=x2*x
            x4=x3*x
            tfunn=tsa0+tspa1*x+tspa2*x2+tspa3*x3+tspa4*x4
            tfund=tsb0+tspb1*x+tspb2*x2+tspb3*x3+tspb4*x4
            tfun=tfunn/tfund
c psat is the saturation temp
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
                  resid=pfun-tl
                  drlp=dpst
                  delp=-resid/drlp
                  x=x+delp
                  if (abs(resid).le.eps) go to 1000
 20            continue
c failed to converge
 1000          continue
               if(k .gt. maxitp) then
                  if(iout.ne.0) 
     &                 write(iout,*) 'failed to converge in psat'
                  if(iptty.ne.0)
     *                 write(iptty,*) 'failed to converge in psat'
               end if
            endif
            psat=x
            dpsatt=0.0
         end if
      endif
 9000 continue

      return
      end
