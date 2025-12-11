      subroutine ss_trans(i,rc_ss,drc_ss)
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
CD1 To compute source/sink terms
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 Initial implementation: ?, Programmer: Hari Viswanathan
CD2 source/sink stuff was in thermc
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/ss_trans.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
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
CD4 nspeci, icns, neq, npt, anl, anv, danl, danv, rovf, rolf, npn,
CD4 drc, rc, nsp, itrc, displx, disply, displz, n, ps, pnx, pny, pnz,
CD4 diffmfl, dispvx, dispvy, dispvz, tclx, tcly, tclz, tcvx, tcvy, tcvz,
CD4 diffmfv, ico2, sk, qh, s, days, t1sk, t2sk, a_henry, dh_henry,
CD4 gas_const, t, temp_conv, phi, avgmolwt, mw_water, cnsk,  rcss
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
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
CD5 fac          real*8      Parameter used in source/sink term of
CD5                             tracer mass balance
CD5 rat          real*8      Parameter used in chemical reaction term
CD5                             of tracer mass balance
CD5 i            int         Do-loop index parameter over all nodes
CD5 mi           int         Index for the concentration arrays
CD5 mim          int         Index for the parameter value arrays
CD5 nsolute      int         Do loop index parameter for looping over
CD5                             all species
CD5 ireg         int         Index denoting current reaction and
CD5                             sorption properties group number
CD5                            rate terms
CD5 anv_subst    real*8      Vapor concentration at current node
CD5 srmiml       real*8      Liquid inlet mass flow rate at node
CD5 srmimv       real*8      Vapor inlet mass flow rate at node
CD5 conc_subst   real*8      Concentration at current node
CD5 drc_subst    real*8      Temporary storage for derivative term
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
CD6 This routine performs the following functions:
CD6 
CD6   Next, within a loop over each node, the code
CD6   
CD6     
CD6     Decides whether it is a liquid or vapor phase solute, and sets
CD6     the flow rate at this node accordingly.
CD6     
CD6     Next, an IF block determines if there is fluid injection at
CD6     this node, and if there is, it either sets the solute source
CD6     term or sets it to 0 if it is not during the time of solute
CD6     injection.  Otherwise, the node may be a sink, and the code
CD6     computes the solute sink term (it will be 0 if there is no
CD6     outflow at the node).
CD6   
CD6     Then, if the constant concentration option is specified the
CD6     code sets the source term and derivative accordingly.
CD6     
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
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
CD9 2.3.7 Sources and Sinks
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
CPS BEGIN ss_trans
CPS
CPS   FOR each node
CPS
CPS
CPS 
CPS     IF this is an isothermal air-water simulation
CPS       Set liquid and vapor source/sink flow rates
CPS     ELSE
CPS       Set liquid and vapor source/sink flow rates
CPS     ENDIF
CPS
CPS     IF liquid is entering the system
CPS
CPS       IF solute can enter with the liquid
CPS         IF the current time is during the solute injection period
CPS           Compute contribution to solute source term
CPS         ENDIF
CPS       ENDIF
CPS   
CPS     ELSE liquid is leaving the system (or the flow rate is 0)
CPS 
CPS       IF this is not a vapor-only solute
CPS         Compute contribution to solute sink term
CPS       ENDIF
CPS   
CPS     ENDIF
CPS 
CPS     IF vapor is entering the system
CPS 
CPS       IF solute can enter with the vapor
CPS         IF the current time is during the solute injection period
CPS           Compute contribution to solute source term
CPS         ENDIF
CPS       ENDIF
CPS   
CPS     ELSE vapor is leaving the system (or the flow rate is 0)
CPS 
CPS       IF this is a vapor-only solute
CPS         Set vapor concentration value
CPS       ELSEIF this is a Henry's Law solute
CPS         Compute exit vapor concentration
CPS       ELSE there is no solute leaving the system via the vapor
CPS         Set vapor concentration value to 0
CPS       ENDIF
CPS   
CPS       Compute contribution to solute sink term
CPS   
CPS     ENDIF
CPS     
CPS     IF the node is to be held at a constant concentration
CPS     
CPS       IF the current time is during the solute injection period
CPS         IF this is not a Henry's Law species with vapor...
CPS         ... concentrations specified
CPS           Compute source term and derivative to keep...
CPS           ... concentration constant
CPS         ELSE
CPS           Compute equilibrium liquid concentration
CPS           Compute source term and derivative to keep...
CPS           ... concentration constant
CPS         ENDIF
CPS       ENDIF
CPS       
CPS     ENDIF
CPS 
CPS   ENDFOR each node
CPS
CPS ENDIF (the previous section is skipped for a solid species)
CPS FOR each node
CPS 
CPS   IF it is a solid tracer
CPS     set ireg = 1
CPS   ELSE it is a liquid or vapor tracer
CPS     set ireg to itrcd
CPS   ENDIF
CPS   Compute integers pointers for rate calculations
CPS   
CPS
CPS ENDFOR each node
CPS END ss_trans
CPS
C**********************************************************************

      use combi
      use comci
      use comdi
      use comgi
      use comrxni
      use comfi
      use comcouple
      use comdti
      use comai
      use comchem, only: ncpnt,cpntnam,co2_couple,pcpnt
      use comco2,only: yc,icarb,wat_prop,carbon_tracer,rate
      use compart
      use davidi
      implicit none

      integer i,ic
      integer mi
      integer mim
      real*8 srmiml
      real*8 srmimv
      real*8 conc_subst
      real*8 drc_subst
      real*8 fac,crit
      real*8 rc_ss
      real*8 drc_ss
      real*8 h_const
      real*8 dvap_conc
      real*8 disco2_flow, mass_rock
      mi=i+npn
      mim = mi - npn
      rcss(mi) = 0.
      rc_ss = 0
      drc_ss = 0
      if (icns(nsp).ne.0)then
c     
c   IF this is an isothermal air-water simualtion
c
         if( ico2 .lt. 0 ) then
c
c     Set liquid and vapor source/sink flow rates
c
            srmiml = sk(mim)
            srmimv = qh(mim)
c   ELSE
c
         elseif( ico2 .gt. 0 ) then
c
c     Set liquid and vapor source/sink flow rates
c
            srmiml = sk(mim)
            srmimv = qc(mim)
         else
            if (irdof .ne. 13 .or. ifree .ne. 0) then              
               srmiml = sk(mim) * s(mim)
            else 
               srmiml = sk(mim)
            end if
            srmimv = sk(mim) - srmiml
c
c   ENDIF
c
         end if
c   
c   IF liquid is entering the system
c   
         if( srmiml .lt. 0. ) then
c
c     IF solute can enter with the liquid
c
            if( icns(nsp) .gt. 0 ) then
c
c       IF the current time is during the solute injection period
c
               if( days .gt. abs(t1sk(mi)) .and.
     2              days .le. abs(t2sk(mi)) ) then
c
c         Compute contribution to solute source term
c

                  rcss(mi) = rcss(mi)+cnsk(mi) * srmiml
                  rc_ss = rc_ss+cnsk(mi)*srmiml
               else
c gaz 031520
                  rcss(mi) = rcss(mi)+cnsk_background(mi) * srmiml
                  rc_ss = rc_ss+cnsk_background(mi)*srmiml       
c
c       ENDIF
c
               end if
c
c     ENDIF
c
            end if
c     
c   ELSE liquid is leaving the system (or the flow rate is 0)
c
         else
c   
c     IF this is not a vapor-only solute
c     
            if( icns(nsp) .ne. -1 ) then
c     
c     Compute contribution to solute sink term
c     

c     Add option to accumulate solute rather than having
c     it leave the system BAR 4-28-99

c     For the accumulation option, this if block is false because
c     the pcnsk flag is set to 1.

               if(pcnsk(mi).le.0.) then
                  rcss(mi) =  rcss(mi) + anl(mi) * srmiml
                  rc_ss = rc_ss + anl(mi)*srmiml
                  drc_ss = drc_ss + srmiml
               end if
c
c     ENDIF
c
            end if
c     
c   ENDIF
c
         end if
c   
c   IF vapor is entering the system
c   
         if( srmimv .lt. 0. ) then
c
c     IF solute can enter with the vapor
c     
            if( icns(nsp) .lt. 0 ) then
c
c       IF the current time is during the solute injection period
c
               if( days .gt. abs(t1sk(mi)) .and. 
     2              days .le. abs(t2sk(mi)) ) then
c
c         Compute contribution to solute source term
c
                  rcss(mi) = rcss(mi) + cnsk(mi) * srmimv
                  rc_ss = rc_ss + cnsk(mi)*srmimv
c
c       ENDIF
c
               end if
c
c     ENDIF
c
            end if
c     
c   ELSE vapor is leaving the system (or the flow rate is 0)
c
         else
c   
c   
c     IF this is a vapor-only solute
c
            if( icns(nsp) .eq. -1 ) then
c
c       Set vapor concentration value
c
               conc_subst = anl(mi)
               drc_subst = srmimv
c
c     ELSEIF this is a Henry's Law solute
c
            else if( abs(icns(nsp)) .eq. 2 ) then
c
c       Compute exit vapor concentration
c
               if(henry_model(nsp).eq.1) then
                  h_const= a_henry(nsp)*
     2                 exp(dh_henry(nsp)/
     3                 (gas_const*0.001)*(1/298.16-1/
     4                 (t(mim)+temp_conv)))
               else if(henry_model(nsp).eq.2) then
                  h_const = 10**(hawwa(nsp,1)+
     2                 hawwa(nsp,2)*(t(mim)+
     3                 temp_conv)+hawwa(nsp,3)
     4                 /(t(mim)+temp_conv)+hawwa
     5                 (nsp,4)*dlog10(t(mim)
     6                 +temp_conv)+hawwa(nsp,5)
     7                 /(t(mim)+temp_conv)**2)
                  h_const= (101325*rolf(mim)*1e-3)/(h_const*1e6*
     2                 mw_water)
               else if(henry_model(nsp).eq.3) then
                  h_const= (phi(mim) - pci(mim)) * dh_henry(nsp)
               endif
               dvap_conc = (mw_water*h_const)/(phi(mim)*
     2              avgmolwt(mim))
               conc_subst = anl(mi)*dvap_conc
               drc_subst = dvap_conc * srmimv
c


c
c     ELSE there is no solute leaving the system via the vapor
c
            else
c
c       Set vapor concentration value to 0
c
               conc_subst = 0.
               drc_subst = 0.
c
c     ENDIF
c     
            end if
c
c     Compute contribution to solute sink term
c
            rcss(mi) = rcss(mi) + conc_subst * srmimv
            rc_ss = rc_ss + conc_subst * srmimv
            drc_ss = drc_ss + drc_subst
c     
c   ENDIF
c
         end if
c CO2 FLOW SOLUTION COUPLING          
               if(carbon_tracer.eq.nsp)then
c check that we want couple co2 flow solubility to the trac macro
c convert mass fraction to moles/kg water 
c this was original (Hari)
                  disco2_flow = (yc(i)/(1-yc(i)))/0.044 
			  
c	         if(disco2_flow.gt.1.e-5.and.disco2_flow.gt.anl(mi)) then
 	          if(disco2_flow.gt.1.e-5) then	  

                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     fac=(sx1(mim)*ps_trac(mim)*s(mim)*denr(mim))/dtotc
                  else
                     fac=(sx1(mim)*ps_trac(mim)*denr(mim))/dtotc
                  end if	                  
                          fac=fac*max(1,month)*1e6

        rcss(mi) = rcss(mi) + fac*(anl(mi) - (disco2_flow) ) 
       rc_ss = rc_ss + fac* (anl(mi) - (disco2_flow) )

                  drc_ss = fac             
c                end do
                endif
                endif                     
c
CPS     IF the node is to be held at a constant concentration
c
         if( pcnsk(mi) .lt. 0. ) then
            call userc(2,i,rc_ss,drc_ss)
c
CPS       IF the current time is during the solute injection period
c
            if( days .gt. abs(t1sk(mi)) .and. 
     2           days .le. abs(t2sk(mi)) ) then

c
CPS       IF this is not a Henry's Law species with vapor...
CPS       ... concentrations specified
c
               if( icns(nsp) .ne. -2 ) then
c
CPS         Compute source term and derivative to keep concentration...
CPS         ... constant
c
c     Changed the boundary condition to consider
c     the tracer time step rather than the heat and
c     mass time step - BAR 12-20-98

                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     fac=(sx1(mim)*ps_trac(mim)*s(mim)*denr(mim))/dtotc
                  else
                     fac=(sx1(mim)*ps_trac(mim)*denr(mim))/dtotc
                  end if
                  fac=fac*max(1,month)*1e6
                  rcss(mi) = rcss(mi) + fac * ( anl(mi) + cnsk(mi) )
                  rc_ss = rc_ss + fac* (anl(mi) + cnsk(mi))
                  drc_ss = fac
c
CPS       ELSE
c
               else
c
CPS         Compute equilibrium liquid concentration
c
                  if(henry_model(nsp).eq.1) then
                     h_const= a_henry(nsp)*
     2                    exp(dh_henry(nsp)/
     3                    (gas_const*0.001)*(1/298.16-1/
     4                    (t(mim)+temp_conv)))
                  else if(henry_model(nsp).eq.2) then
                     h_const = 10**(hawwa(nsp,1)+
     2                    hawwa(nsp,2)*(t(mim)+
     3                    temp_conv)+hawwa(nsp,3)
     4                    /(t(mim)+temp_conv)+hawwa
     5                    (nsp,4)*dlog10(t(mim)
     6                    +temp_conv)+hawwa(nsp,5)
     7                    /(t(mim)+temp_conv)**2)
                     h_const= (101325*rolf(mim)*1e-3)/(h_const*1e6*
     2                       mw_water)
                  else if(henry_model(nsp).eq.3) then
                     h_const= (phi(mim) - pci(mim)) * dh_henry(nsp)
                  endif
                  dvap_conc = (mw_water*h_const)/(phi(mim)*
     2                 avgmolwt(mim))
                  conc_subst = cnsk(mi)/dvap_conc

                  fac=(sx1(mim)*ps_trac(mim)*s(mim)*denr(mim))/dtotc
                  fac=fac*max(1,month)*1e6
                  rcss(mi) = rcss(mi) + fac * (anl(mi) + conc_subst)
                  rc_ss = rc_ss + fac * (anl(mi) + conc_subst)
                  drc_ss = fac
c
               endif
c
            end if
c
CPS     ENDIF
c
         end if
      end if
      call userc(2,i,rc_ss,drc_ss)
      rcss(mi) = rc_ss
      end
