      subroutine solstore( mi )
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
CD1 To compute the tracer mass storage and sorption
CD1 terms of the residuals equations and their derivatives.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 07-26-93      B. Robinson   22      Extracted from old version of
CD2                                     thermc to make code more modular
CD2
CD2 $Log:   /pvcs.config/fehm90/src/solstore.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 02 11:44:46 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   01/28/95 14:20:58   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.1   03/18/94 16:16:00   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:54   pvcs
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
CD3 mi           int      I      Position in concentration arrays for
CD3                              current unknown
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
CD5 conc_min     real*8      Minimum concentration to be used in
CD5                          sorption related floating point
CD5                          calculations
CD5 rtol         real*8      Minimum value of density x saturation to
CD5                          be used in subsequent floating point
CD5                          calculations
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5
CD5 dtin         real*8      Reciprocal of time step
CD5 dencd        real*8      Intermediate variable in calculation
CD5 mim          int         Current node number
CD5 vap_conc     real*8      Concentration of solute in vapor
CD5 dvap_conc    real*8      Derivative of vap_conc with liquid
CD5                              concentration
CD5 betam1l      real*8      Sorption exponent minus 1
CD5 betam1v      real*8      Sorption exponent minus 1
CD5 danld        real*8      Factor used in calculation
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
CD7 It is assumed that upon entry the sorption and reaction parameters
CD7 are properly specified so that no error checks are required to
CD7 catch invalid floating point operations.
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
CD9 2.4.5 Adsorbing solutes
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
CPS BEGIN solstore
CPS
CPS   IF a solid species
CPS      Set storage terms
CPS   ELSE  the species is a liquid, vapor 
CPS      zero out the other phase
CPS     calculate storage terms
CPS     IF current species is sorbing
CPS       calculate sorption terms
CPS     ENDIF
CPS   ENDIF
CPS   
CPS END solstore
CPS
C**********************************************************************

      use comci
      use comdi
      use comrxni
      use comfi
      use comcouple
      use comdti
      use comai
      implicit none

      real*8 dtin
      real*8 dencd
      real*8 vap_conc
      real*8 h_const
      real*8 betam1l
      real*8 betam1v
      real*8 danld
      real*8 dvap_conc
      real*8 conc_subst
      integer mim
      integer mi
      dtin=1./dtotc
      mim=mi-npt(nsp)
      an(mi) = sign(max(abs(an(mi)),1.0d-90),an(mi))
      anl(mi)= sign(max(abs(anl(mi)),1.0d-90),anl(mi))
c      an(mi) = max(1.0d-90,an(mi))
c      anl(mi) = max(1.0d-90,anl(mi))
      if(icns(nsp).eq.0)then
         dencd=anl(mi)*denr(mim)
         denci(mi)=(dencd-dench(mi))*dtin
         akc(mi)=dtin*denr(mim)
         scl(mi) = 0.
         dsccl(mi) = 0.
         scv(mi) = 0.
         dsccv(mi) = 0.
CPS   ELSE the species is a liquid, vapor 
      else
CPS   IF a liquid or vapor species
         if(abs(icns(nsp)).eq.1)then
            vap_conc= 0
            dvap_conc = 0
CPS   ELSE current species is a Henry's species
         else
CPS   calculate conc in vapor phase using Henry's Law
            if(henry_model(nsp).eq.1)then
               h_const= a_henry(nsp)*
     2              exp(dh_henry(nsp)/
     3              (gas_const*0.001)*(1/298.16-1/
     4              (t(mim)+temp_conv)))
            elseif(henry_model(nsp).eq.2)then
               h_const = 10**(hawwa(nsp,1)+
     2              hawwa(nsp,2)*(t(mim)+temp_conv)+
     3              hawwa(nsp,3)/(t(mim)+temp_conv)+
     4              hawwa(nsp,4)*dlog10(t(mim)+temp_conv)+
     5              hawwa(nsp,5)/(t(mim)+temp_conv)**2)
               h_const= (101325*rolf(mim)*1e-3)/(h_const*1e6*
     2              mw_water)
            elseif(henry_model(nsp).eq.3)then
               h_const= (phi(mim) - pci(mim)) * dh_henry(nsp)
            endif
            dvap_conc = (mw_water*h_const)/(phi(mim)*avgmolwt(mim))
            vap_conc = anl(mi)*dvap_conc
            anv(mi) = vap_conc
CPS   ENDIF
         endif
CPS   calculate storage terms
         dencd=danl(mi)*ps_trac(mim)*anl(mi)+
     2        danv(mi)*ps_trac(mim)*vap_conc
         akc(mi) = dtin * (danl(mi)*ps_trac(mim) +
     2        danv(mi)*ps_trac(mim)*dvap_conc )
         scl(mi) = 0.
         scv(mi) = 0.
         dsccl(mi) = 0.
         dsccv(mi) = 0.
         
CPS   IF species is sorbing        
         if(iadsfl(nsp,itrc(mi)).ne.0.or.iadsfv(nsp,itrc(mi)).ne.0)then
CPS   calculate sorption terms
            conc_subst = max(anl(mi),1.d-90)
            scl(mi) = (denr(mim) * a1adfl(nsp,itrc(mi)) * conc_subst
     2           **betadfl(nsp,itrc(mi)) / (1.0 +
     3           a2adfl(nsp,itrc(mi)) *
     4           conc_subst**betadfl(nsp,itrc(mi)) ))
            conc_subst = max(vap_conc,1.d-90)
            scv(mi) = (denr(mim) * a1adfv(nsp,itrc(mi)) * conc_subst
     2           **betadfv(nsp,itrc(mi)) / (1.0 + a2adfv(nsp,itrc(mi))
     3           * conc_subst**betadfv(nsp,itrc(mi)) ))
            betam1l = betadfl(nsp,itrc(mi)) - 1
            betam1v = betadfv(nsp,itrc(mi)) - 1
            danld = 1.0
            conc_subst = max(anl(mi),1.d-90)
            dsccl(mi) = denr(mim) * a1adfl(nsp,itrc(mi)) *
     2           betadfl(nsp,itrc(mi)) *
     3           conc_subst**betam1l / ( 1.0 + a2adfl(nsp,itrc(mi)) *
     4           conc_subst**betadfl(nsp,itrc(mi)))**2
            conc_subst = max(vap_conc,1.d-90)
            dsccv(mi) = dvap_conc*
     2           denr(mim)*a1adfv(nsp,itrc(mi))*betadfv(nsp,itrc(mi))*
     3           conc_subst**betam1v/ ( 1.0 +
     4           a2adfv(nsp,itrc(mi)) * 
     5           conc_subst**betadfl(nsp,itrc(mi)))**2
            dencd = dencd + scl(mi) + scv(mi)
            akc(mi) = akc(mi) + dtin * ( dsccl(mi) + dsccv(mi) )
CPS   ENDIF
         endif
         denci(mi) = ( dencd - dench(mi) ) * dtin
CPS   ENDIF
      end if
      return
      end

