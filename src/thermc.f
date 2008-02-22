      subroutine thermc(ndummy)
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
CD1 To set the tracer mass storage, sorption, and chemical reaction
CD1 terms of the residuals equations and their derivatives.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 ?            G. Zyvoloski   N/A     Initial implementation -
CD2                                     however, previous non-YMP
CD2                                     versions of FEHM exist, and the
CD2                                     current version may differ from
CD2                                     these
CD2 04-07-93      B. Robinson   N/A     Revised to implement the
CD2                                     chemical reaction models
CD2                                     between solutes
CD2
CD2 $Log:   /pvcs.config/fehm90/src/thermc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Fri May 24 09:55:52 1996   hend
CD2 Updated fac parameter manipulations
CD2 
CD2    Rev 1.5   Fri Feb 16 13:01:52 1996   zvd
CD2 Added requirements.
CD2 
CD2    Rev 1.4   04/10/95 11:30:46   robinson
CD2 Solute user sub now called for setting source/sink values
CD2 
CD2    Rev 1.3   04/03/95 08:48:18   robinson
CD2 Corrected source/sink term for ngas problems
CD2 
CD2    Rev 1.2   01/28/95 14:21:04   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.1   03/18/94 16:16:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:38   pvcs
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
CD3 ndummy       int      I      Flag denoting whether main nodes or
CD3                              dual porosity matrix nodes are being
CD3                              used
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
CD3 File with number iatty O    File used to write warning and error
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
CD4 solstore        N/A      Computes the storage of solute at a given
CD4                             node
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
CD5                             sorption related floating point
CD5                             calculations
CD5 rtol         real*8      Minimum value of density x saturation to
CD5                             be used in subsequent floating point
CD5                             calculations
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
CD6   Looping through each solute, it checks to see if it is a liquid
CD6   or vapor solute, and loops through each node, setting variables
CD6   based on fluid density and saturation accordingly (danset).
CD6   
CD6   Looping through each node, the code computes the storage and
CD6   sorption terms based on whether the solute is sorbing of
CD6   conservative (solstore).
CD6   
CD6   It next sets a flag needed to identify the correct flow velocity
CD6   based on whether the solute is liquid or vapor.
CD6   
CD6   Then, it sets a flag that serves the purpose of zeroing out
CD6   advective-dspersion if a matrix node is being handled.
CD6   
CD6   Next, within a loop over each node, the code
CD6   
CD6     Computes the dispersion coefficients;
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
CD6   Finally, a loop over each solute is used to reset temporary
CD6   storage arrays to their values upon entry of this routine.
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
CD9 2.3.7 Sources and sinks
CD9 2.4.5 Adsorbing solutes
CD9 2.4.6 Multiple, interacting solutes
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
CPS BEGIN thermc
CPS
CPS Initialize time step parameter
CPS
CPS FOR each solute
CPS
CPS   IF a liquid,  solid  or a henry's law tracer
CPS   
CPS     FOR each node
CPS       Set terms needed in future calculations
CPS     ENDFOR 
CPS     
CPS   ELSE this is a gas phase tracer
CPS   
CPS     FOR each node
CPS       Set terms needed in future calculations
CPS     ENDFOR 
CPS     
CPS   ENDIF
CPS
CPS ENDFOR 
CPS 
CPS FOR each node
CPS   solstore - compute solute storage in fluid and on rock surface
CPS ENDFOR each node
CPS
CPS FOR each node
CPS   Initialize arrays and indexes
CPS ENDFOR
CPS 
CPS END thermc
CPS
C**********************************************************************

      use combi
      use comchem
      use comci
      use comdi
      use comgi
      use comrxni
      use comfi
      use comcouple
      use comdti
      use comai
      use davidi, only : irdof
      implicit none

      integer ndummy
      integer nsolute
      integer i
      integer mi
      integer mim
      real*8 anv_subst, h_const

      do nsolute = 1, nspeci
**** Add Henry's law, do it like liquid part
CPS   IF a liquid,  solid  or a henry's law tracer
         if(icns(nsolute).ge.0.or.abs(icns(nsolute)).eq.2) then
            do i = 1, neq
               mi=i+ndummy+npt(nsolute)
               mim=mi-npt(nsolute)
               anl(mi)=an(mi)
               if(abs(icns(nsolute)).ne.2) anv(mi)=0.0
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  danv(mi) = max( rovf(mim) * ( 1. - s(mim) ), rtol )
                  danl(mi) = max( rolf(mim) * s(mim), rtol )
               else
                  danv(mi) = 0.
                  danl(mi) = max( rolf(mim), rtol )
               end if
            end do
         else
            do i = 1, neq
               mi=i+ndummy+npt(nsolute)
               anv_subst = anv(mi)
               mim=mi-npt(nsolute)
               anv(mi) = an(mi)
               anl(mi) = anv(mi)
               danv(mi) = max( rolf(mim) * s(mim), rtol )
               danl(mi) = max( rovf(mim) * ( 1. - s(mim) ), rtol )
            end do
         end if
      end do

      do i = 1, neq
         mi = i + ndummy + npt(nsp)
         call solstore( mi )
      end do

      do nsolute = 1, nspeci
         if(icns(nsolute).eq.1) then
            do i = 1, neq
               mi=i+ndummy+npt(nsolute)
               danl(mi) = 1.
               danv(mi) = 0.
            end do
      
         else if(icns(nsolute).eq.-1) then
            do i = 1, neq
               mi=i+ndummy+npt(nsolute)
               danl(mi) = 0.
               danv(mi) = 1.
            end do
            
         else if(abs(icns(nsolute)).eq.2)then
	        if(henry_model(nsolute).eq.1) then
                  h_const= a_henry(nsolute)*
     2                 exp(dh_henry(nsolute)/
     3                 (gas_const*0.001)*(1/298.16-1/
     4                 (t(mim)+temp_conv)))
               else if(henry_model(nsolute).eq.2) then
                  h_const = 10**(hawwa(nsolute,1)+
     2                 hawwa(nsolute,2)*(t(mim)+
     3                 temp_conv)+hawwa(nsolute,3)
     4                 /(t(mim)+temp_conv)+hawwa
     5                 (nsolute,4)*dlog10(t(mim)
     6                 +temp_conv)+hawwa(nsolute,5)
     7                 /(t(mim)+temp_conv)**2)
                  h_const= (101325*rolf(mim)*1e-3)/(h_const*1e6*
     2                 mw_water)
               else if(henry_model(nsolute).eq.3) then
                  h_const= (phi(mim) - pci(mim)) * dh_henry(nsolute)
               endif
               do i = 1, neq
                  mi=i+ndummy+npt(nsolute)
                  mim = i + ndummy
                  danl(mi) = 1.
                  danv(mi) = (mw_water*h_const)/(phi(mim)*
     2                 avgmolwt(mim))
               end do
         endif
      end do

      return
      end




