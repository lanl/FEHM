      subroutine dpdp(iflg)
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
CD1 To provide overall control for a dual porosity/dual permeability
CD1 solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/dpdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 15:07:40 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:47:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:04   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use      Description
CD3 
CD3 iflg         int      I       Control parameter to decide the
CD3                                  purpose for calling this routine
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
CD4 Identifier   Type        Description
CD4 
CD4 
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 Identifier   Type        Description
CD4 
CD4 idpdp, neq, ico2, n, m, nskw, m2, nskw2
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
CD4 rddpdp    N/A        Reads input
CD4 indpdp    N/A        Initializes dpdp variables
CD4 varchk    N/A        Set thermodynamic property values
CD4 gensdp3   N/A        Generates equations for noncondensible gas
CD4                         solution with heat
CD4 gensdp    N/A        Generates equations for pure water or
CD4                         isothermal air-water solution
CD4 crdpdp    N/A        Updates solution variable values
CD4 ctdpdp    N/A        Updates tracer concentration values
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
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
CD9 2.4.9 Double-porosity/double-permeability formulation
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
CPS BEGIN dpdp
CPS 
CPS IF this is a dpdp solution
CPS 
CPS   IF this call is for reading input
CPS     rddpdp - read input
CPS   ELSEIF this call is for initializing dpdp variables
CPS     indpdp - initialize dpdp variables
CPS   ELSEIF this call is for setting thermodynamic properties
CPS     varchk - set thermodynamic properties
CPS   ELSEIF this call is for generating equations
CPS     IF this is a noncondensible gas solution with heat
CPS       gensdp3 - generate eqautions
CPS     ELSE
CPS       gensdp - generate equations
CPS     ENDIF
CPS   ELSEIF this call is for updating solution variable values
CPS     crdpdp - update solution variable values
CPS   ELSEIF this call is to set up the arrays indicating the values...
CPS   ... to be written
CPS     Compute total number of unknowns
CPS     FOR each point to be written
CPS       Compute matrix node number to be written
CPS     ENDFOR
CPS     FOR each point to be written to other file
CPS       Compute matrix node number to be written
CPS     ENDFOR
CPS     Compute total number of nodes written
CPS   ELSEIF this call is for updating tracer concentrations
CPS     ctdpdp - updated tracer concentrations
CPS   ENDIF
CPS 
CPS ENDIF this is a dpdp solution
CPS 
CPS END dpdp
CPS 
C**********************************************************************
C version FEHM5.1J changes
c 17-sept-91
c put in zone capability in input
c 11-oct-91
c added dpdp for tracers
c 16-oct-91
c corrected eqivalence statements in gentdp
c 16-oct-91
c corrected a(i+nmat(4)) to a(i+nmat(2))=0 iv dpdpta
c 17-oct-91
c modified nskw and nskw2 to account for dpdp solution
c 22-nov-91
c modification(major) for symmetry in dpdp formulation
c corrected ...nmat(7)...nmat(3)
c 23-nov-91
c corrected bp(i+nmat.. to bp(i+nrhs..
c 25-nov-91
c changed dpdp for tracer to symmetry
c in dpdpta it8(jm)>idl
c changed to generic intrisic functions
c 29-nov-91
c set volf1 to rp17 in indpdp
c corrected mistaKE in calculation of fractional volumes
c 2-dec-91
c some extra terms deleted in dpdpfa(vapor equations)
c 3-dec-91
c put statement iter(maxit)=2*north before solve4
c 5-dec-91
c equivalenced piv4(16*n0) to piv(n0*16) in gensdp
c changed some calls to solve4 in dpdp
c 18-dec-91
c added bullivants dpdp solver
c worked it under irdof
c testing under air only
c 19-dec-91
c added sol2dp(2 dof solver for dpdp)
c changed meaning of EPN in NLZ4dp
c turned off irdof for tracer dpdp
c got rid of mink and itert count in gentdp
c 19 mar-92
c defined swi in dpdpta
c 26-may-92
c made t9=fid always
c 31-july-92
c zeroed out t6 and t7 and commented sx4d and sx4h
c 7-aug-92
c loaded nmat and nrhs for ico2>0 in indpdp
c 27-aug-92
c changed 20.0(*fdum1) tp parameter variable fdm
c 07-sep-92
c  changed radi in dpdpfh to if block of different dimensions for 2-d or 3-d
c 22-nov-92
c modified gensdp,took out Bullivants programs,added modified solve2
c still normalizing 4 by 4
c 30-nov-92
c temp got rid of call to normalize routine in gensdp
c 2-dec-92
c put maximumon coef1 in routines for numerics
c 30-dec-92
c took out max perm(isotropic =max) loop in INDPDP
c 1-4-93
c put sor iteration in rdof_dp2a
c 27-jan-93
c added new dpdp3 and solve3_dp routines
c 27-jan-93
c put in mod solve routines,set iter=3*north
c 5-feb-93
c changed eqivalences in gensdp and gensdp3
c called solven from gensdp3
c 11-feb-93
c called nornal_uc from gensdp3
c 12-feb-93
c modified geneqc (a little ) so it would work for dpdp
c added air-water vapor diffusion to transfer terms
c 23-feb-93 llt
c changed equivalences to pointers
c 19-mar-93 llt
c replaced gotos embedded in code with a goto at end of subroutine
C***********************************************************************

      use comdi
      use comdti
      use comai
      implicit none

      integer iflg,i

      if(idpdp.ne.0) then
c     
c     organize dpdp solution
c     
         if(iflg.eq.0) then
c     call input routine
            call rddpdp
         else if(iflg.eq.1) then
c     initialize some dpdp variables
            call indpdp
         else if(iflg.eq.2) then
c     call to thermodynamics routines
            call varchk(0,neq)
         else if(iflg.eq.3) then
c     generate algebraic equations,load matrices
            if(ico2.gt.0) then
               call gensdp3
            else
	        if(jswitch.ne.0) then
               call gensdp_switch
	        else
               call gensdp
	        endif
            endif
         else if(iflg.eq.4) then
c     get variable correction
            call crdpdp
         else if(iflg.eq.5) then
c     organize dpdp printout
         else if(iflg.eq.6) then
c     set n=2*neq
            n=neq+neq
c     
c     add nodes to nskw and nskw2
c     
            do i=1,m
               nskw(m+i)=nskw(i)+neq
            enddo
            m=2*m
            do i=1,m2
               nskw2(m2+i)=nskw2(i)+neq
            enddo
            m2=2*m2

c     generate algebraic equations(tracer equations),load matrices
c     routine dpdpta is called by gentdp
c     routine solve2 is called by gentdp
c     
         else if(iflg.eq.8) then
c     get varible correction(tracer)
            call ctdpdp
         endif
      endif

      return
      end
