      subroutine bnswer
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
CD1 To call routines to assemble finite element equations and solve for
CD1 the Newton-Raphson equations. 
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
CD2 $Log:   /pvcs.config/fehm90/src/bnswer.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Fri Mar 01 14:41:52 1996   gaz
CD2 took out " strd=1. 'strd=1. in N_R loop
CD2 
CD2    Rev 1.5   Mon Jan 29 13:15:24 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.4   Mon Jan 29 12:40:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   06/01/95 17:30:30   gaz
CD2 minor mods for minimum iterations
CD2 
CD2    Rev 1.2   03/24/95 00:12:18   gaz
CD2 gaz deleted reference to tmch near end
CD2 
CD2    Rev 1.1   03/18/94 15:45:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:30   pvcs
CD2 original version in process of being certified
CD2 
c 3/20/95 gaz took out if(fdum.le.tmch)... covered in gensl1,2 etc:
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 None
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
CD4 itert, strd, iad, maxit, mlz, idof, idpdp, ico2, fdum, f0, iad,
CD4 itotal, tmch, epe
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
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
CD4 varchk  N/A      Updates variable states, apply corrections to
CD4                     solution vectors
CD4 outbnd  N/A      Determines if any variable values are out of bounds
CD4 gensl3  N/A      Obtain heat transfer solution
CD4 gensl4  N/A      Obtain single porosity solution
CD4 dpdp    N/A      Obtain dpdp solution
CD4 airctr  N/A      Obtain single porosity, isothermal air-water
CD4                     solution
CD4 icectr  N/A      Obtain single porosity, solid-liquid-gas       
CD4                     solution
CD4 gensl1  N/A       Obtain single phase, single porosity solution
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
CD9 2.5.2 Solve nonlinear equation set at each time step  
CD9 2.3.1 Heat-conduction equations
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations
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
CPS BEGIN bnswer
CPS 
CPS REPEAT this subroutine to obtain solution
CPS 
CPS EXITIF the maximum number of iterations is reached after setting...
CPS ... flag
CPS 
CPS   varchk - update variable states
CPS   outbnd - check if variables are out of bounds
CPS   
CPS EXITIF any variables are out of bounds
CPS 
CPS   IF this is a heat transfer only solution
CPS     gensl3 - obtain solution for heat transfer
CPS   ELSE
CPS   
CPS     IF this is a nonisothermal two phase calculation
CPS     
CPS       IF it is single porosity
CPS         gensl4 - obtain solution
CPS       ELSE it is dpdp problem
CPS         dpdp - obtain solution
CPS       ENDIF
CPS       
CPS     ENDIF
CPS 
CPS   
CPS     IF this is an isothermal two phase calculation
CPS     
CPS       IF it is single porosity
CPS         airctr - obtain solution
CPS       ELSE it is dpdp problem
CPS         dpdp - obtain solution
CPS       ENDIF
CPS       
CPS     ENDIF
CPS       
CPS     IF this is a solid-liquid-gas calculation
CPS     
CPS       IF it is single porosity
CPS         icectr - obtain solution
CPS       ELSE it is dpdp problem
CPS
CPS       ENDIF
CPS       
CPS     ENDIF
CPS 
CPS     IF this is a single phase calculation
CPS     
CPS       IF it is single porosity
CPS         gensl1 - obtain solution
CPS       ELSE it is dpdp problem
CPS         dpdp - obtain solution
CPS       ENDIF
CPS       
CPS     ENDIF
CPS     
CPS   ENDIF
CPS   
CPS EXITIF the residual is below the convergence criterion value
CPS 
CPS   Add 1 to the number of iterations and total iteration count
CPS 
CPS EXITIF the residual is below machine precision
CPS 
CPS   varchk - apply corrections to solution vectors
CPS 
CPS UNTIL the residual can go no lower
CPS 
CPS END bnswer
CPS 
C**********************************************************************

      use combi
      use comdi
      use comgi
      use comei
      use comdti
      use comai
      use comco2
      use comsi
      use comsplitts
      use davidi
      use comfem
      implicit none
c gaz 110919 moved iad_min to comai (global variable now)
      integer iad_mult, i
      real*8 fdum_mult
      parameter (fdum_mult=1.d02,iad_mult= 100)

      itert=0 
      minkt=0
      strd=1.
      iad=0
      mlz_save= 0
c gaz 110715
      mlz = 0
      if(g1.lt.0.0.and.jswitch.eq.0) then
	 iad_min = 0
      else if(maxit.ge.0) then
         iad_min=1
      else
         iad_min=2
      endif
      if(iflux_ts.ne.0.and.ico2.lt.0) then
         call cascade_sat(0)
      endif
      fdum = 0.0
      fdum_last = 1.

      if(istrs.gt.0) then
        delta_u = 0.0d0
        delta_v = 0.0d0
        delta_w = 0.0d0
      endif
c     
c     check iterations against maximum
c     
 1000 continue
      if(fdum.gt.fdum_last*fdum_mult.and.iad.ge.iad_mult) then
         mlz = -1
         bp = 0.0
         goto 2000
      else 
         fdum_last = fdum
      endif
      if(iad.gt.abs(maxit)) then
         mlz=-1
         goto 2000
      endif
      if(iad.gt.abs(maxit)/2) then
c     strd=0.95 
      endif
c     
c     update variable states
c   
      if(istrs_coupl.le.-99) then	 
c  istrs_coupl.le.-99 means only a stress solution (like calculating initial lithoststic stress)  
c (4,0) is nonlinear (if any) material properties)   
         call stressctr(4,0)        
      elseif(ice.eq.0) then

c     s kelkar 11feb2011 ihms=-15 is a special option
c     allow uncoupled porosity changes, ihms is reset to -3
c     ps calculated below is a factor so that 
c     pore volume = ps()*initial bulk volume (=sx1d)
c     this factor ps() is used in therm w to calculate deni()=ps()*roho
c     and the sx1d is multiplied in geneq1_.._..
         if(pore_factor.gt.0.0.and.ihms.eq.-3) then
            call porosity_wrt_displacements
         endif

         call varchk(0,0)
c the following is for fully coupled model         
         if(idof_stress.ge.4) then
	    call stressctr(4,0)
         endif
c     RJP 04/10/07 added 
         if(icarb.eq.1) then
            call icectrco2(-1,0)
            call icectrco2(1,0)
            call icectrco2(-34,0)
c     call icectrco2(22,0)
         endif
      else
         call icectr(-1,0)
         call icectr(1,0)
c     added this to check the state of hydrate for the next iteration
c     gaz-commented out 9-2-2003
         call icectr(-6,0)
c     if (mlz.eq.-1) goto 2000
      endif
c
c  very new call outbound here gaz 120709
c      call outbnd
c
c     call appropriate sub to generate equations

c     set up permeability variations with displacements (allocate memory)     
      if(istrs_coupl.gt.-99) then
         flag_pstrain_perm_coupling = 0
         if(flag_permmodel.eq.1) then
            if(flag_element_perm.eq.1) then
           ! finite element option
           ! Setup connectivity list of which elements each node belongs to
               
               if(iPlastic.eq.1) flag_pstrain_perm_coupling=1
               
               if(.not. allocated(NodeElems)) then
                  call Setup_NodeElems()
               endif                        
           ! Setup pointers to edge numbers
               if(.not. allocated(edgeNum1)) then
                  call setup_edgePointers_3D()
               endif
             ! allocate memory for permeability update if necessay
               call stress_perm(-1,0)
             ! update edge permeability factors for current state
               call update_permfactors()
            else if (iPlastic.eq.1) then
c     update permeabilities based on plastic strain, von Mises (if flag is called)
c     Hard-wired the flag, need to input
c     this is for sequential coupling
               flag_pstrain_perm_coupling = 1 
               if (flag_pstrain_perm_coupling.eq.1 .and. 
     &              ifem.eq.1) then
c     Setup connectivity list of which elements each node belongs to
                  if(.not. allocated(NodeElems)) then
                     call Setup_NodeElems()
                  endif
c     Setup pointers to edge numbers
                  if(.not. allocated(edgeNum1)) then
                     call setup_edgePointers_3D()
                  endif
                  
                  call stress_perm(-1,0)
                  
                  call update_permfactors()
c.............................................................
               else
c     cant do plasticity without ifem
                  write(iptty,*)'error. must have fem'
                  write(iptty,*)'for plasticity'
                  stop
               endif
            else
               
c     control volume approach 
c     allocate memory for permeability update if necessay
               call stress_perm(-1,0)
c     update nodal permeabilities (explicit)                  	
               call stress_perm(1,0)
c     deallocate memory for permeability update if necessay
               call stress_perm(-2,0)
            endif
         endif 
      endif

      if(idof_stress.ge.5) then
c 3d coupled THM
         call gensl_stress_coupled_3D
         do i = 1,neq
           delta_u(i) = delta_u(i) - bp(i + nrhs(1))
           delta_v(i) = delta_v(i) - bp(i + nrhs(2))
           delta_w(i) = delta_w(i) - bp(i + nrhs(3))
         enddo

      else if(istrs_coupl.le.-99) then	
c add stress equations and derivatives
c uncoupled stress equations
         call  stressctr(8,0)
         do i = 1,neq
           if(icnl.eq.0) then
            delta_u(i) = delta_u(i) - bp(i + nrhs(1))
            delta_v(i) = delta_v(i) - bp(i + nrhs(2))
            delta_w(i) = delta_w(i) - bp(i + nrhs(3))
           else 
c gaz 052222 bp only allocated 2*neq for 2 dof HM 
c stress needs bp by dimensions                
            delta_u(i) = delta_u(i) - bp(i + nrhs(1))
            delta_v(i) = delta_v(i) - bp(i + nrhs(2))               
           endif
         enddo

      else if(idoff.eq.-1) then
         call gensl3
      else
         if(ico2.gt.0) then
            if(idpdp.eq.0) then
               call gensl4
            else
               call dpdp(3)
            endif
         endif

         if(ico2.lt.0.and.ice.eq.0) then
c     single porosity
            if(idpdp.eq.0) then
               if(isplitts.ge.-1) then
                  call airctr(4,0)
               else if(isplitts.eq.-2) then
                  call airctr(-4,0)
c     explicit iterations finished in airctr so return
                  go to 2000
               endif
            else
c     dpdp enabled
               call dpdp(3)
            endif
         endif

c     gaz 10-15-2001
         if(ico2.lt.0.and.ice.ne.0) then
c     single porosity
            if(idpdp.eq.0) then
               call icectr(4,0)
            else
c     dpdp not enabled
c     
            endif
         endif

         if(ico2.eq.0) then
            if(idpdp.eq.0) then
c     RJP 04/17/07 changed following
               if(icarb.eq.1) then
                  call icectrco2(4,0)
               else
                  call gensl1
               end if
            else
               call dpdp(3)
            endif
         endif

      endif
      if(fdum.le.f0.and.iad.ge.iad_min) goto 2000
      if(mlz.lt.0) goto 2000
      iad=iad+1
      itotal=itotal+1
c     
c     apply nr corrections
c     
      if(idof_stress.ge.5) then
c     3d coupled THM 
         call stressctr(9,0)
      else if(istrs_coupl.le.-99) then
	 call stressctr(9,0) 
      else if(ice.eq.0) then
c     RJP 04/10/07 added following
         if(icarb.eq.1) then
            call icectrco2(2,0)
         else
            call varchk(1,0)
            if(idof_stress.ge.4) then
               call stressctr(9,0)            
            endif
         end if
      else
         call icectr(2,0)
      endif
      if(epe.le.0.) goto 2000
c     
c     check if varibles are out of bounds
c      outbnd moved higher
      call outbnd 
      if(mlz.ne.0) goto 2000
      goto 1000
 2000 continue
      end
