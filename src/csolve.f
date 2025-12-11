      subroutine csolve(hmon)
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
CD1 To loop through all solute time steps, computing the new
CD1 concentrations of each specie at each time.
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
CD2 04-07-93      B. Robinson   22      Revised to implement the
CD2                                     chemical reaction models
CD2                                     between solutes
CD2 03-08-94      B. Robinson   N/A     Revised to allow time steps to
CD2                                     be cut when solution is not
CD2                                     found after maximum number of
CD2                                     iterations
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/csolve.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.26   Wed May 29 15:07:40 1996   hend
CD2 Added variable diffusion with water content
CD2
CD2    Rev 1.25   Tue May 28 16:11:40 1996   hend
CD2 Fixed missing it10 assignment
CD2 
CD2    Rev 1.24   Fri May 24 09:53:56 1996   hend
CD2 Updated trac for mdnodes
CD2 
CD2    Rev 1.23   Tue May 21 10:37:34 1996   hend
CD2 Updated for use with both ldsp and dspl,dspv,or dspb
CD2 
CD2    Rev 1.22   Wed May 08 14:06:50 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.21   Wed May 08 13:34:22 1996   hend
CD2 Updated trac time step manipulations
CD2 
CD2    Rev 1.20   Fri May 03 14:19:48 1996   hend
CD2 Updated for GAZ mdnodes changes
CD2 Vel and Disp. set to 0 for parent connections
CD2 
CD2    Rev 1.19   Mon Apr 29 10:21:50 1996   hend
CD2 Updated to reflect GAZ changes -- all 3 sx components used
CD2 
CD2    Rev 1.18   Thu Apr 25 13:33:06 1996   hend
CD2 Updated for use in long/trans dispersion option
CD2 
CD2    Rev 1.17   Mon Mar 25 11:01:02 1996   hend
CD2 Fixed Indexing for Vapor Phase Velocities
CD2 
CD2    Rev 1.16   Mon Mar 25 09:26:46 1996   hend
CD2 Fixed missing increment in sehindex for
CD2 vapor phase
CD2 
CD2    Rev 1.15   Thu Mar 21 13:18:02 1996   hend
CD2 Fixed to use indexing independent of istrw
CD2 
CD2    Rev 1.14   Tue Mar 19 15:10:32 1996   hend
CD2 kz should be integer, not real
CD2 
CD2    Rev 1.13   Mon Mar 04 16:08:46 1996   hend
CD2 Removed uneccesary calculations from coneq1 and added trac input option
CD2 
CD2    Rev 1.12   Thu Feb 15 09:25:14 1996   zvd
CD2 Added requirement.
CD2    Rev 1.11   08/18/95 10:35:02   llt
CD2
CD2 rxn_interval already defined, removed for cray
CD2 
CD2    Rev 1.10   08/16/95 16:25:16   robinson
CD2 Corrected iteration problem and added write of solute mass balance data
CD2 
CD2    Rev 1.9   08/09/95 15:52:30   zvd
CD2 Corrected write to iatty when unassigned.
CD2 
CD2    Rev 1.8   04/25/95 09:19:12   llt
CD2 retrieved lost log history informatin
CD2 
CD2    Rev 1.7   04/21/95 14:43:52   robinson
CD2 Corrected mechanism for cutting solute time step
CD2 
CD2    Rev 1.6   04/03/95 08:43:58   robinson
CD2 Correction to solute mass balance calculation
CD2 
CD2    Rev 1.5   01/28/95 14:20:08   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.4   06/20/94 11:14:24   zvd
CD2 Added ierr unit number for error output.
CD2
CD2    Rev 1.3   03/28/94 16:36:14   robinson
CD2 Removed unneeded array.
CD2
CD2    Rev 1.2   03/23/94 14:46:52   robinson
CD2 Corrected procedure for cutting solute time step
CD2
CD2    Rev 1.1   03/18/94 16:15:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2
CD2    Rev 1.0   01/20/94 10:22:40   pvcs
CD2 original version in process of being certified
CD2
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
CD4 max_species, residual, diagonal_element, ipbn, nbd, lenreal,
CD4 ipt1, ipt2, ipt3, ipt4, ipt5, ipt6, ipt7, ipt8, ipt9, ipt10,
CD4 ipt11, ipt12, ipt5v, days, daysi, nts, npt, neq, cnsk, t2sk, rc,
CD4 datcmm, daycmx, daycs, iout, t1sk, dtotc, nsp, iaddt, iadd, idpdp,
CD4 idualp, icns, anl, anv, an, danv, danl, rolf, stored_derivative,
CD4 stored_residual, a_henry, dh_henry, gas_const, t, temp_conv, phi,
CD4 avgmolwt, mw_water, cm, n, volume, qcin, qcout, denci, dench, 
CD4 dencj, anlo, dtotdm, iaccmx, daycm, rcss, qcin, qcout, qcrxn
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
CD4 daysi        real*8  faar    Time in days used to determine time
CD4                              step
CD4 days         real*8  faar    Time in days used to determine time
CD4                              step
CD4 day          real*8  faar    Time step in days
CD4 daycmm       real*8  faar    Minimum time step in days
CD4 daycmx       real*8  faar    Maximum time step in days
CD4 dtotc        real*8  faar    Current time step size for tracer
CD4                              solution
CD4 dtotdm       real*8  faar    Heat and mass transfer time step
CD4 daycm        real*8  faar    Solute time step multiplier
CD4 
CD4 
CD4 
CD4 nsp          int     fdd1i   Current solute being computed
CD4 nspeci       int     fdd1i   Total number of solutes being simulated
CD4 iaddt        int     fdd1i   Total number of iterations for each
CD4                              solute
CD4 npn          int     fdd1i   Index directing code to the correct
CD4                              position in solute arrays
CD4 npt          int     fdd1i   Array of starting positions in
CD4                              the concentration array for
CD4                              each specie
CD4 iadd         int     fdd1i   Number of iterations for each solute
CD4                              at the current time step
CD4 
CD4 
CD4 cnsk         real*8  fdd1    Inlet concentration at each node due to
CD4                              source term if it exists
CD4 t1sk         real*8  fdd1    Initial time of tracer injection for
CD4                              each tracer
CD4 t2sk         real*8  fdd1    Final time of tracer injection for
CD4                              each tracer
CD4 rc           real*8  fdd1    Source/sink and chemical reaction terms
CD4                              of the tracer mass balance
CD4 cm           real*8  fdd1    Total solute mass stored for each
CD4                              solute
CD4 dench        real*8  fdd1    Solute mass storage at previous time
CD4                              for each node
CD4 anl          real*8  fdd1    Current concentration at each node
CD4                              for the current tracer
CD4 anlo         real*8  fdd1    Current concentration at each node
CD4                              for the current tracer
CD4 dencj        real*8  fdd1    Solute mass storage at current time
CD4                              for each node
CD4 
CD4 
CD4 
CD4 
CD4 idualp       int     faai    Flag denoting if dual porosity
CD4                              solution is being computed
CD4 iaccmx       int     faai    Maximum number of iterations for
CD4                              which the time step is to be increased
CD4 
CD4 
CD4 
CD4 
CD4 volume       real*8  fdd     Total volume associated with each node
CD4 sk           real*8  fdd     Source/sink flow rate array
CD4 
CD4 
CD4 
CD4 
CD4 denci        real*8  fcc     Solute mass storage at current time
CD4                              for each node
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 cnswer          N/A      Computes the new concentrations at this
CD4                          time step for a given solute
CD4 plotc1          N/A      Writes the solute information to the .trc
CD4                          file
CD4 resettrc        N/A      Resets solute concentration terms, cuts
CD4                          time step
CD4 tyming          real*8   Computes cpu time for a portion of the run
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
CD5 new_concentration
CD5              real*8      Array of new concentrations for each
CD5                             solute at a node
CD5 time1        real*8      Time parameter used in time step
CD5                             determination
CD5 time2        real*8      Time parameter used in time step
CD5                             determination
CD5 tajj         real*8      Cpu time parameter
CD5 caz          real*8      Cpu time parameter
CD5 n_node_sets  int         Number of set of nodes in simulation
CD5 i_node_set   int         Current node set
CD5 ispecies     int         Current species number
CD5 mim          int         Current node number
CD5 mi           int         Current concentration unknown number
CD5 residual     real*8      Residual for concentration of each
CD5                             species at a node
CD5 ndummy       int         Parameter for determining space in
CD5                             concentration array
CD5 dvap_conc    real*8      Derivative of vapor concentration with
CD5                             respect to liquid concentration
CD5 dayst        real*8      Time in days used to determine time step
CD5 daytr        real*8      Time step in days
CD5 nts          int         Number of solute time steps taken
CD5 icfin        int         Flag denoting whether this is the last
CD5                          time step
CD5 id           int         Do loop index parameter
CD5 i            int         Index of current concentration unknown
CD5 rcd          real*8      Term used in setting source/sink terms
CD5 daytrm       real*8      Time step used in solute calculation
CD5 ja           int         Index used in loop over all nodes
CD5 jad          int         Index used in loop over all nodes
CD5 vdum         real*8      Total volume of the current node
CD5 dencht       real*8      Solute mass at previous time step for the
CD5                          current node
CD5 tol_value     int         For each specie, if 0, only a tolerance
CD5                          check calculation was made in cnswer, if
CD5                             1 new concentrations were computed
CD5 tolerance_check
CD5              int         If 0, tolerance checks are not needed, if
CD5                             1, they are needed
CD5 iapos        int         Array of indexes for the A matrix
CD5 ibpos        int         Array of indexes for the right hand
CD5                             side
CD5 irow         int         Do loop index over all rows of matrix
CD5 icolumn      int         Do loop index over all columns of matrix
CD5 reset_tracer logical     Flag denoting if the tracer values have
CD5                          been reset and the time step cut
CD5 iadd_max     int         The maximum number of Newton iteration
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
CD6   Initialize time parameters, initialize number of iterations for
CD6   each solute, looping over each solute.
CD6   
CD6   In an outer loop over all times:
CD6   
CD6     A check is performed to see if the time step needs to be
CD6     altered to exactly hit the final time, and the time step is
CD6     then taken.
CD6     
CD6     Then, in a loop over each solute, a loop over each node is
CD6     carried out to perform the following steps:
CD6     
CD6       The code checks to see if this node is a source of tracer,
CD6       and if it is, the code determines whether the time step
CD6       needs to be reset based on the times that injection starts
CD6       and stops for this solute and node.
CD6       
CD6       At the end of this loop, the time step has been determined.
CD6       
CD6       Then, a loop for each solute is used to set the convergence
CD6       parameter.
CD6     
CD6     Next, a loop is executed until all concentrations are
CD6     converged, by doing the following:
CD6     
CD6       Loop through each solute, calling cnswer to compute the new
CD6       concentrations at this time step.
CD6       
CD6       Also in this loop on each solute, adjust the convergence
CD6       flag, add to the number of iterations, and call thermc to
CD6       compute new storage and reaction terms (for dual porosity,
CD6       call thermc three times, once for each set of nodes).
CD6       
CD6     We then exit if all solutes have successfully converged and a
CD6     second check on convergence shows that ocnvergence is still
CD6     achieved.
CD6     
CD6     Then, after convergence at the current time:
CD6     
CD6       Loop over all solutes and all nodes, computing the overall
CD6       mass balance terms.
CD6       
CD6       Call plotc1 for each solute to write plot information to file.
CD6       Exit the entire outer loop if the final time has been reached.
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
CD9 2.3.4 Solute-transport equations
CD9 2.5.1 Implement time-step mechanism
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
CPS BEGIN csolve
CPS 
CPS Set times and time step
CPS Initialize number of time steps
CPS 
c changed 8/3/94 to incorporate coupling
CPS FOR each grouping
CPS   Initialize the total number of iterations to zero
CPS ENDFOR
c end change
CPS 
CPS LOOP through all time steps
CPS 
CPS   Increase current time
CPS   Set flag indicating the current status of time
CPS   
CPS   IF the time is greater than the final time is supposed to be
CPS     Reset time and time step to hit final time exactly
CPS     Set flag to indicate this is the final time step
CPS   ENDIF
CPS   
CPS   Take actual time step
CPS   
CPS   FOR each aqueous component
CPS     FOR each node of this tracer
CPS   
CPS       IF there is a tracer source at this node for this tracer
CPS         IF the initial time of tracer injection is positive
CPS           IF we have passed this initial time
CPS             IF the initial time of tracer injection is equal to ...
CPS             ... the time the tracer solution is enabled
CPS                 Reset time step to minimum 
CPS             ELSE 
CPS               Reset the time step
CPS               IF this time step is less than the minimum time step
CPS                 Reset time step to minimum
CPS                 Reset the initial time of tracer injection
CPS               ENDIF
CPS             ENDIF
CPS             Write tracer injection information to file
CPS             Set flag indicating current status of time
CPS             Set initial tracer injection time negative so we do...
CPS             ... not turn tracer on again later
CPS           ENDIF
CPS         ENDIF
CPS       
CPS         IF the final time of tracer injection is positive
CPS           IF we have passed this final time
CPS             Reset the time step
CPS             IF this time step is less than the minimum time step
CPS               Reset time step to minimum
CPS               Reset the final time of tracer injection
CPS             ENDIF
CPS             Write tracer injection information to file
CPS             Set flag indicating current status of time
CPS             Set final tracer injection time negative so we do...
CPS             ... not turn tracer on again later
CPS           ENDIF
CPS         ENDIF
CPS       ENDIF
CPS     
CPS     ENDFOR each node of this tracer
CPS   ENDFOR each aqueous component
CPS   
CPS   Set time step in seconds and other time parameters
CPS   
c changed to incorporate components 9/20/95
CPS   Set parameter to indicate that concentrations must be solved for
c end change
CPS   Set flag to indicate tracer values have not yet been reset
CPS   
CPS   Set parameter to indicate that tolerance checks are not yet...
CPS   ... determined to be necessary
CPS     
c changed to incorporate components 9/28/95
CPS    
CPS   cnswer - compute all new concentrations at this time step
CPS       
CPS   IF a new set of concentrations was just solved for
CPS     IF convergence was not achieved
CPS       resettrc - reset values, cut solute time step
CPS       Set flag to denote tracer values have been reset
CPS     ENDIF
CPS     Set parameter indicating tolerance check still ...
CPS     ... must be made
CPS   ENDIF
CPS       
CPS   EXITIF convergence was not achieved
CPS       
CPS       Add to running total of SIA iterations
CPS   
CPS       IF the storage and rate routine needs to be called
CPS         thermc - compute storage and rate terms
CPS         IF this is a dpdp solution
CPS           thermc - compute storage terms for second set of nodes
CPS         
CPS       ELSEIF this is a dual porosity solution
CPS           thermc - compute storage terms for second set of nodes
CPS           thermc - compute storage terms for third set of nodes
CPS         ENDIF
CPS       ENDIF
CPS       
CPS       
c end change
CPS     
CPS   IF the values have not been reset
CPS   EXITIF all tolerance checks have passed
CPS   
CPS       IF this was the maximum number of outer iterations to be taken
CPS         resettrc - reset values, cut solute time step
CPS         Set flag denoting values have been reset
CPS       ENDIF
CPS     ENDIF
CPS   EXITIF maximum number of outer iterations was taken
CPS   ENDLOOP
CPS   
CPS   IF convergence was achieved
CPS     FOR each reaction
CPS       IF the current reaction is an equilibrium reaction THEN
CPS         IF the reactions chosen to be in equilibrium are not ...
CPS         ... in equilibrium THEN
CPS           IF simple multiplier is used to adjust to forward ...
CPS           ... rate constant THEN
CPS              Multiply current forward rate constant by multiplier
CPS           ELSE 
CPS             use more complicated formula to adjust forward ...
CPS             ... rate constant
CPS           ENDIF
CPS           Reset outer iteration counter
CPS           Set flag indicating that another outer iteration ...
CPS           ... must be set
CPS           FOR each node
CPS             Calculate reverse rate constant from forward rate ...
CPS             ... constant and equilibrium constant
CPS           ENDFOR
CPS         ENDIF
CPS       ELSE 
CPS         the current reaction is a kinetic reaction
CPS       ENDIF 
CPS     ENDFOR each reaction
CPS     IF another outer iteration must be performed with higher ...
CPS     ... rate constants to achieve equilibrium
CPS       Add to the equilibrium iteration counter
CPS       Write equilibrium iteration number
CPS       Return to perform another outer iteration
CPS     ENDIF
CPS     Reset equilibrium iteration counter to zero
CPS     Write convergence message to screen
CPS     
CPS     FOR each species
CPS     
CPS       IF this solute is a Henry's Law species with vapor...
CPS       ... concentration specified
CPS     
CPS         IF this is a dpdp simulation
CPS           Set number of node sets to 2
CPS         ELSEIF it is dual porosity
CPS           Set number of node sets to 3
CPS         ELSE
CPS           Set number of node sets to 1
CPS         ENDIF
CPS   
CPS         FOR each node set
CPS           FOR each node
CPS             Compute liquid concentration corresponding to the...
CPS             ... input vapor concentration
CPS           ENDFOR each node
CPS         ENDFOR each node set
CPS       ENDIF
CPS     ENDFOR
CPS   
CPS     FOR each tracer species
CPS   
CPS       FOR each tracer unknown
CPS         Compute tracer mass balance terms
CPS       ENDFOR
CPS       bcon - handle boundary conditions for special case
CPS   
CPS       plotc1 - call plot output
CPS     ENDFOR each tracer species
CPS     Increase number of time steps by 1
CPS     
CPS   ENDIF
CPS   
CPS   IF we are not resetting the tracer
CPS     FOR each group
CPS       Find the maximum number of Newton iterations taken ...
CPS       ... in each grouping
CPS     ENDFOR each group
CPS     IF we are leaving csolve
CPS       Compute average time stepp
CPS       EXIT
CPS     ENDIF
CPS     Adjust the time step
CPS   ENDIF
CPS   
CPS ENDLOOP
CPS
CPS END CSOLVE
CPS 
C**********************************************************************
c****--------------------------------------------------------------****c
c**** set time step for tracer solution , call solution            ****c

      use combi
      use comchem
      use comci
      use comdi
      use comgi
      use comji
      use comrxni
      use comrxnb
      use comcouple
      use comflow
      use comdti
      use comai
      use comuserc, only : usroption
      use davidi
      implicit none

      integer tol_value, nfinal
      real*8 time1
      real*8 time2
      integer iter_counter
      real*8 tajj
      real*4 caz(2)
      real*8 dayst
      real*8 daytr
      integer nts
      integer icfin
      integer id
      integer i
      real*8 rcd
      real*8 daytrm
      integer ja
      integer jad
      real*8 vdum, dumv
      real*8 dencht
      real*8 tyming
      logical reset_tracer
      integer hmon
      real*8 tempx,tempy,tempz,dis,area_t,templength
      real*8 dilfi,dilfkb,axyf,axy,divfi,radi,radkb
      real*8 prodrolf,prodrovf,sx2c,sx3c,sxzc,toldil
      integer j,jj,nmatavw,icd
      integer nmatadd,ii,ii1,ii2,iq,jmi,jm,ij1,ij2,ij,iz,kb,kz
      integer neighc,neqp1,iau,idg,jmia,sehindexl,sehindexv
      real*8 fid,fid1,vxy,divfkb,vxyf,vmag,sehmindays,thetav
      real*8 cord1x,cord1y,cord2x,cord2y,cord2xp,cord2yp
      real*8 cord1z,cord2z,newx,newy,newz,dispzavw,dispyavw,dispxavw
      real*8 newdiff, concadiff, satr, ptime
      real*8 :: last_time = 0.
      integer sia_iter, sia_iter_tot_old
      integer :: sia_iter_tot = 0
      integer isolute
      integer ic
      integer im
      integer iv
      integer ix
      integer :: iprttrc = 0, last_step = 0
      integer :: idebug = 0
      logical :: time2print, istop_flag
      parameter(toldil = 1.d-20)
      real*8 ps_min, ps_max, s_min_salt 
      parameter(ps_min = 1.d-05 ,ps_max = 0.9999, s_min_salt = 1.e-4)
      save daytr, iprttrc, icfin, last_step, last_time
c seh
c set velocities here instead of in coneq1
c gaz 071016
      sia_iter_tot_old = iaddt(1)
      if (daytr .le. 0) daytr = daycmm
      if (ianpe.ne.0) then
         continue
      else if ((ldsp.eq.1).and.(dispsame.eq.0)) then
         continue
      else if (ldsp.eq.1) then
         if ((hmon.eq.1).or.(sehdonevel.eq.0)) then
            if (hmon.eq.0) then
               sehdonevel=1
            else
               sehdonevel=0
            endif
            sehindexl=1
            sehindexv=1
            neqp1=neq+1

            if (hvliquid.eq.1) then
               if (icnl.eq.0) then
                  if(idualp.eq.0) then
                     nfinal = n0
                  else
                     nfinal = neq
                     do i = neq+1, 3*neq
c----------------- phs 9/26/2001 - added s(i) correction
                        if (irdof .ne. 13 .or. ifree .ne. 0) then
                           displx(i) = sehdiff(itrcdsp(i))*ps(i)*s(i)
                        else
                           displx(i) = sehdiff(itrcdsp(i))*ps(i)
                        end if
                     end do
                  end if
                  do i=1,nfinal
                     vmag=sqrt(pnx(n+i)*pnx(n+i)+
     +                    pny(n+i)*pny(n+i)+pnz(n+i)*pnz(n+i))
                     if (irdof .ne. 13 .or. ifree .ne. 0) then
                        satr = s(i)
                     else
                        satr = 1.0d0
                     end if
                     if (vmag.eq.0.) vmag=1e-30
                     if(i.gt.neq.and.idualp.eq.0) then
                        icd=neq
                     else
                        icd=0
                     endif
                     iz=i-icd
                     ii1=nelm(iz)+1
                     ii2=nelm(iz+1)
                     jmi=nelmdg(iz)
                     do jm=jmi+1,ii2
c     compute alpha
c check for multiply defined nodes -- don't compute terms for here
                        iw=istrw(jm-neqp1)
                        if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz).ne.0.)
     2                       then
                           kz=nelm(jm)
                           call rotate(cord(iz,1),cord(iz,2),cord(iz,3),
     &                          pnx(n+i)+cord(iz,1),pny(n+i)+cord(iz,2),
     &                          pnz(n+i)+cord(iz,3),cord(kz,1),
     &                          cord(kz,2),cord(kz,3),newx,newy,newz)
                           cord1x=0.
                           cord1y=0.
                           cord1z=0.
                           cord2x=newx
                           cord2y=newy
                           cord2z=newz
                           ii=itrcdsp(i)
                           dispxavw=tcly(1,ii)
                           dispyavw=tcly(1,ii)
                           dispzavw=tclx(1,ii)
                           tempx = (cord1x-cord2x)**2
                           tempy = (cord1y-cord2y)**2
                           tempz = (cord1z-cord2z)**2
                           templength = tempx+tempy+tempz
                           alphaconl(sehindexl) =
     +                          sqrt((dispxavw*dispyavw*dispzavw*
     +                          templength)/(dispyavw*dispzavw*tempx +
     3                          dispxavw*dispzavw*tempy +
     4                          dispxavw*dispyavw*tempz))
c     compute entire term
                           iw=istrw(jm-neqp1)
                           sx2c=sx(iw,isox)
                           sx3c=sx(iw,isoy)
                           sxzc=sx(iw,isoz)
                           newdiff = concadiff(1,mflagl(1,itrcdsp(i)),
     &                          sehdiff(itrcdsp(i)),ps(i),satr,
     &                          phi(i),t(i))
                           sehvell(sehindexl)=(sx2c+sx3c+sxzc)*
     &                          (vmag*alphaconl(sehindexl)+
     &                          newdiff*ps(i)*satr) 
c------------------ added s(i) to the newdiff 
                           sehindexl=sehindexl+1
                        endif
                     enddo
                  enddo
               else
                  if(idualp.eq.0) then
                     nfinal = n0
                  else
                     nfinal = neq
                     do i = neq+1, 3*neq
c------------------ added s(i) to displx  liquid  
                        if (irdof .ne. 13 .or. ifree .ne. 0) then
                           displx(i) = sehdiff(itrcdsp(i))*ps(i)*s(i)
                        else
                           displx(i) = sehdiff(itrcdsp(i))*ps(i)
                        end if
                     end do
                  end if
                  do i=1,nfinal
                     vmag=sqrt(pnx(n+i)*pnx(n+i)+pny(n+i)*pny(n+i))
                     if (irdof .ne. 13 .or. ifree .ne. 0) then
                        satr = s(i)
                     else
                        satr = 1.0d0
                     end if
                     if (vmag.eq.0.) vmag=1e-30
                     thetav=acos(pnx(n+i)/vmag)
                     if(i.gt.neq.and.idualp.eq.0) then
                        icd=neq
                     else
                        icd=0
                     endif
                     iz=i-icd
                     radi=cord(iz,3)
                     ii1=nelm(iz)+1
                     ii2=nelm(iz+1)
                     jmi=nelmdg(iz)
                     do jm=jmi+1,ii2
c compute alpha
                        iw=istrw(jm-neqp1)
                        if (sx(iw,isox)+sx(iw,isoy).ne.0.) then
                           kz=nelm(jm)
                           cord1x=0.
                           cord1y=0.
                           cord2xp=cord(kz,1)-cord(iz,1)
                           cord2yp=cord(kz,2)-cord(iz,2)
                           cord2x=
     &                          cos(thetav)*cord2xp-sin(thetav)*cord2yp
                           cord2y=
     &                          sin(thetav)*cord2xp+cos(thetav)*cord2yp
                           tempx = (cord1x-cord2x)**2
                           tempy = (cord1y-cord2y)**2
                           ii=itrcdsp(i)
                           alphaconl(sehindexl)=
     &                          sqrt((tclx(1,ii)*tcly(1,ii)*
     &                          (tempx+tempy))/
     &                          (tclx(1,ii)*tempy + tcly(1,ii)*tempx ))
c     compute whole term                        
                           iw=istrw(jm-neqp1)
                           radkb=0.5*(radi+cord(kz,3))
                           sx2c=radkb*sx(iw,isox)
                           sx3c=radkb*sx(iw,isoy)
                           vmag=sqrt(pnx(n+i)*pnx(n+i)+pny(n+i)*
     &                             pny(n+i))
                           newdiff = concadiff(1,mflagl(1,itrcdsp(i)),
     &                          sehdiff(itrcdsp(i)),ps(i),satr,
     &                          phi(i),t(i))
                           sehvell(sehindexl)=(sx2c+sx3c)*
     &                          (vmag*alphaconl(sehindexl)+
     &                          satr*newdiff*ps(i))
c----------------- phs 9/26/2001 - added s(i) to the newdiff 
                           sehindexl=sehindexl+1
                        endif
                     enddo
                  enddo
               endif
            endif
            
            if (hvvapor.eq.1) then
               if(icnl.eq.0) then
                  if(idualp.eq.0) then
                     nfinal = n0
                  else
                     nfinal = neq
                     do i = neq+1, 3*neq
                        dispvx(i) = sehdiff(itrcdsp(i))*ps(i)*(1-s(i))
c----------------- phs 9/26/2001 - added (1-s(i)) correction 
                     end do
                  end if
                  do i=1,nfinal
                     vmag=sqrt(pnx(n*2+i)*pnx(n*2+i)+
     +                    pny(n*2+i)*pny(n*2+i)+pnz(n*2+i)*pnz(n*2+i))
                     if (vmag.eq.0.) vmag=1e-30
                     if(i.gt.neq.and.idualp.eq.0) then
                        icd=neq
                     else
                        icd=0
                     endif
                     iz=i-icd
                     ii1=nelm(iz)+1
                     ii2=nelm(iz+1)
                     jmi=nelmdg(iz)
                     do jm=jmi+1,ii2
c     compute alpha
                        iw=istrw(jm-neqp1)
                        if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz).ne.0.)
     2                       then
                           kz=nelm(jm)
                           call rotate(cord(iz,1),cord(iz,2),cord(iz,3),
     &                          pnx(n*2+i)+cord(iz,1),pny(n*2+i)+
     &                          cord(iz,2),pnz(n*2+i)+cord(iz,3),
     &                          cord(kz,1),cord(kz,2),
     &                          cord(kz,3),newx,newy,newz)
                           cord1x=0.
                           cord1y=0.
                           cord1z=0.
                           cord2x=newx
                           cord2y=newy
                           cord2z=newz
                           ii=itrcdsp(i)
                           dispxavw=tcvy(1,ii)
                           dispyavw=tcvy(1,ii)
                           dispzavw=tcvx(1,ii)
                           tempx = (cord1x-cord2x)**2
                           tempy = (cord1y-cord2y)**2
                           tempz = (cord1z-cord2z)**2
                           templength = tempx+tempy+tempz
                           alphaconv(sehindexv) =
     +                          sqrt((dispxavw*dispyavw*dispzavw*
     +                          templength)/(dispyavw*dispzavw*tempx +
     3                          dispxavw*dispzavw*tempy +
     4                          dispxavw*dispyavw*tempz))
c     compute entire tem
                           iw=istrw(jm-neqp1)
                           sx2c=sx(iw,isox)
                           sx3c=sx(iw,isoy)
                           sxzc=sx(iw,isoz)
                           newdiff = concadiff(2,mflagv(1,itrcdsp(i)),
     &                          sehdiffv(itrcdsp(i)),ps(i),s(i),
     &                          phi(i),t(i))
                           sehvell(sehindexv)=(sx2c+sx3c+sxzc)*
     &                          (vmag*alphaconv(sehindexv)+
     &                          newdiff*ps(i)*(1-s(i)))
c----------------- phs 9/26/2001 - added (1-s(i)) correction
                           sehindexv=sehindexv+1
                        endif
                     enddo
                  enddo
               else
                  if(idualp.eq.0) then
                     nfinal = n0
                  else
                     nfinal = neq
                     do i = neq+1, 3*neq
                        dispvx(i) = sehdiff(itrcdsp(i))*ps(i)*(1-s(i))
c----------------- phs 9/26/2001 - added (1-s(i)) correction
                     end do
                  end if
                  do i=1,nfinal
                     vmag=sqrt(pnx(n*2+i)*pnx(n*2+i)+
     &                    pny(n*2+i)*pny(n*2+i))
                     if (vmag.eq.0.) vmag=1e-30
                     thetav=acos(pnx(n*2+i)/vmag)
                     if(i.gt.neq.and.idualp.eq.0) then
                        icd=neq
                     else
                        icd=0
                     endif
                     iz=i-icd
                     radi=cord(iz,3)
                     ii1=nelm(iz)+1
                     ii2=nelm(iz+1)
                     jmi=nelmdg(iz)
                     do jm=jmi+1,ii2
c compute alpha
                        iw=istrw(jm-neqp1)
                        if (sx(iw,isox)+sx(iw,isoy).ne.0.) then
                           kz=nelm(jm)
                           cord1x=0.
                           cord1y=0.
                           cord2xp=cord(kz,1)-cord(iz,1)
                           cord2yp=cord(kz,2)-cord(iz,2)
                           cord2x=
     &                          cos(thetav)*cord2xp-sin(thetav)*cord2yp
                           cord2y=
     &                          sin(thetav)*cord2xp+cos(thetav)*cord2yp
                           tempx = (cord1x-cord2x)**2
                           tempy = (cord1y-cord2y)**2
                           ii=itrcdsp(i)
                           alphaconl(sehindexv)=
     &                          sqrt((tcvx(1,ii)*tcvy(1,ii)*
     &                          (tempx+tempy))/
     &                          (tcvx(1,ii)*tempy + tcvy(1,ii)*tempx ))
c     compute whole term                        
                           iw=istrw(jm-neqp1)
                           radkb=0.5*(radi+cord(kz,3))
                           sx2c=radkb*sx(iw,isox)
                           sx3c=radkb*sx(iw,isoy)
                           newdiff = concadiff(2,mflagv(1,itrcdsp(i)),
     &                          sehdiffv(itrcdsp(i)),ps(i),s(i),
     &                          phi(i),t(i))
                           sehvelv(sehindexv)=(sx2c+sx3c)*
     &                          (vmag*alphaconv(sehindexv)+
     &                          newdiff*ps(i)*(1-s(i)))
c----------------- phs 9/26/2001 - added (1-s(i)) correction
                           sehindexv=sehindexv+1
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
      else
c ldsp=0, dispsame either 0 or 1
         if ((hmon.eq.1).or.(sehdonevel.eq.0)) then
            if (hmon.eq.0) then
               sehdonevel=1
            else
               sehdonevel=0
            endif
            sehindexl=1
            sehindexv=1
            neqp1=neq+1
            nmatavw = nelm(neqp1)-neqp1
c gaz 100917            
            if(gdkm_flag.ne.0) then
              nfinal = neq
            else if(idualp.eq.0) then
               nfinal = n0
            else
               nfinal = neq
               do i = neq+1, 3*neq
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     displx(i) = diffmfl(1,itrc(i))*ps(i)*s(i)
c----------------- phs 9/26/2001 - added s(i) correction
                  else
                     displx(i) = diffmfl(1,itrc(i))*ps(i)
                  end if
               end do
            end if
               
            do i=1,nfinal
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  satr = s(i)
               else
                  satr = 1.0d0
               end if
               if(i.gt.neq.and.idualp.eq.0) then
                  icd=neq
                  nmatadd = nmatavw
               else
                  icd=0
                  nmatadd = 0
               endif
               iz=i-icd
               
               ii1=nelm(iz)+1
               ii2=nelm(iz+1)
               idg=nelmdg(iz)-ii1+1
               iq=0
               jmi=nelmdg(iz)
               jmia=jmi-neqp1
               do jm=jmi+1,ii2
                  iq=iq+1
                  iw=istrw(jm-neqp1)
                  if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz).eq.
     &                 0.) then
                     iq=iq-1
                  else
                     it8(iq)=nelm(jm)+icd
                     it9(iq)=jm-ii1+1
                     it10(iq)=iw
                     it11(iq)=jm-neqp1
                     ij1=nelm(nelm(jm))+1
                     ij2=nelmdg(nelm(jm))-1
                     do ij=ij1,ij2
                        if(nelm(ij).eq.iz) then
                           it12(iq)=ij-neqp1
                        endif
                     enddo
                  endif
               enddo
               
c     liquid velocities        
               if (hvliquid.eq.1) then
                  do jm = 1, iq
                     kb = it8(jm)
                     iau=it11(jm)
                     neighc=it9(jm)
                     t8(neighc) = a_axy(iau+nmatadd)
                     fid = 0.5
                     if(t8(neighc).lt.0.0) then
                        fid=dnwgta
                        fid1 = 1.-fid
                        t7(neighc) = fid*dil(kb)
                        t9(neighc) = fid1*dil(i)
                     elseif(t8(neighc).gt.0.) then
                        fid=upwgta
                        fid1 = 1.-fid
                        t7(neighc) = fid*dil(kb)
                        t9(neighc) = fid1*dil(i)
                     else
                        t7(neighc) = toldil
                        t9(neighc) = toldil
                     end if
                  enddo
                  
                  if (icnl.eq.0) then
                     do jm=1,iq
                        kb=it8(jm)
                        kz=kb-icd
                        neighc=it9(jm)
                        iw=it10(jm)
                        sx2c=sx(iw,isox)
                        sx3c=sx(iw,isoy)
                        sxzc=sx(iw,isoz)
                        jj = it11(jm) 
                        j = kb
                        axy = t8(neighc)
                        
                        tempx = (cord(iz,1)-cord(kz,1))**2
                        tempy = (cord(iz,2)-cord(kz,2))**2
                        tempz = (cord(iz,3)-cord(kz,3))**2
                        dis=sqrt(tempx + tempy + tempz)
                        area_t=-(sx2c+sx3c+sxzc)*dis
                        dilfi=t9(neighc)
                        dilfkb=t7(neighc)
                        axyf=max(toldil,(dilfkb+dilfi))
                        prodrolf = rolf(i) * rolf(kb)
                        sehvell(sehindexl)=0
                        if(axy.ne.0.) then
                           if(prodrolf.gt.0.) then
                              if (area_t.gt.0.) then
                                 sehvell(sehindexl)=abs(axy/axyf/area_t)
     &                                *((rolf(kb)*dilfi+rolf(i)*dilfkb)/
     +                                prodrolf)
                              endif
                           endif
                        endif
                        if (dispsame.eq.1) then
                           newdiff =  concadiff(1,mflagl(1,itrcdsp(i)),
     &                          sehdiff(itrcdsp(i)),ps(i),satr,
     &                          phi(i),t(i))
                           sehvell(sehindexl)=(sx2c+sx3c+sxzc)*
     &                          (sehvell(sehindexl)*
     +                          alphaconl(sehindexl)+newdiff*ps(i)*satr)
c----------------- phs 9/26/2001 - added s(i) correction
                        endif
                        sehindexl=sehindexl+1
                     enddo
                  else
                     radi=cord(iz,3)
                     do jm=1,iq
                        kb=it8(jm)
                        kz=kb-icd
                        neighc=it9(jm)
                        iw=it10(jm)
                        radkb=0.5*(radi+cord(kz,3))
                        sx2c=radkb*sx(iw,isox)
                        sx3c=radkb*sx(iw,isoy)
                        jj = it11(jm) 
                        j = kb
                        axy = t8(neighc)
                        
                        tempx = (cord(iz,1)-cord(kz,1))**2
                        tempy = (cord(iz,2)-cord(kz,2))**2
                        dis=sqrt(tempx + tempy)
                        area_t=-(sx2c+sx3c)*dis*
     &                       0.5*(cord(iz,3)+cord(kz,3))
                        dilfi=t9(neighc)
                        dilfkb=t7(neighc)
                        axyf=max(toldil,(dilfkb+dilfi))
                        prodrolf = rolf(i) * rolf(kb)
                        sehvell(sehindexl)=0
                        if(axy.ne.0.) then
                           if(prodrolf.gt.0.) then
                              if(area_t.gt.0.) then
                                 sehvell(sehindexl)=abs(axy/axyf/area_t)
     &                                *((rolf(kb)*dilfi+rolf(i)*dilfkb)/
     +                                prodrolf)
                              endif
                           endif
                        endif
                        if (dispsame.eq.1) then
                           newdiff = concadiff(1,mflagl(1,itrcdsp(i)),
     &                          sehdiff(itrcdsp(i)),ps(i),satr,
     &                          phi(i),t(i))
                           sehvell(sehindexl)=(sx2c+sx3c)*
     &                          (sehvell(sehindexl)*
     +                          alphaconl(sehindexl)+newdiff*ps(i)*satr)
c----------------- phs 9/26/2001 - added s(i) correction
                        endif
                        sehindexl=sehindexl+1
                     enddo
                  endif
               endif
               
c     for vapor
               if (hvvapor.eq.1) then
                  do jm = 1, iq
                     kb = it8(jm)
                     iau=it11(jm)
                     neighc=it9(jm)
                     t8(neighc) = a_vxy(iau+nmatadd)
                     fid = 0.5
                     if(t8(neighc).lt.0.0) then
                        fid=dnwgta
                        fid1 = 1.-fid
                        t7(neighc) = fid*div(kb)
                        t9(neighc) = fid1*div(i)
                     elseif(t8(neighc).gt.0.) then
                        fid=upwgta
                        fid1 = 1.-fid
                        t7(neighc) = fid*div(kb)
                        t9(neighc) = fid1*div(i)
                     else
                        t7(neighc) = toldil
                        t9(neighc) = toldil
                     end if
                  enddo
                  
                  if(icnl.eq.0) then
                     do jm=1,iq
                        kb=it8(jm)
                        kz=kb-icd
                        neighc=it9(jm)
                        iw=it10(jm)
                        sx2c=sx(iw,isox)
                        sx3c=sx(iw,isoy)
                        sxzc=sx(iw,isoz)
                        jj = it11(jm) 
                        j = kb
                        vxy = t8(neighc)
                        
                        tempx = (cord(iz,1)-cord(kz,1))**2
                        tempy = (cord(iz,2)-cord(kz,2))**2
                        tempz = (cord(iz,3)-cord(kz,3))**2
                        dis=sqrt(tempx + tempy + tempz)
                        area_t=-(sx2c+sx3c+sxzc)*dis
                        divfi=t9(neighc)
                        divfkb=t7(neighc)
                        vxyf=max(toldil,(divfkb+divfi))
                        prodrovf = rovf(i) * rovf(kb)
                        sehvelv(sehindexv)=0
                        if(vxy.ne.0.) then
                           if(prodrovf.gt.0.) then
                              if(area_t.gt.0.) then
                                 sehvelv(sehindexv)=abs(vxy/vxyf/area_t)
     &                                *((rovf(kb)*divfi+rovf(i)*divfkb)/
     +                                prodrovf)
                              endif
                           endif
                        endif
                        if (dispsame.eq.1) then
                           newdiff = concadiff(2,mflagv(1,itrcdsp(i)),
     &                          sehdiffv(itrcdsp(i)),ps(i),s(i),
     &                          phi(i),t(i))
                           sehvelv(sehindexv)=(sx2c+sx3c+sxzc)
     &                          *(sehvelv(sehindexv)*
     +                          alphaconv(sehindexv)+
     &                          newdiff*ps(i)*(1-s(i)))
c----------------- phs 9/26/2001 - added (1-s(i)) correction
                        endif
                        sehindexv=sehindexv+1
                     enddo
                  else
                     radi=cord(iz,3)
                     do jm=1,iq
                        kb=it8(jm)
                        kz=kb-icd
                        neighc=it9(jm)
                        iw=it10(jm)
                        radkb=0.5*(radi+cord(kz,3))
                        sx2c=radkb*sx(iw,isox)
                        sx3c=radkb*sx(iw,isoy)
                        jj = it11(jm) 
                        j = kb
                        vxy = t8(neighc)
                        
                        tempx = (cord(iz,1)-cord(kz,1))**2
                        tempy = (cord(iz,2)-cord(kz,2))**2
                        dis=sqrt(tempx + tempy)
                        area_t=-(sx2c+sx3c)*dis*
     &                       0.5*(cord(iz,3)+cord(kz,3))
                        divfi=t9(neighc)
                        divfkb=t7(neighc)
                        vxyf=max(toldil,(divfkb+divfi))
                        prodrovf = rovf(i) * rovf(kb)
                        sehvelv(sehindexv)=0
                        if(vxy.ne.0.) then
                           if(prodrovf.gt.0.) then
                              if (area_t.gt.0.) then
                                 sehvelv(sehindexv)=abs(vxy/vxyf/area_t)
     &                                *((rovf(kb)*divfi+rovf(i)*divfkb)/
     +                                prodrovf)
                              endif
                           endif
                        endif
                        if (dispsame.eq.1) then
                           newdiff = concadiff(2,mflagv(1,itrcdsp(i)),
     &                          sehdiffv(itrcdsp(i)),ps(i),s(i),
     &                          phi(i),t(i))
                           sehvelv(sehindexv)=(sx2c+sx3c)*
     &                          (sehvelv(sehindexv)*
     +                          alphaconv(sehindexv)+
     &                          newdiff*ps(i)*(1-s(i)))
c----------------- phs 9/26/2001 - added (1-s(i)) correction
                        endif
                        sehindexv=sehindexv+1
                     enddo
                  endif
               endif
            enddo
         endif
      endif

      dayst=days
      daysi=days-day
      days=daysi
      daytr=min(daytr,day,daycmx/2.0)
      sehmindays=min(daymin,daycmm)
      nts=0
c following causes output whenever a transport time step is initiated
c      iprttrc = nprttrc-1
**** Begin Loop through time   
                      
 1000 continue
c     Add counter for printout of trc output BAR 11-18-98
      iprttrc = iprttrc + 1
      days = days + daytr
      if(days.lt.dayst) icfin=1
      if(days.eq.dayst) icfin=0
      if(days.gt.dayst) then
         days = days - daytr
         daytr = dayst - days
         days=dayst
         icfin = -1
      endif
      days=daysi+daytr
c
c adjust time step for tracer injection start or finish
c
      daytrm=daytr

c Call userc to determine recirculation, the rate is based
c  on previous time step concentrations so only needs to be called
c  once for each tracer timestep
      if (usroption .eq. 5) then
         call userc(3,0,dumv,dumv)
      end if

      do isolute = 1, ncpnt
         npn=npt(isolute)
         istop_flag = .false.
         do id=1,neq
            i=id+npn
c     
            if(abs(cnsk(i)).gt.0.0) then
               time1=t1sk(i)
               time2=t2sk(i)
               rcd=abs(rc(i))
               if(time1.ge.0.0) then
                  if(time1.le.days.and.rcd.le.1.e-7) then
                     daytr=time1-daysi
                     if ((daytr.lt.0.).or.(time1.eq.daycs)) then
                        daytr=daycmm
                     else
                        if(daytr.lt.sehmindays) then
                           t1sk(i)=daysi+sehmindays
                           daytr=sehmindays
                        endif
                     endif
                     if (iout .ne. 0) write(iout,400) id,t1sk(i)
 400                 format(1x,/,'tracer injection started for node ',
     2                    i7,' at days=',g14.6)
                     icfin=1
                     if (t1sk(i).ne.0.) then
                        t1sk(i)=-t1sk(i)
                     else
                        t1sk(i)=-1e-30
                     endif
                     days=daysi+daytrm
                  endif
               endif
               if(time2.gt.0.0) then
                  if(time2.le.days.and.time1.lt.0.0) then
                     daytr=time2-daysi
                     if(daytr.lt.sehmindays) then
                        t2sk(i)=daysi+sehmindays
                        daytr=sehmindays
                     endif
                     if (iout .ne. 0) write(iout,401) i,t2sk(i)
                     istop_flag = .true.
 401                 format(1x,/,'tracer injection stopped for node ',
     2                    i7,' at days=',g14.6)
                     icfin=1
                     t2sk(i)=-t2sk(i)
                  endif
                  daytrm=min(daytrm,daytr)
               endif
            endif
         end do
         if (istop_flag) call wrtcon(1)
      end do
      dtotc=daytrm*86400.
      daytr=daytrm
      days=daysi+daytrm
      daysi=days
      if(days.ge.dayst) icfin=0
c end change
      iter_counter = 0
      tajj = tyming(caz)
      reset_tracer = .FALSE.
 3000 continue
      call cnswer(tol_value, sia_iter)

c     Add counter for total SIA iterations

      sia_iter_tot = sia_iter_tot + sia_iter
      if( tol_value .gt. 0 ) then
***   start over with cut timestep
         call resettrc(daytr)
         reset_tracer = .TRUE.
         if (tol_value .eq. 1) then
            if (iout .ne. 0) then
               write(iout, 3010)
               write(iout, 3020) daytr
            end if
            if( iptty .gt. 0 ) then
               write(iptty, 3010)
               write(iptty, 3020) daytr
            endif
         elseif (tol_value .eq. 2) then
            if (iout .ne. 0) then
               write(iout, 3030)
               write(iout, 3020) daytr
            end if
            if( iptty .gt. 0 ) then
               write(iptty, 3030)
               write(iptty, 3020) daytr
            endif
         else
            if (iout .ne. 0) then
               write(iout, 3040)
               write(iout, 3020) daytr
            end if
            if( iptty .gt. 0 ) then
               write(iptty, 3040) 
               write(iptty, 3020) daytr
            endif
         endif
         tol_value = 0
 3010    format ('Number of SIA iterations exceeded')
 3020    format ('Time step reduced to ', g16.9, ' days')
 3030    format ('Negative concentration encountered')
 3040    format ('Newton-Raphson iteration limit exceeded in ',
     2        'speciation subroutine.')
      else
c*** convergence
         if(iprttrc.ge.abs(nprttrc).or.icfin.le.0) then
            if (nprttrc .gt. 0) then
               if (iout .ne. 0 .and. idebug .eq. 1) 
     1              write(iout, 3050) sia_iter, sia_iter_tot
            endif
         endif
         if(iptty .ne. 0 .and. idebug .eq. 1) then
            write(iptty, 3050) sia_iter, sia_iter_tot
         end if
 3050    format ('# SIA Iterations ', i5, ' (total = ', i5, ')')
         reset_tracer = .FALSE.

c check for porosity changes if salt simulation
         if(isalt.EQ.1) then
           call saltctr(2,0,0.0d00,0.0d00)
         else
c check for porosity changes for non salt simulation
          if(allocated(ps_delta_rxn)) then
           ps_delta_rxn = 0.0d0
          endif

         do im = 1, nimm
            nsp = pimm(im)
            npn = npt(nsp)

            if(mw_mineral(im).ne.0)then
               do i = 1, n0
                  ja = i + npn
                 ps_delta_rxn(i) = rc(ja)*dtotc*mw_mineral(im)
     &                                 /rho_mineral(im)

                 ps_delta_rxn(i) = ps_delta_rxn(i)/sx1(i)
                 an(ja) = anlo(ja)
 
                 if(ps(i).LE.ps_min) ps_delta_rxn(i) = 0.0
                 if(ps(i).GT.ps_max) ps_delta_rxn(i) = 0.0
                 if(s(i).LT.s_min_salt) ps_delta_rxn(i) = 0.0

               enddo
            endif                   
            
            call porosi(1)

            ps_delta_rxn_s = ps_delta_rxn


         enddo
        end if ! isalt.NE.1
      end if  ! convergence check
 3011 format ('*****************************************************')

      if( .not. reset_tracer ) then
         tajj = tyming(caz) - tajj            
         if(iptty .ne. 0) then
            write(iptty, *)
            write(iptty, 3011) 
            write(iptty, *) 'Convergence of tracer solution at time = ',
     2          days, ' days'
            write(iptty ,*) 'cpu seconds for tracer solution at this',
     2            ' time = ', tajj
            write(iptty, 3011) 
            write(iptty, *)
         end if
c gaz debug 073114
         vdum = anl(1)
         vdum = anv(1)
c
c     Compute mass balance terms
c
         ic = 0
         im = 0
         iv = 0
         do nsp = 1, nspeci
            call thermc(0)
            if(idpdp.ne.0)call thermc(neq)
            npn=npt(nsp)
            cm(nsp)=0.
            do jad=1,n
               ja=jad+npn
               vdum=volume(ja-npn)
               if(rcss(ja) .le. 0.) qcin(nsp)=qcin(nsp)+rcss(ja)*dtotc
               if(rcss(ja) .gt. 0.) qcout(nsp)=qcout(nsp)+rcss(ja)*dtotc
               qcrxn(nsp) = qcrxn(nsp) + ( rc(ja) ) * dtotc
               dencht=dench(ja)
               dench(ja)=dencht+denci(ja)*dtotc
               dencj(ja) = denci(ja)
               denci(ja) = 0.
c Do we need to do this here? This is done at start of cnswer?
c We will need to have old value for call to cfluxz so if we have to do
c this it should be after call to cfluxz
c               anlo(ja) = an(ja)
c      ------------------------------------------------------------------
c      PHS 3/28/05 adding integration of sink at each node for Don Neeper
               sinkint(ja) = sinkint(ja) + rcss(ja)*dtotc
               if (ps_trac(jad) .gt. 0) then
                  cm(nsp)=cm(nsp)+dench(ja)*vdum
               end if
            end do
            continue
            call bcon(3)
c
c call plot output
c
         end do
c     Only write to trc if it is time BAR 11-18-98
c     Modify check for time to print ZVD 16-Sep08
c nhist = 1   : Output every heat and mass time step
c nprttrc = 1 : Output every tracer time step

         if (time_flag .eq. 1) then
            ptime = abs(days / 365.25d00)
         else if (time_flag .eq. 3) then
            ptime = abs(days * 86400.d00)
         else if (time_flag .eq. 4) then
            ptime = abs(days * 24.d00)
         else
            ptime = abs(days)
         end if
c
         if (l .eq. 1 .and. iprttrc .eq. 1 .and. last_time .eq. 0.) then
            time2print = .TRUE.
         else if (ifinsh .ne. 0) then
            if (ptime .ne. last_time) then
               time2print = .TRUE.
            else
               time2print = .FALSE.
            end if
         else            
            if (ifinsh .ne. 2 .and. ptime .lt. (last_time + histime)
     &           .and. l .lt. (last_step)
c    &           .and. l .lt. (last_step + nhist)
     &           .and. iprttrc .lt. nprttrc) then
               time2print = .FALSE.
            else
               time2print = .TRUE.
            end if
         end if

         if (time2print) then
c     if(iprttrc.ge.abs(nprttrc).or.icfin.le.0) then
c     iprttrc = 0
            do ii = 1,ncpntprt
               ic = cpntprt(ii)
               nsp = pcpnt(ic)
               npn = npt(nsp)
               call plotc1(2,ic)
            enddo
            do ii = 1,nimmprt
               im = immprt(ii)
               nsp = pimm(im)
               npn = npt(nsp)
               call plotc1(3,im)
            enddo
            do ii = 1,nvapprt
               iv = vapprt(ii)
               nsp = pvap(iv)
               npn = npt(nsp)
               call plotc1(4,iv)
            enddo
            if(ncplx.ne.0)then
               do ii = 1,ncpntprt
                  ic = cpntprt(ii)
                  call plotc1(5,ic)
               enddo
               do ii = 101, ncplxprt+100
                  ix = cplxprt(ii) 
                  call plotc1(6,ix)
               enddo
            endif
            if (cflxz .ne. 0) then
               call cfluxz(days)
            end if
!            if (nprttrc .gt. 0) then
!               call wrtcon(0)
!            end if
            last_time = ptime
            last_step = l
            iprttrc = 0
         end if
         nts=nts+1
      end if
c gaz 071016 iadd does not seem to be used, put sia info in iadd for printout
      iadd(1) = sia_iter_tot
      iaddt(1) = sia_iter_tot_old+sia_iter_tot  
      if( .not. reset_tracer ) then
         if(icfin.le.0) then
            dtotc=dtotdm/nts
            go to 4000
         endif
         if(sia_iter.le.iaccmx) daytr = daytr * daycm
         if(daytr.gt.daycmx) daytr = daycmx          
      end if
      go to 1000
 4000 continue
      end
