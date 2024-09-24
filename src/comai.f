      module comai
!     comai                                              
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global include file for scalar and character variables (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 30-SEP-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comai.f_a  $
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   Identifier      Type     Description
!D4   
!D4   nmacros         int      Number of input macros
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block faai variables *****
!D4   iab             INT      faai   Iteration counter in stress routines
!D4   iac             INT      faai   Counter for print-out interval
!D4   iaccmx          INT      faai   Maximum iterations allowed for time step
!D4                                     increase (tracer solution)
!D4   iad             INT      faai   Current iteration number in flow solution
!D4   iad_up          INT      faai   Parameter for macro itup
!D4   iadif           INT      faai   ?
!D4   iamm            INT      faai   Maximum iterations allowed for time step
!D4                                     increase (heat and mass solution)
!D4   iamx            INT      faai   Iteration count after which the time step 
!D4                                     will be halved
!D4   iatty           INT      faai   Unit number for all tty output
!D4   icapp           INT      faai   Indicates capillary pressure model
!D4   iccen           INT      faai   Parameter which indicates if the tracer 
!D4                                     solution is enabled
!D4   ice             INT      faai   Parameter which indicates if the ice 
!D4                                     solution is enabled
!D4   icf             INT      faai   Parameter indicating status of tracer 
!D4                                     solution
!D4   icgts           INT      faai   Parameter controlling the time of solution
!D4                                     parameter changes
!D4   ichng           INT      faai   Number of thermodynamic region changes
!D4   icnl            INT      faai   Problem dimension
!D4   ico2            INT      faai   Indicates if noncondensible gas solution 
!D4                                     is enabled (co2i)
!D4   icontr          INT      faai   Parameter used in contour plot management
!D4   ics             INT      faai   Parameter indicating status of tracer 
!D4                                     solution
!D4   idof            INT      faai   Number of degrees of freedom per node for 
!D4                                     the current problem
!D4   idpdp           INT      faai   Parameter which indicates if the double
!D4                                     porosity / double permeability
!D4                                     solution is enabled (dpdp)
!D4   idualp          INT      faai   Parameter which indicates if the dual 
!D4                                     porosity solution is enabled (dual)
!D4   ierr            INT      faai   Unit number for error output
!D4   ifinsh          INT      faai   Indicates if the finishing criteria for 
!D4                                     the simulation is achieved
!D4   iflag           INT      faai   Flag used in input subroutine
!D4   igauss          INT      faai   Flag used form maximum gauss 
!D4                                     elimination order
!D4   igrav           INT      faai   Direction of gravity in problem
!D4   ihf             INT      faai   Parameter indicating status of flow 
!D4                                     solution (trac)
!D4   ihs             INT      faai   Parameter indicating status of flow 
!D4                                     solution (trac)
!D4   iprtout         INT      faai   Print-out interval, number of time steps 
!D4                                     (time)
!D4   ilt             INT      faai   Parameter used in time step control
!D4   imdno           INT      faai   Indicates if md_nodes subroutine will be 
!D4                                     called
!D4   impf            INT      faai   ?
!D4   incoor          INT      faai   Unit number for coordinate input file
!D4   inpt            INT      faai   Unit number for input file
!D4   intg            INT      faai   Indicates integration type used
!D4   inzone          INT      faai   Unit number for zone input file
!D4   iocntl          INT      faai   Unit number for control file
!D4   iout            INT      faai   Unit number for output file
!D4   iporos          INT      faai   Indicates if deformation model is enabled 
!D4                                     (ppor)
!D4   iptty           INT      faai   Unit number for selected tty output
!D4   ipqz            INT      faai   The number of nodes used for the contour 
!D4                                     plot output
!D4   ipsat           INT      faai   Parameter regulating call to subroutine 
!D4                                     psat (eos)
!D4   iread           INT      faai   Unit number for restart file (to read)
!D4   ireord          INT      faai   Parameter to reorder system of linear 
!D4                                     equations
!D4   irpd            INT      faai   Indicates relative permeability model
!D4   isave           INT      faai   Unit number for restart file (to write)
!D4   ischk           INT      faai   Unit number for input data check file
!D4   iscon           INT      faai   Unit number for contour data file
!D4   iscon1          INT      faai   Unit number for dual porosity or dpdp 
!D4                                     contour data file
!D4   ishis           INT      faai   Unit numbers for history data file
!D4    ishisp ishist ishishd ishiss ishisf ishise ishishm ishisfz ishisg
!D4    ishisc, ishiswc
!D4   ismis           INT      faai   Unit number for miscellaneous save file
!D4   ispest          INT      faai   Unit number for pest output file 
!D4   ispst1          INT      faai   Unit number for aux. pest output file 
!D4   isptr[1-9]      INT      faai   Unit numbers for streamline particle 
!D4                                     tracking output files
!D4   isstor          INT      faai   Unit number for element coefficient file
!D4   istrc           INT      faai   Unit number for tracer history data file
!D4   istrs           INT      faai   Parameter indicating if the stress 
!D4                                     solution is enabled
!D4   iter            INT      faai   Number of iterations of last call to 
!D4                                     solve(n) subroutine
!D4   itert           INT      faai   Intermediate iteration counter
!D4   ithic           INT      faai   Indicates if thickness subroutine will be 
!D4                                     called
!D4   itotal          INT      faai   Total iteration count
!D4   itsat           INT      faai   Parameter specifying the setting of 
!D4                                     saturation temperature
!D4   iupk            INT      faai
!D4   ivapl           INT      faai   ?
!D4   ivf             INT      faai   ?
!D4   ivfcal          INT      faai   Indicates if vfcal subroutine will be 
!D4                                     called
!D4   iw              INT      faai   Number of storage locations needed to 
!D4                                     store geometric input types
!D4   iwelb           INT      faai   Indicates if wellbore solution is enabled
!D4   iyear           INT      faai   Current year in simulation
!D4   l               INT      faai   Current time step number
!D4   lda             INT      faai   Parameter which specifies if the geometric
!D4                                     coefficients are saved (ctrl)
!D4   m               INT      faai   Total number of nodes used for output 
!D4                                     information (node)
!D4   m2              INT      faai   ?
!D4   maxit           INT      faai   Maximum number of iterations allowed 
!D4                                     before time step is halved (ctrl)
!D4   maxsolve        INT      faai   Maximum solver iterations
!D4   mink            INT      faai   Number of active variables
!D4   mlz             INT      faai   Out of bounds node
!D4   month           INT      faai   Current month in simulation
!D4   n               INT      faai   Total number of nodes
!D4   nbnd            INT      faai   Maximum number of nonzeros per row in the 
!D4                                     incidence matrix
!D4   ncntr           INT      faai   Contour plot interval (cont)
!D4   ndem            INT      faai   Coordinate direction perpendicular to the 
!D4                                     contour plane
!D4   nei             INT      faai   Total number of elements in the problem 
!D4                                     (elem)
!D4   neigh           INT      faai   Maximum number of neighbors occur in 
!D4                                     tracer solution
!D4   nemx            INT      faai   Number of unique (geometrically) elements
!D4   neq             INT      faai   Number of nodes, not including dual 
!D4                                     porosity nodes (cord)
!D4   nflx            INT      faai   
!D4   nflxz           INT      faai   
!D4   ni              INT      faai   Number of integration points per element
!D4   nicg            INT      faai   Parameter used in time step control
!D4   north           INT      faai   Maximum number of orthogonalizations 
!D4                                     allowing in (ctrl)
!D4   ns              INT      faai   Number of nodes per element (elem)
!D4   nsave           INT      faai   Indicates if a restart file will be 
!D4                                     created in the current problem
!D4   nsnk            INT      faai   Number of source or sink terms
!D4   nstep           INT      faai   Maximum number of time steps
!D4   ntty            INT      faai   Flag for file and terminal output
!D4   nptrtrc         INT      faai   Flag for trc file output
!D4
!D4   ***** COMMON Block faar variables *****
!D4   aener           REAL*8   faar   Current energy in  problem
!D4   aiaa            REAL*8   faar   Time step multiplication factor (ctrl)
!D4   aiar            REAL*8   faar   Time step reduction factor (ctrl)
!D4   am0             REAL*8   faar   Initial mass in problem
!D4   amass           REAL*8   faar   Current mass in problem
!D4   ame             REAL*8   faar   Initial energy in problem
!D4   an0             REAL*8   faar   Initial tracer concentration (trac)
!D4   asteam          REAL*8   faar   Current steam mass in problem
!D4   astmo           REAL*8   faar   Initial steam mass in problem
!D4   aw              REAL*8   faar   Time step weighting parameter for heat and
!D4                                     mass solution (ctrl)
!D4   awc             REAL*8   faar   Time step weighting parameter (trac)
!D4   awt             REAL*8   faar   Value of implicitness factor
!D4   ay              REAL*8   faar   Time step weighting parameter for tracer 
!D4                                     (trac)
!D4   ayc             REAL*8   faar   Time step weighting parameter for tracer
!D4   contim          REAL*8   faar   Interval (days) for contour plot output
!D4   day             REAL*8   faar   Current time step size in days
!D4   daycf           REAL*8   faar   Time at which tracer solution stops (trac)
!D4   daycm           REAL*8   faar   ?
!D4   daycmm          REAL*8   faar   Minimum time step size for tracer solution
!D4   daycmx          REAL*8   faar   Maximum time step size for tracer solution
!D4   daycs           REAL*8   faar   Time at which tracer solution starts 
!D4                                     (trac)
!D4   dayhf           REAL*8   faar   Time at which flow solution stops (trac)
!D4   dayhs           REAL*8   faar   Time at which flow solution starts (trac)
!D4   daymax          REAL*8   faar   Maximum time step allowed (ctrl)
!D4   daymin          REAL*8   faar   Minimum time step allowed (ctrl)
!D4   daynew          REAL*8   faar   Parameter used in time step control
!D4   days            REAL*8   faar   Current simulation time
!D4   dayscn          REAL*8   faar   Next time for contour plot
!D4   daysi           REAL*8   faar   Simulation time at last time step
!D4   daysp           REAL*8   faar   Time at next time step
!D4   delat           REAL*8   faar   Given tolerance for gas pressure
!D4   delpt           REAL*8   faar   Given tolerance for pressure
!D4   delst           REAL*8   faar   Given tolerance for saturation
!D4   deltt           REAL*8   faar   Given tolerance for temperature
!D4   den             REAL*8   faar   Intermediate value of mass accumulation 
!D4                                     term
!D4   dene            REAL*8   faar   Intermediate value of energy accumulation 
!D4                                     term
!D4   depcng          REAL*8   faar   Depth at which temperature gradient 
!D4                                     changes
!D4   dife            REAL*8   faar   Energy balance error
!D4   difm            REAL*8   faar   Mass balance error
!D4   ditnd           REAL*8   faar   Next time step change time
!D4   dnwgt           REAL*8   faar   Upwind weighting parameter
!D4   dnwgta          REAL*8   faar   Upwind weighting parameter in tracer 
!D4                                     solution
!D4   dtot            REAL*8   faar   Current time step size in seconds
!D4   dtotc           REAL*8   faar   Tracer time step size in seconds
!D4   dtotdm          REAL*8   faar   Last time step size in seconds
!D4   emiss           REAL*8   faar   ?                               
!D4   epc             REAL*8   faar   Specified solution tolerance for tracer 
!D4                                     solution
!D4   epe             REAL*8   faar   Tolerance for newton-raphson iteration
!D4   eps             REAL*8   faar   Tolerance for linear equation solver
!D4   f0              REAL*8   faar   Initial time step residual
!D4   fdum            REAL*8   faar   Current sum squared of residuals
!D4   fdum1           REAL*8   faar   ?
!D4   fimp            REAL*8   faar   Fraction of variables over a given
!D4                                     tolerance
!D4   g1              REAL*8   faar   Iteration accuracy control parameter 
!D4                                     (iter)
!D4   g2              REAL*8   faar   Iteration accuracy control parameter 
!D4                                     (iter)
!D4   g3              REAL*8   faar   Iteration accuracy control parameter 
!D4                                     (iter)
!D4   grad2           REAL*8   faar   Parameter in description of temperature 
!D4                                     gradient
!D4   gradnt          REAL*8   faar   Parameter used in description of 
!D4                                     temperature gradient
!D4   grav            REAL*8   faar   Value of gravity
!D4   overf           REAL*8   faar   Over relaxation factor for sor equations 
!D4                                     (iter)
!D4   pein            REAL*8   faar   Initial pressure of problem (init)
!D4   pow             REAL*8   faar   Power output for a given time step
!D4   qt              REAL*8   faar   Total outflow for time step
!D4   qte             REAL*8   faar   Total energy outflow for time step
!D4   qtot            REAL*8   faar   Total outflow for problem
!D4   qtote           REAL*8   faar   Total energy outflow for problem
!D4   qtotei          REAL*8   faar   Intermediate energy flow total
!D4   qtoti           REAL*8   faar   Intermediate flow total
!D4   quad            REAL*8   faar   Parameter used in temperature gradient
!D4   rnmax           REAL*8   faar   Maximum run time allowed
!D4   str             REAL*8   faar   Multiplier for newton-raphson corrections
!D4   strd            REAL*8   faar   Multiplier for newton-raphson corrections
!D4   sv              REAL*8   faar   Vapor saturation of a node
!D4   teoutf          REAL*8   faar   Parameter used in calculation of energy 
!D4                                     output
!D4   tims            REAL*8   faar   Ending simulation time (time) 
!D4   tin             REAL*8   faar   Parameter used in temperature gradient
!D4   tin0            REAL*8   faar   Initial problem temperature (init)
!D4   tin1            REAL*8   faar   Parameter used in temperature gradient
!D4   tmch            REAL*8   faar   Machine tolerance (iter)
!D4   tort            REAL*8   faar   ?
!D4   toutfl          REAL*8   faar   Parameter used in calculation of flow 
!D4                                     output
!D4   upwgt           REAL*8   faar   Upwind weighting parameter
!D4   upwgta          REAL*8   faar   Upwind weighting factor for tracer 
!D4                                     solution (trac)
!D4   vernum          REAL*4          Code version number used by GoldSim
!D4   vlmax           REAL*8   faar   Maximum liquid phase velocity
!D4   vtot            REAL*8   faar   Total volume in problem
!D4   vvmax           REAL*8   faar   Maximum vapor phase velocity
!D4
!D4   ***** COMMON Block faac variables *****
!D4   accm            CHAR     faac   Acceleration method for solve
!D4   altc            CHAR     faac   String used for denoting type of output
!D4                                     for contour plots
!D4   sssol           CHAR     faac   Indicates if initial steady state solution
!D4                                     is needed
!D4   verno           CHAR     faac   Contains version identifier of FEHMN code
!D4   wdd             CHAR     faac   Character input string
!D4   wdd1            CHAR     faac   Alternate character input string
!D4
!D4   ***** COMMON Block faac1 variables *****
!D4   jdate           CHAR     faac1  Contains the date (mm/dd/yr)
!D4   jtime           CHAR     faac1  Contains the time (hr:mn:sc)
!D4
!D4   ***** COMMON Block macro variables *****
!D4   macroread(1)    LOGICAL  macro  Flag denoting if macro cap has been read
!D4   macroread(2)    LOGICAL  macro  Flag denoting if macro ngas has been read
!D4   macroread(3)    LOGICAL  macro  Flag denoting if macro dual has been read
!D4   macroread(4)    LOGICAL  macro  Flag denoting if macro ppor has been read
!D4   macroread(5)    LOGICAL  macro  Flag denoting if macro trac has been read
!D4   macroread(6)    LOGICAL  macro  Flag denoting if macro dpdp has been read
!D4   macroread(7)    LOGICAL  macro  Flag denoting if macro rlp has been read
!D4   macroread(8)    LOGICAL  macro  Flag denoting if macro ice has been read
!D4   macroread(9)    LOGICAL  macro  Flag denoting if macro vcon has been read
!D4   macroread(10)   LOGICAL  macro  Flag denoting if macro cond has been read
!D4   macroread(11)   LOGICAL  macro  Flag denoting if macro ctrl has been read
!D4   macroread(12)   LOGICAL  macro  Flag denoting if macro flow has been read
!D4   macroread(13)   LOGICAL  macro  Flag denoting if macro hflx has been read
!D4   macroread(14)   LOGICAL  macro  Flag denoting if macro perm has been read
!D4   macroread(15)   LOGICAL  macro  Flag denoting if macro pres has been read
!D4   macroread(16)   LOGICAL  macro  Flag denoting if macro renm has been read
!D4   macroread(17)   LOGICAL  macro  Flag denoting if macro rock has been read
!D4   macroread(18)   LOGICAL  macro  Flag denoting if macro zone has been read
!D4   macroread(19)   LOGICAL  macro  Flag denoting if macro thic has been read
!D4   macroread(20)   LOGICAL  macro  Flag denoting if macro mdno has been read
!D4   macroread(21)   LOGICAL  macro  Flag denoting if macro ptrk has been read
!D4   macroread(22)   LOGICAL  macro  Flag denoting if macro fper has been read
!D4   macroread(23)   LOGICAL  macro  Flag denoting if macro perm_olivella has been read
!D4   macroread(24)   LOGICAL  macro  Flag denoting if macro  den(spatial density) has been read      
!D4   tpor_flag       LOGICAL         Flag denoting if tracer porosity has been read
!D4 
!D4 Global Subprograms
!D4
!D4   None
!D4
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 N/A
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS None
!PS
!***********************************************************************
c gaz 102919 added iad_min(calculated in bnswer)
      integer             iad_min
      integer             iab, iac, iaccmx, iad, iad_up, iadif, iamm
      integer             iamx, iatty, icapp, iccen, ice, icf, icgts 
      integer             ichng, icnl, ico2, icontr, ics, idof, idpdp 
      integer             idualp, ierr, ifinsh, iflag, igauss, igrav 
      integer             ihf, ihs, iprtout, ilt, imdnode, impf, incoor
      integer             inpt, intg, inzone, iocntl, iout, iporos, ipqz
      integer             ipsat, iptty, iread, ireord, irpd, isave
      integer             ischk, iscon, iscon1, ishis, ismis, isstor
      integer             istrc, istrs, iter, itert, ithic, itotal
      integer             itsat, iupk, ivapl, ivcond, ivelo, ivf
      integer             ivfcal, iw, icons, inapl, iboun, ihead, ichead
      integer             iwelb, iyear, l, lda, m, m2, m3, maxit, mink 
      integer             minkt, mlz, month, n, nbnd, ncntr, nei 
      integer             neigh, nemx, neq, nflx, nflxz, ni, nicg
      integer             north, ns, itotals, nflxt, iflxc, irestart_ts
      integer             nsave, nstep, ntpp, ntty, irsttime 
      integer             ispest, ispst1, ipest, mpest, iexrlp, ivar
c gaz 052322  added global varible ivar_switch      
      integer             ivar_switch
c gaz_082622 added global var nr_completed indicates in NR update just finished     
      integer             nr_completed
      integer             isptr1, isptr2, isptr3, isptr4, isptr5, isptr6
      integer             isptr7, isptr8, isptr9, idoff, rlp_flag
      integer             nprttrc, irlp_fac, interface_flag
      integer             neq_primary, ngdpm_actual, gdpm_flag
      integer             gdkm_flag, nr_stop, neq_gdkm
      integer             ngdpmnodes, maxgdpmlayers, ianpe, ipara
      integer             maxsolve, iflux_ts, inobr, isubbc, isubm
      integer             bin_flag, app_flag, ldna, iflxn, pres_flag
      integer             ishisp, ishist, ishishd, ishiss, ishisf
      integer             ishisfa
      integer             ishise, ishishm, ishisg, ishisfz, ishiswt 
      integer             time_flag, nhist, glob_flag, icont, istea_pest
      integer             icoef_replace
c RJP 12/13/06 added nriver below
      integer             nstep_save, nriver, nrlp, nei_in, ns_in, i_rlp
c ZVD 05/01/07 added form_flag, ishisc, ishiswc
c ZVD 08/05/09 added ishisp2, ishisp3, ishiscsl, ishiscsg
      integer             form_flag, ishisc, ishiswc, ishisp2, ishisp3
c ZVD 01/05/11 added ihisef
      integer             ishisef
c ZVD 07/21/11 added for viscosity and density history output
      integer             ishisd, ishisd2, ishisd3, den_flag
      integer             ishisv, ishisv2, ishisv3, vis_flag
c RJP added for history of carbon mass
      integer             ishiscm, ishiscmf, ishiscmd
      integer             ishiscfz, ishiscsl, ishiscsg 
      integer             imbe_ctr, imbet_ctr
c GAZ 11/02/08 
      integer             ishisstr, ishisstrx, ishisstry, ishisstrz
      integer             ishisstrxy, ishisstrxz, ishisstryz
      integer             ishisdisx, ishisdisy, ishisdisz
      integer             ich_max, ich_m1, ich_m2
      integer             jswitch, joff
      integer             iwellp_flag, icoef_neg, numcoef_neg
      integer             istrs_coupl, iax_b, ipini
      integer             ivboun, itert_s,itotal_s,itotals_s
      integer             ifdm_elem, ienrich, enri_flag
      integer             nenrinodes
      integer             ibcfar,i_subsid,i_init
c GAZ 061513
      integer             ich_pebi
c GAZ 072513      
      integer             isalt,ivaprsalt,isalt_pnx,isalt_ps,isalt_read
      integer             isalt_write
      integer             iactive
c GAZ 011014
c GAZ 012214
      integer             iden_vis, ideng_vis, iphase3
c gaz 110823 added itype_h2
      integer             itype_air,itype_meth,itype_co2,itype_h2
      real*8              pchng, schng 
c gaz 121314
      integer              nrestart_ts
c gaz 081415
      integer              mlz_save   
      integer              nphase_liq, nphase_2, nphase_gas, nphase_sc 
      integer              nphase_liq_0, nphase_2_0, nphase_gas_0
      integer              nphase_sc_0
      integer              dnphase_liq, dnphase_2, dnphase_gas
      integer              dnphase_sc
c gaz 110715, 070720, 081921 
      integer              iwater_table, iair_table, num_eos_table
      integer              ico2wh_table
      integer              itable_files 
      integer :: irun = 0
c gaz 013116
      integer              iaprf
      logical              gdkm_new
c gaz 111216      
      integer              izone_sv_cnt, num_sv_zones, icflux, icconc
      logical              sv_combine
c gaz 081117              
      integer             ivrock 
c gaz 100318      
      integer             initdata_pad
c gaz 042119              neq_active added to report active varables when porosity <0 nodes are disabled      
      integer             neq_active
c gaz 062920              added igdpm_add to facilitate adding gdpm nodes in zone call (after gdpm nodes read)
      integer             igdpm_add
c gaz 111520
      integer             iter_intrvl
      real*8              aener, aiaa, am0, amass, ame, an0, asteam
      real*8              astmo, aw, awc, awt, ay, ayc, contim, day
      real*8              daycf, daycm, daycmm, daycmx, daycs, dayhf
      real*8              dayhs, daymax, daymin, daynew, days, dayscn
      real*8              daysi, daysp, delat, delpt, delst, deltt, den
      real*8              dene, depcng, dife, difm, ditnd, dnwgt
      real*8              dnwgta, dtot, dtotc, dtotdm, emiss, epc, epe
      real*8              eps, f0, fdum, fdum1, fimp, g1, g2, g3, grad2
      real*8              gradnt, grav, overf, pein, pow, qt, qte, qtot
      real*8              qtote, qtotei, qtoti, quad, rnmax, str, strd
      real*8              sv, teoutf, tims, tin, tin0, tin1, tmch, tort
      real*8              toutfl, upwgt, upwgta, vlmax, vtot, vvmax
      real*8              rsttime, sx_mult, enday, time_ss, phi_inc
      real*8              contim_rip, histime, aiar, tims_trans
      real*8              daycs_save, daycf_save, dayhf_save, dayhs_save
      real*8              day_save_ss, fdum_last, aiaa_save, daymax_save
      real*8              tmch_save
c gaz 012521 add some boun macro variables     
      real*8              tmchboun, daymaxboun
c gaz 11-06-2007
      real*8              rho1grav, h_to_p, sat_ich
c gaz 01-06-2014
      real*8              roc0, visc_gas
      real*4              vernum
c gaz 11-29-2009   
      real*8              p_stop,t_stop,s_stop,pa_stop,s2_stop
      real*8              co2f_stop,h_stop    
c gaz 070821
      real*8              p_tol, t_tol, pc_tol
c gaz 032622 nr stopping criteria on raw or fractional equation balance errors
      real*8              eqwm_stop,eqen_stop,eqnm_stop
      real*8              water_mass_err, energy_err, ngas_mass_err
      real*8              water_mass_err_max, energy_err_max 
      real*8              ngas_mass_err_max
      integer             node_water_err_max, node_energy_err_max
      integer             node_ngas_err_max
      logical             var_stop, con_eq_stop
c gaz 090113
      real*8              permavg_salt, poravg_salt

      real*8              strd_iter, strd_rich, tol_phase
c gaz 122623
      real*8              strd_mass 
      character* 3        header_flag                    
      character* 4        altc, accm
      character* 5        sssol
      character* 8        jtime
      character*11        jdate, flux_flag
      character*30        verno
c gaz 021521 overwrite verno in dated with verno read in data file (scanin checked also)    
      character*30        verno_fehmid
      integer             idfehm
c gaz 100318 added wdd2      
      character*80        wdd, wdd1,wdd2
      character*80        salt_read_file, salt_write_file
c gaz 092319 h2o crit T and P moved here
c gaz 110819 pref, tref moved here (now global)  
      real*8 pcrit_h2o, tcrit_h2o, pref, tref
c gaz 112721 removed parameter line because pcrit_h2o=22.064d0, tcrit_h2o=373.946  
c can be changed (slightly) by SC water EOS table    
c      parameter(pcrit_h2o=22.064d0, tcrit_h2o=373.946) 
      real*8 pcrit_h2o_true, tcrit_h2o_true
      integer nmacros
      parameter( nmacros = 60 )
      logical macroread(nmacros) 

      logical tpor_flag, compute_flow, hist_flag, reverse_flow, boun_out
      logical out_zones, spercent, mass_read, pres_read, co2_read
      logical connect_out, resid_out, wflux_flag, vflux_flag
      logical eflux_flag,fperm_flag
c gaz 120421 moved eval_test_h2o from fluid control routines    
      integer eval_test_h2o
c gaz 020522
      integer i_vtk

      end module comai

