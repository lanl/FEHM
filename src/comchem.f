      module comchem
!        comchem
!***********************************************************************
!  Copyright, 1996, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Include file containing passed parameters and pointers related to
!D1 solute module.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 Initial implementation: 15-JUN-96, Programmer: Hari Viswanathan
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comchem.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:20   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2
!**********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!**********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   Identifier    Type     Description
!D4   
!D4   maxrxn     INT     maximum # kinetic reactions
!D4   maxspc     INT     maximum # aqueous components
!D4   
!D4   ***** COMMON Block chemnum variables *****
!D4   ncpnt      INT     # of aqueous components
!D4   ncplx      INT     # of aqueous complexes
!D4   nimm       INT     # of immobile components
!D4   nvap       INT     # of vapor components
!D4   ***** COMMON Block nprints variables *****
!D4   ncpntprt   INT     # of aqueous components to be printed
!D4   ncplxprt   INT     # of aqueous complexes to be printed
!D4   nimmprt    INT     # of immobile components to be printed
!D4   nvapprt    INT     # of vapor components to be printed  
!D4   ***** COMMON Block  idnrxn variables *****
!D4   naqsp      INT     # of aqueous components in a 
!D4                      kinetic reaction
!D4   nimsp      INT     # of immobile components in a 
!D4                      kinetic reaction
!D4   nivsp      INT     # of vapor components in a 
!D4                      kinetic reaction
!D4   irxnic     INT     #'s of the aqueous components
!D4                      in a kinetic reaction
!D4   irxnim     INT     #'s of the immobile components
!D4                      in a kinetic reaction
!D4   irxniv     INT     #'s of the vapor components
!D4                      in a kinetic reaction
!D4   ***** COMMON Block  rdsp variables *****
!D4   rdsmax     REAL*8  tolerance for the speciation solver
!D4   ckeq       REAL*8  equilibrium constants for speciation
!D4   heq        REAL*8  enthalphy's of formation associated
!D4                      with speciation reactions (for Van 
!D4                      Hoff equation temperature dependency
!D4                      model)
!D4   spstoic    REAL*8  stoichiometry for the speciation 
!D4                      reactions
!D4   cpntgs     REAL*8  guesses used for each total aqueous 
!D4                      concentration in the Newton-Raphson
!D4                      speciation solver
!D4   idrxn         INT  type of kinetic reaction 
!D4                      1 - Linear Kinetic Adsorption Model
!D4                      2 - Langmuir Kinetic Adsorption Model
!D4                      3 - General Kinetic Reaction
!D4                      4 - Biodegradation Model
!D4                      5 - Radioactive Decay Model
!D4                      6 - Henry's Law Equilibrium Reaction
!D4   ifxconc       INT  Options for specifying concentration
!D4                      0 - concentration is not fixed
!D4                      1 - free ion concentration corresponding
!D4                          to this total aqueous component is
!D4                          fixed
!D4                      2 - used to specify that component
!D4                          represents the H+ ion (for pH calcs)
!D4   ***** COMMON Block  init_chem variables *****
!D4   cpntsv     REAL*8  contains all of the free ion concentrations
!D4                      of each total aqueous component
!D4   calcderiv  LOGICAL indicates whether speciation derivatives
!D4                      should be taken
!D4   ***** COMMON Block place variables *****
!D4   pcpnt         INT  pointer array which translates a
!D4                      position in a NCPNT array like RRCPNT
!D4                      to a position is NSPECI array like AN
!D4   pimm          INT  pointer array which translates a
!D4                      position in a NIMM array like RRIMM
!D4                      to a position is NSPECI array like AN
!D4   pvap          INT  pointer array which translates a
!D4                      position in a NVAP array like RRVAP
!D4                      to a position is NSPECI array like AN
!D4   ***** COMMON Block couple_index variables *****
!D4   drcpos        INT  pointer array indicating which reaction
!D4                      derivatives are being used in SIA iteration
!D4   nderivs       INT  number of derivative submatrices that
!D4                      need to be stored
!D4   matpos        INT  pointer array indicating which reaction
!D4                      derivatives are being used in SIA iteration
!D4   dimdrc        INT  used to calculate the length of drdctaq
!D4   ***** COMMON Block print_flag variables *****
!D4   cpntprt   INT      aqueous components to be printed
!D4   cplxprt   INT      aqueous complexes to be printed
!D4   immprt    INT      immobile components to be printed
!D4   vapprt    INT      vapor components to be printed  
!D4   ***** COMMON Block chem_name variables *****
!D4   cpntnam   CHAR*20  stores the names of the aqueous components
!D4   cplxnam   CHAR*20  stores the names of the aqueous complexes
!D4   immnam    CHAR*20  stores the names of the immobile components
!D4   vapnam    CHAR*20  stores the names of the vapor components
!D4   ***** COMMON Block rxnpara variables *****
!D4   ckeqlb    REAL*8   the distribution coefficient
!D4                      for the kinetic adsorption reactions
!D4   ckmtrn    REAL*8   the mass transfer coefficient for
!D4                      the kinetic adsorption reactions
!D4   simmmx    REAL*8   maximum sorption capacity for 
!D4                      Langmuir model
!D4   ***** COMMON Block irrpar variables *****
!D4   kfor      REAL*8   forward rate constant for general
!D4                      kinetic reaction
!D4   krev      REAL*8   reverse rate constant for general
!D4                      kinetic reaction
!D4   sticirrv  REAL*8   stoiciometry of aqueous components/
!D4                      complexes for the general kinetic
!D4                      reaction
!D4   stimirrv  REAL*8   stoiciometry of immobile components
!D4                      for the general kinetic reaction
!D4   stivirrv  REAL*8   stoiciometry of vapor components
!D4                      for the general kinetic reaction
!D4   ***** COMMON Block biopara variables *****
!D4   ckc       REAL*8   biological reaction parameter
!D4   cka       REAL*8   biological reaction parameter
!D4   decay     REAL*8    decay rate of biomass
!D4   biofac       REAL*8   biological reaction parameter
!D4   phfac     REAL*8   biological reaction parameter
!D4   phthresh  REAL*8   threshold ph below which biodegredation
!D4                      ceases
!D4   qm        REAL*8   biological reaction parameter
!D4   yield     REAL*8   biological reaction parameter
!D4   ***** COMMON Block concs variables *****
!D4   totaq     REAL*8   total aqueous concentration
!D4   cpnt      REAL*8   free ion concentration
!D4   cplx      REAL*8   complex concentration
!D4   rrcplx    REAL*8   reaction rate associated with the
!D4                      complexes
!D4   rrcpnt    REAL*8   reaction rate associated with the 
!D4                      aqueous component concentrations
!D4   rrimm     REAL*8   reaction rate associated with the
!D4                      immobile component concentrations
!D4   rrvap     REAL*8   reaction rate associated with the
!D4                      vapor component concentrations
!D4   ***** COMMON Block concs variables *****
!D4   drdctaq   REAL*8   the derivative of the reaction rate
!D4                      associated with the aqueous components
!D4                      This array include reaction cross
!D4                      derivative terms between aqueous 
!D4                      components
!D4   drdcimm   REAL*8   the derivative of the reaction rate
!D4                      associated with the immobile components
!D4   drdcvap   REAL*8   the derivative of the reaction rate
!D4                      associated with the vapor components
!D4   drcpnt    REAL*8   the derivative of the reaction rate
!D4                      associated with the aqueous components
!D4   drcplx    REAL*8   the derivative of the reaction rate
!D4                      associated with the aqueous complexes
!D4   drimm     REAL*8   the derivative of the reaction rate
!D4                      associated with the immobile components
!D4   drvap     REAL*8   the derivative of the reaction rate
!D4                      associated with the vapor components
!D4   ***** COMMON Block  jacobian *****
!D4   xjac      REAL*8   jacobian used by speciation solver
!D4   dxct      REAL*8   derivatives used by speciation solver
!D4   resid     REAL*8   residual used by speciation solver
!D4   ipiv         INT   pivots used by speciation solver
!D4   ***** COMMON Block scaling  *******
!D4   ipiv2        INT   pivots used by the scaled speciation
!D4                      solver
!D4   nterms     REAL*8  number of terms used by scaled speciation
!D4   sclmtx     REAL*8  used by scaled speciation
!D4   sclfactr   REAL*8  used by scaled speciation
!D4                    
!D4   ***** COMMON Block rxn_switch *****
!D4   rxnon         INT  nodes at which a kinetic reaction occurs
!D4   numrxn_nodes  INT  # nodes at which a kinetic reaction occurs
!D4   
!**********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!**********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!**********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!**********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!**********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 None
!D9
!**********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!**********************************************************************
!PS
!PS PSEUDOCODE
!PS 
!PS None
!PS
!**********************************************************************

!C     ***** Global constants *****

      integer maxspc
      parameter (maxspc = 40)
      parameter (name_max = 40)

!C ***** Variables in COMMON Block chemnum *****
      integer ncpnt, ph_species 

      integer ncplx

      integer nimm

      integer nvap

      integer numrxn

!C ***** Variables COMMON Block nprints *****
      integer ncpntprt

      integer ncplxprt

      integer nimmprt

      integer nvapprt
      
      integer co2_couple

      integer neg_conc_flag

!C ***** Pointers in COMMON Block idntrxn *****
      integer, allocatable :: naqsp(:) 
      real*8, allocatable :: ps_delta_rxn_s(:)

      integer, allocatable :: nimsp(:)

      integer, allocatable :: nivsp(:)

      integer, allocatable :: irxnic(:,:)

      integer, allocatable :: irxnim(:,:)

      integer, allocatable :: irxniv(:,:)

!C ***** Variables in COMMON Block rdsp *****
      real*8 rsdmax

      real*8, allocatable :: ckeq(:)

      real*8, allocatable :: heq(:,:)

      real*8, allocatable :: spstoic(:,:) 

      real*8, allocatable :: cpntgs(:)

      integer, allocatable :: idrxn(:)

      integer, allocatable :: ifxconc(:)

      integer, allocatable :: neg_conc_possible(:)

!C ***** Variables in COMMON Block initchem *****

      real*8, allocatable :: cpntsv(:,:)

      logical, allocatable :: calcderiv(:)

!C ***** Variables in COMMON Block place *****
      integer, allocatable :: pcpnt(:)

      integer, allocatable :: pimm(:)

      integer, allocatable :: pvap(:)

!C ***** Variables in COMMON Block couple_index *****
      integer drcpos(maxspc,maxspc)

      integer nderivs(maxspc)

      integer matpos(2000)

      integer dimdrc
!C ***** Variables in COMMON Block print_flag *****

      integer, allocatable :: cpntprt(:)

      integer, allocatable :: cplxprt(:)

      integer, allocatable :: immprt(:)

      integer, allocatable :: vapprt(:)

!C ***** Variables in COMMON Block chem_name *****

      character*20, allocatable :: cpntnam(:)

      character*20, allocatable :: cplxnam(:)

      character*20, allocatable :: immnam(:)

      character*20, allocatable :: vapnam(:)
      character(NAME_MAX), allocatable :: species(:)

!C ***** Variables in COMMON Block rxnpara *****

      real*8, allocatable :: ckeqlb(:)

      real*8, allocatable :: ckmtrn(:)

      real*8, allocatable :: simmmx(:)

!C ***** Variables in COMMON Block irrpara *****

      real*8, allocatable :: kfor(:)

      real*8, allocatable :: krev(:)

      real*8, allocatable :: sticirrv(:,:)

      real*8, allocatable :: stimirrv(:,:)

      real*8, allocatable :: stivirrv(:,:)

!C ***** Variables in COMMON Block biopara *****

      real*8, allocatable :: ckc(:)

      real*8, allocatable :: cka(:)

      real*8, allocatable :: decay(:)

      real*8, allocatable :: biofac(:)

      real*8, allocatable :: hfac(:)

      real*8, allocatable :: carbfac(:)

      real*8, allocatable :: ammfac(:)

      real*8, allocatable :: xminit(:)

      integer, allocatable :: nbiofrm(:)

      integer, allocatable :: icbio(:,:) 

      real*8, allocatable :: phthresh(:)

      real*8, allocatable :: qm(:)

      real*8, allocatable :: yield(:)

!C ***** Variables in COMMON Block precdis *****

      real*8, allocatable :: sarea(:)

      real*8, allocatable :: mw_mineral(:), rho_mineral(:)

      real*8, allocatable :: ps_rxn(:,:), ps_delta_rxn(:)

      real*8, allocatable :: pdstic(:,:)

      real*8 pdstim

      integer, allocatable :: pd_flag(:,:)
      

!C ***** Variables in COMMON Block concs *****

      real*8, allocatable :: totaq(:)

      real*8, allocatable :: cplx(:)

      real*8, allocatable :: cpnt(:)

!C ***** Variables in COMMON Block rates *****

      real*8, allocatable :: rrcplx(:)

      real*8, allocatable :: rrcpnt(:)

      real*8, allocatable :: rrimm(:)

      real*8, allocatable :: rrvap(:)

!C ***** Variables in COMMON Block deriv_rates *****

      real*8, allocatable :: drdctaq(:)

      real*8, allocatable :: drdcimm(:,:)

      real*8, allocatable :: drdcvap(:,:)

      real*8, allocatable :: drcpnt(:,:)

      real*8, allocatable :: drtaq(:)

      real*8, allocatable :: drcplx(:)

      real*8, allocatable :: drimm(:)

      real*8, allocatable :: drvap(:)

!C ***** Variables in COMMON Block jacobian *****

      real*8, allocatable :: xjac(:,:)

      real*8, allocatable :: dxct(:,:)

      integer, allocatable :: ipiv(:)

      real*8, allocatable :: resid(:)

!C ***** Variables in COMMON Block scaling *****

      integer, allocatable :: ipiv2(:)

      real*8, allocatable :: nterms(:)

      real*8, allocatable :: sclmtx(:,:)

      real*8, allocatable :: sclfactr(:)

      integer, allocatable :: rxnon(:,:)

      integer, allocatable :: ndconv(:)

      integer iskip

!C ***** Variables in COMMON Block temperature *****
      character*1, allocatable :: temp_model(:)

      character*1, allocatable :: temp_model_kin(:)

      real*8, allocatable :: tcoeff(:,:)

! zvd 26-Jan-07 Added for concentration rate output
      integer cflxz
      integer, allocatable ::  izoncflxz(:),  icflxz(:)
      logical :: cflx_var(5)

! zvd 23-Jul-09 Added for input of total moles
      logical, allocatable :: conc_read(:)

      end module comchem

