      module comdi
!     comdi
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
!D1 Global include file for array variables and pointers (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 28-SEP-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comdi.f_a  $
!D2
!c 12/14/94 gaz added arrays thic,mdnodes
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
!D4   None
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
!D4   ***** COMMON Block fdd pointers and associated variables *****
!D4   ipcpr           pointer  fdd    pointer to array variable cpr    
!D4   ipdeneh         pointer  fdd    pointer to array variable deneh  
!D4   ipdenej         pointer  fdd    pointer to array variable denej  
!D4   ipdenh          pointer  fdd    pointer to array variable denh   
!D4   ipdenj          pointer  fdd    pointer to array variable denj   
!D4   ipdenr          pointer  fdd    pointer to array variable denr   
!D4   ipeflow         pointer  fdd    pointer to array variable eflow  
!D4   ipesk           pointer  fdd    pointer to array variable esk    
!D4   ippcp           pointer  fdd    pointer to array variable pcp    
!D4   ippflow         pointer  fdd    pointer to array variable pflow  
!D4   ipphi           pointer  fdd    pointer to array variable phi    
!D4   ippho           pointer  fdd    pointer to array variable pho    
!D4   ippnx           pointer  fdd    pointer to array variable pnx 
!D4   ippny           pointer  fdd    pointer to array variable pny 
!D4   ippnz           pointer  fdd    pointer to array variable pnz
!D4   ipps            pointer  fdd    pointer to array variable ps     
!D4   ipqflux         pointer  fdd    pointer to array variable qflux  
!D4   ipqflxm         pointer  fdd    pointer to array variable qflxm  
!D4   ipqh            pointer  fdd    pointer to array variable qh     
!D4   ips             pointer  fdd    pointer to array variable s      
!D4   ipsk            pointer  fdd    pointer to array variable sk     
!D4   ipso            pointer  fdd    pointer to array variable so     
!D4   ipt             pointer  fdd    pointer to array variable t      
!D4   ipthx           pointer  fdd    pointer to array variable thx
!D4   ipthy           pointer  fdd    pointer to array variable thy
!D4   ipthz           pointer  fdd    pointer to array variable thz
!D4   ipto            pointer  fdd    pointer to array variable to     
!D4   ipvf            pointer  fdd    pointer to array variable vf     
!D4   ipvolume        pointer  fdd    pointer to array variable volume 
!D4   ipwellim        pointer  fdd    pointer to array variable wellim  
!D4   ipthic          pointer  fdd    pointer to array variable thic    
!D4   ipmdnodes       pointer  fdd    pointer to array variable mdnodes 
!D4   ipmdnode       pointer  fdd    pointer to array variable mdnode 
!D4 
!D4   cpr             REAL*8   fdd    Rock specific heat at each node 
!D4   deneh           REAL*8   fdd    Last time step energy accumulation term
!D4                                     at each node  
!D4   denej           REAL*8   fdd    Last time step energy accumulation time
!D4                                     derivative at each node 
!D4   denh            REAL*8   fdd    Last time step mass accumulation term at
!D4                                     each node 
!D4   denj            REAL*8   fdd    Last time step mass accumulation time
!D4                                     derivative at each node 
!D4   denr            REAL*8   fdd    Rock density at each node 
!D4   eflow           REAL*8   fdd    Energy flow at each source node 
!D4   esk             REAL*8   fdd    Inlet enthalpy associated with a source 
!D4   pcp             REAL*8   fdd    Capillary pressure at each node 
!D4   pflow           REAL*8   fdd    Flowing pressure at each source node 
!D4   phi             REAL*8   fdd    Pressure at each node 
!D4   pho             REAL*8   fdd    Last time step pressure at each node
!D4   pnx             REAL*8   fdd    Permeability in the x-direction, liquid
!D4                                     velocity in the x-direction, vapor
!D4                                     velocity in the x-direction 
!D4   pny             REAL*8   fdd    Permeability in the y-direction liquid
!D4                                     velocity in the y direction, vapor
!D4                                     velocity in the y-direction 
!D4   pnz             REAL*8   fdd    Permeability in the z-direction, liquid
!D4                                     velocity in the z-direction, vapor
!D4                                     velocity in the z-direction 
!D4   ps              REAL*8   fdd    Porosity at each node 
!D4   qflux           REAL*8   fdd    Heat flux at each node 
!D4   qflxm           REAL*8   fdd    Heat flux impedance at each node 
!D4   qh              REAL*8   fdd    Energy source term at each node 
!D4   s               REAL*8   fdd    Liquid saturation at each node 
!D4   sk              REAL*8   fdd    Source strength of each node 
!D4   so              REAL*8   fdd    Last time step saturation at each node 
!D4   t               REAL*8   fdd    Temperature at each node  
!D4   thx             REAL*8   fdd    Thermal conductivity x-direction 
!D4   thy             REAL*8   fdd    Thermal conductivity y-direction 
!D4   thz             REAL*8   fdd    Thermal conductivity z-direction 
!D4   to              REAL*8   fdd    Last time step temperature at each node 
!D4   vf              REAL*8   fdd    Volume factor at each node 
!D4   volume          REAL*8   fdd    Volume associated at each node 
!D4   wellim          REAL*8   fdd    Well impedance at each source node 
!D4 
!D4   ***** COMMON Block fddi pointers and associated variables *****
!D4   ipnskw          pointer  fddi   pointer to array variable nskw 
!D4   ipnskw2         pointer  fddi   pointer to array variable nskw2
!D4
!D4   nskw            INT      fddi   Contains nodes for print-out 
!D4   nskw2           INT      fddi   Contains nodes for print-out 
!D4 
!D4   ***** COMMON Block fdd1 pointers and associated variables *****
!D4   ipa1adfl        pointer  fdd1   pointer to variable array a1adfl
!D4   ipa1adfv        pointer  fdd1   pointer to variable array a1adfv 
!D4   ipa2adfl        pointer  fdd1   pointer to variable array a2adfl 
!D4   ipa2adfv        pointer  fdd1   pointer to variable array a2adfv
!D4   ipan            pointer  fdd1   pointer to variable array an     
!D4   ipanl           pointer  fdd1   pointer to variable array anl    
!D4   ipanlo          pointer  fdd1   pointer to variable array anlo   
!D4   ipanv           pointer  fdd1   pointer to variable array anv    
!D4   ipbetadfl       pointer  fdd1   pointer to variable array betadfl
!D4   ipbetadfv       pointer  fdd1   pointer to variable array betadfv
!D4   ipcm            pointer  fdd1   pointer to variable array cm     
!D4   ipcm0           pointer  fdd1   pointer to variable array cm0    
!D4   ipcnsk          pointer  fdd1   pointer to variable array cnsk   
!D4   ipcp1f          pointer  fdd1   pointer to variable array cp1f   
!D4   ipcp2f          pointer  fdd1   pointer to variable array cp2f   
!D4   ipcp3f          pointer  fdd1   pointer to variable array cp3f   
!D4   ipcp4f          pointer  fdd1   pointer to variable array cp4f   
!D4   ipdench         pointer  fdd1   pointer to variable array dench  
!D4   ipdencj         pointer  fdd1   pointer to variable array dencj  
!D4   ipdiffmfl       pointer  fdd1   pointer to variable array diffmfl
!D4   ipdiffmfv       pointer  fdd1   pointer to variable array diffmfv
!D4   ipdit           pointer  fdd1   pointer to variable array dit    
!D4   ipfc            pointer  fdd1   pointer to variable array fc     
!D4   ipflx12l        pointer  fdd1   pointer to variable array flx12l 
!D4   ipflx12v        pointer  fdd1   pointer to variable array flx12v 
!D4   ipqcin          pointer  fdd1   pointer to variable array qcin   
!D4   ipqcout         pointer  fdd1   pointer to variable array qcout   
!D4   ipqcrxn         pointer  fdd1   pointer to variable array qcrxn  
!D4   iprc            pointer  fdd1   pointer to variable array rc     
!D4   iprcss          pointer  fdd1   pointer to variable array rcss
!D4   iprp1f          pointer  fdd1   pointer to variable array rp1f   
!D4   iprp2f          pointer  fdd1   pointer to variable array rp2f   
!D4   iprp3f          pointer  fdd1   pointer to variable array rp3f   
!D4   iprp4f          pointer  fdd1   pointer to variable array rp4f   
!D4   iprp5f          pointer  fdd1   pointer to variable array rp5f   
!D4   iprp6f          pointer  fdd1   pointer to variable array rp6f   
!D4   iprp7f          pointer  fdd1   pointer to variable array rp7f   
!D4   iprp8f          pointer  fdd1   pointer to variable array rp8f   
!D4   iprp9f          pointer  fdd1   pointer to variable array rp9f   
!D4   iprp10f         pointer  fdd1   pointer to variable array rp10f  
!D4   iprp11f         pointer  fdd1   pointer to variable array rp11f  
!D4   iprp12f         pointer  fdd1   pointer to variable array rp12f  
!D4   iprp13f         pointer  fdd1   pointer to variable array rp13f  
!D4   iprp14f         pointer  fdd1   pointer to variable array rp14f  
!D4   iprp15f         pointer  fdd1   pointer to variable array rp15f  
!D4   iprp16f         pointer  fdd1   pointer to variable array rp16f  
!D4   iprp17f         pointer  fdd1   pointer to variable array rp17f  
!D4   iprp18f         pointer  fdd1   pointer to variable array rp18f  
!D4   iprp19f         pointer  fdd1   pointer to variable array rp19f  
!D4   iprp20f         pointer  fdd1   pointer to variable array rp20f  
!D4   iprp21f         pointer  fdd1   pointer to variable array rp21f  
!D4   iprp22f         pointer  fdd1   pointer to variable array rp22f  
!D4   iprp23f         pointer  fdd1   pointer to variable array rp23f  
!D4   ipt1sk          pointer  fdd1   pointer to variable array t1sk   
!D4   ipt2sk          pointer  fdd1   pointer to variable array t2sk   
!D4   iptclx          pointer  fdd1   pointer to variable array tclx
!D4   iptcly          pointer  fdd1   pointer to variable array tcly
!D4   iptclz          pointer  fdd1   pointer to variable array tclz
!D4   iptcvx          pointer  fdd1   pointer to variable array tcvx
!D4   iptcvy          pointer  fdd1   pointer to variable array tcvy 
!D4   iptcvz          pointer  fdd1   pointer to variable array tcvz
!D4   ipvc1f          pointer  fdd1   pointer to variable array vc1f   
!D4   ipvc2f          pointer  fdd1   pointer to variable array vc2f   
!D4   ipvc3f          pointer  fdd1   pointer to variable array vc3f   
!D4   ipvc4f          pointer  fdd1   pointer to variable array vc4f   
!D4   ipvc5f          pointer  fdd1   pointer to variable array vc5f   
!D4   ipvc6f          pointer  fdd1   pointer to variable array vc6f   
!D4   ipvc7f          pointer  fdd1   pointer to variable array vc7f   
!D4   ipvc8f          pointer  fdd1   pointer to variable array vc8f   
!D4
!D4   a1adfl          REAL*8   fdd1   Alpha1 parameter for each species for
!D4                                     liquid (nonlinear adsorption)  
!D4   a1adfv          REAL*8   fdd1   Alpha1 parameter for each species for
!D4                                     vapor (nonlinear adsorption)  
!D4   a2adfl          REAL*8   fdd1   Alpha2 parameter for each species for
!D4                                     liquid (nonlinear adsorption)  
!D4   a2adfv          REAL*8   fdd1   Alpha2 parameter for each species for
!D4                                     vapor (nonlinear adsorption)  
!D4   an              REAL*8   fdd1   Tracer concentration at each node 
!D4   anl             REAL*8   fdd1   Liquid tracer concentration at each node 
!D4   anlo            REAL*8   fdd1   Last time step tracer concentration at
!D4                                     each node  
!D4   anv             REAL*8   fdd1   Vapor tracer concentration at each node 
!D4   betadfl         REAL*8   fdd1   Beta parameter for each species for
!D4                                     liquid (nonlinear adsorption)  
!D4   betadfv         REAL*8   fdd1   Beta parameter for each species for
!D4                                     vapor (nonlinear adsorption)  
!D4   cm              REAL*8   fdd1   Total tracer mass for each species 
!D4   cm0             REAL*8   fdd1   Initial total tracer mass for each
!D4                                     species 
!D4   cnsk            REAL*8   fdd1   Tracer concentration source term at each
!D4                                     node 
!D4   cp1f            REAL*8   fdd1   Parameter in capillary pressure model 
!D4   cp2f            REAL*8   fdd1   Parameter in capillary pressure model 
!D4   cp3f            REAL*8   fdd1   Parameter in capillary pressure model 
!D4   cp4f            REAL*8   fdd1   Parameter in capillary pressure model 
!D4   dench           REAL*8   fdd1   Last time step tracer accumulation term
!D4                                     at each node
!D4   dencj           REAL*8   fdd1   Last time step tracer accumulation
!D4                                     derivative term at each node
!D4   diffmfl         REAL*8   fdd1   Liquid molecular diffusion coefficient
!D4   diffmfv         REAL*8   fdd1   Vapor molecular diffusion coefficient
!D4   dit             REAL*8   fdd1   Array containing time step changes  
!D4   fc              REAL*8   fdd1   Tracer equation residual at each node 
!D4   flx12l          REAL*8   fdd1   ?
!D4   flx12v          REAL*8   fdd1   ?
!D4   qcin            REAL*8   fdd1   Total injected tracer mass for each
!D4                                     species 
!D4   qcout           REAL*8   fdd1   Total produced tracer mass for each
!D4                                     species 
!D4   qcrxn           REAL*8   fdd1   Total net tracer mass produced by
!D4                                     reaction for each species 
!D4   rc              REAL*8   fdd1   Tracer source/reaction term at each node 
!D4   rcss            REAL*8   fdd1   Tracer source term at each node 
!D4   rp1f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp2f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp3f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp4f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp5f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp6f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp7f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp8f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp9f            REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp10f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp11f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp12f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp13f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp14f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp15f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp16f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp17f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp18f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp19f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp20f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp21f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp22f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp23f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   rp24f           REAL*8   fdd1   Parameter in relative permeability model 
!D4   t1sk            REAL*8   fdd1   Time when the tracer source term is
!D4                                     activated at each node  
!D4   t2sk            REAL*8   fdd1   Time when the tracer source term is
!D4                                     terminated at each node  
!D4   tclx            REAL*8   fdd1   Tracer liquid dispersion coefficient in
!D4                                     the x-direction 
!D4   tcly            REAL*8   fdd1   Tracer liquid dispersion coefficient in
!D4                                     the y-direction 
!D4   tclz            REAL*8   fdd1   Tracer liquid dispersion coefficient in
!D4                                     the z-direction 
!D4   tcvx            REAL*8   fdd1   Tracer vapor dispersion coefficient in
!D4                                     the x-direction 
!D4   tcvy            REAL*8   fdd1   Tracer vapor dispersion coefficient in
!D4                                     the y-direction 
!D4   tcvz            REAL*8   fdd1   Tracer vapor dispersion coefficient in
!D4                                     the z-direction 
!D4   vc1f            REAL*8   fdd1   ? 
!D4   vc2f            REAL*8   fdd1   ?
!D4   vc3f            REAL*8   fdd1   ? 
!D4   vc4f            REAL*8   fdd1   ? 
!D4   vc5f            REAL*8   fdd1   ? 
!D4   vc6f            REAL*8   fdd1   ? 
!D4   vc7f            REAL*8   fdd1   ? 
!D4   vc8f            REAL*8   fdd1   ? 
!D4
!D4   ***** COMMON Block fdd1i pointers and associated variables *****
!D4   ipiadd          pointer  fdd1i  pointer to variable array iadd   
!D4   ipiaddt         pointer  fdd1i  pointer to variable array iaddt  
!D4   ipiadsfl        pointer  fdd1i  pointer to variable array iadsfl 
!D4   ipiadsfv        pointer  fdd1i  pointer to variable array iadsfv 
!D4   ipicapt         pointer  fdd1i  pointer to variable array icapt  
!D4   ipicns          pointer  fdd1i  pointer to variable array icns   
!D4   ipiflx1         pointer  fdd1i  pointer to variable array iflx1  
!D4   ipiflx2         pointer  fdd1i  pointer to variable array iflx2  
!D4   ipirlpt         pointer  fdd1i  pointer to variable array irlpt  
!D4   ipitc           pointer  fdd1i  pointer to variable array itc    
!D4   ipivcn          pointer  fdd1i  pointer to variable array ivcn   
!D4   ipivcon         pointer  fdd1i  pointer to variable array ivcon  
!D4   ipnpt           pointer  fdd1i  pointer to variable array npt    
!D4
!D4   iadd            INT      fdd1i  Iteration count for tracer solution 
!D4   iaddt           INT      fdd1i  Iteration count used for linear equation
!D4                                     solver in tracer solution 
!D4   iadsfl          INT      fdd1i  Adsorption model flag used in tracer
!D4                                     solution for liquid
!D4   iadsfv          INT      fdd1i  Adsorption model flag used in tracer
!D4                                     solution for vapor
!D4   icapt           INT      fdd1i  Capillary pressure model 
!D4   icns            INT      fdd1i  Flag denoting whether a liquid, vapor,
!D4                                     or Henry's law species
!D4   iflx1           INT      fdd1i  ?
!D4   iflx2           INT      fdd1i  ?
!D4   irlpt           INT      fdd1i  Relative permeability model 
!D4   itc             INT      fdd1i  Array containing information used in
!D4                                     time step changes 
!D4   ivcn            INT      fdd1i  ?
!D4   ivcon           INT      fdd1i  ?
!D4   npn             INT      fdd1i  Parameter used in storing tracer data 
!D4   npt             INT      fdd1i  Storage parameter used in tracer solution
!D4   nsp             INT      fdd1i  Current species number 
!D4   nspeci          INT      fdd1i  Number of species for tracer solution 
!D4
!D4   ***** !OMMON Block fdd2 pointers and associated variables *****
!D4   ipamgang        pointer  fdd2   pointer to variable array amgang
!D4   ipdporp         pointer  fdd2   pointer to variable array dporp
!D4   ipdport         pointer  fdd2   pointer to variable array dport
!D4   ippgangi        pointer  fdd2   pointer to variable array pgangi
!D4   ipphini         pointer  fdd2   pointer to variable array phini
!D4   ippsini         pointer  fdd2   pointer to variable array psini
!D4   iptini          pointer  fdd2   pointer to variable array tini
!D4
!D4   amgang          REAL*8   fdd2   Parameter used in gangi model at each
!D4                                     node
!D4   dporp           REAL*8   fdd2   Derivative of porosity with respect to
!D4                                     pressure at each node 
!D4   dport           REAL*8   fdd2   Derivative of porosity with respect to
!D4                                     temperature at each node 
!D4   pgangi          REAL*8   fdd2   Parameter used in gangi model at each
!D4                                     node
!D4   phini           REAL*8   fdd2   Initial pressure at each node 
!D4   psini           REAL*8   fdd2   Initial porosity at each node 
!D4   psdelta         REAL*8   fdd2   Total change in porosity
!D4   psvol           REAL*8   fdd2   Volume of nodes with ppor model
!D4                                     that allows porosity change
!D4   tini            REAL*8   fdd2   Initial temperature at each node 
!D4
!D4   porTemp1--4     REAL*8   fdd2   Parameters used to define the temperature
!D4                                   -dependent porosity at each node
!D4
!D4   ***** COMMON Block fdd2r *****
!D4   phydro          REAL*8   fdd2r  Parameter used in rock deformation model 
!D4   sigini          REAL*8   fdd2r  Parameter used in rock deformation model 
!D4   thexp           REAL*8   fdd2r  Parameter used in rock deformation model 
!D4   young           REAL*8   fdd2r  Parameter used in rock deformation model 
!D4
!D4   ***** COMMON Block fice pointers and associated variables *****
!D4   ipsii           pointer  fice   pointer to variable array sii
!D4   ipsio           pointer  fice   pointer to variable array sio
!D4
!D4   sii             REAL*8   fice   Ice saturation at each node 
!D4   sio             REAL*8   fice   Last time step value of ice saturation 
!D4
!D4   ***** COMMON Block ficer *****
!D4   tmelt           REAL*8   ficer  Freezing temperature of water 
!D4
!D4   ***** COMMON Block iice pointer and associated variable *****
!D4   ipices          pointer  iice   pointer to variable array ices
!D4
!D4   ices            INT      iice   State of ice at each node 
!D4
!D4   ***** COMMON Block fddi1 pointers and associated variables *****
!D4   ipicap          pointer  fddi1  pointer to array variable icap
!D4   ipieos          pointer  fddi1  pointer to array variable ieos
!D4   ipiieos         pointer  fddi1  pointer to array variable iieos
!D4   ipiporf         pointer  fddi1  pointer to array variable iporf
!D4   ipirlp          pointer  fddi1  pointer to array variable irlp
!D4
!D4   icap            INT      fddi1  Capillary pressure model at each node  
!D4   ieos            INT      fddi1  Phase state of fluid at each node 
!D4   ieosc           INT      fddi1  Phase change parameter 
!D4   iieos           INT      fddi1  Thermodynamics set at each node 
!D4   iporf           INT      fddi1  Deformation model at each node 
!D4   irlp            INT      fddi1  Relative permeability model at each node 
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

      integer ieosc, npn, nsp, nspeci, numd, numsorp, numvcon
      real*8 phydro, sigini, thexp, tmelt, young

c KCL 5-2-11, for DeltaPoros subroutine
      real*8, allocatable :: MaxDrainPart(:)
      real*8, allocatable :: InitPoros(:)
      real*8, allocatable :: InitPerm(:)
      real*8, allocatable :: InitPres(:)
      integer :: DrainInitFlag

!     ***** pointers in COMMON Block fdd *****
      real*8, allocatable ::  cpr(:)
      real*8, allocatable ::  deneh(:)
      real*8, allocatable ::  denej(:) 
      real*8, allocatable ::  denh(:)
      real*8, allocatable ::  denj(:) 
      real*8, allocatable ::  denr(:) 
      real*8, allocatable ::  eflow(:) 
      real*8, allocatable ::  esk(:) 
      real*8, allocatable ::  head(:) 
      real*8, allocatable ::  pcp(:)
      real*8, allocatable ::  pflow(:) 
      real*8, allocatable ::  pflowa(:)
      real*8, allocatable ::  xairfl(:)
      real*8, allocatable ::  phi(:) 
      real*8, allocatable ::  pho(:)
      real*8, allocatable ::  pnx(:) 
      real*8, allocatable ::  pny(:)
      real*8, allocatable ::  pnz(:)
      real*8, allocatable ::  pnxi(:) 
      real*8, allocatable ::  pnyi(:)
      real*8, allocatable ::  pnzi(:)
c gaz 090113
      real*8, allocatable ::  pnx_old(:)
      real*8, allocatable ::  pnx_save(:)
      real*8, allocatable ::  pny_save(:)
      real*8, allocatable ::  pnz_save(:)
      real*8, allocatable ::  anxy(:) 
      real*8, allocatable ::  anxz(:) 
      real*8, allocatable ::  anyz(:) 
!     GAZ 3-19-00 rlp_fac residual rel perm
      real*8, allocatable ::  rlp_fac(:)
!     ZVD 09-Apr-2010 fperm for use with rlp
      real*8, allocatable ::  xfperm(:)
      real*8, allocatable ::  yfperm(:)
      real*8, allocatable ::  zfperm(:)
      real*8, allocatable ::  ps(:)
      real*8, allocatable ::  ps_trac(:)
      real*8, allocatable ::  qflux(:) 
      real*8, allocatable ::  qflxm(:) 
      real*8, allocatable ::  qh(:) 
      real*8, allocatable ::  s(:) 
      real*8, allocatable ::  sk(:)
c gaz  071819 added sko(:)  last time step source/sink    
      real*8, allocatable ::  sko(:)
      real*8, allocatable ::  recharge(:)
      real*8, allocatable ::  so(:) 
      real*8, allocatable ::  t(:)
      real*8, allocatable ::  thx(:) 
      real*8, allocatable ::  thy(:) 
      real*8, allocatable ::  thz(:) 
      real*8, allocatable ::  to(:) 
      real*8, allocatable ::  vf(:) 
      real*8, allocatable ::  volume(:) 
      real*8, allocatable ::  wellim(:) 
      real*8, allocatable ::  wellima(:) 
      real*8, allocatable ::  thic(:) 
      real*8, allocatable ::  strgan(:) 
      real*8, allocatable ::  alphae(:) 
      real*8, allocatable ::  var1(:) 
      real*8, allocatable ::  var2(:) 
      real*8, allocatable ::  dvar11(:) 
      real*8, allocatable ::  dvar12(:) 
      real*8, allocatable ::  dvar21(:) 
      real*8, allocatable ::  dvar22(:) 
      integer, allocatable ::  mdnodes(:) 
      integer, allocatable ::  mdnode(:,:) 
      integer, allocatable ::  npest(:) 

!     ***** pointers in COMMON Block fddi *****
      integer, allocatable :: nskw(:) 
      integer, allocatable :: nskw2(:)
!     ****   TENMA   ****
      integer, allocatable :: nskw3(:)
!     ****   TENMA   ****

c pjjohnson change for Leverett Pc function
      real*8, allocatable :: pjki(:)
      real*8, allocatable :: pjk(:) 

c end pjjohnson edit      
      
!     ***** pointers in COMMON Block fdd1 *****
      real*8, allocatable ::  a1adfl(:,:) 
      real*8, allocatable ::  a1adfv(:,:) 
      real*8, allocatable ::  a2adfl(:,:) 
      real*8, allocatable ::  a2adfv(:,:) 
      real*8, allocatable ::  an(:) 
      real*8, allocatable ::  anl(:) 
      real*8, allocatable ::  anlo(:)
      real*8, allocatable ::  anv(:)
      real*8, allocatable ::  betadfl(:,:)
      real*8, allocatable ::  betadfv(:,:)
      real*8, allocatable ::  cm(:) 
      real*8, allocatable ::  cm0(:)
      real*8, allocatable ::  cnsk(:)
      real*8, allocatable ::  cnsk_background(:)
      real*8, allocatable ::  cp1f(:) 
      real*8, allocatable ::  cp2f(:) 
      real*8, allocatable ::  cp3f(:) 
      real*8, allocatable ::  cp4f(:)
      real*8, allocatable ::  dench(:) 
      real*8, allocatable ::  dencj(:) 
      real*8, allocatable ::  diffmfl(:,:) 
      real*8, allocatable ::  diffmfv(:,:) 
      real*8, allocatable ::  dit(:)
      real*8, allocatable ::  fc(:) 
      real*8, allocatable ::  flx12l(:) 
      real*8, allocatable ::  flx12v(:)
      real*8, allocatable ::  qcin(:) 
      real*8, allocatable ::  qcout(:)
      real*8, allocatable ::  qcrxn(:)
      real*8, allocatable ::  rc(:)
      real*8, allocatable ::  rcss(:)
      real*8, allocatable ::  rp1f(:) 
      real*8, allocatable ::  rp2f(:) 
      real*8, allocatable ::  rp3f(:)
      real*8, allocatable ::  rp4f(:) 
      real*8, allocatable ::  rp5f(:) 
      real*8, allocatable ::  rp6f(:) 
      real*8, allocatable ::  rp7f(:) 
      real*8, allocatable ::  rp8f(:)
      real*8, allocatable ::  rp9f(:) 
      real*8, allocatable ::  rp10f(:) 
      real*8, allocatable ::  rp11f(:) 
      real*8, allocatable ::  rp12f(:)
      real*8, allocatable ::  rp13f(:) 
      real*8, allocatable ::  rp14f(:) 
      real*8, allocatable ::  rp15f(:) 
      real*8, allocatable ::  rp16f(:)
      real*8, allocatable ::  rp17f(:) 
      real*8, allocatable ::  rp18f(:) 
      real*8, allocatable ::  rp19f(:) 
      real*8, allocatable ::  rp20f(:) 
      real*8, allocatable ::  rp21f(:)
      real*8, allocatable ::  rp22f(:) 
      real*8, allocatable ::  rp23f(:) 
      real*8, allocatable ::  rp24f(:) 
      real*8, allocatable ::  fmvf(:) 
      real*8, allocatable ::  dfmvf(:) 
      real*8, allocatable ::  fmlf(:) 
      real*8, allocatable ::  dfmlf(:) 
      real*8, allocatable ::  t1sk(:)
      real*8, allocatable ::  t2sk(:) 
      real*8, allocatable ::  tclx(:,:) 
      real*8, allocatable ::  tcly(:,:) 
      real*8, allocatable ::  tclz(:,:)
      real*8, allocatable ::  tcvx(:,:) 
      real*8, allocatable ::  tcvy(:,:) 
      real*8, allocatable ::  tcvz(:,:)
      real*8, allocatable ::  vc1f(:) 
      real*8, allocatable ::  vc2f(:) 
      real*8, allocatable ::  vc3f(:) 
      real*8, allocatable ::  vc4f(:) 
      real*8, allocatable ::  vc5f(:) 
      real*8, allocatable ::  vc6f(:) 
      real*8, allocatable ::  vc7f(:) 
      real*8, allocatable ::  vc8f(:) 
      real*8, allocatable ::  k0f(:) 
      real*8, allocatable ::  bkf(:) 
      real*8, allocatable ::  por0f(:) 

!     ***** pointers in COMMON Block fdd1i *****
      integer, allocatable :: iadd(:) 
      integer, allocatable :: iaddt(:)
      integer, allocatable :: iadsfl(:,:) 
      integer, allocatable :: iadsfv(:,:) 
      integer, allocatable :: icapt(:) 
      integer, allocatable :: icns(:)
      integer, allocatable :: iflx1(:)
      integer, allocatable :: iflx2(:)  
      integer, allocatable :: iflxz(:)  
      integer, allocatable :: irlpt(:) 
      integer, allocatable :: itc(:) 
      integer, allocatable :: ivcn(:) 
      integer, allocatable :: ivcon(:) 
      integer, allocatable :: ivbounf(:) 
      integer, allocatable :: mflagl(:,:)
      integer, allocatable :: mflagv(:,:)
      integer, allocatable :: npt(:) 
 
!     ***** pointers in COMMON Block fdd2 *****
      real*8, allocatable ::  amgang(:)
      real*8, allocatable ::  dporp(:) 
      real*8, allocatable ::  dport(:)
      real*8, allocatable ::  pgangi(:)
!     temperature-dependent porosity
      real*8, allocatable ::  porTemp1(:)
      real*8, allocatable ::  porTemp2(:)
      real*8, allocatable ::  porTemp3(:)
      real*8, allocatable ::  porTemp4(:)
!     ****   TENMA   ****
      real*8, allocatable ::  wgangi(:)
      real*8, allocatable ::  sgangi(:)
      real*8, allocatable ::  agangi(:)
      real*8, allocatable ::  tenma_ww(:)
!     ****   TENMA   ****
      real*8, allocatable ::  phini(:) 
      real*8, allocatable ::  psini(:) 
      real*8 psdelta, psvol 
      real*8, allocatable ::  tini(:) 
c gaz 090113 array for last TS porosity
      real*8, allocatable ::  ps_old(:)
      real*8, allocatable ::  ps_save(:)
 
!     ***** pointers in COMMON Block fice *****
      real*8, allocatable ::  sii(:) 
      real*8, allocatable ::  sio(:)
 
!     ***** pointers in COMMON Block iice *****
      integer, allocatable :: ices(:)
 
!     ***** pointers in COMMON Block fddi1 *****
      integer, allocatable :: icap(:) 
      integer, allocatable :: ieos(:) 
      integer, allocatable :: iieos(:) 
      integer, allocatable :: iporf(:) 
      integer, allocatable :: irlp(:)
      integer, allocatable :: itrc(:)
      integer, allocatable :: itrcdsp(:)
      integer ldsp

!     ***** pointers in flow_boundary_conditions *****
      integer maxmodel,maxtimes,modmin,mmodel
      integer, allocatable ::  time_type(:)
      integer, allocatable ::  time_interpolate(:)
      integer, allocatable ::  tunit_type(:)
      integer, allocatable ::  sourcea_type(:)
      integer, allocatable ::  sourcew_type(:)
      integer, allocatable ::  sourcef_type(:)
      integer, allocatable ::  sourcee_type(:)
      integer, allocatable ::  sourceco2_type(:)
      integer, allocatable ::  esourceco2_type(:)
      integer, allocatable ::  seepfac_type(:)
      integer, allocatable ::  drainar_type(:)
      integer, allocatable ::  humid_type(:)
      integer, allocatable ::  phumid_type(:)
      integer, allocatable ::  thumid_type(:)
      integer, allocatable ::  saturation_type(:)
      integer, allocatable ::  pressurew_type(:)
      integer, allocatable ::  pressurea_type(:)
      integer, allocatable ::  impedance_type(:)
      integer, allocatable ::  enthalpy_type(:)
      integer, allocatable ::  timestep_type(:)
      integer, allocatable ::  eqtol_type(:)
      integer, allocatable ::  tsmax_type(:)
      integer, allocatable ::  temperature_type(:)
      integer, allocatable ::  node_model(:)
      integer, allocatable ::  min_model(:)
      integer, allocatable ::  node_ch_model_type(:)
      integer, allocatable ::  steady_type(:)
      integer, allocatable ::  weight_type(:)
      integer, allocatable ::  permx_type(:)
      integer, allocatable ::  permy_type(:)
      integer, allocatable ::  permz_type(:)

      integer, allocatable ::  saturation_ini_type(:)
      integer, allocatable ::  pressurew_ini_type(:)
      integer, allocatable ::  pressurea_ini_type(:)
      integer, allocatable ::  temperature_ini_type(:)

      integer iqa,ixa,iqw,iqenth,isatb,ienth,itempb,itempb2,ipresa
      integer imped,its,imod,isubmod,iqf,icm,isf,ifd,iqco2,ipresw
      integer isatb_ini,ipresa_ini,ipresw_ini,itempb_ini
      integer iha,ipha,itha,itsmax,ieqtol
      integer isteady, istdy, isty, iwght, lchange
      integer ixperm,iyperm,izperm
      real*8 fac_sec_days, fac_min_days, fac_year_days
      real*8, allocatable ::  qa(:) 
c gaz 111418 added qaxf for air fraction in water      
      real*8, allocatable ::  qaxf(:)       
      real*8, allocatable ::  qw(:) 
      real*8, allocatable ::  qw0(:) 
      real*8, allocatable ::  qco2b(:) 
      real*8, allocatable ::  eco2b(:)
      real*8, allocatable ::  qenth(:) 
      real*8, allocatable ::  satb(:) 
      real*8, allocatable ::  enth(:) 
      real*8, allocatable ::  tempb(:) 
      real*8, allocatable ::  presa(:) 
      real*8, allocatable ::  presw(:) 
      real*8, allocatable ::  huma(:) 
      real*8, allocatable ::  phuma(:) 
      real*8, allocatable ::  thuma(:) 
c gaz debug  033021    
      real*8, allocatable ::  pchuma(:) 
      real*8, allocatable ::  xnva(:) 
      real*8, allocatable ::  entha(:) 
      real*8, allocatable ::  sp(:) 
      real*8, allocatable ::  drain(:)

      real*8, allocatable ::  humida(:) 
      real*8, allocatable ::  phumida(:)
      real*8, allocatable ::  thumida(:)
      real*8, allocatable ::  enth_humid(:)  

      real*8, allocatable ::  satb_ini(:) 
      real*8, allocatable ::  presw_ini(:) 
      real*8, allocatable ::  presa_ini(:) 
      real*8, allocatable ::  tempb_ini(:) 
       
      real*8, allocatable :: time(:,:)
      real*8, allocatable :: time_cycle(:)
      real*8, allocatable :: sourcea(:,:)
      real*8, allocatable :: sourcew(:,:)
      real*8, allocatable :: sourcef(:,:)
      real*8, allocatable :: sourcee(:,:)
      real*8, allocatable :: sourceco2(:,:)
      real*8, allocatable :: esourceco2(:,:)
      real*8, allocatable :: seepfac(:,:)
      real*8, allocatable :: drainar(:,:)
      real*8, allocatable :: saturation(:,:)
      real*8, allocatable :: pressurew(:,:)
      real*8, allocatable :: pressurea(:,:)
      real*8, allocatable :: impedance(:,:)
      real*8, allocatable :: enthalpy(:,:)
      real*8, allocatable :: timestep(:,:)
      real*8, allocatable :: timestepmax(:,:)
      real*8, allocatable :: eqtolerance(:,:)
      real*8, allocatable :: temperature(:,:)
c gaz 123115
      real*8, allocatable ::  humid(:,:) 
      real*8, allocatable ::  phumid(:,:)
      real*8, allocatable ::  thumid(:,:)
       
      real*8, allocatable :: permx(:,:)
      real*8, allocatable :: permy(:,:)
      real*8, allocatable :: permz(:,:)
      real*8, allocatable :: vtotw(:),vtota(:),vtote(:),vtotco2(:)
      real*8, allocatable :: atotd(:)
      real*8, allocatable :: steady_time(:)
      integer, allocatable ::  node_ch_model(:,:)

      real*8, allocatable :: saturation_ini(:,:)
      real*8, allocatable :: pressurew_ini(:,:)
      real*8, allocatable :: pressurea_ini(:,:)
      real*8, allocatable :: temperature_ini(:,:)

      integer nsurf,isurf
      real*8 head0, temp0, pres0, rol0, headconv
      integer, allocatable :: izone_surf(:)
      integer, allocatable :: izone_surf_nodes(:)
      real*8, allocatable :: headconv_val(:)
      real*8, allocatable :: headconv_val_temp(:)

      integer nfree, ifree, ifree1, isw_term
      real*8 head_tol, head_ck
      integer, allocatable :: ifreef(:)  
      integer, allocatable :: ifree_im_ex(:)  
      integer, allocatable :: izone_free_nodes(:)  
      integer, allocatable :: izone_free(:)  
      real*8, allocatable :: head12(:,:)
      real*8, allocatable :: t91(:)
      real*8, allocatable :: dfidf(:)
      real*8, allocatable :: dfid1f(:)
      real*8, allocatable :: rlxyf(:)
      real*8, allocatable :: drlxyf(:)
      real*8, allocatable :: rlzf(:)
      real*8, allocatable :: drlzf(:)
      real*8, allocatable :: dzrg(:)
      real*8, allocatable :: dzrg_new(:)      
      real*8, allocatable :: dxrg(:)
      real*8, allocatable :: dyrg(:)
 
      integer ngrad, igrad, igrad0
      integer, allocatable :: igradf_temp(:)  
      integer, allocatable :: idirg_temp(:)  
      integer, allocatable :: izone_grad_nodes_temp(:)  
      integer, allocatable :: izone_grad_temp(:)  
      real*8, allocatable :: cordg_temp(:)
      real*8, allocatable :: grad1_temp(:)
      real*8, allocatable :: var0_temp(:)
      integer, allocatable :: igradf(:)  
      integer, allocatable :: idirg(:)  
      integer, allocatable :: izone_grad_nodes(:)  
      integer, allocatable :: izone_grad(:)  
      real*8, allocatable :: cordg(:)
      real*8, allocatable :: grad1(:)
      real*8, allocatable :: var0(:)
      character*80, allocatable :: gradmod_filename(:)
      integer, allocatable :: igradmodelfile(:)
      integer, allocatable :: igradmodnamlen(:)
 
      integer igradmd,idgradmc

      integer nconv, iconv, nconv0, nall
      integer, allocatable :: isconv(:)  
      integer, allocatable :: ifconv(:)  
      integer, allocatable :: iconvf_temp(:)  
      integer, allocatable :: idirc_temp(:)  
      integer, allocatable :: izone_conv_nodes_temp(:)  
      integer, allocatable :: izone_conv_temp(:)  
      real*8, allocatable :: cordc_temp(:)
      real*8, allocatable :: conv1_temp(:)
      real*8, allocatable :: conv2_temp(:)
      real*8, allocatable :: conv3_temp(:)
      real*8, allocatable :: varc_temp(:)
      integer, allocatable :: iconvf(:)  
      integer, allocatable :: idirc(:)  
      integer, allocatable :: izone_conv_nodes(:)  
      integer, allocatable :: izone_conv(:)  
      real*8, allocatable :: cordc(:)
      real*8, allocatable :: conv1(:)
      real*8, allocatable :: conv2(:)
      real*8, allocatable :: conv3(:)
      real*8, allocatable :: varc(:)

      integer narea, iarea, narea0
      integer, allocatable :: isarea(:)  
      integer, allocatable :: ifarea(:)  
      integer, allocatable :: iareaf_temp(:)  
      integer, allocatable :: izone_area_nodes_temp(:)  
      integer, allocatable :: izone_area_temp(:)  
      real*8, allocatable :: area01_temp(:)
      real*8, allocatable :: area02_temp(:)
      integer, allocatable :: iareaf(:)  
      integer, allocatable :: izone_area_nodes(:)  
      integer, allocatable :: izone_area(:)  
      real*8, allocatable :: area01(:)
      real*8, allocatable :: area02(:)
      real*8, allocatable :: wgt_area(:)
      real*8, allocatable :: wgt_impd(:)


      integer nwell, iwell, nwell0, nallw
      integer, allocatable :: iswell(:)  
      integer, allocatable :: ifwell(:)  
      integer, allocatable :: iwellf_temp(:)  
      integer, allocatable :: iwellconf_temp(:)  
      integer, allocatable :: izone_well_nodes_temp(:,:)  
      integer, allocatable :: izone_well_temp(:)  
      real*8, allocatable :: cordw_temp(:,:)
      real*8, allocatable :: twell1_temp(:)
      real*8, allocatable :: twell2_temp(:)
      integer, allocatable :: iwellf(:)  
      integer, allocatable :: iwellconf(:)  
      integer, allocatable :: izone_well_nodes(:,:)  
      integer, allocatable :: izone_well(:)  
      real*8, allocatable :: cordw(:,:)
      real*8, allocatable :: twell1(:)
      real*8, allocatable :: twell2(:)

      real*8, allocatable :: totboun(:)

      integer ngh
      real*8, allocatable ::  pflow_gh(:) 
      real*8, allocatable ::  wellim_gh(:) 
      integer, allocatable :: node_gh(:)  
      integer, allocatable :: idir_gh(:)  

      integer iwt_uz
      real*8, allocatable :: dhead12s(:)

      integer ibp_maxw,ibp_maxa,ibpmaxp,ibpmaxs
      real*8 bp_maxw,bp_maxa,bp_maxp,bp_maxs

      real*8, allocatable :: sinkint(:)
 
      integer, allocatable :: ieos_ch(:)
      real*8 time_ch
      real*8, allocatable :: time_ieos(:) 

! gaz 11-Jan-08 Pressure from mass conversion
      real*8, allocatable :: mass_var(:)
c  phase state variable for simplified thermodynamics     
      integer, allocatable :: ieos_aux(:) 
      integer, allocatable :: i_chk(:) 
      integer, allocatable :: izone_bot(:)
      integer, allocatable :: izone_top(:)
      integer, allocatable :: izone_drain(:)
      integer, allocatable :: izone_main(:)
      real*8, allocatable :: z_plot(:)   
      real*8, allocatable :: dzrg_sub(:)
      real*8, allocatable :: z_plot_old(:)
c some simple thermo       
      real*8 LiqEndTemp, VapEndTemp, IceEndTemp   
      
      integer ivol_cnt, iarea_cnt,max_replace
      parameter (max_replace = 1000000)
      real*8,  allocatable :: vol_pri(:)
      real*8,  allocatable :: vol_sec(:)
      real*8,  allocatable :: area_pri(:)
      real*8,  allocatable :: area_sec(:)
      integer,  allocatable :: ii_vol(:)
      integer,  allocatable :: icon_area1(:)
      integer,  allocatable :: icon_area2(:)
c gaz 090113 
      real*8 por_salt_min
      real*8 pressure_std,temperature_std
      parameter(pressure_std = 0.1,temperature_std = 20.0)
c gaz 080817 
      real*8 energy_conv
      parameter(energy_conv = 1.d-6)
c gaz 071819
      real*8 source_ratio_out, source_ratio_in
      integer ntable_roc
      character*200, allocatable :: table_vroc(:)
      integer, allocatable :: ivrn(:) 
      integer, allocatable :: ivroc(:) 
      integer, allocatable :: itroc(:) 
      integer, allocatable :: ivrov(:)
      integer, allocatable :: ntable_vroc(:)
      integer, allocatable :: tblindx_roc(:,:)
      real*8, allocatable ::  roc_table(:,:)       
      real*8, allocatable ::  temp_table(:,:) 
      
      real*8, allocatable ::  vroc1f(:) 
      real*8, allocatable ::  vroc2f(:) 
      real*8, allocatable ::  vroc3f(:) 
      real*8, allocatable ::  vroc4f(:) 
      real*8, allocatable ::  vroc5f(:) 
      real*8, allocatable ::  vroc6f(:) 
      real*8, allocatable ::  vroc7f(:) 
      real*8, allocatable ::  vroc8f(:) 
      real*8, allocatable ::  vroc9f(:)      
      real*8, allocatable ::  dcprt(:) 
      real*8, allocatable ::  ddenrt(:) 
      
      real*8, allocatable ::  urock(:)
      real*8, allocatable ::  durockt(:)
c gaz 111118   
c sk_temp() used in sub thrmwc      
      real*8, allocatable ::  sk_temp(:) 
c gaz 112920  
c sk0() used in sub thrmwc      
      real*8, allocatable ::  sk0(:)       
c gaz 081918
      integer, allocatable :: izone_renum(:)
      integer, allocatable :: nrenu_list(:)
      integer, allocatable :: irb_renum_act(:)
      integer, allocatable :: ncon_renum(:) 
      integer n_renu_zone, min_r_z, max_r_z, neq_act, neq_nonact
c gaz 090819
c arrays and global variables associated with mass equality with
c  phase change
      real*8, allocatable ::  phi_prev(:)
      real*8, allocatable ::  s_prev(:)
      integer, allocatable :: n_phase_nodes(:)
      integer, allocatable :: ieos_prev(:)
      integer n_phase_ch
c gaz 072120 array to track sc phase and T,P
c ieos_sc(i) = 0, no T or P ge than Critical value
c ieos_sc(i) = 0, T and P ge than Critical value      
c ieos_sc(i) = 2, T ge than Critical value
c ieos_sc(i) = 3, P ge than Critical value      
      integer, allocatable :: ieos_sc(:)     
      end module comdi
