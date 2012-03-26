subroutine readtracrxn
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read input for trxn macro.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 Initial implementation: 2011,  Programmer: M. Schauer
!D2
!***********************************************************************

  use combi
  use comchem
  use comrxni
  use comdi
  use comcouple
  use comdti
  use comai
  use comflow
  use compart
  implicit none
  integer NAME_MAX, SPEC_MAX, LINE_MAX, INT_LEN			! Parameters
  character*4 aINT_LEN, aSPEC_MAX, aNAME_MAX			! ...and their character counterparts for use in format strings
  parameter ( NAME_MAX = 20 )
  parameter ( SPEC_MAX = 100 )
  parameter ( LINE_MAX = 200 )
  parameter ( INT_LEN = 3 )
  parameter ( aSPEC_MAX = '100' )	
  parameter ( aNAME_MAX = '20' )
  parameter ( aINT_LEN = '3' )
  character*10 keyword						! Keyword of the current block being read
  character(LINE_MAX) line, line2, null, e, f, g			! Current line and generic variables
  integer linen							! Number of the current line
  integer i, j, k, o, a, c, d					! Generic variables
  real*8 r
  logical flag, flag2, reading					! Generic flags, whether we are reading input
  logical fcomp, fwater, frock, fsorp, fdisp, fvap, fheader, fassign, fhparam, fhenry, fdiff, fspec, &
       fgroup, fprint, fsol, fdist, fequi			! Flags denoting whether each section has been specified.  matrix and lookup are taken care of by stoichmethod.
  integer stoichmethod						! 1 for lookup, 2 for matrix, 0 if nothing has been specified.
  integer nreactions						! The number of reactions that have been specified.
  character(NAME_MAX), dimension(SPEC_MAX) :: array, array2	! Generic, temporary
  character(NAME_MAX), dimension(SPEC_MAX, SPEC_MAX) :: array3	! Genreic, temporary
  integer, dimension(SPEC_MAX) :: narray			! Generic, temporary
  !	integer, dimension(SPEC_MAX, SPEC_MAX) :: narray2	! Generic, temporary
  !	character(NAME_MAX), allocatable :: cpnames(:)		! Component names
  !	real*8, allocatable :: cpguess(:)			! Guesses for initial component concentrations
  !	integer ncp						! Number of components
  character(NAME_MAX), allocatable :: wtnames(:)		! Information on water types will be held internally until the water types are applied to zones.
  character(NAME_MAX), allocatable :: wtspecies(:)		! List of species according to water
  real*8, allocatable :: wtgrid(:, :)				! wtspecies go across the top, and wtnames go down the left side.  wtgrid(i, j) is water type i and species j.
  integer nwt, nwtspecies					! Number of water types and water type species, respectively
  character(NAME_MAX), allocatable :: rtnames(:), rtspecies(:)	! Solid species names and component species.
  real*8, allocatable :: rtgrid(:, :)				! rtspecies go across the top, and rtnames go down the left side.  rtgrid(i, j) is rock type i and species j.
  integer nrt, nrtspecies					! Number of rock types and rock species, respectively.
  character(NAME_MAX), allocatable :: vtnames(:), vtspecies(:)	! Vapor species names and component species.
  real*8, allocatable :: vtgrid(:, :)				! vtgrid(i, j) is vapor type i and species j.
  integer nvt, nvtspecies					! Number of vapor types and vapor species.
  !	real*8, allocatable :: stoichgrid(:, :)			! Stores the stoichiometry data if a matrix is given.  stoichgrid(i, j) is complex i and compnent j, where names of complexes and components are stored in wtnames(nwt) + rtnames(nrt) and cpntnam(ncpnt), respectively.
  !	character(NAME_MAX), allocatable :: scnames(:)		! The component names as given in the stoichiometry matrix section.
  !	character(NAME_MAX), allocatable :: ssnames(:)		! The species names as given in the stoichiometry matrix section.
  character(NAME_MAX), allocatable :: species(:)		! The list of species according to spec
  character(NAME_MAX), allocatable :: masters(:)		! Master species for each component in species(:)
  real*8, allocatable :: guesses(:)				! Concentration guesses for each component in species(:)
  character, allocatable :: states(:)				! The state of each species
  !	character(NAME_MAX), allocatable :: ssids(:)		! The id of each non-aqueous reaction
  !	integer snspecies, sncpnt, naq, nso, nh, nv		! Number of species and components given in the stoichiometry matrix, and number of species of each state
  integer nspecies, nmasters, naq, nso, nh, nv			! Number of all, aqueous, solid, henry, and vapor species
  integer anspecies						! Number of species according to sorp
  character(NAME_MAX), allocatable :: sorpnames(:), sorpspecs(:)	! Names of adsorption models and species
  integer, allocatable :: lsorptypes(:, :), vsorptypes(:, :)	! Types of adsorption models
  logical fltype, fvtype, fa1l, fa2l, fbl, fa1v, fa2v, fbv	! Flags for whether each column has been specified
  real*8, allocatable :: a1l(:, :), a2l(:, :), bl(:, :), a1v(:, :), a2v(:, :), bv(:, :)	! Parameters used for sorp
  integer nsorp							! Number of adsorption models
  integer ndisp							! Number of dispersivity models
  character(NAME_MAX), allocatable :: dispnames(:)		! Names of dispersivity types
  character(NAME_MAX) :: dispparams(SPEC_MAX)			! Dispersivity parameters
  integer dispmode						! 0 for xyz, 1 for ldsp
  real*8, allocatable :: lx(:), ly(:), lz(:), ll(:), lt(:), vx(:), vy(:), vz(:), vl(:), vt(:)	! X, Y, Z, longitudinal, and transverse dispersivity for liquids and vapors
  integer ndiff							! Number of diffusion species
  character(NAME_MAX), allocatable :: diffspecs(:)		! Species as given in diff
  real*8, allocatable :: ldiff(:), vdiff(:)			! Diffusion coefficients for liquid and vapor
  integer nht, nhtspecies					! Number of Henry's Law types and number of species given in henry
  character(NAME_MAX), allocatable :: htnames(:), htspecies(:)	! Names of Henry models
  real*8, allocatable :: htgrid(:, :)				! Grid of Henry's Law types and their constituent species
  integer, allocatable :: hmodels(:)				! Types of Henry models
  character(NAME_MAX), allocatable :: hhtspecies(:)		! Henry's Law species in hparam
  integer hnhtspecies						! Number of Henry's Law species according to hparam
  real, allocatable :: ah(:), ah1(:), ah2(:), ah3(:), ah4(:), ah5(:), dhh(:), hh(:)	! Parameters for Henry models
  integer gngroups, gnspecies					! Number of groups and species according to group
  character(NAME_MAX), allocatable :: gspecies(:)		! List of species according to group
  integer, allocatable :: ggroups(:)				! List of groups
  integer nprint						! Number of printable species
  character(NAME_MAX), allocatable :: printspecies(:)		! List of species to print
  integer nstoichcomponents, nstoichspecies			! Number of components and species in stoich
  integer ncomplexes						! Number of complexes
  character(NAME_MAX), allocatable :: complexes(:)		! List of complex names
  character(NAME_MAX), allocatable :: complexcontents(:, :)	! List of master species that make up each complex
  integer, allocatable :: complexstoich(:, :)			! Stoichiometry of the master species in each complex
  integer ndistmodels						! Number of distribution models
  character(NAME_MAX), allocatable :: distmodelnames(:)		! Names of distribution models
  real*8, allocatable :: distmodels(:, :, :)			! distmodels(i, j) is model i, datapoint j, temperature at k = 1 and distribution coefficient at k = 2
  integer nsolmodels						! Number of solubility models
  character(NAME_MAX), allocatable :: solmodelnames(:)		! Names of solubility models
  real*8, allocatable :: solmodels(:, :, :)			! solmodels(i, j) is model i, datapoint j, temperature at k = 1 and solubility constant at k = 2
  integer encomplexes						! Number of complexes according to equi
  character(NAME_MAX), allocatable :: ecomplexes(:)		! Complexes according to equi
  real*8, allocatable :: equconstants(:), enthalpies(:)		! Equilibrium constants and heats of formation for aqueous complexes
  logical logten						! Whether equi values are log 10 values
  integer nrxns							! Number of reactions
  integer, allocatable :: rxntypes(:)				! Type of each reaction
  !	character(NAME_MAX), allocatable :: rxnnames(:)		! Name of each reaction (defunct)
  character(NAME_MAX), allocatable :: reactants(:, :)		! Reactants for all reactions ("extra" chemicals for type 4)
  character(NAME_MAX), allocatable :: products(:, :)		! Products for all reactions
  integer, allocatable :: stoichiometries(:, :, :)		! Stoichiometry for all reaction types.  Some are tricky...
  real*8, allocatable :: rates(:)				! Rate coefficient for reaction types 1, 2, 7, 8; forward rate coefficient for type 3; qm for type 4; halflife for reaction type 5
  real*8, allocatable :: rates2(:)				! Reverse rate coefficient for reaction type 3; decay coefficient for type 4
  character(NAME_MAX), allocatable :: dscoefs(:)		! Distribution coefficient for reaction types 1, 2; solubility coefficient for types 7, 8
  real*8, allocatable :: xcoef(:)				! Extra coefficients:  maxconc for reaction type 2; surface area for types 7, 8
  character(NAME_MAX), allocatable :: bioparams(:, :)		! Parameters for reaction type 4.  bioparams(i, j) is reaction i, substrate name at j = 1, electron acceptor name at j = 2, biomass name at j = 3, ks at j = 4, ka at j = 5, phthresh at j = 6, yield at j = 7, xminit at j = 8
  character(NAME_MAX), allocatable :: icbioholder(:, :)		! Holds component names for ICBIO for reaction type 4
  real*8, allocatable :: porchange(:, :)			! Hold porosity change information for reaction type 8.  porchange(i, j) is reaction i, weight at j = 1 and density at j = 2
  integer nzones						! Number of zones
  character(NAME_MAX), allocatable :: zones(:)			! List of zones specified in the zone section
  character(NAME_MAX), allocatable :: zonespecs(:)		! List of keywords specified in the zone section (necessary for ordering the data in the grid)
  character(NAME_MAX), allocatable :: zonegrid(:, :)		! The assignments for each zone.  zonegrid(i, j) is zone i and keyword j.
  integer nzonespecs, znzones					! The number specifications for each zone, and the number of zones according to the zone macro.
  character(NAME_MAX), allocatable :: rxns(:, :)		! The set of reactions for each zone
  integer, allocatable :: znrxns(:)				! The number of reactions for each zone
  character(NAME_MAX), allocatable :: spnam(:)			! We may be able to remove this later for 
  integer zinit, zboun, zrock, zrxn, zdisp, zldiff, zvdiff, zsorp, ztpor, ztime, zsehdiff, zhparam, zhenry, zvap	! Column numbers for each specification in assign
  logical inflag, bflag						! Flags for whether we're using the default ("*") init and bound conditions.  False means we are.
  integer rxnnsolid, rxnnaqueous, rxnnvapor
  character(NAME_MAX), allocatable :: rxnsolid(:), rxnaqueous(:), rxnvapor(:)
  ! Variables provided for the copy-and-pasted portions of rdcon and read_rxn
  real*8 h_const, dvap_conc, tempx,tempy,tempz,templength, cord1x,cord1y,cord2x,cord2y, dispxavw,dispyavw,dispzavw,cord1z, &
       cord2z,x, rdum1, rdum2, mvol, psdum, roldum, sdum, sctmp, frac, ctmp, ctol, andum, antmp, anvtmp, mvolv, anqo
  real*8 xmsg(16), lkeq(20),eqtemp(20)
  integer mi,mim,ndummy,iret,iz,jj,ij,ii,ispecies,npt_subst, n_node_sets,i_node_set, ja, jb, jc, kz,iz2,iz2p,i1,i2, &
       ipivt,sehindexl,sehindexv,neqp1, io_stat, open_file, mfile, logkeq, ic,im,iv,ix,idum, nwds, igrp, &
       junk, a1,a2,a3, num_temps, itypex, zonen, iimm, ivap, ncp
  logical null1, readflag
  character*5 user_macro
  character*3 nstring
  character*80 input_msg
  integer msg(16), imsg(16)
  character*32 cmsg(16)
  integer, allocatable :: hflag(:)
  character*4 por_change
  integer icpnt

!!!!!
  !
  ! readtracrxn 1.8 by Matthew Schauer <mschauer@lanl.gov>
  ! Property of the Los Alamos National Laboratory.
  !
!!!!!

  linen = 0
  stoichmethod = 0
  fcomp = .false.
  fwater = .false.
  frock = .false.
  fdisp = .false.
  fassign = .false.
  fvap = .false.
  fheader = .false.
  fsorp = .false.
  fassign = .false.
  fhenry = .false.
  fhparam = .false.
  fdiff = .false.
  fspec = .false.
  fgroup = .false.
  fprint = .false.
  fdist = .false.
  fsol = .false.
  fequi = .false.
  stoichmethod = 0
  nrxns = 0
  nzonespecs = 0	! Zone preprocessing
  znzones = 0
  hnhtspecies = 0
  naq = 0
  nso = 0
  nh = 0
  nv = 0
  array = '*'
  an = 0
  a = 0
  znzones = 0
  do i = 1, n0
     if(izonef(i) .gt. znzones) then
        znzones = izonef(i)
     endif
  enddo
  if(znzones .eq. 0) then
     write(ierr, *) 'A zone macro must precede trxn.'
     write(ierr, *) 'Error in trxn in preprocessing.'
     stop
  endif
  allocate(zonespecs(SPEC_MAX))
  allocate(zones(znzones))
  allocate(zonegrid(znzones, SPEC_MAX))
  iskip = 0	! Header preprocessing
  rsdmax = 1e-9
  strac_max = 0.99
  logkeq = 0
  logten = .false.
  an0 = 0
  awc = 1
  upwgta = 1
  daycs = 0
  daycf = tims
  dayhs = 0
  dayhf = tims
  iaccmx = 10
  daycm = 1.2
  daycmm = 1
  daycmx = 1e20
  nprttrc = 1
  ayc = 1 - awc
  nrt = 0		! Presetting everything to zero.
  nrtspecies = 0
  nwt = 0
  nwtspecies = 0
  ncp = 0
  nvt = 0
  nvtspecies = 0
  nht = 0
  nhtspecies = 0
  nspecies = 0
  naq = 0
  nso = 0
  nh = 0
  nv = 0
  anspecies = 0
  nsorp = 0
  ndisp = 0
  ndiff = 0
  nht = 0
  nhtspecies = 0
  nzones = 0
  ncomplexes = 0
  dispmode = 0
  allocate(dispnames(1))
  allocate(lx(1))
  allocate(ly(1))
  allocate(lz(1))
  allocate(ll(1))
  allocate(lt(1))
  allocate(vx(1))
  allocate(vy(1))
  allocate(vz(1))
  allocate(vl(1))
  allocate(vt(1))
  gngroups = 0
  ! Read main grids
  if(debug) write(iptty, *) 'Reading data for trxn.'
  reading = .true.
  do
100  read(inpt, '(a200)', end=4000) line	! Read the line.
     linen = linen + 1
     ! First, strip off comments.
     if(index(line, '#') .ne. 0) then
        line = line(1:index(line, '#') - 1)
     endif
     ! Then determine the contents of the line.  It could be...
     ! ...a blank line, in which case it is ignored.
     if((len_trim(line) .eq. 0) .or. (line(1:1) .eq. '	')) then
        goto 2000
     endif
     ! ...a line beginning with a keyword, in which case different modes are entered depending on the keyword:
     read(line, *, err=667) keyword
     ! ..."header," in which case we are reading boilerplate information.
     if(keyword .eq. 'ctrl') then
        if(debug) write(iptty, *) 'Reading control information...'
     elseif(keyword .eq. 'header') then
        if(debug) write(iptty, '(a31)') 'Reading header information...  '
        if(fheader .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  header has already been specified.  ', &
                'The new values will overwrite the old ones.'
        endif
        fheader = .true.
        linen = linen + 1
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        read(line, *, end=118, err=118) an0, awc, epc, upwgta
        linen = linen + 1
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        read(line, *, end=118, err=118) daycs, daycf, dayhf, dayhs
        linen = linen + 1
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        read(line, *, end=118, err=118) iaccmx, daycm, daycmm, daycmx, nprttrc
        do
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           array = '*'
           read(line, *, end=138, err=118) (array(i), i = 1, SPEC_MAX)
138        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .eq. 0) then
                 write(ierr, *) "Unknown token in header:  ", array(i)
                 goto 666
              endif
              if(array(i)(1:index(array(i), '=') - 1) .eq. 'iskip') then
                 read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *) iskip
              elseif(array(i)(1:index(array(i), '=') - 1) .eq. 'rsdmax') then
                 read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *) rsdmax
              elseif(array(i)(1:index(array(i), '=') - 1) .eq. 'strac_max') then
                 read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *) strac_max
              else
                 write(ierr, *) 'Error:  Variable "', array(i)(1:index(array(i), '=') - 1), &
                      '" in header not recognized.'
                 goto 666
              endif
           enddo
        enddo
        daycmx = min(daycmx, daymax)
        dtotc=daycmm * 86400.
        ayc = 1 - awc
        goto 2000
118     write(ierr, *) 'Error reading header information.'
        goto 666
        ! Component definitions
     elseif(keyword .eq. 'comp') then
        if(debug) write(iptty, '(a27)', advance='no') 'Reading component data...  '
        if(fcomp .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  comp has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(species)
           deallocate(states)
           deallocate(masters)
           deallocate(guesses)
        endif
        fcomp = .true.
        nspecies = 0
        naq = 0
        nso = 0
        nv = 0
        nh = 0
        backspace inpt
        read(inpt, '(a200)') line
        array = '*'
        read(line, *, end=198) keyword, (array(i), i = 1, SPEC_MAX)
198     do
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           nspecies = nspecies + 1
        enddo
        do i = 1, nspecies + 1
           backspace inpt
        enddo
        allocate(species(nspecies))
        allocate(states(nspecies))
        allocate(masters(nspecies))
        allocate(guesses(nspecies))
        masters = '*'
        guesses = 1e-9
        do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              exit
           elseif((array(i)(1:6) .ne. 'master') .and. (array(i)(1:5) .ne. 'guess')) then
              write(ierr, *) 'Column header "', array(i), '" not recognized in comp.'
              goto 666
           endif
        enddo
        do i = 1, nspecies
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              write(ierr, *) 'Internal error reading components.'
              goto 667
           endif
           if(array(1)(1:6) .eq. 'master') then
              if(array(2)(1:5) .eq. 'guess') then
                 read(line, *, err=167, end=167) states(i), species(i), masters(i), e
                 if(e .eq. '*') then
                    guesses(i) = 1e-9
                 else
                    read(e, *, err=167) guesses(i)
                 endif
              else
                 read(line, *, err=167, end=167) states(i), species(i), masters(i)
                 guesses(i) = 1e-9
              endif
           elseif(array(1)(1:5) .eq. 'guess') then
              if(array(2)(1:6) .eq.'master') then
                 read(line, *, err=167, end=167) states(i), species(i), e, masters(i)
              else
                 read(line, *, err=167, end=167) states(i), species(i), e
              endif
              if(e .eq. '*') then
                 guesses(i) = 1e-9
              else
                 read(e, *, err=167) guesses(i)
              endif
           else
              read(line, *, err=167, end=167) states(i), species(i)
           endif
           if(states(i) .eq. 'a') then
              naq = naq + 1
           elseif(states(i) .eq. 's') then
              nso = nso + 1
           elseif(states(i) .eq. 'h') then
              nh = nh + 1
           elseif(states(i) .eq. 'v') then
              nv = nv + 1
           else
              write(ierr, *) 'Error:  component state "', states(i), '" not recognized.'
              goto 666
           endif
        enddo
        nmasters = 0
        do i = 1, nspecies
           if(masters(i) .ne. '*') then
              nmasters = nmasters + 1
           endif
        enddo
        goto 168
167     write(ierr, *) 'Error reading components in comp.'
        goto 666
168     ncp = nspecies
        if(debug) write(iptty, '(a5, i'//aINT_LEN//', a12)') 'Read ', nspecies, ' components.'
        ! Water type mode
     elseif(keyword .eq. 'water') then
        if(debug) write(iptty, '(a26)', advance='no') 'Reading water type data...'
        if(fwater .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  water has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(wtspecies)
           deallocate(wtnames)
           deallocate(wtgrid)
        endif
        fwater = .true.
        array = '*'
        read(line, *, end=106, err=131) keyword, (array(i), i = 1, SPEC_MAX)
106     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nwtspecies = i - 1
              exit
           endif
        enddo
        allocate(wtspecies(nwtspecies))
        do i = 1, nwtspecies
           wtspecies(i) = array(i)
        enddo
        nwt = 0
        do
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
              goto 105		
           endif
           nwt = nwt + 1
        enddo
105     do i = 1, nwt + 1
           backspace inpt
        enddo
        linen = linen - 1
        allocate(wtnames(nwt))
        allocate(wtgrid(nwt, nwtspecies))
        do i = 1, nwt
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=131, err=131) wtnames(i), (wtgrid(i, j), j = 1, nwtspecies)
        enddo
        goto 132
131     write(ierr, *) 'Error reading values in water.'
        goto 666
132     if(debug) write(iptty, '(i'//aINT_LEN//', a16, i'//aINT_LEN//', a14)') nwt, &
             ' water type(s), ', nwtspecies, ' species read.'
        ! Rock mode -- Note that aside from the variable names, this is identical to the water reader.
     elseif(keyword .eq. 'rock') then
        if(debug) write(iptty, '(a29)', advance='no') 'Reading solid species data...'
        if(frock .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  rock has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(rtnames)
           deallocate(rtspecies)
           deallocate(rtgrid)
        endif
        frock = .true.
        array = '*'
        read(line, *, end=107, err=129) keyword, (array(i), i = 1, SPEC_MAX)
107     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nrtspecies = i - 1
              exit
           endif
        enddo
        allocate(rtspecies(nrtspecies))
        do i = 1, nrtspecies
           rtspecies(i) = array(i)
        enddo
        nrt = 0
        do
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
              goto 109		
           endif
           nrt = nrt + 1
        enddo
109     do i = 1, nrt + 1
           backspace inpt
        enddo
        linen = linen - 1
        allocate(rtnames(nrt))
        allocate(rtgrid(nrt, nrtspecies))
        do i = 1, nrt
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=129, err=129) rtnames(i), (rtgrid(i, j), j = 1, nrtspecies)
        enddo
        goto 130
129     write(ierr, *) 'Error reading values in rock.'
        goto 666
130     if(debug) write(iptty, '(i'//aINT_LEN//', a15, i'//aINT_LEN//', a14)') nrt, &
             ' rock type(s), ', nrtspecies, ' species read.'
        ! Vapor species mode -- this is also very similar to the water reader
     elseif(keyword .eq. 'vap') then
        if(debug) write(iptty, '(a29)', advance='no') 'Reading vapor species data...'
        if(fvap .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  vapor has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(vtnames)
           deallocate(vtspecies)
           deallocate(vtgrid)
        endif
        fvap = .true.
        array = '*'
        read(line, *, end=160, err=153) keyword, (array(i), i = 1, SPEC_MAX)
160     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nvtspecies = i - 1
              exit
           endif
        enddo
        allocate(vtspecies(nvtspecies))
        do i = 1, nvtspecies
           vtspecies(i) = array(i)
        enddo
        nvt = 0
        do
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
              goto 155		
           endif
           nvt = nvt + 1
        enddo
155     do i = 1, nvt + 1
           backspace inpt
        enddo
        linen = linen - 1
        allocate(vtnames(nvt))
        allocate(vtgrid(nvt, nvtspecies))
        do i = 1, nvt
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=153, err=153) vtnames(i), (vtgrid(i, j), j = 1, nvtspecies)
        enddo
        goto 154
153     write(ierr, *) 'Error reading values in vapor.'
        goto 666
154     if(debug) write(iptty, '(i'//aINT_LEN//', a16, i'//aINT_LEN//', a14)') nvt, &
             ' vapor type(s), ', nvtspecies, ' species read.'
        ! Henry's Law mode
     elseif(keyword .eq. 'henry') then
        if(debug) write(iptty, '(a35)', advance='no') 'Reading Henry''s Law species data...'
        if(fhenry .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  henry has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(htnames)
           deallocate(htspecies)
           deallocate(htgrid)
        endif
        fhenry = .true.
        array = '*'
        read(line, *, end=156, err=158) keyword, (array(i), i = 1, SPEC_MAX)
156     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nhtspecies = i - 1
              exit
           endif
        enddo
        allocate(htspecies(nhtspecies))
        do i = 1, nhtspecies
           htspecies(i) = array(i)
        enddo
        nht = 0
        do
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
              goto 157		
           endif
           nht = nht + 1
        enddo
157     do i = 1, nht + 1
           backspace inpt
        enddo
        linen = linen - 1
        allocate(htnames(nht))
        allocate(htgrid(nht, nhtspecies))
        do i = 1, nht
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=158, err=158) htnames(i), (htgrid(i, j), j = 1, nhtspecies)
        enddo
        goto 159
158     write(ierr, *) 'Error reading values in henry.'
        goto 666
159     if(debug) write(iptty, '(i'//aINT_LEN//', a22, i'//aINT_LEN//', a14)') nvt, &
             ' Henry''s Law type(s), ', nhtspecies, ' species read.'
        ! Henry's Law parameter mode
     elseif(keyword .eq. 'hparam') then
        if(debug) write(iptty, '(a35)', advance='no') 'Reading Henry''s Law parameters...  '
        if(fhparam .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  hparam has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(hhtspecies)
           deallocate(hmodels)
           deallocate(ah)
           deallocate(ah1)
           deallocate(ah2)
           deallocate(ah3)
           deallocate(ah4)
           deallocate(ah5)
           deallocate(dhh)
           deallocate(hh)
        endif
        fhparam = .true.
        hnhtspecies = 0
        do
           read(inpt, '(a200)', end=145, err=145) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           hnhtspecies = hnhtspecies + 1
        enddo
        allocate(hhtspecies(hnhtspecies))
        allocate(hmodels(hnhtspecies))
        allocate(ah(hnhtspecies))
        allocate(ah1(hnhtspecies))
        allocate(ah2(hnhtspecies))
        allocate(ah3(hnhtspecies))
        allocate(ah4(hnhtspecies))
        allocate(ah5(hnhtspecies))
        allocate(dhh(hnhtspecies))
        allocate(hh(hnhtspecies))
        do i = 0, hnhtspecies
           backspace inpt
        enddo
        do i = 1, hnhtspecies
           array = '*'
           linen = linen + 1
           read(inpt, '(a200)', end=145, err=145) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           read(line, *, err=145, end=147) (array(j), j = 1, SPEC_MAX)
147        hhtspecies(i) = array(1)
           if((array(2)(1:1) .eq. 'h') .or. (array(2) .eq. '1') .or. (array(2) .eq. '*')) then
              hmodels(i) = 1
           elseif((array(2)(1:1) .eq. 'm') .or. (array(2) .eq. '2')) then
              hmodels(i) = 2
           elseif((array(2)(1:1) .eq. 'w') .or. (array(2) .eq. '3')) then
              hmodels(i) = 3
           else
              write(ierr, *) 'Error:  Unrecognized Henry''s model "', array(2), '" in hparam.'
              goto 666
           endif
           do j = 3, SPEC_MAX
              if(array(j) .eq. '*') then
                 exit
              endif
              e = array(j)(1: index(array(j), '=') - 1)
              f = array(j)(index(array(j), '=') + 1:len_trim(array(j)))
              if(hmodels(i) .eq. 1) then
                 if(e .eq. 'ah') then
                    read(f, *, err=148, end=148) ah(i)
                 elseif(e .eq. 'dhh') then
                    read(f, *, err=148, end=148) dhh(i)
                 elseif((e(1:2) .eq. 'ah') .or. (e .eq. 'hh')) then

                 else
                    write(ierr, *) 'Error:  Unknown hparam variable "', e, '."'
                    goto 666
                 endif
              elseif(hmodels(i) .eq. 2) then
                 if(e .eq. 'ah1') then
                    read(f, *, err=148, end=148) ah1(i)
                 elseif(e .eq. 'ah2') then
                    read(f, *, err=148, end=148) ah2(i)
                 elseif(e .eq. 'ah3') then
                    read(f, *, err=148, end=148) ah3(i)
                 elseif(e .eq. 'ah4') then
                    read(f, *, err=148, end=148) ah4(i)
                 elseif(e .eq. 'ah5') then
                    read(f, *, err=148, end=148) ah5(i)
                 elseif((e .eq. 'ah') .or. (e .eq. 'dhh') .or. (e .eq. 'hh')) then

                 else
                    write(ierr, *) 'Error:  Unknown hparam variable "', e, '."'
                    goto 666
                 endif
              elseif(hmodels(i) .eq. 3) then
                 if(e .eq. 'ah') then
                    read(f, *, err=148, end=148) ah(i)
                 elseif(e .eq. 'hh') then
                    read(f, *, err=148, end=148) hh(i)
                 elseif((e(1:2) .eq. 'ah') .or. (e .eq. 'dhh')) then

                 else
                    write(ierr, *) 'Error:  Unknown hparam variable "', e, '."'
                    goto 666
                 endif
              else
                 write(ierr, *) 'Internal error processing hparam arguments.'
                 goto 667
              endif
           enddo
        enddo
        goto 146
145     write(ierr, *) 'Error reading Henry''s Law data.'
        goto 666
148     write(ierr, *) 'Error:  Invalid number "', f, '" in henry.'
        goto 666
146     if(debug) write(iptty, '(i'//aINT_LEN//', a13)') hnhtspecies, ' models read.'
        ! Aqueous reaction mode - database lookup
     elseif(keyword .eq. 'lookup') then
        write(ierr, *) 'Error:  Database lookup not yet supported.'
        goto 666
        ! Adsorption mode
     elseif(keyword .eq. 'sorp') then
        if(debug) write(iptty, '(a28)', advance='no') 'Reading adsorption data...  '
        if(fsorp .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  sorp has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(sorpnames)
           deallocate(sorpspecs)
           deallocate(lsorptypes)
           deallocate(vsorptypes)
           deallocate(a1l)
           deallocate(a2l)
           deallocate(bl)
           deallocate(a1v)
           deallocate(a2v)
           deallocate(bv)
        endif
        fsorp = .true.
        fltype = .false.
        fvtype = .false.
        fa1l = .false.
        fa2l = .false.
        fbl = .false.
        fa1v = .false.
        fa2v= .false.
        fbv = .false.
        backspace inpt
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        array = '*'
        read(line, *, end=122, err=124) keyword, (array(i), i = 1, SPEC_MAX)
122     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              a = i - 1
              exit
           endif
        enddo
        anspecies = 0
        nsorp = 0
        a = 0
        do
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              anspecies = a
              exit
           endif
           linen = linen + 1
           if(line(1:1) .eq. '.') then
              nsorp = nsorp + 1
              if((anspecies .ne. 0) .and. (anspecies .ne. a)) then
                 write(ierr, *) 'Error:  Different number of species specified in different', &
                      'adsorption models:  ', a, 'and', anspecies, '.'
                 goto 666
              endif
              anspecies = a
              a = 0
           else
              a = a + 1
           endif
        enddo
        if(nsorp .eq. 1) then
           anspecies = a
        endif
        allocate(sorpspecs(anspecies))
        allocate(lsorptypes(nsorp+1, anspecies))
        allocate(vsorptypes(nsorp+1, anspecies))
        allocate(a1l(nsorp+1, anspecies))
        allocate(a2l(nsorp+1, anspecies))
        allocate(bl(nsorp+1, anspecies))
        allocate(a1v(nsorp+1, anspecies))
        allocate(a2v(nsorp+1, anspecies))
        allocate(bv(nsorp+1, anspecies))
        allocate(sorpnames(nsorp+1))
        sorpnames = '*'
        do i = 1, nsorp * anspecies + nsorp + 1
           backspace inpt
           linen = linen - 1
        enddo
        linen = linen + 6
        do i = 1, nsorp
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           linen = linen + 1
           if(line(1:1) .ne. '.') then
              write(ierr, *) 'Error reading adsorption models.'
              goto 666
           endif
           sorpnames(i) = line(2:len_trim(line))
           if(index(line, ' ') .ne. 0) then
              sorpnames(i) = sorpnames(i)(1:index(sorpnames(i), ' '))
           endif
           do k = 1, anspecies
              array2 = '*'
              read(inpt, '(a200)') line
              if(index(line, '#') .ne. 0) then
                 line = line(1:index(line, '#') - 1)
              endif
              read(line, *, end=123, err=124) array3(i, k), (array2(j), j = 1, SPEC_MAX)
123           do j = 1, SPEC_MAX
                 if(array(j) .eq. 'ltype') then
                    fltype = .true.
                    if((array2(j) .eq. '0') .or. (array2(j)(1:2) .eq. 'co') .or. &
                         (array2(j) .eq. '*')) then
                       lsorptypes(i, k) = 0
                    elseif((array2(j) .eq. '1') .or. (array2(j)(1:2) .eq. 'li')) then
                       lsorptypes(i, k) = 1
                    elseif((array2(j) .eq. '2') .or. (array2(j)(1:2) .eq. 'fr')) then
                       lsorptypes(i, k) = 2
                    elseif((array2(j) .eq. '3') .or. (array2(j)(1:2) .eq. 'mf')) then
                       lsorptypes(i, k) = 3
                    elseif((array2(j) .eq. '4') .or. (array2(j)(1:2) .eq. 'la')) then
                       lsorptypes(i, k) = 4
                    else
                       write(ierr, *) 'Error:  Liquid adsorption type "', array2(j), &
                            '" in sorp not recognized.'
                       goto 666
                    endif
                 elseif(array(j) .eq. 'vtype') then
                    fvtype = .true.
                    if((array2(j) .eq. '0') .or. (array2(j)(1:2) .eq. 'co') .or. &
                         (array2(j) .eq. '*')) then
                       vsorptypes(i, k) = 0
                    elseif((array2(j) .eq. '1') .or. (array2(j)(1:2) .eq. 'li')) then
                       vsorptypes(i, k) = 1
                    elseif((array2(j) .eq. '2') .or. (array2(j)(1:2) .eq. 'fr')) then
                       vsorptypes(i, k) = 2
                    elseif((array2(j) .eq. '3') .or. (array2(j)(1:2) .eq. 'mf')) then
                       vsorptypes(i, k) = 3
                    elseif((array2(j) .eq. '4') .or. (array2(j)(1:2) .eq. 'la')) then
                       vsorptypes(i, k) = 4
                    else
                       write(ierr, *) 'Error:  Vapor adsorption type "', array2(j), &
                            '" in sorp not recognized.'
                       goto 666
                    endif
                 elseif(array(j) .eq. 'a1l') then
                    fa1l = .true.
                    if(array2(j) .eq. '*') then
                       a1l(i, k) = 0
                    else
                       read(array2(j), *, err=124) a1l(i, k)
                    endif
                 elseif(array(j) .eq. 'a2l') then
                    fa2l = .true.
                    if(array2(j) .eq. '*') then
                       a2l(i, k) = 0
                    else
                       read(array2(j), *, err=124) a2l(i, k)
                    endif
                 elseif(array(j) .eq. 'bl') then
                    fbl = .true.
                    if(array2(j) .eq. '*') then
                       bl(i, k) = 0
                    else
                       read(array2(j), *, err=124) bl(i, k)
                    endif
                 elseif(array(j) .eq. 'a1v') then
                    fa1v = .true.
                    if(array2(j) .eq. '*') then
                       a1v(i, k) = 0
                    else
                       read(array2(j), *, err=124) a1v(i, k)
                    endif
                 elseif(array(j) .eq. 'a2v') then
                    fa2v = .true.
                    if(array2(j) .eq. '*') then
                       a2v(i, k) = 0
                    else
                       read(array2(j), *, err=124) a2v(i, k)
                    endif
                 elseif(array(j) .eq. 'bv') then
                    fbv = .true.
                    if(array2(j) .eq. '*') then
                       bv(i, k) = 0
                    else
                       read(array2(j), *, err=124) bv(i, k)
                    endif
                 elseif(array(j) .eq. '*') then
                    exit
                 else
                    write(ierr, *) 'Error:  Column "', trim(array(j)), '" in sorp not recognized.'
                    goto 666
                 endif
              enddo
           enddo
        enddo
        linen = linen + 1
        if(fltype .eqv. .false.) then
           lsorptypes = 0
        endif
        if(fvtype .eqv. .false.) then
           vsorptypes = 0
        endif
        if(fa1l .eqv. .false.) then
           a1l = 0
        endif
        if(fa2l .eqv. .false.) then
           a2l = 0
        endif
        if(fbl .eqv. .false.) then
           bl = 0
        endif
        if(fa1v .eqv. .false.) then
           a1v = 0
        endif
        if(fa2v .eqv. .false.) then
           a2v = 0
        endif
        if(fbv .eqv. .false.) then
           bv = 0
        endif
        do i = 1, anspecies
           sorpspecs(i) = array3(1, i)
           do j = 1, nsorp
              if(sorpspecs(i) .ne. array3(j, i)) then
                 write(ierr, *) 'Error:  Species do not match across all models in sorp.'
                 goto 666
              endif
           enddo
        enddo
        if(debug) write(iptty, '(i'//aINT_LEN//', a13, i'//aINT_LEN//', a26)') anspecies, ' species for ', &
             nsorp, ' adsorption model(s) read.'
        goto 2000
124     write(ierr, *) 'Error reading a number in sorp.'
        goto 666
        ! Diffusivity mode
     elseif(keyword .eq. 'diff') then
        if(debug) write(iptty, '(a34)', advance='no') 'Reading diffusion information...  '
        if(fdiff .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  diff has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(diffspecs)
           deallocate(ldiff)
           deallocate(vdiff)
        endif
        fdiff = .true.
        backspace inpt
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        e = ''
        f = ''
        read(line, *, end=162) keyword, e, f
162     ndiff = 0
        do
           read(inpt, '(a200)', end=161) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           ndiff = ndiff + 1
        enddo
161     allocate(diffspecs(ndiff))
        allocate(ldiff(ndiff))
        allocate(vdiff(ndiff))
        do i = 1, ndiff + 1
           backspace inpt
        enddo
        do i = 1, ndiff
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, err=152, end=163) diffspecs(i), array(1), array(2)
163        if(array(1) .eq. '*') then
              array(1) = '0'
           endif
           if(array(2) .eq. '*') then
              array(2) = '0'
           endif
           if((e(1:1) .eq. 'l') .and. (f(1:1) .eq. 'v')) then
              read(array(1), *, err=152, end=152) ldiff(i)
              read(array(2), *, err=152, end=152) vdiff(i)
           elseif((e(1:1) .eq. 'v') .and. (f(1:1) .eq. 'l')) then
              read(array(1), *, err=152, end=152) vdiff(i)
              read(array(2), *, err=152, end=152) ldiff(i)
           elseif((e(1:1) .eq. 'l') .and. (f .eq. '')) then
              read(array(1), *, err=152, end=152) ldiff(i)
              vdiff(i) = 0
           elseif((e(1:1) .eq. 'v') .and. (f .eq. '')) then
              read(array(1), *, err=152, end=152) vdiff(i)
              ldiff(i) = 0
           elseif(e .eq. '') then
              read(inpt, *, err=152, end=152) diffspecs(i)
              vdiff(i) = 0
              ldiff(i) = 0
           else
              write(ierr, *) 'Error:  Problem with headers in diff.'
              goto 666
           endif
        enddo
        goto 151
152     write(ierr, *) 'Error reading numbers in diff.'
        goto 666
151     if(debug) write(iptty, '(a11, i'//aINT_LEN//', a14)') 'Values for ', ndiff, ' species read.'
        ! Dispersivity mode
     elseif(keyword .eq. 'londisp' .or. keyword .eq. 'xyzdisp') then
        if(keyword .eq. 'londisp') then
           if(debug) write(iptty, '(a54)', advance='no') 'Reading longitudinal/transverse dispersivity data...  '
           dispmode = 1
        else
           if(debug) write(iptty, '(a36)', advance='no') 'Reading X/Y/Z dispersivity data...  '
           dispmode = 0
        endif
        if(fdisp .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  xyzdisp or londisp has already been specified.  ', &
                'The new values will overwrite the old ones.'
        endif
        fdisp = .true.
        backspace inpt
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        dispparams = '*'
        read(line, *, end=164, err=165) keyword, (dispparams(i), i = 1, SPEC_MAX)
164     ndisp = 0
        do
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              goto 137
           endif
           ndisp = ndisp + 1
        enddo
137     do i = 0, ndisp
           backspace inpt
        enddo
        deallocate(dispnames)
        deallocate(lx)
        deallocate(ly)
        deallocate(lz)
        deallocate(ll)
        deallocate(lt)
        deallocate(vx)
        deallocate(vy)
        deallocate(vz)
        deallocate(vl)
        deallocate(vt)
        allocate(dispnames(ndisp+1))
        allocate(lx(ndisp+1))
        allocate(ly(ndisp+1))
        allocate(lz(ndisp+1))
        allocate(ll(ndisp+1))
        allocate(lt(ndisp+1))
        allocate(vx(ndisp+1))
        allocate(vy(ndisp+1))
        allocate(vz(ndisp+1))
        allocate(vl(ndisp+1))
        allocate(vt(ndisp+1))
        lx = 0
        ly = 0
        lz = 0
        vx = 0
        vy = 0
        vz = 0
        ll = 0
        lt = 0
        vl = 0
        vt = 0
        do i = 1, ndisp
           read(inpt, '(a200)') line
           linen = linen + 1
           array = '*'
           if(dispmode .eq. 0) then
              read(line, *, end=139, err=140) dispnames(i), (array(j), j = 1, SPEC_MAX)	!lx(i), ly(i), lz(i), vx(i), vy(i), vz(i)
139           do j = 1, SPEC_MAX
                 if((dispparams(j) .eq. '*')) then
                    exit
                 endif
                 if(array(j) .eq. '*') then
                    array(j) = '0'
                 endif
                 if(dispparams(j) .eq. 'lx') then
                    read(array(j), *) lx(i)
                 elseif(dispparams(j) .eq. 'ly') then
                    read(array(j), *) ly(i)
                 elseif(dispparams(j) .eq. 'lz') then
                    read(array(j), *) lz(i)
                 elseif(dispparams(j) .eq. 'vx') then
                    read(array(j), *) vx(i)
                 elseif(dispparams(j) .eq. 'vy') then
                    read(array(j), *) vy(i)
                 elseif(dispparams(j) .eq. 'vz') then
                    read(array(j), *) vz(i)
                 else
                    write(ierr, *) 'Unknown xyzdisp header ', trim(dispparams(j)), '.'
                    goto 666
                 endif
              enddo
           else
              read(line, *, end=166, err=140) dispnames(i), (array(j), j = 1, SPEC_MAX)	!ll(i), lt(i), vl(i), vt(i)
166           do j = 1, SPEC_MAX
                 if((dispparams(j) .eq. '*') .or. (array(j) .eq. '*')) then
                    exit
                 elseif(dispparams(j) .eq. 'll') then
                    read(array(j), *) ll(i)
                 elseif(dispparams(j) .eq. 'lt') then
                    read(array(j), *) lt(i)
                 elseif(dispparams(j) .eq. 'vl') then
                    read(array(j), *) vl(i)
                 elseif(dispparams(j) .eq. 'vt') then
                    read(array(j), *) vt(i)
                 else
                    write(ierr, *) 'Unknown londisp header ', trim(dispparams(j)), '.'
                    goto 666
                 endif
              enddo
           endif
        enddo
        goto 141
140     write(ierr, *) 'Error reading numbers in disp.'
        goto 666
165     write(ierr, *) 'Error reading column headers in disp.'
        goto 666
141     if(debug) write(iptty, '(a19, i'//aINT_LEN//', a10)') 'Read parameters for ', ndisp, ' model(s).'
        ! Component grouping mode
     elseif(keyword .eq. 'group') then
        if(debug) write(iptty, '(a36)', advance='no') 'Reading component grouping data...  '
        if(fgroup .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  group has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(gspecies)
           deallocate(ggroups)
        endif
        fgroup = .true.
        gngroups = 0
        gnspecies = 0
        do
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           array = '*'
           read(line, *, end=178, err=179) (array(i), i = 1, SPEC_MAX)
178        i = 0
           do
              i = i + 1
              if(i .gt. SPEC_MAX) then
                 exit
              endif
              if(array(i) .eq. '*') then
                 gnspecies = gnspecies + i - 1
                 exit
              endif
           enddo
           gngroups = gngroups + 1
        enddo
        do i = 1, gngroups + 1
           backspace inpt
        enddo
        allocate(gspecies(gnspecies))
        allocate(ggroups(gnspecies))
        k = 0
        do i = 1, gngroups
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           linen = linen + 1
           array = '*'
           read(line, *, end=177, err=179) (array(j), j = 1, SPEC_MAX)
177        do j = 1, SPEC_MAX
              if(array(j) .eq. '*') then
                 exit
              endif
              k = k + 1
              gspecies(k) = array(j)
              ggroups(k) = i
           enddo
        enddo
        if(debug) write(iptty, '(a5, i'//aINT_LEN//', a15, i'//aINT_LEN//', a8)') 'Read ', &
             gnspecies, ' components in ', gngroups, ' groups.'		
        goto 2000
179     write(ierr, *) 'Error reading grouping information.'
        goto 666
        ! Information printing selection mode
     elseif(keyword .eq. 'print') then
        if(debug) write(iptty, *) 'Reading component printing data...  '
        if(fprint .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  print has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(printspecies)
        endif
        fprint = .true.
        array = '*'
        backspace inpt
        read(inpt, *, end=201) keyword, (array(i), i = 1, SPEC_MAX)
201     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nprint = i - 1
              exit
           endif
        enddo
        allocate(printspecies(nprint))
        do i = 1, nprint
           printspecies(i) = array(i)
           if(printspecies(i) .eq. 'all') then 
              if(debug) write(iptty, *) 'All species will print.'
              deallocate(printspecies)
              allocate(printspecies(1))
              nprint = 1
              printspecies(1) = 'all'
              goto 2000
           elseif(printspecies(i) .eq. 'none') then
              if(debug) write(iptty, *) 'No species will print.'
              deallocate(printspecies)
              allocate(printspecies(1))
              nprint = 1
              printspecies(1) = 'none'
              goto 2000
           endif
        enddo
        !			do i = 1, SPEC_MAX
        !				read(inpt, '(a200)') line
        !				if(index(line, '#') .ne. 0) then
        !					line = line(1:index(line, '#') - 1)
        !				endif
        !				if(len_trim(line) .eq. 0) then
        !					exit
        !				endif
        !				nprint = nprint + 1
        !			enddo
        !			allocate(printspecies(nprint))
        !			do i = 1, nprint
        !				backspace(inpt)
        !			enddo
        !			do i = 1, nprint
        !				read(inpt, '(a200)') line
        !				if(index(line, '#') .ne. 0) then
        !					line = line(1:index(line, '#') - 1)
        !				endif
        !				linen = linen + 1
        !				read(line, *, err=180) printspecies(i), null
        !			enddo
        if(debug) write(iptty, '(a5, i'//aINT_LEN//', a19)') 'Read ', nprint, ' printable species.'
        goto 2000
180     write(ierr, *) 'Error reading printable species.'
        goto 666
        ! Stoichiometry input mode
     elseif(keyword .eq. 'stoich') then
        if(debug) write(iptty, *) 'Reading stoichiometry data...  '
        if(stoichmethod .ne. 0) then
           if(debug) write(iptty, *) 'Warning:  Stoichiometry has already been handled.  ', &
                'The new values will overwrite the old ones.'
           deallocate(complexes)
           deallocate(complexcontents)
           deallocate(complexstoich)
        endif
        stoichmethod = 2
        ncomplexes = 0
        do i = 1, SPEC_MAX
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              ncomplexes = i - 1
              exit
           endif
        enddo
        do i = 1, ncomplexes + 1
           backspace inpt
        enddo
        allocate(complexes(ncomplexes))
        allocate(complexcontents(ncomplexes, SPEC_MAX))
        allocate(complexstoich(ncomplexes, SPEC_MAX))
        complexcontents = '*'
        complexstoich = 0
        do i = 1, ncomplexes
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=239) complexes(i), null, (complexstoich(i, j), complexcontents(i, j), null, j = 1, SPEC_MAX)
239        continue
        enddo
        ! Distribution coefficient model mode
     elseif(keyword .eq. 'dist') then
        if(debug) write(iptty, '(a32)', advance='no') 'Reading distribution models...  '
        if(fdist .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  dist has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(distmodels)
           deallocate(distmodelnames)
        endif
        fdist = .true.
        ndistmodels = 0
        i = 0
        j = 0
        do k = 1, SPEC_MAX
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           if(line(1:1) .eq. '.') then
              ndistmodels = ndistmodels + 1
              if(j .gt. i) then
                 i = j
              endif
              j = 0
           else
              j = j + 1
           endif
        enddo
        k = k - 1
        allocate(distmodels(ndistmodels, i, 2))
        allocate(distmodelnames(ndistmodels))
        distmodels = -20000.4	! An absurd value
        do i = 1, k + 1
           backspace inpt
        enddo
        j = 0
        o = -1	! Throw an error if the user doesn't start with a model name
        do i = 1, k
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(line(1:1) .eq. '.') then
              j = j + 1
              distmodelnames(j) = line(2:len_trim(line))
              o = 0
           else
              o = o + 1
              read(line, *, err=181) distmodels(j, o, 1), distmodels(j, o, 2)
           endif
        enddo
        if(debug) write(iptty, '(a5, i'//aINT_LEN//', a8)') 'Read ', ndistmodels, ' models.'
        goto 2000
181     write(ierr, *) 'Error reading models in dist.'
        goto 666
        ! Solubility coefficient model mode
     elseif(keyword .eq. 'sol') then
        if(debug) write(iptty, '(a30)', advance='no') 'Reading solubility models...  '
        if(fsol .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  sol has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(solmodels)
           deallocate(solmodelnames)
        endif
        fsol = .true.
        nsolmodels = 0
        i = 0
        j = 0
        do k = 1, SPEC_MAX
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              if(j .gt. i) then
                 i = j
              endif
              exit
           endif
           if(line(1:1) .eq. '.') then
              nsolmodels = nsolmodels + 1
              if(j .gt. i) then
                 i = j
              endif
              j = 0
           else
              j = j + 1
           endif
        enddo
        k = k - 1
        allocate(solmodels(nsolmodels, i, 2))
        allocate(solmodelnames(nsolmodels))
        solmodels = -20000.4	! An absurd value
        do i = 1, k + 1
           backspace inpt
        enddo
        j = 0
        o = -1	! Throw an error if the user doesn't start with a model name
        do i = 1, k
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(line(1:1) .eq. '.') then
              j = j + 1
              solmodelnames(j) = line(2:len_trim(line))
              o = 0
           else
              o = o + 1
              read(line, *, err=222) solmodels(j, o, 1), solmodels(j, o, 2)
           endif
        enddo
        if(debug) write(iptty, '(a5, i'//aINT_LEN//', a8)') 'Read ', nsolmodels, ' models.'
        goto 2000
222     write(ierr, *) 'Error reading models in sol.'
        goto 666
        ! Equilibrium constant input mode
     elseif(keyword .eq. 'equi') then
        ! TODO Make this less picky:  allow the user to specify headers "CKEQ" and "HEQ" to put them in different orders or leave the columns out to use default values.
        if(debug) write(iptty, *) 'Reading equilibrium constants...'
        if(fequi .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  equi has already been specified.  ', &
                'The new values will overwrite the old ones.'
           deallocate(ecomplexes)
           deallocate(equconstants)
           deallocate(enthalpies)
        endif
        fequi = .true.
        encomplexes = 0
        logten = .false.
        backspace inpt
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        e = '*'
        read(line, *, end=248) keyword, e
248     if(e(1:3) .eq. 'log') then
           logten = .true.
        endif
        do i = 1, SPEC_MAX
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              encomplexes = i
              exit
           endif
        enddo
        do i = 1, encomplexes + 1
           backspace inpt
        enddo
        allocate(ecomplexes(encomplexes))
        allocate(equconstants(encomplexes))
        allocate(enthalpies(encomplexes))
        equconstants = 0
        enthalpies = 0
        do i = 1, encomplexes
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           linen = linen + 1
           read(line, *, end=246, err=247) ecomplexes(i), equconstants(i), enthalpies(i)
246        if(logten .eqv. .true.) then
              equconstants(i) = 10**(equconstants(i))
              enthalpies(i) = 10**(enthalpies(i))
           endif
           continue
        enddo
        goto 2000
247     write(ierr, *) 'Error reading values in equi.'
        goto 666
        ! Reaction input mode
     elseif(keyword .eq. 'rxn') then
        ! TODO If the user doesn't include the key/value pair line, use defaults or throw an error
        ! TODO If the user leaves one or more values on a line out, use defaults or throw an error
        nrxns = nrxns + 1
        if(debug) write(iptty, '(a17, i'//aINT_LEN//', a5)', advance='no') 'Reading reaction ', nrxns, '...  '
        if(nrxns .eq. 1) then
           allocate(rxntypes(numrxn))
           !allocate(rxnnames(numrxn))
           allocate(reactants(numrxn, SPEC_MAX))
           allocate(products(numrxn, SPEC_MAX))
           allocate(stoichiometries(numrxn, 2, SPEC_MAX * 2))
           allocate(rates(numrxn))
           allocate(dscoefs(numrxn))
           allocate(rates2(numrxn))
           allocate(xcoef(numrxn))
           allocate(bioparams(numrxn, 8))
           allocate(icbioholder(numrxn, SPEC_MAX))
           allocate(porchange(numrxn, 2))
           reactants = '*'
           products = '*'
           dscoefs = '*'
           stoichiometries = 0
           icbioholder = '*'
           porchange = 0
        endif
        backspace inpt
        read(inpt, '(a200)') line
        if(index(line, '#') .ne. 0) then
           line = line(1:index(line, '#') - 1)
        endif
        read(line, *, err=182) keyword, rxntypes(nrxns)!, rxnnames(nrxns)
        if(rxntypes(nrxns) .eq. 1) then
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=186, err=182) (array(i), i = 1, SPEC_MAX)
186        reactants(nrxns, 1) = array(1)
           if(index(array(2), '=') .eq. 0) then
              goto 183
           endif
           products(nrxns, 1) = array(3)
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=233) (array(i), i = 1, SPEC_MAX)
233        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0) then
                 e = array(i)(1:index(array(i), '=') - 1)
                 f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:4) .eq. 'rate') then
                    read(f, *) rates(nrxns)
                 elseif(e(1:4) .eq. 'dist') then
                    read(f, *) dscoefs(nrxns)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read linear kinetic reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        elseif(rxntypes(nrxns) .eq. 2) then
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=187, err=182) (array(i), i = 1, SPEC_MAX)
187        reactants(nrxns, 1) = array(1)
           if(index(array(2), '=') .eq. 0) then
              goto 183
           endif
           products(nrxns, 1) = array(3)
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=234) (array(i), i = 1, SPEC_MAX)
234        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0) then
                 e = array(i)(1:index(array(i), '=') - 1)
                 f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:4) .eq. 'rate') then
                    read(f, *) rates(nrxns)
                 elseif(e(1:4) .eq. 'dist') then
                    read(f, *) dscoefs(nrxns)
                 elseif(e(1:3) .eq. 'max') then
                    read(f, *) xcoef(nrxns)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read Langmuir kinetic reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        elseif(rxntypes(nrxns) .eq. 3) then
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=188, err=182) (array(i), i = 1, SPEC_MAX)
188        flag = .false.
           j = 1
           do i = 1, SPEC_MAX
              if(index(array(i), '=') .ne. 0) then
                 if(flag .eqv. .true.) then
                    goto 183
                 endif
                 flag = .true.
                 j = 1
              elseif(array(i) .eq. '+') then
                 continue
              else
                 if(flag .eqv. .false.) then
                    read(array(i), *, err=189) stoichiometries(nrxns, 1, j)
                 else
                    read(array(i), *, err=190) stoichiometries(nrxns, 2, j)
                 endif
                 goto 191
189              read(array(i), *, err=183) reactants(nrxns, j)
                 stoichiometries(nrxns, 1, j) = 1
                 j = j + 1
                 goto 191
190              read(array(i), *, err=183) products(nrxns, j)
                 stoichiometries(nrxns, 2, j) = 1
                 j = j + 1
191              continue
              endif
           enddo
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=235) (array(i), i = 1, SPEC_MAX)
235        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0) then
                 e = array(i)(1:index(array(i), '=') - 1)
                 f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:3) .eq. 'for') then
                    read(f, *) rates(nrxns)
                 elseif(e(1:3) .eq. 'rev') then
                    read(f, *) rates2(nrxns)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read general kinetic reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        elseif(rxntypes(nrxns) .eq. 4) then
           do
              read(inpt, '(a200)') line
              if(index(line, '#') .ne. 0) then
                 line = line(1:index(line, '#') - 1)
              endif
              if(index(line, ':') .eq. 0) then
                 exit
              endif
              e = line(1:index(line, ':') - 1)
              f = line(index(line, ':') + 1:len_trim(line))
              if(e(1:4) .eq. 'subs') then
                 bioparams(nrxns, 1) = e
              elseif(e(1:4) .eq. 'elec') then
                 bioparams(nrxns, 2) = e
              elseif(e(1:7) .eq. 'biomass') then
                 bioparams(nrxns, 3) = e
              elseif(e(1:4) .eq. 'reac') then
                 array = '*'
                 read(e, *, end=214) (array(j), j = 1, SPEC_MAX)
                 k = 1
214              do j = 1, SPEC_MAX
                    if(array(k) .eq. '*') then
                       exit
                    endif
                    read(array(k), *, err=215) stoichiometries(i, 1, j)
                    k = k + 1
                    if(array(k) .eq. '*') then
                       goto 183
                    endif
                    goto 220
215                 stoichiometries(i, 1, j) = 1
220                 read(array(k), *, err=183, end=183) reactants(i, j)
                    k = k + 1
                 enddo
              elseif(e(1:4) .eq. 'prod') then
                 array = '*'
                 read(e, *, end=216) (array(j), j = 1, SPEC_MAX)
                 k = 1
216              do j = 1, SPEC_MAX
                    if(array(k) .eq. '*') then
                       exit
                    endif
                    read(array(k), *, err=217) stoichiometries(i, 2, j)
                    k = k + 1
                    if(array(k) .eq. '*') then
                       goto 183
                    endif
                    goto 221
217                 stoichiometries(i, 2, j) = 1
221                 read(array(k), *, err=183, end=183) products(i, j)
                    k = k + 1
                 enddo
              elseif(e(1:6) .eq. 'biodeg') then
                 read(f, *, end=192) (icbioholder(nrxns, i), i = 1, SPEC_MAX)
              else
                 goto 185
              endif
192           continue
           enddo
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=236) (array(i), i = 1, SPEC_MAX)
236        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0) then
                 e = array(i)(1:index(array(i), '=') - 1)
                 f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:2) .eq. 'ks') then
                    read(f, *) bioparams(nrxns, 4)
                 elseif(e(1:2) .eq. 'ka') then
                    read(f, *) bioparams(nrxns, 5)
                 elseif(e(1:3) .eq. 'dec') then
                    read(f, *) rates2(nrxns)
                 elseif(e(1:2) .eq. 'ph') then
                    read(f, *) bioparams(nrxns, 6)
                 elseif(e(1:2) .eq. 'qm') then
                    read(f, *) rates(nrxns)
                 elseif(e(1:5) .eq. 'yield') then
                    read(f, *) bioparams(nrxns, 7)
                 elseif(e(1:2) .eq. 'xm') then
                    read(f, *) bioparams(nrxns, 8)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read biodegradation reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        elseif(rxntypes(nrxns) .eq. 5) then
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=193) (array(i), i = 1, SPEC_MAX)
193        reactants(nrxns, 1) = array(1)
           if(index(array(2), '=') .eq. 0) then
              goto 183
           endif
           products(nrxns, 1) = array(3)
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=237) (array(i), i = 1, SPEC_MAX)
237        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0)then
                 e = array(1)(1:index(array(i), '=') - 1)
                 f = array(1)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:4) .eq. 'half') then
                    read(f, *) rates(nrxns)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read radioactive decay reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        elseif(rxntypes(nrxns) .eq. 6) then
           write(ierr, *) 'Error: reaction type 6 no longer supported.'
           goto 666
        elseif((rxntypes(nrxns) .eq. 7) .or.(rxntypes(nrxns) .eq. 8)) then
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=194, err=182) (array(i), i = 1, SPEC_MAX)
194        flag = .false.
           j = 1
           do i = 1, SPEC_MAX
              if(index(array(i), '=') .ne. 0) then
                 if(flag .eqv. .true.) then
                    goto 183
                 endif
                 flag = .true.
                 j = 1
              elseif(array(i) .eq. '+') then
                 continue
              else
                 if(flag .eqv. .false.) then
                    read(array(i), *, err=195) stoichiometries(nrxns, 1, j)
                 else
                    read(array(i), *, err=196) stoichiometries(nrxns, 2, j)
                 endif
                 goto 197
195              read(array(i), *, err=183) reactants(nrxns, j)
                 j = j + 1
                 goto 197
196              read(array(i), *, err=183) products(nrxns, j)
                 j = j + 1
197              continue
              endif
           enddo
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           array = '*'
           read(line, *, end=238) (array(i), i = 1, SPEC_MAX)
238        do i = 1, SPEC_MAX
              if(array(i) .eq. '*') then
                 exit
              endif
              if(index(array(i), '=') .ne. 0) then
                 e = array(i)(1:index(array(i), '=') - 1)
                 f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
                 if(e(1:3) .eq. 'sol') then
                    read(f, *) dscoefs(nrxns)
                 elseif(e(1:4) .eq. 'rate') then
                    read(f, *) xcoef(nrxns)
                 elseif(e(1:7) .eq. 'surface') then
                    read(f, *) dscoefs(nrxns)
                 elseif(e(1:3) .eq. 'mol') then
                    if(rxntypes(i) .eq. 7) then
                       write(ierr, *) 'Error:  Porosity change only available for ', &
                            'reaction type 8.'
                       goto 666
                    endif
                    read(f, *) porchange(nrxns, 1)
                 elseif(e(1:4) .eq. 'dens') then
                    if(rxntypes(i) .eq. 7) then
                       write(ierr, *) 'Error:  Porosity change only available for ', &
                            'reaction type 8.'
                       goto 666
                    endif
                    read(f, *) porchange(nrxns, 2)
                 else
                    goto 185
                 endif
              else
                 goto 184
              endif
           enddo
           if(debug) write(iptty, *) 'Read precipitation/dissolution reaction.' ! ', nrxns, '.' !"', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
           goto 2000
        else
           write(ierr, *) 'Error:  reaction type ', rxntypes(nrxns), ' not known for reaction ', nrxns, '.'
        endif
182     write(ierr, *) 'Error reading reaction.'
        goto 666
183     write(ierr, *) 'Error:  Check reaction format.'
        goto 666
184     write(ierr, *) 'Error:  Check key/value list syntax.'
        goto 666
185     write(ierr, *) 'Error:  Unknown key "', e(1:len_trim(e)), '" in key/value list.'
        goto 666
        ! Assignment mode
     elseif(keyword .eq. 'assign') then
        if(debug) write(iptty, '(a43)', advance='no') 'Reading zone attribute assignment data...  '
        if(fassign .eqv. .true.) then
           if(debug) write(iptty, *) 'Warning:  assign has already been specified.  ', &
                'The new values will overwrite the old ones.'
        endif
        fassign = .true.
        array = '*'
        zonespecs = '*'
        zones = '*'
        zonegrid = '*'
        read(line, *, end=117, err=144) keyword, (array(i), i = 1, SPEC_MAX)
117     do i = 1, SPEC_MAX
           if(array(i) .eq. '*') then
              nzonespecs = i - 1
              exit
           endif
        enddo
        do i = 1, nzonespecs
           zonespecs(i) = array(i)
        enddo
        nzones = 0
        do
           read(inpt, '(a200)') line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           nzones = nzones + 1
        enddo
        a = 0
        if(nzones .gt. znzones) then
           if(debug) write(iptty, *) 'Warning:  More zones were specified in the assign grid'// & 
                ' than in the zone macro.  Only the first ', znzones, ' zone(s) will be used.'
           a = znzones
        elseif(nzones .lt. znzones) then
           if(debug) write(iptty, *) 'Warning:  Fewer zones were specified in the assign grid'// &
                ' than in the zone macro.  Extra zones will receive the default settings.'
           a = nzones
        else
           a = znzones
        endif
        do i = 1, nzones + 1
           backspace inpt
        enddo
        linen = linen - 1
        do i = 1, a
           read(inpt, '(a200)') line
           linen = linen + 1
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           read(line, *, end=144, err=144) zones(i), (zonegrid(i, j), j = 1, nzonespecs)
        enddo
        do i = a + 1, nzones
           read(inpt, '(a200)') line
           linen = linen + 1
        enddo
        nzones = a
        goto 143
144     write(ierr, *) 'Error reading zone data.'
        goto 666
143     if(debug) write(iptty, '(a5, i'//aINT_LEN//', a18, i'//aINT_LEN//', a8)', advance='no') 'Read ', nzonespecs, &
             ' parameter(s) for ', a, ' zone(s)'
        if(nzonespecs .eq. 0) then
           if(debug) write(iptty, '(a1)') '.'
        else
           if(debug) write(iptty, '(a3)', advance='no') ':  '
           if(debug) write(iptty, '('//aSPEC_MAX//'a8)') (zonespecs(i), i = 1, nzonespecs)
        endif
        ! ..."endtrxn," in which case we know that we are all done.
     elseif(keyword(1:3) .eq. 'end') then
        goto 4000
        ! ...a nonsense line, in which case an error message is printed and an error status returned to the calling program.
     else
        write(ierr, *) 'Error:  keyword "'//trim(keyword)//'" not recognized.'
        goto 666
     endif
     ! Wrap up the processing of each line
2000 continue
  enddo
4000 reading = .false.
  ! Process the data and put them in a format that FEHM will like.
  if(debug) write(iptty, '(i'//aINT_LEN//', a12)') linen, ' lines read.'
  if(debug) write(iptty, *) 'Checking data...'
  ! Sanity checks!
  ! ...Ensure that the bare minimums (minima?) are present
  if(fcomp .eqv. .false.) then
     write(ierr, *) 'Error:  comp was not specified.'
     goto 666
  endif
  if(ncp .lt. 1) then
     write(ierr, *) 'Error:  At least one component must be specified.'
     goto 666
  endif
  if(fheader .eqv. .false.) then
     write(ierr, *) 'Error:  header was not specified.'
     goto 666
  endif
  ! ...Check that all the columns in assign are good and fill in values for missing columns.
  zinit = 0
  zboun = 0
  zrock = 0
  zdisp = 0
  zsorp = 0
  ztpor = 0
  zldiff = 0
  zvdiff = 0
  ztime = 0
  zvap = 0
  zhenry = 0
  zhparam = 0
  do i = 1, nzonespecs
     if(zonespecs(i) .eq. 'init') then
        zinit = i
     elseif(zonespecs(i) .eq. 'boun') then
        zboun = i
     elseif(zonespecs(i) .eq. 'rock') then
        zrock = i
     elseif(zonespecs(i) .eq. 'disp') then
        zdisp = i
     elseif(zonespecs(i) .eq. 'ldiff') then
        zldiff = i
     elseif(zonespecs(i) .eq. 'vdiff') then
        zvdiff = i
     elseif(zonespecs(i) .eq. 'sorp') then
        zsorp = i
     elseif(zonespecs(i) .eq. 'tpor') then
        ztpor = i
     elseif(zonespecs(i) .eq. 'time') then
        ztime = i
     elseif(zonespecs(i) .eq. 'vap') then
        zvap = i
     elseif(zonespecs(i) .eq. 'henry') then
        zhenry = i
     else
        write(ierr, *) 'Error:  Column "', zonespecs(i), '" in assign not recognized.'
        goto 666
     endif
  enddo
  ! ...Ensure that the appropriate section has been included for each column specified in assign.
  if(((zinit .ne. 0) .or. (zboun .ne. 0)) .and. (fwater .eqv. .false.)) then
     write(ierr, *) 'Error:  init or boun column specified in assign, but no water section specified.'
     goto 666
  elseif((zrock .ne. 0) .and. (frock .eqv. .false.)) then
     write(ierr, *) 'Error:  rock column specified in assign, but no rock section specified.'
     goto 666
  elseif((zsorp .ne. 0) .and. (fsorp .eqv. .false.)) then
     write(ierr, *) 'Error:  sorp column specified in assign, but no sorp section specified.'
     goto 666
  elseif((zdisp .ne. 0) .and. (fdisp .eqv. .false.)) then
     write(ierr, *) 'Error:  disp column specified in assign, but no xyzdisp or londisp section specified.'
     goto 666
  endif
  if(zinit .eq. 0) then
     nzonespecs = nzonespecs + 1
     zinit = nzonespecs
     zonespecs(zinit) = 'init'
     do i = 1, nzones
        zonegrid(i, zinit) = '*'
     enddo
  endif
  if(zboun .eq. 0) then
     nzonespecs = nzonespecs + 1
     zboun = nzonespecs
     zonespecs(zboun) = 'boun'
     do i = 1, nzones
        zonegrid(i, zboun) = '*'
     enddo
  endif
  if(zrock .eq. 0) then
     nzonespecs = nzonespecs + 1
     zrock = nzonespecs
     zonespecs(zrock) = 'rock'
     do i = 1, nzones
        zonegrid(i, zrock) = '*'
     enddo
  endif
  if(zdisp .eq. 0) then
     nzonespecs = nzonespecs + 1
     zdisp = nzonespecs
     zonespecs(zdisp) = 'disp'
     do i = 1, nzones
        zonegrid(i, zdisp) = '*'
     enddo
  endif
  if(zldiff .eq. 0) then
     nzonespecs = nzonespecs + 1
     zldiff = nzonespecs
     zonespecs(zldiff) = 'ldiff'
     do i = 1, nzones
        zonegrid(i, zldiff) = '0'
     enddo
  endif
  if(zvdiff .eq. 0) then
     nzonespecs = nzonespecs + 1
     zvdiff = nzonespecs
     zonespecs(zvdiff) = 'vdiff'
     do i = 1, nzones
        zonegrid(i, zvdiff) = '0'
     enddo
  endif
  if(zsorp .eq. 0) then
     nzonespecs = nzonespecs + 1
     zsorp = nzonespecs
     zonespecs(zsorp) = 'sorp'
     do i = 1, nzones
        zonegrid(i, zsorp) = '*'
     enddo
  endif
  if(ztime .eq. 0) then
     nzonespecs = nzonespecs + 1
     ztime = nzonespecs
     zonespecs(ztime) = 'time'
     do i = 1, nzones
        zonegrid(i, ztime) = '*'
     enddo
  endif
  if(zvap .eq. 0) then
     nzonespecs = nzonespecs + 1
     zvap = nzonespecs
     zonespecs(zvap) = 'vap'
     do i = 1, nzones
        zonegrid(i, zvap) = '*'
     enddo
  endif
  if(zhenry .eq. 0) then
     nzonespecs = nzonespecs + 1
     zhenry = nzonespecs
     zonespecs(zhenry) = 'henry'
     do i = 1, nzones
        zonegrid(i, zhenry) = '*'
     enddo
  endif
  if(ztpor.eq. 0) then
     nzonespecs = nzonespecs + 1
     ztpor = nzonespecs
     zonespecs(ztpor) = 'tpor'
     do i = 1, nzones
        zonegrid(i, ztpor) = '*'
     enddo
  endif
  if(znzones .gt. nzones) then
     do i = nzones + 1, znzones
        write(zones(i), *) i
        zonegrid(i, zinit) = '*'
        zonegrid(i, zboun) = '*'
        zonegrid(i, zrock) = '*'
        zonegrid(i, zdisp) = '*'
        zonegrid(i, zsorp) = '*'
        zonegrid(i, ztpor) = '0'
        zonegrid(i, zldiff) = '0'
        zonegrid(i, zvdiff) = '0'
        zonegrid(i, ztime) = '*'
        zonegrid(i, zhenry) = '*'
     enddo
  endif
  ! ...Ensure that for each name given in zone, there is a matching name in the appropriate section.
  do i = 1, nzonespecs
     if(zonespecs(i) .eq. 'init') then
        do j = 1, nzones
           if(zonegrid(j, i) .eq. '*') then
              goto 114
           endif
           if(nwt .eq. 0) then
              write(ierr, *) 'Error:  Water type "', trim(zonegrid(j, i)), '" specified in assign, ', &
                   'but no water types defined in water.'
              goto 666
           endif
           do k = 1, nwt
              if(zonegrid(j, i) .eq. wtnames(k)) then
                 goto 114
              endif
           enddo
           write(ierr, *) 'Error:  Water type "', trim(zonegrid(j, i)), '" specified in assign but not in water.'
           goto 666
114        continue
        enddo
     elseif(zonespecs(i) .eq. 'boun') then
        do j = 1, nzones
           if(zonegrid(j, i) .eq. '*') then
              goto 176
           endif
           e = zonegrid(j, i)//'.'
           do while(index(e, '.') .ne. 0)
              f = e(1:index(e, '.') - 1)
              e = e(index(e, '.') + 1:len(e))
              do k = 1, nwt
                 if(f .eq. wtnames(k)) then
                    goto 175
                 endif
              enddo
              do k = 1, nrt
                 if(f .eq. rtnames(k)) then
                    goto 175
                 endif
              enddo
              do k = 1, nvt
                 if(f .eq. vtnames(k)) then
                    goto 175
                 endif
              enddo
              do k = 1, nht
                 if(f .eq. htnames(k)) then
                    goto 175
                 endif
              enddo
              write(ierr, *) 'Error:  Type "', trim(f), &
                   '" specified in assign but not in water, rock, vap, or henry.'
              goto 666
175           continue
           enddo
176        continue
        enddo
     elseif(zonespecs(i) .eq. 'rock') then
        do j = 1, nzones
           if(zonegrid(j, i) .eq. '*') then
              goto 115
           endif
           if(nrt .eq. 0) then
              write(ierr, *) 'Error:  Rock type "', trim(zonegrid(j, i)), '" specified in assign, ', &
                   'but no rock types defined in rock.'
              goto 666
           endif
           do k = 1, nrt
              if((zonegrid(j, i) .eq. rtnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
                 goto 115
              endif
           enddo
           write(ierr, *) 'Error:  Rock type "', trim(zonegrid(j, i)), '" specified in assign but not in rock.'
           goto 666
115        continue
        enddo
     elseif(zonespecs(i) .eq. 'vap') then
        do j = 1, nzones
           if(zonegrid(j, i) .eq. '*') then
              goto 149
           endif
           if(nrt .eq. 0) then
              write(ierr, *) 'Error:  Vapor type "', trim(zonegrid(j, i)), '" specified in assign, ', &
                   'but no vapor types defined in vap.'
              goto 666
           endif
           do k = 1, nvt
              if((zonegrid(j, i) .eq. vtnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
                 goto 149
              endif
           enddo
           write(ierr, *) 'Error:  Vapor type "', trim(zonegrid(j, i)), '" specified in assign but not in vap.'
           goto 666
149        continue
        enddo
     elseif(zonespecs(i) .eq. 'henry') then
        do j = 1, nzones
           if(zonegrid(j, i) .eq. '*') then
              goto 150
           endif
           if(nht .eq. 0) then
              write(ierr, *) 'Error:  Henry''s Law type "', trim(zonegrid(j, i)), &
                   '" specified in assign, but no Henry''s Law types defined in henry.'
              goto 666
           endif
           do k = 1, nht
              if((zonegrid(j, i) .eq. htnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
                 goto 150
              endif
           enddo
           write(ierr, *) 'Error:  henry type "', trim(zonegrid(j, i)), '" specified in assign but not in henry.'
           goto 666
150        continue
        enddo
     elseif(zonespecs(i) .eq. 'sorp') then
        do j = 1, nzones
           if(nrt .eq. 0) then
              goto 136
           endif
           do k = 1, nrt
              if((zonegrid(j, i) .eq. sorpnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
                 goto 136
              endif
           enddo
           write(ierr, *) 'Error:  Adsorption type "', trim(zonegrid(j, i)), &
                '" specified in assign but not in sorp.'
           goto 666
136        continue
        enddo
     elseif((zonespecs(i) .eq. 'ldisp') .or. (zonespecs(i) .eq. 'vdisp')) then
        do j = 1, nzones
           do k = 1, ndisp
              if(zonegrid(j, i) .eq. dispnames(k)) then
                 goto 125
              endif
           enddo
           write(ierr, *) 'Error:  Dispersivity model "', trim(zonegrid(j, i)), &
                '" specified in assign but not in disp.'
125        continue
        enddo
     endif
  enddo
  ! ...Convert diff columns to FEHM-compatible model numbers.
  do i = 1, nzones
     if((zonegrid(i, zldiff) .eq. '0') .or. (zonegrid(i, zldiff) .eq. '*') .or. (zonegrid(i, zldiff) .eq. 'con')) then
        zonegrid(i, zldiff) = '0'
     elseif((zonegrid(i, zldiff) .eq. '1') .or. (zonegrid(i, zldiff) .eq. 'mq')) then
        zonegrid(i, zldiff) = '1'
     elseif((zonegrid(i, zldiff) .eq. '2') .or. (zonegrid(i, zldiff) .eq. 'cw')) then
        zonegrid(i, zldiff) = '2'
     else
        write(ierr, *) 'Liquid diffusion model "', zonegrid(i, zldiff), '" not recognized in assign.'
        goto 666
     endif
     if((zonegrid(i, zvdiff) .eq. '0') .or. (zonegrid(i, zvdiff) .eq. '*') .or. (zonegrid(i, zvdiff) .eq. 'con')) then
        zonegrid(i, zvdiff) = '0'
     elseif((zonegrid(i, zvdiff) .eq. '1') .or. (zonegrid(i, zvdiff) .eq. 'mq')) then
        zonegrid(i, zvdiff) = '1'
     elseif((zonegrid(i, zvdiff) .eq. '2') .or. (zonegrid(i, zvdiff) .eq. 'cw')) then
        zonegrid(i, zvdiff) = '2'
     else
        write(ierr, *) 'Vapor diffusion model "', trim(zonegrid(i, zvdiff)), '" not recognized in assign.'
        goto 666
     endif
  enddo
  ! ...If H+ concentration is given in terms of pH, convert to [H+]
  do i = 1, nwtspecies
     if(wtspecies(i) .eq. 'pH') then
        wtspecies(i) = 'H+'
        do j = 1, nwt
           wtgrid(j, i) = 10 ** (-1.0 * wtgrid(j, i))
        enddo
     endif
  enddo
  ! ...Add an extra adsorption model with all parameters 0 for unspecified zones.
  do i = 1, anspecies
     lsorptypes(nsorp+1, i) = 0
     vsorptypes(nsorp+1, i) = 0
     a1l(nsorp+1, i) = 0
     a2l(nsorp+1, i) = 0
     bl(nsorp+1, i) = 0
     a1v(nsorp+1, i) = 0
     a2v(nsorp+1, i) = 0
     bv(nsorp+1, i) = 0
     sorpnames(nsorp+1) = '[default]'
  enddo
  ! ...Add an extra dispersion model with all parameters 0 for unspecified zones.
  lx(ndisp+1) = 0
  ly(ndisp+1) = 0
  lz(ndisp+1) = 0
  ll(ndisp+1) = 0
  lt(ndisp+1) = 0
  vx(ndisp+1) = 0
  vy(ndisp+1) = 0
  vz(ndisp+1) = 0
  vl(ndisp+1) = 0
  vt(ndisp+1) = 0
  dispnames(ndisp+1) = '[default]'
  ! ...Check that all components in water/rock/vap/henry appear in comp and are of the appropriate state.
  if(fwater .eqv. .true.) then
     do i = 1, nwtspecies
        do j = 1, nspecies
           if(species(j) .eq. wtspecies(i)) then
              if(states(j) .ne. 'a') then
                 write(ierr, *) 'Component "', trim(wtspecies(i)), '" in water is not aqueous.'
                 goto 666
              else
                 goto 169
              endif
           endif
        enddo
        write(ierr, *) 'Component "', trim(wtspecies(i)), '" specified in water but not in comp.'
        goto 666
169     continue
     enddo
  endif
  if(frock .eqv. .true.) then
     do i = 1, nrtspecies
        do j = 1, nspecies
           if(species(j) .eq. rtspecies(i)) then
              if(states(j) .ne. 's') then
                 write(ierr, *) 'Component "', trim(rtspecies(i)), '" in rock is not solid.'
                 goto 666
              else
                 goto 170
              endif
           endif
        enddo
        write(ierr, *) 'Component "', trim(rtspecies(i)), '" specified in rock but not in comp.'
        goto 666
170     continue
     enddo
  endif
  if(fvap .eqv. .true.) then
     do i = 1, nvtspecies
        do j = 1, nspecies
           if(species(j) .eq. vtspecies(i)) then
              if(states(j) .ne. 'v') then
                 write(ierr, *) 'Component "', trim(vtspecies(i)), '" in vap is not vapor.'
                 goto 666
              else
                 goto 171
              endif
           endif
        enddo
        write(ierr, *) 'Component "', trim(vtspecies(i)), '" specified in vap but not in comp.'
        goto 666
171     continue
     enddo
  endif
  if(fhenry .eqv. .true.) then
     do i = 1, nhtspecies
        do j = 1, nspecies
           if(species(j) .eq. htspecies(i)) then
              if(states(j) .ne. 'h') then
                 write(ierr, *) 'Component "', trim(htspecies(i)), &
                      '" in henry is not a Henry''s Law component.'
                 goto 666
              else
                 goto 172
              endif
           endif
        enddo
        write(ierr, *) 'Component "', trim(htspecies(i)), '" specified in henry but not in comp.'
        goto 666
172     continue
     enddo
  endif
  ! ...Ensure that hparam exists if there are Henry's Law components, and that hparam contains all Henry's Law components.
  if(fhparam .eqv. .false.) then
     do i = 1, nspecies
        if(states(i) .eq. 'h') then
           write(ierr, *) 'Error:  Henry''s Law components have been specified in comp, but hparam is not present.'
           goto 666
        endif
     enddo
  else
     do i = 1, nspecies
        if(states(i) .eq. 'h') then
           do j = 1, hnhtspecies
              if(hhtspecies(j) .eq. species(i)) then
                 goto 173
              endif
           enddo
           write(ierr, *) 'Error:  Henry''s Law component "', trim(species(i)), '" does not appear in hparam.'
           goto 666
        endif
173     continue
     enddo
  endif
  ! ...Ensure that everything in group is a component.
  if(fgroup .eqv. .true.) then
     do i = 1, gnspecies
        do j = 1, nspecies
           if(gspecies(i) .eq. species(j)) then
              goto 199
           endif
        enddo
        write(ierr, *) 'Error:  species "', trim(gspecies(i)), '" in group not found in comp.'
        goto 666
199     continue
     enddo
  endif
  ! ...Ensure that everything in print is a component or master species.
  if(fprint .eqv. .true.) then
     do i = 1, nprint
        if(printspecies(i) .eq. 'all') then
           nprint = nspecies + nmasters
           deallocate(printspecies)
           allocate(printspecies(nprint))
           do j = 1, nspecies
              printspecies(i) = species(i)
           enddo
           k = 1
           do j = 1, nspecies
              if(masters(i) .ne. '*') then
                 printspecies(nspecies + k) = masters(i)
                 k = k + 1
              endif
           enddo
           goto 202
        elseif(printspecies(i) .eq. 'none') then
           nprint = 0
           goto 202
        endif
        do j = 1, nspecies
           if((printspecies(i) .eq. species(j)) .or. (printspecies(i) .eq. masters(j))) then
              goto 200
           endif
        enddo
        write(ierr, *) 'Error:  Print specification "', printspecies(i), &
             '" in print is not a component or master species.'
200     continue
     enddo
202  continue
  else
     nprint = 0
  endif
  ! ...Ensure that all reactants and products (except the daughter in radioactive decay), substrate, electron acceptor, biomass, and icbioholder are all components.  TODO ...right?
  ! ...Ensure that all dscoefs are names of models in dist and sol.
  do i = 1, nrxns
     if(dscoefs(i) .eq. '*') then
        continue
     endif
     r = -20000.4
     read(dscoefs(i), *, err=204) r
204  if(r .ne. -20000.4) then
        goto 203
     endif
     if((rxntypes(i) .eq. 1) .or. (rxntypes(i) .eq. 2)) then
        do j = 1, ndistmodels
           if(dscoefs(i) .eq. distmodelnames(j)) then
              goto 203
           endif
        enddo
        write(ierr, *) 'Error:  Distribution model "', dscoefs(i), '" in reaction ', i, &
             ' not found in dist.'
        goto 666
     elseif((rxntypes(i) .eq. 7) .or. (rxntypes(i) .eq. 8)) then
        do j = 1, nsolmodels
           if(dscoefs(i) .eq. solmodelnames(j)) then
              goto 203
           endif
        enddo
        write(ierr, *) 'Error:  Solubility model "', dscoefs(i), '" in reaction ', i, &
             '" not found in sol.'
        goto 666
     endif
203  continue
  enddo
  ! If lookup is being used, use the database to create a matrix.

  ! Map all internal variables to FEHM's variables.
  ! Set values for trac.
  if(debug) write(iptty, *) 'Applying data for trac...'
  ctol = 1.e-06
  iret = 0
  nspeci = nspecies
  if(nspeci .ne. 0) then
     ntpp = n7 / nspeci
  else
     ntpp = 0
  endif
  allocate(spnam(nspeci))
  do nsp = 1, nspeci
     npt(nsp) = (nsp-1)*ntpp
     qcout(nsp) = 0.0
     qcin(nsp) = 0.0
     qcrxn(nsp) = 0.0
  enddo
  ! Groups 1-3 have already been set in header.  Groups 4-7 are set here.
  ! zvd 06/02/03 Set default transport porosity to bogus value, assign all unassigned
  ! nodes the value set for rock macro after all data has been read
  ps_trac = 10. ! ps_trac = 0.1
  if(ztpor .ne. 0) then
     do i = 1, n0
        if(zonegrid(izonef(i), ztpor) .eq. '*') then
           ! Use default value
        else
           read(zonegrid(izonef(i), ztpor), *) ps_trac(i)
        endif
     enddo
  endif
  do i = 1, nspecies
     spnam(i) = species(i)
  enddo
  if(dispmode .eq. 0) then
     ldsp = 0
  else
     ldsp = 1
  endif
  ! Set group 8	
  dispsame = 0	! 1
  ! Set groups 9-10
  do i = 1, nsorp + 1 !ndisp + 1  <-  FIXME
     do j = 1, nspeci
        if(dispmode .eq. 0) then
           tclx(j, i) = max(zero_t, lx(i))**2	! FIXME These variables are allocated by numsorp!
           tcly(j, i) = max(zero_t, ly(i))**2
           tclz(j, i) = max(zero_t, lz(i))**2
           tcvx(j, i) = max(zero_t, vx(i))**2
           tcvy(j, i) = max(zero_t, vy(i))**2
           tcvz(j, i) = max(zero_t, vz(i))**2
        else
           tclx(j, i) = max(zero_t, ll(i))**2
           tcly(j, i) = max(zero_t, lt(i))**2
           tcvx(j, i) = max(zero_t, vl(i))**2
           tcvy(j, i) = max(zero_t, vt(i))**2
        endif
     enddo
  enddo
  do i = 1, n0
     a = 0
     do j = 1, ndisp
        if(zonegrid(izonef(i), zdisp) .eq. dispnames(j)) then
           a = j
           exit
        endif
     enddo
     if((a .eq. 0) .and. (zonegrid(izonef(i), zdisp) .ne. '*')) then
        write(ierr, *) 'Internal error setting dispersivity data.'
        goto 667
     endif
     if(zonegrid(izonef(i), zdisp) .eq. '*') then
        itrcdsp(i) = ndisp + 1
     else
        itrcdsp(i) = a
     endif
  enddo
  icpnt = 1
  iimm = 1
  ivap = 1
  do nsp = 1, nspeci
     if((states(nsp) .eq. 'a') .or. (states(nsp) .eq. 'h')) then
        pcpnt(icpnt) = nsp
        cpntnam(icpnt) = species(nsp)
        icpnt = icpnt + 1
     elseif(states(nsp) .eq. 's') then
        pimm(iimm) = nsp
        immnam(iimm) = species(nsp)
        iimm = iimm + 1
     elseif(states(nsp) .eq. 'v') then
        pvap(ivap) = nsp
        vapnam(ivap) = species(nsp)
        ivap = ivap + 1
     else
        write(ierr, *) 'Internal error processing phase data.'
        goto 667
     endif
  enddo
  c = 1
  do i = 1, nspecies
     npn = npt(i)
     conc_read(i) = .true.
     ! Set group 11 values
     if(states(i) .eq. 'h') then
        hvliquid = 1
        hvvapor = 1
        icns(i) = 2
     elseif(states(i) .eq. 'v') then
        hvvapor = 1
        icns(i) = -1
     elseif(states(i) .eq. 's') then
        icns(i) = 0
     elseif(states(i) .eq. 'a') then
        hvliquid = 1
        icns(i) = 1
     else
        write(ierr, *) 'Internal error setting state data.'
        goto 667
     endif
     ! Set diffm
     do k = 1, ndiff
        if(diffspecs(k) .eq. species(i)) then
           do j = 1, nsorp + 1
              diffmfl(i, j) = ldiff(k)
              diffmfv(i, j) = vdiff(k)
           enddo
           goto 174
        endif
     enddo
     do j = 1, nsorp + 1
        diffmfl(i, j) = 0
        diffmfv(i, j) = 0
     enddo
174  continue
     ! Set group 15, 16 values
     an = an0
     do j = 1, n0
        k = (i - 1) * n0 + j
        inflag = .true.
        bflag = .true.
        if(zonegrid(izonef(j), zinit) .eq. '*') then
           inflag = .false.
           an((i - 1) * n0 + j) = an0
        endif
        if(zonegrid(izonef(j), zboun) .eq. '*') then
           bflag = .false.
           cnsk(k) = 0
           t1sk(k) = 0
           t2sk(k) = 0
        endif
        if(states(i) .eq. 'a') then
           if(inflag .eqv. .true.) then
              a = 0
              do o = 1, nwtspecies
                 if(wtspecies(o) .eq. species(i)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              do o = 1, nwt
                 if(wtnames(o) .eq. zonegrid(izonef(j), zinit)) then
                    d = o
                    exit
                 endif
              enddo
              if(a .eq. 0) then
                 an((i - 1) * n0 + j) = 0
              elseif(d .eq. 0) then
                 write(ierr, *) 'Error:  Could not find a match for species "', &
                      trim(species(i)), '" in water.'
                 goto 666
              else
                 an((i - 1) * n0 + j) = wtgrid(d, a)
              endif
           endif
           if(bflag .eqv. .true.) then
              a = 0
              do o = 1, nwtspecies
                 if(wtspecies(o) .eq. species(i)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              e = zonegrid(izonef(j), zboun)//'.'
              do o = 1, nwt
                 if(index(e, wtnames(o)//'.') .ne. 0) then
                    d = o
                    exit
                 endif
              enddo
              if(a .eq. 0) then
                 cnsk(k) = 0
              elseif(d .eq. 0) then
                 write(ierr, *) 'Internal error setting boundary water concentration data.'
                 goto 667
              else
                 cnsk(k) = wtgrid(d, a)
              endif
              e = zonegrid(izonef(j), ztime)
              if(e .eq. '*') then
                 t1sk(k) = 0
                 t2sk(k) = tims
              elseif(e .eq. '0') then
                 t1sk(k) = 0
                 t2sk(k) = 0
              elseif(index(e, '>') .eq. 0) then
                 read(e, *, err=250) t1sk(k)
                 t2sk(k) = t1sk(k) + 1
              else
                 f = e(1:index(e, '>') - 1)
                 g = e(index(e, '>') + 1:len(e))
                 if(len_trim(f) .eq. 0) then
                    t1sk(k) = 0
                 else
                    read(f, *, err=250, end=250) t1sk(k)
                 endif
                 if(len_trim(g) .eq. 0) then
                    t2sk(k) = tims
                 else
                    read(g, *, err=250, end=250) t2sk(k)
                 endif
              endif
              goto 251
250           write(ierr, *) 'Error setting aqueous tracer injection start/stop times.'
              goto 666
251           continue
           endif
        elseif(states(i) .eq. 's') then
           an((i - 1) * n0 + j) = an0
           if(bflag .eqv. .true.) then
              cnsk(k) = 0
              t1sk(k) = 0
              t2sk(k) = 0
           endif
        elseif(states(i) .eq. 'h') then
           if(inflag .eqv. .true.) then
              a = 0
              do o = 1, nht
                 if(htnames(o) .eq. zonegrid(izonef(j), zhenry)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              do o = 1, nhtspecies
                 if(htspecies(o) .eq. species(i)) then
                    d = o
                    exit
                 endif
              enddo
              if(d .eq. 0) then
                 an((i - 1) * n0 + j) = 0
              elseif(a .eq. 0) then
                 write(ierr, *) 'Internal error setting Henry''s Law components data.'
                 goto 667
              else
                 an((i - 1) * n0 + j) = htgrid(a, d)
              endif
           endif
           if(bflag .eqv. .true.) then
              a = 0
              do o = 1, nhtspecies
                 if(htspecies(o) .eq. species(i)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              e = zonegrid(izonef(j), zboun)//'.'
              do o = 1, nht
                 if(index(e, htnames(o)//'.') .ne. 0) then
                    d = o
                    exit
                 endif
              enddo
              if(a .eq. 0) then
                 cnsk(k) = 0
              elseif(d .eq. 0) then
                 write(ierr, *) 'Internal error setting boundary Henry''s concentration data.'
                 goto 667
              else
                 cnsk(k) = htgrid(d, a)
              endif
              e = zonegrid(izonef(j), ztime)
              if(e .eq. '*') then
                 t1sk(k) = 0
                 t2sk(k) = tims
              elseif(e .eq. '0') then
                 t1sk(k) = 0
                 t2sk(k) = 0
              elseif(index(e, '>') .eq. 0) then
                 read(e, *, err=252) t1sk(k)
                 t2sk(k) = t1sk(k) + 1
              else
                 f = e(1:index(e, '>') - 1)
                 g = e(index(e, '>') + 1:len(e))
                 if(f .eq. '') then
                    t1sk(k) = 0
                 else
                    read(f, *, err=252, end=252) t1sk(k)
                 endif
                 if(g .eq. '') then
                    t2sk(k) = tims
                 else
                    read(g, *, err=252, end=252) t2sk(k)
                 endif
              endif
              goto 253
252           write(ierr, *) 'Error setting Henry''s Law tracer injection start/stop times.'
              goto 666
253           continue
           endif
        elseif(states(i) .eq. 'v') then
           if(inflag .eqv. .true.) then
              a = 0
              do o = 1, nvt
                 if(vtnames(o) .eq. zonegrid(izonef(j), zvap)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              do o = 1, nvtspecies
                 if(vtspecies(o) .eq. species(i)) then
                    d = o
                    exit
                 endif
              enddo
              if(d .eq. 0) then
                 an((i - 1) * n0 + j) = 0
              elseif(a .eq. 0) then
                 write(ierr, *) 'Internal error setting vapor data.'
                 goto 667
              else
                 an((i - 1) * n0 + j) = vtgrid(a, d)
              endif
           endif
           if(bflag .eqv. .true.) then
              a = 0
              do o = 1, nvtspecies
                 if(vtspecies(o) .eq. species(i)) then
                    a = o
                    exit
                 endif
              enddo
              d = 0
              e = zonegrid(izonef(j), zboun)//'.'
              do o = 1, nvt
                 if(index(e, vtnames(o)//'.') .ne. 0) then
                    d = o
                    exit
                 endif
              enddo
              if(a .eq. 0) then
                 cnsk(k) = 0
              elseif(d .eq. 0) then
                 write(ierr, *) 'Internal error setting boundary vapor concentration data.'
                 goto 667
              else
                 cnsk(k) = vtgrid(d, a)
              endif
              e = zonegrid(izonef(j), ztime)
              if(e .eq. '*') then
                 t1sk(k) = 0
                 t2sk(k) = tims
              elseif(e .eq. '0') then
                 t1sk(k) = 0
                 t2sk(k) = 0
              elseif(index(e, '>') .eq. 0) then
                 read(e, *, err=254) t1sk(k)
                 t2sk(k) = t1sk(k) + 1
              else
                 f = e(1:index(e, '>') - 1)
                 g = e(index(e, '>') + 1:len(e))
                 if(f .eq. '') then
                    t1sk(k) = 0
                 else
                    read(f, *, err=254, end=254) t1sk(k)
                 endif
                 if(g .eq. '') then
                    t2sk(k) = tims
                 else
                    read(g, *, err=254, end=254) t2sk(k)
                 endif
              endif
              goto 255
254           write(ierr, *) 'Error setting vapor tracer injection start/stop times.'
              goto 666
255           continue
           endif
        endif
        if(cnsk(k) .lt. 0) then
           if(t2sk(k) .lt. 0) then
              write(ierr, *) 'Error:  Solute accumulation option specified by end time', t2sk(k), &
                   'cannot be used with constant concentration', cnsk(k), '.'
              goto 666
           endif
           pcnsk(k) = -1.0
        elseif(t2sk(k) .lt. 0) then
           t2sk(k) = abs(t2sk(k))
           pcnsk(k) = 1.0
        endif
        anlo(j) = an(j)
        e = zonegrid(izonef(j), zsorp)
        if(e .eq. '*') then
           itrc((i-1)*n0 + j) = nsorp + 1
        else
           a = 0
           do k = 1, nsorp
              if(sorpnames(k) .eq. e) then
                 a = k
                 exit
              endif
           enddo
           if(a .eq. 0) then
              write(ierr, *) 'Internal error setting adsorption data.'
              goto 667
           endif
           itrc((i-1)*n0 + j) = a
        endif
     enddo
     if(states(i) .eq. 'h') then
        ! Set group 14 values
        d = 0
        do o = 1, hnhtspecies
           if(hhtspecies(o) .eq. species(i)) then
              d = o
              exit
           endif
        enddo
        if(d .eq. 0) then
           write(ierr, *) 'Internal error processing Henry''s Law models.'
           goto 667
        endif
        henry_model(i) = hmodels(d)
        a = 0
        do o = 1, hnhtspecies
           if(hhtspecies(o) .eq. species(i)) then
              a = o
              exit
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Internal error setting Henry''s Law parameters.'
           goto 667
        endif
        if(henry_model(i) .eq. 1) then
           hawwa(c, 1) = ah(a)
           hawwa(c, 2) = dhh(a)
           hawwa(c, 3) = 0
           hawwa(c, 4) = 0
           hawwa(c, 5) = 0
        elseif(henry_model(i) .eq. 2) then
           hawwa(c, 1) = ah1(a)
           hawwa(c, 2) = ah2(a)
           hawwa(c, 3) = ah3(a)
           hawwa(c, 4) = ah4(a)
           hawwa(c, 5) = ah5(a)
        elseif(henry_model(i) .eq. 3) then
           hawwa(c, 1) = ah(a)
           hawwa(c, 2) = hh(a)
           hawwa(c, 3) = 0
           hawwa(c, 4) = 0
           hawwa(c, 5) = 0
        else
           write(ierr, *) 'Internal error setting Henry''s Law parameters.'
           goto 667
        endif
        c = c + 1
     endif
     ! Set group 12, 13 values
     if(states(i) .eq. 's') then
        do j = 1, nsorp + 1
           iadsfl(i, j) = 0
           a1adfl(i, j) = 0
           a2adfl(i, j) = 0
           betadfl(i, j) = 1
           iadsfv(i, j) = 0
           a1adfv(i, j) = 0
           a2adfv(i, j) = 0
           betadfv(i, j) = 1
        enddo
     elseif((states(i) .eq. 'a') .or. (states(i) .eq. 'v') .or. (states(i) .eq. 'h')) then
        ! Set k to be the index of the current species in anspecies.
        k = 0
        do a = 1, anspecies
           if(species(i) .eq. sorpspecs(a)) then
              k = a
              exit
           endif
        enddo
        if(k .ne. 0) then
           do j = 1, nsorp + 1
              iadsfl(i, j) = lsorptypes(j, k)
              iadsfv(i, j) = vsorptypes(j, k)
              a1adfl(i, j) = a1l(j, k)
              a2adfl(i, j) = a2l(j, k)
              betadfl(i, j) = bl(j, k)
              a1adfv(i, j) = a1v(j, k)
              a2adfv(i, j) = a2v(j, k)
              betadfv(i, j) = bv(j, k)
           enddo
        else
           do j = 1, nsorp + 1
              iadsfl(i, j) = 0
              a1adfl(i, j) = 0
              a2adfl(i, j) = 0
              betadfl(i, j) = 1
              iadsfv(i, j) = 0
              a1adfv(i, j) = 0
              a2adfv(i, j) = 0
              betadfv(i, j) = 1
           enddo
        endif
     else
        write(ierr, *) 'Internal error setting state data.'
        goto 667
     endif
  enddo
  ! End of trac value-setting.
  ! Set values for rxn.
  if(rxn_flag .eq. 0) then
     ncpntprt = ncpnt
     do i = 1, ncpnt
        cpntprt(i) = i
     enddo
     !		nimmprt = nimm		! TODO If reactions aren't disabled, this code kills the simulation.  Why?
     !		do i = 1, nimm
     !			immprt(i) = i
     !		enddo
     !		nvapprt = nvap
     !		do i = 1, nvap
     !			vapprt(i) = i
     !		enddo
     !		ngroups = 0
     goto 6000
  endif
  if(debug) write(iptty, *) 'Applying data for rxn...'
  rxnon = 1
  rxnnaqueous = 0
  rxnnsolid = 0
  rxnnvapor = 0
  allocate(rxnaqueous(SPEC_MAX * 2))
  allocate(rxnsolid(SPEC_MAX))
  allocate(rxnvapor(SPEC_MAX))
  rxnaqueous = '*'
  rxnsolid = '*'
  rxnvapor = '*'
  do i = 1, nspecies
     if((states(i) .eq. 'a') .or. (states(i) .eq. 'h')) then
        rxnnaqueous = rxnnaqueous + 1
        rxnaqueous(rxnnaqueous) = species(i)
     elseif(states(i) .eq. 'v') then
        rxnnvapor = rxnnvapor + 1
        rxnvapor(rxnnvapor) = species(i)
     elseif(states(i) .eq. 's') then
        rxnnsolid = rxnnsolid + 1
        rxnsolid(rxnnsolid) = species(i)
     endif
  enddo
  ngroups = gngroups
  ! Set group 3 - 6 values
  ncpnt = 0
  nimm = 0
  nvap = 0
  cpntprt = 0
  immprt = 0
  vapprt = 0
  ncpntprt = nprint
  do i = 1, nspecies
     if((states(i) .eq. 'a') .or. (states(i) .eq. 'h')) then
        ncpnt = ncpnt + 1
        cpntnam(ncpnt) = species(i)
        ifxconc(ncpnt) = 0
        cpntgs(ncpnt) = guesses(i)
        do j = 1, nprint
           if(printspecies(j) .eq. species(i)) then
              cpntprt(ncpnt) = 1
           endif
        enddo
        !idcpnt(ncpnt) = ncpnt
     elseif(states(i) .eq. 'v') then
        nvap = nvap + 1
        vapnam(nvap) = species(i)
        do j = 1, nprint
           if(printspecies(j) .eq. species(i)) then
              vapprt(nvap) = 1
           endif
        enddo
        !idvap(nvap) = nimm
     elseif(states(i) .eq. 's') then
        nimm = nimm + 1
        immnam(nimm) = species(i)
        do j = 1, nprint
           if(printspecies(j) .eq. species(i)) then
              immprt(nimm) = 1
           endif
        enddo
        !idimm(nimm) = nimm
     else
        write(ierr, *) 'Internal error setting species names.'
        goto 667
     endif
241  continue
  enddo
  ! Set group 9 and 10 values
  if(logten .eqv. .true.) then
     logkeq = 1
  else
     logkeq = 0
  endif
  if(fequi .eqv. .true.) then
     do i = 1, ncomplexes
        do j = 1, encomplexes
           if(complexes(i) .eq. ecomplexes(j)) then
              ckeq(i + 100) = equconstants(j)
              heq(i + 100, 1) = enthalpies(j)
              goto 249
           endif
        enddo
        ckeq(i + 100) = 1	! TODO Need a good default value for equilibrium constant
        heq(i + 100, 1) = 0
249     continue
     enddo
  else
     do i = 1, ncomplexes
        ckeq(i + 100) = 1
        heq(i + 100, 1) = 0
     enddo
  endif
  ! Set group 11 values
  cplxprt = 0
  ncplx = ncomplexes
  do i = 1, ncomplexes
     cplxnam(i + 100) = complexes(i)
     do j = 1, SPEC_MAX
        if(complexcontents(i, j) .eq. '*') then
           exit
        endif
        do k = 1, nspecies
           if((masters(k) .eq. '*') .or. (masters(k) .eq. '')) then
              continue
           endif
           if(complexcontents(i, j) .eq. masters(k)) then
              goto 240
           endif
        enddo
        write(ierr, *) 'Error:  Master species "', complexcontents(i, j), '" not found in comp.'
        goto 666
240     continue
     enddo
     do j = 1, nprint
        if(complexes(i) .eq. printspecies(j)) then
           cplxprt(i) = 1
           exit
        endif
     enddo
  enddo
  spstoic = 0
  do i = 1, ncomplexes
     o = 0
     do j = 1, nspecies
        if((states(j) .eq. 'a') .or. (states(j) .eq. 'h')) then
           o = o + 1
        endif
        do k = 1, SPEC_MAX
           if(complexcontents(i, k) .eq. '*') then
              exit
           endif
           if(masters(j) .eq. complexcontents(i, k)) then
              spstoic(i + 100, o) = complexstoich(i, k)
              exit
           endif
        enddo
     enddo
  enddo
  ! Copied from read_rxn
  ii=0
  do ic = 1,ncpnt
     if(cpntprt(ic).eq.0)then
        ii = ii+1
        cpntprt(ii)=ic
     endif
  enddo
  ncpntprt=ii
  ii=100
  do ix=101,ncplx+100
     if(cplxprt(ix).eq.0)then
        ii= ii+1
        cplxprt(ii)=ix
     endif
  enddo
  ncplxprt=ii-100
  ii=0
  do im=1,nimm
     if(immprt(im).eq.0)then
        ii=ii+1
        immprt(ii)=im
     endif
  enddo
  nimmprt=ii
  ii=0
  do iv=1,nvap
     if(vapprt(iv).eq.0)then
        ii=ii+1
        vapprt(ii)=iv
     endif
  enddo
  nvapprt=ii
  ! Set group 12 for each reaction
  numrxn = nrxns
  do i = 1, numrxn
     idrxn(i) = rxntypes(i)
  enddo
  ! Set reaction-specific groups
  ! GAAH I'M DROWNING IN SEQUENTIAL SEARCHES!
  do i = 1, nrxns
     if(rxntypes(i) .eq. 1) then
        naqsp(i) = 1
        nimsp(i) = 1
        nivsp(i) = 0
        do j = 1, nspecies
           if((reactants(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
              do o = 1, nimm
                 if(immnam(o) .eq. species(j)) then
                    irxnim(i, 1) = o
                 endif
              enddo
              do k = 1, nspecies
                 if(products(i, 1) .eq. masters(k)) then
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, 1) = o
                       endif
                    enddo
                    goto 243
                 endif
              enddo
              do k = 1, ncomplexes
                 if(products(i, 1) .eq. complexes(k)) then
                    irxnic(i, 1) = k + 100
                    goto 243
                 endif
              enddo
              goto 242
           else if((products(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
              do o = 1, nimm
                 if(immnam(o) .eq. species(j)) then
                    irxnim(i, 1) = o
                 endif
              enddo
              do k = 1, nspecies
                 if(reactants(i, 1) .eq. masters(k)) then
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, 1) = o
                       endif
                    enddo
                    goto 243
                    goto 243
                 endif
              enddo
              do k = 1, ncomplexes
                 if(reactants(i, 1) .eq. complexes(k)) then
                    irxnic(i, 1) = k + 100
                    goto 243
                 endif
              enddo
              goto 242
           endif
        enddo
242     write(ierr, *) 'Error:  Type 1 reaction ', i, ' requires exactly one aqueous complex or master ', &
             'species and one solid component.'
        goto 666
243     read(dscoefs(i), *, err=206) r
        ckeqlb(i) = r
        goto 207
206     a = 0
        do j = 1, ndistmodels
           if(dscoefs(i) .eq. distmodelnames(j)) then
              a = j
              exit
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Internal error setting distribution coefficient information.'
           goto 667
        endif
        j = 0
        do
           if(distmodels(a, j + 1, 1) .eq. -20000.4) then
              exit
           endif
           j = j + 1
           eqtemp(j) = distmodels(a, j, 1);
           lkeq(j) = distmodels(a, j, 2);
        enddo
        call lstsq(j, lkeq, eqtemp)
        do j = 1, 3
           tcoeff(i, j) = lkeq(i)
        enddo
        ckeqlb(i)=tcoeff(i,1)+tcoeff(i,2)*25+tcoeff(i,3)*25**2
207     ckmtrn(i) = rates(i)
     elseif(rxntypes(i) .eq. 2) then
        naqsp(i) = 1
        nimsp(i) = 1
        nivsp(i) = 0
        do j = 1, nspecies
           if((reactants(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
              do o = 1, nimm
                 if(immnam(o) .eq. species(j)) then
                    irxnim(i, 1) = o
                 endif
              enddo
              do k = 1, nspecies
                 if(products(i, 1) .eq. masters(k)) then
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, 1) = o
                       endif
                    enddo
                    goto 243
                    goto 245
                 endif
              enddo
              do k = 1, ncomplexes
                 if(products(i, 1) .eq. complexes(k)) then
                    irxnic(i, 1) = k + 100
                    goto 245
                 endif
              enddo
              goto 244
           else if((products(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
              do o = 1, nimm
                 if(immnam(o) .eq. species(j)) then
                    irxnim(i, 1) = o
                 endif
              enddo
              do k = 1, nspecies
                 if(reactants(i, 1) .eq. masters(k)) then
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, 1) = o
                       endif
                    enddo
                    goto 243
                    goto 245
                 endif
              enddo
              do k = 1, ncomplexes
                 if(reactants(i, 1) .eq. complexes(k)) then
                    irxnic(i, 1) = k + 100
                    goto 245
                 endif
              enddo
              goto 244
           endif
        enddo
244     write(ierr, *) 'Error:  Type 2 reaction ', i, ' requires exactly one aqueous complex or master ', &
             'species and one solid component.'
        goto 666
245     read(dscoefs(i), *, err=208) r
        ckeqlb(i) = r
        goto 209
208     a = 0
        do j = 1, ndistmodels
           if(dscoefs(i) .eq. distmodelnames(j)) then
              a = j
              exit
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Internal error setting distribution coefficient information.'
           goto 667
        endif
        j = 0
        do
           if(distmodels(a, j + 1, 1) .eq. -20000.4) then
              exit
           endif
           j = j + 1
           eqtemp(j) = distmodels(a, j, 1);
           lkeq(j) = distmodels(a, j, 2);
        enddo
        call lstsq(j, lkeq, eqtemp)
        do j = 1, 3
           tcoeff(i, j) = lkeq(i)
        enddo
        ckeqlb(i)=tcoeff(i,1)+tcoeff(i,2)*25+tcoeff(i,3)*25**2
209     ckmtrn(i) = rates(i)
        simmmx(i) = xcoef(i)
     elseif(rxntypes(i) .eq. 3) then
        nimsp(i) = 0
        naqsp(i) = 0
        nivsp(i) = 0
        do j = 1, SPEC_MAX
           if(reactants(i, j) .eq. '*') then
              exit
           endif
           do k = 1, nspecies
              if((species(k) .eq. reactants(i, j)) .or. (masters(k) .eq. reactants(i, j))) then
                 if((states(k) .eq. 'a') .or. (states(k) .eq. 'h')) then
                    naqsp(i) = naqsp(i) + 1
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, naqsp(i)) = o
                       endif
                    enddo
                    sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
                 elseif(states(k) .eq. 's') then
                    nimsp(i) = nimsp(i) + 1
                    do o = 1, nimm
                       if(immnam(o) .eq. species(k)) then
                          irxnim(i, nimsp(i)) = o
                       endif
                    enddo
                    stimirrv(i, nimsp(i)) = stoichiometries(i, 1, j)

                 elseif(states(k) .eq. 'v') then
                    nivsp(i) = nivsp(i) + 1
                    do o = 1, nvap
                       if(vapnam(o) .eq. species(k)) then
                          irxniv(i, nivsp(i)) = o
                       endif
                    enddo
                    stivirrv(i, nivsp(i)) = stoichiometries(i, 1, j)
                 endif
                 goto 210
              endif
           enddo
           do k = 1, ncomplexes
              if(complexes(k) .eq. reactants(i, j)) then
                 naqsp(i) = naqsp(i) + 1
                 irxnic(i, naqsp(i)) = k + 100
                 sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
                 goto 210
              endif
           enddo
           write(ierr, *) 'Error setting reactant information.'
           goto 666
210        continue
        enddo
        do j = 1, SPEC_MAX
           if(products(i, j) .eq. '*') then
              exit
           endif
           do k = 1, nspecies
              if((species(k) .eq. products(i, j)) .or. (masters(k) .eq. products(i, j))) then
                 if((states(k) .eq. 'a') .or. (states(k) .eq. 'h')) then
                    naqsp(i) = naqsp(i) + 1
                    do o = 1, ncpnt
                       if(cpntnam(o) .eq. species(k)) then
                          irxnic(i, naqsp(i)) = o
                       endif
                    enddo
                    sticirrv(i, naqsp(i)) = - stoichiometries(i, 2, j)
                 elseif(states(k) .eq. 's') then
                    nimsp(i) = nimsp(i) + 1
                    do o = 1, nimm
                       if(immnam(o) .eq. species(k)) then
                          irxnim(i, nimsp(i)) = o
                       endif
                    enddo
                    stimirrv(i, nimsp(i)) = - stoichiometries(i, 2, j)
                 elseif(states(k) .eq. 'v') then
                    nivsp(i) = nivsp(i) + 1
                    do o = 1, nvap
                       if(vapnam(o) .eq. species(k)) then
                          irxniv(i, nivsp(i)) = o
                       endif
                    enddo
                    stivirrv(i, nivsp(i)) = - stoichiometries(i, 2, j)
                 endif
                 goto 211
              endif
           enddo
           do k = 1, ncomplexes
              if(complexes(k) .eq. reactants(i, j)) then
                 naqsp(i) = naqsp(i) + 1
                 irxnic(i, naqsp(i)) = k + 100
                 sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
                 goto 211
              endif
           enddo
           write(ierr, *) 'Error setting product information.'
           goto 666
211        continue
        enddo
        kfor(i) = rates(i)
        krev(i) = rates2(i)
     elseif(rxntypes(i) .eq. 4) then
        naqsp(i) = 0
        nimsp(i) = 0
        nivsp(i) = 0
        read(bioparams(i, 1), *, end=218, err=230) f, e
        goto 219
218     e = f
        f = '1'
        read(f, *) o
219     do j = 1, nspecies
           if(e .eq. species(j)) then
              if((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
                 write(ierr, *) 'Error:  The substrate "', e, '" in type 4 reaction ', &
                      i, ' must be aqueous or Henry''s Law.'
                 goto 666
              endif
              if(o .ne. 1) then
                 write(ierr, *) 'Error:  The stoichiometry of the substrate "', e, &
                      '" in type 4 reaction ', i, ' must be 1.'
              endif
              naqsp(i) = naqsp(i) + 1
              irxnic(i, naqsp(i)) = j
           endif
        enddo
        if(naqsp(i) .eq. 0) then
           write(*, *) 'Error:  Substrate "', e, '" in type 4 reaction ', i, ' not found in comp.'
           goto 666
        endif
        read(bioparams(i, 2), *, end=223, err=230) f, e
        goto 224
223     e = f
        f = '1'
        read(f, *) o
224     do j = 1, nspecies
           if(e .eq. species(j)) then
              if((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
                 write(ierr, *) 'Error:  The electron acceptor "', e, '" in type 4 reaction ', &
                      i, ' must be aqueous or Henry''s Law.'
                 goto 666
              endif
              naqsp(i) = naqsp(i) + 1
              irxnic(i, naqsp(i)) = j
              biofac(j) = o
           endif
        enddo
        if(naqsp(i) .eq. 1) then
           write(*, *) 'Error:  Electron acceptor "', e, '" in type 4 reaction ', i, &
                ' not found in comp.'
           goto 666
        endif
        read(bioparams(i, 3), *, end=225, err=230) f, e
        goto 226
225     e = f
        f = '1'
        read(f, *) o
226     do j = 1, nspecies
           if(e .eq. species(j)) then
              if((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
                 write(ierr, *) 'Error:  The biomass "', e, '" in type 4 reaction ', &
                      i, ' must be solid.'
                 goto 666
              endif
              nimsp(i) = nimsp(i) + 1
              irxnim(i, nimsp(i)) = j
           endif
        enddo
        if(nimsp(i) .eq. 0) then
           write(*, *) 'Error:  Biomass "', e, '" in type 4 reaction ', i, ' not found in comp.'
           goto 666
        endif
        do j = 1, SPEC_MAX
           if(reactants(i, j) .eq. '*') then
              exit
           endif
           naqsp(i) = naqsp(i) + 1
           if(naqsp(i) .gt. 5) then
              write(ierr, *) 'Error:  Type 4 reaction ', i, &
                   ' can only have a total of three extra reactants and products.'
              goto 666
           endif
           do k = 1, nspecies
              if(reactants(i, j) .eq. species(k)) then
                 if((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
                    write(ierr, *) 'Error:  Extra reactant "', reactants(i, j), '" in type ', &
                         '4 reaction ', i, ' must be aqueous or Henry''s Law.'
                    goto 666
                 endif
                 irxnic(i, naqsp(i)) = k
                 if(naqsp(i) .eq. 3) then
                    hfac(i) = stoichiometries(i, 1, j)
                 elseif(naqsp(i) .eq. 4) then
                    carbfac(i) = stoichiometries(i, 1, j)
                 elseif(naqsp(i) .eq. 5) then
                    ammfac(i) = stoichiometries(i, 1, j)
                 else
                    write(ierr, *) 'Internal error setting biodegradation extra species data.'
                    goto 667
                 endif
              endif
           enddo
        enddo
        goto 231
230     write(ierr, *) 'Error:  Bad format for species in type 4 reaction ', i, '.'
        goto 666
231     do j = 1, SPEC_MAX
           if(products(i, j) .eq. '*') then
              exit
           endif
           naqsp(i) = naqsp(i) + 1
           if(naqsp(i) .gt. 5) then
              write(ierr, *) 'Error:  Type 4 reaction ', i, ' can only have a total of three ', &
                   'extra reactants and products.'
              goto 666
           endif
           do k = 1, nspecies
              if(products(i, j) .eq. species(k)) then
                 if((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
                    write(ierr, *) 'Error:  Extra product "', reactants(i, j), '" in type ', &
                         '4 reaction ', i, ' must be aqueous or Henry''s Law.'
                    goto 666
                 endif
                 irxnic(i, naqsp(i)) = k
                 if(naqsp(i) .eq. 3) then
                    hfac(i) = - stoichiometries(i, 1, j)
                 elseif(naqsp(i) .eq. 4) then
                    carbfac(i) = - stoichiometries(i, 1, j)
                 elseif(naqsp(i) .eq. 5) then
                    ammfac(i) = - stoichiometries(i, 1, j)
                 else
                    write(ierr, *) 'Internal error setting biodegradation extra species data.'
                    goto 667
                 endif
              endif
           enddo
        enddo
        read(bioparams(i, 4), *) ckc(i)
        read(bioparams(i, 5), *) cka(i)
        decay(i) = rates2(i)
        read(bioparams(i, 6), *) phthresh(i)
        qm(i) = rates(i)
        read(bioparams(i, 7), *) yield(i)
        read(bioparams(i, 8), *) xminit(i)
        nbiofrm(i) = 0
        do j = 1, SPEC_MAX
           if(icbioholder(i, j) .eq. '*') then
              exit
           endif
           do k = 1, nspecies
              if(icbioholder(i, j) .eq. species(k)) then
                 if((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
                    write(ierr, *) 'Error:  Biodegradable form "', icbioholder(i, j), &
                         '" in type 4 reaction ', i, ' must be aqueous or ', &
                         'Henry''s Law.'
                    goto 666
                 endif
                 nbiofrm(i) = nbiofrm(i) + 1
                 icbio(i, nbiofrm(i)) = k
              endif
           enddo
           if(nbiofrm(i) .ne. j) then
              write(ierr, *) 'Error:  Biodegradable form "', icbioholder(i, j), '" in type 4 reaction ', &
                   i, ' not found in comp.'
           endif
        enddo
     elseif(rxntypes(i) .eq. 5) then
        ckmtrn(i) = rates(i)
        ckmtrn(j) = log(2.0)/(ckmtrn(j)*(365.25*24))
        a = 0
        do j = 1, nspecies
           if(reactants(i, 1) .eq. species(j)) then
              a = j
              exit	
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Error:  Component "', reactants(i, 1), '" in reaction ', i, &
                ' not found in comp.'
           goto 666
        endif
        c = 0
        do j = 1, nspecies
           if(products(i, 1) .eq. '*') then
              exit
           elseif(products(i, 1) .eq. species(j)) then
              c = j
              exit	
           endif
        enddo
        if(states(a) .ne. states(c)) then
           write(ierr, *) 'The parent and daughter in type 5 reaction ', i, &
                ' must be in the same state.'
           goto 666
        endif
        if((states(a) .eq. 'a') .or. (states(a) .eq. 'h')) then
           naqsp(i) = 2
           irxnic(i, 1) = a
           if(c .eq. 0) then
              naqsp(i) = 1
           else
              irxnic(i, 2) = c
           endif
        elseif(states(a) .eq. 's') then
           nimsp(i) = 2
           irxnim(i, 1) = a
           if(c .eq. 0) then
              nimsp(i) = 1
           else
              irxnim(i, 2) = c
           endif
        elseif(states(a) .eq. 'v') then
           nivsp(i) = 2
           irxniv(i, 1) = a
           if(c .eq. 0) then
              nivsp(i) = 1
           else
              irxniv(i, 2) = c
           endif
        else
           write(ierr, *) 'Internal error setting decay reaction state.'
           goto 667
        endif
     elseif((rxntypes(i) .eq. 7) .or. (rxntypes(i) .eq. 8)) then
        nimsp(i) = 1
        naqsp(i) = 0
        nivsp(i) = 0
        a = 0
        do j = 1, nspecies
           if(reactants(i, 1) .eq. species(j)) then
              a = j
              exit
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Error:  Reactant "', reactants(i, 1), '" in type ', rxntypes(i), &
                ' reaction ', i, ' not found in comp.'
           goto 666
        endif
        if(states(a) .ne. 's') then
           if (reactants(i, 2) .ne. '*') then
              goto 212
           endif
           irxnim(i, 1) = a
           pdstim = stoichiometries(i, 1, 1)
           do j = 1, SPEC_MAX
              if(products(i, j) .eq. '*') then
                 exit
              endif
              do k = 1, nspecies
                 if(products(i, j) .eq. masters(k)) then
                    if(states(k) .ne. 'a') then
                       write(ierr, *) 'Error:  One side of type ', rxntypes(i), &
                            ' reaction ', i, ' must consist entirely ', &
                            'of aqueous master species.'
                       goto 666
                    endif
                    naqsp(i) = naqsp(i) + 1
                    irxnic(i, naqsp(i)) = k
                    pdstic(i, naqsp(i)) = - stoichiometries(i, 2, j)
                    goto 213
                 endif
              enddo
              write(ierr, *) 'Error:  Product "', products(i, j), '" in type ', rxntypes(i), &
                   ' reaction ', i, ' not found in comp.'
213           continue
           enddo
        else
           c = 0
           do j = 1, nspecies
              if(products(i, 1) .eq. species(j)) then
                 c = j
                 exit
              endif
           enddo
           if(c .eq. 0) then
              write(ierr, *) 'Error:  Product "', reactants(i, 1), '" in type ', rxntypes(i), &
                   ' reaction ', i, ' not found.'
              goto 666
           endif
           if(states(c) .eq. 's') then
              if(products(i, 2) .ne. '*') then
                 goto 212
              endif
              irxnim(i, 1) = c
              pdstim = stoichiometries(i, 2, 1) ! pdstim(i, 1)
              do j = 1, SPEC_MAX
                 if(reactants(i, j) .eq. '*') then
                    exit
                 endif
                 do k = 1, nspecies
                    if(reactants(i, j) .eq. masters(k)) then
                       if(states(k) .ne. 'a') then
                          write(ierr, *) 'Error:  One side of type ', rxntypes(i), &
                               ' reaction ', i, ' must consist ', &
                               'entirely of aqueous master species.'
                          goto 666
                       endif
                       naqsp(i) = naqsp(i) + 1
                       irxnic(i, naqsp(i)) = k
                       pdstic(i, naqsp(i)) = stoichiometries(i, 1, j)
                       goto 227
                    endif
                 enddo
                 write(ierr, *) 'Error:  Reactant "', reactants(i, j), '" in type ', rxntypes(i), &
                      ' reaction ', i, ' not found in comp.'
227              continue
              enddo
           else
              goto 212
           endif
        endif
        goto 232
212     write(ierr, *) 'Error:  One side of type ', rxntypes(i), ' reaction ', i, &
             ' must contain only a single solid species.'
        goto 666
232     read(dscoefs(i), *, err=228) r
        ckeqlb(i) = r
        goto 229
228     a = 0
        do j = 1, nsolmodels
           if(dscoefs(i) .eq. solmodelnames(j)) then
              a = j
              exit
           endif
        enddo
        if(a .eq. 0) then
           write(ierr, *) 'Internal error setting distribution coefficient information.'
           goto 667
        endif
        j = 0
        do
           if(solmodels(a, j + 1, 1) .eq. -20000.4) then
              exit
           endif
           j = j + 1
           eqtemp(j) = solmodels(a, j, 1);
           lkeq(j) = solmodels(a, j, 2);
        enddo
        call lstsq(j, lkeq, eqtemp)
        do j = 1, 3
           tcoeff(i, j) = lkeq(i)
        enddo
        ckeqlb(i)=tcoeff(i,1)+tcoeff(i,2)*25+tcoeff(i,3)*25**2
229     ckmtrn(i) = rates(i)
        if(rxntypes(i) .eq. 8) then
           mw_mineral(irxnim(i, 1)) = porchange(i, 1)
           rho_mineral(irxnim(i, 1)) = porchange(i, 2)
        endif
        sarea(i) = xcoef(i)
     else
        write(ierr, *) 'Internal error setting reaction types.'
        goto 667
     endif
  enddo
  ! End of rxn value-setting.
  ! End of FEHM compatibility block.
  goto 6000

  ! Go to 666 for user errors
666 if(reading .eqv. .true.) then
     write(ierr, '(a14, i3, a9)') 'Error on line ', linen, ' of trxn.'
     write(ierr, '(i3, a5, a, a3)') linen, ': >> ', line(1:len_trim(line)), ' <<'
     write(iptty, '(a14, i3, a9)') 'Error on line ', linen, ' of trxn.'
     write(iptty, '(i3, a5, a, a3)') linen, ': >> ', line(1:len_trim(line)), ' <<'
  else
     write(ierr, *) 'Error in trxn in postprocessing.'
     write(iptty, *) 'Error in trxn in postprocessing.'
  endif
  stop

  ! Go to 667 for internal errors
667 write(iptty, *) 'An internal error has occurred in trxn.  This should never happen.'
  write(iptty, *) 'Please contact mschauer@lanl.gov with information about the error.'
  stop

  ! Go to 800 to perform cleanup if no errors have occurred
6000 macroread(5) = .true.
  if(debug) then
     write(iptty, *) 'trxn read successfully.'
     !	call trxn_varcheck
  endif

end subroutine readtracrxn
