module lookupvars ! Variables for database lookup subroutine in trxn

	use trxnvars

	logical flookup							! Whether lookup has been called
	character(LINE_MAX), dimension(SPEC_MAX) :: options		! Options
	character(LINE_MAX) database					! Database file name
	integer dat							! Filehandle for database
	integer nlocomp, nlocompaq					! Number of components lookup needs to add
	character(NAME_MAX), dimension(SPEC_MAX, 3) :: locomp		! Components to be added by lookup (1 = state, 2 = component, 3 = master)
	integer ndatcomp						! Number of complexes in the database
	character(NAME_MAX), allocatable :: datcomp(:), datmaster(:)	! Complexes in the database and their master species
	integer ndatcplx, ndatmin, maxdatcplx, maxdatmin		! Number of complexes and minerals in the database, and the maximum number of reactants and products for each
	real*8, allocatable :: datcden(:)				! Molecular weights in the database
	character(NAME_MAX), allocatable :: datcplxmain(:)		! Complex formulae
	character(NAME_MAX), allocatable :: datcplx(:, :)		! Species for complexes in the database.  datcplx(i, j) is reaction i, species j
	real*8, allocatable :: datcplxstoic(:, :)			! Stoichiometry for complexes in the database.  Same indexing as datcplx
	real*8, allocatable :: datcplxlkeq(:)				! log_k for complexes in the database
	real*8, allocatable :: datcplxheq(:)				! delta_h for complexes in the database
	real*8, allocatable :: datcplxtemp(:, :)			! Temperature dependency coefficients for complexes in the database
	character(NAME_MAX), allocatable :: datminnam(:)		! Names of minerals in the database
	character(NAME_MAX), allocatable :: datminmain(:)		! Mineral formulae
	character(NAME_MAX), allocatable :: datmin(:, :)		! Species of minerals in the database.  Same indexing as datcplx
	real*8, allocatable :: datminstoic(:, :)			! Stoichiometry for minerals in the database. Same indexing as datcplx
	real*8, allocatable :: datminlkeq(:)				! log_k for minerals in the database
	real*8, allocatable :: datminheq(:)				! delta_h for minerals in the database
	real*8, allocatable :: datmintemp(:, :)				! Temperature dependency coefficients for minerals in the database

end module lookupvars
