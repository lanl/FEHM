module trxnvars ! Basic variables shared among trxn's subroutines

	integer NAME_MAX, SPEC_MAX, LINE_MAX		! Parameters for the maximum space to allocate for different fields.  Probably should get rid of them
	character*4 aSPEC_MAX, aNAME_MAX
	parameter ( NAME_MAX = 40 )
	parameter ( SPEC_MAX = 100 )
	parameter ( LINE_MAX = 200 )
	parameter ( aSPEC_MAX = '100' )	
	parameter ( aNAME_MAX = '40' )
	integer inpttmp					! Holder for input file handle
	integer linen					! Number of the current line
	character(LINE_MAX) line			! The current line being read (scrapped)
	logical debug             			! Whether we are in debug mode
	logical debug_stop				! Whether to halt execution after reading the macro
	integer trxn_flag         			! Whether we are calling tracrxn or trac
	integer numzones, zonemax			! Number of zones and maximum zone number (persistent)
	integer, allocatable :: zonenums(:)		! Numbers corresponding to zones
	character(20), allocatable :: zonenames(:)	! Zone names
	logical molein					! Whether we are using the total mole initial water option
	logical interactive				! Whether we are reading interactively

end module trxnvars
