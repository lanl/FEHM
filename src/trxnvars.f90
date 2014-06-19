module trxnvars ! Basic variables shared among trxn's subroutines
!***********************************************************************
! Copyright 2012. Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los
! Alamos National Laboratory (LANL), which is operated by Los Alamos
! National Security, LLC for the U.S. Department of Energy. The U.S.
! Government has rights to use, reproduce, and distribute this software.
! Neither the U.S. Government nor Los Alamos National Security, LLC or
! persons acting on their behalf, make any warranty, express or implied,
! or assumes any liability for the accuracy, completeness, or usefulness
! of the software, any information pertaining to the software, or
! represents that its use would not infringe privately owned rights.
!
! The software being licensed may be Export Controlled.   It may not be
! distributed or used by individuals or entities prohibited from having
! access to the software package, pursuant to United States export
! control laws and regulations.  An export control review and
! determination must be completed before LANS will provide access to the
! identified Software.
!***********************************************************************

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
