subroutine trxn_varcheck
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
!D1 Check input for trxn macro.
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

  integer i

  !if(debug .eqv. .false.) then
  !	return
  !endif

  if(rxn_flag .eq. 1) then
     write(*, *)
     write(*, *) 'RXN DUMP:'
     write(*, *) ncpnt, ' aqueous components:'
     write(*, *) (i, ':', cpntnam(i), '	', i = 1, ncpnt)
     write(*, *)
     write(*, *) nimm, ' immobile components:'
     write(*, *) (i, ':', immnam(i), '	', i = 1, nimm)
     write(*, *)
     write(*, *) nvap, ' vapor components:'
     write(*, *) (i, ':', vapnam(i), '	', i = 1, nvap)
     write(*, *)
     write(*, *) ncplx, ' aqueous complexes:'
     write(*, *) (i, ':', cplxnam(i), '	', i = 101, ncplx + 100)
  endif

  write(*, *) 'EXTRA INFORMATION:'
  write(*, *) numsorp, ldsp, dispsame, strac_max
  write(*, *) hvliquid, hvvapor
  write(*, *) tclx
  write(*, *) tcly
  write(*, *) tclz
  write(*, *) sehdiff
  write(*, *) diffmfl
  write(*, *) iadsfl
  write(*, *) a1adfl
  write(*, *) a2adfl
  write(*, *) betadfl
  write(*, *) qcout
  write(*, *) qcin
  write(*, *) qcrxn
  write(*, *) hawwa
  write(*, *) pimm
  write(*, *) icns
  write(*, *) '*****'
  write(*, *) ncplx, numrxn, ngroups, iskip, rsdmax
  write(*, *) ifxconc
  write(*, *) cpntgs
  write(*, *) group
  write(*, *) '1', neg_conc_possible
  write(*, *) '2', temp_model_kin
  write(*, *) '3', simmmx
  write(*, *) ckeq
  write(*, *) heq
  write(*, *) idrxn
  write(*, *) irxnic
  write(*, *) irxnim
  write(*, *) kd
  write(*, *) pdstic
  write(*, *) pdstim
  write(*, *) ckeqlb
  write(*, *) ckmtrn
  write(*, *) naqsp
  write(*, *) nimsp
  write(*, *) spstoic
  write(*, *) tcoeff
  write(*, *) sticirrv
  write(*, *) stimirrv
  write(*, *) cplx
  write(*, *) cpnt

  write(*, *) 'Press enter to continue.'
  read(*, *)

end subroutine trxn_varcheck
