subroutine tracrxn_init
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
!D1 Initialize variables for trxn macro.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 Initial implementation: 2011,  Programmer: M. Schauer
!D2
!***********************************************************************

  use comchem
  use comdi
  use comai
  use comcouple
  use comrxni

  implicit none

  character*20, dimension(100) :: array, species, array3
  character*20, dimension(100, 100) :: array2
  character*200 line
  character*20 keyword, null
  character state, master
  integer i, j, k

  ! Set tpor_flag true or false
  ! Read or calculate nspeci,ncpnt,nimm,nvap,isorp,and numsorp
  debug = .false.
  tpor_flag = .false.
  nspeci = 0
  ncpnt = 0
  nimm = 0
  nvap = 0
  numrxn = 0
  ncplx = 0
  ngroups = 0
  rxn_flag = 0
  numsorp = 1
  do
1200 read(inpt, '(a200)', end=1800, err=1666) line
     if(index(line, '#') .ne. 0) then
        line = line(1:index(line, '#') - 1)
     endif
     if(len_trim(line) .eq. 0) then
        goto 1200
     endif
     read(line, *, end=1666) keyword
     if(keyword(1:3) .eq. 'end') then
        goto 1800
     elseif(keyword .eq. 'ctrl') then
        array = '*'
        read(line, *, end=1104) keyword, (array(i), i = 1, 100)
1104    do i = 1, 100
           if(array(i) .eq. '*') then
              exit
           elseif(array(i) .eq. 'debug') then
              debug = .true.
              write(iptty, *) 'trxn debug on.'
           elseif(array(i) .eq. 'rxnon') then
              rxn_flag = 1
              if(debug) write(iptty, *) 'Reactions enabled for this simulation.'
           else
              write(ierr, *) 'Error:  keyword "', array(i), '" in ctrl not known.'
              goto 1666
           endif
        enddo

     elseif(keyword .eq. 'sorp') then
        do
           read(inpt, '(a200)', end=1666, err=1666) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           if(line(1:1) .eq. '.') then
              numsorp = numsorp + 1
           endif
        enddo
     elseif(keyword .eq. 'comp') then
        master='*'
        array3 = '*'
        do
           read(inpt, '(a200)', end=1666, err=1666) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           nspeci = nspeci + 1
           read(line, *, end=1102) state, species(nspeci), master
1102       if(state(1:1) .eq. 's') then
              nimm = nimm + 1
           elseif(state(1:1) .eq. 'v') then
              nvap = nvap + 1
           elseif((state(1:1) .eq. 'a') .or. (state(1:1) .eq. 'h')) then
              ncpnt = ncpnt + 1
              array3(ncpnt) = species(nspeci)
           else
              write(ierr, *) 'Unknown phase "', state, '" in comp.'
              goto 1666
           endif
           !				if((master .ne. '*') .and. (master .ne. '')) then
           !					ncplx = ncplx + 1
           !				endif
        enddo
     elseif(keyword .eq. 'assign') then
        array = '*'
        read(line, *, err=1666, end=1101) keyword, (array(i), i = 1, 100)
1101    do i = 1, 100
           if(array(i) .eq. 'tpor') then
              tpor_flag = .true.
              exit
              !				elseif(array(i) .eq. 'rxn') then
              !					rxn_flag = 1
              !					exit
           endif
        enddo
     elseif(keyword .eq. 'rxn') then
        numrxn = numrxn + 1
        do
           read(inpt, '(a200)', end=1666, err=1666) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
        enddo
     elseif(keyword .eq. 'group') then
        array2 = '*'
        do
           read(inpt, '(a200)', end=1666, err=1666) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           ngroups = ngroups + 1
           read(line, *, end=1103) (array2(ngroups, i), i = 1, 100)
1103       continue
        enddo
     elseif(keyword .eq. 'stoich') then
        ncplx = 0
        do
           read(inpt, '(a200)', end=1666, err=1666) line
           if(index(line, '#') .ne. 0) then
              line = line(1:index(line, '#') - 1)
           endif
           if(len_trim(line) .eq. 0) then
              exit
           endif
           ncplx = ncplx + 1
        enddo
     elseif(keyword .eq. 'debug') then
        write(*, *) 'trxn debug on.'
        debug = .true.
     else
        continue
     endif
  enddo

1800 if(rxn_flag .eq. 1) then
     if(irun.eq.1) then
        if (ngroups .ne. 0) then
           allocate(group(ngroups,ncpnt),pos(ngroups,ncpnt))
           allocate(n_couple_species(ngroups),fzero(ngroups))
           group = 0
           pos = 0
           n_couple_species = 0
        end if
     end if
     do i = 1, ngroups
        do j = 1, 100
           if(array2(i, j) .eq. '*') then
              exit
           endif
           do k = 1, ncpnt
              if(array2(i, j) .eq. array3(k)) then
                 group(i, k) = 1
                 exit
              endif
           enddo
           ! TODO Perhaps we should implement a check here to ensure that everything in group is in comp...
        enddo
     enddo
  endif

  if(debug) write(iptty, *) 'trxn preprocessing complete.'
  return

1666 write(ierr, *) 'Error in trxn in preprocessing.'
  write(iptty, *) 'Error in trxn in preprocessing.'
  stop

end subroutine tracrxn_init
