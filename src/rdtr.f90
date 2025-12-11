!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC   All rights reserved
! Unless otherwise indicated,  this information  has been authored by an
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
!D1 Read and process input for trxn macro.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 Initial implementation: August 2011,  Programmer: Matthew Schauer
!D2 
!D2 August 2012  M. Schauer
!D2 Major bug fixes and added functionality
!D2
!***********************************************************************

subroutine readline ! Reads a line of input, strips off comments and whitespace, and updates the line number.

	use comai
	use trxnvars
	use comdi, only:nspeci

	implicit none

	read(inpt, '(a200)', end=2101, err=2102) line
!	write(*,*) 'inreadline ', line
	if (index(line, '#') .ne. 0) then
		line = line(1:index(line, '#') - 1)
	endif
	line = trim(adjustl(line))
	do
		if (line(1:1) .eq. '	') then
			line = line(2:len_trim(line))
		else
			exit
		endif
	enddo
	goto 3000

2101	write(ierr, '(a)') 'Error:  Unexpected end of input.'
	goto 2666

2102	write(ierr, '(a, i3, a)') 'Error reading input on line ', linen + 1
2666	write(ierr, '(a)') 'Error in trxn in readline.'
	write(iptty, '(a)') 'Error in trxn in readline.'
	stop
3000	linen = linen + 1
	return

end subroutine readline

subroutine trxninit ! Called by scanin.f to preset some basic variables.

      	use comchem
      	use comdi
      	use comai
	use comcouple
	use comrxni
	use comco2
!	use lookupvars
	use trxnvars

	implicit none

	character*20, dimension(100) :: array, array3,  masters, states, groupspecs
	character*20, dimension(100, 100) :: array2
	integer, dimension(100) :: iarray
	character*20 keyword, null, e, f
	character*200 word, libfname
	character state
	integer blocknum, i, j, k, o, p, libread, libnum
	logical flag, fctrl

	! trxninit's jobs:
	! Set tpor_flag true or false
	! Read or calculate nspeci, ncpnt, nimm, nvap, isorp, numd, and numsorp
	debug = .false.
	debug_stop = .false.
	tpor_flag = .false.
	strong=.false.
	blocknum = 0
	flookup = .false.
	fctrl = .false.
	libread = 0
	nspeci = 0
	ncpnt = 0
	nimm = 0
	nvap = 0
	numrxn = 0
	ncplx = 0
	ngroups = 0
	rxn_flag = 0
	numsorp = 0
	numd = 0
	nlocomp = 0
	nlocompaq = 0
	do
1200		call readline
		if (len_trim(line) .eq. 0) goto 1200
		read(line, *, end=1666, err=1600) keyword
		if (keyword(1:3) .eq. 'end') then
			if (libread .gt. 0) goto 1801
			goto 1800
		elseif (keyword .eq. 'include') then
			blocknum = blocknum + 1
			libfname = '*'
			if (index(line, ' ') .eq. 0) then
				write(ierr, '(a)') 'Error:  Missing include filename or bad input format in include.  '// &
					'Please note that TAB cannot be used to separate the keyword from the filename.'
				goto 1666
			endif
			libfname = line(index(line, ' ') + 1:len_trim(line))
			if (debug) write(iptty, '(a)') 'Reading from library "'//trim(libfname)//'".'
			libread = libread + 1
			libnum = 23 + libread
			open(libnum, file=libfname, err=1109)
			if (libread .eq. 1) inpttmp = inpt
			inpt = libnum
			rewind(libnum)
			call readline
			do
				if (len_trim(line) .ne. 0) exit
				call readline
			enddo
			read(line, *, end=340, err=1600) e, f, null
340			if ((e .ne. 'trxn') .or. (f(1:3) .ne. 'lib')) then
				write(ierr, '(a)') 'Error:  "'//trim(libfname)//'" is not a valid trxn library.'
				goto 1666
			endif
			goto 1200
1801			if (libread .eq. 1) then
				inpt = inpttmp
			else
				inpt = inpt - 1
			endif
			close(libnum)
			if (debug) write(iptty, '(a)') 'Library read done.'
			libread = libread - 1
			goto 1200
1109			write(ierr, '(a)') 'Error opening include file "'//trim(libfname)//'".'
		elseif (keyword .eq. 'ctrl') then
			if (debug) write(iptty, '(a)') 'Prereading ctrl.'
			blocknum = blocknum + 1
			if (blocknum .gt. 1) then
				write(ierr, '(a)') 'Error:  If a ctrl block is provided, '// &
					'it must be the first block in trxn.'
				goto 1666
			endif
			fctrl = .true.
			array = '*'
			read(line, *, end=1104, err=1600) keyword, (array(i), i = 1, 100)
1104			do i = 1, 100
				if (array(i) .eq. '*') then
					exit
				elseif (array(i) .eq. 'debug') then
					debug = .true.
					write(iptty, '(a)') 'trxn debug on.'
				elseif (array(i) .eq. 'rxnon') then
					rxn_flag = 1
					if (debug) write(iptty, '(a)') 'Reactions enabled for this simulation.'
				elseif (array(i) .eq. 'stop') then
					debug_stop = .true.
				elseif(array(i).eq.'strong') then
					if(icarb.eq.1) then
					strong=.true.  ! want to modify this to search for calcite
					endif
				elseif (array(i) .eq. 'co2_couple') then
					if (icarb .eq. 1) then
						if (debug) write(iptty, '(a)') &
							'CO2 coupling enabled for this simulation.'
						co2_couple = 1
					else
						write(ierr, '(a)') 'Warning:  Attempted to enable CO2 coupling '// &
							'without a carb macro.  CO2 coupling will remain disabled.'
					endif
				else
					write(ierr, '(a)') 'Error:  keyword "'//trim(array(i))// &
						'" in ctrl not recognized.'
					goto 1666
				endif
			enddo
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
		elseif (keyword .eq. 'sorp') then
			if (debug) write(iptty, '(a)') 'Prereading sorp.'
			blocknum = blocknum + 1
			do
				call readline
				if (len_trim(line) .eq. 0) exit
				if (line(1:1) .eq. '.') then
					numsorp = numsorp + 1
				endif
			enddo
		elseif (keyword .eq. 'disp') then
			if (debug) write(iptty, '(a)') 'Prereading disp.'
			blocknum = blocknum + 1
			do
				call readline
				if (len_trim(line) .eq. 0) exit
				numd = numd + 1
			enddo
		elseif (keyword .eq. 'lookup') then
			if (debug) write(iptty, '(a)') 'Prereading lookup.'
			blocknum = blocknum + 1
			if (((blocknum .gt. 1) .and. (.not. fctrl)) .or. ((blocknum .gt. 2) .and. (fctrl))) then
				write(ierr, '(a)') 'Error:  If a lookup block is provided, it must be the '// &
					'first block in trxn (excepting the ctrl block, if provided).'
				goto 1666
			endif
			if (debug) write(iptty, '(a)') 'Reading database information...'
			if (flookup) then
				write(ierr, '(a)') 'Warning:  lookup has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(datcomp)
				deallocate(datmaster)
				deallocate(datcden)
				deallocate(datcplxmain)
				deallocate(datcplx)
				deallocate(datcplxstoic)
				deallocate(datcplxlkeq)
				deallocate(datcplxheq)
				deallocate(datcplxtemp)
				deallocate(datminnam)
				deallocate(datminmain)
				deallocate(datmin)
				deallocate(datminstoic)
				deallocate(datminlkeq)
				deallocate(datminheq)
				deallocate(datmintemp)
			endif
			flookup = .true.
			options = '*'
			database = '*'
			ndatcomp = 0
			ndatcplx = 0
			ndatmin = 0
			maxdatcplx = 0
			maxdatmin = 0
			read(line, *, end=282, err=283) keyword, (options(i), i = 1, SPEC_MAX)
282			do i = 1, SPEC_MAX - 1
				if (options(i) .eq. '*') then
					exit
				elseif (options(i)(1:5) .eq. 'inter') then
					interactive = .true.
				else
					write(ierr, '(a)') 'Error:  Unrecognized option "'//trim(options(i))// &
						'" in lookup.'
					goto 1666
				endif
			enddo
			if (index(line, ' ') .eq. 0) then
				write(ierr, '(a)') 'Missing database filename or bad input format in lookup.  '// &
				'Please note that TAB cannot be used to separate the keyword from the database name.'
				goto 1666
			endif
			database = trim(adjustl(line(index(line, ' ') + 1: len_trim(line))))
			if (debug) write(iptty, '(a)', advance='no')  'Reading from database file "'// &
				trim(database)//'"...  '
			open(dat, file=database, status='old', err=281)
			do
300				read(dat, '(a200)', end=293, err=292) line ! Preread to determine how to allocate arrays
				if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
				line = trim(line)
				if (len_trim(line) .eq. 0) then
					continue
				elseif (line(1:1) .eq. '%') then
					if (line .eq. '% MASTER') then
						do
325							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 300
							elseif (len_trim(line) .eq. 0) then
								goto 325
							endif
							ndatcomp = ndatcomp + 1
						enddo
					elseif (line .eq. '% CPLX') then
						do
301							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 300
							elseif (len_trim(line) .eq. 0) then
								goto 301
							endif
							ndatcplx = ndatcplx + 1
							backspace dat
							array = '*'
							array3 = '*'
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=295, err=292) (array(i), array3(i), i = 1, SPEC_MAX) ! Reactants
295							do i = 1, SPEC_MAX
								if (array3(i) .eq. '*') then
									maxdatcplx = max(maxdatcplx, i)
									exit
								endif
							enddo
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=296, err=292) (array(i), array3(i), i = 1, SPEC_MAX) ! Products
296							do i = 1, SPEC_MAX
								if (array3(i) .eq. '*') then
									maxdatcplx = max(maxdatcplx, i)
									exit
								endif
							enddo
							read(dat, *, end=292, err=292) ! log_k
							read(dat, *, end=292, err=292) ! delta_h
							read(dat, *, end=292, err=292) ! temp
						enddo
					elseif (line .eq. '% IMM') then
						do
302							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 300
							elseif (len_trim(line) .eq. 0) then
								goto 302
							endif
							ndatmin = ndatmin + 1
							backspace dat
							array = '*'
							array3 ='*'
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=297, err=292) (array(i), array3(i), i = 1, SPEC_MAX) ! Reactants
297							do i = 1, SPEC_MAX
								if (array3(i) .eq. '*') then
									maxdatmin = max(maxdatmin, i)
									exit
								endif
							enddo
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=298, err=292) (array(i), array3(i), i = 1, SPEC_MAX) ! Products
298							do i = 1, SPEC_MAX
								if (array3(i) .eq. '*') then
									maxdatmin = max(maxdatmin, i)
									exit
								endif
							enddo
							read(dat, *, end=292, err=292) ! log_k
							read(dat, *, end=292, err=292) ! delta_h
							read(dat, *, end=292, err=292) ! temp
						enddo
					elseif (line .eq. '% END') then
						exit
					else
						write(ierr, '(a)') 'Error:  Unrecognized section header "'// &
							trim(line)//'" in database.'
						goto 1666
					endif
				else
					write(ierr, '(a)') 'Error:  Unrecognized token in database.'
					goto 1666
				endif		
			enddo
293			allocate(datcomp(ndatcomp))
			allocate(datmaster(ndatcomp))
			allocate(datcden(ndatcomp))
			allocate(datcplxmain(ndatcplx))
			allocate(datcplx(ndatcplx, maxdatcplx))
			allocate(datcplxstoic(ndatcplx, maxdatcplx))
			allocate(datcplxlkeq(ndatcplx))
			allocate(datcplxheq(ndatcplx))
			allocate(datcplxtemp(ndatcplx, 5))
			allocate(datminnam(ndatmin))
			allocate(datminmain(ndatmin))
			allocate(datmin(ndatmin, maxdatmin))
			allocate(datminstoic(ndatmin, maxdatmin))
			allocate(datminlkeq(ndatmin))
			allocate(datminheq(ndatmin))
			allocate(datmintemp(ndatmin, 5))
			datcplx = '*'
			datmin = '*'
			datcden = 0
			datcplxtemp = 0
			datmintemp = 0
			j = 0
			k = 0
			o = 0
			rewind dat
			do
303				read(dat, '(a200)', end=324, err=292) line ! The real read
				if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
				line = trim(line)
				if (len_trim(line) .eq. 0) then
					continue
				elseif (line(1:1) .eq. '%') then
					if (line .eq. '% MASTER') then
						do
326							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 303
							elseif (len_trim(line) .eq. 0) then
								goto 326
							endif
							o = o + 1
							read(line, *, end=292, err=292) datcomp(o), datmaster(o), datcden(o)
						enddo
					elseif (line .eq. '% CPLX') then
						do
304							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 303
							elseif (len_trim(line) .eq. 0) then
								goto 304
							endif
							backspace dat
							j = j + 1
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=306, err=292) datcplxmain(j) !(datcplxstoic(j, 1, i), datcplx(j, 1, i), i = 1, maxdatcplx) ! Reactants
306							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=307, err=292) &
								(datcplxstoic(j, i), datcplx(j, i), i = 1, maxdatcplx) ! Products
307							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=310, err=292) datcplxlkeq(j) ! log_k
							goto 311
310							datcplxlkeq(j) = 0 ! log_k not specified
311							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=312, err=292) datcplxheq(j) ! delta_h
							goto 313
312							datcplxheq(j) = 0 ! delta_h not specified
313							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=314, err=292) (datcplxtemp(j, i), i = 1, 5) ! temp
							goto 315
314							datcplxtemp(j, 1:5) = 0 ! temp not (completely) specified
315							continue
						enddo
					elseif (line .eq. '% IMM') then
						do
305							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							line = trim(line)
							if (line(1:1) .eq. '%') then
								backspace dat
								goto 303
							elseif (len_trim(line) .eq. 0) then
								goto 305
							endif
							backspace dat
							k = k + 1
							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=299, err=292) datminnam(k), datminmain(k)
308							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=309, err=292) &
								(datminstoic(k, i), datmin(k, i), i = 1, maxdatmin) ! Products
309							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=316, err=292) datminlkeq(k) ! log_k
							goto 317
316							datminlkeq(k) = 0 ! log_k not specified
317							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=318, err=292) datminheq(k) ! delta_h
							goto 319
318							datminheq(k) = 0 ! delta_h not specified
319							read(dat, '(a200)', end=294, err=292) line
							if (index(line, '#') .ne. 0) line = line(1:index(line, '#') - 1)
							read(line, *, end=320, err=292) (datmintemp(k, i), i = 1, 5) ! temp
							goto 321
320							datmintemp(k, 1:5) = 0 ! temp not (completely) specified
321							continue
						enddo
					elseif (line .eq. '% END') then
						exit
					else
						write(ierr, '(a)') 'Internal error reading headers in database.'
						goto 1667
					endif
				else
					write(ierr, '(a)') 'Internal error parsing tokens in database.'
					goto 1667
				endif		
			enddo
324			close(dat)
			if (debug) write(iptty, '(i3, a, i3, a)') ndatcplx, ' aqueous complexation reactions and ', &
				ndatmin, ' mineral dissolution/precipitation reactions read.'
			goto 1200
280			write(ierr, '(a)') 'Error reading path to lookup database.'
			goto 1666
281			write(ierr, '(a)') 'Error:  Could not open lookup database.'
			goto 1666
283			write(ierr, '(a)') 'Error reading options in lookup.'
			goto 1666
292			write(ierr, '(a)') 'Syntax error reading database in lookup.'
			goto 1666
294			write(ierr, '(a)') 'Error:  Premature end of database file in lookup.'
			goto 1666
299			write(ierr, '(a)') 'Internal error reading database in lookup.'
			goto 1667
		elseif (keyword .eq. 'comp') then
			if (debug) write(iptty, '(a)') 'Prereading comp.'
			blocknum = blocknum + 1
			masters = '*'
			states = '*'
			array3 = '*'
			do
1105			call readline
				if(len_trim(line).eq.0) exit
				nspeci=nspeci+1
				goto 1105
			end do	
			do i=1,nspeci+1
					backspace inpt
				end do
                                if(.not.allocated(species)) allocate (species(nspeci))
				nspeci=0
			do

 1117 				call readline
				if (len_trim(line) .eq. 0) exit
				nspeci=nspeci+1
				read(line, *, end=1102, err=1600) states(nspeci), species(nspeci), masters(nspeci)
1102				do i = 1, nlocomp
					if (locomp(i, 2) .eq. species(nspeci)) then
						species(nspeci) = '*'
						nspeci = nspeci - 1
						goto 1117
					endif	
				enddo
				if (states(nspeci)(1:1) .eq. 's') then
					nimm = nimm + 1
				elseif (states(nspeci)(1:1) .eq. 'g') then
					nvap = nvap + 1
				elseif ((states(nspeci)(1:1) .eq. 'a') .or. (states(nspeci)(1:1) .eq. 'h')) then
					ncpnt = ncpnt + 1
					array3(ncpnt) = species(nspeci)
				else
					write(ierr, '(a)') 'Error:  unknown phase "'//trim(states(i))//'" in comp.'
					goto 1666
				endif
				!if ((master .ne. '*') .and. (master .ne. '')) then
				!	ncplx = ncplx + 1
				!endif
			enddo
	
		elseif (keyword .eq. 'assign') then
			if (debug) write(iptty, '(a)') 'Prereading assign.'
			blocknum = blocknum + 1
			array = '*'
			read(line, *, err=1666, end=1101) keyword, (array(i), i = 1, 100)
1101			do i = 1, 100
				if (array(i) .eq. 'tpor') then
					tpor_flag = .true.
					exit

				!elseif (array(i) .eq. 'rxn') then
				!	rxn_flag = 1
				!	exit
				endif
			enddo
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
		elseif (keyword .eq. 'rxn') then
			if (debug) write(iptty, '(a, i3, a)') 'Prereading rxn ', numrxn + 1, '.'
			blocknum = blocknum + 1
			numrxn = numrxn + 1
			call readline
			e = '*'
			f = '*'
			read(line, *, end=289, err=1600) e, f
289			if (f .eq. '*') then
				if (.not. flookup) then
						write(ierr, '(a)') 'Error:  A lookup block must be provided before '// &
							'rxn if dynamic reaction lookup is to be used.'
						goto 1666
					endif
					do i = 1, ndatmin
						if ((datminmain(i) .eq. e) .or. (datminnam(i) .eq. e)) then
							do k = 1, nspeci
								if (datminmain(i) .eq. species(k)) goto 327
							enddo
							nlocomp = nlocomp + 1
							nspeci = nspeci + 1
							nimm = nimm + 1
							locomp(nlocomp, 1) = 's'
							locomp(nlocomp, 2) = datminmain(i)
							locomp(nlocomp, 3) = '*'
327							do j = 1, maxdatmin
								if (datmin(i, j) .eq. '*') exit
								do k = 1, nlocomp
									if (locomp(k, 3) .eq. datmin(i, j)) goto 290
								enddo
								do k = 1, nspeci
									if (masters(k) .eq. datmin(i, j)) goto 290
								enddo
								nlocomp = nlocomp + 1
								nlocompaq = nlocompaq + 1
								nspeci = nspeci + 1
								ncpnt = ncpnt + 1
								locomp(nlocomp, 1) = 'a'
								locomp(nlocomp, 3) = datmin(i, j)
								do k = 1, ndatcomp
									if (datmaster(k) .eq. datmin(i, j)) then
										locomp(nlocomp, 2) = datcomp(k)
										array3(ncpnt) = datcomp(k)
										goto 290
									endif
								enddo
								write(ierr, '(a)') 'Error:  Master species "'// &
									trim(datcplx(i, j))//'" in mineral included '// &
									'from lookup database not found in MASTER '// &
									'section of lookup database.'
								goto 1666
290								continue
							enddo
						goto 1200
						endif
					enddo
					write(ierr, '(a)') 'Error:  Complex "'//trim(e)// &
						'" in cplx not found in lookup database.'
					goto 1666
				endif
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
		elseif (keyword .eq. 'group') then
			if (debug) write(iptty, '(a)') 'Prereading group.'
			blocknum = blocknum + 1
			array2 = '*'
			groupspecs = '*'
			j = 0
			do
				call readline
				if (len_trim(line) .eq. 0) exit
				ngroups = ngroups + 1
				read(line, *, end=1103, err=1600) (array2(ngroups, i), i = 1, 100)
1103				do i = 1, 100
					if (array2(ngroups, i) .eq. '*') exit
					j = j + 1
					groupspecs(j) = array2(ngroups, i)
				enddo
			enddo
			do i = 1, 100
				if (groupspecs(i) .eq. '*') exit
				do j = 1, i - 1
					if (groupspecs(j) .eq. groupspecs(i)) then
						write(ierr, '(a)') 'Error:  Multiple definition of component "'// &
							trim(groupspecs(i))//'" in group.'
						goto 1666
					endif
				enddo
			enddo
		elseif (keyword .eq. 'cplx') then
			if (debug) write(iptty, '(a)') 'Prereading cplx.'
			blocknum = blocknum + 1
			ncplx = 0
			do
291				call readline
				if (len_trim(line) .eq. 0) exit
				ncplx = ncplx + 1
				e = '*'
				f = '*'
				read(line, *, end=288, err=1600) e, f
288				if (f .ne. '=') then
					if (.not. flookup) then
						write(ierr, '(a)') 'Error:  A lookup block must be provided before '// &
							'cplx if dynamic complex lookup is to be used.'
						goto 1666
					endif
					do i = 1, ndatcplx
						if (datcplxmain(i) .eq. e) then
							do j = 1, maxdatcplx
								if (datcplx(i, j) .eq. '*') exit
								do k = 1, nlocomp
									if (locomp(k, 3) .eq. datcplx(i, j)) goto 287
								enddo
								do k = 1, nspeci
									if (datcplx(i, j) .eq. masters(k)) goto 287
								enddo
								nlocomp = nlocomp + 1
								nlocompaq = nlocompaq + 1
								nspeci = nspeci + 1
								ncpnt = ncpnt + 1
								locomp(nlocomp, 1) = 'a'
								locomp(nlocomp, 3) = datcplx(i, j)
								do k = 1, ndatcomp
									if (datmaster(k) .eq. datcplx(i, j)) then
										locomp(nlocomp, 2) = datcomp(k)
										array3(ncpnt) = datcomp(k)
										goto 287
									endif
								enddo
								write(ierr, '(a)') 'Error:  Master species "'// &
									trim(datcplx(i, j))//'" in complex included '// &
									'from lookup database not found in MASTER '// &
									'section of lookup database.'
								goto 1666
287								continue
							enddo
						goto 291
						endif
					enddo
					write(ierr, '(a)') 'Error:  Complex "'//trim(e)// &
						'" in cplx not found in lookup database.'
					goto 1666
				endif
			enddo
		elseif (keyword .eq. 'null') then
			read(line, *, end=1110, err=1600) keyword, e
1110			if (debug) write(iptty, '(a)') 'Skipping null statement "'//trim(e)//'".'
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
		else
			blocknum = blocknum + 1
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
		endif
	enddo

	write(ierr, '(a)') 'Warning:  Missing "end" statement may cause unexpected behavior.'
1800	if (debug) write(iptty, '(a)') 'Done prereading.'
	numsorp = numsorp + 1 ! Leave room for a default model
	numd = numd + 1
	if (nspeci .gt. 40) then
		write(ierr, '(a)') 'Error:  Due to arbitrary limits in FEHM''s solver code, no more than 40 '// &
			'species may be specified.'
		goto 1666
	endif
	if (rxn_flag .eq. 1) then
		if (ngroups .gt. 0) then
			do i = 1, nspeci
				if ((states(i)(1:1) .ne. 'a') .and. (states(i)(1:1) .ne. 'h')) goto 1108
				do j = 1, 100
					if (groupspecs(j) .eq. '*') exit
					if (species(i) .eq. groupspecs(j)) goto 1108
				enddo
				write(ierr, '(a)') 'Error:  Aqueous or Henry''s Law component "'// &
					trim(species(i))//'" in comp not found in group.'
				goto 1666
1108				continue
			enddo
			o = ngroups
			if (flookup) then
				ngroups = nlocomp
				do i = 1, o
					k = 0
					do j = 1, ncpnt
						if (array2(i, j) .eq. '*') exit
						k = k + 1
					enddo
					k = k - 1
					if (k .gt. 0) ngroups = ngroups - k
				enddo
			endif
			allocate(group(ngroups,ncpnt),pos(ngroups,ncpnt))
			allocate(n_couple_species(ngroups),fzero(ngroups))
			group = 0
			pos = 0
			n_couple_species = 0
			do i = 1, o
				do j = 1, 100
					if (array2(i, j) .eq. '*') then
						exit
					endif
					flag = .false.
					do k = 1, ncpnt
						if (array2(i, j) .eq. array3(k)) then
							group(i, k) = 1
							flag = .true.
							exit
						endif
					enddo
					if (.not. flag) then
						write(ierr, '(a)') 'Error:  Aqueous or Henry''s Law component "'// &
							trim(array2(i, j))//'" in group not found in comp.'
						goto 1666
					endif
				enddo
			enddo
			if (flookup) then
				i = o + 1
				do
					if (i .gt.ngroups) exit
					do j = 1, ncpnt
						do k = 1, i - 1
							if (group(k, j) .eq. 1) goto 1106
						enddo
						group(i, j) = 1
						goto 1107
1106						continue
					enddo
1107					i = i + 1
				enddo
			endif
		else
			ngroups = ncpnt
			allocate(group(ngroups,ncpnt),pos(ngroups,ncpnt))
			allocate(n_couple_species(ngroups),fzero(ngroups))
			group = 0
			do i = 1, ncpnt
				group(i, i) = 1
			enddo
			pos = 0
			n_couple_species = 0
		endif
	else
		ngroups = 0
	endif

	if (debug) write(iptty, '(a)') 'trxn preprocessing complete.'
	return

1600	write(ierr, '(a)') 'Error:  Generic read error in trxn preprocessing.'
	goto 1666
1666	write(ierr, '(a)') 'Error in trxn in preprocessing.'
	write(iptty, '(a)') 'Error in trxn in preprocessing.'
	stop
1667	write(iptty, '(a)') 'An internal error has occurred in trxn.  This should never happen, '// &
		'and indicates a bug in the program.  Please contact fehm-dev@lanl.gov with information about the error.'
	stop

	end subroutine

subroutine rdtr ! The main input reader

	use comai
      	use combi
	use comci
      	use comchem
      	use comcouple
	use comdi
      	use comdti
      	use comflow
      	use compart
      	use comrxni
!	use lookupvars
	use trxnvars

	implicit none

	character*10 keyword						! Keyword of the current block being read
	character(LINE_MAX) line2, null, e, f, g, answer		! Generic variables
	integer i, j, k, o, a, c, d					! Generic variables
	real*8 r, u							! Generic variables
	logical flag, flag2, reading					! Generic flags, whether we are reading input
	logical fph							! Whether [H+] is specified as pH
	logical fcomp, fwater, frock, fsorp, fdisp, fvap, fheader, &	! 
		fassign, fhparam, fdiff, fspec, fgroup, fprint, fsol, &	! 
		fdist, fequi, fcplx, fmoles, fuserc			! Flags denoting whether each section has been specified.
	integer libread, libnum						! Whether we are currently reading from an included file (number of levels deep), filehandle for library file
	character(LINE_MAX) :: libfname					! Filename of current library
	integer ncdenspecs						! Number of species according to cden
	character(NAME_MAX), allocatable :: cdenspecs(:)		! Species names according to cden
	real*8, allocatable :: molweights(:)				! Molecular weights given in cden
	character(NAME_MAX), dimension(SPEC_MAX) :: array, array2	! Generic, temporary
	character(NAME_MAX), dimension(SPEC_MAX, SPEC_MAX) :: array3	! Generic, temporary
	integer, dimension(SPEC_MAX) :: narray				! Generic, temporary
	real*8, dimension(100) :: rarray1, rarray2			! Generic, temporary
	character(NAME_MAX), allocatable :: wtnames(:)			! Information on water types will be held internally until the water types are applied to zones.
	character(NAME_MAX), allocatable :: wtspecies(:)		! List of species according to water
	real*8, allocatable :: wtgrid(:, :)				! wtspecies go across the top, and wtnames go down the left side.  wtgrid(i, j) is water type i and species j.
	integer nwt, nwtspecies						! Number of water types and water type species, respectively
	character(NAME_MAX), allocatable :: rtnames(:), rtspecies(:)	! Solid species names and component species.
	real*8, allocatable :: rtgrid(:, :)				! rtspecies go across the top, and rtnames go down the left side.  rtgrid(i, j) is rock type i and species j.
	integer nrt, nrtspecies						! Number of rock types and rock species, respectively.
	character(NAME_MAX), allocatable :: vtnames(:), vtspecies(:)	! Vapor species names and component species.
	real*8, allocatable :: vtgrid(:, :)				! vtgrid(i, j) is vapor type i and species j.
	integer nvt, nvtspecies						! Number of vapor types and vapor species.
	character(NAME_MAX), allocatable :: masters(:)			! Master species for each component in species(:)
	real*8, allocatable :: guesses(:)				! Concentration guesses for each component in species(:)
	character, allocatable :: states(:)				! The state of each species
	integer nspecies, nmasters, naq, nso, nh, nv			! Number of all, aqueous, solid, henry, and vapor species
	integer anspecies						! Number of species according to sorp
	character(NAME_MAX), allocatable :: sorpnames(:), sorpspecs(:)	! Names of adsorption models and species
	integer, allocatable :: lsorptypes(:, :), vsorptypes(:, :)	! Types of adsorption models
	logical fltype, fvtype, fa1l, fa2l, fbl, fa1v, fa2v, fbv	! Flags for whether each column has been specified
	real*8, allocatable :: a1l(:, :), a2l(:, :), bl(:, :), &	! 
		a1v(:, :), a2v(:, :), bv(:, :)				! Parameters used for sorp
	integer nsorp							! Number of adsorption models
	integer ndisp							! Number of dispersivity models
	character(NAME_MAX), allocatable :: dispnames(:)		! Names of dispersivity types
	character(NAME_MAX), dimension(SPEC_MAX) :: dispparams		! Dispersivity parameters
	integer dispmode						! 0 for xyz, 1 for ldsp
	real*8, allocatable :: lx(:), ly(:), lz(:), ll(:), lt(:), &
		vx(:), vy(:), vz(:), vl(:), vt(:)			! X, Y, Z, longitudinal, and transverse dispersivity for liquids and vapors
	real*8 ldiff, vdiff						! Diffusion coefficients for liquid and vapor
	integer ldiffm, vdiffm						! Diffusion models for liquid and vapor
	integer, allocatable :: hmodels(:)				! Types of Henry models
	character(NAME_MAX), allocatable :: hhtspecies(:)		! Henry's Law species in hparam
	integer hnhtspecies						! Number of Henry's Law species according to hparam
	real, allocatable :: ah(:), ah1(:), ah2(:), ah3(:), ah4(:), &	! 
		ah5(:), dhh(:), hh(:)					! Parameters for Henry models
!	integer gngroups, gnspecies					! Number of groups and species according to group
!	character(NAME_MAX), allocatable :: gspecies(:)			! List of species according to group
!	integer, allocatable :: ggroups(:)				! List of groups
!	integer, allocatable :: ngroupdefs(:)				! The number of times each component is specified in group.  Should be 1 for every component.
	integer nprint							! Number of printable species
	integer nstoichcomponents, nstoichspecies			! Number of components and species in stoich
	integer ncomplexes						! Number of complexes
	character(NAME_MAX), allocatable :: complexes(:)		! List of complex names
	character(NAME_MAX), allocatable :: complexcontents(:, :)	! List of master species that make up each complex
	real*8, allocatable :: complexstoich(:, :)			! Stoichiometry of the master species in each complex
	integer ndistmodels						! Number of distribution models
	character(NAME_MAX), allocatable :: distmodelnames(:)		! Names of distribution models
	real*8, allocatable :: distmodels(:, :, :)			! distmodels(i, j) is model i, datapoint j, temperature at k = 1 and distribution coefficient at k = 2
	integer nsolmodels						! Number of solubility models
	character(NAME_MAX), allocatable :: solmodelnames(:)		! Names of solubility models
	real*8, allocatable :: solmodels(:, :, :)			! solmodels(i, j) is model i, datapoint j, temperature at k = 1 and solubility constant at k = 2
	real*8, allocatable :: equconstants(:), enthalpies(:)		! Equilibrium constants and heats of formation for aqueous complexes
	logical logten, cequi						! Whether cplx values are log 10 values, whether we are using equi for our constants
	integer necomplexes, neqtemps					! The number of species according to equi, and the maximum number of equilibrium temperatures in equi
	character(NAME_MAX), allocatable :: ecomplexes(:)		! Species according to equi
	character(NAME_MAX), allocatable :: printspecies(:)		! List of species to print

	real*8, allocatable :: eqtemps(:, :), eqconsts(:, :)		! Temperatures and constants for equilibrium constants in equi.  eqtemps(i, j) is species i and temperature j
	integer nrxns							! Number of reactions
	integer, allocatable :: rxntypes(:)				! Type of each reaction
	character(NAME_MAX), allocatable :: reactants(:, :)		! Reactants for all reactions ("extra" chemicals for type 4)
	character(NAME_MAX), allocatable :: products(:, :)		! Products for all reactions
	real*8, allocatable :: stoichiometries(:, :, :)			! Stoichiometry for all reaction types.  Some are tricky...
	real*8, allocatable :: rates(:)					! Rate coefficient for reaction types 1, 2, 7, 8; forward rate coefficient for type 3; qm for type 4; halflife for reaction type 5
	real*8, allocatable :: rates2(:)				! Reverse rate coefficient for reaction type 3; decay coefficient for type 4
	character(NAME_MAX), allocatable :: dscoefs(:)			! Distribution coefficient for reaction types 1, 2; solubility coefficient for types 7, 8
	real*8, allocatable :: xcoef(:)					! Extra coefficients:  maxconc for reaction type 2; surface area for types 7, 8
	character(NAME_MAX), allocatable :: bioparams(:, :)		! Parameters for reaction type 4.  bioparams(i, j) is reaction i, substrate name at j = 1, electron acceptor name at j = 2, biomass name at j = 3, ks at j = 4, ka at j = 5, phthresh at j = 6, yield at j = 7, xminit at j = 8
	logical fsubs, felec, fbiomass					! Flags for whether the required fields of reaction type for have been provided
	character(NAME_MAX), allocatable :: icbioholder(:, :)		! Holds component names for ICBIO for reaction type 4
	real*8, allocatable :: porchange(:, :)				! Hold porosity change information for reaction type 8.  porchange(i, j) is reaction i, weight at j = 1 and density at j = 2
	integer nzones							! Number of zones
	character(NAME_MAX), allocatable :: zones(:)			! List of zones specified in the zone section
	character(NAME_MAX), allocatable :: zonespecs(:)		! List of keywords specified in the zone section (necessary for ordering the data in the grid)
	integer, allocatable :: zoneresolv(:)				! List of the row in zonegrid corresponding to each zone number (izone); list of zone numbers
	character(NAME_MAX), allocatable :: zonegrid(:, :)		! The assignments for each zone.  zonegrid(i, j) is zone i and keyword j
	integer nzonespecs						! The number specifications for each zone, the maximum zone number
	character(NAME_MAX), allocatable :: rxns(:, :)			! The set of reactions for each zone
	integer, allocatable :: znrxns(:)				! The number of reactions for each zone
	character(NAME_MAX), allocatable :: spnam(:)			! We may be able to remove this later for 
	integer zinit, zboun, zrock, zrxn, zdisp, zsorp, ztpor, &	! 
		ztime, zsehdiff, zhparam, zhenry, zvap, zopt		! Column numbers for each specification in assign
	logical zfaccum, zfconst					! Flags for solute accumulation and constant concentration
	logical inflag, bflag						! Flags for whether we're using the default ("*") init and bound conditions.  False means we are
	integer rxnnsolid, rxnnaqueous, rxnnvapor			! Number of solid, aqueous, and vapor species for reactions
	character(NAME_MAX), allocatable :: rxnsolid(:), &		! 
		rxnaqueous(:), rxnvapor(:)				! List of species in reactions according to type
	real*8, allocatable :: heqfit(:, :)				! The data for the PHREEQC fit
	character(LINE_MAX) :: zonefile					! The file containing the zones for moles
	integer nmspecs, nmzones					! Number of species and zones in moles
	character(NAME_MAX), allocatable :: mspecs(:)			! Species listed in moles
	integer, allocatable :: mzones(:)				! Zones listed in moles
	real*8, allocatable :: molegrid(:, :)				! Grid in moles.  molegrid(i, j) is zone i and species k
	character(LINE_MAX) :: usercfname				! Userc input filename
	
	! Variables provided for the copy-and-pasted portions of rdcon and read_rxn
	real*8 ctol
	integer ic, ii, iimm, im, iret, ivap, iv, ix, ncp, icpnt

		!!!!!
		!
		! readtracrxn by Matthew Schauer <mschauer@lanl.gov>
		! Initial implementation June 2010 - August 2012
		! Property of Los Alamos National Security, LLC
		!
		!!!!!

	! Preset everything to 0 or default values
	linen = 0
	fcomp = .false.
	fwater = .false.
	frock = .false.
	fdisp = .false.
	fassign = .false.
	fvap = .false.
	fheader = .false.
	fsorp = .false.
	fassign = .false.
	fhparam = .false.
	fdiff = .false.
	fspec = .false.
	fgroup = .false.
	fprint = .false.
	fdist = .false.
	fsol = .false.
	fcplx = .false.
	fequi = .false.
	fmoles = .false.
	fuserc = .false.
	libread = 0
	nrxns = 0
	hnhtspecies = 0
	ncdenspecs = 0
	naq = 0
	nso = 0
	nh = 0
	nv = 0
	array = '*'
	an = 0
	ldiff = 0
	vdiff = 0
	ldiffm = 0
	vdiffm = 0
	iskip = 0
	rsdmax = 1e-9
	strac_max = 0.99
!	co2_couple = 0
	logten = .false.
	cequi = .false.
	molein = .false.
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
	nrt = 0
	nrtspecies = 0
	nwt = 0
	nwtspecies = 0
	ncp = 0
	nvt = 0
	nvtspecies = 0
	nspecies = 0
	naq = 0
	nso = 0
	nh = 0
	nv = 0
	anspecies = 0
	nsorp = 0
	ndisp = 0
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
	dat = 20
	allocate(heqfit(ncplx, 5))
	! Allocate reaction blocks
	allocate(rxntypes(numrxn))
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
	! Check for a zone macro
	nzonespecs = 0
	if (numzones .le. 0) then
		write(ierr, '(a)') 'A zone macro with at least one valid zone must precede trxn.'
		goto 269
	endif
	if (debug) then
		write(iptty, '(i3, a19)', advance='no') numzones, ' available zones (zonenames):  '
		do i = 1, numzones
			write(iptty, '(a)', advance='no') trim(zonenames(i))
			write(iptty, '(a)', advance='no') ' '
		enddo
		write(iptty, *)
	endif
	allocate(zonespecs(SPEC_MAX))
	allocate(zoneresolv(0:zonemax))
	! Check for a time macro
	if (tims .le. 0) then
		write(ierr, '(a)') 'A time macro with a stop time greater than zero must precede trxn.'
		goto 269
	endif
	goto 270
269	write(ierr, '(a)') 'Error in trxn in preprocessing.'
	write(iptty, '(a)') 'Error in trxn in preprocessing.'
	stop
	! Read main grids
270	if (debug) write(iptty, '(a)') 'Reading data for trxn.'
	reading = .true.
	do
100		call readline ! First, read the line.
		! Then determine the contents of the line.  It could be...
		! a blank line, in which case it is ignored.
		if (len_trim(line) .eq. 0) goto 2000
		! a line beginning with a keyword, in which case different modes are entered depending on the keyword:
		read(line, *, err=667) keyword
		
		! Header information
		if (keyword .eq. 'ctrl') then
			if (debug) write(iptty, '(a)') 'Reading control information...' ! Already taken care of by trxninit
			goto 2000
		elseif (keyword .eq. 'include') then
			libfname = '*'
			if (index(line, ' ') .eq. 0) then
				write(ierr, '(a)') 'Error:  Missing include filename or bad input format in include.  '// &
					'Please note that TAB cannot be used to separate the keyword from the filename.'
				goto 666
			endif
			libfname = line(index(line, ' ') + 1:len_trim(line))
			if (debug) write(iptty, '(a)') 'Reading from library "'//trim(libfname)//'".'
			libread = libread + 1
			libnum = 23 + libread
			open(libnum, file=libfname, err=339)
			if (libread .eq. 1) inpttmp = inpt
			inpt = libnum
			rewind(libnum)
			call readline
			do
				if (len_trim(line) .ne. 0) exit
				call readline
			enddo
			goto 100
801			if (libread .eq. 1) then
				inpt = inpttmp
			else
				inpt = inpt - 1
			endif
			close(libnum)
			if (debug) write(iptty, '(a)') 'Library read done.'
			libread = libread - 1
			goto 2000
339			write(ierr, '(a)') 'Error opening include file "'//trim(libfname)//'".'
			goto 666
		elseif (keyword .eq. 'header') then
			if (debug) write(iptty, '(a)') 'Reading header information...  '
			if (fheader) then
				write(ierr, '(a)') 'Warning:  header has already been specified.  '// &
					'The new values will overwrite the old ones.'
			endif
			fheader = .true.
			call readline
			read(line, *, end=118, err=118) an0, awc, epc, upwgta ! Verbatim from trac
			call readline
			read(line, *, end=118, err=118) daycs, daycf, dayhf, dayhs
			call readline
			read(line, *, end=118, err=118) iaccmx, daycm, daycmm, daycmx, nprttrc
			if (awc .lt. 1) then
				awc = 1
			else if (awc .gt. 1) then
				awc = 0.5
			endif
			if (upwgta .lt. 0.5) then
				write(ierr, '(a)') 'Warning:  upwgta must be from 0.5 to 1.  Value contrained to 0.5.'
				upwgta = 0.5
			elseif (upwgta .gt. 1) then
				write(ierr, '(a)') 'Warning:  upwgta must be from 0.5 to 1.  Value constrained to 1.'
				upwgta = 1
			endif
			if (dayhs .lt. dayhf) then
				write(ierr, '(a)') 'Warning:  dayhs cannot be earlier than dayhf.  Value constrained to dayhf.'
				dayhs = dayhf
			endif
			call readline
			if (len_trim(line) .eq. 0) goto 2000
			array = '*'
			read(line, *, end=138, err=118) (array(i), i = 1, SPEC_MAX) ! Next, read optional variables.
138				do i = 1, SPEC_MAX
				if (array(i) .eq. '*') exit
				if (index(array(i), '=') .eq. 0) then
					write(ierr, '(a)') 'Unknown token in header:  "'//trim(array(i))//'"'
					goto 666
				endif
				if (array(i)(1:index(array(i), '=') - 1) .eq. 'iskip') then
					read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *, err=118) iskip
				elseif (array(i)(1:index(array(i), '=') - 1) .eq. 'rsdmax') then
					read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *,err=118) rsdmax
				elseif (array(i)(1:index(array(i), '=') - 1) .eq. 'strac_max') then
					read(array(i)(index(array(i), '=') + 1:len_trim(array(i))), *, err=118) strac_max
				else
					write(ierr, '(a)') 'Error:  Variable "'//array(i)(1:index(array(i), '=') - 1)// &
						'" in header not recognized.'
					goto 666
				endif
			enddo
			if (iskip .lt. 0) then
				write(ierr, '(a)') 'Error:  iskip cannot be negative.'
				goto 666
			endif
			if (rsdmax .lt. 0) then
				write(ierr, '(a)') 'Error:  rsdmax cannot be negative.'
				goto 666
			endif
			if (strac_max .lt. 0) then
				write(ierr, '(a)') 'Error:  strac_max cannot be negative.'
				goto 666
			endif
			daycmx = min(daycmx, daymax)
			dtotc = daycmm * 86400
			ayc = 1 - awc
			goto 2000
118			write(ierr, '(a)') 'Error reading header information.'
			goto 666
		! Userc call
		elseif (keyword .eq. 'userc') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading user subroutine filename...  '
			if (fuserc) write(ierr, '(a)') 'Warning:  userc has already been specified.  '// &
				'The new filename will be used rather than the old one.'
			fuserc = .true.
			if (index(line, ' ') .eq. 0) then
				write(ierr, '(a)') 'Error:  userc input filename missing.  '// &
					'Please note that TAB cannot be used to separate the keyword from the file name.'
				goto 666
			endif
			usercfname = trim(adjustl(line(index(line, ' ') + 1:len_trim(line))))
			if (debug) write(iptty, '(a)') 'User subroutine will be called after reading trxn.'
			goto 2000
		! Component definitions
		elseif (keyword .eq. 'comp') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading component data...  '
			if (fcomp) then
				write(ierr, '(a)') 'Warning:  comp has already been specified.  '// &
					'The new values will overwrite the old ones.'
!				deallocate(species)
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
			linen = linen - 1
			call readline
			array = '*'
			read(line, *, end=198, err=600) keyword, (array(i), i = 1, SPEC_MAX)
198			do
				call readline
				if (len_trim(line) .eq. 0) exit
				nspecies = nspecies + 1
			enddo
			do i = 1, nspecies + 1
				backspace inpt
				linen = linen - 1
			enddo
			if (flookup) nspecies = nspecies + nlocomp
!			allocate(species(nspecies))
			allocate(states(nspecies))
			allocate(masters(nspecies))
			allocate(guesses(nspecies))
			masters = '*'
			guesses = 1e-9
			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
					exit
				elseif ((array(i)(1:6) .ne. 'master') .and. (array(i)(1:5) .ne. 'guess')) then
					write(ierr, '(a)') 'Column header "'//trim(array(i))//'" not recognized in comp.'
					goto 666
				endif
			enddo
			j = 0
			do i = 1, nspecies
				call readline
				if (len_trim(line) .eq. 0) then
					j = i - 1
					exit
				endif
				if (array(1)(1:6) .eq. 'master') then
					if (array(2)(1:5) .eq. 'guess') then
						read(line, *, err=167, end=167) states(i), species(i), masters(i), e
						if(states(i).eq.'a'.and.species(i)(1:1).eq.'H') ph_species=i
						if (e .eq. '*') then
							guesses(i) = 1e-9
						else
							read(e, *, err=167) guesses(i)
						endif
					else
						read(line, *, err=167, end=167) states(i), species(i), masters(i)
						guesses(i) = 1e-9
					endif
				elseif (array(1)(1:5) .eq. 'guess') then
					if (array(2)(1:6) .eq.'master') then
						read(line, *, err=167, end=167) states(i), species(i), e, masters(i)
					else
						read(line, *, err=167, end=167) states(i), species(i), e
					endif
					if (e .eq. '*') then
						guesses(i) = 1e-9
					else
						read(e, *, err=167) guesses(i)
					endif
				else
					read(line, *, err=167, end=167) states(i), species(i)
				endif
				if (states(i) .eq. 'a') then
					naq = naq + 1
				elseif (states(i) .eq. 's') then
					nso = nso + 1
				elseif (states(i) .eq. 'h') then
					nh = nh + 1
				elseif (states(i) .eq. 'g') then
					nv = nv + 1
				else
					write(ierr, '(a)') 'Error:  Component state "'//trim(states(i))//'" not recognized.'
					goto 666
				endif
				if (guesses(i) .lt. 0) guesses(i) = 10 ** guesses(i)
				if ((states(i) .ne. 'a') .and. (states(i) .ne. 'h') .and. (masters(i) .ne. '*')) then
					write(ierr, '(a)') 'Error:  Component "'//trim(species(i))//'" is not '// &
						'aqueous or Henry''s Law and therefore cannot have a master species.'
					goto 666
				endif
			enddo
			if (flookup) then
				do i = 1, nlocomp
					states(j + i) = locomp(i, 1)
					if (locomp(i, 1) .eq. 'a') then
						naq = naq + 1
					elseif (locomp(i, 1) .eq. 's') then
						nso = nso + 1
					elseif (locomp(i, 1) .eq. 'v') then
						nv = nv + 1
					elseif (locomp(i, 1) .eq. 'h') then
						nh = nh + 1
					else
						write(ierr, '(a)') 'Internal error setting states of components '// &
							'added by lookup.'
						goto 667
					endif
					species(j + i) = locomp(i, 2)
					masters(j + i) = locomp(i ,3)
				enddo
			endif
			nmasters = 0
			do i = 1, nspecies
				if (masters(i) .ne. '*') then
					nmasters = nmasters + 1
				endif
			enddo
			array = '*' ! Don't let a component be defined more than once
			do i = 1, nspecies
				do j = 1, i - 1
					if (species(i) .eq. species(j)) then
						write(ierr, '(a)') 'Error:  Multiple definition of component "'// &
							trim(species(i))//'" in comp.'
						goto 666
					endif
				enddo
			enddo
			array = '*' ! Master species must be unique
			do i = 1, nspecies
				do j = 1, i - 1
					if ((masters(i) .eq. masters(j)) .and. (masters(i) .ne. '*')) then
						write(ierr, '(a)') 'Error:  Master species "'//trim(masters(i))// &
							'" is used for two components, "'//trim(species(i))// &
							'" and "'//trim(species(j))//'".'
						goto 666
					endif
				enddo
			enddo
			ncp = nspecies
			if (debug) write(iptty, '(a, i3, a)') 'Read ', nspecies, ' components.'
			goto 2000
167			write(ierr, '(a)') 'Error reading components in comp.'
			goto 666
		! Read molecular weights
		elseif (keyword .eq. 'cden') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading molecular weight table...  '
			if (cden_flag .ne. 0) then
				write(ierr, '(a)') 'Warning:  A cden macro may have been previously specified.  '// &
					'Use of the cden block in trxn may overwrite these values.'
			endif
			if (cden) then
				write(ierr, '(a)') 'Warning:  cden has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(cdenspecs)
				deallocate(molweights)
			endif
			
			cden = .true.
			cden_flag = 1
			ncdenspecs = 0
			flag = .false.
			j = 0
			do
				call readline
				if (len_trim(line) .eq. 0) exit
				if (line .eq. '*') then
					if (.not. flookup) then
						write(ierr, '(a)') 'Error:  A lookup block must be provided '// &
							'before cden if dynamic molecular weight lookup is to be used.'
						goto 666
					endif
					j = j + 1
					flag = .true.
				else
					ncdenspecs = ncdenspecs + 1
				endif
			enddo
			do i = 1, ncdenspecs + j
				backspace inpt
				linen = linen - 1
			enddo
			if (flag) then
				allocate(cdenspecs(ncdenspecs + nlocompaq))
				allocate(molweights(ncdenspecs + nlocompaq))
			else
				backspace inpt
				linen = linen - 1
				allocate(cdenspecs(ncdenspecs))
				allocate(molweights(ncdenspecs))
			endif
			
			do i = 1, ncdenspecs
328				call readline
				cdenspecs(i) = '*'
				e = '*'
				read(line, *, end=286, err=271) cdenspecs(i), e
286				if (cdenspecs(i) .eq. '*') then
					goto 328
				else
					read(e, *, end=271, err=271) molweights(i)
				endif
322				if (molweights(i) .lt. 0) then
					write(ierr, '(a)') 'Error:  Molecular weight for component "'//cdenspecs(i)// &
						'" cannot be negative.'
					goto 666
				endif
			enddo
			if (flag) then
				k = 0
				do i = 1, nlocomp
					if (locomp(i, 1) .eq. 'a') then
						k = k + 1
						cdenspecs(ncdenspecs + k) = locomp(i, 2)
						do j = 1, ndatcomp
							if (datcomp(j) .eq. locomp(i, 2)) then
								molweights(ncdenspecs + k) = datcden(j)
								goto 272
							endif
						enddo
						write(ierr, '(a)') 'Internal error looking up component molecular weights.'
						goto 667
272						continue
					endif
				enddo
				ncdenspecs = ncdenspecs + nlocompaq
			endif
			if (debug) write(iptty, '(a, i3, a)') 'Read molecular weights for ', ncdenspecs, ' aqueous components.'
			goto 2000
271			write(ierr, '(a)') 'Error reading molecular weights in cden.'
			goto 666
		! Water type mode
		elseif (keyword .eq. 'water') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading water type data...  '
			if (fwater) then
				write(ierr, '(a)') 'Warning:  water has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(wtspecies)
				deallocate(wtnames)
				deallocate(wtgrid)
			endif
			fwater = .true.
			array = '*'
			read(line, *, end=106, err=131) keyword, (array(i), i = 1, SPEC_MAX)
106			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
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
				call readline
				if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
					goto 105		
				endif
				nwt = nwt + 1
			enddo
105			do i = 1, nwt + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(wtnames(nwt))
			allocate(wtgrid(nwt, nwtspecies))
			do i = 1, nwt
				call readline
				read(line, *, end=131, err=131) wtnames(i), (wtgrid(i, j), j = 1, nwtspecies)
				do j = 1, nwtspecies
					if (wtgrid(i, j) .lt. 0) wtgrid(i, j) = 10 ** wtgrid(i, j)
				enddo
			enddo
			if (debug) write(iptty, '(i3, a, i3, a)') nwt, ' water type(s), ', nwtspecies, &
				' species read.'
			goto 2000
131			write(ierr, '(a)') 'Error reading values in water.'
			goto 666
		! Rock mode -- Note that aside from the variable names, this is identical to the water reader.
		elseif (keyword .eq. 'rock') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading solid species data...  '
			if (frock) then
				write(ierr, '(a)') 'Warning:  rock has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(rtnames)
				deallocate(rtspecies)
				deallocate(rtgrid)
			endif
			frock = .true.
			array = '*'
			read(line, *, end=107, err=129) keyword, (array(i), i = 1, SPEC_MAX)
107			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
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
				call readline
				if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
					goto 109		
				endif
				nrt = nrt + 1
			enddo
109			do i = 1, nrt + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(rtnames(nrt))
			allocate(rtgrid(nrt, nrtspecies))
			do i = 1, nrt
				call readline
				read(line, *, end=129, err=129) rtnames(i), (rtgrid(i, j), j = 1, nrtspecies)
				do j = 1, nrtspecies
					if (rtgrid(i, j) .lt. 0) rtgrid(i, j) = 10 ** rtgrid(i, j)
				enddo
			enddo
			if (debug) write(iptty, '(i3, a15, i3, a14)') nrt, ' rock type(s), ', nrtspecies, &
				' species read.'
			goto 2000
129			write(ierr, '(a)') 'Error reading values in rock.'
			goto 666
		! Vapor species mode -- This is also very similar to the water reader
		elseif (keyword .eq. 'gas') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading gaseous species data...  '
			if (fvap) then
				write(ierr, '(a)') 'Warning:  gas has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(vtnames)
				deallocate(vtspecies)
				deallocate(vtgrid)
			endif
			fvap = .true.
			array = '*'
			read(line, *, end=160, err=153) keyword, (array(i), i = 1, SPEC_MAX)
160			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
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
				call readline
				if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
					goto 155		
				endif
				nvt = nvt + 1
			enddo
155			do i = 1, nvt + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(vtnames(nvt))
			allocate(vtgrid(nvt, nvtspecies))
			do i = 1, nvt
				call readline
				read(line, *, end=153, err=153) vtnames(i), (vtgrid(i, j), j = 1, nvtspecies)
				do j = 1, nvtspecies
					if (vtgrid(i, j) .lt. 0) vtgrid(i, j) = 10 ** vtgrid(i, j)
				enddo
			enddo
			if (debug) write(iptty, '(i3, a, i3, a)') nvt, ' gas type(s), ', nvtspecies, &
				' species read.'
			goto 2000
153			write(ierr, '(a)') 'Error reading values in gas.'
			goto 666
		! Mole input
		elseif (keyword .eq. 'moles') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading contant-mole species data...  '
			if (fmoles) then
				write(ierr, '(a)') 'Warning:  moles has already been specified.  '// &
					'The new values will overwrite the old ones.'
				open(20, file='moles_in.tmp')
				close(20, status='delete')
				deallocate(mspecs)
				deallocate(mzones)
				deallocate(molegrid)
			endif
			fmoles = .true.
			molein = .true.
			zonefile = '*'
			if (index(line, ' ') .ne. 0) zonefile = trim(adjustl(line(index(line, ' ') + 1:len_trim(line))))
			array = '*'
			read(line, *, end=333, err=332) (array(i), i = 1, SPEC_MAX)
333			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
					nmspecs = i - 1
					exit
				endif
			enddo
			allocate(mspecs(nmspecs))
			do i = 1, nmspecs
				mspecs(i) = array(i)
			enddo
			nmzones = 0
			do
				call readline
				if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '	') then
					goto 334		
				endif
				nmzones = nmzones + 1
			enddo
334			do i = 1, nmzones + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(mzones(nmzones))
			allocate(molegrid(nmzones, nmspecs))
			do i = 1, nmzones
				call readline
				read(line, *, end=332, err=332) mzones(i), (molegrid(i, j), j = 1, nmspecs)
				do j = 1, nmspecs
					if (molegrid(i, j) .lt. 0) molegrid(i, j) = 10 ** molegrid(i, j)
				enddo
			enddo
			if (debug) write(iptty, '(i3, a, i3, a)') nmzones, ' zone(s), ', nmspecs, &
				' species read.'
			goto 2000
332			write(ierr, '(a)') 'Error reading values in moles.'
			goto 666
		! Henry's Law parameter mode
		elseif (keyword .eq. 'hparam') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading Henry''s Law parameters...  '
			if (fhparam) then
				write(ierr, '(a)') 'Warning:  hparam has already been specified.  '// &
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
				call readline
				if (len_trim(line) .eq. 0) exit
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
				call readline
				if (len_trim(line) .eq. 0) exit
				read(line, *, err=145, end=147) (array(j), j = 1, SPEC_MAX)
147				hhtspecies(i) = array(1)
				if ((array(2)(1:1) .eq. 'h') .or. (array(2) .eq. '1') .or. (array(2) .eq. '*')) then
					hmodels(i) = 1
				elseif ((array(2)(1:1) .eq. 'm') .or. (array(2) .eq. '2')) then
					hmodels(i) = 2
				elseif ((array(2)(1:1) .eq. 'w') .or. (array(2) .eq. '3')) then
					hmodels(i) = 3
				else
					write(ierr, '(a)') 'Error:  Unrecognized Henry''s model "'//trim(array(2))// &
						'" in hparam.'
					goto 666
				endif
				do j = 3, SPEC_MAX
					if (array(j) .eq. '*') then
						exit
					endif
					e = array(j)(1: index(array(j), '=') - 1)
					f = array(j)(index(array(j), '=') + 1:len_trim(array(j)))
					if (hmodels(i) .eq. 1) then
						if (e .eq. 'ah') then
							read(f, *, err=148, end=148) ah(i)
						elseif (e .eq. 'dhh') then
							read(f, *, err=148, end=148) dhh(i)
						elseif ((e(1:2) .eq. 'ah') .or. (e .eq. 'hh')) then
							
						else
							write(ierr, '(a)') 'Error:  Unknown hparam variable "'//trim(e)//'".'
							goto 666
						endif
					elseif (hmodels(i) .eq. 2) then
						if (e .eq. 'ah1') then
							read(f, *, err=148, end=148) ah1(i)
						elseif (e .eq. 'ah2') then
							read(f, *, err=148, end=148) ah2(i)
						elseif (e .eq. 'ah3') then
							read(f, *, err=148, end=148) ah3(i)
						elseif (e .eq. 'ah4') then
							read(f, *, err=148, end=148) ah4(i)
						elseif (e .eq. 'ah5') then
							read(f, *, err=148, end=148) ah5(i)
						elseif ((e .eq. 'ah') .or. (e .eq. 'dhh') .or. (e .eq. 'hh')) then
							
						else
							write(ierr, '(a)') 'Error:  Unknown hparam variable "'//trim(e)//'".'
							goto 666
						endif
					elseif (hmodels(i) .eq. 3) then
						if (e .eq. 'ah') then
							read(f, *, err=148, end=148) ah(i)
						elseif (e .eq. 'hh') then
							read(f, *, err=148, end=148) hh(i)
						elseif ((e(1:2) .eq. 'ah') .or. (e .eq. 'dhh')) then
							
						else
							write(ierr, '(a)') 'Error:  Unknown hparam variable "'//trim(e)//'".'
							goto 666
						endif
					else
						write(ierr, '(a)') 'Internal error processing hparam arguments.'
						goto 667
					endif
				enddo
			enddo
			if (debug) write(iptty, '(i3, a)') hnhtspecies, ' models read.'
			goto 2000
145			write(ierr, '(a)') 'Error reading Henry''s Law data.'
			goto 666
148			write(ierr, '(a)') 'Error:  Invalid number "'//trim(f)//'" in henry.'
			goto 666
		! Aqueous reaction mode - database lookup
		elseif (keyword .eq. 'lookup') then
			if (debug) write(iptty, '(a)') 'Using lookup database '// &	! All processing done in trxninit
				trim(adjustl(line(index(line, ' ') + 1:len_trim(line))))//'.'
			goto 2000
		! Adsorption mode
		elseif (keyword .eq. 'sorp') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading adsorption data...  '
			if (fsorp) then
				write(ierr, '(a)') 'Warning:  sorp has already been specified.  '// &
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
			linen = linen - 1
			call readline
			array = '*'
			read(line, *, end=122, err=124) keyword, (array(i), i = 1, SPEC_MAX)
122			anspecies = 0
			nsorp = 0
			a = 0
			do
				call readline
				if (len_trim(line) .eq. 0) then
					anspecies = a
					exit
				endif
				if (line(1:1) .eq. '.') then
					nsorp = nsorp + 1
					if ((anspecies .ne. 0) .and. (anspecies .ne. a)) then
						write(ierr, '(a, i3, a, i3, a)') 'Error:  Different number of species '// &
							'specified in different adsorption models:  ', a, ' and ', &
							anspecies, '.'
						goto 666
					endif
					anspecies = a
					a = 0
				else
					a = a + 1
				endif
			enddo
			if (nsorp .eq. 1) then
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
			do i = 1, nsorp
				call readline
				if (line(1:1) .ne. '.') then
					write(ierr, '(a)') 'Error reading adsorption models.'
					goto 666
				endif
				sorpnames(i) = line(2:len_trim(line))
				if (index(line, ' ') .ne. 0) then
					sorpnames(i) = sorpnames(i)(1:index(sorpnames(i), ' '))
				endif
				do k = 1, anspecies
					array2 = '*'
					call readline
					read(line, *, end=123, err=124) array3(i, k), (array2(j), j = 1, SPEC_MAX)
123 					do j = 1, SPEC_MAX
						if (array(j) .eq. 'ltype') then
							fltype = .true.
							if ((array2(j) .eq. '0') .or. (array2(j)(1:2) .eq. 'co') .or. &
								(array2(j) .eq. '*')) then
								lsorptypes(i, k) = 0
							elseif ((array2(j) .eq. '1') .or. (array2(j)(1:2) .eq. 'li')) then
								lsorptypes(i, k) = 1
							elseif ((array2(j) .eq. '2') .or. (array2(j)(1:2) .eq. 'fr')) then
								lsorptypes(i, k) = 2
							elseif ((array2(j) .eq. '3') .or. (array2(j)(1:2) .eq. 'mf')) then
								lsorptypes(i, k) = 3
							elseif ((array2(j) .eq. '4') .or. (array2(j)(1:2) .eq. 'la')) then
								lsorptypes(i, k) = 4
							else
								write(ierr, '(a)') 'Error:  Liquid adsorption type "'// & 
									trim(array2(j))//'" in sorp not recognized.'
								goto 666
							endif
						elseif (array(j) .eq. 'vtype') then
							fvtype = .true.
							if ((array2(j) .eq. '0') .or. (array2(j)(1:2) .eq. 'co') .or. &
								(array2(j) .eq. '*')) then
								vsorptypes(i, k) = 0
							elseif ((array2(j) .eq. '1') .or. (array2(j)(1:2) .eq. 'li')) then
								vsorptypes(i, k) = 1
							elseif ((array2(j) .eq. '2') .or. (array2(j)(1:2) .eq. 'fr')) then
								vsorptypes(i, k) = 2
							elseif ((array2(j) .eq. '3') .or. (array2(j)(1:2) .eq. 'mf')) then
								vsorptypes(i, k) = 3
							elseif ((array2(j) .eq. '4') .or. (array2(j)(1:2) .eq. 'la')) then
								vsorptypes(i, k) = 4
							else
								write(ierr, '(a)') 'Error:  Vapor adsorption type "'// &
									trim(array2(j))//'" in sorp not recognized.'
								goto 666
							endif
						elseif (array(j) .eq. 'a1l') then
							fa1l = .true.
							if (array2(j) .eq. '*') then
								a1l(i, k) = 0
							else
								read(array2(j), *, err=124) a1l(i, k)
							endif
						elseif (array(j) .eq. 'a2l') then
							fa2l = .true.
							if (array2(j) .eq. '*') then
								a2l(i, k) = 0
							else
								read(array2(j), *, err=124) a2l(i, k)
							endif
						elseif (array(j) .eq. 'bl') then
							fbl = .true.
							if (array2(j) .eq. '*') then
								bl(i, k) = 0
							else
								read(array2(j), *, err=124) bl(i, k)
							endif
						elseif (array(j) .eq. 'a1v') then
							fa1v = .true.
							if (array2(j) .eq. '*') then
								a1v(i, k) = 0
							else
								read(array2(j), *, err=124) a1v(i, k)
							endif
						elseif (array(j) .eq. 'a2v') then
							fa2v = .true.
							if (array2(j) .eq. '*') then
								a2v(i, k) = 0
							else
								read(array2(j), *, err=124) a2v(i, k)
							endif
						elseif (array(j) .eq. 'bv') then
							fbv = .true.
							if (array2(j) .eq. '*') then
								bv(i, k) = 0
							else
								read(array2(j), *, err=124) bv(i, k)
							endif
						elseif (array(j) .eq. '*') then
							exit
						else
							write(ierr, '(a)') 'Error:  Column "'//trim(array(j))// &
								'" in sorp not recognized.'
							goto 666
						endif
					enddo
				enddo
			enddo
			if (.not. fltype) lsorptypes = 0
			if (.not. fvtype) vsorptypes = 0
			if (.not. fa1l) a1l = 0
			if (.not. fa2l) a2l = 0
			if (.not. fbl) bl = 0
			if (.not. fa1v) a1v = 0
			if (.not. fa2v) a2v = 0
			if (.not. fbv) bv = 0
			do i = 1, anspecies
				sorpspecs(i) = array3(1, i)
				do j = 1, nsorp
					if (sorpspecs(i) .ne. array3(j, i)) then
						write(ierr, '(a)') 'Error:  Species do not match across all models in sorp.'
						goto 666
					endif
				enddo
			enddo
			if (debug) write(iptty, '(i3, a13, i3, a26)') anspecies, ' species for ', &
				nsorp, ' adsorption model(s) read.'
			goto 2000
124			write(ierr, '(a)') 'Error reading a number in sorp.'
			goto 666
		! Diffusivity mode
		elseif (keyword .eq. 'diff') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading diffusion information...  '
			if (fdiff) then
				write(ierr, '(a)') 'Warning:  diff has already been specified.  '// &
					'The new values will overwrite the old ones.'
			endif
			fdiff = .true.
			array = '*'
			read(line, *, end=162, err=151) keyword, (array(i), i = 1, SPEC_MAX)
162			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') exit
				e = array(i)
				f = e(1:index(e, '=') - 1)
				g = e(index(e, '=') + 1:len_trim(e))
				if ((f .eq. '') .or. (g .eq. '')) then
					write(ierr, '(a)') 'Error:  Token "'//trim(e)// &
						'" does not have format "variable=value".'
					goto 666
				endif
				if (f .eq. 'l') then
					read(g, *, end=151, err=151) ldiff
				elseif (f .eq. 'v') then
					read(g, *, end=151, err=151) vdiff
				elseif (f .eq. 'lm') then
					if ((g .eq. '0') .or. (g .eq. '*') .or. (g .eq. 'con')) then
						ldiffm = 0
					elseif ((g .eq. '1') .or. (g .eq. 'mq')) then
						ldiffm = 1
					elseif ((g .eq. '2') .or. (g .eq. 'cw')) then
						ldiffm = 2
					else
						write(ierr, '(a)') 'Error:  Liquid diffusion model "'//trim(g)// &
							'" not recognized.'
						goto 666
					endif
				elseif (f .eq. 'vm') then
					if ((g .eq. '0') .or. (g .eq. '*') .or. (g .eq. 'con')) then
						vdiffm = 0
					elseif ((g .eq. '1') .or. (g .eq. 'mq')) then
						vdiffm = 1
					elseif ((g .eq. '2') .or. (g .eq. 'cw')) then
						vdiffm = 2
					else
						write(ierr, '(a)') 'Error:  Vapor diffusion model "'//trim(g)// &
							'" not recognized.'
						goto 666
					endif
				else
					write(ierr, '(a)') 'Error:  Parameter "'//trim(f)//'" not recognized.'
					goto 666
				endif
			enddo
			if (debug) write(iptty, '(a)') 'Read diffusion values.'
			goto 2000
151			write(ierr, '(a)') 'Error reading diffusion values.'
			goto 666
		! Dispersivity mode
		elseif (keyword .eq. 'disp') then
			if (fdisp) then
				write(ierr, '(a)') 'Warning:  disp has already been specified.  '// &
					'The new values will overwrite the old ones.'
			endif
			fdisp = .true.
			dispmode = -1 ! 0 is XYZ, 1 is LT
			backspace inpt
			call readline
			dispparams = '*'
			read(line, *, end=164, err=165) keyword, (dispparams(i), i = 1, SPEC_MAX)
164			do i = 1, SPEC_MAX
				if (dispparams(i) .eq. '*') exit
				if ((dispparams(i)(2:2) .eq. 'x') .or. (dispparams(i)(2:2) .eq. 'y') .or. &
					(dispparams(i)(2:2) .eq. 'z')) then
					if (dispmode .eq. 1) then
						write(ierr, '(a)') 'Error:  Dispersivity column header "' &
							//trim(dispparams(i))//'" conflicts with previous '// &
							'longitudinal/transverse column headers.'
						goto 666
					endif
					dispmode = 0
				elseif ((dispparams(i)(2:2) .eq. 'l') .or. (dispparams(i)(2:2) .eq. 't')) then
					if (dispmode .eq. 0) then
						write(ierr, '(a)') 'Error:  Dispersivity column header "' &
							//trim(dispparams(i))//'" conflicts with previous XYZ '// &
							'column headers.'
						goto 666
					endif
					dispmode = 1
				else
					write(ierr, '(a)') 'Error:  Unrecognized header "'// trim(dispparams(i))// &
						'" in disp.'
					goto 666
				endif
				if ((dispparams(i)(1:1) .ne. 'l') .and. (dispparams(i)(1:1) .ne. 'v')) then
					write(ierr, '(a)') 'Error:  Unrecognized header "'//trim(dispparams(i))// &
						'" in disp.'
					goto 666
				endif
			enddo
			if (dispmode .eq. 0) then
				if (debug) write(iptty, '(a)', advance='no') 'Reading XYZ dispersion models...  '
			elseif (dispmode .eq. 1) then
				if (debug) write(iptty, '(a)', advance='no') 'Reading longitudinal/'// &
					'transverse dispersivity models...  '
			else
				write(ierr, '(a)') 'Error:  At least one column must be provided in disp.'
				goto 666
			endif
			ndisp = 0
			do
				call readline
				if (len_trim(line) .eq. 0) goto 137
				ndisp = ndisp + 1
			enddo
137			do i = 0, ndisp
				backspace inpt
				linen = linen - 1
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
				call readline
				array = '*'
				if (dispmode .eq. 0) then
					read(line, *, end=139, err=140) dispnames(i), (array(j), j = 1, SPEC_MAX)	!lx(i), ly(i), lz(i), vx(i), vy(i), vz(i)
139					do j = 1, SPEC_MAX
						if ((dispparams(j) .eq. '*')) then
							exit
						endif
						if (array(j) .eq. '*') then
							array(j) = '0'
						endif
						if (dispparams(j) .eq. 'lx') then
							read(array(j), *, err=140) lx(i)
							if (lx(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  X liquid dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'ly') then
							read(array(j), *, err=140) ly(i)
							if (ly(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Y liquid dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'lz') then
							read(array(j), *, err=140) lz(i)
							if (lz(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Z liquid dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'vx') then
							read(array(j), *, err=140) vx(i)
							if (vx(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  X gaseous dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'vy') then
							read(array(j), *, err=140) vy(i)
							if (vy(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Y gaseous dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'vz') then
							read(array(j), *, err=140) vz(i)
							if (vz(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Z gaseous dispersivity '// &
									'for dispersivity model "'//dispnames(i)// &
									'" cannot be negative.'
								goto 666
							endif
						else
							write(ierr, '(a)') 'Unknown XYZ dispsersivity header "'// &
								trim(dispparams(j))//'".'
							goto 666
						endif
					enddo
				else
					read(line, *, end=166, err=140) dispnames(i), (array(j), j = 1, SPEC_MAX)	!ll(i), lt(i), vl(i), vt(i)
166					do j = 1, SPEC_MAX
						if ((dispparams(j) .eq. '*') .or. (array(j) .eq. '*')) then
							exit
						elseif (dispparams(j) .eq. 'll') then
							read(array(j), *, err=140) ll(i)
							if (ll(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Longitudinal liquid'// &
									'dispersivity for dispersivity model "'// &
									dispnames(i)//'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'lt') then
							read(array(j), *, err=140) lt(i)
							if (lt(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Transverse liquid '// &
									'dispersivity for dispersivity model "'// &
									dispnames(i)//'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'vl') then
							read(array(j), *, err=140) vl(i)
							if (vl(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Longitudinal gaseous'// &
									'dispersivity for dispersivity model "'// &
									dispnames(i)//'" cannot be negative.'
								goto 666
							endif
						elseif (dispparams(j) .eq. 'vt') then
							read(array(j), *, err=140) vt(i)
							if (vt(i) .lt. 0) then
								write(ierr, '(a)') 'Error:  Transverse gaseous '// &
									'dispersivity for dispersivity model "'// &
									dispnames(i)//'" cannot be negative.'
								goto 666
							endif
						else
							write(ierr, '(a)') 'Unknown longitudinal/transverse '// &
								'dispersivity header "'//trim(dispparams(j))//'".'
							goto 666
						endif
					enddo
				endif
			enddo			
			if (debug) write(iptty, '(a, i3, a)') 'Read parameters for ', ndisp, ' model(s).'
			goto 2000
140			write(ierr, '(a)') 'Error reading numbers in disp.'
			goto 666
165			write(ierr, '(a)') 'Error reading column headers in disp.'
			goto 666
		! Information printing selection mode
		elseif (keyword .eq. 'print') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading component printing data...  '
			if (fprint) then
				write(ierr, '(a)') 'Warning:  print has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(printspecies)
			endif
			fprint = .true.
			array = '*'
			read(line, *, end=201, err=600) keyword, (array(i), i = 1, SPEC_MAX)
201			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
					nprint = i - 1
					exit
				endif
			enddo
			allocate(printspecies(nprint))
			do i = 1, nprint
				printspecies(i) = array(i)
				if (printspecies(i) .eq. 'all') then 
					if (debug) write(iptty, '(a)') 'All species will print.'
					deallocate(printspecies)
					allocate(printspecies(1))
					nprint = 1
					printspecies(1) = 'all'
					goto 2000
				elseif (printspecies(i) .eq. 'none') then
					if (debug) write(iptty, '(a)') 'No species will print.'
					deallocate(printspecies)
					allocate(printspecies(1))
					nprint = 1
					printspecies(1) = 'none'
					goto 2000
				endif
			enddo
			if (debug) write(iptty, '(a5, i3, a19)') 'Read ', nprint, ' printable species.'
			goto 2000
180			write(ierr, '(a)') 'Error reading printable species.'
			goto 666
		! Component grouping mode
		elseif (keyword .eq. 'group') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading component grouping data...  '
			if (fgroup) then
				write(ierr, '(a)') 'Warning:  group has already been specified.  '// &
					'The new values will overwrite the old ones.'
			endif
			do ! Already done in trxninit
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
			goto 2000
		! Aqueous complex input mode
		elseif (keyword .eq. 'cplx') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading aqueous complexes...  '
			if (fcplx) then
				write(ierr, '(a)') 'Warning:  Aqueous complexes have already been input.  '// &
					'The new values will overwrite the old ones.'
				deallocate(complexes)
				deallocate(complexcontents)
				deallocate(complexstoich)
				deallocate(equconstants)
				deallocate(enthalpies)
			endif
			fcplx = .true.
			ncomplexes = 0
			logten = .false.
			cequi = .false.
			backspace inpt
			linen = linen - 1
			call readline
			e = '*'
			f = '*'
			read(line, *, end=248, err=600) keyword, e, f
248			if ((e(1:3) .eq. 'log') .or. (f(1:3) .eq. 'log')) then
				logten = .true.
			endif
			if ((e(1:4) .eq. 'equi') .or. (f(1:4) .eq. 'equi')) then
				if (flookup) then
					write(ierr, '(a)') 'Error:  equi cannot be used if the lookup block is specified.'
					goto 666
				endif
				cequi = .true.
			endif
			do i = 1, SPEC_MAX
				call readline
				if (len_trim(line) .eq. 0) then
					ncomplexes = i - 1
					exit
				endif
			enddo
			do i = 1, ncomplexes + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(complexes(ncomplexes))
			allocate(complexcontents(ncomplexes, SPEC_MAX))
			allocate(complexstoich(ncomplexes, SPEC_MAX))
			allocate(equconstants(ncomplexes))
			allocate(enthalpies(ncomplexes))
			equconstants = 0
			enthalpies = 0
			complexcontents = '*'
			complexstoich = 0
			do i = 1, ncomplexes
				call readline
				array = '*'
				a = 1
				e = '*'
				read(line, *, end=261, err=600) complexes(i), e, (array(j), j = 1, SPEC_MAX)
261				if (e .ne. '=') then ! Lookup mode
					if (.not. flookup) then
						write(ierr, '(a)') 'Error:  A lookup block must be provided before '// &
							'cplx if dynamic complex lookup is to be used.'
						goto 666
					endif
					do j = 1, ndatcplx
						if (datcplxmain(j) .eq. complexes(i)) then
							do k = 1, maxdatcplx
								if (datcplx(j, k) .eq. '*') exit
								complexcontents(i, k) = datcplx(j, k);
								complexstoich(i, k) = datcplxstoic(j, k);
							enddo
							if (datcplxtemp(j, 1) .eq. 0) then
								equconstants(i) = datcplxlkeq(j)
								enthalpies(i) = datcplxheq(j)
								heqfit(i, 1:5) = 0
							else
								equconstants(i) = 0
								enthalpies(i) = 0
								heqfit(i, 1:5) = datcplxtemp(j, 1:5)
							endif
							goto 330
						endif
					enddo
					write(ierr, '(a)') 'Error:  Complex "'//trim(complexes(i))// &
						'" not found in lookup database.'
					goto 666
330					if (e .ne. '*') then
						if (e(1:5) .eq. 'ckeq=') then
							e = e(6:len_trim(e))
							if (e .eq. '*') then
								equconstants(i) = 0
							else
								read(e, *, err=263, end=263) equconstants(i)
								if (logten) equconstants(i) = 10 ** equconstants(i)
							endif
						elseif (e(1:4) .eq. 'heq=') then
							e = e(5:len_trim(e))
							if (e .eq. '*') then
								enthalpies(i) = 0
							else
								read(e, *, err=263, end=263) enthalpies(i)
								if (logten) equconstants(i) = 10 ** enthalpies(i)
							endif
						endif
					endif
					if (array(2) .ne. '*') then
						e = array(2)
						if (e(1:5) .eq. 'ckeq=') then
							e = e(6:len_trim(e))
							if (e .eq. '*') then
								equconstants(i) = 0
							else
								read(e, *, err=263, end=263) equconstants(i)
								if (logten) equconstants(i) = 10 ** equconstants(i)
							endif
						elseif (e(1:4) .eq. 'heq=') then
							e = e(5:len_trim(e))
							if (e .eq. '*') then
								enthalpies(i) = 0
							else
								read(e, *, err=263, end=263) enthalpies(i)
								if (logten) equconstants(i) = 10 ** enthalpies(i)
							endif
						endif
					endif
					goto 284
				endif
				do j = 1, SPEC_MAX
					if (array(j) .eq. '*') exit
					if (array(j) .eq. '+') a = a + 1
				enddo
				read(line, *, err=262, end=239) complexes(i), null, (complexstoich(i, j), &
					complexcontents(i, j), null, j = 1, a - 1), complexstoich(i, a), &
					complexcontents(i, a), e, f
239				if (len_trim(e) .eq. 0) then
					equconstants(i) = 0
					enthalpies(i) = 0
				elseif (e(1:5) .eq. 'ckeq=') then
					e = e(6:len_trim(e))
					if (e .eq. '*') then
						equconstants(i) = 0
					else
						read(e, *, err=263, end=263) equconstants(i)
						if (logten) equconstants(i) = 10 ** equconstants(i)
					endif
					if (f(1:4) .eq. 'heq=') then
						f = f(5:len_trim(f))
						if (f .eq. '*') then
							enthalpies(i) = 0
						else
							read(f, *, err=263, end=263) enthalpies(i)
							if (logten) enthalpies(i) = 10 ** enthalpies(i)
						endif
					elseif ((len_trim(f) .eq. 0) .or. (f .eq. '*')) then
						enthalpies(i) = 0
					else
						goto 263
					endif
				elseif (e(1:4) .eq. 'heq=') then
					e = e(5:len_trim(e))
					if (e .eq. '*') then
						enthalpies(i) = 0
					else
						read(e, *, err=263, end=263) enthalpies(i)
						if (logten) enthalpies(i) = 10 ** enthalpies(i)
					endif
					if (f(1:5) .eq. 'ckeq=') then
						f = f(6:len_trim(f))
						if (f .eq. '*') then
							equconstants(i) = 0
						else 
							read(f, *, err=263, end=263) equconstants(i)
							if (logten) equconstants(i) = 10 ** equconstants(i)
						endif
					elseif ((len_trim(f) .eq. 0) .or. (f .eq. '*')) then
						equconstants(i) = 0
					else
						goto 263
					endif
				else
					goto 263
				endif
284				continue
			enddo
			if (debug) write(iptty, '(a, i3, a)') 'Read ', ncomplexes, ' complexes.'
			goto 2000
262			write(ierr, '(a)') "Error:  Bad reaction format in cplx."
			write(ierr,*) 'Make sure there are ',6+3*(a-1),' words in this line'
			goto 666
263			write(ierr, '(a)') "Error:  Bad parameter after equation in cplx."
			goto 666
		! Equilibrium constant model mode
		elseif(keyword .eq. 'equi') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading equilibrium constants...  '
			if (fequi) then
				write(ierr, '(a)') 'Warning:  equi has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(ecomplexes)
				deallocate(eqtemps)
				deallocate(eqconsts)
			endif
			fequi = .true.
			logten = .false.
			e = '*'
			read(line, *, end=279, err=600) keyword, e
279			if (e(1:3) .eq. 'log') logten = .true.
			necomplexes = 0
			neqtemps = 0
			i = 0
			j = 0
			do
				call readline
				if (len_trim(line) .eq. 0) then
					neqtemps = max(i, neqtemps)
					exit
				elseif (line(1:1) .eq. '.') then
					if ((j .ne. 0) .and. (i .eq. 0)) then
						write(ierr, '(a)') 'Error:  Each complex must have at least '// &
							'one temperature entry.'
						goto 666
					endif
					neqtemps = max(i, neqtemps)
					i = 0
					necomplexes = necomplexes + 1
				else
					if (j .eq. 0) then
						write(ierr, '(a)') 'Error:  The first item in equi must be a '// &
							'species name preceded by a period.'
						goto 666
					endif
					i = i + 1
				endif
				j = j + 1
			enddo
			do i = 1, j + 1
				backspace inpt
				linen = linen - 1
			enddo
			neqtemps = max(neqtemps, 3) ! LSTSQ doesn't like less than 3 parameters.  Thus:  If we are given one parameter, duplicate it for params 2 and 3.  If we are given 2 params, interpolate for param 3.
			allocate(ecomplexes(necomplexes))
			allocate(eqtemps(necomplexes, neqtemps))
			allocate(eqconsts(necomplexes, neqtemps))
			eqtemps = 0
			eqconsts = -1
			i = 0	! Counter for temps
			j = 0	! Counter for species
			do
				call readline
				if (len_trim(line) .eq. 0) then
					exit
				elseif (line(1:1) .eq. '.') then
					if (j .ne. 0) then
						if (i .le. 0) then
							write(ierr, '(a)') 'Internal error processing equilibrium constants.'
							goto 666
						elseif (i .eq. 1) then
							eqtemps(j, 2) = eqtemps(j, 1) + 1
							eqconsts(j, 2) = eqconsts(j, 1)
							eqtemps(j, 3) = eqtemps(j, 1) + 2
							eqconsts(j, 3) = eqconsts(j, 1)
						elseif (i .eq. 2) then
							if (eqtemps(j, 1) .eq. eqtemps(j, 2)) then
								eqtemps(j, 3) = eqtemps(j, 1) + 1
							else
								eqtemps(j, 3) = (eqtemps(j, 1) + eqtemps(j, 2)) / 2
							endif
							eqconsts(j, 3) = (eqconsts(j, 1) + eqconsts(j, 2)) / 2
						endif
					endif
					i = 0
					j = j + 1
					ecomplexes(j) = line(2:len_trim(line))
				else
					i = i + 1
					read(line, *, err=274, end=274) eqtemps(j, i), eqconsts(j, i)
					if (logten) eqconsts(j, i) = 10 ** eqconsts(j, i)
					if (eqconsts(j, i) .lt. 0) then
						write(ierr, '(a)') 'Error:  Equilibrium constant cannot be less than zero.'
						goto 666
					endif
				endif
			enddo
			if (debug) write(iptty, '(a, i3, a)') 'Read constants for ', necomplexes, ' species.'
			goto 2000
274			write(ierr, '(a)') 'Error reading equilibrium constants in equi.'
			goto 666
		! Distribution coefficient model mode
		elseif (keyword .eq. 'dist') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading distribution models...  '
			if (fdist) then
				write(ierr, '(a)') 'Warning:  dist has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(distmodels)
				deallocate(distmodelnames)
			endif
			fdist = .true.
			ndistmodels = 0
			i = 0
			j = 0
			do k = 1, SPEC_MAX
				call readline
				if (len_trim(line) .eq. 0) exit
				if (line(1:1) .eq. '.') then
					if ((j .eq. 0) .and. (ndistmodels .ne. 0)) then
						write(ierr, '(a)') 'Error:  Each distribution model must have '// &
							'at least one temperature entry.'
						goto 666
					endif
					ndistmodels = ndistmodels + 1
					if (j .gt. i) then
						i = j
					endif
					j = 0
				else
					j = j + 1
				endif
			enddo
			k = k - 1
			i = max(i, 3)
			allocate(distmodels(ndistmodels, i, 2))
			allocate(distmodelnames(ndistmodels))
			distmodels = -1
			do i = 1, k + 1
				backspace inpt
				linen = linen - 1
			enddo
			j = 0
			o = -1	! Throw an error if the user doesn't start with a model name
			do i = 1, k
				call readline
				if (line(1:1) .eq. '.') then
					if (j .ne. 0) then
						if (o .le. 0) then
							write(ierr, '(a)') 'Internal error processing distribution constants.'
							goto 667
						elseif (o .eq. 1) then
							distmodels(j, 2, 1) = distmodels(j, 1, 1) + 1
							distmodels(j, 2, 2) = distmodels(j, 1, 2)
							distmodels(j, 3, 1) = distmodels(j, 3, 1) + 2
							distmodels(j, 3, 2) = distmodels(j, 3, 2)
						elseif (o .eq. 2) then
							if (distmodels(j, 1, 1) .eq. distmodels(j, 2, 1)) then
								distmodels(j, 3, 1) = distmodels(j, 1, 1) + 1
							else
								distmodels(j, 3, 1) = (distmodels(j, 1, 1) + &
									distmodels(j, 2, 1)) / 2
							endif
							distmodels(j, 3, 2) = (distmodels(j, 1, 2) + distmodels(j, 2, 2)) / 2
						endif
					endif
					j = j + 1
					distmodelnames(j) = line(2:len_trim(line))
					o = 0
				else
					if (o .lt. 0) then
						write(ierr, '(a)') 'Error:  The second line of dist must be a model name.'
						goto 666
					endif
					o = o + 1
					read(line, *, end=181, err=181) distmodels(j, o, 1), distmodels(j, o, 2)
					if (distmodels(j, o, 2) .lt. 0) then
						write(ierr, '(a)') 'Error:  Distribution coefficients in distribution '// &
							'model "'//distmodelnames(j)//'" cannot be negative.'
						goto 666
					endif
				endif
			enddo
			if (debug) write(iptty, '(a, i3, a)') 'Read ', ndistmodels, ' models.'
			goto 2000
181			write(ierr, '(a)') 'Error reading models in dist.'
			goto 666
		! Solubility coefficient model mode
		elseif (keyword .eq. 'sol') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading solubility models...  '
			if (fsol) then
				write(ierr, '(a)') 'Warning:  sol has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(solmodels)
				deallocate(solmodelnames)
			endif
			fsol = .true.
			nsolmodels = 0
			i = 0
			j = 0
			do k = 1, SPEC_MAX
				call readline
				if (len_trim(line) .eq. 0) then
					if (j .gt. i) then
						i = j
					endif
					exit
				endif
				if (line(1:1) .eq. '.') then
					if ((j .eq. 0) .and. (nsolmodels .ne. 0)) then
						write(ierr, '(a)') 'Error:  Each solubility model must have '// &
							'at least one temperature entry.'
						goto 666
					endif
					nsolmodels = nsolmodels + 1
					if (j .gt. i) then
						i = j
					endif
					j = 0
				else
					j = j + 1
				endif
			enddo
			k = k - 1
			i = max(i, 3)
			allocate(solmodels(nsolmodels, i, 2))
			allocate(solmodelnames(nsolmodels))
			solmodels = -1
			do i = 1, k + 1
				backspace inpt
				linen = linen - 1
			enddo
			j = 0
			o = -1	! Throw an error if the user doesn't start with a model name
			do i = 1, k
				call readline
				if (line(1:1) .eq. '.') then
					if (j .ne. 0) then
						if (o .le. 0) then
							write(ierr, '(a)') 'Internal error processing distribution constants.'
							goto 667
						elseif (o .eq. 1) then
							solmodels(j, 2, 1) = solmodels(j, 1, 1) + 1
							solmodels(j, 2, 2) = solmodels(j, 1, 2)
							solmodels(j, 3, 1) = solmodels(j, 3, 1) + 2
							solmodels(j, 3, 2) = solmodels(j, 3, 2)
						elseif (o .eq. 2) then
							if (solmodels(j, 1, 1) .eq. solmodels(j, 2, 1)) then
								solmodels(j, 3, 1) = solmodels(j, 1, 1) + 1
							else
								solmodels(j, 3, 1) = (solmodels(j, 1, 1) + &
									solmodels(j, 2, 1)) / 2
							endif
							solmodels(j, 3, 2) = (solmodels(j, 1, 2) + solmodels(j, 2, 2)) / 2
						endif
					endif
					j = j + 1
					solmodelnames(j) = line(2:len_trim(line))
					o = 0
				else
					if (o .lt. 0) then
						write(ierr, '(a)') 'Error:  The second line of sol must be a model name.'
						goto 666
					endif
					o = o + 1
					read(line, *, end=222, err=222) solmodels(j, o, 1), solmodels(j, o, 2)
					!if (solmodels(j, o, 2) .lt. 0) then
					!	write(ierr, '(a)') 'Error:  Solubility coefficients in solubility '// &
					!		'model "'//solmodelnames(j)//'" cannot be negative.'
					!	goto 666
					!endif
				endif
			enddo
			if (debug) write(iptty, '(a, i3, a8)') 'Read ', nsolmodels, ' models.'
			goto 2000
222			write(ierr, '(a)') 'Error reading models in sol.'
			goto 666
		! Reaction input mode
		elseif (keyword .eq. 'rxn') then
			nrxns = nrxns + 1
			if (debug) write(iptty, '(a, i3, a)', advance='no') 'Reading reaction ', nrxns, '...  '
			backspace inpt
			linen = linen - 1
			call readline
			read(line, *, err=182) keyword, rxntypes(nrxns)!, rxnnames(nrxns)
			if (rxntypes(nrxns) .eq. 1) then
				call readline
				array = '*'
				read(line, *, end=186, err=182) (array(i), i = 1, SPEC_MAX)
186				reactants(nrxns, 1) = array(1)
				if (index(array(2), '=') .eq. 0) then
					goto 183
				endif
				products(nrxns, 1) = array(3)
				stoichiometries(nrxns, 1, 1) = 1
				stoichiometries(nrxns, 2, 1) = 1
				rates(nrxns) = -1
				dscoefs(nrxns) = '*'
				call readline
				array = '*'
				read(line, *, end=233, err=600) (array(i), i = 1, SPEC_MAX)
233				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0) then
						e = array(i)(1:index(array(i), '=') - 1)
						f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:4) .eq. 'rate') then
							read(f, *, err=184) rates(nrxns)
						elseif (e(1:4) .eq. 'dist') then
							read(f, *, err=184) dscoefs(nrxns)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Rate constant negative or missing.'
					goto 666
				endif
				if (dscoefs(nrxns) .eq. '*') then
					write(ierr, '(a)') 'Error:  Distribution coefficient omitted or invalid.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read linear kinetic reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			elseif (rxntypes(nrxns) .eq. 2) then
				call readline
				array = '*'
				read(line, *, end=187, err=182) (array(i), i = 1, SPEC_MAX)
187				reactants(nrxns, 1) = array(1)
				if (index(array(2), '=') .eq. 0) then
					goto 183
				endif
				products(nrxns, 1) = array(3)
				stoichiometries(nrxns, 1, 1) = 1
				stoichiometries(nrxns, 2, 1) = 1
				rates(nrxns) = -1
				dscoefs(nrxns) = '*'
				xcoef(nrxns) = -1
				call readline
				array = '*'
				read(line, *, end=234, err=600) (array(i), i = 1, SPEC_MAX)
234				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0) then
						e = array(i)(1:index(array(i), '=') - 1)
						f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:4) .eq. 'rate') then
							read(f, *, err=184) rates(nrxns)
						elseif (e(1:4) .eq. 'dist') then
							read(f, *, err=184) dscoefs(nrxns)
						elseif (e(1:3) .eq. 'max') then
							read(f, *, err=184) xcoef(nrxns)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Rate constant negative or missing.'
					goto 666
				endif
				if (dscoefs(nrxns) .eq. '*') then
					write(ierr, '(a)') 'Error:  Distribution coefficient invalid or missing.'
					goto 666
				endif
				if (xcoef(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Maximum sorption concentration negative or missing.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read Langmuir kinetic reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			elseif (rxntypes(nrxns) .eq. 3) then
				call readline
				array = '*'
				read(line, *, end=188, err=182) (array(i), i = 1, SPEC_MAX)
188				flag = .false.
				i = 0
				j = 1
				do
					i = i + 1
					if (i .gt. SPEC_MAX) then
						exit
					elseif (array(i) .eq. '+') then
						continue
					elseif (index(array(i), '=') .ne. 0) then
						if (flag) then
							goto 183
						endif
						flag = .true.
						j = 1
					elseif (array(i) .eq. '+') then
						continue
					else
						if (.not. flag) then
							stoichiometries(nrxns, 1, j) = 1
							read(array(i), *, err=189) stoichiometries(nrxns, 1, j)
							i = i + 1
							goto 189
						else
							stoichiometries(nrxns, 2, j) = 1
							read(array(i), *, err=190) stoichiometries(nrxns, 2, j)
							i = i + 1
							goto 190
						endif
						goto 191
189						read(array(i), *, err=183) reactants(nrxns, j)
						j = j + 1
						goto 191
190						read(array(i), *, err=183) products(nrxns, j)
						j = j + 1
191						continue
					endif
				enddo
				rates(nrxns) = -1
				rates2(nrxns) = -1
				call readline
				array = '*'
				read(line, *, end=235, err=600) (array(i), i = 1, SPEC_MAX)
235				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0) then
						e = array(i)(1:index(array(i), '=') - 1)
						f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:3) .eq. 'for') then
							read(f, *, err=184) rates(nrxns)
						elseif (e(1:3) .eq. 'rev') then
							read(f, *, err=184) rates2(nrxns)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Forward rate constant negative or missing.'
					goto 666
				endif
				if (rates2(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Reverse rate constant negative or missing.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read general kinetic reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			elseif (rxntypes(nrxns) .eq. 4) then
				fsubs = .false.
				felec = .false.
				fbiomass = .false.
				do
					call readline
					if (index(line, ':') .eq. 0) exit
					e = line(1:index(line, ':') - 1)
					f = line(index(line, ':') + 1:len_trim(line))
					if (e(1:4) .eq. 'subs') then
						fsubs = .true.
						bioparams(nrxns, 1) = e
					elseif (e(1:4) .eq. 'elec') then
						felec = .true.
						bioparams(nrxns, 2) = e
					elseif (e(1:7) .eq. 'biomass') then
						fbiomass = .true.
						bioparams(nrxns, 3) = e
					elseif (e(1:4) .eq. 'reac') then
						array = '*'
						read(e, *, end=214, err=600) (array(j), j = 1, SPEC_MAX)
						k = 1
214						do j = 1, SPEC_MAX
							if (array(k) .eq. '*') then
								exit
							endif
							read(array(k), *, err=215) stoichiometries(i, 1, j)
							k = k + 1
							if (array(k) .eq. '*') then
								goto 183
							endif
							goto 220
215							stoichiometries(i, 1, j) = 1
220							read(array(k), *, err=183, end=183) reactants(i, j)
							k = k + 1
						enddo
					elseif (e(1:4) .eq. 'prod') then
						array = '*'
						read(e, *, end=216, err=600) (array(j), j = 1, SPEC_MAX)
						k = 1
216						do j = 1, SPEC_MAX
							if (array(k) .eq. '*') then
								exit
							endif
							read(array(k), *, err=217) stoichiometries(i, 2, j)
							k = k + 1
							if (array(k) .eq. '*') then
								goto 183
							endif
							goto 221
217							stoichiometries(i, 2, j) = 1
221							read(array(k), *, err=183, end=183) products(i, j)
							k = k + 1
						enddo
					elseif (e(1:6) .eq. 'biodeg') then
						read(f, *, end=192, err=184) (icbioholder(nrxns, i), i = 1, SPEC_MAX)
					else
						goto 185
					endif
192					continue
				enddo
				if (.not. fsubs) then
					write(ierr, '(a, i3, a)') 'Error:  No substrate provided for type 4 '// &
						'reaction ', nrxns, '.'
					goto 666
				endif
				if (.not. felec) then
					write(ierr, '(a, i3, a)') 'Error:  No electron acceptor provided for type 4 '// &
						'reaction ', nrxns, '.'
					goto 666
				endif
				if (.not. fbiomass) then
					write(ierr, '(a, i3, a)') 'Error:  No biomass provided for type 4 '// &
						'reaction ', nrxns, '.'
					goto 666
				endif
				bioparams(nrxns, 4) = '*'
				bioparams(nrxns, 5) = '*'
				bioparams(nrxns, 6) = '*'
				bioparams(nrxns, 7) = '*'
				bioparams(nrxns, 8) = '*'
				rates(nrxns) = 0
				rates2(nrxns) = 0
				call readline
				array = '*'
				read(line, *, end=236, err=600) (array(i), i = 1, SPEC_MAX)
236				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0) then
						e = array(i)(1:index(array(i), '=') - 1)
						f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:2) .eq. 'ks') then
							read(f, *, err=184) bioparams(nrxns, 4)
						elseif (e(1:2) .eq. 'ka') then
							read(f, *, err=184) bioparams(nrxns, 5)
						elseif (e(1:3) .eq. 'dec') then
							read(f, *, err=184) rates2(nrxns)
						elseif (e(1:2) .eq. 'ph') then
							read(f, *, err=184) bioparams(nrxns, 6)
						elseif (e(1:2) .eq. 'qm') then
							read(f, *, err=184) rates(nrxns)
						elseif (e(1:5) .eq. 'yield') then
							read(f, *, err=184) bioparams(nrxns, 7)
						elseif (e(1:2) .eq. 'xm') then
							read(f, *, err=184) bioparams(nrxns, 8)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (bioparams(nrxns, 4) .eq. '*') then
					write(ierr, '(a)') 'Error:  Substrate half maximum rate concentration '// &
						'invalid or missing.'
					goto 666
				endif
				if (bioparams(nrxns, 5) .eq. '*') then
					write(ierr, '(a)') 'Error:  Electron acceptor half maximum rate '// &
						'concentration invalid or missing.'
					goto 666
				endif
				if (bioparams(nrxns, 6) .eq. '*') then
					write(ierr, '(a)') 'Error:  Biodegradation pH threshold invalid or missing.'
					goto 666
				endif
				if (bioparams(nrxns, 7) .eq. '*') then
					write(ierr, '(a)') 'Error:  Microbial yield coefficient invalid or missing.'
					goto 666
				endif
				if (bioparams(nrxns, 8) .eq. '*') then
					write(ierr, '(a)') 'Error:  Minimum biomass concentration invalid or missing.'
					goto 666
				endif
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Maximum substrate utilization rate negative or missing.'
					goto 666
				endif
				if (rates2(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Microbial decay rate negative or missing.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read biodegradation reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			elseif (rxntypes(nrxns) .eq. 5) then
				call readline
				array = '*'
				read(line, *, end=193, err=600) (array(i), i = 1, SPEC_MAX)
193				reactants(nrxns, 1) = array(1)
				if (index(array(2), '=') .eq. 0) then
					goto 183
				endif
				products(nrxns, 1) = array(3)
				stoichiometries(nrxns, 1, 1) = 1
				stoichiometries(nrxns, 2, 1) = 1
				rates(nrxns) = -1
				call readline
				array = '*'
				read(line, *, end=237, err=600) (array(i), i = 1, SPEC_MAX)
237				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0)then
						e = array(1)(1:index(array(i), '=') - 1)
						f = array(1)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:4) .eq. 'half') then
							read(f, *, err=184) rates(nrxns)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Halflife negative or missing.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read radioactive decay reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			elseif (rxntypes(nrxns) .eq. 6) then
				write(ierr, '(a)') 'Error: reaction type 6 not supported by trxn.'
				goto 666
			elseif ((rxntypes(nrxns) .eq. 7) .or.(rxntypes(nrxns) .eq. 8)) then
				call readline
				array = '*'
				dscoefs(nrxns) = '*'
				read(line, *, end=194, err=182) (array(i), i = 1, SPEC_MAX)
194				if (array(2) .eq. '*') then ! Lookup option
					if (.not. flookup) then
						write(ierr, '(a)') 'Error:  A lookup block must be provided '// &
							'before reaction blocks if the dynamic lookup option is to be used.'
						goto 666
					endif
					do i = 1, ndatmin
						if ((array(1) .eq. datminnam(i)) .or. (array(1) .eq. datminmain(1))) then
							reactants(nrxns, 1) = datminmain(i)
							stoichiometries(nrxns, 1, 1) = 1
							do j = 1, maxdatmin
								if (datmin(i, j) .eq. '*') exit
								products(nrxns, j) = datmin(i, j)
								stoichiometries(nrxns, 2, j) = datminstoic(i, j)
							enddo
							write(dscoefs(nrxns), *) datminlkeq(i)
							goto 285
						endif
					enddo
					write(ierr, '(a)') 'Error:  Mineral "'//trim(reactants(nrxns, 1))// &
						'" not found in lookup database.'
					goto 666
				endif
				flag = .false.
				i = 0
				j = 1
				do
					i = i + 1
					if (i .gt. SPEC_MAX) then
						exit
					elseif (array(i) .eq. '+') then
						continue
					elseif (index(array(i), '=') .ne. 0) then
						if (flag) then
							goto 183
						endif
						flag = .true.
						j = 1
					elseif (array(i) .eq. '*') then
						exit
					else
						if (.not. flag) then
							stoichiometries(nrxns, 1, j) = 1
							read(array(i), *, err=195) stoichiometries(nrxns, 1, j)
							i = i + 1
							goto 195
						else
							stoichiometries(nrxns, 2, j) = 1
							read(array(i), *, err=196) stoichiometries(nrxns, 2, j)
							i = i + 1
							goto 196
						endif
						goto 197
195						read(array(i), *, err=183) reactants(nrxns, j)
						j = j + 1
						goto 197
196						read(array(i), *, err=183) products(nrxns, j)
						j = j + 1
197						continue
					endif
				enddo
285				call readline
				array = '*'
				rates(nrxns) = -1
				xcoef(nrxns) = -1
				porchange(nrxns, 1) = -1
				porchange(nrxns, 2) = -1
				read(line, *, end=238, err=600) (array(i), i = 1, SPEC_MAX)
238				do i = 1, SPEC_MAX
					if (array(i) .eq. '*') then
						exit
					endif
					if (index(array(i), '=') .ne. 0) then
						e = array(i)(1:index(array(i), '=') - 1)
						f = array(i)(index(array(i), '=') + 1:len_trim(array(i)))
						if (e(1:3) .eq. 'sol') then
							read(f, *, err=184) dscoefs(nrxns)
						elseif (e(1:4) .eq. 'rate') then
							read(f, *, err=184) rates(nrxns)
						elseif (e(1:5) .eq. 'sarea') then
							read(f, *, err=184) xcoef(nrxns)
						elseif (e(1:3) .eq. 'mol') then
							if (rxntypes(nrxns) .eq. 7) then
								write(ierr, '(a)') 'Error:  Molecular weight only '// &
									'available for reaction type 8.'
								goto 666
							endif
							read(f, *, err=184) porchange(nrxns, 1)
						elseif (e(1:4) .eq. 'dens') then
							if (rxntypes(nrxns) .eq. 7) then
								write(ierr, '(a)') 'Error:  Mineral density only '// &
									'available for reaction type 8.'
								goto 666
							endif
							read(f, *, err=184) porchange(nrxns, 2)
						else
							goto 185
						endif
					else
						goto 184
					endif
				enddo
				if (dscoefs(nrxns) .eq. '*') then
					write(ierr, '(a)') 'Error:  Solubility product invalid or missing.'
					goto 666
				endif
				if (rates(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Reaction rate negative or missing.'
					goto 666
				endif
				if (xcoef(nrxns) .lt. 0) then
					write(ierr, '(a)') 'Error:  Mineral surface area negative or missing.'
					goto 666
				endif
				if ((porchange(nrxns, 1) .lt. 0) .and. (rxntypes(nrxns) .eq. 8)) then
					write(ierr, '(a)') 'Error:  Mineral molecular weight negative or missing.'
					goto 666
				endif
				if ((porchange(nrxns, 2) .lt. 0) .and. (rxntypes(nrxns) .eq. 8)) then
					write(ierr, '(a)') 'Error:  Mineral density negative or missing.'
					goto 666
				endif
				if (debug) write(iptty, '(a)') 'Read precipitation/dissolution reaction.' !', rxnnames(nrxns)(1:len_trim(rxnnames(nrxns))), '".'
				goto 2000
			else
				write(ierr, '(a, i3, a, i3)') 'Error:  reaction type ', rxntypes(nrxns), &
					' not known for reaction ', nrxns, '.'
				goto 666
			endif
182			write(ierr, '(a)') 'Error reading reaction.'
			goto 666
183			write(ierr, '(a)') 'Error:  Invalid reaction format.'
			goto 666
184			write(ierr, '(a)') 'Error:  Invalid syntax in key/value list.'
			goto 666
185			write(ierr, '(a)') 'Error:  Unknown key "'//trim(e)//'" in key/value list.'
			goto 666
		! Assignment mode
		elseif (keyword .eq. 'assign') then
			if (debug) write(iptty, '(a)', advance='no') 'Reading zone attribute assignment data...  '
			if (fassign) then
				write(ierr, '(a)') 'Warning:  assign has already been specified.  '// &
					'The new values will overwrite the old ones.'
				deallocate(zones)
				deallocate(zonegrid)
			endif
			fassign = .true.
			array = '*'
			zonespecs = '*'
			read(line, *, end=117, err=144) keyword, (array(i), i = 1, SPEC_MAX)
117			do i = 1, SPEC_MAX
				if (array(i) .eq. '*') then
					nzonespecs = i - 1
					exit
				endif
			enddo
			do i = 1, nzonespecs
				zonespecs(i) = array(i)
			enddo
			nzones = 0
			do
				call readline
				if (len_trim(line) .eq. 0) exit
				nzones = nzones + 1
			enddo
			do i = 1, nzones + 1
				backspace inpt
				linen = linen - 1
			enddo
			allocate(zones(max(numzones, nzones) + 1))
			allocate(zonegrid(max(numzones, nzones) + 1, SPEC_MAX))
			
			zones = '*'
			zonegrid = '*'
			do i = 1, nzones
				call readline
				read(line, *, end=144, err=144) zones(i), (zonegrid(i, j), j = 1, nzonespecs)
			enddo
			
			
			
			
			
			goto 143
144			write(ierr, '(a)') 'Error reading zone data.'
			goto 666
143			if (debug) write(iptty, '(a, i3, a, i3, a)', advance='no') 'Read ', nzonespecs, &
				' parameter(s) for ', nzones, ' zone(s)'
			if (nzonespecs .eq. 0) then
				if (debug) write(iptty, '(a)') '.'
			else
				if (debug) write(iptty, '(a)', advance='no') ':  '
				if (debug) write(iptty, '('//aSPEC_MAX//'a20)') (trim(zones(i)), i = 1, nzones)
			endif
			goto 2000
		! "end," in which case we know that we are all done.
		elseif (keyword(1:3) .eq. 'end') then
			if (libread .gt. 0) goto 801
			reading = .false.
			goto 4000
		! "null," which we ignore.
		elseif (keyword .eq. 'null') then
			if (debug) then
				e = '*'
				read(line, *, end=256, err=600) keyword, e
256				if (e .ne. '*') write(iptty, '(a)') 'Skipping null statement "'//trim(e)//'".'
			endif
			do
				call readline
				if (len_trim(line) .eq. 0) exit
			enddo
			goto 2000
		! a nonsense line, in which case an error message is printed and an error status returned to the calling program.
		else
			write(ierr, '(a)') 'Error:  keyword "'//trim(keyword)//'" not recognized.'
			goto 666
		endif
		! Wrap up the processing of each line
2000		continue
	enddo
4000	if (reading) then
		write(ierr, '(a)') 'Warning:  Missing "end trxn" keyword may cause unexpected behavior.'
		reading = .false.
	endif
	! Process the data and put them in a format that FEHM will like.
	!if (debug) write(iptty, '(i3, a)') linen, ' lines read.'
	if (debug) write(iptty, '(a)') 'Checking data...'
	! Sanity checks!
	! Ensure that the bare minimums (minima?) are present
	if (.not. fcomp) then
		write(ierr, '(a)') 'Error:  comp was not specified.'
		goto 666
	endif
	if (ncp .lt. 1) then
		write(ierr, '(a)') 'Error:  At least one component must be specified in comp.'
		goto 666
	endif
	if (.not. fheader) then
		write(ierr, '(a)') 'Error:  header was not specified.'
		goto 666
	endif
	! Check that selected set in trxninit match the correct values (more checks are done when setting the values)
	if (nspeci .ne. nspecies) then
		write(ierr, '(a, i3, a, i3, a)') 'Internal error:  Components do not pass presetting check.  (', &
			nspeci, '/', nspecies, ')'
		goto 667
	endif
	if (numrxn .ne. nrxns) then
		write(ierr, '(a, i3, a, i3, a)') 'Internal error:  Reactions do not pass presetting check.  (', &
			numrxn, '/', nrxns, ')'
		goto 667
	endif
	if (ncplx .ne. ncomplexes) then
		write(ierr, '(a, i3, a, i3, a)') 'Internal error:  Complexes do not pass presetting check.  (', &
			ncplx, '/', ncomplexes, ')'
		goto 667
	endif
	! Check that all the columns in assign are good and fill in values for missing columns.
	zinit = 0
	zboun = 0
	zrock = 0
	zdisp = 0
	zsorp = 0
	ztpor = 0
	ztime = 0
	zvap = 0
	zhparam = 0
	zopt = 0
	do i = 1, nzonespecs
		if (zonespecs(i) .eq. 'water') then
			zinit = i
		elseif (zonespecs(i) .eq. 'boun') then
			zboun = i
		elseif (zonespecs(i) .eq. 'rock') then
			zrock = i
		elseif (zonespecs(i) .eq. 'disp') then
			zdisp = i
		elseif (zonespecs(i) .eq. 'sorp') then
			zsorp = i
		elseif (zonespecs(i) .eq. 'tpor') then
			ztpor = i
		elseif (zonespecs(i) .eq. 'time') then
			ztime = i
		elseif (zonespecs(i) .eq. 'gas') then
			zvap = i
		elseif (zonespecs(i) .eq. 'opt') then
			zopt = i
		else
			write(ierr, '(a)') 'Error:  Column "'//trim(zonespecs(i))//'" in assign not recognized.'
			goto 666
		endif
	enddo
	! Ensure that the appropriate section has been included for each column specified in assign.
	if ((zboun .ne. 0) .and. (.not. fwater) .and. (.not. fvap)) then
		write(ierr, '(a)') 'Error:   boun column specified in assign, but no water or gas section specified.'
		goto 666
	elseif ((zinit .ne. 0) .and. (.not. fwater)) then
		write(ierr, '(a)') 'Error:  water column specified in assign, but no water section specified.'
		goto 666
	elseif ((zrock .ne. 0) .and. (.not. frock)) then
		write(ierr, '(a)') 'Error:  rock column specified in assign, but no rock section specified.'
		goto 666
	elseif ((zvap .ne. 0) .and. (.not. fvap)) then
		write(ierr, '(a)') 'Error:  gas column specified in assign, but no gas section specified.'
		goto 666
	elseif ((zsorp .ne. 0) .and. (.not. fsorp)) then
		write(ierr, '(a)') 'Error:  sorp column specified in assign, but no sorp section specified.'
		goto 666
	elseif ((zdisp .ne. 0) .and. (.not. fdisp)) then
		write(ierr, '(a)') 'Error:  disp column specified in assign, but no xyzdisp or londisp section specified.'
		goto 666
	elseif ((zinit + zrock + zvap .ne. 0) .and. (fmoles)) then
		write(ierr, '(a)') 'Error:  The water, rock, and gas columns must be omitted if the moles block is used.'
		goto 666
	endif
	if (zinit .eq. 0) then
		nzonespecs = nzonespecs + 1
		zinit = nzonespecs
		zonespecs(zinit) = 'water'
		do i = 1, nzones
			zonegrid(i, zinit) = '*'
		enddo
	endif
	if (zboun .eq. 0) then
		nzonespecs = nzonespecs + 1
		zboun = nzonespecs
		zonespecs(zboun) = 'boun'
		do i = 1, nzones
			zonegrid(i, zboun) = '*'
		enddo
	endif
	if (zrock .eq. 0) then
		nzonespecs = nzonespecs + 1
		zrock = nzonespecs
		zonespecs(zrock) = 'rock'
		do i = 1, nzones
			zonegrid(i, zrock) = '*'
		enddo
	endif
	if (zdisp .eq. 0) then
		nzonespecs = nzonespecs + 1
		zdisp = nzonespecs
		zonespecs(zdisp) = 'disp'
		do i = 1, nzones
			zonegrid(i, zdisp) = '*'
		enddo
	endif
	if (zsorp .eq. 0) then
		nzonespecs = nzonespecs + 1
		zsorp = nzonespecs
		zonespecs(zsorp) = 'sorp'
		do i = 1, nzones
			zonegrid(i, zsorp) = '*'
		enddo
	endif
	if (ztime .eq. 0) then
		nzonespecs = nzonespecs + 1
		ztime = nzonespecs
		zonespecs(ztime) = 'time'
		do i = 1, nzones
			zonegrid(i, ztime) = '*'
		enddo
	endif
	if (zvap .eq. 0) then
		nzonespecs = nzonespecs + 1
		zvap = nzonespecs
		zonespecs(zvap) = 'gas'
		do i = 1, nzones
			zonegrid(i, zvap) = '*'
		enddo
	endif
	if (ztpor .eq. 0) then
		nzonespecs = nzonespecs + 1
		ztpor = nzonespecs
		zonespecs(ztpor) = 'tpor'
		do i = 1, nzones
			zonegrid(i, ztpor) = '0.32'
		enddo
	endif
	if (zopt .eq. 0) then
		nzonespecs = nzonespecs + 1
		zopt = nzonespecs
		zonespecs(zopt) = 'opt'
		do i = 1, nzones
			zonegrid(i, zopt) = '*'
		enddo
	endif
	! Check zone names, set up zone number resolution, and assign default properties to extra zones.
	zoneresolv = 0
! first check that every zone in the assign macro was predefined, if not, report warning
	do i = 1, nzones   ! this the number of zones defined in the assign macro
		do j = 1, numzones  ! this is the number of zones defined before trxn macro
			if (zonenames(j) .eq. zones(i)) then  ! zonenames are predefined, zones are defined in assign macro
				zoneresolv(zonenums(j)) = i
				goto 259
			endif
		end do
		write(ierr, '(a)') 'Warning:  Zone "'//trim(zones(i))//'" in assign macro not '// &
					'included defined as a zone before trxn call'
 259 end do
! now check that every node was assigned something in the assign block
		do k=1,n0
			do i=1,nzones
				if(zonenames(izonef(k)).eq.zones(i)) goto 2164
			end do
		write(ierr, *) 'Error:  Node ', k ,' belongs to ',izonef(k),'**',trim(zonenames(izonef(k))),' ,not present in assign macro'
!# this logic expects zones to be named 1,2,3,4 	gaz 033121 (put '!' in front of '#'	)
		write(ierr,*) (zones(i),i=1,nzones)
		goto 666   
! would be good to assign this node default values and move on		
2164		end do
				
	! Ensure that for each name given in zone, there is a matching name in the appropriate section.
	
	do i = 1, nzonespecs
		if (zonespecs(i) .eq. 'water') then
			do j = 1, nzones
				if (zonegrid(j, i) .eq. '*') then
					goto 114
				endif
				if (nwt .eq. 0) then
					write(ierr, '(a)') 'Error:  Water type "'//trim(zonegrid(j, i))// &
						'" specified in assign, but no water types defined in water.'
					goto 666
				endif
				do k = 1, nwt
					if (zonegrid(j, i) .eq. wtnames(k)) then
						goto 114
					endif
				enddo
				write(ierr, '(a)') 'Error:  Water type "'//trim(zonegrid(j, i))// &
					'" specified in assign but not in water.'
				goto 666
114				continue
			enddo
		elseif (zonespecs(i) .eq. 'boun') then
			do j = 1, nzones
				if (zonegrid(j, i) .eq. '*') then
					goto 176
				endif
				e = zonegrid(j, i)//'.'
				do while(index(e, '.') .ne. 0)
					f = e(1:index(e, '.') - 1)
					e = e(index(e, '.') + 1:len(e))
					do k = 1, nwt
						if (f .eq. wtnames(k)) then
							goto 175
						endif
					enddo
					do k = 1, nrt
						if (f .eq. rtnames(k)) then
							goto 175
						endif
					enddo
					do k = 1, nvt
						if (f .eq. vtnames(k)) then
							goto 175
						endif
					enddo
					write(ierr, '(a)') 'Error:  Type "'//trim(f)// &
						'" specified in assign but not in water, rock, or gas.'
					goto 666
175					continue
				enddo
176				continue
			enddo
		elseif (zonespecs(i) .eq. 'rock') then
			do j = 1, nzones
				if (zonegrid(j, i) .eq. '*') then
					goto 115
				endif
				if (nrt .eq. 0) then
					write(ierr, '(a)') 'Error:  Rock type "'//trim(zonegrid(j, i))// &
						'" specified in assign, but no rock types defined in rock.'
					goto 666
				endif
				do k = 1, nrt
					if ((zonegrid(j, i) .eq. rtnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
						goto 115
					endif
				enddo
				write(ierr, '(a)') 'Error:  Rock type "'//trim(zonegrid(j, i))// &
					'" specified in assign but not in rock.'
				goto 666
115				continue
			enddo
		elseif (zonespecs(i) .eq. 'gas') then
			do j = 1, nzones
				if (zonegrid(j, i) .eq. '*') then
					goto 149
				endif
				if (nvt .eq. 0) then
					write(ierr, '(a)') 'Error:  Gas type "'//trim(zonegrid(j, i))// &
						'" specified in assign, but no gas types defined in gas.'
					goto 666
				endif
				do k = 1, nvt
					if ((zonegrid(j, i) .eq. vtnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
						goto 149
					endif
				enddo
				write(ierr, '(a)') 'Error:  Gas type "'//trim(zonegrid(j, i))// &
						'" specified in assign but not in gas.'
				goto 666
149				continue
			enddo
		elseif (zonespecs(i) .eq. 'sorp') then
			do j = 1, nzones
				if ((nsorp .eq. 0) .and. (zonegrid(j, i) .eq. '*')) then
					goto 136
				endif
				do k = 1, nsorp
					if ((zonegrid(j, i) .eq. sorpnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
						goto 136
					endif
				enddo
				write(ierr, '(a)') 'Error:  Adsorption type "'//trim(zonegrid(j, i))// &
					'" specified in assign but not in sorp.'
				goto 666
136				continue
			enddo
		elseif (zonespecs(i) .eq. 'disp') then
			do j = 1, nzones
				if ((ndisp .eq. 0) .and. (zonegrid(j, i) .eq. '*')) then
					goto 125
				endif
				do k = 1, ndisp
					if ((zonegrid(j, i) .eq. dispnames(k)) .or. (zonegrid(j, i) .eq. '*')) then
						goto 125
					endif
				enddo
				write(ierr, '(a)') 'Error:  Dispersivity model "'//trim(zonegrid(j, i))// &
					'" specified in assign but not in disp.'
					goto 666
125				continue
			enddo
		endif
	enddo
	! If H+ concentration is given in terms of pH, convert to [H+]
	fph = .false.
	do i = 1, nwtspecies
		if (wtspecies(i) .eq. 'pH') then
			fph = .true.
			wtspecies(i) = 'H' ! H is component, H+ is master species
			!do j = 1, nwt
			!	wtgrid(j, i) = 10 ** (-1.0 * wtgrid(j, i))
			!enddo
		endif
	enddo
	! Add an extra adsorption model with all parameters 0 for unspecified zones.
	do i = 1, anspecies
		lsorptypes(nsorp + 1, i) = 0
		vsorptypes(nsorp + 1, i) = 0
		a1l(nsorp + 1, i) = 0
		a2l(nsorp + 1, i) = 0
		bl(nsorp + 1, i) = 1
		a1v(nsorp + 1, i) = 0
		a2v(nsorp + 1, i) = 0
		bv(nsorp + 1, i) = 1
		sorpnames(nsorp + 1) = '*'
	enddo
	! Add an extra dispersion model with all parameters 0 for unspecified zones.
	lx(ndisp + 1) = 0
	ly(ndisp + 1) = 0
	lz(ndisp + 1) = 0
	ll(ndisp + 1) = 0
	lt(ndisp + 1) = 0
	vx(ndisp + 1) = 0
	vy(ndisp + 1) = 0
	vz(ndisp + 1) = 0
	vl(ndisp + 1) = 0
	vt(ndisp + 1) = 0
	dispnames(ndisp+1) = '*'
	! Check that all components in water/rock/gas/moles appear in comp and are of the appropriate state.
	if (fwater) then
		do i = 1, nwtspecies
			do j = 1, nspecies
				if (species(j) .eq. wtspecies(i)) then
					if ((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
						write(ierr, '(a)') 'Component "'//trim(wtspecies(i))// &
							'" in water is not aqueous or Henry''s Law.'
						goto 666
					else
						goto 169
					endif
				endif
			enddo
			write(ierr, '(a)') 'Component "'//trim(wtspecies(i))//'" specified in water but not in comp.'
			if (wtspecies(i) .eq. 'H') then
				write(ierr, '(a)') 'Remember that if you use "pH" as a column header in water, '// &
					'"H" must be defined in comp.'
			endif
			goto 666
169			continue
		enddo
	endif
	
	if (frock) then
		do i = 1, nrtspecies
			do j = 1, nspecies
				if (species(j) .eq. rtspecies(i)) then
					if (states(j) .ne. 's') then
						write(ierr, '(a)') 'Component "'//trim(rtspecies(i))// &
							'" in rock is not solid.'
						goto 666
					else
						goto 170
					endif
				endif
			enddo
			write(ierr, '(a)') 'Component "'//trim(rtspecies(i))//'" specified in rock but not in comp.'
			goto 666
170			continue
		enddo
	endif
	if (fvap) then
		do i = 1, nvtspecies
			do j = 1, nspecies
				if (species(j) .eq. vtspecies(i)) then
					if ((states(j) .ne. 'g') .and. (states(j) .ne. 'h')) then
						write(ierr, '(a)') 'Component "'//trim(vtspecies(i))// &
							'" in gas is not gaseous or Henry''s Law.'
						goto 666
					else
						goto 171
					endif
				endif
			enddo
			write(ierr, '(a)') 'Component "'//trim(vtspecies(i))//'" specified in gas but not in comp.'
			goto 666
171			continue
		enddo
	endif
	if (fmoles) then
		do i = 1, nmspecs
			do j = 1, nspecies
				if (species(j) .eq. mspecs(i)) then
					goto 335
				endif
			enddo
			write(ierr, '(a)') 'Component "'//trim(mspecs(i))//'" specified in moles but not in comp.'
			goto 666
335			continue
		enddo
	endif
	! Ensure that hparam exists if there are Henry's Law components, and that hparam contains all Henry's Law components.
	if (.not. fhparam) then
		do i = 1, nspecies
			if (states(i) .eq. 'h') then
				write(ierr, '(a)') 'Error:  Henry''s Law components have been specified in comp, '// &
					'but hparam is not present.'
				goto 666
			endif
		enddo
	else
		do i = 1, nspecies
			if (states(i) .eq. 'h') then
				do j = 1, hnhtspecies
					if (hhtspecies(j) .eq. species(i)) then
						goto 173
					endif
				enddo
				write(ierr, '(a)') 'Error:  Henry''s Law component "'//trim(species(i))// &
					'" does not appear in hparam.'
				goto 666
			endif
173		continue
		enddo
	endif
	! Ensure that everything in print is a component or master species.
	if (fprint) then
		do i = 1, nprint  ! at this point nprint is == 1
			if (printspecies(i) .eq. 'all') then
				nprint = nspecies + nmasters + ncomplexes
! print total component, free ion, 
				deallocate(printspecies)
				allocate(printspecies(nprint))
				do j = 1, nspecies
					printspecies(j) = species(j)  ! for each component
				enddo
				do j = 1, ncomplexes
					printspecies(nspecies + j) = complexes(j)  ! each complex
				enddo
				k = 1
				do j = 1, nspecies
					if (masters(j) .ne. '*') then
						printspecies(nspecies + ncomplexes + k) = masters(j)  ! each free ion
						k = k + 1
					endif
				enddo
				goto 202
			elseif (printspecies(i) .eq. 'none') then
				nprint = 0
				goto 202
			endif
			do j = 1, nspecies
				if ((printspecies(i) .eq. species(j)) .or. (printspecies(i) .eq. masters(j))) then
					goto 200
				endif
			enddo
			write(ierr, '(a)') 'Error:  Print specification "'//trim(printspecies(i))// &
				'" in print is not a component or master species.'
			goto 666
200			continue
		enddo
202		continue
		do j=1,nprint
		if(debug) write(*,*) printspecies(j), 'is going to print'
		end do
	else
		nprint = nspecies + nmasters + ncomplexes
!		write(*,*) '3866: ',nspecies,nmasters,ncomplexes (gaz 041921 debug line)
		allocate(printspecies(nprint))
		do j = 1, nspecies
			printspecies(j) = species(j)
		enddo
		do j = 1, ncomplexes
			printspecies(nspecies + j) = complexes(j)
		enddo
		k = 1
		do j = 1, nspecies
			if (masters(j) .ne. '*') then
				printspecies(nspecies + ncomplexes + k) = masters(j)
				k = k + 1
			endif
		enddo	
	endif
	! Ensure that all dscoefs are names of models in dist and sol.
	do i = 1, nrxns
		if (dscoefs(i) .eq. '*') then
			goto 203
		endif
		r = 0
		read(dscoefs(i), *, err=204) r
		goto 203
204		if ((rxntypes(i) .eq. 1) .or. (rxntypes(i) .eq. 2)) then
			do j = 1, ndistmodels
				if (dscoefs(i) .eq. distmodelnames(j)) then
					goto 203
				endif
			enddo
			write(ierr, '(a, i3, a)') 'Error:  Distribution model "'//trim(dscoefs(i))//'" in reaction ', i, &
				' not found in dist.'
			goto 666
		elseif ((rxntypes(i) .eq. 7) .or. (rxntypes(i) .eq. 8)) then
			do j = 1, nsolmodels
				if (dscoefs(i) .eq. solmodelnames(j)) then
					goto 203
				endif
			enddo
			write(ierr, '(a, i3, a)') 'Error:  Solubility model "'//trim(dscoefs(i))//'" in reaction ', i, &
				'" not found in sol.'
			goto 666
		endif
203		continue
	enddo
	! Check that if equi keyword is provided in cplx, an equi block appears.
	if ((cequi) .and. (.not. fequi)) then
		write(ierr, '(a)') 'Error:  "equi" keyword supplied in cplx, but no equi block provided.'
		goto 666
	endif
	! Check that all complexes in equi are in cplx and vice versa.
	if (cequi) then
		do i = 1, ncomplexes
			do j = 1, necomplexes
				if (complexes(i) .eq. ecomplexes(j)) goto 276
			enddo
			write(ierr, '(a)') 'Error:  Complex "'//trim(complexes(i))//'" in cplx does not appear in equi.'
			goto 666
276			continue
		enddo
		do i = 1, necomplexes
			do j = 1, ncomplexes
				if (ecomplexes(i) .eq. complexes(j)) goto 277
			enddo
			write(ierr, '(a)') 'Error:  Complex "'//trim(ecomplexes(i))//'" in equi does not appear in cplx.'
			goto 666
277			continue
		enddo
	endif
	! Create pseudo-states 'i' and 'j' for Henry's Law components that are being treated as gases and Henry's Law components whose primary state has not been specified, respectively.
	do i = 1, nspecies
		if (states(i) .eq. 'h') then
			a = 0
			c = 0
			do j = 1, nwtspecies
				if (species(i) .eq. wtspecies(j)) then
					a = j
					exit
				endif
			enddo
			do j = 1, nvtspecies
				if (species(i) .eq. vtspecies(j)) then
					c = j
					exit
				endif
			enddo
			if ((a .ne. 0) .and. (c .ne. 0)) then
				write (ierr, '(a)') 'Error:  Henry''s Law species "'//trim(species(i))// &
					'" cannot be included in both water and gas.'
				goto 666
			elseif ((a .eq. 0) .and. (c .eq. 0)) then
				states(i) = 'j'
			elseif (a .ne. 0) then
				states(i) = 'h'
			elseif (c .ne. 0) then
				states(i) = 'i'
			else
				write (ierr, '(a)') 'Internal error computing appropriate types for Henry''s', &
					'Law species.'
				goto 667
			endif
		endif
	enddo
	! Check tracer injection start/stop times.
	if (ztime .ne. 0) then
		do i = 1, nzones
			if (index(zonegrid(i, ztime), '>') .eq. 0) then
				if (zonegrid(i, ztime) .eq. '*') goto 267
				read(zonegrid(i, ztime), *, err=260, end=260) r
				if (r .lt. 0) then
					goto 265
				elseif (r .gt. tims) then
					goto 266
				endif
			else
				e = zonegrid(i, ztime)(1:index(zonegrid(i, ztime), '>') - 1)
				f = zonegrid(i, ztime)(index(zonegrid(i, ztime), '>') + 1:len_trim(zonegrid(i, ztime)))
				r = -1
				u = -1
				if (len_trim(e) .ne. 0) then
					read(e, *, end=260, err=260) r
					if (r .lt. 0) then
						goto 265
					elseif (r .gt. tims) then
						goto 266
					endif
				endif
				if (len_trim(f) .ne. 0) then
					read(f, *, end=260, err=260) u
					if (u .lt. 0) then
						goto 265
					elseif (u .gt. tims) then
						goto 266
					endif
				endif
				if ((r .ne. -1) .and. (u .ne. -1) .and. (abs(u) .lt. abs(r))) then
					write(ierr, '(a)') 'Error:  Stop time in zone "'//trim(zones(i))// &
						'" must occur after the start time.'
					goto 666
				endif
			endif
267			continue
		enddo
		goto 264
260		write(ierr, '(a)') 'Error:  Invalid start/stop time in zone "'//trim(zones(i))//'".'
		goto 666
265		write(ierr, '(a)') 'Error:  Start time in zone "'//trim(zones(i))//'" may not be negative.'
		goto 666
266		write(ierr, '(a)') 'Error:  Start/stop time in zone "'//trim(zones(i))//'" may not occur after the '// &
			'end of the simulation.'
		goto 666
264		continue
	endif

	! If we are in debug mode, print a short report to iout of what has been enabled.
	if (debug) then
		write(iout, *)
		write(iout, '(a)') '----- TRXN parameter report -----'
		write(iout, *)
		write(iout, '(a)') 'Components:'
		write(iout, '(a)') '(State) (Component Name)    (Master Species)    (Guess)   (Weight)'
		do i = 1, nspecies
			write(iout, '(a2, a6, a20, a20, es8.2, a2)', advance='no') states(i), '      ', &
				species(i), masters(i), guesses(i), '  '
				do j = 1, ncdenspecs
					if (species(i) .eq. cdenspecs(j)) then
						write(iout, '(f7.2)', advance='no') molweights(j)
						goto 329
					endif
				enddo
				write(iout, '(a)', advance='no') '      *'
329			write(iout, *)
		enddo
		write(iout, *)
		write(iout, '(a)') 'Water types:'
		if (nwt .eq. 0) write(iout, '(a)', advance='no') '(none)'
		do i = 1, nwt
			write(iout, '(a)', advance='no') trim(wtnames(i))//' = '
			flag = .false.
			do j = 1, nwtspecies
				if (wtgrid(i, j) .ne. 0) then
					if (flag) write(iout, '(a)', advance='no') ' + '
					write(iout, '(es10.2, a)', advance='no') wtgrid(i, j), ' m '//trim(wtspecies(j))
					flag = .true.
				endif
			enddo
		enddo
		write(iout, *)
		write(iout, *)
		write(iout, '(a)') 'Rock types:'
		if (nrt .eq. 0) write(iout, '(a)', advance='no') '(none)'
		do i = 1, nrt
			write(iout, '(a)', advance='no') trim(rtnames(i))//' = '
			flag = .false.
			do j = 1, nrtspecies
				if (rtgrid(i, j) .ne. 0) then
					if (flag) write(iout, '(a)', advance='no') ' + '
					write(iout, '(es10.2, a)', advance='no') rtgrid(i, j), ' m '//trim(rtspecies(j))
					flag = .true.
				endif
			enddo
		enddo
		write(iout, *)
		write(iout, *)
		write(iout, '(a)') 'Gas types:'
		if (nvt .eq. 0) write(iout, '(a)', advance='no') '(none)'
		do i = 1, nvt
			write(iout, '(a)', advance='no') trim(vtnames(i))//' = '
			flag = .false.
			do j = 1, nvtspecies
				if (vtgrid(i, j) .ne. 0) then
					if (flag) write(iout, '(a)', advance='no') ' + '
					write(iout, '(es10.2, a)', advance='no') vtgrid(i, j), ' m '//trim(vtspecies(j))
					flag = .true.
				endif
			enddo
		enddo
		write(iout, *)
		write(iout, *)
		write(iout, '(a)') 'Aqueous complexes:'
		if (ncomplexes .eq. 0) write(iout, '(a)') '(none)'
		do i = 1, ncomplexes
			write(iout, '(a, f6.2, a)', advance='no') trim(complexes(i))//' = ', complexstoich(i, 1), &
				' '//trim(complexcontents(i, 1))
			do j = 2, SPEC_MAX
				if (complexcontents(i, j) .eq. '*') exit
				write(iout, '(a, f6.2, a)', advance='no') ' + ', complexstoich(i, j), &
					' '//trim(complexcontents(i, j))
			enddo
			if ((flookup) .and. (equconstants(i) .ne. 0)) then
				write(iout, '(a, es10.2, a)') '    (constant K = ', equconstants(i), ')'
			elseif (flookup) then
				write(iout, '(a)') '    (5-parameter fit K)'
			elseif (cequi) then
				write(iout, '(a)') '    (piecewise K)'
			else
				write(iout, '(a, es10.2, a)') '    (constant K = ', equconstants(i), ')'
			endif
		enddo
		write(iout, *)
		write(iout, '(a)') 'Component groups:'
		if (ngroups .eq. 0) write(iout, '(a)') '(none)'
		do i = 1, ngroups
			write(iout, '(i3, a)', advance='no') i, ':  '
			k = 0
			do j = 1, nspecies
				if (states(j) .eq. 'a') then
					k = k + 1
					if (group(i, k) .ne. 0) write(iout, '(a)', advance='no') trim(species(j))//' '
				endif
			enddo
			write(iout, *)
		enddo
		write(iout, *)
		write(iout, '(a)') 'Reactions:'
		if (nrxns .eq. 0) write(iout, '(a)') '(none)'
		do i = 1, nrxns
			write(iout, '(f6.2, a)', advance='no') stoichiometries(i, 1, 1), ' '//trim(reactants(i, 1))
			do j = 2, SPEC_MAX
				if (reactants(i, j) .eq. '*') exit
				write(iout, '(a, f6.2, a)', advance='no') ' + ', stoichiometries(i, 1, j), &
					' '//trim(reactants(i, j))
			enddo
			write(iout, '(a)', advance='no') ' = '
			write(iout, '(f6.2, a)', advance='no') stoichiometries(i, 2, 1), ' '//trim(products(i, 1))
			do j = 2, SPEC_MAX
				if (products(i, j) .eq. '*') exit
				write(iout, '(a, f6.2, a)', advance='no') ' + ', stoichiometries(i, 2, j), &
					' '//trim(products(i, j))
			enddo
			write(iout, *)
		enddo
		write(iout, *)
		write(iout, '(100a20)') 'Zones:              ', (zonespecs(i), i = 1, nzonespecs)
		do i = 1, nzones + 1
			write(iout, '(100a20)') zones(i), (zonegrid(i, j), j = 1, nzonespecs)
		enddo
		write(iout, *)
		write(iout, '(a)') '----- --------------------- -----'
		write(iout, *)
	endif

	! Map all internal variables to FEHM's variables.
	! Set values for trac.
	if (debug) write(iptty, '(a)') 'Applying data for trac...'
	ctol = 1.e-06
      	iret = 0
	if (nspeci .ne. 0) then
	        ntpp = n7 / nspeci
	else
		ntpp = 0
	endif
	allocate(spnam(nspeci))
	do nsp = 1, nspeci
		npt(nsp) = (nsp - 1) * ntpp
		qcout(nsp) = 0.0
		qcin(nsp) = 0.0
		qcrxn(nsp) = 0.0
	enddo
	! Groups 1-3 have already been set in header.  Groups 4-7 are set here.
	if (tpor_flag) then
		if (ztpor .eq. 0) then
			write(ierr, '(a)') 'Internal error setting tracer porosity.'
			goto 667
		endif
		ps_trac = 0.32
		do i = 1, n0
			if (zonegrid(zoneresolv(izonef(i)), ztpor) .eq. '*') then
				ps_trac(i) = ps(i)
			else
				read(zonegrid(zoneresolv(izonef(i)), ztpor), *, err=341) ps_trac(i)
				if (ps_trac(i) .le. 0) then
					write(ierr, '(a)') 'Error:  Tracer porosity must be positive.'
					goto 666
				endif
			endif
		enddo
	endif
	goto 342
341	write(ierr, '(a)') 'Error:  Tracer porosity is not a valid number.'
	goto 666
342	if ((dispmode .ne. 0) .and. (dispmode .ne. 1)) then
		write(ierr, '(a)') 'Internal error determining dispersivity mode.'
		goto 667
	endif
	ldsp = dispmode
	! Set groups 8 and 9
	dispsame = 1
	if (ndisp + 1 .ne. numd) then
		write(ierr, '(a)') 'Internal error determining number of dispersivity models.'
		goto 667
	endif
	if (nsorp + 1 .ne. numsorp) then
		write(ierr, '(a)') 'Internal error determining number of sorption models.'
		goto 667
	endif
	if ((ldiffm .lt. 0) .or. (ldiffm .gt. 3) .or. (vdiffm .lt. 0) .or. (vdiffm .gt. 3)) then
		write(ierr, '(a)') 'Internal error setting diffusion models.'
		goto 666
	endif
	do i = 1, numd
		mflagl(1, i) = ldiffm
		mflagv(1, i) = vdiffm
		sehdiff(i) = ldiff
		sehdiffv(i) = vdiff
		diffmfl = ldiff
		diffmfv = vdiff
		if (ldsp .eq. 0) then
			tclx(1, i) = lx(i)
			tcly(1, i) = ly(i)
			tclz(1, i) = lz(i)
			tcvx(1, i) = vx(i)
			tcvy(1, i) = vy(i)
			tcvz(1, i) = vz(i)
		else
			tclx(1, i) = ll(i)
			tcly(1, i) = lt(i)
			tclz(1, i) = 0
			tcvx(1, i) = vl(i)
			tcvy(1, i) = vt(i)
			tcvz(1, i) = 0
		endif
	enddo
	! Set group 10
	do i = 1, n0
		a = 0
		if (zonegrid(zoneresolv(izonef(i)), zdisp) .eq. '*') then
			a = ndisp + 1
			goto 268
		endif
		do j = 1, ndisp
			if (dispnames(j) .eq. zonegrid(zoneresolv(izonef(i)), zdisp)) then
				a = j
				goto 268
			endif
		enddo
268		if (a .eq. 0) then
			write(ierr, '(a)') 'Internal error setting dispersivity models.'
			goto 667
		endif
		itrcdsp(i) = a
	enddo
	! Set group 11
	do i = 1, nspecies
		spnam(i) = species(i)
	enddo
	! Set groups 12, 13 (mflag)
!	do i = 1, nspecies
!		do j = 1, numsorp
!			mflagl(i, j) = 0
!			mflagv(i, j) = 0
!		enddo
!	enddo
	if (molein) then
		open(20, file='moles_in.tmp', err=331)
		goto 336
331		write(ierr, '(a)') 'Error:  Could not open temporary mole input storage file.  '// &
			'Ensure that write access is available in the current directory.'
		goto 666
	endif
336	c = 1
	do i = 1, nspecies
		conc_read(i) = .true.
		! Set group 11 values
		if (states(i) .eq. 'h') then
			hvliquid = 1
			hvvapor = 1
			icns(i) = 2
		elseif (states(i) .eq. 'i') then
			hvliquid = 1
			hvvapor = 1
			icns(i) = -2
		elseif (states(i) .eq. 'g') then
			hvvapor = 1
			icns(i) = -1
		elseif (states(i) .eq. 's') then
			icns(i) = 0
		elseif (states(i) .eq. 'a') then
			hvliquid = 1
			icns(i) = 1
		else
			write(ierr, '(a)') 'Internal error setting state data.'
			goto 667
		endif
		! Set group 12 values (sorption)
		if ((states(i) .eq. 's') .or. (states(i) .eq. 'j')) then
			do j = 1, numsorp
				iadsfl(i, j) = 0
		                a1adfl(i, j) = 0
		                a2adfl(i, j) = 0
		                betadfl(i, j) = 1
				iadsfv(i, j) = 0
		                a1adfv(i, j) = 0
		                a2adfv(i, j) = 0
		                betadfv(i, j) = 1
			enddo
		elseif ((states(i) .eq. 'a') .or. (states(i) .eq. 'g') .or. (states(i) .eq. 'h') &
			.or. (states(i) .eq. 'i')) then
			k = 0 ! Set k to be the index of the current species in anspecies.
			do a = 1, anspecies
				if ((species(i) .eq. sorpspecs(a)) .or. (sorpspecs(a) .eq. '*')) then
					k = a
					exit
				endif
			enddo
			if (k .ne. 0) then
				do j = 1, nsorp
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
				do j = 1, nsorp
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
			write(ierr, '(a)') 'Internal error setting state data.'
			goto 667
		endif
		iadsfl(i, nsorp + 1) = 0 ! Default model
		a1adfl(i, nsorp + 1) = 0
                a2adfl(i, nsorp + 1) = 0
                betadfl(i, nsorp + 1) = 1
		iadsfv(i, nsorp + 1) = 0
                a1adfv(i, nsorp + 1) = 0
                a2adfv(i, nsorp + 1) = 0
                betadfv(i, nsorp + 1) = 1
		! Group 13
		do j = 1, n0
			e = zonegrid(zoneresolv(izonef(j)), zsorp)
			if (e .eq. '*') then
				itrc((i-1)*n0 + j) = nsorp + 1
			else
				a = 0
				do k = 1, nsorp + 1
					if (sorpnames(k) .eq. e) then
						a = k
						exit
					endif
				enddo
				if (a .eq. 0) then
					write(ierr, '(a)') 'Internal error setting adsorption data.'
					goto 667
				endif
				itrc((i-1)*n0 + j) = a
			endif
		enddo
		! Set group 14 values if we're Henry's Law
		if ((states(i) .eq. 'h') .or. (states(i) .eq. 'i') .or. (states(i) .eq. 'j')) then
			d = 0
			do o = 1, hnhtspecies
				if (hhtspecies(o) .eq. species(i)) then
					d = o
					exit
				endif
			enddo
			if (d .eq. 0) then
				write(ierr, '(a)') 'Internal error processing Henry''s Law models.'
				goto 667
			endif
			henry_model(i) = hmodels(d)
			a = 0
			do o = 1, hnhtspecies
				if (hhtspecies(o) .eq. species(i)) then
					a = o
					exit
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a)') 'Internal error setting Henry''s Law parameters.'
				goto 667
			endif
			if (henry_model(i) .eq. 1) then
				a_henry(c) = ah(a)
				dh_henry(c) = dhh(a)
				hawwa(c, 1) = 0
				hawwa(c, 2) = 0
				hawwa(c, 3) = 0
				hawwa(c, 4) = 0
				hawwa(c, 5) = 0
			elseif (henry_model(i) .eq. 2) then
				hawwa(c, 1) = ah1(a)
				hawwa(c, 2) = ah2(a)
				hawwa(c, 3) = ah3(a)
				hawwa(c, 4) = ah4(a)
				hawwa(c, 5) = ah5(a)
			elseif (henry_model(i) .eq. 3) then
				a_henry(c) = ah(a)
				dh_henry(c) = hh(a)
				hawwa(c, 1) = 0
				hawwa(c, 2) = 0
				hawwa(c, 3) = 0
				hawwa(c, 4) = 0
				hawwa(c, 5) = 0
			else
				write(ierr, '(a)') 'Internal error setting Henry''s Law parameters.'
				goto 667
			endif
			c = c + 1
		endif
		! Set group 15, 16 values
		do k = (i - 1) * n0 + 1, i * n0
			an(k) = an0
		enddo
		if (molein) then ! Mole input option
			write(20, *) i
			if ((zonefile .ne. '*') .and. (i .eq. 1)) then
				write(20, *) 'file'
				write(20, *) zonefile
			endif
			write(20, *) '1 0 0 0'
			o = 0
			do k = 1, nmspecs
				if (mspecs(k) .eq. species(i)) then
					o = k
					exit
				endif
			enddo
			if (o .ne. 0) then
				do j = 1, nmzones
					write(20, *) -1 * mzones(j), 0, 0, molegrid(j, o)
				enddo
			endif
			write(20, *)
		endif
		do j = 1, n0
			k = (i - 1) * n0 + j
			inflag = .true.
			if (zonegrid(zoneresolv(izonef(j)), zboun) .eq. '*') then
				cnsk(k) = 0
			endif
			if ((states(i) .eq. 'a') .or. (states(i) .eq. 'h')) then
				if ((zonegrid(zoneresolv(izonef(j)), zinit) .eq. '*') .or. (molein)) then
					inflag = .false.
					an((i - 1) * n0 + j) = an0
				endif
				if (inflag) then
					a = 0
					do o = 1, nwtspecies
						if (wtspecies(o) .eq. species(i)) then
							a = o
							exit
						endif
					enddo
					d = 0
					do o = 1, nwt
						if (wtnames(o) .eq. zonegrid(zoneresolv(izonef(j)), zinit)) then
							d = o
							exit
						endif
					enddo
					if (a .eq. 0) then
						an((i - 1) * n0 + j) = 0
					elseif (d .eq. 0) then
						write(ierr, '(a)') 'Internal error finding a match for inflow '// &
							'species in water.'
						goto 666
					else
						an((i - 1) * n0 + j) = wtgrid(d, a)
					endif
				endif
				if (zonegrid(zoneresolv(izonef(j)), zboun) .eq. '*') then
					cnsk(k) = 0
				else
					a = 0
					do o = 1, nwtspecies
						if (wtspecies(o) .eq. species(i)) then
							a = o
							exit
						endif
					enddo
					d = 0
					e = zonegrid(zoneresolv(izonef(j)), zboun)//'.'
					do o = 1, nwt
						if (index(e, wtnames(o)//'.') .ne. 0) then
							d = o
							exit
						endif
					enddo
					if (a .eq. 0) then
						cnsk(k) = 0
					elseif (d .eq. 0) then
						write(ierr, '(a)') 'Internal error setting inflow concentration data.'
						goto 667
					else
						cnsk(k) = wtgrid(d, a)
					endif
				endif
				e = zonegrid(zoneresolv(izonef(j)), ztime)
				if (e .eq. '*') then
					if (zonegrid(zoneresolv(izonef(j)), zboun) .eq. '*') then
						t1sk(k) = 0
						t2sk(k) = 0
					else
						t1sk(k) = 0
						t2sk(k) = tims
					endif
				elseif (e .eq. '0') then
					t1sk(k) = 0
					t2sk(k) = 0
				elseif (index(e, '>') .eq. 0) then
					read(e, *, err=250) t1sk(k)
					t2sk(k) = t1sk(k) + 1
				else
					f = e(1:index(e, '>') - 1)
					g = e(index(e, '>') + 1:len(e))
					if (len_trim(f) .eq. 0) then
						t1sk(k) = 0
					else
						read(f, *, err=250, end=250) t1sk(k)
					endif
					if (len_trim(g) .eq. 0) then
						t2sk(k) = tims
					else
						read(g, *, err=250, end=250) t2sk(k)
					endif
				endif
				goto 251
250					write(ierr, '(a)') 'Error setting aqueous tracer injection start/stop times.'
				goto 666
251				continue
			elseif (states(i) .eq. 's') then
				if (zonegrid(zoneresolv(izonef(j)), zrock) .eq. '*') then
					inflag = .false.
					an((i - 1) * n0 + j) = an0
				endif
				inflag = .true.
				if ((zonegrid(zoneresolv(izonef(j)), zrock) .eq. '*') .or. (molein)) then
					inflag = .false.
					an((i - 1) * n0 + j) = an0
				endif
				if (inflag) then
					a = 0
					do o = 1, nrtspecies
						if (rtspecies(o) .eq. species(i)) then
							a = o
							exit
						endif
					enddo
					d = 0
					do o = 1, nrt
						if (rtnames(o) .eq. zonegrid(zoneresolv(izonef(j)), zrock)) then
							d = o
							exit
						endif
					enddo
					if (a .eq. 0) then
						an((i - 1) * n0 + j) = 0
					elseif (d .eq. 0) then
						write(ierr, '(a)') 'Internal error finding a match for species in rock.'
						goto 666
					else
						an((i - 1) * n0 + j) = rtgrid(d, a)
					endif
				endif
!				stop
				cnsk(k) = 0	! Solids can't flow!
				t1sk(k) = 0
				t2sk(k) = 0
			elseif ((states(i) .eq. 'g') .or. (states(i) .eq. 'i')) then
				if ((zonegrid(zoneresolv(izonef(j)), zvap) .eq. '*') .or. (molein)) then
					inflag = .false.
					an((i - 1) * n0 + j) = an0
				endif
				if (inflag) then
					a = 0
					do o = 1, nvt
						if (vtnames(o) .eq. zonegrid(zoneresolv(izonef(j)), zvap)) then
							a = o
							exit
						endif
					enddo
					d = 0
					do o = 1, nvtspecies
						if (vtspecies(o) .eq. species(i)) then
							d = o
							exit
						endif
					enddo
					if (d .eq. 0) then
						an((i - 1) * n0 + j) = 0
					elseif (a .eq. 0) then
						write(ierr, '(a)') 'Internal error setting gas data.'
						goto 667
					else
						an((i - 1) * n0 + j) = vtgrid(a, d)
					endif
				endif
				if (zonegrid(zoneresolv(izonef(j)), zboun) .eq. '*') then
					cnsk(k) = 0
				else
					a = 0
					do o = 1, nvtspecies
						if (vtspecies(o) .eq. species(i)) then
							a = o
							exit
						endif
					enddo
					d = 0
					e = zonegrid(zoneresolv(izonef(j)), zboun)//'.'
					do o = 1, nvt
						if (index(e, vtnames(o)//'.') .ne. 0) then
							d = o
							exit
						endif
					enddo
					if (a .eq. 0) then
						cnsk(k) = 0
					elseif (d .eq. 0) then
						write(ierr, '(a)') 'Internal error setting boundary gas concentration data.'
						goto 667
					else
						cnsk(k) = vtgrid(d, a)
					endif
				endif
				e = zonegrid(zoneresolv(izonef(j)), ztime)
				if (e .eq. '*') then
					if (zonegrid(zoneresolv(izonef(j)), zboun) .eq. '*') then
						t1sk(k) = 0
						t2sk(k) = 0
					else
						t1sk(k) = 0
						t2sk(k) = tims
					endif
				elseif (e .eq. '0') then
					t1sk(k) = 0
					t2sk(k) = 0
				elseif (index(e, '>') .eq. 0) then
					read(e, *, err=254) t1sk(k)
					t2sk(k) = t1sk(k) + 1
				else
					f = e(1:index(e, '>') - 1)
					g = e(index(e, '>') + 1:len(e))
					if (f .eq. '') then
						t1sk(k) = 0
					else
						read(f, *, err=254, end=254) t1sk(k)
					endif
					if (g .eq. '') then
						t2sk(k) = tims
					else
						read(g, *, err=254, end=254) t2sk(k)
					endif
				endif
				goto 255
254					write(ierr, '(a)') 'Error setting gaseous tracer injection start/stop times.'
				goto 666
255					continue
			endif
			! Set solute accumulation and constant concentration options
			zfaccum = .false.
			zfconst = .false.
			e = zonegrid(zoneresolv(izonef(j)), zopt)
			if (e .ne. '*') then
				e = trim(e)//'.'
				do
					do
						if (index(e, '.') .eq. 1) then
							e = e(2:len_trim(e))
						else
							exit
						endif
					enddo
					if (index(e, '.') .eq. 0) exit
					f = e(1:index(e, '.') - 1)
					e = e(index(e, '.') + 1:len_trim(e))
					if (f(1:5) .eq. 'accum') then
						if (zfconst) then
							write(ierr, '(a)') 'Error:  Solute accumulation and '// &
								'constant concentration at inflow nodes cannot '// &
								'both be enabled.'
							goto 666
						endif
						zfaccum = .true.
						pcnsk(k) = 1
					elseif (f(1:5) .eq. 'const') then
						if (zfaccum) then
							write(ierr, '(a)') 'Error:  Solute accumulation and '// &
								'constant concentration at inflow nodes cannot '// &
								'both be enabled.'
							goto 666
						endif
						zfconst = .true.
						cnsk(k) = -1 * cnsk(k)
						pcnsk(k) = -1
					else
						write(ierr, '(a)') 'Error:  Zone option "'//trim(f)//'" not recognized.'
						goto 666
					endif
				enddo
			endif
			anlo(j) = an(j)
		enddo
	enddo
	if (molein) close(20)
	! Set values for cden.
	if (cden) then
		allocate(mw_speci(ncpnt))
		mw_speci = 0
		do i = 1, ncdenspecs
			do j = 1, nspecies
				if ((cdenspecs(i) .eq. masters(j)) .or. (cdenspecs(i) .eq. species(j)) .and. &
					(cdenspecs(i) .ne. '*')) then
					if ((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
						write(ierr, '(a)') 'Error:  Master species "'//trim(cdenspecs(i))// &
							'" in cden is not aqueous or aqueous Henry''s Law.'
						goto 666
					endif
					do k = 1, ncpnt
!						if (trim(cdenspecs(i)) .eq. trim(cpntnam(k))) then
						if (trim(cdenspecs(i)) .eq. trim(species(k))) then

							mw_speci(k) = molweights(i)
							goto 273
						endif
					enddo
					write(ierr, '(a)') 'Error:  Master species "'//trim(cdenspecs(i))//'" in cden not found in comp.'
					write(ierr,*) ncpnt,(cpntnam(k),k=1,ncpnt)
					goto 666
				endif
			enddo
	273		continue
		enddo
	endif
	! Call userc
	if (fuserc) then
		if (debug) write(iptty, '(a)') 'Calling user transport subroutine in '//trim(usercfname)//'.'
		open(20, status='scratch', err=337)
		write(20, '(a4)') 'file'
		write(20, '(a100)') usercfname
		rewind 20
		inpttmp = inpt
		inpt = 20
		call userc(0, i, r, u)
		inpt = inpttmp
		close(20)
		if (debug) write(iptty, '(a)') 'User transport subroutine finished.'
		goto 338
337		write(ierr, '(a)') 'Error:  Could not open temporary file for userc.  '// &
			'Ensure that write access is enabled in the current directory.'
		goto 666
	endif
338	continue
	! Printing, guesses, names...
	ncpnt = 0
	nimm = 0
	nvap = 0
	cpntprt = 0
	immprt = 0
	vapprt = 0
	ncpntprt = 0
	nimmprt = 0
	nvapprt = 0
	ncplxprt = 0
	cplxprt = 0
! nspecies is number of components
	do i = 1, nspecies
		if ((states(i) .eq. 'a') .or. (states(i) .eq. 'h') .or. (states(i) .eq. 'i')) then
			ncpnt = ncpnt + 1
			cpntnam(ncpnt) = species(i)
            if (cpntnam(ncpnt) .eq. cden_spnam) ispcden = i
			pcpnt(ncpnt) = i
			if (rxn_flag .eq. 1) then
				cpntgs(ncpnt) = guesses(i)
				if ((cpntnam(ncpnt) .eq. 'H') .and. fph) then
					ifxconc(ncpnt) = 2
				else
					ifxconc(ncpnt) = 0 ! Removed support for free-ion concentration using ifxconc = 1
				endif
			endif
			do j = 1, nprint
				if (printspecies(j) .eq. species(i)) then
					ncpntprt = ncpntprt + 1
					cpntprt(ncpntprt) = ncpnt
				endif
			enddo
			!idcpnt(ncpnt) = ncpnt
		elseif (states(i) .eq. 'g') then
			nvap = nvap + 1
			vapnam(nvap) = species(i)
			pvap(nvap) = i
			do j = 1, nprint
				if (printspecies(j) .eq. species(i)) then
					nvapprt = nvapprt + 1
					vapprt(nvapprt) = nvap
				endif
			enddo
			!idvap(nvap) = nimm
		elseif (states(i) .eq. 's') then
			nimm = nimm + 1
			immnam(nimm) = species(i)
			pimm(nimm) = i
			do j = 1, nprint
				if (printspecies(j) .eq. species(i)) then
					nimmprt = nimmprt + 1
					immprt(nimmprt) = nimm
				endif
			enddo
			!idimm(nimm) = nimm
		else
			write(ierr, '(a)') 'Internal error setting species names.'
			goto 667
		endif
241		continue
	enddo
	! End of trac value-setting.
	if (rxn_flag .eq. 0) goto 6000
	! Set values for rxn.
	if (debug) write(iptty, '(a)') 'Applying data for rxn...'
	rxnon = 1
	rxnnaqueous = 0
	rxnnsolid = 0
	rxnnvapor = 0
	allocate(rxnaqueous(SPEC_MAX * 2))
	allocate(rxnsolid(SPEC_MAX * 2))
	allocate(rxnvapor(SPEC_MAX))
	rxnaqueous = '*'
	rxnsolid = '*'
	rxnvapor = '*'
	do i = 1, nspecies
		if ((states(i) .eq. 'a') .or. (states(i) .eq. 'h') .or. (states(i) .eq. 'i')) then
			rxnnaqueous = rxnnaqueous + 1
			rxnaqueous(rxnnaqueous) = species(i)
		elseif (states(i) .eq. 'g') then
			rxnnvapor = rxnnvapor + 1
			rxnvapor(rxnnvapor) = species(i)
		elseif (states(i) .eq. 's') then
			rxnnsolid = rxnnsolid + 1
			rxnsolid(rxnnsolid) = species(i)
		endif
	enddo
	! Printing of complexes
	do i = 1, ncomplexes
		do j = 1, nprint
			if (complexes(i) .eq. printspecies(j)) then
				ncplxprt = ncplxprt + 1
				cplxprt(ncplxprt + 100) = i + 100
			endif
		enddo
	enddo
	! Group 1 and 2 values have already been set in trxninit.
	! Groups 3 - 6 were set in trac, as they apply to all species.
	! Groups 7 and 8 were taken care of in header.
	! Set group 9, 10, and 11 values
	temp_model = ' '
	do i = 1, ncomplexes
		cplxnam(i + 100) = complexes(i)
		if (flookup) then ! PHREEQC-style multi-parameter fit
			if (equconstants(i) .ne. 0) then
				ckeq(i + 100) = equconstants(i)
				heq(i + 100, 1) = enthalpies(i)
			else
				temp_model(i + 100) = 't'
				ckeq(i + 100) = 0
				heq(i + 100, 1:5) = heqfit(i, 1:5)
			endif
		elseif (cequi) then ! If we used equi to set ckeq as a function of temperature
			temp_model(i + 100) = 'l'
			do j = 1, necomplexes
				if (complexes(i) .eq. ecomplexes(j)) then
					a = 0
					do k = 1, neqtemps
						if (eqtemps(j, k) .eq. -1) then
							a = k - 1
							exit
						endif
					enddo
					if (a .eq. 0) a = neqtemps
					rarray1 = 0
					rarray2 = 0
					do k = 1, a
						rarray1(k) = eqtemps(j, k)
						rarray2(k) = eqconsts(j, k)
					enddo
					call lstsq(a, rarray2, rarray1)
					do k = 1, 3
						heq(i + 100, k) = rarray2(k)
					enddo
					ckeq(i + 100) = heq(i + 100, 1) + heq(i + 100, 2) * 25 + heq(i + 100, 3) * 625
					ckeq(i + 100) = 10 ** ckeq(i + 100)
					goto 278
				endif
			enddo
			write(ierr, '(a)') 'Internal error setting temperature-variant equilibrium constants.'
			goto 667
278			continue
		else
			ckeq(i + 100) = equconstants(i)
			heq(i + 100, 1) = enthalpies(i)
		endif
		do j = 1, SPEC_MAX
			if (complexcontents(i, j) .eq. '*') then
				exit
			endif
			do k = 1, nspecies
				if ((masters(k) .eq. '*') .or. (masters(k) .eq. '')) then
					continue
				endif
				if (complexcontents(i, j) .eq. masters(k)) then
					goto 240
				endif
			enddo
			write(ierr, '(a)') 'Error:  Master species "'//trim(complexcontents(i, j))//'" not found in comp.'
			goto 666
240			continue
		enddo
	enddo
	spstoic = 0
	neg_conc_possible = 0
	do i = 1, ncomplexes
		o = 0
		do j = 1, nspecies
			if ((states(j) .eq. 'a') .or. (states(j) .eq. 'h')) then
				o = o + 1
			endif
			do k = 1, SPEC_MAX
				if (complexcontents(i, k) .eq. '*') then
					exit
				endif
				if (masters(j) .eq. complexcontents(i, k)) then
					spstoic(i + 100, o) = complexstoich(i, k)
					if (complexstoich(i, k) .lt. 0) neg_conc_possible(o) = 1
					exit
				endif
			enddo
		enddo
	enddo
	! Set group 12 for each reaction
	do i = 1, numrxn
		idrxn(i) = rxntypes(i)
	enddo
	temp_model_kin = ' '
	! Set reaction-specific groups
	do i = 1, nrxns
		if (rxntypes(i) .eq. 1) then
			naqsp(i) = 1
			nimsp(i) = 1
			nivsp(i) = 0
			do j = 1, nspecies
				if ((reactants(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
					do o = 1, nimm
						if (immnam(o) .eq. species(j)) irxnim(i, 1) = o
					enddo
					do k = 1, nspecies
						if (products(i, 1) .eq. masters(k)) then
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, 1) = o
							enddo
							goto 243
						endif
					enddo
					do k = 1, ncomplexes
						if (products(i, 1) .eq. complexes(k)) then
							irxnic(i, 1) = k + 100
							goto 243
						endif
					enddo
					goto 242
				else if ((products(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
					do o = 1, nimm
						if (immnam(o) .eq. species(j)) irxnim(i, 1) = o
					enddo
					do k = 1, nspecies
						if (reactants(i, 1) .eq. masters(k)) then
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, 1) = o
							enddo
							goto 243
						endif
					enddo
					do k = 1, ncomplexes
						if (reactants(i, 1) .eq. complexes(k)) then
							irxnic(i, 1) = k + 100
							goto 243
						endif
					enddo
					goto 242
				endif
			enddo
242			write(ierr, '(a, i3, a)') 'Error:  Type 1 reaction ', i, ' requires exactly one aqueous '// &
				'complex or master species and one solid component.'
			goto 666
243			read(dscoefs(i), *, err=206) r
			ckeqlb(i) = r
			goto 207
206			if (dscoefs(i) .eq. '*') then
				write(ierr, '(a)') 'Internal error setting distribution coefficients.'
				goto 667
			endif
			a = 0
			do j = 1, ndistmodels
				if (dscoefs(i) .eq. distmodelnames(j)) then
					temp_model_kin(i) = 'l'
					a = j
					exit
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a)') 'Internal error setting distribution coefficient information.'
				goto 667
			endif
			j = 0
			do
				if (distmodels(a, j + 1, 1) .eq. -1) exit
				j = j + 1
				rarray1(j) = distmodels(a, j, 1);
				rarray2(j) = distmodels(a, j, 2);
			enddo
			call lstsq(j, rarray2, rarray1)
			do j = 1, 3
				tcoeff(i, j) = rarray2(i)
			enddo
			ckeqlb(i) = tcoeff(i,1) + tcoeff(i, 2) * 25 + tcoeff(i, 3) * 625
207			ckmtrn(i) = rates(i)
		elseif (rxntypes(i) .eq. 2) then
			naqsp(i) = 1
			nimsp(i) = 1
			nivsp(i) = 0
			do j = 1, nspecies
				if ((reactants(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
					do o = 1, nimm
						if (immnam(o) .eq. species(j)) irxnim(i, 1) = o
					enddo
					do k = 1, nspecies
						if (products(i, 1) .eq. masters(k)) then
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, 1) = o
							enddo
							goto 245
						endif
					enddo
					do k = 1, ncomplexes
						if (products(i, 1) .eq. complexes(k)) then
							irxnic(i, 1) = k + 100
							goto 245
						endif
					enddo
					goto 244
				else if ((products(i, 1) .eq. species(j)) .and. (states(j) .eq. 's')) then
					do o = 1, nimm
						if (immnam(o) .eq. species(j)) irxnim(i, 1) = o
					enddo
					do k = 1, nspecies
						if (reactants(i, 1) .eq. masters(k)) then
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, 1) = o
							enddo
							goto 245
						endif
					enddo
					do k = 1, ncomplexes
						if (reactants(i, 1) .eq. complexes(k)) then
							irxnic(i, 1) = k + 100
							goto 245
						endif
					enddo
					goto 244
				endif
			enddo
244			write(ierr, '(a, i3, a)') 'Error:  Type 2 reaction ', i, ' requires exactly one aqueous '// &
				'complex or master species and one solid component.'
			goto 666
245			read(dscoefs(i), *, err=208) r
			ckeqlb(i) = r
			goto 209
208			a = 0
			do j = 1, ndistmodels
				if (dscoefs(i) .eq. distmodelnames(j)) then
					temp_model_kin(i) = 'l'
					a = j
					exit
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a)') 'Internal error setting distribution coefficient information.'
				goto 667
			endif
			j = 0
			do
				if (distmodels(a, j + 1, 1) .eq. -1) exit
				j = j + 1
				rarray1(j) = distmodels(a, j, 1);
				rarray2(j) = distmodels(a, j, 2);
			enddo
			call lstsq(j, rarray2, rarray1)
			do j = 1, 3
				tcoeff(i, j) = rarray2(i)
			enddo
			ckeqlb(i) = tcoeff(i, 1) + tcoeff(i, 2) * 25 + tcoeff(i, 3) * 625
209			ckmtrn(i) = rates(i)
			simmmx(i) = xcoef(i)
		elseif (rxntypes(i) .eq. 3) then
			nimsp(i) = 0
			naqsp(i) = 0
			nivsp(i) = 0
			do j = 1, SPEC_MAX
				if (reactants(i, j) .eq. '*') then
					exit
				endif
				do k = 1, nspecies
					if ((species(k) .eq. reactants(i, j)) .or. (masters(k) .eq. reactants(i, j))) then
						if ((states(k) .eq. 'a') .or. (states(k) .eq. 'h') .or. &
							(states(k) .eq. 'i')) then
							naqsp(i) = naqsp(i) + 1
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
							enddo
							sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
						elseif (states(k) .eq. 's') then
							nimsp(i) = nimsp(i) + 1
							do o = 1, nimm
								if (immnam(o) .eq. species(k)) irxnim(i, nimsp(i)) = o
							enddo
							stimirrv(i, nimsp(i)) = stoichiometries(i, 1, j)

						elseif (states(k) .eq. 'g') then
							nivsp(i) = nivsp(i) + 1
							do o = 1, nvap
								if (vapnam(o) .eq. species(k)) irxniv(i, nivsp(i)) = o
							enddo
							stivirrv(i, nivsp(i)) = stoichiometries(i, 1, j)
						endif
						goto 210
					endif
				enddo
				do k = 1, ncomplexes
					if (complexes(k) .eq. reactants(i, j)) then
						naqsp(i) = naqsp(i) + 1
						irxnic(i, naqsp(i)) = k + 100
						sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
						goto 210
					endif
				enddo
				write(ierr, '(a)') 'Error setting reactant information.'
				goto 666
210				continue
			enddo
			do j = 1, SPEC_MAX
				if (products(i, j) .eq. '*') then
					exit
				endif
				do k = 1, nspecies
					if ((species(k) .eq. products(i, j)) .or. (masters(k) .eq. products(i, j))) then
						if ((states(k) .eq. 'a') .or. (states(k) .eq. 'h') .or. &
							(states(k) .eq. 'i')) then
							naqsp(i) = naqsp(i) + 1
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
							enddo
							sticirrv(i, naqsp(i)) = - stoichiometries(i, 2, j)
						elseif (states(k) .eq. 's') then
							nimsp(i) = nimsp(i) + 1
							do o = 1, nimm
								if (immnam(o) .eq. species(k)) irxnim(i, nimsp(i)) = o
							enddo
							stimirrv(i, nimsp(i)) = - stoichiometries(i, 2, j)
						elseif ((states(k) .eq. 'g')) then
							nivsp(i) = nivsp(i) + 1
							do o = 1, nvap
								if (vapnam(o) .eq. species(k)) irxniv(i, nivsp(i)) = o
							enddo
							stivirrv(i, nivsp(i)) = - stoichiometries(i, 2, j)
						endif
						goto 211
					endif
				enddo
				do k = 1, ncomplexes
					if (complexes(k) .eq. reactants(i, j)) then
						naqsp(i) = naqsp(i) + 1
						irxnic(i, naqsp(i)) = k + 100
						sticirrv(i, naqsp(i)) = stoichiometries(i, 1, j)
						goto 211
					endif
				enddo
				write(ierr, '(a)') 'Error setting product information.'
				goto 666
211				continue
			enddo
			kfor(i) = rates(i)
			krev(i) = rates2(i)
		elseif (rxntypes(i) .eq. 4) then
			naqsp(i) = 0
			nimsp(i) = 0
			nivsp(i) = 0
			read(bioparams(i, 1), *, end=218, err=230) f, e
			goto 219
218			e = f
			f = '1'
			read(f, *, err=602) o
219			do j = 1, nspecies
				if (e .eq. species(j)) then
					if ((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
						write(ierr, '(a, i3, a)') 'Error:  The substrate "'//trim(e)// &
							'" in type '// '4 reaction ', i, ' must be aqueous or '// &
							'aqueous Henry''s Law.'
						goto 666
					endif
					if (o .ne. 1) then
						write(ierr, '(a, i3, a)') 'Error:  The stoichiometry of the '// &
							'substrate "'//trim(e)//'" in type 4 reaction ', i, ' must be 1.'
						goto 666
					endif
					naqsp(i) = naqsp(i) + 1
					do k = 1, ncpnt
						if (cpntnam(k) .eq. species(j)) irxnic(i, naqsp(i)) = k
					enddo
				endif
			enddo
			if (naqsp(i) .eq. 0) then
				write(ierr, '(a, i3, a)') 'Error:  Substrate "'//e//'" in type 4 reaction ', i, &
					' not found in comp.'
				goto 666
			endif
			read(bioparams(i, 2), *, end=223, err=230) f, e
			goto 224
223			e = f
			f = '1'
			read(f, *, err=602) o
224			do j = 1, nspecies
				if (e .eq. species(j)) then
					if ((states(j) .ne. 'a') .and. (states(j) .ne. 'h')) then
						write(ierr, '(a, i3, a)') 'Error:  The electron acceptor "'//trim(e)// &
							'" in type 4 reaction ', i, ' must be aqueous or '// &
							'aqueous Henry''s Law.'
						goto 666
					endif
					naqsp(i) = naqsp(i) + 1
					do k = 1, ncpnt
						if (cpntnam(k) .eq. species(j)) irxnic(i, naqsp(i)) = k
					enddo
					biofac(j) = o
				endif
			enddo
			if (naqsp(i) .eq. 1) then
				write(ierr, '(a, i3, a)') 'Error:  Electron acceptor "'//trim(e)// &
					'" in type 4 reaction ', i, ' not found in comp.'
				goto 666
			endif
			read(bioparams(i, 3), *, end=225, err=230) f, e
			goto 226
225			e = f
			f = '1'
			read(f, *, err=602) o
226			do j = 1, nspecies
				if (e .eq. species(j)) then
					if (states(j) .ne. 's') then
						write(ierr, '(a, i3, a)') 'Error:  The biomass "'//trim(e)// &
							'" in type 4 reaction ', i, ' must be solid.'
						goto 666
					endif
					nimsp(i) = nimsp(i) + 1
					do k = 1, nimm
						if (immnam(k) .eq. species(j)) irxnim(i, nimsp(i)) = k
					enddo
				endif
			enddo
			if (nimsp(i) .eq. 0) then
				write(ierr, '(a, i3, a)') 'Error:  Biomass "'//trim(e)//'" in type 4 reaction ', i, &
					' not found in comp.'
				goto 666
			endif
			do j = 1, SPEC_MAX
				if (reactants(i, j) .eq. '*') then
					exit
				endif
				naqsp(i) = naqsp(i) + 1
				if (naqsp(i) .gt. 5) then
					write(ierr, '(a, i3, a)') 'Error:  Type 4 reaction ', i, &
						' can only have a total of three extra reactants and products.'
					goto 666
				endif
				do k = 1, nspecies
					if (reactants(i, j) .eq. species(k)) then
						if ((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
							write(ierr, '(a, i3, a)') 'Error:  Extra reactant "'// &
								trim(reactants(i, j))//'" in type 4 reaction ', i, &
								' must be aqueous or Henry''s Law.'
							goto 666
						endif
						do o = 1, ncpnt
							if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
						enddo
						if (naqsp(i) .eq. 3) then
							hfac(i) = stoichiometries(i, 1, j)
						elseif (naqsp(i) .eq. 4) then
							carbfac(i) = stoichiometries(i, 1, j)
						elseif (naqsp(i) .eq. 5) then
							ammfac(i) = stoichiometries(i, 1, j)
						else
							write(ierr, '(a)') 'Internal error setting biodegradation '// &
								'extra species data.'
							goto 667
						endif
					endif
				enddo
			enddo
			goto 231
230			write(ierr, '(a, i3, a)') 'Error:  Bad format for species in type 4 reaction ', i, '.'
			goto 666
231			do j = 1, SPEC_MAX
				if (products(i, j) .eq. '*') then
					exit
				endif
				naqsp(i) = naqsp(i) + 1
				if (naqsp(i) .gt. 5) then
					write(ierr, '(a, i3, a)') 'Error:  Type 4 reaction ', i, ' can only have a '// &
						'total of three extra reactants and products.'
					goto 666
				endif
				do k = 1, nspecies
					if (products(i, j) .eq. species(k)) then
						if ((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
							write(ierr, '(a, i3, a)') 'Error:  Extra product "'// &
								trim(reactants(i, j))//'" in type 4 reaction ', i, &
								' must be aqueous or aqueous Henry''s Law.'
							goto 666
						endif
						do o = 1, ncpnt
							if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
						enddo
						if (naqsp(i) .eq. 3) then
							hfac(i) = - stoichiometries(i, 1, j)
						elseif (naqsp(i) .eq. 4) then
							carbfac(i) = - stoichiometries(i, 1, j)
						elseif (naqsp(i) .eq. 5) then
							ammfac(i) = - stoichiometries(i, 1, j)
						else
							write(ierr, '(a)') 'Internal error setting biodegradation '// &
								'extra species data.'
							goto 667
						endif
					endif
				enddo
			enddo
			read(bioparams(i, 4), *, err=602) ckc(i)
			read(bioparams(i, 5), *, err=602) cka(i)
			decay(i) = rates2(i)
			read(bioparams(i, 6), *, err=602) phthresh(i)
			qm(i) = rates(i)
			read(bioparams(i, 7), *, err=602) yield(i)
			read(bioparams(i, 8), *, err=602) xminit(i)
			nbiofrm(i) = 0
			do j = 1, SPEC_MAX
				if (icbioholder(i, j) .eq. '*') then
					exit
				endif
				do k = 1, nspecies
					if (icbioholder(i, j) .eq. species(k)) then
						if ((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
							write(ierr, '(a, i3, a)') 'Error:  Biodegradable form "'// &
								trim(icbioholder(i, j))//'" in type 4 reaction ', i, &
								' must be aqueous or aqueous Henry''s Law.'
							goto 666
						endif
						nbiofrm(i) = nbiofrm(i) + 1
						icbio(i, nbiofrm(i)) = k
					endif
				enddo
				if (nbiofrm(i) .ne. j) then
					write(ierr, '(a, i3, a)') 'Error:  Biodegradable form "'//trim(icbioholder(i, j))// &
						'" in type 4 reaction ', i, ' not found in comp.'
					goto 666
				endif
			enddo
		elseif (rxntypes(i) .eq. 5) then
			ckmtrn(i) = rates(i)
			ckmtrn(i) = log(2.0)/(ckmtrn(i) * 365.25 * 24)
			a = 0
			do j = 1, nspecies
				if (reactants(i, 1) .eq. species(j)) then
					a = j
					exit	
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a, i3, a)') 'Error:  Component "'//trim(reactants(i, 1))//'" in reaction ', &
					i, ' not found in comp.'
				goto 666
			endif
			c = 0
			do j = 1, nspecies
				if (products(i, 1) .eq. '*') then
					exit
				elseif (products(i, 1) .eq. species(j)) then
					c = j
					exit	
				endif
			enddo
			if (states(a) .ne. states(c)) then
				write(ierr, '(a, i3, a)') 'The parent and daughter species in type 5 reaction ', i, &
					' must be in the same state.'
				goto 666
			endif
			if ((states(a) .eq. 'a') .or. (states(a) .eq. 'h') .or. (states(k) .eq. 'i')) then
				naqsp(i) = 2
				do j = 1, ncpnt
					if (cpntnam(j) .eq. species(a)) irxnic(i, 1) = j
				enddo
				if (c .eq. 0) then
					naqsp(i) = 1
				else
					do j = 1, ncpnt
						if (cpntnam(j) .eq. species(c)) irxnic(i, 2) = j
					enddo
				endif
			elseif (states(a) .eq. 's') then
				nimsp(i) = 2
				do j = 1, nimm
					if (immnam(j) .eq. species(a)) irxnim(i, 1) = j
				enddo
				if (c .eq. 0) then
					nimsp(i) = 1
				else
					do j = 1, nimm
						if (immnam(j) .eq. species(c)) irxnim(i, 2) = j
					enddo
				endif
			elseif (states(k) .eq. 'g') then
				nivsp(i) = 2
				do j = 1, nvap
					if (vapnam(j) .eq. species(a)) irxniv(i, 1) = j
				enddo
				if (c .eq. 0) then
					nivsp(i) = 1
				else
					do j = 1, nvap
						if (vapnam(j) .eq. species(c)) irxniv(i, 2) = j
					enddo
				endif
			else
				write(ierr, '(a)') 'Internal error setting decay reaction state.'
				goto 667
			endif
		elseif ((rxntypes(i) .eq. 7) .or. (rxntypes(i) .eq. 8)) then
			nimsp(i) = 1
			naqsp(i) = 0
			nivsp(i) = 0
			a = 0
			do j = 1, nspecies
				if ((reactants(i, 1) .eq. species(j)) .or. (reactants(i, 1) .eq. masters(j)) .and. &
					(reactants(i, 1) .ne. '*')) then
					a = j
					exit
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a, i3, a, i1, a)') 'Error:  Reactant "'//trim(reactants(i, 1))// &
					'" in type ', rxntypes(i), ' reaction ', i, ' not found in comp.'
				goto 666
			endif
			if (states(a) .eq. 's') then
				if (reactants(i, 2) .ne. '*') then
					goto 212
				endif
				do j = 1, nimm
					if (pimm(j) .eq. a) then
						irxnim(i, 1) = j
						goto 257
					endif
				enddo
				write(ierr, '(a)') 'Internal error setting immoble species for reaction type 7/8.'
				goto 667
257				pdstim = stoichiometries(i, 1, 1)
				do j = 1, SPEC_MAX
					if (products(i, j) .eq. '*') then
						exit
					endif
					do k = 1, nspecies
						if (products(i, j) .eq. masters(k)) then
							if ((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
								write(ierr, '(a, i1, a, i3, a)') 'Error:  One side of type ', &
									rxntypes(i), ' reaction ', i, ' must consist '// &
									'entirely of aqueous master species.'
								goto 666
							endif
							naqsp(i) = naqsp(i) + 1
							do o = 1, ncpnt
								if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
							enddo
							pdstic(i, naqsp(i)) = stoichiometries(i, 2, j)
							goto 213
						endif
					enddo
					write(ierr, '(a, i1, a, i3, a)') 'Error:  Product "'//trim(products(i, j))// &
						'" in type ', rxntypes(i), ' reaction ', i, ' not found in comp.'
					goto 666
213					continue
				enddo
			else
				c = 0
				do j = 1, nspecies
					if ((products(i, 1) .eq. species(j)) .or. (products(i, 1) .eq. masters(j)) .and. &
						(products(i, 1) .ne. '*')) then
						c = j
						exit
					endif
				enddo
				if (c .eq. 0) then
					write(ierr, '(a, i1, a, i3, a)') 'Error:  Product "'//trim(products(i, 1))// &
						'" in type ', rxntypes(i), ' reaction ', i, ' not found in comp.'
					goto 666
				endif
				if (states(c) .eq. 's') then
					if (products(i, 2) .ne. '*') then
						goto 212
					endif
					do j = 1, nimm
						if (pimm(j) .eq. c) then
							irxnim(i, 1) = j
							goto 258
						endif
					enddo
					write(ierr, '(a)') 'Internal error setting immoble species for reaction type 7/8.'
					goto 667
258					pdstim = stoichiometries(i, 2, 1) ! pdstim(i, 1)
					do j = 1, SPEC_MAX
						if (reactants(i, j) .eq. '*') then
							exit
						endif
						do k = 1, nspecies
							if (reactants(i, j) .eq. masters(k)) then
								if ((states(k) .ne. 'a') .and. (states(k) .ne. 'h')) then
									write(ierr, '(a, i1, a, i3, a)') 'Error:  One '// &
										'side of type ', rxntypes(i), &
										' reaction ', i, ' must consist '// &
										'entirely of aqueous master species.'
									goto 666
								endif
								naqsp(i) = naqsp(i) + 1
								do o = 1, ncpnt
									if (cpntnam(o) .eq. species(k)) irxnic(i, naqsp(i)) = o
								enddo
								pdstic(i, naqsp(i)) = stoichiometries(i, 1, j)
								goto 227
							endif
						enddo
						write(ierr, '(a, i1, a, i3, a)') 'Error:  Reactant "'// &
							trim(reactants(i, j))//'" in type ', rxntypes(i), &
							' reaction ', i, ' not found in comp.'
						goto 666
227						continue
					enddo
				else
					goto 212
				endif
			endif
			goto 232
212			write(ierr, '(a, i1, a, i3, a)') 'Error:  One side of type ', rxntypes(i), ' reaction ', i, &
				' must contain only a single solid species.'
			goto 666
232			read(dscoefs(i), *, err=228) r
			ckeqlb(i) = r
			goto 229
228			a = 0
			do j = 1, nsolmodels
				if (dscoefs(i) .eq. solmodelnames(j)) then
					temp_model_kin(i) = 'l'
					a = j
					exit
				endif
			enddo
			if (a .eq. 0) then
				write(ierr, '(a)') 'Internal error setting distribution coefficient information.'
				goto 667
			endif
			j = 0
			do
				if (solmodels(a, j + 1, 1) .eq. -1) exit
				j = j + 1
				rarray1(j) = solmodels(a, j, 1);
				rarray2(j) = solmodels(a, j, 2);
			enddo
			call lstsq(j, rarray2, rarray1)
			do j = 1, 3
				tcoeff(i, j) = rarray2(i)
			enddo
			ckeqlb(i) = tcoeff(i, 1) + tcoeff(i, 2) * 25 + tcoeff(i, 3) * 625
229			ckmtrn(i) = rates(i)
			if (rxntypes(i) .eq. 8) then
				mw_mineral(irxnim(i, 1)) = porchange(i, 1)
				rho_mineral(irxnim(i, 1)) = porchange(i, 2)
			endif
			sarea(i) = xcoef(i)
		else
			write(ierr, '(a)') 'Internal error setting reaction types.'
			goto 667
		endif
	enddo
	! End of rxn value-setting.
	! End of FEHM compatibility block.
	goto 6000

600	write(ierr, '(a)') 'Error:  Generic read error while reading input in trxn.'
	goto 666
602	write(ierr, '(a)') 'Error:  Generic read error in trxn postprocessing.'
	goto 666

	! Go to 666 for user errors
666	if (reading) then
		write(ierr, '(a)') 'Error in trxn during main read.'
		write(ierr, '(a)') ' >> '//trim(line)//' <<'
		write(iptty, '(a)') 'Error in trxn during main read.'
		write(iptty, '(a)') ' >> '//trim(line)//' <<'
	else
		write(ierr, '(a)') 'Error in trxn in postprocessing.'
		write(iptty, '(a)') 'Error in trxn in postprocessing.'
	endif
	stop
	! Go to 667 for internal errors
667	write(iptty, '(a)') 'An internal error has occurred in trxn.  This should never happen, '// &
		'and indicates a bug in the program.  Please contact fehm-dev@lanl.gov with information about the error.'
	stop

	! Perform cleanup if no errors have occurred
	! Deallocate local arrays
6000	if (allocated(a1l)) deallocate(a1l)
	if (allocated(a1v)) deallocate(a1v)
	if (allocated(a2l)) deallocate(a2l)
	if (allocated(a2v)) deallocate(a2v)
	if (allocated(ah)) deallocate(ah)
	if (allocated(ah1)) deallocate(ah1)
	if (allocated(ah2)) deallocate(ah2)
	if (allocated(ah3)) deallocate(ah3)
	if (allocated(ah4)) deallocate(ah4)
	if (allocated(ah5)) deallocate(ah5)
	if (allocated(bioparams)) deallocate(bioparams)
	if (allocated(bl)) deallocate(bl)
	if (allocated(bv)) deallocate(bv)
	if (allocated(cdenspecs)) deallocate(cdenspecs)
	if (allocated(complexcontents)) deallocate(complexcontents)
	if (allocated(complexes)) deallocate(complexes)
	if (allocated(complexstoich)) deallocate(complexstoich)
	if (allocated(datcden)) deallocate(datcden)
	if (allocated(datcomp)) deallocate(datcomp)
	if (allocated(datcplx)) deallocate(datcplx)
	if (allocated(datcplxheq)) deallocate(datcplxheq)
	if (allocated(datcplxlkeq)) deallocate(datcplxlkeq)
	if (allocated(datcplxmain)) deallocate(datcplxmain)
	if (allocated(datcplxstoic)) deallocate(datcplxstoic)
	if (allocated(datcplxtemp)) deallocate(datcplxtemp)
	if (allocated(datmaster)) deallocate(datmaster)
	if (allocated(datmin)) deallocate(datmin)
	if (allocated(datminheq)) deallocate(datminheq)
	if (allocated(datminlkeq)) deallocate(datminlkeq)
	if (allocated(datminmain)) deallocate(datminmain)
	if (allocated(datminnam)) deallocate(datminnam)
	if (allocated(datminstoic)) deallocate(datminstoic)
	if (allocated(datmintemp)) deallocate(datmintemp)
	if (allocated(dhh)) deallocate(dhh)
	if (allocated(dispnames)) deallocate(dispnames)
	if (allocated(distmodelnames)) deallocate(distmodelnames)
	if (allocated(distmodels)) deallocate(distmodels)
	if (allocated(dscoefs)) deallocate(dscoefs)
	if (allocated(ecomplexes)) deallocate(ecomplexes)
	if (allocated(enthalpies)) deallocate(enthalpies)
	if (allocated(eqconsts)) deallocate(eqconsts)
	if (allocated(eqtemps)) deallocate(eqtemps)
	if (allocated(equconstants)) deallocate(equconstants)
	if (allocated(guesses)) deallocate(guesses)
	if (allocated(heqfit)) deallocate(heqfit)
	if (allocated(hh)) deallocate(hh)
	if (allocated(hhtspecies)) deallocate(hhtspecies)
	if (allocated(hmodels)) deallocate(hmodels)
	if (allocated(icbioholder)) deallocate(icbioholder)
	if (allocated(ll)) deallocate(ll)
	if (allocated(lsorptypes)) deallocate(lsorptypes)
	if (allocated(lt)) deallocate(lt)
	if (allocated(lx)) deallocate(lx)
	if (allocated(ly)) deallocate(ly)
	if (allocated(lz)) deallocate(lz)
	if (allocated(masters)) deallocate(masters)
	if (allocated(molegrid)) deallocate(molegrid)
	if (allocated(molweights)) deallocate(molweights)
	if (allocated(mspecs)) deallocate(mspecs)
	if (allocated(mzones)) deallocate(mzones)
	if (allocated(porchange)) deallocate(porchange)
	if (allocated(printspecies)) deallocate(printspecies)
	if (allocated(products)) deallocate(products)
	if (allocated(rates)) deallocate(rates)
	if (allocated(rates2)) deallocate(rates2)
	if (allocated(reactants)) deallocate(reactants)
	if (allocated(rtgrid)) deallocate(rtgrid)
	if (allocated(rtnames)) deallocate(rtnames)
	if (allocated(rtspecies)) deallocate(rtspecies)
	if (allocated(rxnaqueous)) deallocate(rxnaqueous)
	if (allocated(rxnsolid)) deallocate(rxnsolid)
	if (allocated(rxntypes)) deallocate(rxntypes)
	if (allocated(rxnvapor)) deallocate(rxnvapor)
	if (allocated(solmodelnames)) deallocate(solmodelnames)
	if (allocated(solmodels)) deallocate(solmodels)
	if (allocated(sorpnames)) deallocate(sorpnames)
	if (allocated(sorpspecs)) deallocate(sorpspecs)
!	if (allocated(species)) deallocate(species)
	if (allocated(spnam)) deallocate(spnam)
	if (allocated(states)) deallocate(states)
	if (allocated(stoichiometries)) deallocate(stoichiometries)
	if (allocated(vl)) deallocate(vl)
	if (allocated(vsorptypes)) deallocate(vsorptypes)
	if (allocated(vt)) deallocate(vt)
	if (allocated(vtgrid)) deallocate(vtgrid)
	if (allocated(vtnames)) deallocate(vtnames)
	if (allocated(vtspecies)) deallocate(vtspecies)
	if (allocated(vx)) deallocate(vx)
	if (allocated(vy)) deallocate(vy)
	if (allocated(vz)) deallocate(vz)
	if (allocated(wtgrid)) deallocate(wtgrid)
	if (allocated(wtnames)) deallocate(wtnames)
	if (allocated(wtspecies)) deallocate(wtspecies)
	if (allocated(xcoef)) deallocate(xcoef)
	if (allocated(zonegrid)) deallocate(zonegrid)
	if (allocated(zoneresolv)) deallocate(zoneresolv)
	if (allocated(zones)) deallocate(zones)
	if (allocated(zonespecs)) deallocate(zonespecs)

	! Finish up
	macroread(5) = .true.
	if (debug) write(iptty, '(a)') 'trxn read successfully.'
	if (debug_stop) then
		write(iptty, '(a)') 'Stop requested:  stopping.'
		stop
	endif
	return ! All done

8000	end subroutine rdtr
