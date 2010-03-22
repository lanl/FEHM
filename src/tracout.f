      subroutine tracout
!***********************************************************************
! Copyright 2007 Los Alamos National Security, LLC  All rights reserved
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
!D1 Setup files for tracer history output.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 15-May-07,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/tracout.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai
      use comchem
      use comdi
      use comxi
      implicit none

      integer i, ic, im, iv, ix, iroot, iunit
      logical null1
      character*3 spnum
      character*4 trcsfx
      character*14 time_string, units
      character*20 tmpname
      character*80 form_string, title_string
      character*14, allocatable :: var_string(:)
      character*150 fname, root

      if (null1(root_name)) then
! Find root of tracer file name or if not present, output or input file name
         if (nmfil(9) .ne. nmfily(3) .and. nmfil(9) .ne. ' ') then
            call file_prefix(nmfil(9), iroot)
            if (iroot .gt. 100) iroot = 100
            root(1:iroot) = nmfil(9)(1:iroot)
         else
            if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
               call file_prefix(nmfil(5), iroot)
               if (iroot .gt. 100) iroot = 100
               root(1:iroot) = nmfil(5)(1:iroot)
            else
               if (nmfil(2)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil2 file: ', nmfil(2),
     .                 ' unable to determine history file prefix'
                  stop
               else
                  call file_prefix(nmfil(2), iroot)
                  if (iroot .gt. 100) iroot = 100
                  root(1:iroot) = nmfil(2)(1:iroot)
               end if
            end if
         endif
      else
         iroot = len_trim (root_name)
         if (iroot .gt. 100) iroot = 100
         root(1:iroot) = root_name(1:iroot)
      end if

      title_string = ''
      time_string = ''

      if (time_flag .eq. 1) then
         time_string = 'Time (years)'
      else if (time_flag .eq. 2) then
         time_string = 'Time (days)'
      else if (time_flag .eq. 3) then
         time_string = 'Time (seconds)'
      else if (time_flag .eq. 4) then
         time_string = 'Time (hours)'
      end if

      allocate (var_string(m))
      do i = 1, m
         var_string(i)(1:4) = 'Node'  
         write (var_string(i)(5:14), '(i9)') nskw(i)
      end do
      form_string = ''
      select case (form_flag)
      case (1)
c Tecplot
         write (form_string, 200) m
         trcsfx = '.dat'
      case (2)
c Surfer (csv)
         write (form_string, 210) m
         trcsfx = '.csv'
      case default
c Standard text
         write (form_string, 220) m
         trcsfx = '.trc'
      end select

      do i = 1, ncpntprt
         fname = ''
         title_string = ''
         ic = cpntprt(i)
         nsp =  pcpnt(ic)
         if (icns(i) .ne. -2) then
            units = 'Moles/kg water'
         else
            units = 'Moles/kg vapor'
            tmpname = cpntnam(ic)
            if (tmpname(1:7) .eq. 'Aqueous') then
               tmpname(1:5) = 'Vapor'
               tmpname(6:17) = tmpname(8:19)
               tmpname(18:20) = ''
               cpntnam(ic) = tmpname
            end if
         end if
         title_string = cpntnam(ic)
         iunit = ishisc + nsp
         write(spnum, '(i3.3)') nsp
         fname =  root(1:iroot) // '_' // trim(cpntnam(ic)) // trcsfx
         open (unit=iunit, file=fname, form='formatted')
         call write_fileheader
      end do
      units = 'Moles/kg rock'
      do i = 1, nimmprt
         fname = ''
         title_string = ''
         im = immprt(i)
         nsp = pimm(im)
         title_string = immnam(im)  
         iunit = ishisc + nsp
         write(spnum, '(i3.3)') nsp
         fname =  root(1:iroot) // '_' // trim(immnam(im)) // trcsfx
         open (unit=iunit, file=fname, form='formatted')
         call write_fileheader
      end do
      units = 'Moles/kg vapor'
      do i = 1, nvapprt
         fname = ''
         title_string = ''
         iv = vapprt(i)
         nsp = pvap(iv)
         title_string = vapnam(iv)     
         iunit = ishisc + nsp
         write(spnum, '(i3.3)') nsp
         fname =  root(1:iroot) // '_' // trim(vapnam(iv)) // trcsfx
         open (unit=iunit, file=fname, form='formatted')
         call write_fileheader
      end do
      if(ncplx.ne.0)then
         units = 'Moles/kg water'
         do i = 1,ncpntprt
            fname = ''
            title_string = ''
            ic = cpntprt(i)
            nsp =  pcpnt(ic)
            title_string = 'Free Ion ' // trim(cpntnam(ic))
            iunit = ishisc + ic + nspeci
            write(spnum, '(i3.3)') nsp
            fname =  root(1:iroot) // '_FreeIon_' // trim(cpntnam(ic)) 
     &            // trcsfx
            open (unit=iunit, file=fname, form='formatted')
            call write_fileheader
         enddo
         do i = 101, ncplxprt+100
            fname = ''
            title_string = ''
            ix = cplxprt(i) 
            title_string = cplxnam(ix)
            nsp =  ix
            iunit = ishisc + nsp
            write(spnum, '(i3.3)') nsp
            fname =  root(1:iroot) // '_' // trim(cplxnam(ix)) // trcsfx
            open (unit=iunit, file=fname, form='formatted')
            call write_fileheader
         enddo
      end if
 200  format ("('variables = ',", "'", '"', "', a, '", '" ', "',", 
     &     i5, "('", '"', "', a, '", '" ', "'))")
 210  format ("( a,", i5, '(", ", a))')
 220  format ("( a,", i5, "(1x, a))")

      contains
      subroutine write_fileheader

      implicit none
      integer j
      character*30 tmp_string

      select case (form_flag)
      case (1)
c Tecplot
         tmp_string = 'Concentration ' // '(' // trim(units) // ')'
         write(iunit, 250) verno, jdate, jtime
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(j)),j = 1, m)
         write(iunit, 230) 50., 95., trim(wdd)
c         write(iunit, 230) 50., 90., trim(tmp_string)
         write(iunit, 260) trim(title_string), tmp_string
      case (2)
c Surfer
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(j)), j = 1, m)
      case default
c Standard text
         write(iunit, 240) verno, jdate, jtime, trim(wdd)
         write(iunit, 270) trim(title_string), trim(units)
         write(iunit, form_string) trim(time_string), 
     &        (trim(var_string(j)), j = 1, m)
      end select

 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')
 240  format(a30, 2x, a10, 2x, a8, /, a)
 250  format('TITLE = "', a30, 2x, a10, 2x, a8,'"')
 260  format('ZONE T = "', a, 1x, a, '"')
 270  format(a, ' Concentration ', '(', a, ')')
      end subroutine write_fileheader

      end subroutine tracout
