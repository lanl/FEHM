      subroutine  cfluxz(ptime)
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
!D1 Output concentration flux for a zone.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 26-Jan-07,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/cfluxz.f_a  $
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
      use combi
      use comci, only : rolf, rovf
      use comdi
      use comdti
      use comflow
      use comchem
      use compart, only : pcnsk
      use comxi
      use davidi, only : irdof

      implicit none
      integer addnode, iconn, idummy, i1, i2, ipr_vapor
      integer i, indexa_axy, inneq, inode, inodec, izone, md
      integer num_aq, num_v
      integer open_file, iroot, is, ie, iname
      integer, allocatable :: icfile(:)
      real*8 ptime, sumfout, sumsink, sumsource, sumboun, sumfin,sum_vap
      real*8 ancell, anlocell, sumcell, sumcello, sumdelta, cellh2o
      character*80 string, formstring
      character*85, allocatable :: flux_string(:)
      character*150 cflx_name, cflx_root, string2
      logical matrix_node, null1
      save icfile

! Fluxes are written to concentration flux history file
      if (.not. allocated (flux_string)) 
     &     allocate (flux_string(cflxz))
      flux_string = ''

      if (.not. allocated(icfile)) then
         allocate (icfile(nspeci))
         if (null1(root_name)) then
! Use trc file root name
            if (nmfil(9) .ne. nmfily(3) .and. nmfil(9) .ne. ' ') 
     &           then
               call file_prefix(nmfil(9), iroot)
               if (iroot .gt. 100) iroot = 100
               cflx_root(1:iroot) = nmfil(9)(1:iroot)
            else
               if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') 
     &              then
                  call file_prefix(nmfil(5), iroot)
                  if (iroot .gt. 100) iroot = 100
                  cflx_root(1:iroot) = nmfil(5)(1:iroot)
               else
                  if (nmfil(2)(1:1) .eq. ' ' ) then
                     write (ierr, *) 'FILE ERROR: nmfil2 file: ', 
     &                    nmfil(2),
     &                    ' unable to determine cflx file prefix'
                     stop
                  else
                     call file_prefix(nmfil(2), iroot)
                     if (iroot .gt. 100) iroot = 100
                     cflx_root(1:iroot) = nmfil(2)(1:iroot)
                  end if
               end if
            endif
         else
            iroot = len_trim (root_name)
            if (iroot .gt. 100) iroot = 100
            cflx_root(1:iroot) = root_name(1:iroot)
         end if
         string = '# Zone Flux (moles/day): '
         is = 26
         ie = 35
         do i = 1, 5
            if (cflx_var(i)) then
               select case (i)
               case (1)
                  write(string(is:ie), '(a, x)') 'Source'
               case (2)
                  write(string(is:ie), '(a, x)') 'Sink'
               case (3)
                  write(string(is:ie), '(a, x)') 'NetIn'
               case (4)
                  write(string(is:ie), '(a, x)') 'NetOut'
               case (5)
                  write(string(is:ie), '(a, x)') 'Boundary'
               end select
               is = ie + 1
               ie = ie + 9
            end if
         end do
         string2 = '# Time (days) Zones:'
         is = 21
         ie = 26
         do i = 1, cflxz
            write (string2(is:ie), '(1x, i5)') icflxz(i)
            is = ie + 1
            ie = ie + 6
         end do

! Start counter for number of aqueous or vapor species ABJ 3/2012
         num_aq=1
         num_v=1
         
         do nsp = 1, nspeci

         cflx_name = ''
         cflx_name(1:iroot) = cflx_root(1:iroot)

! Check if an aqueous or vapor species filename is needed ABJ 3/2012
            if (icns(nsp).ge.0.or.abs(icns(nsp)).eq.2) then
                iname = len_trim (cpntnam(num_aq))
                cflx_name(iroot+1:iroot+iname) = cpntnam(num_aq)
                num_aq = num_aq + 1
            else
                iname = len_trim (vapnam(num_v))
                cflx_name(iroot+1:iroot+iname) = vapnam(num_v)
                num_v = num_v + 1
            end if
            cflx_name(iroot+iname+1:iroot+iname+5) = '.cflx'
            icfile(nsp) = open_file(cflx_name, 'unknown')
            write(icfile(nsp), '(a)') trim(string)
            write(icfile(nsp), '(a)') trim(string2)
         end do
         formstring = ''
         write (formstring, 200) cflxz
      end if

 100  format ('("#Time (days) Species# Zones:",', i3, '(1x, i3))')
 104  format ('(a19,', i3, '(1x, i3))')
 200  format ("(g16.9, ",i4,"(1x, a))")

c     Compute fluxes out of zone (>0), leaving out flues
c     into other parts of the zone

      do nsp = 1, nspeci
         npn = npt(nsp)
         do izone = 1, cflxz
            sumfout = 0.
            sumfin  = 0.
            sumsink = 0.
            sumsource = 0.
            sumboun = 0.0
            sum_vap = 0.0
            sumcell = 0.
            sumcello = 0.
            md=0
c     Loop over all nodes
            do inode = 1, n0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
               if(inode.gt.neq) then
                  inneq = inode-neq
                  matrix_node = .true.
                  addnode = nelm(neq+1)-neq-1
                  idummy = neq
               else
                  matrix_node = .false.
                  inneq = inode
                  addnode = 0
                  idummy = 0
               end if
c     Determine if node is part of the zone being summed
               if(izoncflxz(inode).eq.izone) then
                  md = md+1
                  inodec = inode + npn
c     Add boundary condition sources
                  if (pcnsk(inodec) .eq. -1) then
                     sumboun=sumboun + rcss(inodec)
                  else if (sk(inode) .lt. 0.) then
! Incoming
                     sumboun=sumboun + sk(inode)*cnsk(inodec)
                  else
                     sumboun=sumboun + sk(inode)*an(inodec)
                  end if
c     Compute mass in cell if liquid (moles/kg H2O *  kg / m^3 H2O * m^3 cell)
                  if (icns(nsp).ge.0.or.abs(icns(nsp)).eq.2) then
                     cellh2o = rolf(inode) * sx1(inode) * ps(inode) 
     &                    * s(inode)
                     ancell = an(inodec) * cellh2o
                     anlocell = anlo(inodec) * cellh2o

!     Compute mass in cell if vapor (moles/kg air * kg/m^3 air * m^3 cell)
!     Still called cellh2o, however  ABJ 3/2012
                  else
                     cellh2o = rovf(inode) * sx1(inode) * ps(inode)
     &                    * (1-s(inode))
                     ancell = an(inodec) * cellh2o
                     anlocell = anlo(inodec) * cellh2o
                  end if

                  sumcell = sumcell + ancell
                  sumcello = sumcello + anlocell
                  if (rcss(inodec) .gt. 0.) then
                     sumsink = sumsink + rcss(inodec)
                  else if (rcss(inodec) .lt. 0.) then
                     sumsource = sumsource + rcss(inodec)
                  end if
c     calculate vapor out (assume zero if in at this time)
!               if(irdof.ne.13.or.ifree.ne.0) then
!                  if(sk(inode).gt.0.0d0) then
!                     sum_vap = sum_vap + (1.0-s(inode))*sk(inode)*cnsk(inodec)
!                  endif
!               endif
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
c                  i1 = nelm(inneq)+1
c                  i2 = nelm(inneq+1)
c Skip this, we are summing total moles in zone to find change
                  i1 = 1
                  i2 = 0
c     loop over every connecting node
                  do iconn = i1, i2
                     indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
                     if(a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoncflxz(idummy+nelm(iconn))
     2                       .ne.izone .or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfout = sumfout + a_axy(indexa_axy)*
     &                          an(inodec)
                           if(nelm(iconn).eq.inneq) then
                              sumsink = sumsink +
     2                             a_axy(indexa_axy)*an(inodec)
                           end if
                        end if
                     elseif(a_axy(indexa_axy).lt.0.) then
                        if(izoncflxz(idummy+nelm(iconn))
     2                       .ne.izone .or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfin = sumfin + a_axy(indexa_axy)*
     &                          an(inodec)
c     add to source sum
                           if(nelm(iconn).eq.inneq) then
                              sumsource = sumsource +
     2                             a_axy(indexa_axy)*cnsk(inodec)
                           end if
                        end if
                     end if
                  end do
               end if
            end do
            sumdelta = (sumcell - sumcello) / dtotc
            sumfin = sumsource
            sumfout = sumsink
            if (sumdelta .gt. 0.) then
               sumfin = sumfin + sumdelta
            else
               sumfout = sumfout +sumdelta
            end if
c     Write results
! Fluxes (moles/day) are written to concentration flux history file
            flux_string(izone) = ''
            is = 1
            ie = 17
            do i = 1, 5
               if (cflx_var(i)) then
                  select case (i)
                  case (1)
                     write (flux_string(izone)(is:ie), 1050) sumsource *
     &                    86400.
                  case (2)
                     write (flux_string(izone)(is:ie), 1050) sumsink *
     &                    86400.
                  case (3)
                     write (flux_string(izone)(is:ie), 1050) sumfin *
     &                    86400.
                  case (4)
                     write (flux_string(izone)(is:ie), 1050) sumfout *
     &                    86400.
                  case (5)
                     write (flux_string(izone)(is:ie), 1050) sumboun *
     &                    86400.
                  end select
                  is = ie + 1
                  ie = ie + 17
               end if
            end do
 1050       format(1x, g16.9)
         end do
         
! Write to flux history file
c         write (icfile(nsp), formstring) ptime, 
c     &        (trim(flux_string(izone)), izone = 1, cflxz)
         write (icfile(nsp), '(g16.9, 2(1x, a))') ptime, 
     &        (trim(flux_string(izone)), izone = 1, cflxz)
         call flush(icfile(nsp))
      end do

      if (days .ge. tims) deallocate (icfile)

      end
