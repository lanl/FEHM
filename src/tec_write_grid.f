      subroutine tec_write_grid(lu)
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Output TECPLOT grid information for FEHM.
CD1
C***********************************************************************
CD2
CD2 Initial implementation: 12-Oct-11,  Programmer: Z. Dash
CD2
C***********************************************************************
CD3
CD3 This format for grid and solution output requires Tecplot 360 2010
CD3 or later.
CD3
C***********************************************************************

      use avsio
      use comai, only : icnl, ierr, ifdm_elem, ivf, neq_primary, 
     &     nei_in, ns_in, wdd
      use combi, only : cord, nelm, elem_geo, nact_elem_geo
      implicit none

      integer i, j, lu, i1, il, open_file
      integer, allocatable :: nelm2(:)
      character(18) :: title(3), typestring
      character(60) :: varstring
      real*8 dummy

c
c     ASCII formatted output
c
c     if(icnl2 .eq. 0) => three dimensional problem
c     else            => two   dimensional problem
c
      title(1) = '"X coordinate (m)"'
      title(2) = '"Y coordinate (m)"'
      title(3) = '"Z coordinate (m)"'

      if (ivf .eq. -1 .and. ifdm_elem .eq. 1) then
c first generate elements      
         call structured(4)
         il = open_file('fdm_elem.macro','old')
         read(il,*) 
         read(il,*)  ns_in , nei_in
         if (lu .eq. 0) close (il)
      end if

      select case (icnl)
      case (0)
         varstring = title(1) // " " // title(2) // " " // title(3)
         if (ns_in .eq. 4) then
            typestring = 'FETETRAHEDRON'
         else if (ns_in .eq. 0) then
            typestring = ''
         else
            typestring = 'FEBRICK'
         end if
      case (1, 4)
         varstring = title(1) // " " // title(2)
      case (2, 5)
         varstring = title(1) // " " // title(3)
      case (3, 6)
         varstring = title(2) // " " // title(3)
      end select

      if (icnl .ne. 0) then
         select case (ns_in)
         case (2)
            typestring = 'FELINESEG'
         case (3)
            typestring = 'FETRIANGLE'
         case (4)
            typestring = 'FEQUADRILATERAL'
         end select
      end if
      
      if (nei_in .eq. 0) then
         gridstring = ''
      else
         write (gridstring, 40) neq_primary, nei_in,
     &        trim(typestring)
      end if

c     If lu =0, we are only defining the gridstring
      if (lu .eq. 0) return

      write (lu, 10) trim(wdd)
      write (lu, 20)
      write (lu, 30) trim(varstring)
      write (lu, 45) trim(gridstring)
      do i = 1, neq_primary
         select case (icnl)
         case (0)
            write (lu, 50) (cord(i, j), j = 1, 3)
         case (1, 4)
            write (lu, 50) (cord(i, j), j = 1, 2)
         case (2, 5)
            write (lu, 50) (cord(i, j), j = 1, 3, 2)
         case (3, 6)
            write (lu, 50) (cord(i, j), j = 2, 3)
         end select
      end do

      if (ivf .eq. -1 .and. ifdm_elem .eq. 1) then
        if(nact_elem_geo.eq.0) then
         allocate (nelm2(ns_in))
         do i = 1, nei_in
            read (il,*) i1, (nelm2(j), j=1,ns_in)
            write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
         end do
         deallocate(nelm2)
         close (il)
        else
         do i = 1, nei_in
            write (lu, 60) (elem_geo((i-1)*ns_in + j), j = 1, ns_in)
         end do   
        endif    
      else if(nact_elem_geo.eq.1) then
         do i = 1, nei_in
            write (lu, 60) (elem_geo((i-1)*ns_in + j), j = 1, ns_in)
         end do
      end if

 10   format ('TITLE = "', a, '"')
 20   format ('FILETYPE = "GRID"')
 30   format ('VARIABLES = ', a)
 40   format (', N = ', i8, ', E = ', i8, 
     &     ', DATAPACKING = POINT', ', ZONETYPE = ', a)
 45   format ('ZONE T = "GRID"', a, ', STRANDID = 0, ',
     &     'SOLUTIONTIME = 0.')
 50   format(3(e16.9,2x))
 60   format(8(i8))
 
      end subroutine tec_write_grid

