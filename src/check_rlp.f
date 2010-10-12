      subroutine check_rlp
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
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
!D1 Test calculation of relative permeabilities and capillary pressures
!D1 and derivatives and output values.
!D1
!***********************************************************************
      
      use comai, only : form_flag, idpdp, ierr, neq, nrlp, wdd
      use comrlp, only : ishisrlp, rlpnew, delta_sat, num_sat, sat_out
      use comci, only : rlf, rvf
      use comdi, only : ieos, irlp, pcp, s

      implicit none
      integer i, j, ndummy, neqtemp
      integer, allocatable :: ieostemp(:)
      integer, allocatable :: irlptemp(:)
      real*8, allocatable  :: stemp(:)
      character*100 form_string, title_string
      
      neqtemp = neq
      if (num_sat .eq. 0) then
         neq = 1 / delta_sat + 1
      else
         neq = num_sat
      end if
      if (idpdp .ne. 0) then
         allocate (stemp(2*neq), ieostemp(2*neq), irlptemp(2*neq))
      else
         allocate (stemp(neq), ieostemp(neq), irlptemp(neq))
      end if
      ndummy = 0

      do i = 1, neq
         stemp(i) = s(i)
         ieostemp(i) = ieos(i)
         irlptemp(i) = irlp(i)
         if (idpdp .ne. 0) then
            stemp(i + neq) = s(i + neqtemp)
            ieostemp(i + neq) = ieos(i + neqtemp)
            irlptemp(i + neq) = irlp(i + neqtemp)
         end if
         if (num_sat .eq. 0) then
            s(i) = (i - 1) * delta_sat
         else
            s(i) = sat_out(i)
         end if
         ieos(i) = 2
         if (idpdp .ne. 0) then
            s(i + neq) = s(i)
            ieos(i + neq) = 2
         end if
      end do

      title_string = "Relative permeability and " //
     &     "Capillary pressure"
      if (form_flag .eq. 1) then
         form_string = 'variables = "Saturation" ' //
     &        '"Liquid" "Vapor" "Capillary pressure"'
         write(ishisrlp, '(a)') trim(form_string)
         write(ishisrlp, 230) 50., 95., trim(title_string)
         write(ishisrlp, 230) 50., 90., trim(wdd)
      else if (form_flag .eq. 2) then
         form_string = 'Saturation, Liquid, Vapor, ' //
     &        'Capillary pressure'
         write(ishisrlp, '(a)') trim(form_string)
      else
         form_string = '"Saturation" "Liquid" "Vapor" ' //
     &        '"Capillary pressure"'
         write(ishisrlp, '(a)') trim(title_string)
         write(ishisrlp, '(a)') trim(wdd)
         write(ishisrlp, '(a)') trim(form_string)
      end if
         
      do j = 1, nrlp
         do i = 1, neq
            irlp(i) = j
            if (idpdp .ne. 0) irlp(i+neq) = j
         end do
         if (form_flag .eq. 1) then
            write (ishisrlp, 240) j
         else if (form_flag .eq. 2) then
         else
         end if              
         if (rlpnew) then
            call rlp_cap(0)
            if (idpdp .ne. 0) call rlp_cap(neq)
         else
            call rlperm(0,1)
            if (idpdp .ne. 0) call rlperm(neq,1) 
         end if
         do i = 1, neq
            write (ishisrlp, '(4(g16.9, 1x))') s(i), rlf(i), rvf(i),
     &           pcp(i)
         end do
         if (idpdp .ne. 0) then
            if (form_flag .eq. 1) then
               write (ishisrlp, 245) j
            else if (form_flag .eq. 2) then
            else
            end if              
            do i = neq + 1, 2*neq
               write (ishisrlp, '(4(g16.9, 1x))') s(i), rlf(i), rvf(i),
     &              pcp(i)
            end do
         end if
      end do

      do i = 1, neq
         s(i) = stemp(i)
         ieos(i) = ieostemp(i)
         irlp(i) = irlptemp(i)
      end do
      neq = neqtemp

      deallocate (stemp, ieostemp, irlptemp)
      close (ishisrlp)

 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')
 240  format ('zone t = "model ', i4, '"')
 245  format ('zone t = "model ', i4, ' matrix"')

      end subroutine check_rlp
