      subroutine check_rlp_carb
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
      use comco2
      use comdi, only : ieos, irlp, pcp, icap

      implicit none
      integer i, j, mi, ndummy, neqtemp
      integer, allocatable :: ieostemp(:)
      integer, allocatable :: irlptemp(:)
      integer, allocatable :: icaptemp(:)
      real*8, allocatable  :: fwtemp(:)
      real*8, allocatable  :: fltemp(:)
      real*8, allocatable  :: fgtemp(:)
      real*8, allocatable  :: pcptemp(:)
      real*8, allocatable  :: pcgtemp(:)
      real*8, allocatable  :: rlwtemp(:)
      real*8, allocatable  :: rlltemp(:)
      real*8, allocatable  :: rlvtemp(:)
      real*8 :: dum1, dum2, dum3
      character*100 form_string, title_string
      
      neqtemp = neq
      if (idpdp .ne. 0) then
         allocate (fwtemp(2*neq), fltemp(2*neq), fgtemp(2*neq))
         allocate (ieostemp(2*neq), irlptemp(2*neq), icaptemp(2*neq))
         allocate (pcptemp(2*neq), pcgtemp(2*neq))
         allocate (rlwtemp(2*neq), rlltemp(2*neq), rlvtemp(2*neq))
      else
         allocate (fwtemp(neq), fltemp(neq), fgtemp(neq))
         allocate (ieostemp(neq), irlptemp(neq), icaptemp(neq))
         allocate (pcptemp(neq), pcgtemp(neq))
         allocate (rlwtemp(neq), rlltemp(neq), rlvtemp(neq))
      end if

      ndummy = 0

      fwtemp = fw
      fltemp = fl
      fgtemp = fg
      ieostemp = ieos
      irlptemp = irlp
      icaptemp = icap
      pcptemp = pcp
      pcgtemp = pcg
      rlwtemp = rl_w
      rlltemp = rl_l
      rlvtemp = rl_v

      if (num_sat .eq. 0) then
         neq = 1 / delta_sat + 1
      else
         neq = num_sat
      end if

      do i = 1, neq
         if (num_sat .eq. 0) then
            fw(i) = (i - 1) * delta_sat
            fl(i) = 1 - fw(i)
         else
            fw(i) = sat_out(i)
            fl(i) = 1 - fw(i)
         end if
         ieos(i) = 2
         if (idpdp .ne. 0) then
            fw(i + neq) = fw(i)
            fl(i + neq) = fl(i)
            ieos(i + neq) = 2
         end if
      end do
      
      title_string = "Relative permeability and " //
     &     "Capillary pressure"
      if (form_flag .eq. 1) then
         form_string = 'variables = "Saturation" ' //
     &        '"Liquid" "Vapor" "Capillary pressure" '
c     &       , '"drl_ww" "drl_lw" "drl_lg"'
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
         write(ishisrlp, '(a)') trim(form_string)
      end if
         
      do j = 1, nrlp
         do i = 1, neq
            irlp(i) = j
            icap(i) = j
            if (idpdp .ne. 0) then
               irlp(i+neq) = j
               icap(i+neq) = j
            end if
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
            do mi = 1, neq
c     calculate multi-phase relative perms.
               call rlperm_co2(0,0,mi,rl_w(mi),
     &              drl_ww(mi),drl_wg(mi),rl_l(mi),drl_lw(mi),
     &              drl_lg(mi),rl_v(mi),drl_vw(mi),drl_vg(mi))
c     calculate multi-phase cap. pres.
               call rlperm_co2(0,1,mi,pcp(mi),
     &              dpcpw(mi),dpcpg(mi),pcg(mi),dpcgw(mi),dpcgg(mi),
     &              dum1,dum2,dum3)
            end do
            if (idpdp .ne. 0) then
c                  call rlperm(neq,1)
c                  call cappr(1,neq)
            end if 
         end if
         do i = 1, neq
            write (ishisrlp, '(7(g16.9, 1x))') fw(i), rl_w(i), 
     &           rl_l(i), pcp(i)
c     &              rl_l(i), pcp(i), drl_ww(i), drl_lw(i), drl_lg(i)
         end do
         if (idpdp .ne. 0) then
            if (form_flag .eq. 1) then
               write (ishisrlp, 245) j
            else if (form_flag .eq. 2) then
            else
            end if              
            do i = neq + 1, 2*neq
               write (ishisrlp, '(7(g16.9, 1x))') fw(i), rl_w(i),
     &              rl_l(i), pcp(i) 
c     &              rl_l(i), pcp(i), drl_ww(i), drl_lw(i), drl_lg(i)
            end do
         end if
      end do

      fw = fwtemp
      fl = fltemp
      fg = fgtemp
      ieostemp = ieos
      irlptemp = irlp
      icaptemp = icap
      pcp = pcptemp
      pcg = pcgtemp
      rl_w = rlwtemp
      rl_l = rlltemp
      rl_v = rlvtemp
      neq = neqtemp

      deallocate (fwtemp, fltemp, fgtemp, ieostemp, irlptemp, icaptemp, 
     &     pcptemp, pcgtemp, rlwtemp, rlltemp, rlvtemp)
      close (ishisrlp)

 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')
 240  format ('zone t = "model ', i4, '"')
 245  format ('zone t = "model ', i4, ' matrix"')

      end subroutine check_rlp_carb
