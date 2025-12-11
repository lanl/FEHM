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
      
      use comai, only : form_flag, idpdp, ierr, neq, nrlp, wdd,
     & gdkm_flag,jdate,jtime,verno
      use comrlp, only : ishisrlp, rlpnew, delta_sat, num_sat, sat_out,
     &     rlp_type2, rlp_fparam,max_rp,max_cp,cap_param
      use comci, only : rlf, rvf
      use comdi, only : ieos, irlp, pcp, s, icap, irlpt
      use comxi, only : nmfil

      implicit none
      integer i, j, ndummy, neqtemp,mm
      integer, allocatable :: ieostemp(:)
      integer, allocatable :: irlptemp(:)
      integer, allocatable :: icaptemp(:)
      real*8, allocatable  :: stemp(:)
      real*8, allocatable  :: pcptemp(:)
      real*8, allocatable  :: rlftemp(:)
      real*8, allocatable  :: rvftemp(:)
c gaz 071423 s_test to deal with real to integer  conversion
      integer ii,jj,ishisrlpn
      real*8 s_test, s_comtol
      logical s_last
      logical op
      character*120  dum_tbl, input_file_string
c gaz 101023 changed character*120 to character*200
      character*200 fname, fname_new, model_string
      logical :: frac_model = .false.
      character*100 form_string, title_string
      parameter (s_comtol = 1.d-8)
      neqtemp = neq      
      if (idpdp .ne. 0) then
         allocate (stemp(2*neq), ieostemp(2*neq), irlptemp(2*neq))
         allocate (pcptemp(2*neq), rlftemp(2*neq), rvftemp(2*neq))
         allocate (icaptemp(2*neq))
      else
         allocate (stemp(neq), ieostemp(neq), irlptemp(neq))
         allocate (pcptemp(neq), rlftemp(neq), rvftemp(neq))
         allocate (icaptemp(neq))
      end if
      ndummy = 0

      stemp(1:neq)    = s(1:neq)
      ieostemp(1:neq) = ieos(1:neq)
      irlptemp(1:neq) = irlp(1:neq)
      icaptemp(1:neq) = icap(1:neq)
      pcptemp(1:neq)  = pcp(1:neq)
      rlftemp(1:neq)  = rlf(1:neq)
      rvftemp(1:neq)  = rvf(1:neq)

      if (num_sat .eq. 0) then
c gaz 071423 made sure there are enough s increments
         neq = int(1 / delta_sat) + 2
      else
         neq = num_sat
         sat_out(num_sat) = 1.d0
      end if
      s_last = .false.
      do i = 1, neq
         if (num_sat .eq. 0) then
            s(i) = (i - 1) * delta_sat
            if(s(i).ge.1.0d0) then
             if(abs(s(i)-s(i-1)).lt.s_comtol) then
              neq = i-1
              s_last = .true.
              go to 100
             else
              s(i) = 1.d0
              neq = i
              go to 100
              s_last = .true.
            endif
          endif
         else
            s(i) = sat_out(i)       
         end if
      enddo
100      continue
         do i = 1, neq
          ieos(i) = 2
c gaz 071423 added gdkm flag
          if (idpdp .ne. 0.or.gdkm_flag.ne.0) then
            s(i + neq) = s(i)
            ieos(i + neq) = 2
          end if
         end do
c get file name
      inquire(unit = ishisrlp, name = fname, opened = op) 
      if(op) then
      fname_new = trim(fname)
      close(unit = ishisrlp,status='delete')
      endif
      do j = 1, nrlp 
c gaz 071623 get file name 
      ishisrlpn = ishisrlp + j
      fname = fname_new
      write(model_string,'(a6,i4)')  '_model',1000+j 
      model_string(7:7)='_' 
      ii = len_trim(fname)
      write(fname(ii+1:ii+10),'(a)') model_string(1:10)
      select case (form_flag)
      case (0)
                fname = trim(fname) // '.tbl'
      case (1)
               fname = trim(fname) // '_tbl.dat'
      case (2)
               fname = trim(fname) // '_tbl.csv'
      end select 
      open (unit=ishisrlpn, file=fname, form='formatted')
      select case (form_flag)
      case (0)
               write(ishisrlpn, 6001) verno, jdate, jtime
      case (1)
               write(ishisrlpn, 6005) verno, jdate, jtime
      end select    
c form_flag = 1 tecplot; 2 csv or sur      
      title_string = "Relative permeability and " //
     &     "Capillary pressure"
      input_file_string ='Table generated from ' // nmfil(2)
      if (form_flag .eq. 1) then
c tecplot style
         form_string = 'variables = "Saturation" ' //
     &        '"Liquid" "Vapor" "Capillary pressure"'
         write(ishisrlpn, '(a)') trim(form_string)
         write(ishisrlpn, 230) 50., 95., trim(title_string)
         write(ishisrlpn, 230) 50., 90., trim(wdd)
      else if (form_flag .eq. 2) then
         form_string = 'Saturation, Liquid, Vapor, ' //
     &        'Capillary pressure'
         write(ishisrlpn, '(a)') trim(form_string)
      else
         form_string = '"Saturation" "Liquid" "Vapor" ' //
     &        '"Capillary pressure"'
         write(ishisrlpn, '(a)') trim(title_string)
         write(ishisrlpn, '(a)') trim(input_file_string)
         write(ishisrlpn, '(a)') trim(form_string)
      end if
         
         do i = 1, neq
            irlp(i) = j
            icap(i) = j
c gaz 071623 add gdkm_flag
            if (idpdp .ne. 0.or.gdkm_flag.ne.0) then
               irlp(i+neq) = j
               icap(i+neq) = j
            end if
         end do
         if (form_flag .eq. 1) then
            write (ishisrlpn, 240) j
         else if (form_flag .eq. 2) then
         else
         end if  
c calculate relperms and cap pressure for each saturation value                     
         if (rlpnew) then
            call rlp_cap(0)
c          ! its a vg model with fractures
            if (idpdp .ne. 0 .and. rlp_fparam(j, 1) .ge. 0.) 
     &         frac_model = .true.  
            if (frac_model) call rlp_cap(neq)
         else
            call rlperm(0,1)
            call cappr(1,0)
            if (abs(irlpt(j)) .eq. 4 .or. irlpt(j) .eq. 6 .or. 
     &           irlpt(j) .eq.7) then
c     This is a Van Genuchten fracture model
               frac_model = .true.
            end if
            if (idpdp .ne. 0 .and. frac_model) then
               call rlperm(neq,1)
               call cappr(1,neq)
            end if 
         end if
c output values         
! fracture   
		if(frac_model) write (ishisrlpn,*) 'Fracture nodes:'      
         do i = 1, neq
c gaz 071523 
            write (dum_tbl, '(4(g16.9, 1x))') s(i), rlf(i), rvf(i),
     &           pcp(i)
            ii = len_trim(dum_tbl)
            do jj = 1,ii
             if(dum_tbl(jj:jj).ne.' ') then
              go to 101
             endif
            enddo
101         write (ishisrlpn, '(t1,a)') dum_tbl(jj:ii)
         end do
! Matrix         
        if (frac_model) then
		  write (ishisrlpn,*) 'Matrix nodes:'        
         do i = neq + 1, 2*neq
               write (dum_tbl, '(4(g16.9, 1x))') s(i), rlf(i), rvf(i),
     &              pcp(i)
             do jj = 1,180
              if(dum_tbl(jj:jj).ne.' ') then
               go to 102
              endif
             enddo
102         write (ishisrlpn, '(t1,a)') dum_tbl(jj:80)
         end do
        end if
c gaz 071723 close each individual rlp table
      close (ishisrlpn)
      end do

      neq = neqtemp

      s(1:neq)    = stemp(1:neq)
      ieos(1:neq) = ieostemp(1:neq)
      irlp(1:neq) = irlptemp(1:neq)
      icap(1:neq) = icaptemp(1:neq)
      pcp(1:neq)  = pcptemp(1:neq)
      rlf(1:neq)  = rlftemp(1:neq)
      rvf(1:neq)  = rvftemp(1:neq)

      deallocate (stemp, ieostemp, irlptemp,  icaptemp,  
     &     pcptemp, rlftemp, rvftemp)
      close (ishisrlpn)
 230  format ("text X=", f4.1, " Y=", f4.1, " AN=center T=", '"',
     &     a, '"')
 240  format ('zone t = "model ', i4, '"')
 245  format ('zone t = "model ', i4, ' matrix"')
 6000 format(a30, 2x, a10, 2x, a8, /, a)
 6001 format(a30, 2x, a10, 2x, a8)
 6005 format('TITLE = "', a30, 2x, a10, 2x, a8,'"')
      end subroutine check_rlp
