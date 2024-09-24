      subroutine write_avs_node_hf(icall,neq,nheatflux,lu,ifdual)
!***********************************************************************
! Copyright 2014. Los Alamos National Security, LLC.  This material was 
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Output AVS heat flux node information for FEHM
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 16-July-2014  Z. Dash
CD2
C***********************************************************************

      use avsio
      use comai, only : altc, days, icnl, jdate, jtime, verno, wdd,
     &     nei_in, ns_in, neq_primary, idoff
      use combi, only : corz, izonef, nelm, elem_geo
      use comdi, only : nsurf, izone_surf, izone_surf_nodes
      use comflow, only : e_adv_nodal, e_cond_nodal
      implicit none
      
      integer maxvector
      parameter (maxvector = 3)
      character*3 dlstring
      character*5 cstring, dual_char
      character*6 share_string
      character*14 tailstring
      character*150 :: tec_string = ''
      character*45 title(3*maxvector), title2(2)
      character*50 fstring
      character*500 tstring
      integer icall, neq, nheatflux, lu, ifdual, length
      integer i, istep, nout, iz, idz, iendz, j, k, iocord_temp
      integer icord1, icord2, icord3, il, i1, i2, open_file
      integer, allocatable :: nelm2(:)
      real*8  write_array(9)
      character*5 char_type
      
      save tec_string

      if(nheatflux .gt. maxvector)then
         write(lu,*)'--------------------------------------------'
         write(lu,*)'ERROR: WRITE_AVS_NODE_HF'
         write(lu,*)'nvector   = ',nheatflux,' is greater than '
         write(lu,*)'maxvector = ',maxvector
         write(lu,*)'Subroutine only able to handle 2 vectors'
         write(lu,*)'--------------------------------------------'
         return
      endif
      
      iocord_temp = iocord
      if(ifdual .eq. 0)then
         dual_char = ''
         tailstring = '_hf_node'
      elseif(ifdual .eq. 1)then
         if (iodual .eq. 1) then
            dual_char ='Dual '
            tailstring = '_hf_dual_node'
         else
            dual_char = 'GDKM '
            tailstring = '_hf_gdkm_node'
            if (icnl .eq. 0) then
               iocord = 3
            else
               iocord = 2
            end if
         end if
      else
c     error
      endif

      if (iozone .ne. 0 ) then
         iendz = nsurf
      else 
         iendz = 1
         idz = iozone
      end if
      
      if (altc(1:3) .eq. 'tec' .or. altc(1:3) .eq. 'sur') then
         if (iocord .eq. 0 .and. iogrid .eq. 0) then
! Coordinates need to be output unless grid is set
            if (icnl .eq. 0) then
               iocord = 3
            else
               iocord = 2
            end if
         end if
      end if
      if (iocord .ne. 0) then
         select case (icnl)
         case (1, 4)
            icord1 = 1
            icord2 = 2
            icord3 = 1
         case (2, 5)
            icord1 = 1
            icord2 = 3
            icord3 = 2
         case(3, 6)
            icord1 = 2
            icord2 = 3
            icord3 = 1
         case default
            icord1 = 1
            icord2 = 3
            icord3 = 1
         end select
      end if

      nout = nheatflux
      if (iocord .ne. 0) nout = nout + 1

      title = ''
      title2 = ''
      if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
         j = 0
         if (iocord .ne. 0 ) then
            j = j + 1
            title(j) = trim(dual_char) // 'XYZ coordinates, (m)'
            icord1 = 1
            icord2 = 3
            icord3 = 1
            iocord = 3
         end if
         if (idoff .ne. -1 ) then
            j = j + 1
            title(j) =  trim(dual_char) // 
     &           'Advective Heat Flux (MW/m2), (MW/m2)'
         end if
         j = j + 1
         title(j) =  trim(dual_char) // 
     &        'Conductive Heat Flux (MW/m2), (MW/m2)'
         dlstring = '" "'
      else
         if (altc(1:4) .eq. 'avsx') then
            dlstring = ' : '
         else if (altc(1:3) .eq. 'tec') then
            dlstring = '" "'
         else if (altc(1:3) .eq. 'sur') then
            dlstring = ', '
         end if
         j = 0
         if (altc(1:3) .eq. 'tec' .or. altc(1:3) .eq. 'sur') then
            if (altc(1:3) .eq. 'tec' .and. iocord .eq. 0) then
               title2(1) = trim(dual_char) // 'Node'         
            else if (altc(1:3) .eq. 'tec') then
               title2(1) = dlstring // trim(dual_char) // 'Node'
            else
               title2(1) = trim(dual_char) // 'Node'
            end if    
            k = 1
            if (iozid .ne. 0) then
               title2(2) = dlstring // trim(dual_char) // 'Zone'
               k = 2
            end if
         end if
         if (iocord .ne. 0 ) then
            j = 0
            if (icnl .ne. 3 .and. icnl .ne. 6) then
               j = j + 1
               if (altc(1:3) .eq. 'tec') then
                  title(j) = 'X Coordinate (m)'
               else
                  title(j) = dlstring // 'X Coordinate (m)'
               end if
            end if
            if (icnl .ne. 2 .and. icnl .ne. 5) then
               j = j + 1
               if (j .eq. 1 .and. altc(1:3) .eq. 'tec') then
                  title(j) = 'Y Coordinate (m)'
               else
                  title(j) = dlstring // 'Y Coordinate (m)'
               end if
            end if
            if (icnl .ne. 1 .and. icnl .ne. 4) then
               j = j + 1
               title(j) = dlstring // 'Z Coordinate (m)'
            end if
            if (iozid .ne. 0) then
               write (share_string, '(i1, "-", i1, ", ", i1)') 1,j,j+2
            else
               write (share_string, '(i1, "-", i1)') 1, j
            end if
         end if
         if (idoff .ne. -1) then
            j = j + 1
            if (icnl .ne. 3 .and. icnl .ne. 6) then
               title(j) = dlstring // trim(dual_char) // 
     &              'Advective X Heat Flux (MW/m2)'
               j = j + 1
            end if
            if (icnl .ne. 2 .and. icnl .ne. 5) then
               title(j) = dlstring // trim(dual_char) //  
     &              'Advective Y Heat Flux (MW/m2)'
               j = j + 1
            end if
            if (icnl .ne. 1 .and. icnl .ne. 4) then
               title(j) = dlstring // trim(dual_char) //  
     &              'Advective Z Heat Flux (MW/m2)'
            end if
         end if
         j = j + 1
         if (icnl .ne. 3 .and. icnl .ne. 6) then
            title(j) = dlstring // trim(dual_char) //  
     &           'Conductive X Heat Flux (MW/m2)'
            j = j + 1
         end if
         if (icnl .ne. 2 .and. icnl .ne. 5) then
            title(j) = dlstring // trim(dual_char) //  
     &           'Conductive Y Heat Flux (MW/m2)'
            j = j + 1
         end if
         if (icnl .ne. 1 .and. icnl .ne. 4) then
            title(j) = dlstring // trim(dual_char) //  
     &           'Conductive Z Heat Flux (MW/m2)'
         end if
      end if
      nout = j

      if (altc(1:3) .ne. 'sur') then
         call namefile2(icall,lu,ioformat,tailstring,0)
! file will be opened in zone loop for surfer
      end if

      if(altc(1:4).eq.'avsx') then
         write(lu,200) days, (trim(title(i)), i = 1, nout)           
      else if(altc(1:3).eq.'avs') then
         write(lu,90) nout, (3, i = 1, nout)
         write(lu, 100) (trim(title(i)), i = 1, nout)
      else if(altc(1:3).eq.'tec') then
         if (icall .eq. 1 .or. (iogrid .eq. 1 .and. iocord .eq. 0)) then
            write(lu,301) verno, jdate, jtime, trim(wdd)
            if (iogrid .eq. 1 .and. iocord .eq. 0) write(lu,304)
c            write(lu,300) (trim(title(i)), i = 1, iocord), 
c     &           (trim(title2(i)), i = 1, k), 
c     &           (trim(title(i)), i = iocord + 1, nout), '"'
            tstring = 'VARIABLES = "'
            do i = 1, iocord
               tstring = trim(tstring) // trim(title(i))
            end do
            do i = 1, k
               tstring = trim(tstring) // trim(title2(i))
            end do
            do i = iocord + 1, nout
               tstring = trim(tstring) // trim(title(i))
            end do
            tstring = trim(tstring) // '"'
            write(lu, '(a)') trim(tstring)
         end if
         if (iozone .ne. 0) write(lu,302) trim(timec_string)
      else if(altc(1:3).eq.'sur') then
c         write(tstring,400) (trim(title2(i)), i = 1, k),
c     &        (trim(title(i)), j = 1, nout)
         tstring = ''
         do i = 1, k
            tstring = trim(tstring) // trim(title2(i))
         end do
         do i = 1, nout
            tstring = trim(tstring) // trim(title(i))
         end do
      end if

      do iz = 1, iendz
! Zone loop
         if (iozone .ne. 0) then
            idz = izone_surf(iz)
         end if
         if (altc(1:3) .eq. 'tec') then
            if (iozone .ne. 0) then
               write (lu, 95) idz, trim(tec_string)
            else
               if (iogrid .eq. 1 .and. iocord .eq. 0) then
                  write (lu, 96) trim(timec_string),
     &              trim(gridstring), trim(times_string)
               else if (iogeo .eq. 1 .or. (ifdual .eq. 1 .and. 
     &              iogdkm .eq. 1)) then
                   write (lu, 96) trim(timec_string), trim(gridstring), 
     &                 trim(tec_string)
               else
                  write (lu, 96) trim(timec_string), trim(tec_string)
               end if
            end if
            if (icall .eq. 1) then
               if (iogeo .eq. 1) then
                  write (tec_string, 305) trim(share_string), iz
               else
                  write (tec_string, 303) trim(share_string), iz
               end if
            end if
         else if (altc(1:3) .eq. 'sur') then
            call namefile2(icall,lu,ioformat,tailstring,idz)
            write(lu, '(a)') trim(tstring)
         end if
         do i = 1,neq
! Node loop
            if (iozone .ne. 0) then
               if (izone_surf_nodes(i).ne.idz) goto 199
            end if
            if (iocord .ne. 0) then
               k = 1
               do j = icord1, icord2, icord3
                  write_array(k) = corz(i, j)
                  k = k + 1
               end do
               j = iocord + 1
            else
               j = 1
            end if
            if (idoff .ne. -1) then
               write_array(j) = e_adv_nodal(i,1)
               j = j + 1
               write_array(j) = e_adv_nodal(i,2)
               j = j + 1
               if (icnl .eq. 0 .or. altc .eq. 'avs') then
                  write_array(j) = e_adv_nodal(i,3)
                  j = j + 1
               end if
            end if
            write_array(j) = e_cond_nodal(i,1)
            j = j + 1
            write_array(j) = e_cond_nodal(i,2)
            j = j + 1
            if (icnl .eq. 0 .or. altc .eq. 'avs') then
               write_array(j) = e_cond_nodal(i,3)
               j = j + 1
            end if
            j = j - 1
            if (altc(1:3) .eq. 'tec') then
               if (icall .eq. 1 .and. iozid .ne. 0) then
                  write (fstring, 130) iocord, j - iocord
                  write(lu, fstring) (write_array(k), k = 1, iocord),
     &                 i, izonef(i), (write_array(k), k = iocord + 1, j)
               else if (icall .eq. 1 .and. iozid .eq. 0 .and.
     &                 iogrid .eq. 0) then
                  write (fstring, 120) iocord, j - iocord
                  write(lu, fstring) (write_array(k), k = 1, iocord),
     &                 i, (write_array(k), k = iocord + 1, j)
               else if (icall .eq. 1 .and. ifdual .eq. 1 .and. 
     &                 iogdkm .eq. 1) then
                  write (fstring, 120) iocord, j - iocord
                  write(lu, fstring) (write_array(k), k = 1, iocord),
     &                 i, (write_array(k), k = iocord + 1, j)
               else
                  write (fstring, 110) j - iocord
                  write(lu, fstring) i, (write_array(k), 
     &                 k = iocord + 1, j)
               end if
            else if (altc(1:3).eq.'sur') then
               if (iozid .ne. 0) then
                  write (fstring, 220) dlstring, j, dlstring
                  write(lu, fstring) i, izonef(i), 
     &                 (write_array(k), k = 1, j)
               else
                  write (fstring, 210) j, dlstring
                  write(lu, fstring) i, (write_array(k), k = 1, j)
               end if
            else if (altc(1:4).eq.'avsx') then
               write (fstring, 210) j, dlstring
               write(lu, fstring) i, (write_array(k), k = 1, j)
            else              
               write (fstring, 110) j
               write(lu, fstring) i, (write_array(k), k = 1, j)
            end if
 199     enddo 
      end do

      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. iogeo .eq. 1)
     &     then
c     Read the element connectivity and write to tec file
         if(ifdual .eq. 1 .and. iogdkm .eq. 1) then
c     Do nothing, no connectivity defined
         else 

c gaz 052222 adding tec and other if blocks            
c gaz 060822 using elem_geo              
            do i = 1, nei_in
               write(lu, '(8(i8))') (elem_geo((i-1)*ns_in+j), j=1,ns_in)
            end do          
          endif
      end if

      iocord = iocord_temp

 90   format(i1,2x,5(i5,2x))
 95   format('ZONE T = "',i4.4, '"', a)
c 96   format('ZONE T = "Simulation time ',1p,g16.9,' days"', a)
 96   format('ZONE T = ', a, a, a, a)   
 500  format(i10.10, 9(1x, g16.9))
 666  format(i10.10, 9(' : ', g16.9))
 667  format(i10.10, 9(', ', g16.9))
 100  format(a)
 110  format("(i10.10, ", i1, "(1x, g16.9))")
 120  format("(", i1, "(g16.9, 1x), i10.10, ", i1, "(1x, g16.9))" )
 130  format("(", i1, "(g16.9, 1x), i10.10, 1x, i4,", i1, 
     &     "(1x, g16.9))") 
 200  format('nodes at ', g16.9,' days', 9(a))
 210  format("(i10.10, ", i1, "('", a, "',g16.9))") 
 220  format("(i10.10, '", a, "', i4, ", i1, "('", a, "',g16.9))") 
 300  format('VARIABLES = "', 12(a))
 301  format ('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
c 302  format('TEXT T = "Simulation time ',1p,g16.9,' days"')
 302  format('TEXT T = ', a)
 303  format(', VARSHARELIST = ([', a,'] = ', i4, ')')
 304  format('FILETYPE = "SOLUTION"')
 305  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')
 400  format(11(a))
 
      return
      end
