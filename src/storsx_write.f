      subroutine storsx_write
!***********************************************************************
! Copyright 2009 Los Alamos National Security, LLC  All rights reserved
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
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University.
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Store element coefficients.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 13-FEB-09,  Programmer: Z. Dash
!D2 
!***********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.2 Finite-Element Coefficient Generation
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
!D4    FEHM Application Version 2.20
!D4
!D2 storsx routine was split to simplify code maintenance
!D4
!***********************************************************************

      use comai
      use combi
      use comdti
      use comsi
      use comxi
      implicit none
c
      integer i, iwtotl, iwtotl_strs, j, ncont, neq_old, ncoef 
      integer ncont_new, neq_save, ncont_primary, neqp1
c     local
      logical opend
      logical exists
      integer ilen, rlen, flen
      integer ityp
      integer :: max_con = 0
      character*1000 filename, tail
      character*72 cline
      character*32 sxformat
      character*3 stat_var

c     Storing coefficients
      stat_var='old'
      inquire(unit=isstor,name=filename,form=sxformat)
c     If the file already contains data we don't want to overwrite it 
      filename = trim(filename)
      if(sxformat(1:9).eq.'FORMATTED') then
         read(isstor,'(a72)',end=900) cline
      else if(sxformat(1:11).eq.'UNFORMATTED') then
         read(isstor,end=900) cline
      endif
c     This file already contains data so a new one should be created
      flen = len_trim (filename)
      do i = flen, 1, -1
         if (filename(i:i) .eq. '.') then
            rlen = i-1
            ilen = flen-rlen
            tail(1:ilen) = filename(i:flen)
            exit
         end if
         if (i .eq. 1) rlen = flen
      end do
      if ((flen + 3) .gt. 100) rlen = rlen-3
      filename(rlen+1:rlen+3) = 'NEW'
      filename(rlen+4:rlen+4+ilen) = tail(1:ilen)

c     Do not overwrite this one until we are sure
      stat_var='new'
      inquire(file=filename,exist=exists)
      if (exists) then
         inquire(file='fehmn_temp.stor',exist=exists)
         if (exists) then
            write (ierr, 902) trim(filename)
            if (iout .ne. 0) write(iout,902) trim(filename)
            if(iptty.ne.0) write(iptty,902) trim(filename)
            stop
         else
            filename = ''
            filename(1:15)='fehmn_temp.stor'
         end if
      end if
      if (iout .ne. 0) write(iout,901) trim(filename)
      if(iptty.ne.0) write(iptty,901) trim(filename)
 901  format(1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     &     /,1x,'>>> Changing name of stor file (old file exists)',
     &     /,1x,'>>> New file name is: ', a)
      inquire(file='filename',exist=exists)
 902  format(1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!', /, 1x, '>>> Files: ', 
     &     a, /, ' and fehmn_temp.stor exist : stopping',
     &     /,1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
 900  continue
      close(isstor)

      if (lda .eq.-1)  then
C**** store coefficients (formatted) ****
         open(isstor,file=filename,status=stat_var,
     &        form='formatted')

 1000    format(a30, 3x, a8, 3x, a8)

c     Eliminated gdpm switch, just save primary node coefficients
c     for non-gdpm case neq = neq_primary
         iwtotl = iw
         neqp1 = neq_primary + 1
         ncont_primary  = nelm(neqp1)

         write(isstor, 1000)  verno, jdate, jtime
         write(isstor, '(a80)') wdd
         if(istrs.ne.0) then
            iwtotl_strs = ncont_primary - neqp1
            if(icnl.eq.0) then
               ncoef=15
            else
               ncoef=8
            endif
         else
            iwtotl_strs  = 0
            if(isoy.eq.1) then
               ncoef=1
            else if(isoy.eq.2) then
               ncoef=3
            endif
         endif

         write(isstor, '(7i10)'  ) iwtotl, neq_primary, 
     &        ncont_primary, ncoef, max_con, iwtotl_strs, intg
         write(isstor, '(5(1pe20.12))') 
     &        (sx1(i), i = 1, neq_primary)
         write(isstor, '(5i10)'  ) 
     &        (nelm(i), i = 1, ncont_primary)
         write(isstor, '(5i10)'  ) 
     &        (istrw(i), i = 1, ncont_primary)
         write(isstor, '(5i10)'  ) 
     &        (nelmdg(i), i = 1, neq_primary)

         if (istrs .eq. 0)  then

            if(isoy.eq.2) then
               write(isstor,'(5(1pe20.12))') 
     &              (sx(j,1)+sx(j,2)+sx(j,3), j = 1, iwtotl)
               write(isstor,'(5(1pe20.12))') (0.0d00 
     &              , j = 1, iwtotl)
               write(isstor,'(5(1pe20.12))') (0.0d00 
     &              , j = 1, iwtotl)

            else if(isoy.eq.1) then
               if(icnl.eq.0) then
                  write(isstor,'(5(1pe20.12))')
     &                 (3.0d00*sx(j,1), j = 1, iwtotl)
               else               
                  write(isstor,'(5(1pe20.12))')
     &                 (2.0d00*sx(j,1), j = 1, iwtotl)
               endif
               
            endif
         else

            do i = 1, 3
               write(isstor,'(5g20.12)') (sx(j,i)
     &              , j = 1, iwtotl)
            end do
c     there are more stress coeficients (ncont_primary) because symmetry is not used  
            write(isstor, '(5i10)'  ) 
     &           (istrws(i), i = 1, iwtotl_strs)
            do i = 1, ncoef
               write(isstor,'(5(1pe20.12))') (sxs(j,i),
     &              j = 1, iwtotl_strs)
c     &                 j = 1, ncont_primary)
            end do
         end if


      else if(lda.eq.-2) then
C**** store coefficients (unformatted) ****
         open(isstor,file=filename,status=stat_var,
     &        form='unformatted')

c     Eliminated gdpm switch, just save primary node coefficients
c     for non-gdpm case neq = neq_primary
         iwtotl = iw
         neqp1 = neq_primary + 1
         ncont_primary  = nelm(neqp1)

         cline = verno // ' ' // jdate // ' ' // jtime
         write(isstor) cline
         write(isstor) wdd(1:72)
         if(istrs.ne.0) then
            iwtotl_strs = ncont_primary - neqp1
            if(icnl.eq.0) then
               ncoef=12
            else
               ncoef=8
            endif
         else
            iwtotl_strs = 0
            if(isoy.eq.1) then
               ncoef=1
            else if(isoy.eq.2) then
               ncoef=3
            endif
         endif
         write(isstor) iwtotl, neq_primary, ncont_primary, ncoef, 
     &        max_con, iwtotl_strs, intg
         write(isstor) 
     &        (sx1(i), i = 1, neq_primary)
         write(isstor) 
     &        (nelm(i), i = 1, ncont_primary)
         write(isstor) 
     &        (istrw(i), i = 1, ncont_primary)
         write(isstor) 
     &        (nelmdg(i), i = 1, neq_primary)

         if (istrs .eq. 0)  then

            if(isoy.eq.2) then
               write(isstor) (sx(j,1)+sx(j,2)+sx(j,3)
     &              , j = 1, iwtotl)
               write(isstor) (0.0d00 
     &              , j = 1, iwtotl)
               write(isstor) (0.0d00 
     &              , j = 1, iwtotl)

            else if(isoy.eq.1) then
               if(icnl.eq.0) then
                  write(isstor)
     &                 (3.0d00*sx(j,1), j = 1, iwtotl)
               else               
                  write(isstor)
     &                 (2.0d00*sx(j,1), j = 1, iwtotl)
               endif

            endif
         else

            do i = 1, 3
               write(isstor) (sx(j,i)
     &              , j = 1, iwtotl)
            end do
c     there are more stress coeficients because symmetry is not used
            write(isstor) 
     &           (istrws(i), i = 1, iwtotl_strs)                
            do i = 1, ncoef
               write(isstor) (sxs(j,i),
     &              j = 1, iwtotl_strs)
c     &                 j = 1, ncont_primary)
            end do
         end if

      endif

      if (iout .ne. 0) then
         write(iout, 2000)
         write(iout, 2010) trim (filename)
         write(iout, 2000)
      end if
      if(iptty.gt.0) then
         write(iptty, 2000)
         write(iptty, 2010) trim (filename)
         write(iptty, 2000)
      endif

 2000 format (' !!!!!!!!!!!!!!!!!!!!!!!!!! ')
 2010 format (' Coefficients written to file: ', a)

      close (isstor)

      end subroutine storsx_write
