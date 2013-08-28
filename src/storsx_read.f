      subroutine storsx_read
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
      integer i, iwtotl, iwtotl_strs, j, ncont,  neq_old,  ncoef
      integer ncont_new, neq_save, ncont_primary, neqp1
      integer iwt, jjt, ict, iwtotl_strs_temp, i1, i2, kb, intg_old
c parser variables
      character*80 input_msg
      integer msg(10)
      integer nwds,sehtemp
      real*8 xmsg(10)
      integer imsg(10)
      character*32 cmsg(10)
c local
      logical opend
      logical exists
      integer ilen, rlen, flen
      integer ityp
      integer :: max_con = 0
      integer, allocatable :: nelm_temp(:)
      integer, allocatable :: istrw_temp(:)
      real*8 tot
c      character*100 filename, tail
      character*120 filename, tail
      character*72 cline
      character*32 sxformat
      character*3 stat_var

c**** retrieve coefficients ****

      if(lda.eq.1 .or. lda.eq.5 .or. lda.eq.6) then
         sxformat(1:5) = 'ascii'
         ityp = 2
      else if(lda.eq.2 .or. lda.eq.7 .or. lda.eq.8) then
         sxformat(1:5) = 'unfor'
         ityp = 3
      else if(lda.eq.3) then
         sxformat(1:5) = 'binar'
         ityp = 1
      else if(lda.eq.4) then
         sxformat(1:5) = 'unkno'
         ityp = 0
      endif

      inquire(unit=isstor,name=filename)
      inquire(unit=isstor,opened=opend)
      if (opend) close(isstor)

 1500 format (1x, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
 1501 format (1x, 'stor file has neq less than data file: quit')
 1502 format (1x, 'error in parsing beginning of stor file')
 1503 format (1x, 'stor file has unrecognized format: quit')
 1504 format (1x, 'Coefficients read from file ', a)

      intg_old = intg
      iwtotl_strs_temp = 0
c-----------READ ASCII  format stor file ------------------------------
      if (sxformat(1:5).eq.'ascii') then

         if(iptty.ne.0) then
            write(iptty,*)'WARNING: ASSUMING NEW ASCII STOR FORMAT'
         end if

         open(isstor,file=filename ,status='old',form='formatted')

         read (isstor,'(a72)', END=5000, err=5000) cline
         read (isstor,'(a72)', END=5000, err=5000) cline
c         read (isstor,'(7i10)', END=4100, ERR=4100) iwtotl,neq_old,
c     &        ncont,sehtemp,max_con, iwtotl_strs_temp, intg
         read (isstor, '(a80)') input_msg
         call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
         iwtotl=imsg(1)
         neq_old=imsg(2)
         ncont=imsg(3)
         if (nwds.ge.4) then
            sehtemp=imsg(4)
         end if
         if (nwds.ge.5) then
            max_con=imsg(5)
         end if
         if (nwds .eq. 7) then
            iwtotl_strs_temp = imsg(6)
            intg = imsg(7)
            if(intg.ne.intg_old) then
               if(iout.ne.0)write(iout,*)
     &              'NOTE >>> read coef. generated with intg = ', intg
               if(iptty.ne.0)write(iptty,*)
     &              'NOTE >>>read coef. generated with intg = ', intg  
            endif 
            if(iwtotl_strs_temp.eq.0.and.istrs.ne.0) then
               lda = -abs(lda)
               return
            endif
         end if
c         go to 4101
c 4100    iwtotl_strs_temp = 0
c  4101    continue
c     check for no stress solution (set number of coefficients to 3)
         if(iwtotl_strs_temp.ne.0.and.istrs.eq.0) then
            sehtemp = 3
         endif

         if(neq_old.lt.neq_primary) then
c     enable mdnodes
            imdnode=-neq_old
         else if(neq_old.gt.neq_primary) then
c     error occurs
            write(ierr, 1500)
            write(ierr, 1501)
            write(ierr, 1500)
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1501)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1501)
               write(iptty, 1500)
            endif
            stop
         endif
         
c     add 1 more conection to be used for BCtoBC and mdnode
         nr=iwtotl+1

c     call mdnode to arrange additional space 
c     for connectivity when we have mdnodes

         if(imdnode.ne.0) then
            ncont_new=ncont
            call md_nodes(4,0,ncont_new)
         else
            ncont_new=ncont
         endif
c     
         deallocate(nelm)
         allocate(nelm(ncont_new),istrw(ncont_new))
c     gaz 11-09-2001
c     allocation of istrw_itfc and istrw_cold done in startup
         nelmd=ncont_new
         neqp1 = neq_old + 1
         read (isstor,*)  (sx1(i), i = 1, neq_old)
         read (isstor,*)  (nelm(i), i = 1, ncont)
         read (isstor,*)  (istrw(i), i = 1, ncont)
         read (isstor,*)  (nelmdg(i), i = 1, neq_old)

         if (istrs .eq. 0)  then
            if (sehtemp.eq.1) then
               isoy=1
               isoz=1
               allocate(sx(nr,1))
               sx=0
               call read_sx(nr,iwtotl,1,isstor,0,ityp)
            else if (sehtemp.eq.3.or.sehtemp.eq.4) then
               allocate(sx(nr,3))
               call read_sx(nr,iwtotl,3,isstor,0,ityp)
            else if (sehtemp.eq.6) then
               allocate(sx(nr,6))
               call read_sx(nr,iwtotl,6,isstor,0,ityp)
            endif
         else
            iwtotl_strs = ncont-neqp1
            allocate (istrws(iwtotl_strs))
            allocate(sx(nr,3))
            call read_sx(nr,iwtotl,3,isstor,0,ityp)
            read (isstor,*)  (istrws(i), i = 1, iwtotl_strs)
            if(icnl.eq.0) then
               ncoef=15
            else
               ncoef=8
            endif
            allocate(sxs(iwtotl_strs,ncoef))
            if(icnl.eq.0) then 
               allocate(dnidnj(iwtotl_strs, 9))
            else
               allocate(dnidnj(iwtotl_strs, 4))
            endif

            call read_sx(iwtotl_strs,iwtotl_strs,ncoef,
     &           isstor,1,ityp)
         end if

c-----------READ UNFORMATTED  format stor file ------------------------
      elseif (sxformat(1:5).eq.'unfor') then
         if(iptty.ne.0) then
            write(iptty,*)'WARNING: ASSUMING NEW UNFORMATTED',
     &           ' STOR FORMAT'
         end if

         open(isstor,file=filename ,status='old',
     &        form='unformatted')

         read (isstor, END=5000, err=5000) cline
         read (isstor, END=5000, err=5000) cline
         read (isstor, END=4002, ERR=4002) iwtotl, neq_old, ncont, 
     &        sehtemp, max_con, iwtotl_strs_temp, intg
         if(intg.ne.intg_old) then
            if(iout.ne.0)write(iout,*)
     &           'NOTE >>> read coef. generated with intg = ',intg
            if(iptty.ne.0)write(iptty,*)
     &           'NOTE >>>read coef. generated with intg = ',intg  
         endif      
         go to 4003
 4002    iwtotl_strs_temp = 0
         if(iwtotl_strs_temp.eq.0.and.istrs.ne.0) then
            lda = -abs(lda) 
            return
         endif
 4003    continue
c     check for no stress solution (set number of coefficients to 3)
         if(iwtotl_strs_temp.ne.0.and.istrs.eq.0) then
            sehtemp = 3
         endif

c     adding possibility that the *.stor file may have less nodes than
c     what was defined in the data file. We assume that means additional
c     nodes were defined for multiply connected situations. Further
c     all the new nodes will be after neq
c     
         if(neq_old.lt.neq_primary) then
c     enable mdnodes
            imdnode=-neq_old
         else if(neq_old.gt.neq_primary) then
c     error occurs
            write(ierr, 1500)
            write(ierr, 1501)
            write(ierr, 1500)
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1501)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1501)
               write(iptty, 1500)
            endif
            stop
         endif
c     
c     always add one more conection to be used for BC to BC
c     and mdnode connections
c     
         nr=iwtotl+1
c     
c     call mdnode to arrange additional space for connectivity
c     when we have mdnodes
c     
         if(imdnode.ne.0) then
            ncont_new=ncont
            call md_nodes(4,0,ncont_new)
         else
            ncont_new=ncont
         endif
c     
         neqp1 = neq_old + 1
         deallocate(nelm)
         allocate(nelm(ncont_new),istrw(ncont_new))
c     gaz 11-09-2001
c     allocation of istrw_itfc and istrw_cold done in startup
         nelmd=ncont_new

         read (isstor)  (sx1(i), i = 1, neq_old)
         read (isstor)  (nelm(i), i = 1, ncont)
         read (isstor)  (istrw(i), i = 1, ncont)
         read (isstor)  (nelmdg(i), i = 1, neq_old)
         if (istrs .eq. 0)  then
            if (sehtemp.eq.1) then
               isoy=1
               isoz=1
               allocate(sx(nr,1))
               sx=0
               call read_sx(nr,iwtotl,1,isstor,0,ityp)
            else if (sehtemp.eq.3.or.sehtemp.eq.4) then
               allocate(sx(nr,3))
               call read_sx(nr,iwtotl,3,isstor,0,ityp)
            else if (sehtemp.eq.6) then
               allocate(sx(nr,6))
               call read_sx(nr,iwtotl,6,isstor,0,ityp)
            endif
         else
            iwtotl_strs = ncont-neqp1
            allocate (istrws(iwtotl_strs))
            allocate(sx(nr,3))
            call read_sx(nr,iwtotl,3,isstor,0,ityp)
            read (isstor)  (istrws(i), i = 1, iwtotl_strs)
            if(icnl.eq.0) then
               ncoef=15
            else
               ncoef=8
            endif
            allocate(sxs(iwtotl_strs,ncoef))
            if(icnl.eq.0) then 
               allocate(dnidnj(iwtotl_strs, 9))
            else
               allocate(dnidnj(iwtotl_strs, 4))
            endif

            call read_sx(iwtotl_strs,iwtotl_strs,ncoef,
     &           isstor,1,ityp)

         end if

c-----------READ OLD ASCII  format stor file --------------------------
c     if format not recognized, try old ascii format
      elseif (sxformat(1:7).eq.'unknown') then

         open(isstor,file=filename ,status='old',form='formatted')

         read (isstor, '(a80)', END=5000, err=5000)     input_msg
         read (isstor, '(a80)', END=5000, err=5000)     input_msg
         read (isstor, '(a80)') input_msg
         call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
         iwtotl=imsg(1)
         neq_old=imsg(2)
         ncont=imsg(3)
         if (nwds.eq.4) then
            sehtemp=imsg(4)
         else if (nwds.eq.5) then
            sehtemp=imsg(4)
         else if (nwds.eq.3) then
            sehtemp=6
         else
            write(ierr, 1500)
            write(ierr, 1502)
            write(ierr, 1500)
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1502)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1502)
               write(iptty, 1500)
            endif
            stop
         endif
c     
c     adding possibility that the *.stor file may have less nodes than
c     what was defined in the data file. We assume that means additional
c     nodes were defined for multiply connected situations. Further
c     all the new nodes will be after neq
c     
         if(neq_old.lt.neq_primary) then
c     enable mdnodes
            imdnode=-neq_old
         else if(neq_old.gt.neq_primary) then
c     error occurs
            write(ierr, 1500)
            write(ierr, 1501)
            write(ierr, 1500)
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1501)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1501)
               write(iptty, 1500)
            endif
            stop
         endif
c     
c     always add one more conection to be used for BC to BC 
c     and mdnode connections
c     
         nr=iwtotl+1
c     
c     call mdnode to arrange additional space for connectivity
c     when we have mdnodes
c     
         if(imdnode.ne.0) then
            ncont_new=ncont
            call md_nodes(4,0,ncont_new)
         else
            ncont_new=ncont
         endif
c     
         neqp1 = neq_old + 1
         deallocate(nelm)
         allocate(nelm(ncont_new),istrw(ncont_new))
c     gaz 11-09-2001
c     allocation of istrw_itfc and istrw_cold done in startup
         nelmd=ncont_new

         read (isstor, *,err=990)  (sx1(i), i = 1, neq_old)
         read (isstor, *,err=990)  (nelm(i), i = 1, ncont)
         read (isstor, *,err=990)  (istrw(i), i = 1, ncont)
         read (isstor, *,err=990)  (nelmdg(i), i = 1, neq_old)
         if (istrs .eq. 0)  then
            if (sehtemp.eq.1) then
               isoy=1
               isoz=1
               allocate(sx(nr,1))
               sx=0
               call read_sx(nr,iwtotl,1,isstor,0,ityp)
            else if (sehtemp.eq.3) then
               allocate(sx(nr,3))
               call read_sx(nr,iwtotl,3,isstor,0,ityp)
            else if (sehtemp.eq.6) then
               allocate(sx(nr,6))
               call read_sx(nr,iwtotl,6,isstor,0,ityp)
            endif
         else 
            iwtotl_strs = ncont-neqp1
            allocate (istrws(iwtotl_strs))
            read (isstor)  (istrws(i), i = 1, iwtotl_strs)
            allocate(sx(nr,3))
            call read_sx(nr,iwtotl,3,isstor,0,ityp)
            if(icnl.eq.0) then
               ncoef=15
            else
               ncoef=8
            endif
            allocate(sxs(iwtotl_strs,ncoef))
            if(icnl.eq.0) then 
               allocate(dnidnj(iwtotl_strs, 9))
            else
               allocate(dnidnj(iwtotl_strs, 4))
            endif

            call read_sx(iwtotl_strs,iwtotl_strs,ncoef,
     &           isstor,1,ityp)
         end if
         
c     made it without encountering read errors
         sxformat='oldascii'
 990     if (sxformat(1:7).eq.'unknown') then
            write(ierr, 1500)
            write(ierr, 1503)
            write(ierr, 1500)
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1503)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1503)
               write(iptty, 1500)
            endif
            stop
         endif

      endif
c-----END switch on reading stor formats--------------------------------
      if(iptty.ne.0) then
         write(iptty,*)'Coefficient file read with format ',
     &        sxformat(1:15)
      end if

c     simplify ncon if stress solution is disabled
      if(iwtotl_strs_temp.ne.0.and.istrs.eq.0) then
         allocate(nelm_temp(ncont),istrw_temp(ncont-neqp1))
         ict = neqp1
         nelm_temp(1) = ict
         istrw_temp(1) = 1
         do i = 1,neq_old
            i1 = nelm(i)+1
            i2 = nelm(i+1)
            do j = i1,i2
               iwt = istrw(j-neqp1)
               kb = nelm(j)
               if(kb.ne.i) then
                  tot = abs(sx(iwt,1)) + abs(sx(iwt,2)) + abs(sx(iwt,3))
               else
c     skip kb = i
                  tot = 1.0d0
               endif
               if(tot.gt.0.0d0) then
                  ict = ict +1
                  nelm_temp(ict) = kb
                  istrw_temp(ict-neqp1) = iwt
               endif
            enddo
            nelm_temp(i+1) = ict
         enddo 
         deallocate(nelm,istrw)
         allocate(nelm(ict),istrw(ict-neqp1))
         nelm(1:ict) = nelm_temp(1:ict)
         istrw(1:ict-neqp1) = istrw_temp(1:ict-neqp1)
         deallocate(nelm_temp,istrw_temp)
      endif

c     
c     insure isotropy for non stress solutions
c     
      if(isoy.eq.2) then
         if(gdpm_flag.ne.0) then
            neq_save = neq
            neq = neq_primary
            call sx_combine(0)
            neq = neq_save
         else
            call sx_combine(0)
         end if
c     done now in startup.f
c     call sx_combine(1)
         call sx_combine(2)
      else if(isoy.eq.1) then
         call sx_combine(-3)
c     done now in startup.f
c     call sx_combine(1)
         call sx_combine(3)
      endif

c**** printout positions required for geometric coefficients ****
c     
      if (iout .ne. 0) then
         write(iout, 1500)
         write(iout, 1504) trim(filename)
         write(iout, 1500)
      end if
      if(iptty.gt.0) then
         write(iptty, 1500)
         write(iptty, 1504) trim(filename)
         write(iptty, 1500)
      endif
c     
      if (iout .ne. 0) write(iout, 2000) nr, nr
      if (iptty .gt. 0) write(iptty, 2000) nr,nr       
 2000 format(/,1x,'storage for geometric coefficients ', i10,
     *     ' in common(nr) ', i10)

      if (iwtotl .gt. nr)  then

         if (iout .ne. 0) write(iout, 3000)
         write(ierr, 3000)
         if (iptty .ne. 0) write(iptty, 3000)
 3000    format(/, 1x, 'program terminated because of ',
     *        'insufficient storage')

         stop

      end if
c     
c     
c     calculate additional storage if necessary for multiply defined nodes
c     
      if(imdnode.ne.0) then
         call md_nodes(5,0,ncont_new)
      endif

      close  (isstor)

      return
         
! Data could not be read
 5000 continue
 5001 format (/, 1x, 'program terminated because ',
     *     'data could not be read from storage file')
      if (iout .ne. 0) write(iout, 5001)
      write(ierr, 5001)
      if (iptty .ne. 0) write(iptty, 5001)
      stop  
 
      end subroutine storsx_read
