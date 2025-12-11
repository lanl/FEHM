      subroutine wrtptrk
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine outputs particle concentration data.
CD1
C***********************************************************************
CD2 
CD2  REVISION HISTORY 
CD2
CD2 Initial implementation:  Summer '94, Programmer: Stephen Henderson
CD2
CD2 $Log:   /pvcs.config/fehm90/src/wrtptrk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2 Made changes for ingrowth.  July 15, 1997. cli
CD2
CD2 Made necessary changes for multiple sepecies.   Apr. 18, 1997, CLI
CD2
CD2    Rev 1.4   Wed May 08 14:18:02 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.3   Thu Jan 11 12:43:04 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   03/15/95 17:05:48   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.1   02/02/95 15:22:54   llt
CD2 added pvcs log info
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.5 Cell-based particle-tracking module
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use comdti
      use comai
      use combi
      use comdi
      use comgi
      use comxi
      use compart
      use comsk

      implicit none

      integer i, i2, ith, mdd, num, mdd1, iroot, is, ie, n1
c bhl_5/15/08
      integer imbl1,ipart
c bhl_5/15/08
      integer open_file, inode, num_current, isptrk, isptrk1
      integer varlen, itmp
      real*8 num_enter, num_leave, num_decay, num_filtered, ret_wght
      real*8, allocatable :: num_in(:,:), num_left(:)
c bhl_5/15/08
      real*8 mass_enter,mass_leave,mass_decay,mass_filtered,mass_current
      real*8 mass_dout, mass_dfilt
c bhl_5/15/08
      character*120 ptrk_name, ptrk_root
      character*21 tstring, fstring
      character*66, allocatable :: string(:)
c bhl_5/15/08
      character*21 tstring2, fstring2
      character*112, allocatable :: string2(:)
c bhl_5/15/08
      logical end_flag, null1, opnd
c bhl_5/15/08
      logical imbl
c bhl_5/15/08
      save ptrk_root, iroot, isptrk1, num_left, varlen, fstring

      dtotc=dtotc/8.64e4

      end_flag  = .false.
	if (ripfehm .ne. 0) then
		itmp = iout
	else
		itmp = iptty
	end if

      if (.not. allocated(num_left)) then
         allocate (num_left(nspeci))
         num_left = 0
         if (abs(prnt_rst) .ge. 10) then
            if (null1(root_name)) then
               if (nmfil(10) .ne. nmfily(3) .and. nmfil(10) .ne. ' ') 
     &              then
                  call file_prefix(nmfil(10), iroot)
                  if (iroot .gt. 100) iroot = 100
                  ptrk_root(1:iroot) = nmfil(10)(1:iroot)
               else
                  if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') 
     &                 then
                     call file_prefix(nmfil(5), iroot)
                     if (iroot .gt. 100) iroot = 100
                     ptrk_root(1:iroot) = nmfil(5)(1:iroot)
                  else
                     if (nmfil(2)(1:1) .eq. ' ' ) then
                        write (ierr, *) 'FILE ERROR: nmfil2 file: ', 
     &                       nmfil(2),
     &                       ' unable to determine contour file prefix'
                        stop
                     else
                        call file_prefix(nmfil(2), iroot)
                        if (iroot .gt. 100) iroot = 100
                        ptrk_root(1:iroot) = nmfil(2)(1:iroot)
                     end if
                  end if
               end if
            else
               iroot = len_trim (root_name)
               if (iroot .gt. 100) iroot = 100
               ptrk_root(1:iroot) = root_name(1:iroot) 
            end if
            ptrk_name = ''
            num_left = 0
            if (abs(prnt_rst) .ge. 20) then
               ptrk_name(1:iroot) = ptrk_root(1:iroot)
               ptrk_name(iroot+1:iroot+5) = '.ptrk'
               opnd = .false.
               inquire (file = ptrk_name, opened=opnd)
               if (.not. opnd) isptrk1 = open_file(ptrk_name,'unknown')
               if (.not. allocated (string)) allocate (string(nspeci))
               do ith = 1,nspeci
                  string(ith) = ''
                  is = 1
                  ie = 11
                  do i = 1, 6
                     if (prnt_var(i)) then
                        write(string(ith)(is:ie), 101) ith, i
                        is = ie + 1
                        ie = ie + 11
                     end if
                  end do
               end do
               fstring = ''
               varlen = len_trim(string(1))
               write (fstring,104) nspeci, varlen
               write (isptrk1,105)
               write (isptrk1,fstring) 'VARIABLES="Time (days)" ',
     &              (string(ith),ith = 1,nspeci) 
            end if
         end if
      end if
 101  format (' "Sp', i3.3, ' V', i1, '"')
 102  format (1x,i10)
 103  format (g21.14)
 104  format ('(1x, a24,', i3, '(a', i2, '))')
 105  format (1x,'TITLE="V1=Number Having Entered System, ',
     &     'V2=Number Currently In System, ',
     &     'V3=Number Having Left System, ',
     &     'V4=Number Having Decayed, ',
     &     'V5=Number Having Been Filtered, ',
     &     'V6=Number That Left This Time"')
      
      if (days .ge. tims .and. abs(prnt_rst) .ge. 10 .and.
     &     abs(prnt_rst) .lt. 30) then
         allocate ( num_in(n0,nspeci) )
         end_flag = .true.
         num_in = 0
         ptrk_name = ''
         ptrk_name(1:iroot) = ptrk_root(1:iroot)
         ptrk_name(iroot+1:iroot+9) = '.ptrk_fin'
         opnd = .false.
         inquire (file = ptrk_name, opened=opnd)
         if (.not. opnd) isptrk = open_file(ptrk_name,'unknown')
      end if

      if (.not. allocated (string)) allocate (string(nspeci))
      write(tstring,103) days
c calculate the number which have entered the system, as well as 
c the number which have already exited
      do ith=1,nspeci
         string(ith) = ''
         num_enter=0.
         num_leave=0.
         num_decay=0.
         num_filtered=0.
         do i=1,num_particles(ith)
             ret_wght = 1.d0
             if(flag_diversity(ith).and.flag_col_irrev(ith))then
                ret_wght = ret_weight(i,irrevs(divs(ith)))
             endif
             if(flag_col_daughter(ith) .ne. 0)then
                ret_wght = ret_weight_daughter(i,divs_d(ith))
             endif

            if (start_time(i,ith).le.86400.*days) 
     &           num_enter=num_enter+ret_wght
            if (box(i,ith).lt.0)then
               if(timeleft(i,ith).ge.0.)then
                  num_leave=num_leave+ret_wght
               else
                  num_filtered=num_filtered+ret_wght
               end if
            else if (box(i,ith) .gt. 0) then
               if (end_flag) then
                  inode = box(i,ith)
                  num_in(inode,ith) =  
     &                 num_in(inode,ith) + ret_wght
               end if
            else
               num_decay=num_decay+ret_wght
            end if
	 
         enddo

         num_current = num_enter - num_leave - num_decay
         is = 1
         ie = 11
         do i = 1, 6
            if (prnt_var(i)) then
               select case (i)
               case (1)
                  write(string(ith)(is:ie), 102) int(num_enter)
               case (2)
                  write(string(ith)(is:ie), 102) int(num_current)
               case (3)
                  write(string(ith)(is:ie), 102) int(num_leave)
               case (4)
                  write(string(ith)(is:ie), 102) int(num_decay)
               case (5)
                  write(string(ith)(is:ie), 102) int(num_filtered)
               case (6)
                  write(string(ith)(is:ie), 102)
     &                 int(num_leave - num_left(ith))
               end select
               is = ie + 1
               ie = ie + 11
            end if
         end do
         if (ntty.eq.2) then
 1          format(/,19x,'*************************',/,23x,
     &           'Particle Tracking ==> Species: ',i3)
 2          format(1x,'Number Having Entered System: ',i10)
 3          format(1x,'Number Currently In System  : ',i10)
 4          format(1x,'Number Having Left System   : ',i10)
 5          format(1x,'Number Having Decayed       : ',i10)
 6          format(1x,'Number Having Been Filtered : ',i10)
 7          format(1x,'Number That Left This Time  : ',i10)
            if (itmp .gt. 0) then
               write(itmp,1) ith 
               if (prnt_var(1)) write(itmp,2) int(num_enter)
               if (prnt_var(2)) write(itmp,3) int(num_current)
               if (prnt_var(3)) write(itmp,4) int(num_leave)
               if (prnt_var(4)) write(itmp,5) int(num_decay)
               if (prnt_var(5)) write(itmp,6) int(num_filtered)
               if (prnt_var(6)) write(itmp,7) 
     &              int(num_leave - num_left(ith))
               if (m .gt. 0) write(itmp,6010)
            endif

c output for all nodes specified in macro node
!            if (iout .ne. 0 .and. ith .eq. 1) write(iout,6009)
            do i=1,m
               mdd=nskw(i)
               num=0
               do i2=1,num_particles(ith)
                  if (box(i2,ith).eq.mdd .and.
     2                 start_time(i2,ith).le.86400.*days) num=num+1
               enddo
               mdd1=mdd+npt(ith)
               if (iout .ne. 0) write(iout,6011) mdd,an(mdd1),num
!               if (iout .ne. 0) write(iout,6012) ith,mdd,an(mdd1),num
               if (iptty .gt. 0) write(iptty,6011) mdd,an(mdd1),num
            enddo
         endif
! Update num_left to new cumulative
         num_left(ith) = num_leave
      end do
      if (abs(prnt_rst) .ge. 20) then
         write (isptrk1, fstring) tstring, (string(ith),ith = 1,nspeci)
         if (days. ge. tims .and. ripfehm .eq. 0) then
! This was the last time step
            close (isptrk1)
            deallocate (num_left)
         end if
      end if

      if (end_flag) then
         write (isptrk, 6004)
         do ith = 1, nspeci
            write (isptrk, 6003) ith
            do i = 1, n0
               if (i .gt. neq_primary) then
                  n1 = i - neq_primary
               else
                  n1 = i
               end if
               if (num_exit(i,ith) .ne. 0) 
     &              write(isptrk,6005) i, cord(n1,1), cord(n1,2), 
     &              cord(n1,3), int(num_exit(i,ith)), izonef(i)
            end do
         end do
         deallocate (num_in)
         close(isptrk)
      end if

c bhl_6/4/08

      if (abs(prnt_rst) .ge. 40 .and. ripfehm .eq. 1) then
         if (.not. allocated (string2)) allocate (string2(nspeci))
      
         imbl = .false.
         inquire (file = "FEHM_GSM_Mass_balance.txt", opened=imbl)
         if (.not. imbl) imbl1 = 
     &        open_file("FEHM_GSM_Mass_balance.txt",'unknown')

         do ith = 1,nspeci
            string2(ith) = ''
            is = 1
            ie = 16
            do i = 1, 7
               write(string2(ith)(is:ie), 201) ith, i
               is = ie + 1
               ie = ie + 16
            end do
         end do
         fstring2 = ''
         write (fstring2,204) nspeci
         if (ipmbal.le.0) then
            write (imbl1,205)
            write (imbl1,fstring2) 'VARIABLES="Time (years)"',
     &           (string2(ith),ith = 1,nspeci)
            ipmbal=ipmbal+1
         end if
             
 201     format ('      "Sp', i3.3, ' V', i1, '"')
 202     format (7(x,e15.8))
 203     format (g21.14)
 204     format ('(1x, a24,', i3, '(a112))')
 205     format (x,'TITLE="V1=Mass Having Entered System, ',
     &        'V2=Mass Currently In System, ',
     &        'V3=Mass Having Left System, ',
     &        'V4=Mass Having Decayed, ',
     &        'V5=Mass Having Been Filtered, ',
     &        'V6=Mass Having Decayed Outside The UZ, ',
     &        'V7=Filtered Mass Having Decayed"')
c bhl_5/15/08

         write(tstring2,203) days/365.25
         do ith=1,nspeci
            string2(ith) = ''
            mass_enter=0
            mass_leave=0
            mass_decay=0
            mass_filtered=0
            mass_current=0
            mass_dout=0
            mass_dfilt=0
            do ipart=1,num_particles(ith)
               if (start_time(ipart,ith).le.86400.*days) 
     &              mass_enter=mass_enter+gmol(ith)*bconf_sav(ipart,ith)
               if (box(ipart,ith).lt.0)then
                  if(timeleft(ipart,ith).ge.0.)then
                     mass_leave=mass_leave+gmol(ith)*
     &                    bconf_sav(ipart,ith)
                     if(abs(box(ipart,ith)).ge.ibox_offset)then
                        mass_dout=mass_dout+gmol(ith)*
     &                       bconf_sav(ipart,ith)
                     endif
                  else
                     mass_filtered=mass_filtered+gmol(ith)*
     &                    bconf_sav(ipart,ith)
                     if(abs(box(ipart,ith)).ge.ibox_offset)then
                        mass_dfilt=mass_dfilt+gmol(ith)*
     &                       bconf_sav(ipart,ith)
                     endif
                  end if
               
               else if (box(ipart,ith) .gt. 0. and.
     &                 start_time(ipart,ith).le.86400.*days) then
                  mass_current=mass_current+gmol(ith)*
     &                 bconf_sav(ipart,ith)
               else
                  mass_decay=mass_decay+gmol(ith)*bconf_sav(ipart,ith)
               end if
            enddo
            write(string2(ith),202) mass_enter, mass_current, 
     &           mass_leave, mass_decay, mass_filtered,
     &           mass_dout, mass_dfilt
         enddo
     
         write (imbl1,fstring2) tstring2,(string2(ith),ith = 1,nspeci) 
          
      end if

c bhl_6/4/08

      return

 6000 format(10x,20(a20,2x))
 6001 format(4x,'Node',2x,20(a20,2x))
 6002 format(1x,i7,2x,20(i9,2x,i9,2x))
 6003 format('ZONE T = "Species', i3, '"')
 6004 format('VARIABLES = ','"Node"', 2x,'"X"',13x,'"Y"',7x,'"Z"',
     &     4x, '"Num_exited"',2x,'"Zone"')
 6005 format(i7,2x,3(g14.5,1x),1x,i10,2x,i5)
 6009 format(/,1x,'Species',6x,'Node',3x,'Concentration',4x,
     &     '# of Particles')
 6010 format(/,2x,'Node',4x,'Concentration',4x,'# of Particles')
 6011 format(i7,3x,e14.6,8x,i10)
 6012 format(3x,i3.3,4x,i7,1x,e14.6,8x,i10)
      end
