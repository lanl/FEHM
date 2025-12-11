      subroutine diskp(read_ptrk)
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
CD1  This subroutine writes or reads particle tracking information to/
CD1  from the restart file.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/diskp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jan 10 11:05:20 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.5   Wed Jan 10 08:31:34 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   05/16/95 09:34:36   robinson
CD2 Added option to write partial restart info to .fin file
CD2 
CD2    Rev 1.3   04/25/95 09:35:10   llt
CD2 retrieved lost log history
CD2 
CD2    Rev 1.2   03/15/95 17:04:52   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.1   02/02/95 15:22:28   llt
CD2 added pvcs log info
CD2
CD2    Rev 1.0   01/28/95 14:00:20   llt
CD2 new particle tracking module
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3 2.7 Provide Restart Capability
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

      use comai
      use combi
      use comdi, only : nspeci
      use comdti
      use compart
      use comrxni
      implicit none

      integer isolflg, i, idum, nwds, j
      real*8 dumm
      character*80 dummy_string
      integer imsg(4),msg(4)
      real*8  xmsg(4)
      character*32 cmsg(4)
      character*11 ptrkword
      logical read_ptrk
            
! zvd 22-Jan-04 Fluxes are no longer read or written in diskp

      if (nsave .eq. 0 ) then
         isolflg = 0
         if (wdd1(5:8) .eq. 'ptrk') then
            if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
               read(iread, end = 10, err = 10) dummy_string
            else
               read(iread, *, end = 10, err = 10) dummy_string
            end if
            isolflg=1
         else if (read_ptrk) then
            if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
               read(iread, end = 10, err = 10) ptrkword
               read(iread, end = 10, err = 10) dummy_string
            else
               read(iread, *, end = 10, err = 10) ptrkword
               read(iread, *, end = 10, err = 10) dummy_string
            end if
            if (ptrkword .eq. 'ptrk') isolflg=1
         end if
 10      continue
         if (isolflg.eq.0) then
c     either eof found or trac data is in the file instead of particle
c     tracking data, so we continue without reading stuff in
            if (iout .ne. 0) write(iout,6000)
            if (iatty.gt.0) write(iatty,6000)
 6000       format('Particle tracking data not found in restart file')
         else if (.not. read_ptrk) then
            write(ierr, 6001)
            if (iout .ne. 0) write(iout, 6001)
            if (iptty .gt. 0 ) write(iptty, 6001)
 6001       format('Particle tracking data found in restart file will',
     &           ' not be used')
         else
c     read in data
            call parse_string(dummy_string,imsg,msg,xmsg,cmsg,nwds)
            num_particles(1) = imsg(1)+xmsg(1)
            rseed = imsg(2)+xmsg(2)
            if (rseed .gt. 0) then
               rseed_release=imsg(3)+xmsg(3)
               nspeci = imsg(4)+xmsg(4)
            else
               rseed_release=abs(rseed)
               nspeci = imsg(3)+xmsg(3) 
            end if
            if (nspeci .gt. 1) read(iread) (num_particles(i), i = 1,
     &           nspeci)
            box = 0
            frac_done = 0.
            theta = 1.
            timeleft = 0.
            if( rseed .gt. 0 ) then
               if (ischk .ne. 0) then
                  write(ischk,*)
                  write(ischk,*)'Cell locations, fractions done, and ',
     &                 'times read for '
                  write(ischk,*) num_particles(1),
     &                 ' particles from restart file.'
                  write(ischk,*)
               end if
               if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
                  do j = 1, nspeci
                     read(iread) (box(i,j),i=1,num_particles(j))
                     read(iread) (frac_done(i,j),i=1,num_particles(j))
                     read(iread) (theta(i,j),i=1,num_particles(j))
                     read(iread) (timeleft(i,j),i=1,num_particles(j))
                  end do
               else
                  do j = 1, nspeci
                     read(iread,*) (box(i,j),i=1,num_particles(j))
                     read(iread,*) (frac_done(i,j),i=1,num_particles(j))
                     read(iread,*) (theta(i,j),i=1,num_particles(j))
                     read(iread,*) (timeleft(i,j),i=1,num_particles(j))
                  end do
               endif
            else
               if (ischk .ne. 0) then
                  write(ischk,*)
                  write(ischk,*)'Warning - restart does not contain all'
                  write(ischk,*)'particle information needed for exact'
                  write(ischk,*)'restart.  Particle tracking simulation'
                  write(ischk,*)'will not be exact. Times spent by each'
                  write(ischk,*)'particle in its restart cell will be'
                  write(ischk,*)'inaccurate.'
                  write(ischk,*)'Only the cell locations and times were'
                  write(ischk,*)'read for ', num_particles(1),
     &                 ' particles from restart file.'
                  write(ischk,*)
               end if
               rseed = -rseed
               if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
                  do j = 1, nspeci
                     read(iread) (box(i,1),i=1,num_particles(1))
                     read(iread) (timeleft(i,1),i=1,num_particles(1))
                  end do
               else
                  do j = 1, nspeci
                     read(iread,*) (box(i,j),i=1,num_particles(j))
                     read(iread,*) (timeleft(i,j),i=1,num_particles(j))
                  end do
               endif
c               do i = 1, num_particles(1)
c                  frac_done(i,1) = 0.
c                  theta(i,1) = 1.
c               end do
            end if
c     fix data
            do j = 1, nspeci
               do i=1,num_particles(j)
                  if (box(i,j).gt.0) then
                     start_time(i,j)=86400.*days-timeleft(i,j)
                  endif
               end do
            end do
            restarting = .TRUE.
         endif
         
      else if ((isave.gt.0)) then
c zvd 6/9/08 prnt_rst values of 41, 42 for GoldSim mass output option
         ptrkword = 'ptrk       '
         if (bin_flag .eq. 1 .or. bin_flag .eq. 3) then
            select case (prnt_rst)
            case (1, 2, 11, 12, 21, 22, 31, 32, 41, 42)
               if (header_flag .eq. 'new') write(isave) ptrkword
               write(dummy_string, *)  num_particles(1), rseed, 
     &              rseed_release, nspeci
               write(isave) dummy_string
               if (nspeci .gt. 1) write(isave) (num_particles(i),
     &              i = 1, nspeci)
               do j = 1, nspeci
                  write(isave) (box(i,1),i=1,num_particles(j))
                  write(isave) (frac_done(i,1),i=1,num_particles(j))
                  write(isave) (theta(i,1),i=1,num_particles(j))
                  write(isave) (timeleft(i,1),i=1,num_particles(j))
               end do
            case (-1, -2, -11, -12, -21, -22, -31, -32, -41, -42)
               if (header_flag .eq. 'new') write(isave,*) ptrkword
               write(dummy_string, *) num_particles(1), -rseed, nspeci
               write(isave) dummy_string
               if (nspeci .gt. 1) write(isave) (num_particles(i),
     &              i = 1, nspeci)
               do j = 1, nspeci
                  write(isave) (box(i,j),i=1,num_particles(j))
                  write(isave) (timeleft(i,j),i=1,num_particles(j))
               end do
            case default
! Do nothing, particle data is not written to restart file
            end select
         else
            select case (prnt_rst)
            case (1, 2, 11, 12, 21, 22, 31, 32, 41, 42)
               if (header_flag .eq. 'new') write(isave,'(a4)') 'ptrk'
               write(isave,*) num_particles(1), rseed, rseed_release,
     &              nspeci
               if (nspeci .gt. 1) write(isave, *) (num_particles(i),
     &              i = 1, nspeci)
               do j = 1, nspeci
                  write(isave, 20) (box(i,j),i=1,num_particles(j))
                  write(isave, 30) (frac_done(i,j),i=1,num_particles(j))
                  write(isave, 30) (theta(i,j),i=1,num_particles(j))
                  write(isave, 30) (timeleft(i,j),i=1,num_particles(j))
               end do
            case (-1, -2, -11, -12, -21, -22, -31, -32 -41, -42)
               if (header_flag .eq. 'new') write(isave,'(a4)') 'ptrk'
               write(isave,*) num_particles(1), -rseed, nspeci
               if (nspeci .gt. 1) write(isave, *) (num_particles(i),
     &              i = 1, nspeci)
               do j = 1, nspeci
                  write(isave, 20) (box(i,j),i=1,num_particles(j))
                  write(isave, 30) (timeleft(i,j),i=1,num_particles(j))
               end do
            case default
! Do nothing, particle data is not written to restart file
            end select
        end if 
      endif
      
 20   format (10i10)
 30   format (4g25.16)
      return
      end
