      subroutine  inrestart
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read control information for reading and writing restart files.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inrestart.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3 2.7 Provide Restart Capability
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
      use combi, only : rstr, rstw
      use comxi
      implicit none

      integer i, imsg(8), msg(8), nwds
      real*8  xmsg(8)
      character*32 cmsg(8)
      character*80 chdum
      logical null1

! Default will be to read or write all variables from/to restart file 
! (based on problem type) using new header style
      header_flag = 'new'

! While loop to read restart control information
      do
         chdum = ''
         read (inpt, '(a80)') chdum
         if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or. 
     &        chdum(1:3) .eq. 'END') exit
         flags : select case (chdum(1:4))
! ascii - both read and write files are ascii
! binary - both read and write files
! rbinary - binary read file
! wbinary - binary write file
         case ('asci', 'ASCI')
            bin_flag = 0
! Default format is formatted
         case ('bina', 'BINA')
            bin_flag = 1
! Restart input
            cform(6) = 'binary'
! Restart output
            cform(7) = 'binary'
         case ('rbin', 'RBIN')
            bin_flag = 2
! Restart input
            cform(6) = 'binary'
         case ('wbin', 'WBIN')
            bin_flag = 3
! Restart output
            cform(7) = 'binary'
! overwrite - only 1 set of restart data is saved
! append - multiple sets of restart data are written (not implemented)
         case ('over')
            app_flag = 0
         case ('appe')
            app_flag = 1
         case ('read')
! List parameters to be read from restart file
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            do i = 2, nwds
               if (msg(i) .eq. 3) rstr(i-1) = cmsg(i)(1:4)
            end do
         case ('writ')
! List parameters to be written to restart file
            call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
            do i = 2, nwds
               if (msg(i) .eq. 3) rstw(i-1) = cmsg(i)(1:4)
            end do
! noflux - no flux output
! flux - both liquid and vapor
! lflux - liquid flux only
! vflux - vapor flux only
         case ('nofl', 'NOFL')
            flux_flag = 'no fluxes  '
         case ('flux', 'FLUX')
            flux_flag = 'all fluxes '
         case ('lflu', 'LFLU')
            flux_flag = 'liquid flux'
         case ('vflu', 'VFLU')
            flux_flag = 'vapor flux '
! fluxes and flow source/sink flow for mixing model
         case ('rtdm', 'RTDM')
            call gen_mixmodel
! new - restart file will be written using new header format
         case ('new ', 'NEW ')
            header_flag = 'new'
! old - restart file will be written using old header format
         case ('old ', 'OLD ')
            header_flag = 'old'
         case default
            write (ierr, 100) chdum
         end select flags
      enddo

 100  format ('WARNING: Unrecognized restart option: ', a, /, 
     &     'Check input deck')
      end
      
