      subroutine  flxn
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
!D1 To write source/sink fluxes by node.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 15-JAN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/flxn.f_a  $
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

      use comai, only : neq, days, icnl, iflxn
      use combi, only : cord, nelmdg
      use comflow, only : a_axy, a_vxy
      use davidi, only : irdof

      integer i, index_axy, iflux
      real*8 a_vxy_value
      logical it_is_open

! Find a unit that is not used
      do 
         iflux = 100
         inquire (unit = iflux, opened = it_is_open)
         if (.not. it_is_open) exit
         iflux = iflux + 1
      end do

      open (unit = iflux, file = "source_sink.flux")

      write (iflux, *) '"Time = ', days, 'days (source < 0, sink > 0)"'
      write (iflux, *) '"X [m]"  "Y [m]"  "Z [m]"  "Node"',
     &     '  "Liquid Flux (kg/s)"  "Vapor Flux (kg/s)"'
      do i = 1, neq
         index_axy = nelmdg (i) - (neq + 1)
         if (irdof .ne. 13) then
            a_vxy_value = a_vxy (index_axy)
         else
            a_vxy_value = 0.d0
         end if
         if (a_axy (index_axy) .ne. 0. .or. a_vxy_value .ne. 0.) then
            select case (icnl)
            case (0)
               write (iflux, 100) cord(i,1), cord(i,2), cord(i,3),
     &              i, a_axy (index_axy), a_vxy_value
            case (1, 4)
               write (iflux, 100) cord(i,1), cord(i,2), 0.,
     &              i, a_axy (index_axy), a_vxy_value
            case (2, 5)
               write (iflux, 100) cord(i,1), 0., cord(i,3),
     &              i, a_axy (index_axy), a_vxy_value
            case (3, 6)
               write (iflux, 100) 0., cord(i,2), cord(i,3),
     &              i, a_axy (index_axy), a_vxy_value
            end select
         end if
      end do
      close (iflux)
 100  format (3(g16.9, 1x), i8, 2(1x, g16.9))

      end
