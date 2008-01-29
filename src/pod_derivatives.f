      subroutine pod_derivatives
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
!D1 Calculate derivative functions needed for POD basis functions.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 09-SEP-04, Programmer: Z. Dash, B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/pod_derivatives.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************
      
      use comai
      use comci
      use comdi
      use comsptr
      use comxi
      implicit none

      integer coord, i, j, ispod, iroot
      real*8, allocatable :: fpod(:,:), gpod(:,:), drho(:), velocity(:)
      real*8, allocatable :: Dx(:), Dy(:), Dz(:)
      real*8 du, d2rho, durho, dphiSD, dispersivity,d(3)
      real*8 grad_array, velocities
      character*120 fname, root
      logical null1
      
      if (null1(root_name)) then
! Find root of output file name or if not present, input file name
         if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
            call file_prefix(nmfil(5), iroot)
            if (iroot .gt. 100) iroot = 100
            root(1:iroot) = nmfil(5)(1:iroot)
         else
            if (nmfil(2)(1:1) .eq. ' ' ) then
               write (ierr, *) 'FILE ERROR: nmfil2 file: ', nmfil(2),
     .              ' unable to determine pod file prefix'
               stop
            else
               call file_prefix(nmfil(2), iroot)
            end if
         end if
      else
         iroot = len_trim (root_name)
         if (iroot .gt. 100) iroot = 100
         root(1:iroot) = root_name(1:iroot)
      end if
      fname =  root(1:iroot) // '.pod'
      ispod = 321
      open (unit=ispod, file=fname, form='formatted')
      write(ispod, 6000)  verno, jdate, jtime, wdd

      allocate (fpod(neq,3), gpod(neq,3), drho(neq), velocity(neq))
      allocate (Dx(neq), Dy(neq), Dz(neq))

! Load dispersivity arrays
      do i = 1, neq
         call dispersion_node(i, d)
         Dx(i) = d(1)
         Dy(i) = d(2)
         Dz(i) = d(3)
      end do

! x = 1, y = 2, z = 3
      do coord = 1, 3
         do i = 1, neq
! Derivative of density
            drho(i) = grad_array(coord, rolf, i)
            velocity(i) = velocities(coord, i)
         end do
         do i = 1, neq
! 2nd derivative of density
            d2rho = grad_array(coord, drho, i)
! derivative of velocity
            du =  grad_array(coord, velocity, i)
            durho =  rolf(i) * du + velocity(i) * drho(i)
!            call ?dispersion_derivative?
! 3 lines below temporary to force zero derivatives
            durho = 0.
            dphiSD = 0.
            dispersivity = 0.

            fpod(i, coord) = (drho(i) *  dphiSD - durho) / 
     &           (ps(i) * s(i) * rolf(i)) +
     &           dispersivity * d2rho / rolf(i)

            gpod(i, coord) = dphiSD / (ps(i) * s(i)) + 
     &           2*dispersivity * drho(i) / rolf(i) -
     &           velocity(i) / (ps(i) * s(i))
            
         end do
      end do

! Write out the derivative functions and dispersivities
      do i = 1, 3
         write (ispod, *) ( fpod(j, i), j = 1, neq)
      end do
      do i = 1, 3
         write (ispod, *) (gpod(j, i), j = 1, neq)
      end do
      write (ispod, *) (Dx(i), i = 1, neq)
      write (ispod, *) (Dy(i), i = 1, neq)
      write (ispod, *) (Dz(i), i = 1, neq)

      deallocate (fpod, gpod, drho, velocity, Dx, Dy, Dz)
      close (ispod)

 6000 format(a30, 3x, a8, 3x, a8, /, a80, /)
      end
