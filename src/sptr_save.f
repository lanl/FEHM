      subroutine sptr_save (sflag)
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
!D1 Perform initial setup functions for the streamline particle 
!D1 tracking calculations. 
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/sptr_save.f_a  $
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
      use comdi
      use comsk
      use comsptr
      use comxi

      integer sflag, i, j, k

      open (isptr7, file = nmfil(23), form = cform(23), 
     .     status = cstats(23))

      if (sflag .eq. 0) then
! Read stored sptr arrays
         if (cform(23) .eq. 'formatted' .or. 
     .        cform(23) .eq. 'FORMATTED') then
            read (isptr7, *) node_count
            if (node_count .ne. 0) then
               read (isptr7, *) (node_omr(i), i = 1, node_count)
               do j = 1, 7
                  read (isptr7, *) (iboulist(i,j), i=1,neq)
               end do
            end if
            do j = -3, 3
               read (isptr7, *) (irray(i,j), i=1,neq)
            end do
            read (isptr7, *) (ddx(i), i=1,neq)
            read (isptr7, *) (ddy(i), i=1,neq)
            read (isptr7, *) (ddz(i), i=1,neq)
            do j = 1, 3
               read (isptr7, *) (corn(i, j), i=1,neq)
            end do
            if (.not. compute_flow) then
               do j = -3, 3
                  read (isptr7, *) (ggg(i, j), i=1,neq)
               end do
            end if
         else
            read (isptr7) node_count
            if (node_count .ne. 0) then
               read (isptr7) (node_omr(i), i = 1, node_count)
               do j = 1, 7
                  read (isptr7) (iboulist(i,j), i=1,neq)
               end do
            end if
            do j = -3, 3
               read (isptr7) (irray(i,j), i=1,neq)
            end do
            read (isptr7) (ddx(i), i=1,neq)
            read (isptr7) (ddy(i), i=1,neq)
            read (isptr7) (ddz(i), i=1,neq)
            do j = 1, 3
               read (isptr7) (corn(i, j), i=1,neq)
            end do
            if (.not. compute_flow) then
               do j = -3, 3
                  read (isptr7) (ggg(i, j), i=1,neq)
               end do
            end if
         end if

      else
! Save sptr arrays
         if (cform(23) .eq. 'formatted' .or. 
     .        cform(23) .eq. 'FORMATTED') then
            write (isptr7, *) node_count
            if (node_count .ne. 0) then
               write (isptr7, *) (node_omr(i), i = 1, node_count)
               do j = 1, 7
                  write (isptr7, *) (iboulist(i,j), i=1,neq)
               end do
            end if
            do j = -3, 3
               write (isptr7, *) (irray(i, j), i=1,neq)
            end do
            write (isptr7, *) (ddx(i), i=1,neq)
            write (isptr7, *) (ddy(i), i=1,neq)
            write (isptr7, *) (ddz(i), i=1,neq)
            do j = 1, 3
               write (isptr7, *) (corn(i, j), i=1,neq)
            end do
            if (.not. compute_flow) then
               do j = -3, 3
                  write (isptr7, *) (ggg(i, j), i=1,neq)
               end do
               write(isptr7,*)(ps_trac(i),i=1,neq)
            end if
         else
            write (isptr7) node_count
            if (node_count .ne. 0) then
               write (isptr7) (node_omr(i), i = 1, node_count)
               do j = 1, 7
                  write (isptr7) (iboulist(i,j), i=1,neq)
               end do
            end if
            do j = -3, 3
               write (isptr7) (irray(i, j), i=1,neq)
            end do
            write (isptr7) (ddx(i), i=1,neq)
            write (isptr7) (ddy(i), i=1,neq)
            write (isptr7) (ddz(i), i=1,neq)
            do j = 1, 3
               write (isptr7) (corn(i, j), i=1,neq)
            end do
            if (.not. compute_flow) then
               do j = -3, 3
                  write (isptr7) (ggg(i, j), i=1,neq)
               end do
               write(isptr7)(ps_trac(i),i=1,neq)
            end if
         end if
      end if

      close (isptr7)
      
      end
