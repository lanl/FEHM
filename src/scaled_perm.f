      subroutine scaled_perm
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Scale permeability data.
!D1 
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: ?, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/scaled_perm.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Ouput Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comai
      use comdi
      use comdti
      use comki

      implicit none
      real*8, allocatable :: scalex(:)
      real*8, allocatable :: scaley(:)
      real*8, allocatable :: scalez(:)

      integer i, ir
 
      allocate(scalex(n0))
      allocate(scaley(n0))
      allocate(scalez(n0))

      scalex = 1.0d00
      scaley = 1.0d00
      scalez = 1.0d00

      macro = 'fper'

c**** read permeability scaling data ****
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = 1.0d00
      default(2) = 1.0d00
      default(3) = 1.0d00
      igroup = 1

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(22), macro, igroup, ireturn,
     2     r8_1=scalex(1:n0),r8_2=scaley(1:n0),
     3     r8_3=scalez(1:n0)) 

      do i=1,n0
         pnx(i) = scalex(i)*pnx(i)       
         pny(i) = scaley(i)*pny(i)       
         pnz(i) = scalez(i)*pnz(i)       
      end do

      macroread(22) = .TRUE.

      if (nrlp .ne. 0) then
c     Assign scaling factors for rlp models
         if (.not. allocated(xfperm)) then
c     use fperm arrays for temporary storage for later model assignment
            allocate(xfperm(n0),yfperm(n0),zfperm(n0))
            xfperm = scalex
            yfperm = scaley
            zfperm = scalez
         else
c     If rlp has already been read, these arrays have been allocated
c     and a model assignment needs to be made
            do i = 1, n0
               ir = irlp(i)
               if (ir  .ne. 0) then
                  xfperm(ir) = scalex(i)
                  yfperm(ir) = scaley(i)
                  zfperm(ir) = scalez(i)
               end if
            end do
         end if
      end if

      deallocate(scalex)
      deallocate(scaley)
      deallocate(scalez)

      end
