      subroutine DoubleDot(A, B, res)
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
C**********************************************************************
!     Computes the tensor product of two tensors A and B
!     res = A:B
!
!     Author : Sai Rapaka
!
      implicit none
      real*8, dimension(6)         :: A, B
      real*8                       :: res

      res = sum(A(1:3)*B(1:3)) + 2.0d0*sum(A(4:6)*B(4:6))

      end subroutine DoubleDot
