      real*8 function interp8 (a1, a2, av, b1, b2, bv, c1, c2, cv,
     .     t1, t2, t3, t4, t5, t6, t7, t8)
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
!D1 To calculate particle time delay using type curves.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 13-NOV-02, Programmer: Zora Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/interp8.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************
! Return the interpolated time value

      implicit none

      real*8 a1, a2, av, b1, b2, bv, c1, c2, cv
      real*8 t1, t2, t3, t4, t5, t6, t7, t8
      real*8 d1, d2, d3, d4, d5, d6, d7, d8, totald

      d1 = dsqrt( (dlog10 (a1) - dlog10(av))**2 +
     .     (dlog10 (b1) - dlog10(bv))**2 + 
     .     (dlog10 (c1) - dlog10(cv))**2)

      d2 = dsqrt( (dlog10 (a1) - dlog10(av))**2 +
     .     (dlog10 (b1) - dlog10(bv))**2 + 
     .     (dlog10 (c2) - dlog10(cv))**2)

      d3 = dsqrt( (dlog10 (a1) - dlog10(av))**2 +
     .     (dlog10 (b2) - dlog10(bv))**2 + 
     .     (dlog10 (c1) - dlog10(cv))**2)

      d4 = dsqrt( (dlog10 (a1) - dlog10(av))**2 +
     .     (dlog10 (b2) - dlog10(bv))**2 + 
     .     (dlog10 (c2) - dlog10(cv))**2)

      d5 = dsqrt( (dlog10 (a2) - dlog10(av))**2 +
     .     (dlog10 (b1) - dlog10(bv))**2 + 
     .     (dlog10 (c1) - dlog10(cv))**2)

      d6 = dsqrt( (dlog10 (a2) - dlog10(av))**2 +
     .     (dlog10 (b1) - dlog10(bv))**2 + 
     .     (dlog10 (c2) - dlog10(cv))**2)

      d7 = dsqrt( (dlog10 (a2) - dlog10(av))**2 +
     .     (dlog10 (b2) - dlog10(bv))**2 + 
     .     (dlog10 (c1) - dlog10(cv))**2)

      d8 = dsqrt( (dlog10 (a2) - dlog10(av))**2 +
     .     (dlog10 (b2) - dlog10(bv))**2 + 
     .     (dlog10 (c2) - dlog10(cv))**2)

      if (d1 .eq. 0) then
         interp8 = t1
         return
      end if

      if (d2 .eq. 0) then
         interp8 = t2
         return
      end if

      if (d3 .eq. 0) then
         interp8 = t3
         return
      end if

      if (d4 .eq. 0) then
         interp8 = t4
         return
      end if

      if (d5 .eq. 0) then
         interp8 = t5
         return
      end if

      if (d6 .eq. 0) then
         interp8 = t6
         return
      end if

      if (d7 .eq. 0) then
         interp8 = t7
         return
      end if

      if (d8 .eq. 0) then
         interp8 = t8
         return
      end if

      d1 = 1./d1
      d2 = 1./d2
      d3 = 1./d3
      d4 = 1./d4
      d5 = 1./d5
      d6 = 1./d6
      d7 = 1./d7
      d8 = 1./d8

      totald = 1. /(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8)

      interp8 = totald * (d1*t1 + d2*t2 + d3*t3 + d4*t4
     .     + d5*t5 + d6*t6 + d7*t7 + d8*t8)

      return
      end
