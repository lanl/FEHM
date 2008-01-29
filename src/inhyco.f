      subroutine inhyco
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
!D1 Read hydraulic conductivity data.                   
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 09-FEB-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inhyco.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
      use comdti
      use comki
      implicit none

      integer i
      real*8 conv 

      macro = 'hyco'
c**** read hydraulic conductivity data (in m/sec)  
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = zero_t
      default(2) = zero_t
      default(3) = zero_t
      igroup = 1

      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(14), macro, igroup, ireturn,
     2     r8_1=pnx(1:n0),r8_2=pny(1:n0),r8_3=pnz(1:n0)) 

      do i=1,n0
         if(pnx(i).ge.0.0d00) then
            pnx(i) = max (zero_t, pnx(i))
         else
            pnx(i) = 1.0d00/10.0d00**(abs(pnx(i)))
         endif
         if(pny(i).ge.0.0d00) then
            pny(i) = max (zero_t, pny(i))
         else
            pny(i) = 1.0d00/10.0d00**(abs(pny(i)))
         endif
         if(pnz(i).ge.0.0d00) then
            pnz(i) = max (zero_t, pnz(i))
         else
            pnz(i) = 1.0d00/10.0d00**(abs(pnz(i)))
         endif
      end do
c convert to intrinsic perms at 20 C and 0.1 Mpa
      conv = 9.89d-13/9.66d-6
      do i=1,n0
         pnx(i) = conv*pnx(i)
         pny(i) = conv*pny(i)
         pnz(i) = conv*pnz(i)
      end do

      macroread(14) = .TRUE.
   

      end
