      subroutine infrlp
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
!D1 To read relative permeability information.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: 26-MAR-2000, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/infrlp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.4 Relative-permeability and capillary-pressure functions
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
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

      use comdi
      use comdti
      use comai
      use comki
      implicit none

      real*8, allocatable ::  tmp(:)
      integer i
      
      macro = 'frlp'
      
      allocate(tmp(n0))

c**** read relative permeability factors ****
      narrays = 2
      itype(1) = 8
      itype(2) = 8
      default(1) = 0.0d00
      default(2) = 0.0d00
      igroup = 1
      
      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(17), macro, igroup, ireturn,
     2     r8_1=rlp_fac(1:n0),r8_2=tmp(1:n0)) 
      do i=1,n0            
         rlp_fac(n0+i) = tmp(i)
      enddo  
      
      macroread(17) = .TRUE.
      
      deallocate(tmp)

      end
