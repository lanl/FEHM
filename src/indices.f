      subroutine indices (low, high, datav, datar, numr) 
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
!D1 Find indices of appropriate sigma and omega curves to be used in 
!D1 interpolation (log value method).
!D1  
!***********************************************************************
!D2 
!D2 REVISION HISTORY
!D2 
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 Initial implementation: 05-FEB-99, Programmer: Z. Dash
!D2 Modification for data revision: 16-AUG-99
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/indices.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:30:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:10   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
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
 
      implicit none

      real*8 datav, datar(*)
      integer low, high, numr, midpoint

C Find the T vs C curve indices

      if (datav .le. datar(1)) then
         low = 1
         high = low
         return
      else if (datav .ge. datar(numr)) then
         low = numr
         high = low
         return
      else
         low = 1
         high = numr
 100     midpoint = (low + high) / 2
         if (datav .eq. datar(midpoint)) then
            low = midpoint
            high = low
            return
         else if (datav .lt. datar(midpoint)) then
            high = midpoint
         else
            low = midpoint
         endif
         if (high - low .eq. 1) return
         goto 100
      endif

      end
