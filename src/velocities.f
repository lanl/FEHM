      real*8 function velocities(coord, i)
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
!D1 Find velocity values needed for POD basis function derivatives.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 08-SEP-04, Programmer: B. Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/dispersion_node.f_a  $
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

      use comdi, only : ps
      use comsptr, only : ggg
      
      implicit none
      real*8 velocity1, velocity2
      integer i,ix1,ixm1
      integer coord
      
c account for boundary nodes by taking the other connection if
c     there is no connection at that location
c     Then, we average the two velocities to do either a centered
c     difference or just get the one that exists

      velocity1 = -ggg(i,-coord)*ps(i)
      velocity2 = +ggg(i,coord)*ps(i)
      velocities = (velocity1 + velocity2)/2.

      return
      
      end

