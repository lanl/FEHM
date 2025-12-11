      subroutine pcp_save(iflg,n,ndummy,j,pcp,dpcps,rlf,drlfs,
     &                    rvf,drvfs,s,pcp0,dpcps0,rlf0,drlfs0,
     &                    rvf0,drvfs0,s0)
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
!D1 To hold transmissibilities constant if phase change occurs.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/pcp_save.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:36   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of thrair 08-Feb-02
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.1 Pressure- and temperature-dependent water properties
!D3 2.4.2 Properties of air and air/water vapor mixtures
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

      implicit none

      integer iflg, n, i, j, ii, ndummy
      real*8 pcp(*),  dpcps(*), s(*), rlf(*), drlfs(*)
      real*8 rvf(*),  drvfs(*)
      real*8 pcp0(*),  dpcps0(*), s0(*)
      real*8 rlf0(*),  drlfs0(*)
      real*8 rvf0(*),  drvfs0(*)
      if(iflg.eq.0) then
       do ii=1,n
        i=ii+ndummy
         rlf0(i)=rlf(i)
         drlfs0(i)=drlfs(i)
         rvf0(i)=rvf(i)
         drvfs0(i)=drvfs(i)
         dpcps0(i)=dpcps(i)
         pcp0(i)=pcp(i)
         dpcps0(i)=dpcps(i)
         s0(i)=s(i)
       enddo
      else if(iflg.eq.1) then
         pcp(j)=pcp0(j)+dpcps0(j)*(s(j)-s0(j))
         dpcps(j)=dpcps0(j)
         rlf(j)=rlf0(j)+drlfs0(j)*(s(j)-s0(j))
         drlfs(j)=drlfs0(j)
         rvf(j)=rvf0(j)+drvfs0(j)*(s(j)-s0(j))
         drvfs(j)=drvfs0(j)
      else if(iflg.eq.2) then
       do ii=1,n
        i=ii+ndummy
         pcp(i)=pcp0(i)+dpcps0(i)*(s(i)-s0(i))
         dpcps(i)=dpcps0(i)
         rlf(i)=rlf0(i)+drlfs0(i)*(s(i)-s0(i))
         drlfs(i)=drlfs0(i)
         rvf(i)=rvf0(i)+drvfs0(i)*(s(i)-s0(i))
         drvfs(i)=drvfs0(i)
       enddo
      endif
c
      return
      end
