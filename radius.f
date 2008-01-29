      subroutine  radius
!***********************************************************************
!  Copyright, 1994, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To determine the radius in radial geometry problems.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 $Log:   /pvcs.config/fehm90/src/radius.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Wed Jan 17 12:41:36 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.2   Thu Jan 11 09:43:32 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:43:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:40   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation 
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5 
CD5 None
CD5
CD5 Interface Tables
CD5
CD5 None
CD5
CD5 Files
CD5 
CD5 None
CD5
C***********************************************************************CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6 
CD6 neq, icnl, cord, nelm, sx1
CD6
CD6 Global Subprograms
CD6
CD6 None
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   i               int      Do loop index
CD7   i1              int      Do loop limit
CD7   i2              int      Do loop limit
CD7   j               int      Do loop index
CD7   kb              int      Node number index
CD7   pi2             real*8   Pi x 2
CD7   radi            real*8   Radius
CD7   radmin          real*8   Minimum radius
CD7   radmax          real*8   Maximum radius
CD7   radkb           real*8   Average radius between i and kb
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CD8
CD8 ASSUMPTIONS AND LIMITATIONS
CD8
CD8 This routine assumes it is not being called for a 3-D problem
CD8 (icnl=0).
CD8
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN radius
CPS
CPS FOR each node
CPS
CPS   IF this is a 2-D cartesian coordinate problem
CPS     Set coordinate value of third coordinate to 1
CPS   ELSE this is a radial problem
CPS
CPS     Set parameter values
CPS
CPS     FOR each node connected to this node
CPS       EXITIF there is no connected node stored in the array ...
CPS       ... at this position
CPS       Computed average radius, new minimum and maximum radii.
CPS     ENDFOR each node connected to this node
CPS
CPS     Compute volume associated with this node
CPS     Store value used later in calculation in coordinate array
CPS
CPS   ENDIF
CPS
CPS ENDFOR each node
CPS
CPS END radius
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      integer i,i1,i2,j,kb,n_n_n,imodel
      real*8 pi2,radi,radmin,radmax,radkb
      real*8 rad_tol
      parameter (rad_tol=1.d-10)

      pi2    =  6.283185
      do i = 1, neq_primary
         if ( icnl .le. 3 ) then
            cord(i,3) =  1.0
         else
            radi   =  cord(i,1)
            radmin =  radi
            radmax =  radi
            i1     =  nelm(i  )+1
            i2     =  nelm(i+1)
            do j = i1, i2
               kb     =  nelm(j)
               if ( kb .eq. 0 )  go  to  1000
               radkb  =  0.5*( radi+cord(kb,1) )
               radmax =  max( radkb,radmax )
               radmin =  min( radkb,radmin )
            end do
 1000       continue
            if((radmin+radmax)/2..lt.rad_tol) then
               write (ierr,900) i,(radmin+radmax)/2.,rad_tol
               if (iout .ne. 0) 
     &              write(iout,900) i,(radmin+radmax)/2.,rad_tol
               if (iptty.ne.0)  
     &              write(iptty,900) i,(radmin+radmax)/2.,rad_tol
 900           format('>>> radius node ',i9,' lt rad_tol (stopping)',
     &              /,'radius = ',g12.4,' rad_tol = ',g12.4)
               stop
            endif
            sx1 (i  ) =  sx1(i)*0.5*( radmin+radmax )*pi2
            cord(i,3) =         0.5*( radmin+radmax )*pi2
         end if
      end do
      
      return
      end
