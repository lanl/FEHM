      subroutine  near3 (x, y, z, node, iflg)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
CD1 Determine the nearest node to a given set of coordinates (x, y, z).
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 18-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/near3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Tue Jan 30 16:45:34 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:55:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:46   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   x               REAL*8   I    X coordinate position
CD3   y               REAL*8   I    Y coordinate position
CD3   z               REAL*8   I    Z coordinate position
CD3   node            INT      O    Number of node nearest to coordinate
CD3                                   position
CD3   iflg            INT      I    Flag to indicate if zones are defined
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   icnl            INT      faai   Problem dimension
CD4   izonef          INT      fbb    Zone in which each node is locatedCD4   n0              INT      param  Maximum number of nodes allowed
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4
CD4 Global Subprograms
CD4
CD4   None
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   dismin          REAL*8   Minimum distance between given coordinates and
CD5                              a node
CD5   dist            REAL*8   Distance between given coordinates and node
CD5   i               INT      Loop index
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This routine is used to allow another type of input specification
CD8 (zone by coordinates), a minor modification solely for user ease
CD8 of use.
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 Not Applicable.  See Special Comments.
CD9
C*******************************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN  near3
CPS 
CPS   initialize minimum distance between nodes to a large value
CPS   
CPS   IF zones are not defined
CPS   
CPS      FOR each node
CPS          IF the problem is 3D
CPS             compute the distance between the node and the given x,y 
CPS              z coordinates
CPS          ELSE
CPS             compute the distance between the node and the given x,y 
CPS              coordinates
CPS          ENDIF
CPS          
CPS          IF the distance between the current node and the 
CPS           coordinates 
CPS           is less than the minimum distance
CPS             set minimum distance to current distance
CPS             set node to the current node
CPS          ENDIF
CPS          
CPS      ENDFOR    
CPS          
CPS   ELSE
CPS   
CPS      FOR each node
CPS          IF the problem is 3D
CPS             compute the distance between the node and the given x,y 
CPS              z coordinates
CPS          ELSE
CPS             compute the distance between the node and the given x,y 
CPS              coordinates
CPS          ENDIF
CPS          
CPS          IF the distance between the current node and the 
CPS           coordinates 
CPS           is less than the minimum distance and the zone is defined
CPS             set minimum distance to current distance
CPS             set node to the current node
CPS          ENDIF
CPS          
CPS      ENDFOR    
CPS          
CPS   END IF
CPS   
CPS END  near3
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      real*8 dismin,dist,x,y,z
      integer i,iflg,node

      dismin = 1.0d+10

      if (iflg .eq. 0) then
c     Changed to neq_primary (used to be neq) BAR - 12-15-99
         do i = 1, neq_primary
            if (icnl .eq. 0)  then
               dist = sqrt((cord(i,1) - x)**2
     *                  +  (cord(i,2) - y)**2
     *                  +  (cord(i,3) - z)**2)
            else
               dist = sqrt((cord(i,1) - x)**2 + (cord(i,2) - y)**2)
            end if
            if (dist .lt. dismin)  then
               dismin = dist
               node   =  i
            end if
         end do
      else
c     Changed to neq_primary (used to be neq) BAR - 12-15-99
         do i = 1, neq_primary
            if (icnl .eq. 0) then
               dist = sqrt((cord(i,1) - x)**2
     *                  +  (cord(i,2) - y)**2
     *                  +  (cord(i,3) - z)**2)
            else
               dist =  sqrt((cord(i,1) - x)**2 + (cord(i,2) - y)**2)
            end if
            if (dist .lt. dismin .and. izonef(i) .ne. 0)  then
               dismin =  dist
               node   =  i
            end if
         end do
      endif

      end
