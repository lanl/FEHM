      subroutine  rarng(iflg)
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
CD1 To rearrange values in coordinate array for 2-d problems.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rarng.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jan 17 12:53:00 1996   zvd
CD2 Added prolog.
CD2 
CD2    Rev 1.4   04/25/95 08:20:46   llt
CD2 retrieved lost log history
CD2 
CD2    Rev 1.3   04/25/95 08:15:08   llt
CD2 added log history information to prolog
CD2
CD2    Rev 1.2   03/20/95 13:31:56   gaz   
CD2 zeroed out cord(i,3) for 2-d problems
CD2
CD2    Rev 1.1   03/18/94 15:43:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2
CD2    Rev 1.0   01/20/94 10:26:44   pvcs
CD2 original version in process of being certified
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable.  See Special Comments.
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4  This is a general utility routine used in the code.
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
CD6 icnl, neq, cord
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
CD7   i               int      do loop index
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN rarng
CPS
CPS IF this is an xz or xz radial simulation
CPS
CPS   FOR each node
CPS     Move z coordinate value into position for y coordinate
CPS     Set z coordinate to zero
CPS   ENDFOR each node
CPS
CPS ELSEIF this is a yz or yz radial simulation
CPS
CPS   FOR each node
CPS     Move y coordinate value into position for x coordinate
CPS     Move z coordinate value into position for y coordinate
CPS     Set z coordinate to zero
CPS   ENDFOR each node
CPS
CPS ENDIF
CPS
CPS END rarng
CPS
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      implicit none

      integer i,iflg,icode, isubst
      if(iflg.eq.0) then
         allocate(corz(neq,3))
c       call mmgetblk ("corz",  "combi", ipcorz,  n0*3,  2, icode)
          do i = 1, neq
             corz(i,1)   =  cord(i,1)
             corz(i,2)   =  cord(i,2)
             corz(i,3)   =  cord(i,3)
          end do
       if ( icnl .eq. 2 .or. icnl .eq. 5 )  then
          do i = 1, neq
             cord(i,2)   =  cord(i,3)
             cord(i,3)   =  0.0
          end do
 
       else if ( icnl .eq. 3 .or. icnl .eq. 6 )  then
          do i = 1, neq
             cord(i,1)   =  cord(i,2)
             cord(i,2)   =  cord(i,3)
             cord(i,3)   =  0.0
          end do
       end if
      else if(iflg.eq.1) then
c release dummy memory
         deallocate(corz)
c       call mmrelblk ("corz",  "combi", ipcorz, icode)
      endif


      return
      end
