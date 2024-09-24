      subroutine elem_type(nelm,nssave,icnl,char_type)
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
CD1 To determine the element type.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/elem_type.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Wed Jan 17 12:05:02 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.3   Wed Jan 10 11:13:40 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.2   Wed Jan 10 09:19:32 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:43:04   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:42   pvcs
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
CD4  This is a general utility routine used in the code whenever 
CD4  necessary.
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5 
CD5 Identifier    Type   Use    Definition
CD5
CD5 nelm          INT     I     Array of nodes comprising the element
CD5 nssave        INT     I     Number of nodes in the element
CD5 icnl          INT     I     Flag to denote dimensionality of problem
CD5 char_type     CHAR    O     Used to denote element type
CD5
CD5 Interface Tables
CD5
CD5 None
CD5
CD5 Files
CD5 
CD5 None
CD5
C***********************************************************************
CD6
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
CD7   ns              INT      Number of nodes in the element
CD7   done            LOGICAL  Flag deonting if element has been 
CD7                            successfully determined
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN elem_type
CPS
CPS LOOP to determine the element type
CPS
CPS   Set flag to denote element has been successfully determined
CPS
CPS   IF the element has 8 nodes
CPS     Set flag to denote hex
CPS   ELSEIF the element has 6 nodes
CPS     Set flag to denote hex
CPS   ELSEIF the element has 5 nodes
CPS     Set flag to denote prism
CPS   ELSEIF the element has 4 nodes
CPS     IF this is a 3-D problem
CPS       Set flag to denote tet
CPS     ELSE this is a 2-D problem
CPS       Set flag to denote quad
CPS     ENDIF
CPS   ELSEIF the element has 3 nodes
CPS     Set flag to denote triangle
CPS   ELSEIF the element has 2 nodes
CPS     Set flag to denote line
CPS   ELSEIF the element has 1 node
CPS     Set flag to denote point
CPS   ENDIF
CPS
CPS   IF the number of nodes in element needs to be reduced
CPS     IF there are still at least 2 nodes in the element
CPS       Reduce number of elements by 1
CPS       Set flag to continue determination of element type
CPS     ENDIF
CPS   ENDIF
CPS
CPS EXITIF element has been successfully determined
CPS
CPS ENDLOOP to determine the element type
CPS 
CPS END elem_type
CPS
C***********************************************************************
c***date written   930910   (yymmdd)
c***keywords element type
c***author   Carl W. Gable   (Los Alamos National Laboratory)
c***purpose  Determine element type.
c***description At present routine cannot tell the
c               difference between quad and tet. Default
c               is tet for 4 node elements when icnl = 0,
c                 quad for 4 node elements when icnl .ne. 0.

      implicit none

      character*5 char_type
      integer nssave,icnl,nelm(nssave),ns
      logical done

      ns = nssave
c
c     if(icnl .eq. 0) => three dimensional problem
c     else            => two   dimensional problem
c

 1000 continue

      done = .TRUE.
      
      if(ns .eq. 8)then
         char_type = 'hex  '
      elseif(ns .eq. 6)then
         char_type = 'prism'
      elseif(ns .eq. 5)then
         char_type = 'pyr  '
      elseif(ns .eq. 4)then
         
c     
c     check for quad or tet
c     
         if(icnl .eq. 0)then
            char_type = 'tet  '
         else
            char_type = 'quad '
         endif
      elseif(ns .eq. 3)then
         char_type = 'tri  '
      elseif(ns .eq. 2)then
         char_type = 'line '
      elseif(ns .eq. 1)then
         char_type = 'pt   '
      endif
c     
c     check to see if the size of nelm is larger than
c     the structure of the element. i.e. a prism written
c     into an 8 node structure. Condition that will flag
c     this is the end nodes are redundent or set to node 0.
c
      if((nelm(ns) .eq. 0) .or. 
     1     (nelm(ns) .eq. nelm(ns-1)))then
         if(ns .ge. 2)then
            ns = ns - 1
            done = .FALSE.
         endif
      endif
      if(done) goto 2000
      goto 1000
 2000 continue
      return
      end
      subroutine elem_geo_ctr(iflg)
c
c gaz 050822 new routine to save element data for geo file
c
      use comai, only : ns_in, nei_in, ns, nei, ivf, ifdm_elem
      use combi
      use avsio
      implicit none
c      
      integer iflg, i, elem_stor_size
c   
      if(iogeo.eq.0) return
      if(iflg.eq.0.and.ivf.ne.-1) then
c allocate storage          
       if(.not.allocated(elem_geo)) then
        if(ns_in.ne.0.and.nei_in.ne.0) then
         allocate(elem_geo(ns_in*nei_in))       
        endif
       endif
      else if(iflg.eq.0.and.ivf.eq.-1) then
c gaz 122622 added element generation for fdm
       ifdm_elem = 1
       call structured(4)
       if(.not.allocated(elem_geo)) then
        if(ns_in.ne.0.and.nei_in.ne.0) then
         allocate(elem_geo(ns_in*nei_in))       
        endif
       endif
      else if(iflg.eq.-1) then 
c  gaz 061522 wait until after vtk calcs to deallocate
c       deallocate(elem_geo)
      else if(iflg.eq.1.and.ivf.ne.-1) then 
       elem_stor_size = ns_in*nei_in
       do i = 1, elem_stor_size
         elem_geo(i) = nelm(i)  
       enddo
       nact_elem_geo = 1
       continue
      else if(iflg.eq.2) then
      else
      endif
      return
      end
