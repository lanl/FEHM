      subroutine avs_write_cord(x,y,z,nnode,lu,ioformat,icnl2,ierr2)
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
CD1
CD1 PURPOSE
CD1
CD1 Output AVS coordinate information for FEHM.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-SEP-93    Carl Gable     22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/avs_write_cord.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Thu Apr 17 12:01:58 1997   llt
CD2 minor format error corrected (gaz)
CD2 
CD2    Rev 1.6   Wed Apr 09 09:38:26 1997   gaz
CD2 changes made to accomodate head output
CD2 
CD2    Rev 1.5   Thu Oct 24 15:11:30 1996   zvd
CD2 Increased output precision
CD2
CD2    Rev 1.4   Wed Feb 14 10:47:16 1996   zvd
CD2 Modified requirements
CD2 
CD2    Rev 1.3   Mon Jan 29 13:10:58 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.2   Mon Jan 29 10:20:18 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   01/26/95 14:22:12   tam
CD2 added icnl2 to argument list
CD2 If 2-dim problem, write 0's for z coordinate
CD2 
CD2    Rev 1.0   08/23/94 15:33:58   llt
CD2 Original version
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
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
CD4 None
CD4
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
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
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
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
c
c----------------------------------------------------------------------------
c  phs 4/27/00   added comai  and   if(avsx) then add header to *_geo
c  phs 4/27/00   needed to rename dummys icnl + ierr  to icnl2 + ierr2
c----------------------------------------------------------------------------
      use comai
      implicit none

      integer i, nnode, lu, ioformat, icnl2, ierr2
      real*8 x(nnode), y(nnode), z(nnode)
      real*8 dummy

c
c     ASCII formatted output
c
c     if(icnl2 .eq. 0) => three dimensional problem
c     else            => two   dimensional problem
c

      if(ioformat .eq. 2)then

	 if(altc(1:4).EQ.'avsx')       
     x      write(lu,110) nnode, nei, 0, 0, 0  

            do i = 1,nnode
                write(lu,100)i,x(i),y(i),z(i)
            enddo
      else
c
c     ERROR
c
      write(ierr2, *)'AVS_WRITE_CORD'
      write(ierr2, *)'Invalid output format'
      write(ierr2, *)'No action'
      return
      endif
 
  100 format(i10.10,5x,3(e16.9,2x))
  110 format(i10.10,5x,4(i10.2,2x))
 
      return
      end

