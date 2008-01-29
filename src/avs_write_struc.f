      subroutine avs_write_struc
     1   (nelm,ns,icnl,nei,x,y,z,nnode,lu,ioformat,ierr)
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
CD1 Output AVS mesh connectivity information with FEHM
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
CD2 $Log:   /pvcs.config/fehm90/src/avs_write_struc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Wed Feb 14 10:47:18 1996   zvd
CD2 Modified requirements
CD2 
CD2    Rev 1.2   Mon Jan 29 13:15:20 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.1   Mon Jan 29 10:27:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.0   08/23/94 15:34:22   llt
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

      use comai, only : altc
      implicit none
      integer ns, icnl, nei, nnode, lu, ioformat, ierr, i, ipropelm, j
      integer nelm(ns*nei)
      real*8  x(nnode), y(nnode), z(nnode)
      character*5 char_type

      ipropelm = 1

      if (ioformat .eq. 2 )then
         do i = 1,nei
            call elem_type(nelm(((i-1)*ns)+1),ns,icnl,char_type)
               write(lu,100)i,ipropelm,char_type,(nelm((i-1)*ns + j), 
     .              j=1,ns)
         enddo
 100     format(i10.10,2x,i8,1x,a5,1x,8i8)
 110     format(8i8)
      elseif(ioformat .eq. 1)then
         write(ierr,*)'AVS_WRITE_STRUC'
         write(ierr,*)'Unformatted IO not implimented'
         write(ierr,*)'No action'
         return
      else
         write(ierr,*)'AVS_WRITE_STRUC'
         write(ierr,*)'Invalid output format'
         write(ierr,*)'No action'
         return
      endif

      return
      end

