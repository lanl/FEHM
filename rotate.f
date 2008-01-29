      subroutine rotate (xxi,yyi,zzi,xxj,yyj,zzj,xix,xiy,xiz,xoa,xob,
     +     xoc)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Driver for Harold's (X3D) rotation matrix code
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 25-APR-96    S. Henderson   22      Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rotate.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:14:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   Wed Jun 26 12:38:48 1996   hend
CD2 Updated Prolog
CD2 
CD2    Rev 1.0   Thu Apr 25 13:36:04 1996   hend
CD2 Initial coding for use in long/trans dispersion
CD2 
C**********************************************************************
CD3
CD3 REQUIREMENTS TRACEABILITY
CD3
CD3 2.3.4 Solute-transport equations
CD3
C**********************************************************************
CD4
CD4 SPECIAL COMMENTS
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      implicit none

      integer itp
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3
      real*8 xxi,yyi,zzi,xxj,yyj,zzj,xix,xiy,xiz,xoa,xob,xoc

c first two define new z axis, last is point to convert
      x1=0.0
      y1=0.0
      z1=0.0
      x2=xxj
      y2=yyj
      z2=zzj
      x3=xix-xxi
      y3=xiy-yyi
      z3=xiz-zzi
      itp=1
      call eullag3(itp,x1,y1,z1,x2,y2,z2,x3,y3,z3,
     *             xoa,xob,xoc)

      return
      end


