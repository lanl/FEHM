      subroutine eullag3(itp,xxi,yyi,zzi,xxj,yyj,zzj,xix,xiy,xiz,
     *     xoa,xob,xoc)
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
CD1 this routine converts between two geometric reference frames
CD1 given a point i(xxi,yyi,zzi) and a nearest neighbor j(xxj,yyj,zzj)
CD1 along with a vector (xix,xiy,xiz) do a rotation of this
CD1 vector form x,y,z space to a,b,c space.
CD1      itp = 1 ==> x,y,z to a,b,c
CD1      itp= 2 ==> a,b,c to x,y,z
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
CD2 $Log:   /pvcs.config/fehm90/src/eullag3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   Wed Jun 26 12:37:50 1996   hend
CD2 Updated Prolog
CD2 
CD2    Rev 1.0   Thu Apr 25 13:34:24 1996   hend
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
CD4 This code is part of X3D
CD4
C**********************************************************************

      implicit none

      integer itp
      real*8  xxi,yyi,zzi,xxj,yyj,zzj,xix,xiy,xiz,xoa,xob,xoc
      real*8  sinph,cosph,sinth,costh
      real*8  a11,a12,a13,a21,a22,a23,a31,a32,a33
      real*8  ai11,ai12,ai13,ai21,ai22,ai23,ai31,ai32,ai33
      character*132 logmess
      
      call angle3(xxi,yyi,zzi,xxj,yyj,zzj,sinph,cosph,sinth,costh)
      goto (100,200) itp
100   continue
      a11=cosph
      a12=sinph
      a13=0.0
      a21=-costh*sinph
      a22=costh*cosph
      a23=sinth
      a31=sinth*sinph
      a32=-sinth*cosph
      a33=costh
      xoa=a11*xix+a12*xiy+a13*xiz
      xob=a21*xix+a22*xiy+a23*xiz
      xoc=a31*xix+a32*xiy+a33*xiz
      goto 9999
200   continue
      ai11=cosph
      ai12=-sinph*costh
      ai13=sinph*sinth
      ai21=sinph
      ai22=cosph*costh
      ai23=-cosph*sinth
      ai31=0.0
      ai32=sinth
      ai33=costh
      xoa=ai11*xix+ai12*xiy+ai13*xiz
      xob=ai21*xix+ai22*xiy+ai23*xiz
      xoc=ai31*xix+ai32*xiy+ai33*xiz
      goto 9999
9999  continue
      return
      end
