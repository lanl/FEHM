      subroutine angle3(xxi,yyi,zzi,xxj,yyj,zzj,sinph,cosph,sinth,costh)
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
CD1 this routine calculates the rotation angles needed to rotate
CD1 the x-axis from i to j and the from here
CD1 to rotate the z-axis from i to j.
CD1
CD1 determine the angles of rotation needed to point the
CD1 z-axis from mass point 'i' to neighbor 'j'.
CD1
CD1 find the cosine and sine of the angle with reference to the z-axis
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
CD2 $Log:   /pvcs.config/fehm90/src/angle3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.0   Thu Apr 25 13:30:18 1996   hend
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
CD4 This code is part of X3D
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
C**********************************************************************

      implicit none

      real*8 xxi, yyi, zzi, xxj, yyj, zzj
      real*8 sinph, cosph, sinth, costh
      real*8 xnoise, ds, dsxy, dsxysq, dsiv
      character*132 logmess
      
      xnoise=1.0e-10
      dsxysq=(xxj-xxi)**2 + (yyj-yyi)**2
      dsxy=sqrt(dsxysq)
      ds=sqrt(dsxysq+(zzj-zzi)**2)
      dsiv=1.0/ds
      costh=(zzj-zzi)*dsiv
      sinth=dsxy*dsiv
c
c find the cosine and sine of the angle with reference to the x-axis
c
      if(dsxy.lt.ds*xnoise) then
         cosph=0.0
         sinph=1.0
      else
         cosph=-(yyj-yyi)/dsxy
         sinph= (xxj-xxi)/dsxy
      endif
      goto 9999
 9999 continue
      return
      end
