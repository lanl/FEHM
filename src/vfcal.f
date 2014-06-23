          function vfcal(mi,rl,dvfp,drlp)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This function calculates the aperture space.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vfcal.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:06   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 16 13:23:40 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.2   Thu Jan 11 12:26:14 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 16:12:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:18   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.4.10 Stress-dependent properties
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      implicit none

      integer mi
      real*8 dvfp, drlp, rl, vfcal
      real*8 a, b, c, d0, pmul, fracw, phd, rat

      vfcal=1.0
      dvfp=0.
      rl=1.
      drlp=0.
      return
c     for aperture change erase return
 100  pmul=1000.
      fracw=.00100*pmul
      if(ps(mi).ne.1.0d0) return
      phd=29.7
      rat=(phi(mi)-10.0)/phd
      a=1.d-4
      b=log(1./a)
      c=a*a*a
      d0=3*b
      vfcal=a*exp(b*rat)*fracw
      dvfp=a*b*exp(b*rat)*fracw/phd
      rl=c*exp(d0*rat)
      drlp=c*d0*exp(d0*rat)/phd
      if(vfcal.gt.2.) rl=8.
      if(vfcal.gt.2.) drlp=0.
      if(ps(mi).gt..01) return
      vfcal=vf(mi)
      
      return
      end
