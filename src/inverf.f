      function inverf(x)
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
CD1  This subroutine calculates the inverse of the error function
CD1  for a value of x between 0 and 1.  It does so using linear 
CD1  interpolation between chosen values of the erf(x) curve.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/inverf.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Wed Jan 10 13:11:52 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/16/95 09:48:46   llt
CD2 added PVCS log history
CD2
CD2    Rev 1.0   03/16/95 09:00:44   robinson
CD2 Initial revision.
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.5  Cell-based particle-tracking module
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

      implicit none

      integer i
      real inverf,ef(8),a(8),b(8), x
      data ef/.2227,.42839,.60386,.79691,.91031,.97635,.99532,1./
      data a/.8980691,.972337,-1.255877,-1.033911,-0.8452046,
     +     -0.6909501,-0.5685161,-0.4299497/
      data b/0.,-1.653945e-2,9.494722e-2,0.1842109,0.3148551,0.476399,
     +     0.6754972,0.9983227/

      if(x.lt.0.) x=-x
      if(x.gt.1.) then
         write(6,110)
 110     format(1x,'value of x greater than 1')
         return
      end if
      do i=1,2
         if(x.lt.ef(i)) then
            inverf=a(i)*x+b(i)
            return
         end if
      enddo
      do i=3,7
         if(x.lt.ef(i)) goto 11
      enddo
      i=8
 11   inverf=a(i)*alog10(1.-x)+b(i)

      return
      end
