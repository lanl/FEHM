      subroutine humidity(hum,alp,beta,sr,smax,t,s)
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
CD1  This subroutine calculates a saturation for a given humidity value.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/humidity.f_a  $ 
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Feb 15 11:07:02 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.2   Wed Jan 10 13:00:14 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   04/12/95 10:57:36   gaz
CD2 changed t to t+273.16(corrected error)
CD2 
CD2    Rev 1.0   03/10/95 11:58:22   llt
CD2 new routine
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
CD4  This is a general utility routine used as needed.
CD4  
C***********************************************************************

      implicit none

      real*8 hum,h,alp,beta,alamda,sr,smax,t,s,dgls
      parameter (dgls=0.461)

      h=-log(hum)*dgls*(t+273.16)/9.8e-3
      alamda=1.0-1.0/beta
      s=1.0/((alp*h)**(1.0/(1.0-alamda))+1.0)**alamda
      s=s*(smax-sr)+sr

      return
      end
