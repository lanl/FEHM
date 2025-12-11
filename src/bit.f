      logical function bit(nbitnum,ibitnum,iword)
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
CD1 PURPOSE
CD1
CD1  This routine sets the specified bit to the input state. The bits
CD1  are arranged from least-significant-bit (0 which is on the right)
CD1  to the most-significant-bit (nbitnum which is on the left).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/bit.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Wed Jan 10 10:48:28 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.1   Tue Jan 09 14:29:00 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.0   04/27/95 17:53:20   llt
CD2 original routine
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable. See special Comments.
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  This is a general purpose utility routine, used throughout the
CD4  code where appropriate.
CD4
C***********************************************************************
C
C      INPUT ARGUMENTS -
C
C         nbitnum - TOTAL NUMBER OF BITS IN THE INPUT WORD. THIS
C                      CAN BE ANY LENGTH.
C         ibitnum - THE BIT THAT IS TO BE CHANGED.
C         iword() - THE STARTING LOCATION OF THE MOST-SIGNIFICANT-BIT.
C         istate  - THE STATE OF THE BIT TO BE SET.
C
C      OUTPUT ARGUMENTS -
C
C         iword() - THE STARTING LOCATION OF THE MOST-SIGNIFICANT-BIT.
C
C#######################################################################
C
      implicit none

      integer nbitnum,ibitnum,iword(*)
      integer itotal, ichange, j, jchange, nwd1, idum, imask

c     Note: when compiling with Digital FORTRAN, the external statement
c     must be commented out

C     external integer ishft
C
C#######################################################################
C
C     DEFINE A STATEMENT FUNCTION FOR CREATING A MASK BY SHIFTING
C        BY A SPECIFIED NUMBER OF BITS.
C
      imask(idum)=2**(idum-1)
C
C#######################################################################
C
C     *** CALCULATE THE TOTAL NUMBER OF 32-BIT ELEMENTS IN THE INPUT
C            VECTOR.
      itotal=1+(nbitnum-1)/32
C
C     *** CALCULATE THE INDEX NUMBER OF THE 32-BIT ELEMENT TO CHANGE.
      ichange=1+(ibitnum-1)/32
C
C     *** CALCULATE THE BIT NUMBER OF THE 32-BIT ELEMENT TO CHANGE.
      j=ibitnum-32*(ichange-1)+1
C
C     *** CALCULATE THE INDEX NUMBER OF THE INPUT VECTOR TO CHANGE.
      jchange=itotal-ichange+1
C
C     *** FIND THE CURRENT VALUE OF THE BIT IN THE INPUT VECTOR.
      nwd1=ishft(iand(iword(jchange),imask(j)),-(j-1))
C
C     ******************************************************************
C     TEST THE CURRENT BIT SETTING. 
C         =  0 ==> FALSE
C        < > 0 ==> TRUE
C
      if(nwd1.eq.0) then
         bit=.false.
      else
         bit=.true.
      endif
C
C     ******************************************************************
C
      goto 9999
 9999 continue
      return
      end
