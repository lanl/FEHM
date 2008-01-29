      subroutine setbit(nbitnum,ibitnum,iword,istate)
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
CD1  This subroutine sets the specified bit to the input state.  The 
CD1  bits are arranged from least-significant-bit (0 which is on the
CD1  right) to the most-significant-bit (nbitnum which is on the left).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/setbit.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   Thu Jan 11 10:42:16 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.0   04/27/95 17:53:08   llt
CD2 original routine
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
      implicit none

      integer nbitnum,ibitnum,istate,iword(*)
      integer itotal, ichange, j, jchange, nwd1, jstate, idum, imask

c     Note: for the Digital FORTRAN version comment out the 
c     external line below

C     external integer ishft

C     DEFINE A STATEMENT FUNCTION FOR CREATING A MASK BY SHIFTING
C        BY A SPECIFIED NUMBER OF BITS.
      imask(idum)=2**(idum-1)

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

C     ******************************************************************
C     SET THE BIT TO THE INPUT STATE BY DOING A TRUTH TABLE WITH THE
C        CURRENT BIT SETTING AND THE STATE VALUE.
C
C        CURRENT BIT  STATE VALUE  OUTPUT BIT
C        ***********  ***********  **********
C                  0            1           1
C                  1            0           0
C                  1            1           1
C                  0            0           0
C
      if(nwd1.eq.0.and.istate.ne.0) then
         jstate=ishft(1,(j-1))
         iword(jchange)=or(iword(jchange),jstate)
      elseif(nwd1.eq.1.and.istate.eq.0) then
         jstate=not(ishft(1,(j-1)))
         iword(jchange)=and(iword(jchange),jstate)
      elseif(nwd1.eq.1.and.istate.ne.0) then
      elseif(nwd1.eq.0.and.istate.eq.0) then
      endif

      goto 9999
 9999 continue

      return
      end
