      function ran_sp( k )
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
CD1 To generate a pseudorandom number with a uniform distribution
CD1 between 0 and 1.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 Oct. 2, 1991 B. Robinson    00023    Initial Implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/ran_sp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Thu Feb 01 15:24:44 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   02/02/95 15:22:42   llt
CD2 added pvcs log info
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier            Type  Use  Description
CD3
CD3 k                      I    I/O  Random number seed used to
CD3                                  generate random number - a new
CD3                                  seed is computed in the process.
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 None
CD3
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 This function computes a random number by first performing a
CD6 remainder calculation using the mod function.  Then, it performs a
CD6 floating point calculation dividing by the real value of the
CD6 argument used in the mod calculation, thereby resulting in a value
CD6 between 0 an 1 dominated by roundoff error.
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 It is assumed that the integer k passed to the function is a
CD7 6-digit integer.
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 This is a general purpose utility routine used in the code 
CD8 wherever necessary.
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 Not Applicable.  See Special Comments.
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA FRACNET SRS
CDA
CDA Carnahan, B., H. A. Luther, and J. O. Wilkes, Applied Numerical
CDA Methods, John Wiley and Sons, Inc., NY (1969).
CDA
C**********************************************************************
CPS PSEUDOCODE
CPS 
CPS BEGIN ran_sp
CPS 
CPS Compute new value for k using remainder function mod
CPS Generate new random number
CPS
CPS END ran_sp
C**********************************************************************

      implicit none

      integer k
      real ran_sp

      k=mod(k*125,2796203)
      ran_sp=real(k)/2796203.0

      return
      end

