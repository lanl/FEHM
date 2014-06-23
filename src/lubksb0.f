	subroutine lubksb0(A, n, np, indx, b)
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
CD1  PURPOSE
CD1
CD1  To solve the set of n linear equations A*x = b.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/lubksb0.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Feb 16 10:04:28 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.3   Thu Jan 18 09:53:50 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Wed Jan 10 14:17:34 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:57:46   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:28   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable. This is a general utility routine used in the code
CD3  as needed for solving a matrix.
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  Taken from Numerical Recipes, section 2.3.
CD4  
CD4  Solves the set of n linear equations A*x=b.  Here A is input, not
CD4  as the matrix A but rather as its LU decomposition, determined by 
CD4  routine ludcmp. indx is input as the permutation vector returned by
CD4  ludcmp. b is input as the right-hnd side vector b, and returns with
CD4  the solution vector x.  A, n, np, and indx are not modified by this
CD4  routine and can be left in place for successive calls with differnt
CD4  right-hand sides b. This routine takes into account the possibility
CD4  that b will begin with many zero elements, so it is efficient for 
CD4  use in matrix inversion.
CD4  
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5   
CD5   A               REAL*8   I    LU decomposition of A
CD5   b               REAL*8   I/O  Vector, input right-hand side, output
CD5                                   solution vector
CD5   indx            INT      I    Permutation vector returned by ludcmp0
CD5   n               INT      I    Number of elements in b
CD5   np              INT      I    Number of elements in A
CD5
CD5 Interface Tables
CD5
CD5   None 
CD5
CD5 Files
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6   None
CD6
CD6 Global Subprograms
CD6
CD6   None
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   i               INT      Loop index
CD7   ii              INT      Index of first nonvanishing element of b
CD7   j               INT      Loop index
CD7   ll              INT      Element index
CD7   sum             REAL*8   Modified element value
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN lubksb0
CPS 
CPS   set the index for the first nonzero element of b to 0
CPS   FOR each element of b [do forward substitution]
CPS       set the index, initial sum, initial value of b
CPS       IF the nonzero element index is set
CPS          FOR each element from the first nonzero b to the previous b
CPS              compute the new sum
CPS          END FOR
CPS       ELSE IF the sum (element value) is not zero
CPS          set the index for the first nonzero element of b
CPS       END IF
CPS       set value of b for this element to sum
CPS   END FOR
CPS   
CPS   FOR each element of b [do backward substitution]
CPS       set the initial sum
CPS       IF this isn't the last element
CPS          FOR each element from the previous element to the last
CPS              compute the new sum
CPS          END FOR
CPS       END IF
CPS       set value of b for this element to sum over diagonal A of 
CPS        that element
CPS   END FOR
CPS   
CPS END lubksb0
CPS
C***********************************************************************

 	implicit none

	integer np,n,i,ii,indx(n),j,ll
        real*8 A(np,np),b(n),sum

c  When ii is set to a positive value, it will become the index of the
c  first nonvanishing element of b.  We now do the forward substitution,
c  equation 2.3.6.  The only new wrinkle is to unscramble the
c  permutation as we go.
	ii = 0
	do i = 1, n
	   ll = indx(i)
	   sum = b(ll)
	   b(ll) = b(i)
	   if (ii .ne. 0) then
	      do j = ii, i-1
		 sum = sum - A(i,j)*b(j)
	      enddo
	   else if (sum .ne. 0) then
c  A nonzero element was encountered, so from now on we will have to do
c  the sums in the loop above.
	      ii = i
	   end if
	   b(i) = sum
	enddo

c  Now we do the backsubstitution, equation 2.3.7.
	do i = n, 1, -1
	   sum = b(i)
	   if (i .lt. n) then
	      do j = i+1, n
		 sum = sum - A(i,j)*b(j)
	      enddo
	   end if
c  Store a component of the solution vector x.
	   b(i) = sum/A(i,i)
	enddo
	
	return
	end
