        subroutine ludcmp0(A,n,np,indx,d)
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
CD1  Given an n by n matrix A, with physical dimension np, this routine
CD1  replaces it by the LU decomposition of a rowwise permutation of
CD1  itself.  
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/ludcmp0.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:18   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Fri Apr 26 15:45:34 1996   gaz
CD2 returns different singular matrix info
CD2 
CD2    Rev 1.4   Fri Feb 16 10:04:30 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.3   Thu Jan 18 09:54:00 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.2   Wed Jan 10 14:33:52 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:57:48   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:30   pvcs
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
CD4 Given an n by n matrix A, with physical dimension np, this routine
CD4 replaces it by the LU decompostion of a row wise permutation of
CD4 itself.  A and n are input.  A is output, arranged as in equation
CD4 2.3.14 above; indx is an output vector which records the row
CD4 permutation effected by the partial pivoting; d is output as +/-1
CD4 depending on whether the number of row interchanges was even or odd,
CD4 respectively.  This routine is used in combination with lubksb0 to
CD4 solve linear equations or invert a matrix.
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5
CD5   A               REAL*8   I/O  LU decomposition of A
CD5   d               REAL*8   O    Flag denoting even or odd number of
CD5                                   row interchanges
CD5   indx            INT      O    Permutation vector 
CD5   n               INT      I    Number of rows/columns in A
CD5   np              INT      I    Dimension of A
CD5
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
CD7   aamax           REAL*8   Maximum value of A in a given row
CD7   dum             REAL*8   Figure of merit for pivot, dummy storage
CD7                              for element intercahnge values
CD7   i               INT      Loop index
CD7   imax            INT      Index of pivot row
CD7   j               INT      Loop index
CD7   k               INT      Loop index
CD7   sum             REAL*8   Modified element value
CD7   vv              REAL*8   Implicit scaling of row
CD7
CD7 Local Subprograms
CD7
CD7   Identifier      Type     Description
CD7
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN ludcmp0
CPS
CPS   initialize row interchanges flag to even
CPS   
CPS   FOR each row
CPS       set max row value of A to zero
CPS       FOR each element in the row
CPS           IF the absolute value of the current element is greater 
CPS            than the current maximum value
CPS              set the maximum value to the absolute value of the 
CPS               current element
CPS           END IF
CPS       END FOR
CPS
CPS       IF the row maximum value is zero
CPS          pause, the matrix is singular
CPS       END IF
CPS       
CPS       set scaling to 1 over maximum row value
CPS   END FOR
CPS       
CPS   FOR each column     
CPS   
CPS       FOR each row down to but excluding diagonal element
CPS           set sum to current element value
CPS           FOR all previous elements in column
CPS               compute updated sum
CPS           END FOR
CPS           set element value to current sum
CPS       END FOR
CPS       
CPS       set max row value of A to zero for pivot element search
CPS       set pivot row to current diagonal
CPS       
CPS       FOR each row from diagonal element to end
CPS           set sum to current element value
CPS           FOR all previous elements in column up to diagonal
CPS               compute updated sum
CPS           END FOR
CPS           set element value to current sum
CPS           
CPS           compute figure of merit for the pivot
CPS           IF this is the best value for pivot
CPS              set the pivot row value to current row number
CPS              set the maximum row value to figure of merit
CPS           END IF
CPS       END FOR
CPS       
CPS       IF rows should be interchanged
CPS          FOR each row element 
CPS              interchange the elements in the two rows
CPS          END FOR
CPS          toggle row interchanges flag
CPS          interchange scale factor
CPS       END IF
CPS       
CPS       set the permutation index
CPS       
CPS       IF the diagonal element is now zero 
CPS          set it to tiny
CPS       END IF
CPS       
CPS       IF this isn't the last column
CPS          compute the pivot scaling factor
CPS          FOR each row below th diagonal
CPS              compute the new element value
CPS          END FOR
CPS       END IF
CPS       
CPS   END FOR
CPS       
CPS END ludcmp0
CPS
C***********************************************************************

      implicit none

      integer n,i,j,k,imax,indx(n),nmax,np

*     Largest expected n, and a small number.
      parameter (nmax = 100)
*     vv stores the implicit scaling of each row.
      real*8 A(np, np), aamax, d, dum, sum, tiny, vv(nmax) 
      parameter (tiny = 1.e-20)
      
c     No row interchanges yet.
      d = 1.
c     Loop over rows to get the implicit scaling information.
      do i = 1, n
         aamax = 0.
         do j = 1, n
	    if (abs(A(i,j)) .gt. aamax)  aamax = abs(A(i,j))
         enddo
c     No nonzero largest element.
	  if (aamax .eq. 0.)  then                       
             indx(1)=-1
             return     
          endif
c     Save the scaling.
	  vv(i) = 1./aamax
       enddo

c     This is the loop over columns of Crout's method.
       do j = 1, n
c     This is equation 2.3.12 except for i=j.
	  do i = 1, j-1
             sum = A(i,j)
             do k = 1, i-1
                sum = sum - A(i,k)*A(k,j)
             enddo
             A(i,j) = sum
          enddo
c     Initialize for the search for largest pivot element.
	  aamax = 0.
          imax = j
c     This is i=j of equation 2.3.12 and i = j+1...n of equation 2.3.13.
	  do i = j, n
             sum = A(i,j)
             do k = 1, j-1
                sum = sum - A(i,k)*A(k,j)
             enddo
             A(i,j) = sum
c     Figure of merit for the pivot.
             dum = vv(i)*abs(sum)
c     Is it better than the best so far?
             if (dum .ge. aamax) then
                imax = i
                aamax = dum
             end if
          enddo
c     Do we need to interchange rows?
	  if (j .ne. imax) then
c     Yes, do so...
             do k = 1, n
                dum = A(imax,k)
                A(imax,k) = A(j,k)
                A(j,k) = dum
             enddo
c     ...and change the parity of d.
             d = -d
c     Also interchange the scale factor.
             vv(imax) = vv(j)
	  end if
	  indx(j) = imax
	  if (A(j,j) .eq. 0.)  A(j,j) = tiny
c     Now, finally, divide by the pivot element.
	  if (j .ne. n) then
c     If the pivot element is zero the matrix is singular (at least to
c     precision of the algorithm).  For some applications on singular
c     matrics, it is desirable to substitute tiny for zero.
             dum = 1./A(j,j)
             do i = j+1, n
                A(i,j) = A(i,j)*dum
             enddo
	  endif
c     Go back for the next column in the reduction.
       enddo

       return
       end
      
