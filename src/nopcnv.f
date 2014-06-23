      subroutine nopcnv  ( b,
     2                   igaus,
     3                   neq,
     4                   ncon,
     5                   nop,
     6                   npvc,
     7                   npvt,
     8                   nopt,
     9                   dum,
     t                   nbd )
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
CD1 To perform the symbolic factorization in natural order.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-12-92     G. Zyvoloski   00097   Initial Implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/nopcnv.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Sep 12 08:25:16 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.4   Thu Feb 01 14:01:44 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   11/15/95 15:14:36   gaz
CD2 faster symbolic factorization
CD2 
CD2    Rev 1.1   03/18/94 16:07:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:56   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier      Type  Use  Description
CD3
CD3 b               int    I   Scratch space used to perform symbolic
CD3                            factorization
CD3 igaus           int    I   Factorization level for each node
CD3 neq             int    I   Number of equations
CD3 ncon            int    I   Connectivity matrix for the equation set
CD3 nop             int    O   Connectivity matrix for factorization
CD3                            matrix
CD3 npvc            int    I   Array of positions in ncon array where
CD3                            ncon(ncon(npvc(i))) = i
CD3 npvt            int    O   Array of pivot positions in nop matrix
CD3 nopt            int    I   Storage space needed for operations
CD3 dum             int    I   Storage space needed for operations
CD3 nbd             int    I   Storage space alloted for nop
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
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 None
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4 None
CD4
CD4 Global Subprograms
CD4
CD4 None
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier  Type        Description
CD5
CD5 i           int         Array index
CD5 i1          int         Do loop index
CD5 i2          int         Do loop index
CD5 i3          int         Do loop index
CD5 i4          int         Do loop index
CD5 icount      int         Index parameter used in computing nop
CD5 idd         int         Index parameter used in computing nop
CD5 igmax       int         Maximum value of igaus array
CD5 k           int         Do loop index
CD5 kb          int         Index parameter used in symbolic factor-...
CD5...                      ization loops
CD5 kc          int         Index parameter used in symbolic factor-...
CD5...                      ization loops
CD5 kbp         int         Index parameter used in symbolic factor-...
CD5...                      ization loops
CD5 kk          int         Do loop index
CD5 kp          int         Do loop index
CD5 minkb       int         Parameter used in symbolic factorization
CD5 nalot       int         Maximum index in b array
CD5 neqp1       int         neq + 1
CD5
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 First, the code loops through all values of igaus to determine the
CD6 maximum value.  For the case of some values of igaus greater than
CD6 1, the code first initializes all value of b.  Next, the first row
CD6 is sybolically factored separately, after which each remaining row
CD6 is symbolically factored within a recursive loop.  In this loop,
CD6 the code first initializes a working array, symbolically factors
CD6 row kp, adds the combinations of row kbp, and then calculates nop.
CD6 Finally, nop is changed as needed for variable factorization before
CD6 returning.
CD6
CD6 If there is no value greater than 1 (i.e., no fill),
CD6 the code loops through all values of ncon and sets the array nop
CD6 equal to ncon.  It then sets every diagonal pointer in nop to that
CD6 in ncon before returning.
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 N/A
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8 N/A
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 3.1.2. Perform Symbolic Factorization
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA GZSOLVE SRS.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN nopcnv
CPS
CPS Initialize maximum igaus value to 0
CPS
CPS FOR each value of igaus
CPS   Update largest value of igaus
CPS ENDFOR
CPS
CPS IF there is some fill (igmax > 1)
CPS
CPS   Compute maximum address in b array
CPS
CPS   FOR each value in b array
CPS     Initialize value in the b array
CPS   ENDFOR
CPS
CPS   Initialize index parameters
CPS
CPS   FOR each unknown connected to node 1
CPS
CPS     Set value of b
CPS     Set value of nop
CPS
CPS     IF this entry in ncon is 1
CPS       Set value of npvt
CPS     ENDIF
CPS
CPS   ENDFOR
CPS
CPS   Set first and second values of nop
CPS
CPS   FOR all remaining rows
CPS
CPS     Determine index values
CPS
CPS     FOR all values of working array
CPS       Initialize value in working array
CPS     ENDFOR
CPS
CPS     FOR each unknown connected to the current node
CPS       Set value in array dum to 1
CPS       Update lowest number of unknown in this section of ncon array
CPS     ENDFOR
CPS
CPS     FOR all nonzero values up to the one previous to the diagonal
CPS       FOR all nonzero values from the diagonal to the last entry
CPS         Update value in dum array
CPS       ENDFOR
CPS     ENDFOR
CPS
CPS     FOR each equation
CPS
CPS       IF fill in level is less than maximum fill in level
CPS         Set value of nop
CPS         IF column row index exists
CPS           Set pivot index
CPS         ENDIF
CPS         Set the value of b
CPS         EXITIF more than neq elements have been set
CPS       ENDIF
CPS
CPS     ENDFOR
CPS
CPS     Set the value of nop for this row
CPS     
CPS   ENDFOR
CPS
CPS   FOR all rows except the first
CPS
CPS     FOR all nonzero elements pointer array for lu factorization...
CPS     ... matrix
CPS       IF the neighbor is above the diagonal
CPS         Set current value of ig
CPS       ENDIF
CPS       IF the next value of nopt needs to be set
CPS         Set next value of nopt
CPS       ENDIF
CPS     ENDFOR
CPS
CPS     FOR all nonzero elements in this row
CPS       Set value of nop
CPS       IF this is the pivot element
CPS         Set value of npvt
CPS       ENDIF
CPS     ENDFOR
CPS
CPS   ENDFOR
CPS
CPS   Set final value of nop
CPS
CPS ELSE there is no fill
CPS
CPS   FOR all values in ncon array
CPS     Set nop equal to ncon
CPS   ENDFOR
CPS
CPS   FOR each equation
CPS     Set npvt equal to npvc
CPS   ENDFOR
CPS ENDIF
CPS
CPS END nopcnv
C**********************************************************************
c
c symbolic factorization in natural  order
c this routine is called from slvesu
c
c b- integer array defined as ib in slvesu
c dum(neq)- integer array defined as idum in slvesu
c nbd-size of array nop(nbd=nnop in slvesu)

      implicit none
      integer b(*)
      integer dum(*)
      integer i
      integer i1
      integer i2
      integer i3
      integer i4
      integer icount
      integer idd
      integer igmax
      integer igmin
      integer k
      integer kb
      integer kc
      integer kbp
      integer kk
      integer kp
      integer minkb
      integer maxkb
      integer nalot
      integer nbd
      integer neq
      integer neqp1
      integer nop(*)
      integer npvc(*)
      integer npvt(*)
      integer ncon(*)
      integer igaus(*)
      integer nopt(*)
      integer itparr
      integer icountc

      neqp1=neq+1
c
c check variable factorization to find max igauss(igmax)
c
      igmax=0
      igmin=0

      do i = 1, neq
         igmax = max ( igmax, igaus(i) )
         igmin = min ( igmin, igaus(i) )
      end do

      if( igmax .gt. 1 ) then
c
c
c initialize matrix b
c
         nalot = ncon(neqp1)    
         do i = 1, nalot
          nop(i) = ncon(i)
         end do
         do i=1,neq
          npvt(i)=npvc(i)
         enddo
c
c do first row separately
c
         i = 1
         i1 = ncon(i) + 1
         i2 = ncon(i+1)
         icount = neqp1 

         do k = i1, i2
            kb = ncon(k)
            icount = icount + 1
            nop(icount) = kb
            if( kb .eq. 1 ) then
               npvt(1) = icount
            end if
         end do

         nop(1) = neqp1
         nop(2) = icount 
         icountc=icount
c
c start parallel iterations              
c
         do itparr = 1,igmax-1
c
c set old nop(b) = new nop
c
         do i=1,nop(neqp1)
          b(i)=nop(i)
         enddo
         do i=1,neq
          nopt(i)=npvt(i)
         enddo
c
         icount = icountc
c
c symbolically factor remaining rows
c 
         do i = 1, neq
           dum(i) = 0
         end do
         do kp = 2, neq
          if(igaus(kp).ge.itparr+1) then
            i1 = b(kp) + 1
            i2 = b(kp+1)

c initialize working vector

            minkb=neq
            maxkb=1
            do i=i1,i2
             dum(b(i))=1
               minkb=min(minkb,b(i))
               maxkb=max(maxkb,b(i))
            enddo

c symbolically factor row kp

c add combination of row kpp
            do kbp = i1,nopt(kp)-1    
               idd = b(kbp)
               i3 = nopt(idd)
               i4 = b(idd+1)
               minkb=min(minkb,idd)
               do kk = i3, i4
                  kc = b(kk)
                  maxkb=max(maxkb,kc)
                  dum(kc) = 1                                     
               end do

            end do
c
c calculate nop
c

            do i = minkb,maxkb
               idd = dum(i)
               dum(i) = 0
               if(idd.ne. 0 ) then
                  icount=icount+1
                  nop(icount) = i
                  if( i .eq. kp ) then
                     npvt(kp)=icount
                  end if

               end if
            end do
            nop(kp+1)=icount
          endif
            if(icount.gt.nbd) then
             nbd = -nbd
             return
            endif
         end do
         end do
c
c set nop=ncon for the case with no fill
c
      else
         i2 = ncon(neq+1)

         do i = 1, i2
            nop(i)=ncon(i)
         end do

c
c set diagonal pointer in nop to that in ncon
c
         do i = 1, neq
            npvt(i) = npvc(i)
         end do

c
c
      end if
      return
      end
