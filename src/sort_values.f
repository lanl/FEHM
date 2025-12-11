      subroutine sort_values(n, nsize, indx, value)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To sort parameter values to facillitate finding of closest transfer
!D1 function curve.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 
!D2 Initial implementation: 15-DEC-02, Programmer: Bruce Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/sort_values.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2
!!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************
c
c     Indexing and Ranking algorithm for sorting, from Numerical Recipes
c     Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T.
c     Vetterling, 1986, Numerical Recipes. The Art of Scientific
c     Computing, Cambridge University Press, Cambridge, pp. 232-234.

      implicit none

      integer nsize
      integer n, j, l, ir, indxt, i
      integer indx(nsize)
      real(8) value(nsize), q

      do j = 1, n
         indx(j) = j
      end do
      l = n/2+1
      ir=n
 10   continue
      if(l.gt.1) then
         l=l-1
         indxt=indx(l)
         q=value(indxt)
      else
         indxt=indx(ir)
         q=value(indxt)
         indx(ir)=indx(1)
         ir=ir-1
         if(ir.eq.1) then
            indx(1)=indxt
            return
         end if
      end if
      i=l
      j=l+l
 20   if(j.le.ir) then
         if(j.lt.ir) then
            if(value(indx(j)).lt.value(indx(j+1))) j=j+1
         end if
         if(q.lt.value(indx(j))) then
            indx(i)=indx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         goto 20
      end if
      indx(i)=indxt
      goto 10
      end
