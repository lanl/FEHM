      subroutine sub_bksub31(neq,a,b,r,na,nrhs,ncon,nop
     *     ,npvt,piv,sto5,icoupl,tollr,overf)         
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
CD1 To perform backsubstitution for a 3-1 degree of freedom (rdof) problem.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 4-6-94       G. Zyvoloski   97      Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/sub_bksub31.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:22   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:04 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:21:10 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:18:48   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:36   pvcs
CD2 original version
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 neq          integer  I      Number of entries in the array for
CD3                                  each degree of freedom
CD3 b            real*8   I      lu factorization matrix
CD3 r            real*8   I/O    Residual array on input, solution
CD3                                  vector on output
CD3 nop          integer  I      Connectivity matrix
CD3 npvt         integer  I      Pivot positions in nop
CD3 piv          real*8   I      Array of pivots
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 NONE
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4 
CD4 NONE
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4 
CD4 NONE
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 NONE
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 neqm1        integer     neq-1
CD5 neqp1        integer     neq+1
CD5 k            integer     Do loop index
CD5 k1           integer     Do loop index
CD5 npivk1       integer     Index in pivot array
CD5 i1           integer     First position in solution array
CD5 jj           integer     Do loop index
CD5 nj           integer     Position in lu pointer array
CD5 idum         integer     Position in lu factorization matrix
CD5 i            integer     Do loop index
CD5 mpiv         integer     Index in pivot array
CD5 i2           integer     Second position in solution array
CD5 jq           integer     Do loop index
CD5 num          integer     Position in lu pointer array
CD5 inum         integer     Position in lu factorization matrix
CD5 pivot        real*8      Pivot value
CD5 sum          real*8      Term in calculation of factorization
CD5 c1           real*8      Term in calculation of factorization
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7 
CD7 N/A
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8 
CD8 3.2 Solve Linear Equation Set
CD8    3.2.3 Form Approximate Solution
CD8 
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN sub_bksub1
CPS 
CPS Define integer parameters
CPS 
CPS FOR each equation
CPS 
CPS   Identify diagonal term, multiply rhs by that term
CPS   
CPS   FOR all nonzero terms below the diagonal
CPS     Multiply lower diagonal matrix for the current row by right...
CPS     ... hand side, add to running sum
CPS   ENDFOR
CPS 
CPS ENDFOR each equation
CPS 
CPS Compute solution at last equation
CPS 
CPS FOR each equation, last-1 to first
CPS 
CPS   Identify pivot position
CPS   FOR all nonzero terms above the diagonal
CPS
CPS     Multiply upper diagonal matrix for the current row by right...
CPS     ... hand side, add to running sum
CPS     
CPS   ENDFOR
CPS 
CPS   Compute solution value for current equation
CPS 
CPS ENDFOR
CPS 
CPS END sub_bksub1
CPS
C**********************************************************************
c
c
c     back substitution routine for 3-1 dof problem

      implicit none

      integer neq
      real*8 a(*), b(*), r(*), piv(neq,*)
      integer ncon(*), nop(*),npvt(*),nrhs(*),na(*)
      integer neqm1, neqp1, k, k1, npivk1, i1, jj, nj, idum, i, mpiv
      integer j,kb,i2, jq, num, inum
      real*8 pivot, sum, c1
      integer icoupl,kk
      real*8 sto5(*)
      real*8 sum1,sum2,sum3,fdum2,tollr,overf
      
c
c     define some parameters
c
      neqm1=neq-1
      neqp1=neq+1
c
      if (icoupl.gt.0) then
        do i=1,neq
          sto5(i+nrhs(1))=r(nrhs(1)+i)
          sto5(i+nrhs(2))=r(nrhs(2)+i)
          sto5(i+nrhs(3))=r(nrhs(3)+i)
        enddo
      endif
c
c     use previously factored matrix
c
      do  k=1,neqm1
         pivot=piv(k,1)
         r(k)=r(k)*pivot
         k1=k+1
         npivk1=npvt(k1)-1
         i1=nop(k1)+1
c     factor row k1
         sum=0.0d00
         do jj=i1,npivk1
            nj=nop(jj)
            idum=jj-neqp1
            c1=b(idum)
            sum=sum-c1*r(nj)
         enddo
         r(k1)=r(k1)+sum
      enddo
c
c     back substitution
c
      r(neq)=r(neq)*piv(neq,1)
      do i=neqm1,1,-1
         mpiv=npvt(i)+1
         i2=nop(i+1)
         sum=0.0d00
         do jq=mpiv,i2
            num=nop(jq)
            inum=jq-neqp1
            sum=sum+b(inum)*r(num)
         enddo 
         r(i)=r(i)-sum
      enddo
c
c       back substitute to find r(nrhs(2)+i)
c
        do i=1,neq
           i1=ncon(i)+1
           i2=ncon(i+1)
           sum=r(nrhs(2)+i)
            do j=i1,i2
             k=j-neqp1
             kb=ncon(j)
             sum=sum-a(na(4)+k)*r(nrhs(1)+kb)
            enddo
           r(nrhs(2)+i)=sum
         enddo        
c
c       back substitute to find r(nrhs(3)+i)
c
        do i=1,neq
          i1=ncon(i)+1
          i2=ncon(i+1)
          sum=r(nrhs(3)+i)
          doj=i1,i2
            k=j-neqp1
            kb=ncon(j)
            sum=sum-a(na(7)+k)*r(nrhs(1)+kb)
          enddo
          r(nrhs(3)+i)=sum
        enddo
c
      if(icoupl.gt.0) then
      do  j=1,icoupl
      fdum2 = 0.0d00
        do  i=1,neq
          i1=ncon(i)+1
          i2=ncon(i+1)
          sum1=sto5(i+nrhs(1))
          sum2=sto5(i+nrhs(2))
          sum3=sto5(i+nrhs(3))
          do  k=i1,i2
            kk=k-neqp1
            kb=ncon(k)
            sum1=sum1-a(na(1)+kk)*r(nrhs(1)+kb)
     &        -a(na(2)+kk)*r(nrhs(2)+kb)
     &        -a(na(3)+kk)*r(nrhs(3)+kb)
            sum2=sum2-a(na(4)+kk)*r(nrhs(1)+kb)
     &        -a(na(5)+kk)*r(nrhs(2)+kb)
     &        -a(na(6)+kk)*r(nrhs(3)+kb)
            sum3=sum3-a(na(7)+kk)*r(nrhs(1)+kb)
     &        -a(na(8)+kk)*r(nrhs(2)+kb)
     &        -a(na(9)+kk)*r(nrhs(3)+kb)
          enddo
         r(nrhs(1)+i)=r(nrhs(1)+i)+overf*sum1
         r(nrhs(2)+i)=r(nrhs(2)+i)+overf*sum2
         r(nrhs(3)+i)=r(nrhs(3)+i)+overf*sum3
         fdum2=fdum2+abs(sum1)+abs(sum2)+abs(sum3)
        enddo
       if(fdum2.le.tollr) go to 500
      enddo
500   continue
      endif
      return
      end
