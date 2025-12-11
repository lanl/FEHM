      subroutine sub_bksub62(neq,a,b,r,na,nb,nrhs,ncon,nop
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
CD1 To perform backsubstitution for a 2 degree of freedom (dof) problem.
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_bksub62.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:08   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:18   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:06 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:22:32 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:18:52   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:26   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:38   pvcs
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
CD3 nb           integer  I      Pointer array for positions in b array
CD3 nrhs         integer  I      Pointer array for positions of
CD3                                  solution vector
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
CD4 Identifier    Type     Description
CD4 
CD4 alm           real*8   Linear combination function
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
CD5 nr1          integer     Pointer for first unknown in solution
CD5                              vector
CD5 nr2          integer     Pointer for second unknown in solution
CD5                              vector
CD5 piv1         real*8      Pivot value
CD5 piv2         real*8      Pivot value
CD5 piv3         real*8      Pivot value
CD5 piv4         real*8      Pivot value
CD5 sum1         real*8      Term in calculation of factorization
CD5 sum2         real*8      Term in calculation of factorization
CD5 r1           real*8      Current residual value, first unknown
CD5 r2           real*8      Current residual value, second unknown
CD5 c1           real*8      Term in calculation of factorization
CD5 c2           real*8      Term in calculation of factorization
CD5 c3           real*8      Term in calculation of factorization
CD5 c4           real*8      Term in calculation of factorization
CD5 a1           real*8      Variable in definition of function alm
CD5 a2           real*8      Variable in definition of function alm
CD5 b1           real*8      Variable in definition of function alm
CD5 b2           real*8      Variable in definition of function alm
CD5 
CD5 Local Subprograms
CD5
CD5 
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
CD8 3.4.2     Solve Nonlinear Equation Set at Each Time Step
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
CPS BEGIN sub_bksub2
CPS 
CPS Define integer parameters
CPS 
CPS FOR each equation
CPS 
CPS   Identify diagonal terms, multiply rhs's by those terms
CPS   
CPS   FOR all nonzero terms below the diagonal
CPS     Multiply lower diagonal matrix for the current row by right...
CPS     ... hand side, add to running sum
CPS   ENDFOR
CPS 
CPS ENDFOR each equation
CPS 
CPS Compute solution at last equation, each unknown
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
CPS   Compute solution values for current equation, each unknown
CPS 
CPS ENDFOR
CPS 
CPS END sub_bksub2
CPS
C**********************************************************************
c
c
c     back substitution routine for 2dof problem

      implicit none

      integer neq
      real*8 a(*), b(*), r(*), piv(neq ,*)
      integer ncon(*), nop(*),npvt(*)
      integer na(*),nb(*),nrhs(*)
      integer neqm1, neqp1, k, k1, npivk1, i1, jj, nj, idum, i, mpiv
      integer j, kb, kk, i2, jq, num, inum, nr1, nr2
      real*8 sum1, sum2, c1, piv1, piv2, piv3, piv4, r1, r2
      real*8 c2, c3, c4, b3, b4
      real*8 a1, a2, b1, b2
      real*8 alm, fdum2
      real*8 sto5(*)
      real*8 sum3,sum4,sum5,sum6
      integer icoupl
      real*8 tollr,overf
c
c     linear algebra
c
      alm(a1,a2,b1,b2)=a1*b1+a2*b2
c
c     define some parameters
c
      neqm1=neq-1
      neqp1=neq+1     
      nr1=nrhs(1)
      nr2=nrhs(2)
c
      if (icoupl.gt.0) then
        do i=1,neq
          sto5(i+nrhs(1))=r(nrhs(1)+i)
          sto5(i+nrhs(2))=r(nrhs(2)+i)
          sto5(i+nrhs(3))=r(nrhs(3)+i)
          sto5(i+nrhs(4))=r(nrhs(4)+i)
          sto5(i+nrhs(5))=r(nrhs(5)+i)
          sto5(i+nrhs(6))=r(nrhs(6)+i)
        enddo
      endif
c
c
c     use previously factored matrix
c
      do k=1,neqm1
         piv1=piv(k,1)
         piv2=piv(k,2)
         piv3=piv(k,3)
         piv4=piv(k,4)
         r1=r(nr1+k)
         r2=r(nr2+k)
         r(nr1+k)=alm(piv1,piv2,r1,r2)
         r(nr2+k)=alm(piv3,piv4,r1,r2)
         k1=k+1
         npivk1=npvt(k1)-1
         i1=nop(k1)+1
c     factor row k1
         sum1=0.0
         sum2=0.0
         do jj=i1,npivk1
            nj=nop(jj)
            idum=jj-neqp1
            c1=b(idum+nb(1))
            c2=b(idum+nb(2))
            c3=b(idum+nb(3))
            c4=b(idum+nb(4))
            r1=r(nr1+nj)
            r2=r(nr2+nj)
            sum1=sum1-alm(c1,c2,r1,r2)
            sum2=sum2-alm(c3,c4,r1,r2)
         enddo
         r(nr1+k1)=r(nr1+k1)+sum1
         r(nr2+k1)=r(nr2+k1)+sum2
      enddo
c     
c     back substitution
c     
      c1=piv(neq,1)
      c2=piv(neq,2)
      c3=piv(neq,3)
      c4=piv(neq,4)
      r1=r(nr1+neq)
      r2=r(nr2+neq)
      r(nr1+neq)=alm(c1,c2,r1,r2)
      r(nr2+neq)=alm(c3,c4,r1,r2)
      do  i=neqm1,1,-1
         mpiv=npvt(i)+1
         i2=nop(i+1)
         sum1=0.0
         sum2=0.0
         do jq=mpiv,i2
            num=nop(jq)
            inum=jq-neqp1
            b1=b(inum+nb(1))
            b2=b(inum+nb(2))
            b3=b(inum+nb(3))
            b4=b(inum+nb(4))
            r1=r(nr1+num)
            r2=r(nr2+num)
            sum1=sum1+alm(b1,b2,r1,r2)
            sum2=sum2+alm(b3,b4,r1,r2)
         enddo
         r(nr1+i)=r(nr1+i)-sum1
         r(nr2+i)=r(nr2+i)-sum2
      enddo
c     
        do i=1,neq
           i1=ncon(i)+1
           i2=ncon(i+1)
           sum1=r(nrhs(3)+i)
           sum2=r(nrhs(4)+i)
           sum3=r(nrhs(5)+i)
           sum4=r(nrhs(6)+i)
            do j=i1,i2
             k=j-neqp1
             kb=ncon(j)
             sum1=sum1-a(na(13)+k)*r(nrhs(1)+kb)
     &                -a(na(14)+k)*r(nrhs(2)+kb)
             sum2=sum2-a(na(19)+k)*r(nrhs(1)+kb)
     &                -a(na(20)+k)*r(nrhs(2)+kb)
             sum3=sum3-a(na(25)+k)*r(nrhs(1)+kb)
     &                -a(na(26)+k)*r(nrhs(2)+kb)
             sum4=sum4-a(na(31)+k)*r(nrhs(1)+kb)
     &                -a(na(32)+k)*r(nrhs(2)+kb)
            enddo
           r(nrhs(3)+i)=sum1
           r(nrhs(4)+i)=sum2
           r(nrhs(5)+i)=sum3
           r(nrhs(6)+i)=sum4
         enddo
c
      do  j=1,icoupl
        fdum2=0.0d00
        do  i=1,neq
          i1=ncon(i)+1
          i2=ncon(i+1)
          sum1=sto5(i+nrhs(1))
          sum2=sto5(i+nrhs(2))
          sum3=sto5(i+nrhs(3))
          sum4=sto5(i+nrhs(4))
          sum5=sto5(i+nrhs(5))
          sum6=sto5(i+nrhs(6))
          do  k=i1,i2
            kk=k-neqp1
            kb=ncon(k)
            sum1=sum1-a(na(1)+kk)*r(nrhs(1)+kb)
     &        -a(na(2)+kk)*r(nrhs(2)+kb)
     &        -a(na(3)+kk)*r(nrhs(3)+kb)
     &        -a(na(4)+kk)*r(nrhs(4)+kb)
     &        -a(na(5)+kk)*r(nrhs(5)+kb)
     &        -a(na(6)+kk)*r(nrhs(6)+kb)
            sum2=sum2-a(na(7)+kk)*r(nrhs(1)+kb)
     &        -a(na(8)+kk)*r(nrhs(2)+kb)
     &        -a(na(9)+kk)*r(nrhs(3)+kb)
     &        -a(na(10)+kk)*r(nrhs(4)+kb)
     &        -a(na(11)+kk)*r(nrhs(5)+kb)
     &        -a(na(12)+kk)*r(nrhs(6)+kb)
            sum3=sum3-a(na(13)+kk)*r(nrhs(1)+kb)
     &        -a(na(14)+kk)*r(nrhs(2)+kb)
     &        -a(na(15)+kk)*r(nrhs(3)+kb)
     &        -a(na(16)+kk)*r(nrhs(4)+kb)
     &        -a(na(17)+kk)*r(nrhs(5)+kb)
     &        -a(na(18)+kk)*r(nrhs(6)+kb)
            sum4=sum4-a(na(19)+kk)*r(nrhs(1)+kb)
     &        -a(na(20)+kk)*r(nrhs(2)+kb)
     &        -a(na(21)+kk)*r(nrhs(3)+kb)
     &        -a(na(22)+kk)*r(nrhs(4)+kb)
     &        -a(na(23)+kk)*r(nrhs(5)+kb)
     &        -a(na(24)+kk)*r(nrhs(6)+kb)
            sum5=sum5-a(na(25)+kk)*r(nrhs(1)+kb)
     &        -a(na(26)+kk)*r(nrhs(2)+kb)
     &        -a(na(27)+kk)*r(nrhs(3)+kb)
     &        -a(na(28)+kk)*r(nrhs(4)+kb)
     &        -a(na(29)+kk)*r(nrhs(5)+kb)
     &        -a(na(30)+kk)*r(nrhs(6)+kb)
            sum6=sum6-a(na(31)+kk)*r(nrhs(1)+kb)
     &        -a(na(32)+kk)*r(nrhs(2)+kb)
     &        -a(na(33)+kk)*r(nrhs(3)+kb)
     &        -a(na(34)+kk)*r(nrhs(4)+kb)
     &        -a(na(35)+kk)*r(nrhs(5)+kb)
     &        -a(na(36)+kk)*r(nrhs(6)+kb)
          enddo
         r(nrhs(1)+i)=r(nrhs(1)+i)+overf*sum1
         r(nrhs(2)+i)=r(nrhs(2)+i)+overf*sum2
         r(nrhs(3)+i)=r(nrhs(3)+i)+overf*sum3
         r(nrhs(4)+i)=r(nrhs(4)+i)+overf*sum4
         r(nrhs(5)+i)=r(nrhs(5)+i)+overf*sum5
         r(nrhs(6)+i)=r(nrhs(6)+i)+overf*sum6
         fdum2=fdum2+abs(sum1)+abs(sum2)+abs(sum3)
     &              +abs(sum4)+abs(sum5)+abs(sum6)
        enddo
        if(fdum2.le.tollr) go to 500
      enddo
500   continue

      return
      end        
