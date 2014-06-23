      subroutine  sub_ilu3(neq,a,b,na,nb,ncon,nop
     *     ,irb,iirb,npvt,sorthm,dum1,dum2,dum3,piv) 
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
CD1 To perform incomplete lu factorization for 3 degree of freedom
CD1 problem.
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilu3.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:42 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:42:22 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:19:22   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:44   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:58   pvcs
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
CD3 a            real*8   I      Solution matrix
CD3 b            real*8   O      lu factorization matrix
CD3 na           integer  I      Pointer array in a matrix
CD3 nb           integer  I      Pointer array in b matrix
CD3 ncon         integer  I      Connectivity matrix for solution matrix
CD3 nop          integer  I      Connectivity matrix for factorization
CD3                                 matrix
CD3 irb          integer  I      Renumber array - inew=irb(iold)
CD3 iirb         integer  I      Renumber array - iold=iirb(inew)
CD3 npvt         integer  I      Pivot positions in nop
CD3 sorthm       real*8   I      Scratch storage array
CD3 dum1         real*8   I      Scratch storage array
CD3 dum2         real*8   I      Scratch storage array
CD3 dum3         real*8   I      Scratch storage array
CD3 piv          real*8   O      Array of pivots
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
CD5 na1          integer     Pointer value in a array
CD5 na2          integer     Pointer value in a array
CD5 na3          integer     Pointer value in a array
CD5 na4          integer     Pointer value in a array
CD5 na5          integer     Pointer value in a array
CD5 na6          integer     Pointer value in a array
CD5 na7          integer     Pointer value in a array
CD5 na8          integer     Pointer value in a array
CD5 na9          integer     Pointer value in a array
CD5 nb1          integer     Pointer value in b array
CD5 nb2          integer     Pointer value in b array
CD5 nb3          integer     Pointer value in b array
CD5 nb4          integer     Pointer value in b array
CD5 nb5          integer     Pointer value in b array
CD5 nb6          integer     Pointer value in b array
CD5 nb7          integer     Pointer value in b array
CD5 nb8          integer     Pointer value in b array
CD5 nb9          integer     Pointer value in b array
CD5 i            integer     Do loop index
CD5 i1           integer     Do loop index
CD5 i2           integer     Do loop index
CD5 i3           integer     Do loop index
CD5 i4           integer     Do loop index
CD5 i5           integer     Do loop index
CD5 i6           integer     Do loop index
CD5 ip           integer     Pointer to reorder array
CD5 ik           integer     Do loop index
CD5 ikb          integer     Integer index
CD5 k            integer     Do loop index
CD5 ka           integer     Pointer in connectivity array
CD5 kb           integer     Pointer to reorder array
CD5 l            integer     Do loop index
CD5 npiv         integer     Pivot position
CD5 k1           integer     k+1
CD5 npivk1       integer     Index in pivot array
CD5 nz           integer     Pointer to multiplier of previous row
CD5 kd           integer     Integer index
CD5 jj           integer     Do loop index
CD5 nj           integer     Position in lu pointer array
CD5 npivnj       integer     Index in pivot array
CD5 ii           integer     Do loop index
CD5 nk           integer     Position in factorization connectivity
CD5                             array
CD5 ijk          integer     Do loop index
CD5 kj           integer     Integer index
CD5 idum         integer     Integer index
CD5 iii          integer     Integer index
CD5 pp           real*8      Scratch storage space
CD5 a1           real*8      Variable in function definition
CD5 a2           real*8      Variable in function definition
CD5 a3           real*8      Variable in function definition
CD5 a4           real*8      Variable in function definition
CD5 a5           real*8      Variable in function definition
CD5 b1           real*8      Value in b array
CD5 b2           real*8      Value in b array
CD5 b3           real*8      Value in b array
CD5 b4           real*8      Value in b array
CD5 b5           real*8      Value in b array
CD5 b6           real*8      Value in b array
CD5 b7           real*8      Value in b array
CD5 b8           real*8      Value in b array
CD5 b9           real*8      Value in b array
CD5 c1           real*8      Value used in calculation
CD5 c2           real*8      Value used in calculation
CD5 c3           real*8      Value used in calculation
CD5 c4           real*8      Value used in calculation
CD5 c5           real*8      Value used in calculation
CD5 c6           real*8      Value used in calculation
CD5 c7           real*8      Value used in calculation
CD5 c8           real*8      Value used in calculation
CD5 c9           real*8      Value used in calculation
CD5 p1           real*8      Value in b array
CD5 p2           real*8      Value in b array
CD5 p3           real*8      Value in b array
CD5 p4           real*8      Value in b array
CD5 p5           real*8      Value in b array
CD5 p6           real*8      Value in b array
CD5 p7           real*8      Value in b array
CD5 p8           real*8      Value in b array
CD5 p9           real*8      Value in b array
CD5 det          real*8      Determinant used in calculation
CD5 
CD5 Local Subprograms
CD5
CD5 Identifier    Type      Description
CD5 
CD5 alm           real*8    Linear combination function
CD5 ali           real*8    Linear combination function
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
CD8    3.2.1 Perform Incomplete Factorization
CD8 3.3 Provide Multiple Degree of Freedom Option
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
CPS BEGIN sub_ilu3
CPS 
CPS Define integer parameters
CPS 
CPS FOR each equation
CPS 
CPS   Determine index parameters
CPS   
CPS   For each position in dummy storage arrays
CPS     Initialize to zero
CPS   ENDFOR
CPS   
CPS   FOR each nonzero position in solution array
CPS     Copy solution matrix values to dummy storage arrays
CPS   ENDFOR
CPS   
CPS   FOR each nonzero position in lu factorization  array
CPS     Copy dummy storage array values into lu factorization array
CPS   ENDFOR
CPS   
CPS ENDFOR each equation
CPS 
CPS FOR each equation except the last
CPS 
CPS   FOR each value in b array for this equation
CPS     Store b value in scratch storage array
CPS   ENDFOR
CPS   
CPS   Identify pivot position
CPS   Invert diagonal element to form pivot
CPS   
CPS   Save pivot value in arrays
CPS   
CPS   FOR each nonzero element in current equation
CPS     Compute values in b array
CPS   ENDFOR
CPS   
CPS   FOR each nonzero element of lu matrix
CPS     Copy value into dummy storage array
CPS   ENDFOR
CPS   
CPS   FOR each nonzero element in lower diagonal part of lu matrix
CPS   
CPS     Save term
CPS     Identify column number of that term
CPS     
CPS     FOR each column number
CPS       Add linear combination (place multiplier in array sorthm)...
CPS       ... of previous row with row number equal to this column...
CPS       ... number (place in array dum)
CPS     ENDFOR
CPS     
CPS   ENDFOR each nonzero element in lower diagonal part of lu matrix
CPS   
CPS   FOR each nonzero position in lu factorization matrix
CPS     Copy value in dum to lu factorization matrix
CPS   ENDFOR
CPS   
CPS   FOR each nonzero value of lu factoriation matrix below diagonal
CPS     Copy value in sorthm to lu factorization matrix
CPS   ENDFOR
CPS   
CPS ENDFOR each equation except the last
CPS 
CPS Compute pivot for last equation
CPS 
CPS END sub_ilu3
CPS
C**********************************************************************
c
c three degree of freedom solver
c
c call ilu factorization at preconditioner stage
c

      implicit none

      integer neq
      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),sorthm(*),piv(neq,9)
      real*8 pp(9)
      real*8 dum1(neq,3),dum2(neq,3),dum3(neq,3)
      integer na(*),nb(*)
      real*8 a1, a2, a3, a4, a5, b1, b2, b3, b4, c1, c2, c3, c4
      real*8 p1, p2, p3, p4, det, b5, b6, b7, b8, b9
      real*8 c5, c6, c7, c8, c9, p5, p6, p7, p8, p9
      integer neqm1, neqp1, na1, na2, na3, na4, nb1, nb2, nb3, nb4
      integer na5, na6, na7, na8, na9, nb5, nb6, nb7, nb8, nb9
      integer i, i1, i2, ip, i3, i4, ik, ikb, k, ka, kb, l, npiv
      integer k1, npivk1, i5, nz, i6, kd, jj, nj, npivnj, ii, nk
      integer kj, ijk, idum, iii
      real*8 ali, alm
c
c***     linear algebra     ***
c
      ali(a1,a2,a3,a4,a5)=(a1*a2-a3*a4)/a5
      alm(a1,a2,a3,b1,b2,b3)=a1*b1+a2*b2+a3*b3
c     
      neqp1=neq+1
      neqm1=neq-1
      na1=na(1)
      na2=na(2)
      na3=na(3)
      na4=na(4)
      na5=na(5)
      na6=na(6)
      na7=na(7)
      na8=na(8)
      na9=na(9)
      nb1=nb(1)
      nb2=nb(2)
      nb3=nb(3)
      nb4=nb(4)
      nb5=nb(5)
      nb6=nb(6)
      nb7=nb(7)
      nb8=nb(8)
      nb9=nb(9)
c
      do i=1,neq
        i1=nop(i)+1
        i2=nop(i+1)
        ip=iirb(i)
        i3=ncon(ip)+1
        i4=ncon(ip+1)
          do ik=i1,i2
            ikb=nop(ik)
            dum1(ikb,1)=0.0d00
            dum1(ikb,2)=0.0d00
            dum1(ikb,3)=0.0d00
            dum2(ikb,1)=0.0d00
            dum2(ikb,2)=0.0d00
            dum2(ikb,3)=0.0d00
            dum3(ikb,1)=0.0d00
            dum3(ikb,2)=0.0d00
            dum3(ikb,3)=0.0d00
          enddo
          do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            dum1(kb,1)=a(k-neqp1+na1)
            dum1(kb,2)=a(k-neqp1+na2)
            dum1(kb,3)=a(k-neqp1+na3)
            dum2(kb,1)=a(k-neqp1+na4)
            dum2(kb,2)=a(k-neqp1+na5)
            dum2(kb,3)=a(k-neqp1+na6)
            dum3(kb,1)=a(k-neqp1+na7)
            dum3(kb,2)=a(k-neqp1+na8)
            dum3(kb,3)=a(k-neqp1+na9)
          enddo
          do l=i1,i2
            kb=nop(l)
            b(l-neqp1+nb1)=dum1(kb,1)
            b(l-neqp1+nb2)=dum1(kb,2)
            b(l-neqp1+nb3)=dum1(kb,3)
            b(l-neqp1+nb4)=dum2(kb,1)
            b(l-neqp1+nb5)=dum2(kb,2)
            b(l-neqp1+nb6)=dum2(kb,3)
            b(l-neqp1+nb7)=dum3(kb,1)
            b(l-neqp1+nb8)=dum3(kb,2)
            b(l-neqp1+nb9)=dum3(kb,3)
          enddo
      enddo
c
c start factorization
c
      do k=1,neqm1
        npiv=npvt(k)
        do i=1,9
          pp(i)=b(npiv-neqp1+nb(i))
        enddo
        det = pp(1)*pp(5)*pp(9) + pp(2)*pp(6)*pp(7) + pp(3)*pp(4)*pp(8)
     *      - pp(3)*pp(5)*pp(7) - pp(1)*pp(6)*pp(8) - pp(2)*pp(4)*pp(9)
        piv(k,1) = ali(pp(5),pp(9),pp(6),pp(8),det)
        piv(k,2) = ali(pp(3),pp(8),pp(2),pp(9),det)
        piv(k,3) = ali(pp(2),pp(6),pp(3),pp(5),det)
        piv(k,4) = ali(pp(6),pp(7),pp(4),pp(9),det)
        piv(k,5) = ali(pp(1),pp(9),pp(3),pp(7),det)
        piv(k,6) = ali(pp(3),pp(4),pp(1),pp(6),det)
        piv(k,7) = ali(pp(4),pp(8),pp(5),pp(7),det)
        piv(k,8) = ali(pp(2),pp(7),pp(1),pp(8),det)
        piv(k,9) = ali(pp(1),pp(5),pp(2),pp(4),det)
c divide row k and rhs by pivot
        i2=nop(k+1)
        do i=npiv,i2
          iii=i-neqp1
          b1=b(iii+nb(1))
          b2=b(iii+nb(2))
          b3=b(iii+nb(3))
          b4=b(iii+nb(4))
          b5=b(iii+nb(5))
          b6=b(iii+nb(6))
          b7=b(iii+nb(7))
          b8=b(iii+nb(8))
          b9=b(iii+nb(9))
          b(iii+nb(1))=alm(piv(k,1),piv(k,2),piv(k,3),b1,b4,b7)
          b(iii+nb(2))=alm(piv(k,1),piv(k,2),piv(k,3),b2,b5,b8)
          b(iii+nb(3))=alm(piv(k,1),piv(k,2),piv(k,3),b3,b6,b9)
          b(iii+nb(4))=alm(piv(k,4),piv(k,5),piv(k,6),b1,b4,b7)
          b(iii+nb(5))=alm(piv(k,4),piv(k,5),piv(k,6),b2,b5,b8)
          b(iii+nb(6))=alm(piv(k,4),piv(k,5),piv(k,6),b3,b6,b9)
          b(iii+nb(7))=alm(piv(k,7),piv(k,8),piv(k,9),b1,b4,b7)
          b(iii+nb(8))=alm(piv(k,7),piv(k,8),piv(k,9),b2,b5,b8)
          b(iii+nb(9))=alm(piv(k,7),piv(k,8),piv(k,9),b3,b6,b9)
        enddo
        k1=k+1
c factor row k1
        npivk1=npvt(k1)-1
        i5=nop(k1)+1
        nz=0
        i6=nop(k1+1)
        do i=i5,i6
          kd=nop(i)
          dum1(kd,1)=b(i-neqp1+nb(1))
          dum1(kd,2)=b(i-neqp1+nb(2))
          dum1(kd,3)=b(i-neqp1+nb(3))
          dum2(kd,1)=b(i-neqp1+nb(4))
          dum2(kd,2)=b(i-neqp1+nb(5))
          dum2(kd,3)=b(i-neqp1+nb(6))
          dum3(kd,1)=b(i-neqp1+nb(7))
          dum3(kd,2)=b(i-neqp1+nb(8))
          dum3(kd,3)=b(i-neqp1+nb(9))
        enddo
        do jj=i5,npivk1
          nj=nop(jj)
          nz=nz+1
          c1=dum1(nj,1)
          c2=dum1(nj,2)
          c3=dum1(nj,3)
          c4=dum2(nj,1)
          c5=dum2(nj,2)
          c6=dum2(nj,3)
          c7=dum3(nj,1)
          c8=dum3(nj,2)
          c9=dum3(nj,3)
          sorthm(nz)=c1
          sorthm(nz+neq)=c2
          sorthm(nz+neq+neq)=c3
          sorthm(nz+neq+neq+neq)=c4
          sorthm(nz+neq+neq+neq+neq)=c5
          sorthm(nz+neq+neq+neq+neq+neq)=c6
          sorthm(nz+neq+neq+neq+neq+neq+neq)=c7
          sorthm(nz+neq+neq+neq+neq+neq+neq+neq)=c8
          sorthm(nz+neq+neq+neq+neq+neq+neq+neq+neq)=c9
c add linear combo of row nj
          i4=nop(nj+1)
          npivnj=npvt(nj)
          do ii=npivnj,i4
            nk=nop(ii)
            p1=b(ii-neqp1+nb(1))
            p2=b(ii-neqp1+nb(2))
            p3=b(ii-neqp1+nb(3))
            p4=b(ii-neqp1+nb(4))
            p5=b(ii-neqp1+nb(5))
            p6=b(ii-neqp1+nb(6))
            p7=b(ii-neqp1+nb(7))
            p8=b(ii-neqp1+nb(8))
            p9=b(ii-neqp1+nb(9))
            dum1(nk,1)=dum1(nk,1)-alm(c1,c2,c3,p1,p4,p7)
            dum1(nk,2)=dum1(nk,2)-alm(c1,c2,c3,p2,p5,p8)
            dum1(nk,3)=dum1(nk,3)-alm(c1,c2,c3,p3,p6,p9)
            dum2(nk,1)=dum2(nk,1)-alm(c4,c5,c6,p1,p4,p7)
            dum2(nk,2)=dum2(nk,2)-alm(c4,c5,c6,p2,p5,p8)
            dum2(nk,3)=dum2(nk,3)-alm(c4,c5,c6,p3,p6,p9)
            dum3(nk,1)=dum3(nk,1)-alm(c7,c8,c9,p1,p4,p7)
            dum3(nk,2)=dum3(nk,2)-alm(c7,c8,c9,p2,p5,p8)
            dum3(nk,3)=dum3(nk,3)-alm(c7,c8,c9,p3,p6,p9)
          enddo
        enddo
        do i=npivk1+1,i6
          kj=nop(i)
          b(i-neqp1+nb(1))=dum1(kj,1)
          b(i-neqp1+nb(2))=dum1(kj,2)
          b(i-neqp1+nb(3))=dum1(kj,3)
          b(i-neqp1+nb(4))=dum2(kj,1)
          b(i-neqp1+nb(5))=dum2(kj,2)
          b(i-neqp1+nb(6))=dum2(kj,3)
          b(i-neqp1+nb(7))=dum3(kj,1)
          b(i-neqp1+nb(8))=dum3(kj,2)
          b(i-neqp1+nb(9))=dum3(kj,3)
        enddo
        nz=0
        do ijk=i5,npivk1
          nz=nz+1
          idum=ijk-neqp1
          b(idum+nb(1))=sorthm(nz)
          b(idum+nb(2))=sorthm(nz+neq)
          b(idum+nb(3))=sorthm(nz+neq+neq)
          b(idum+nb(4))=sorthm(nz+neq+neq+neq)
          b(idum+nb(5))=sorthm(nz+neq+neq+neq+neq)
          b(idum+nb(6))=sorthm(nz+neq+neq+neq+neq+neq)
          b(idum+nb(7))=sorthm(nz+neq+neq+neq+neq+neq+neq)
          b(idum+nb(8))=sorthm(nz+neq+neq+neq+neq+neq+neq+neq)
          b(idum+nb(9))=sorthm(nz+neq+neq+neq+neq+neq+neq+neq+neq)
        enddo
      enddo
      idum=nop(neqp1)-neqp1
      do i=1,9
        pp(i)=b(idum+nb(i))
      enddo
      det = pp(1)*pp(5)*pp(9) + pp(2)*pp(6)*pp(7) + pp(3)*pp(4)*pp(8)
     *    - pp(3)*pp(5)*pp(7) - pp(1)*pp(6)*pp(8) - pp(2)*pp(4)*pp(9)
      piv(neq,1)=ali(pp(5),pp(9),pp(6),pp(8),det)
      piv(neq,2)=ali(pp(3),pp(8),pp(2),pp(9),det)
      piv(neq,3)=ali(pp(2),pp(6),pp(3),pp(5),det)
      piv(neq,4)=ali(pp(6),pp(7),pp(4),pp(9),det)
      piv(neq,5)=ali(pp(1),pp(9),pp(3),pp(7),det)
      piv(neq,6)=ali(pp(3),pp(4),pp(1),pp(6),det)
      piv(neq,7)=ali(pp(4),pp(8),pp(5),pp(7),det)
      piv(neq,8)=ali(pp(2),pp(7),pp(1),pp(8),det)
      piv(neq,9)=ali(pp(1),pp(5),pp(2),pp(4),det)
c
        return
        end
