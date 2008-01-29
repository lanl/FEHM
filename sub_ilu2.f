      subroutine  sub_ilu2(neq,a,b,na,nb,ncon,nop
     *     ,irb,iirb,npvt,sorthm,dum1,dum2,piv) 
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
CD1 To perform incomplete lu factorization for 2 degree of freedom
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilu2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:40   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:40 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:41:26 1996   hend
CD2 Added Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:19:18   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:42   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:56   pvcs
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
CD5 nb1          integer     Pointer value in b array
CD5 nb2          integer     Pointer value in b array
CD5 nb3          integer     Pointer value in b array
CD5 nb4          integer     Pointer value in b array
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
CD5 pp           real*8      Scratch storage space
CD5 a1           real*8      Variable in function definition
CD5 a2           real*8      Variable in function definition
CD5 a3           real*8      Variable in function definition
CD5 a4           real*8      Variable in function definition
CD5 b1           real*8      Value in b array
CD5 b2           real*8      Value in b array
CD5 b3           real*8      Value in b array
CD5 b4           real*8      Value in b array
CD5 c1           real*8      Value used in calculation
CD5 c2           real*8      Value used in calculation
CD5 c3           real*8      Value used in calculation
CD5 c4           real*8      Value used in calculation
CD5 p1           real*8      Value in b array
CD5 p2           real*8      Value in b array
CD5 p3           real*8      Value in b array
CD5 p4           real*8      Value in b array
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
CPS BEGIN sub_ilu2
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
CPS END sub_ilu2
CPS
C**********************************************************************
c
c two degree of freedom solver
c
c     call ilu factorization at preconditioner stage

      implicit none
 
      integer neq
      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),sorthm(*),piv(neq,4)
      real*8 pp(4)
      real*8 dum1(neq,2),dum2(neq,2)
      integer na(*),nb(*)
      real*8 a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4
      real*8 p1, p2, p3, p4
      integer neqm1, neqp1, na1, na2, na3, na4, nb1, nb2, nb3, nb4
      integer i, i1, i2, ip, i3, i4, ik, ikb, k, ka, kb, l, npiv
      integer k1, npivk1, i5, nz, i6, kd, jj, nj, npivnj, ii, nk
      integer kj, ijk, idum
      real*8 alm, ali
c
c***     linear algebra     ***
c
      alm(a1,a2,b1,b2)=a1*b1+a2*b2
      ali(a1,a2,a3,a4,b1)=b1/(a1*a4-a2*a3)
c
c     define some parameters
c
      neqm1=neq-1
      neqp1=neq+1
      na1=na(1)
      na2=na(2)
      na3=na(3)
      na4=na(4)
      nb1=nb(1)
      nb2=nb(2)
      nb3=nb(3)
      nb4=nb(4)
c
      do i=1,neq
         i1=nop(i)+1
         i2=nop(i+1)
         ip=irb(i)
         i3=ncon(ip)+1
         i4=ncon(ip+1)
         do ik=i1,i2
            ikb=nop(ik)
            dum1(ikb,1)=0.0d00
            dum1(ikb,2)=0.0d00
            dum2(ikb,1)=0.0d00
            dum2(ikb,2)=0.0d00
         enddo
         do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            dum1(kb,1)=a(k-neqp1+na1)
            dum1(kb,2)=a(k-neqp1+na2)
            dum2(kb,1)=a(k-neqp1+na3)
            dum2(kb,2)=a(k-neqp1+na4)
         enddo
         do l=i1,i2
            kb=nop(l)
            b(l-neqp1+nb1)=dum1(kb,1)
            b(l-neqp1+nb2)=dum1(kb,2)
            b(l-neqp1+nb3)=dum2(kb,1)
            b(l-neqp1+nb4)=dum2(kb,2)
         enddo
      enddo
c     
c     start factorization
c     
      do  k=1,neqm1
         npiv=npvt(k)
         do i=1,4
            pp(i)=b(npiv-neqp1+nb(i))
         enddo
         piv(k,1)=ali(pp(1),pp(2),pp(3),pp(4),pp(4))
         piv(k,2)=-ali(pp(1),pp(2),pp(3),pp(4),pp(2))
         piv(k,3)=-ali(pp(1),pp(2),pp(3),pp(4),pp(3))
         piv(k,4)=ali(pp(1),pp(2),pp(3),pp(4),pp(1))
c     divide row k and rhs by pivot
         i2=nop(k+1)
         do  i=npiv,i2
            b1=b(i-neqp1+nb(1))
            b2=b(i-neqp1+nb(2))
            b3=b(i-neqp1+nb(3))
            b4=b(i-neqp1+nb(4))
            b(i-neqp1+nb(1))=alm(piv(k,1),piv(k,2),b1,b3)
            b(i-neqp1+nb(2))=alm(piv(k,1),piv(k,2),b2,b4)
            b(i-neqp1+nb(3))=alm(piv(k,3),piv(k,4),b1,b3)
            b(i-neqp1+nb(4))=alm(piv(k,3),piv(k,4),b2,b4)
         enddo
         k1=k+1
c     factor row k1
         npivk1=npvt(k1)-1
         i5=nop(k1)+1
         nz=0
         i6=nop(k1+1)
         do  i=i5,i6
            kd=nop(i)
            dum1(kd,1)=b(i-neqp1+nb(1))
            dum1(kd,2)=b(i-neqp1+nb(2))
            dum2(kd,1)=b(i-neqp1+nb(3))
            dum2(kd,2)=b(i-neqp1+nb(4))
         enddo
         do jj=i5,npivk1
            nj=nop(jj)
            nz=nz+1
            c1=dum1(nj,1)
            c2=dum1(nj,2)
            c3=dum2(nj,1)
            c4=dum2(nj,2)
            sorthm(nz)=c1
            sorthm(nz+neq)=c2
            sorthm(nz+neq+neq)=c3
            sorthm(nz+neq+neq+neq)=c4
c     add linear combo of row nj
            i4=nop(nj+1)
            npivnj=npvt(nj)
            do ii=npivnj,i4
               nk=nop(ii)
               p1=b(ii-neqp1+nb(1))
               p2=b(ii-neqp1+nb(2))
               p3=b(ii-neqp1+nb(3))
               p4=b(ii-neqp1+nb(4))
               dum1(nk,1)=dum1(nk,1)-alm(c1,c2,p1,p3)
               dum1(nk,2)=dum1(nk,2)-alm(c1,c2,p2,p4)
               dum2(nk,1)=dum2(nk,1)-alm(c3,c4,p1,p3)
               dum2(nk,2)=dum2(nk,2)-alm(c3,c4,p2,p4)
            enddo
         enddo
         do i=npivk1+1,i6
            kj=nop(i)
            b(i-neqp1+nb(1))=dum1(kj,1)
            b(i-neqp1+nb(2))=dum1(kj,2)
            b(i-neqp1+nb(3))=dum2(kj,1)
            b(i-neqp1+nb(4))=dum2(kj,2)
         enddo
         nz=0
         do  ijk=i5,npivk1
            nz=nz+1
            idum=ijk-neqp1
            b(idum+nb(1))=sorthm(nz)
            b(idum+nb(2))=sorthm(nz+neq)
            b(idum+nb(3))=sorthm(nz+neq+neq)
            b(idum+nb(4))=sorthm(nz+neq+neq+neq)
         enddo
      enddo 
      idum=nop(neqp1)-neqp1
      b1=b(idum+nb(1))
      b2=b(idum+nb(2))
      b3=b(idum+nb(3))
      b4=b(idum+nb(4))
      piv(neq,1)=ali(b1,b2,b3,b4,b4)
      piv(neq,2)=-ali(b1,b2,b3,b4,b2)
      piv(neq,3)=-ali(b1,b2,b3,b4,b3)
      piv(neq,4)=ali(b1,b2,b3,b4,b1)
c     
      return
      end
