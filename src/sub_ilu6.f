      subroutine  sub_ilu6(nel,a,b,na,nb,ncon,nop
     *     ,irb,iirb,npvt,dum1,dum2,dum3,dum4,dum5,dum6,piv)
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
!C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To perform incomplete lu factorization for 6 degree of freedom
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
CD2 6-24-94      B. Robinson    97      Made 4 dof routine into 6 dof
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilu6.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Thu Sep 19 12:22:08 1996   llt
CD2 added log history
CD2
CD2    Rev 1.2   09/12/96 08:26:52   robinson
CD2 Prolog Changes
CD2
CD2    Rev 1.1   02/02/96 12:44:16   hend
CD2 Updated Requirements Traceability
CD2
CD2    Rev 1.0   01/16/96 14:48:42   llt
CD2 Initial revision.
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 nel          integer  I      Number of entries in the array for
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
CD3 dum1         real*8   I      Scratch storage used in calculation
CD3 dum2         real*8   I      Scratch storage used in calculation
CD3 dum3         real*8   I      Scratch storage used in calculation
CD3 dum4         real*8   I      Scratch storage used in calculation
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
CD5 nelp1        integer     nel+1
CD5 i1           integer     Do loop limit parameter
CD5 i2           integer     Do loop limit parameter
CD5 i3           integer     Do loop limit parameter
CD5 i4           integer     Do loop limit parameter
CD5 na11         integer     Pointer position in a array
CD5 na12         integer     Pointer position in a array
CD5 na13         integer     Pointer position in a array
CD5 na14         integer     Pointer position in a array
CD5 na21         integer     Pointer position in a array
CD5 na22         integer     Pointer position in a array
CD5 na23         integer     Pointer position in a array
CD5 na24         integer     Pointer position in a array
CD5 na31         integer     Pointer position in a array
CD5 na32         integer     Pointer position in a array
CD5 na33         integer     Pointer position in a array
CD5 na34         integer     Pointer position in a array
CD5 na41         integer     Pointer position in a array
CD5 na42         integer     Pointer position in a array
CD5 na43         integer     Pointer position in a array
CD5 na44         integer     Pointer position in a array
CD5 i            integer     Do loop index
CD5 ifin         integer     Do loop limit parameter
CD5 ja           integer     Integer index
CD5 ij           integer     Do loop index
CD5 ijind        integer     Integer index
CD5 im1nq2       integer     Integer index
CD5 ipvt         integer     Integer index
CD5 ista         integer     Do loop limit parameter
CD5 ipvtp1       integer     Integer index
CD5 ipvtm1       integer     Integer index
CD5 j            integer     Integer index
CD5 ik           integer     Do loop index
CD5 k            integer     Do loop index
CD5 kj           integer     Integer index
CD5 l            integer     Do loop index
CD5 j2           integer     Integer index
CD5 kpvt         integer     Integer index
CD5 ikind        integer     Integer index
CD5 kjind        integer     Integer index
CD5 ik           integer     Integer index
CD5 ipvind       integer     Integer index
CD5 a1           real*8      Variable in function definition
CD5 a2           real*8      Variable in function definition
CD5 a3           real*8      Variable in function definition
CD5 a4           real*8      Variable in function definition
CD5 b1           real*8      Variable in function definition
CD5 b2           real*8      Variable in function definition
CD5 dd11         real*8      Value in b array
CD5 dd12         real*8      Value in b array
CD5 dd13         real*8      Value in b array
CD5 dd14         real*8      Value in b array
CD5 dd21         real*8      Value in b array
CD5 dd22         real*8      Value in b array
CD5 dd23         real*8      Value in b array
CD5 dd24         real*8      Value in b array
CD5 dd31         real*8      Value in b array
CD5 dd32         real*8      Value in b array
CD5 dd33         real*8      Value in b array
CD5 dd34         real*8      Value in b array
CD5 dd41         real*8      Value in b array
CD5 dd42         real*8      Value in b array
CD5 dd43         real*8      Value in b array
CD5 dd44         real*8      Value in b array
CD5 bksb1        real*8      Term used in calculation
CD5 bksb2        real*8      Term used in calculation
CD5 bksb3        real*8      Term used in calculation
CD5 bksb4        real*8      Term used in calculation
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
CD7 Note that although this routine uses essentially the same algorithm
CD7 as sub_ilu1, sub_ilu2, and sub_ilu3, the actual
CD7 implementation is somewhat different for the sake of computational
CD7 efficiency.  Thus the code structure is somewhat different than
CD7 these other routines.
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
CPS BEGIN sub_ilu6
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
CPS FOR each equation
CPS 
CPS   Identify pivot position
CPS   
CPS   FOR each nonzero element of lu factorization matrix below the...
CPS   ... diagonal (column is j, term is l(i,j))
CPS   
CPS     FOR each nonzero element of lu factorization matrix up to...
CPS     ... the current one (position k, term is l(i,k))
CPS     
CPS       LOOP to find position k,j for u(k,j)
CPS       EXITIF we have found it
CPS         Decrease search index in lu factorization matrix by 1
CPS       ENDLOOP
CPS       
CPS       IF we found a match
CPS         Subtract products of l(i,k)*u(k,j) from l(i,j)
CPS       ENDIF
CPS     ENDFOR
CPS   ENDFOR
CPS   
CPS   Calculate inverse of l(i,i)
CPS   
CPS   FOR each nonzero element of lu factorization matrix above the...
CPS   ... diagonal (column is j, term is u(i,j))
CPS   
CPS     FOR each nonzero element of lu factorization matrix up to...
CPS     ... the current one (position k, term is l(i,k))
CPS     
CPS       LOOP to find position k,j for u(k,j)
CPS       EXITIF we have found it
CPS         Decrease search index in lu factorization matrix by 1
CPS       ENDLOOP
CPS       
CPS       IF we found a match
CPS         Subtract products of l(i,k)*u(k,j) from u(i,j)
CPS       ENDIF
CPS     ENDFOR
CPS     
CPS     Divide term u(i,j) by diagonal value (l(i,i))
CPS     
CPS   ENDFOR
CPS   
CPS ENDFOR
CPS 
CPS END sub_ilu6
CPS
C**********************************************************************
c
c six degree of freedom solver
c
c     6 degree of freedom ilu preconditioner
   
      implicit none

      integer nel
      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),piv(*)
      real*8 dum1(nel,6),dum2(nel,6),dum3(nel,6),dum4(nel,6)
      real*8 dum5(nel,6),dum6(nel,6)
      integer na(*),nb(*)                    
      integer nelp1
      integer i1, i2, i3, i4
      integer na11, na12, na13, na14, na15, na16
      integer na21, na22, na23, na24, na25, na26
      integer na31, na32, na33, na34, na35, na36
      integer na41, na42, na43, na44, na45, na46
      integer na51, na52, na53, na54, na55, na56
      integer na61, na62, na63, na64, na65, na66
      integer nb11, nb12, nb13, nb14, nb15, nb16
      integer nb21, nb22, nb23, nb24, nb25, nb26
      integer nb31, nb32, nb33, nb34, nb35, nb36
      integer nb41, nb42, nb43, nb44, nb45, nb46
      integer nb51, nb52, nb53, nb54, nb55, nb56
      integer nb61, nb62, nb63, nb64, nb65, nb66
      integer i, ifin, ij, ijind, ikb, ka, kb
      integer im1nq2, ipvt, ista, ipvtp1, ipvtm1, j, ik, k, kj, l
      integer j2, kpvt, ikind, kjind
      integer ipvind, ip
      real*8 dd11, dd12, dd13, dd14, dd15, dd16
      real*8 dd21, dd22, dd23, dd24, dd25, dd26
      real*8 dd31, dd32, dd33, dd34, dd35, dd36
      real*8 dd41, dd42, dd43, dd44, dd45, dd46
      real*8 dd51, dd52, dd53, dd54, dd55, dd56
      real*8 dd61, dd62, dd63, dd64, dd65, dd66
      real*8 bksb1, bksb2, bksb3, bksb4, bksb5, bksb6
c       
      nelp1=nel+1
      na11=na( 1)
      na12=na( 2)
      na13=na( 3)
      na14=na( 4)
      na15=na( 5)
      na16=na( 6)
      na21=na( 7)
      na22=na( 8)
      na23=na( 9)
      na24=na(10)
      na25=na(11)
      na26=na(12)
      na31=na(13)
      na32=na(14)
      na33=na(15)
      na34=na(16)
      na35=na(17)
      na36=na(18)
      na41=na(19)
      na42=na(20)
      na43=na(21)
      na44=na(22)
      na45=na(23)
      na46=na(24)
      na51=na(25)
      na52=na(26)
      na53=na(27)
      na54=na(28)
      na55=na(29)
      na56=na(30)
      na61=na(31)
      na62=na(32)
      na63=na(33)
      na64=na(34)
      na65=na(35)
      na66=na(36)
      nb11=nb( 1)
      nb12=nb( 2)
      nb13=nb( 3)
      nb14=nb( 4)
      nb15=nb( 5)
      nb16=nb( 6)
      nb21=nb( 7)
      nb22=nb( 8)
      nb23=nb( 9)
      nb24=nb(10)
      nb25=nb(11)
      nb26=nb(12)
      nb31=nb(13)
      nb32=nb(14)
      nb33=nb(15)
      nb34=nb(16)
      nb35=nb(17)
      nb36=nb(18)
      nb41=nb(19)
      nb42=nb(20)
      nb43=nb(21)
      nb44=nb(22)
      nb45=nb(23)
      nb46=nb(24)
      nb51=nb(25)
      nb52=nb(26)
      nb53=nb(27)
      nb54=nb(28)
      nb55=nb(29)
      nb56=nb(30)
      nb61=nb(31)
      nb62=nb(32)
      nb63=nb(33)
      nb64=nb(34)
      nb65=nb(35)
      nb66=nb(36)
c     initialize b, which is equal to a but has extra storage for
c       storing the lu fill-in
      do i=1,nel
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
            dum1(ikb,4)=0.0d00
            dum1(ikb,5)=0.0d00
            dum1(ikb,6)=0.0d00
            dum2(ikb,1)=0.0d00
            dum2(ikb,2)=0.0d00
            dum2(ikb,3)=0.0d00
            dum2(ikb,4)=0.0d00
            dum2(ikb,5)=0.0d00
            dum2(ikb,6)=0.0d00
            dum3(ikb,1)=0.0d00
            dum3(ikb,2)=0.0d00
            dum3(ikb,3)=0.0d00
            dum3(ikb,4)=0.0d00
            dum3(ikb,5)=0.0d00
            dum3(ikb,6)=0.0d00
            dum4(ikb,1)=0.0d00
            dum4(ikb,2)=0.0d00
            dum4(ikb,3)=0.0d00
            dum4(ikb,4)=0.0d00
            dum4(ikb,5)=0.0d00
            dum4(ikb,6)=0.0d00
            dum5(ikb,1)=0.0d00
            dum5(ikb,2)=0.0d00
            dum5(ikb,3)=0.0d00
            dum5(ikb,4)=0.0d00
            dum5(ikb,5)=0.0d00
            dum5(ikb,6)=0.0d00
            dum6(ikb,1)=0.0d00
            dum6(ikb,2)=0.0d00
            dum6(ikb,3)=0.0d00
            dum6(ikb,4)=0.0d00
            dum6(ikb,5)=0.0d00
            dum6(ikb,6)=0.0d00
          enddo
          do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            dum1(kb,1)=a(k-nelp1+na11)
            dum1(kb,2)=a(k-nelp1+na12)
            dum1(kb,3)=a(k-nelp1+na13)
            dum1(kb,4)=a(k-nelp1+na14)
            dum1(kb,5)=a(k-nelp1+na15)
            dum1(kb,6)=a(k-nelp1+na16)
            dum2(kb,1)=a(k-nelp1+na21)
            dum2(kb,2)=a(k-nelp1+na22)
            dum2(kb,3)=a(k-nelp1+na23)
            dum2(kb,4)=a(k-nelp1+na24)
            dum2(kb,5)=a(k-nelp1+na25)
            dum2(kb,6)=a(k-nelp1+na26)
            dum3(kb,1)=a(k-nelp1+na31)
            dum3(kb,2)=a(k-nelp1+na32)
            dum3(kb,3)=a(k-nelp1+na33)
            dum3(kb,4)=a(k-nelp1+na34)
            dum3(kb,5)=a(k-nelp1+na35)
            dum3(kb,6)=a(k-nelp1+na36)
            dum4(kb,1)=a(k-nelp1+na41)
            dum4(kb,2)=a(k-nelp1+na42)
            dum4(kb,3)=a(k-nelp1+na43)
            dum4(kb,4)=a(k-nelp1+na44)
            dum4(kb,5)=a(k-nelp1+na45)
            dum4(kb,6)=a(k-nelp1+na46)
            dum5(kb,1)=a(k-nelp1+na51)
            dum5(kb,2)=a(k-nelp1+na52)
            dum5(kb,3)=a(k-nelp1+na53)
            dum5(kb,4)=a(k-nelp1+na54)
            dum5(kb,5)=a(k-nelp1+na55)
            dum5(kb,6)=a(k-nelp1+na56)
            dum6(kb,1)=a(k-nelp1+na61)
            dum6(kb,2)=a(k-nelp1+na62)
            dum6(kb,3)=a(k-nelp1+na63)
            dum6(kb,4)=a(k-nelp1+na64)
            dum6(kb,5)=a(k-nelp1+na65)
            dum6(kb,6)=a(k-nelp1+na66)
          enddo
          do l=i1,i2
            kb=nop(l)
            b(l-nelp1+nb11)=dum1(kb,1)
            b(l-nelp1+nb12)=dum1(kb,2)
            b(l-nelp1+nb13)=dum1(kb,3)
            b(l-nelp1+nb14)=dum1(kb,4)
            b(l-nelp1+nb15)=dum1(kb,5)
            b(l-nelp1+nb16)=dum1(kb,6)
            b(l-nelp1+nb21)=dum2(kb,1)
            b(l-nelp1+nb22)=dum2(kb,2)
            b(l-nelp1+nb23)=dum2(kb,3)
            b(l-nelp1+nb24)=dum2(kb,4)
            b(l-nelp1+nb25)=dum2(kb,5)
            b(l-nelp1+nb26)=dum2(kb,6)
            b(l-nelp1+nb31)=dum3(kb,1)
            b(l-nelp1+nb32)=dum3(kb,2)
            b(l-nelp1+nb33)=dum3(kb,3)
            b(l-nelp1+nb34)=dum3(kb,4)
            b(l-nelp1+nb35)=dum3(kb,5)
            b(l-nelp1+nb36)=dum3(kb,6)
            b(l-nelp1+nb41)=dum4(kb,1)
            b(l-nelp1+nb42)=dum4(kb,2)
            b(l-nelp1+nb43)=dum4(kb,3)
            b(l-nelp1+nb44)=dum4(kb,4)
            b(l-nelp1+nb45)=dum4(kb,5)
            b(l-nelp1+nb46)=dum4(kb,6)
            b(l-nelp1+nb51)=dum5(kb,1)
            b(l-nelp1+nb52)=dum5(kb,2)
            b(l-nelp1+nb53)=dum5(kb,3)
            b(l-nelp1+nb54)=dum5(kb,4)
            b(l-nelp1+nb55)=dum5(kb,5)
            b(l-nelp1+nb56)=dum5(kb,6)
            b(l-nelp1+nb61)=dum6(kb,1)
            b(l-nelp1+nb62)=dum6(kb,2)
            b(l-nelp1+nb63)=dum6(kb,3)
            b(l-nelp1+nb64)=dum6(kb,4)
            b(l-nelp1+nb65)=dum6(kb,5)
            b(l-nelp1+nb66)=dum6(kb,6)
          enddo
      enddo
c     8 december 1990.
c     partial lu factorization of a (using storage of b)
c     assumptions
c     1 a has identity blocks on the diagonal (?? makes sense for
c       >1 equations per block [decoupling] but maybe not for 1)
c     basic algorithm
c     for i=1 to n
c       for j=1 to i
c         for k=1 to j-1
c           l(i,j)=l(i,j)-l(i,k)*u(k,j)
c       for j=i+1 to n
c         for k=1 to i-1
c           u(i,j)=u(i,j)-l(i,k)*u(k,j)
c         u(i,j)=u(i,j)/l(i,i)
c     calculation order - this is the way it was done in the original
c       routines
c       - more natural for sparse storage method
c       - faster on the one test problem that i tried
c     1 1 1 ...    1
c     2 3 3 ...    3
c     4 4 5 ...    5
c     .
c     .
c     .
c
c     for i=1 to n
      im1nq2=0
      do 360 i=1,nel
        ista=nop(i)+1
        ipvt=npvt(i)
        ipvtp1=ipvt+1
        ipvtm1=ipvt-1
        ifin=nop(i+1)
c       for j=1 to i
        do 270 ij=ista,ipvt
          j=nop(ij)
          ijind=ij-nelp1
c         for k=1 to j-1
          do 260 ik=ista,ij-1
            k=nop(ik)
            kj=nop(k+1)
            j2=nop(kj)
            kpvt=npvt(k)
  210       if ((j2.le.j).or.(kj.le.kpvt)) goto 220
              kj=kj-1
              j2=nop(kj)
              goto 210
  220       continue
            if (j2.eq.j) then
c             l(i,j)=l(i,j)-l(i,k)*u(k,j)
              ikind=ik-nelp1
              kjind=kj-nelp1
              b(ijind+nb11)=b(ijind+nb11)
     &          -b(ikind+nb11)*b(kjind+nb11)
     &          -b(ikind+nb12)*b(kjind+nb21)
     &          -b(ikind+nb13)*b(kjind+nb31)
     &          -b(ikind+nb14)*b(kjind+nb41)
     &          -b(ikind+nb15)*b(kjind+nb51)
     &          -b(ikind+nb16)*b(kjind+nb61)
              b(ijind+nb12)=b(ijind+nb12)
     &          -b(ikind+nb11)*b(kjind+nb12)
     &          -b(ikind+nb12)*b(kjind+nb22)
     &          -b(ikind+nb13)*b(kjind+nb32)
     &          -b(ikind+nb14)*b(kjind+nb42)
     &          -b(ikind+nb15)*b(kjind+nb52)
     &          -b(ikind+nb16)*b(kjind+nb62)
              b(ijind+nb13)=b(ijind+nb13)
     &          -b(ikind+nb11)*b(kjind+nb13)
     &          -b(ikind+nb12)*b(kjind+nb23)
     &          -b(ikind+nb13)*b(kjind+nb33)
     &          -b(ikind+nb14)*b(kjind+nb43)
     &          -b(ikind+nb15)*b(kjind+nb53)
     &          -b(ikind+nb16)*b(kjind+nb63)
              b(ijind+nb14)=b(ijind+nb14)
     &          -b(ikind+nb11)*b(kjind+nb14)
     &          -b(ikind+nb12)*b(kjind+nb24)
     &          -b(ikind+nb13)*b(kjind+nb34)
     &          -b(ikind+nb14)*b(kjind+nb44)
     &          -b(ikind+nb15)*b(kjind+nb54)
     &          -b(ikind+nb16)*b(kjind+nb64)
              b(ijind+nb15)=b(ijind+nb15)
     &          -b(ikind+nb11)*b(kjind+nb15)
     &          -b(ikind+nb12)*b(kjind+nb25)
     &          -b(ikind+nb13)*b(kjind+nb35)
     &          -b(ikind+nb14)*b(kjind+nb45)
     &          -b(ikind+nb15)*b(kjind+nb55)
     &          -b(ikind+nb16)*b(kjind+nb65)
              b(ijind+nb16)=b(ijind+nb16)
     &          -b(ikind+nb11)*b(kjind+nb16)
     &          -b(ikind+nb12)*b(kjind+nb26)
     &          -b(ikind+nb13)*b(kjind+nb36)
     &          -b(ikind+nb14)*b(kjind+nb46)
     &          -b(ikind+nb15)*b(kjind+nb56)
     &          -b(ikind+nb16)*b(kjind+nb66)




              b(ijind+nb21)=b(ijind+nb21)
     &          -b(ikind+nb21)*b(kjind+nb11)
     &          -b(ikind+nb22)*b(kjind+nb21)
     &          -b(ikind+nb23)*b(kjind+nb31)
     &          -b(ikind+nb24)*b(kjind+nb41)
     &          -b(ikind+nb25)*b(kjind+nb51)
     &          -b(ikind+nb26)*b(kjind+nb61)
              b(ijind+nb22)=b(ijind+nb22)
     &          -b(ikind+nb21)*b(kjind+nb12)
     &          -b(ikind+nb22)*b(kjind+nb22)
     &          -b(ikind+nb23)*b(kjind+nb32)
     &          -b(ikind+nb24)*b(kjind+nb42)
     &          -b(ikind+nb25)*b(kjind+nb52)
     &          -b(ikind+nb26)*b(kjind+nb62)
              b(ijind+nb23)=b(ijind+nb23)
     &          -b(ikind+nb21)*b(kjind+nb13)
     &          -b(ikind+nb22)*b(kjind+nb23)
     &          -b(ikind+nb23)*b(kjind+nb33)
     &          -b(ikind+nb24)*b(kjind+nb43)
     &          -b(ikind+nb25)*b(kjind+nb53)
     &          -b(ikind+nb26)*b(kjind+nb63)
              b(ijind+nb24)=b(ijind+nb24)
     &          -b(ikind+nb21)*b(kjind+nb14)
     &          -b(ikind+nb22)*b(kjind+nb24)
     &          -b(ikind+nb23)*b(kjind+nb34)
     &          -b(ikind+nb24)*b(kjind+nb44)
     &          -b(ikind+nb25)*b(kjind+nb54)
     &          -b(ikind+nb26)*b(kjind+nb64)
              b(ijind+nb25)=b(ijind+nb25)
     &          -b(ikind+nb21)*b(kjind+nb15)
     &          -b(ikind+nb22)*b(kjind+nb25)
     &          -b(ikind+nb23)*b(kjind+nb35)
     &          -b(ikind+nb24)*b(kjind+nb45)
     &          -b(ikind+nb25)*b(kjind+nb55)
     &          -b(ikind+nb26)*b(kjind+nb65)
              b(ijind+nb26)=b(ijind+nb26)
     &          -b(ikind+nb21)*b(kjind+nb16)
     &          -b(ikind+nb22)*b(kjind+nb26)
     &          -b(ikind+nb23)*b(kjind+nb36)
     &          -b(ikind+nb24)*b(kjind+nb46)
     &          -b(ikind+nb25)*b(kjind+nb56)
     &          -b(ikind+nb26)*b(kjind+nb66)




              b(ijind+nb31)=b(ijind+nb31)
     &          -b(ikind+nb31)*b(kjind+nb11)
     &          -b(ikind+nb32)*b(kjind+nb21)
     &          -b(ikind+nb33)*b(kjind+nb31)
     &          -b(ikind+nb34)*b(kjind+nb41)
     &          -b(ikind+nb35)*b(kjind+nb51)
     &          -b(ikind+nb36)*b(kjind+nb61)
              b(ijind+nb32)=b(ijind+nb32)
     &          -b(ikind+nb31)*b(kjind+nb12)
     &          -b(ikind+nb32)*b(kjind+nb22)
     &          -b(ikind+nb33)*b(kjind+nb32)
     &          -b(ikind+nb34)*b(kjind+nb42)
     &          -b(ikind+nb35)*b(kjind+nb52)
     &          -b(ikind+nb36)*b(kjind+nb62)
              b(ijind+nb33)=b(ijind+nb33)
     &          -b(ikind+nb31)*b(kjind+nb13)
     &          -b(ikind+nb32)*b(kjind+nb23)
     &          -b(ikind+nb33)*b(kjind+nb33)
     &          -b(ikind+nb34)*b(kjind+nb43)
     &          -b(ikind+nb35)*b(kjind+nb53)
     &          -b(ikind+nb36)*b(kjind+nb63)
              b(ijind+nb34)=b(ijind+nb34)
     &          -b(ikind+nb31)*b(kjind+nb14)
     &          -b(ikind+nb32)*b(kjind+nb24)
     &          -b(ikind+nb33)*b(kjind+nb34)
     &          -b(ikind+nb34)*b(kjind+nb44)
     &          -b(ikind+nb35)*b(kjind+nb54)
     &          -b(ikind+nb36)*b(kjind+nb64)
              b(ijind+nb35)=b(ijind+nb35)
     &          -b(ikind+nb31)*b(kjind+nb15)
     &          -b(ikind+nb32)*b(kjind+nb25)
     &          -b(ikind+nb33)*b(kjind+nb35)
     &          -b(ikind+nb34)*b(kjind+nb45)
     &          -b(ikind+nb35)*b(kjind+nb55)
     &          -b(ikind+nb36)*b(kjind+nb65)
              b(ijind+nb36)=b(ijind+nb36)
     &          -b(ikind+nb31)*b(kjind+nb16)
     &          -b(ikind+nb32)*b(kjind+nb26)
     &          -b(ikind+nb33)*b(kjind+nb36)
     &          -b(ikind+nb34)*b(kjind+nb46)
     &          -b(ikind+nb35)*b(kjind+nb56)
     &          -b(ikind+nb36)*b(kjind+nb66)







              b(ijind+nb41)=b(ijind+nb41)
     &          -b(ikind+nb41)*b(kjind+nb11)
     &          -b(ikind+nb42)*b(kjind+nb21)
     &          -b(ikind+nb43)*b(kjind+nb31)
     &          -b(ikind+nb44)*b(kjind+nb41)
     &          -b(ikind+nb45)*b(kjind+nb51)
     &          -b(ikind+nb46)*b(kjind+nb61)
              b(ijind+nb42)=b(ijind+nb42)
     &          -b(ikind+nb41)*b(kjind+nb12)
     &          -b(ikind+nb42)*b(kjind+nb22)
     &          -b(ikind+nb43)*b(kjind+nb32)
     &          -b(ikind+nb44)*b(kjind+nb42)
     &          -b(ikind+nb45)*b(kjind+nb52)
     &          -b(ikind+nb46)*b(kjind+nb62)
              b(ijind+nb43)=b(ijind+nb43)
     &          -b(ikind+nb41)*b(kjind+nb13)
     &          -b(ikind+nb42)*b(kjind+nb23)
     &          -b(ikind+nb43)*b(kjind+nb33)
     &          -b(ikind+nb44)*b(kjind+nb43)
     &          -b(ikind+nb45)*b(kjind+nb53)
     &          -b(ikind+nb46)*b(kjind+nb63)
              b(ijind+nb44)=b(ijind+nb44)
     &          -b(ikind+nb41)*b(kjind+nb14)
     &          -b(ikind+nb42)*b(kjind+nb24)
     &          -b(ikind+nb43)*b(kjind+nb34)
     &          -b(ikind+nb44)*b(kjind+nb44)
     &          -b(ikind+nb45)*b(kjind+nb54)
     &          -b(ikind+nb46)*b(kjind+nb64)
              b(ijind+nb45)=b(ijind+nb45)
     &          -b(ikind+nb41)*b(kjind+nb15)
     &          -b(ikind+nb42)*b(kjind+nb25)
     &          -b(ikind+nb43)*b(kjind+nb35)
     &          -b(ikind+nb44)*b(kjind+nb45)
     &          -b(ikind+nb45)*b(kjind+nb55)
     &          -b(ikind+nb46)*b(kjind+nb65)
              b(ijind+nb46)=b(ijind+nb46)
     &          -b(ikind+nb41)*b(kjind+nb16)
     &          -b(ikind+nb42)*b(kjind+nb26)
     &          -b(ikind+nb43)*b(kjind+nb36)
     &          -b(ikind+nb44)*b(kjind+nb46)
     &          -b(ikind+nb45)*b(kjind+nb56)
     &          -b(ikind+nb46)*b(kjind+nb66)







              b(ijind+nb51)=b(ijind+nb51)
     &          -b(ikind+nb51)*b(kjind+nb11)
     &          -b(ikind+nb52)*b(kjind+nb21)
     &          -b(ikind+nb53)*b(kjind+nb31)
     &          -b(ikind+nb54)*b(kjind+nb41)
     &          -b(ikind+nb55)*b(kjind+nb51)
     &          -b(ikind+nb56)*b(kjind+nb61)
              b(ijind+nb52)=b(ijind+nb52)
     &          -b(ikind+nb51)*b(kjind+nb12)
     &          -b(ikind+nb52)*b(kjind+nb22)
     &          -b(ikind+nb53)*b(kjind+nb32)
     &          -b(ikind+nb54)*b(kjind+nb42)
     &          -b(ikind+nb55)*b(kjind+nb52)
     &          -b(ikind+nb56)*b(kjind+nb62)
              b(ijind+nb53)=b(ijind+nb53)
     &          -b(ikind+nb51)*b(kjind+nb13)
     &          -b(ikind+nb52)*b(kjind+nb23)
     &          -b(ikind+nb53)*b(kjind+nb33)
     &          -b(ikind+nb54)*b(kjind+nb43)
     &          -b(ikind+nb55)*b(kjind+nb53)
     &          -b(ikind+nb56)*b(kjind+nb63)
              b(ijind+nb54)=b(ijind+nb54)
     &          -b(ikind+nb51)*b(kjind+nb14)
     &          -b(ikind+nb52)*b(kjind+nb24)
     &          -b(ikind+nb53)*b(kjind+nb34)
     &          -b(ikind+nb54)*b(kjind+nb44)
     &          -b(ikind+nb55)*b(kjind+nb54)
     &          -b(ikind+nb56)*b(kjind+nb64)
              b(ijind+nb55)=b(ijind+nb55)
     &          -b(ikind+nb51)*b(kjind+nb15)
     &          -b(ikind+nb52)*b(kjind+nb25)
     &          -b(ikind+nb53)*b(kjind+nb35)
     &          -b(ikind+nb54)*b(kjind+nb45)
     &          -b(ikind+nb55)*b(kjind+nb55)
     &          -b(ikind+nb56)*b(kjind+nb65)
              b(ijind+nb56)=b(ijind+nb56)
     &          -b(ikind+nb51)*b(kjind+nb16)
     &          -b(ikind+nb52)*b(kjind+nb26)
     &          -b(ikind+nb53)*b(kjind+nb36)
     &          -b(ikind+nb54)*b(kjind+nb46)
     &          -b(ikind+nb55)*b(kjind+nb56)
     &          -b(ikind+nb56)*b(kjind+nb66)






              b(ijind+nb61)=b(ijind+nb61)
     &          -b(ikind+nb61)*b(kjind+nb11)
     &          -b(ikind+nb62)*b(kjind+nb21)
     &          -b(ikind+nb63)*b(kjind+nb31)
     &          -b(ikind+nb64)*b(kjind+nb41)
     &          -b(ikind+nb65)*b(kjind+nb51)
     &          -b(ikind+nb66)*b(kjind+nb61)
              b(ijind+nb62)=b(ijind+nb62)
     &          -b(ikind+nb61)*b(kjind+nb12)
     &          -b(ikind+nb62)*b(kjind+nb22)
     &          -b(ikind+nb63)*b(kjind+nb32)
     &          -b(ikind+nb64)*b(kjind+nb42)
     &          -b(ikind+nb65)*b(kjind+nb52)
     &          -b(ikind+nb66)*b(kjind+nb62)
              b(ijind+nb63)=b(ijind+nb63)
     &          -b(ikind+nb61)*b(kjind+nb13)
     &          -b(ikind+nb62)*b(kjind+nb23)
     &          -b(ikind+nb63)*b(kjind+nb33)
     &          -b(ikind+nb64)*b(kjind+nb43)
     &          -b(ikind+nb65)*b(kjind+nb53)
     &          -b(ikind+nb66)*b(kjind+nb63)
              b(ijind+nb64)=b(ijind+nb64)
     &          -b(ikind+nb61)*b(kjind+nb14)
     &          -b(ikind+nb62)*b(kjind+nb24)
     &          -b(ikind+nb63)*b(kjind+nb34)
     &          -b(ikind+nb64)*b(kjind+nb44)
     &          -b(ikind+nb65)*b(kjind+nb54)
     &          -b(ikind+nb66)*b(kjind+nb64)
              b(ijind+nb65)=b(ijind+nb65)
     &          -b(ikind+nb61)*b(kjind+nb15)
     &          -b(ikind+nb62)*b(kjind+nb25)
     &          -b(ikind+nb63)*b(kjind+nb35)
     &          -b(ikind+nb64)*b(kjind+nb45)
     &          -b(ikind+nb65)*b(kjind+nb55)
     &          -b(ikind+nb66)*b(kjind+nb65)
              b(ijind+nb66)=b(ijind+nb66)
     &          -b(ikind+nb61)*b(kjind+nb16)
     &          -b(ikind+nb62)*b(kjind+nb26)
     &          -b(ikind+nb63)*b(kjind+nb36)
     &          -b(ikind+nb64)*b(kjind+nb46)
     &          -b(ikind+nb65)*b(kjind+nb56)
     &          -b(ikind+nb66)*b(kjind+nb66)



            endif
  260     continue
  270   continue
c       calculate the inverse of l(i,i)
        ipvind=ipvt-nelp1
        dd11=b(ipvind+nb11)
        dd12=b(ipvind+nb12)
        dd13=b(ipvind+nb13)
        dd14=b(ipvind+nb14)
        dd15=b(ipvind+nb15)
        dd16=b(ipvind+nb16)
        dd21=b(ipvind+nb21)
        dd22=b(ipvind+nb22)
        dd23=b(ipvind+nb23)
        dd24=b(ipvind+nb24)
        dd25=b(ipvind+nb25)
        dd26=b(ipvind+nb26)
        dd31=b(ipvind+nb31)
        dd32=b(ipvind+nb32)
        dd33=b(ipvind+nb33)
        dd34=b(ipvind+nb34)
        dd35=b(ipvind+nb35)
        dd36=b(ipvind+nb36)
        dd41=b(ipvind+nb41)
        dd42=b(ipvind+nb42)
        dd43=b(ipvind+nb43)
        dd44=b(ipvind+nb44)
        dd45=b(ipvind+nb45)
        dd46=b(ipvind+nb46)
        dd51=b(ipvind+nb51)
        dd52=b(ipvind+nb52)
        dd53=b(ipvind+nb53)
        dd54=b(ipvind+nb54)
        dd55=b(ipvind+nb55)
        dd56=b(ipvind+nb56)
        dd61=b(ipvind+nb61)
        dd62=b(ipvind+nb62)
        dd63=b(ipvind+nb63)
        dd64=b(ipvind+nb64)
        dd65=b(ipvind+nb65)
        dd66=b(ipvind+nb66)





        dd21=dd21/dd11
        dd31=dd31/dd11
        dd41=dd41/dd11
        dd51=dd51/dd11
        dd61=dd61/dd11
        dd22=dd22
     &    -dd21*dd12
        dd32=dd32
     &    -dd31*dd12
        dd42=dd42
     &    -dd41*dd12
        dd52=dd52
     &    -dd51*dd12
        dd62=dd62
     &    -dd61*dd12


        dd32=dd32/dd22
        dd42=dd42/dd22
        dd52=dd52/dd22
        dd62=dd62/dd22


        dd23=dd23
     &    -dd21*dd13
        dd33=dd33
     &    -dd31*dd13
     &    -dd32*dd23
        dd43=dd43
     &    -dd41*dd13
     &    -dd42*dd23
        dd53=dd53
     &    -dd51*dd13
     &    -dd52*dd23
        dd63=dd63
     &    -dd61*dd13
     &    -dd62*dd23
        dd43=dd43/dd33
        dd53=dd53/dd33
        dd63=dd63/dd33

        dd24=dd24
     &    -dd21*dd14
        dd34=dd34
     &    -dd31*dd14
     &    -dd32*dd24
        dd44=dd44
     &    -dd41*dd14
     &    -dd42*dd24
     &    -dd43*dd34
        dd54=dd54
     &    -dd51*dd14
     &    -dd52*dd24
     &    -dd53*dd34
        dd64=dd64
     &    -dd61*dd14
     &    -dd62*dd24
     &    -dd63*dd34
        dd54=dd54/dd44
        dd64=dd64/dd44
        

        dd25=dd25
     &    -dd21*dd15
        dd35=dd35
     &    -dd31*dd15
     &    -dd32*dd25
        dd45=dd45
     &    -dd41*dd15
     &    -dd42*dd25
     &    -dd43*dd35
        dd55=dd55
     &    -dd51*dd15
     &    -dd52*dd25
     &    -dd53*dd35
     &    -dd54*dd45
        dd65=dd65
     &    -dd61*dd15
     &    -dd62*dd25
     &    -dd63*dd35
     &    -dd64*dd45
        dd65=dd65/dd55

        dd26=dd26
     &    -dd21*dd16
        dd36=dd36
     &    -dd31*dd16
     &    -dd32*dd26
        dd46=dd46
     &    -dd41*dd16
     &    -dd42*dd26
     &    -dd43*dd36
        dd56=dd56
     &    -dd51*dd16
     &    -dd52*dd26
     &    -dd53*dd36
     &    -dd54*dd46
        dd66=dd66
     &    -dd61*dd16
     &    -dd62*dd26
     &    -dd63*dd36
     &    -dd64*dd46
     &    -dd65*dd56



        piv(im1nq2+ 1)=dd11
        piv(im1nq2+ 2)=dd12
        piv(im1nq2+ 3)=dd13
        piv(im1nq2+ 4)=dd14
        piv(im1nq2+ 5)=dd15
        piv(im1nq2+ 6)=dd16
        piv(im1nq2+ 7)=dd21
        piv(im1nq2+ 8)=dd22
        piv(im1nq2+ 9)=dd23
        piv(im1nq2+10)=dd24
        piv(im1nq2+11)=dd25
        piv(im1nq2+12)=dd26
        piv(im1nq2+13)=dd31
        piv(im1nq2+14)=dd32
        piv(im1nq2+15)=dd33
        piv(im1nq2+16)=dd34
        piv(im1nq2+17)=dd35
        piv(im1nq2+18)=dd36
        piv(im1nq2+19)=dd41
        piv(im1nq2+20)=dd42
        piv(im1nq2+21)=dd43
        piv(im1nq2+22)=dd44
        piv(im1nq2+23)=dd45
        piv(im1nq2+24)=dd46
        piv(im1nq2+25)=dd51
        piv(im1nq2+26)=dd52
        piv(im1nq2+27)=dd53
        piv(im1nq2+28)=dd54
        piv(im1nq2+29)=dd55
        piv(im1nq2+30)=dd56
        piv(im1nq2+31)=dd61
        piv(im1nq2+32)=dd62
        piv(im1nq2+33)=dd63
        piv(im1nq2+34)=dd64
        piv(im1nq2+35)=dd65
        piv(im1nq2+36)=dd66
c       for j=i+1 to n
        do 350 ij=ipvtp1,ifin
          ijind=ij-nelp1
          j=nop(ij)
c         for k=1 to i-1
          do 340 ik=ista,ipvtm1
            k=nop(ik)
            kj=nop(k+1)
            j2=nop(kj)
            kpvt=npvt(k)
  290       if ((j2.le.j).or.(kj.le.kpvt)) goto 300
              kj=kj-1
              j2=nop(kj)
              goto 290
  300       continue
            if (j2.eq.j) then
c             u(i,j)=u(i,j)-l(i,k)*u(k,j)
              ikind=ik-nelp1
              kjind=kj-nelp1
              b(ijind+nb11)=b(ijind+nb11)
     &          -b(ikind+nb11)*b(kjind+nb11)
     &          -b(ikind+nb12)*b(kjind+nb21)
     &          -b(ikind+nb13)*b(kjind+nb31)
     &          -b(ikind+nb14)*b(kjind+nb41)
     &          -b(ikind+nb15)*b(kjind+nb51)
     &          -b(ikind+nb16)*b(kjind+nb61)
              b(ijind+nb12)=b(ijind+nb12)
     &          -b(ikind+nb11)*b(kjind+nb12)
     &          -b(ikind+nb12)*b(kjind+nb22)
     &          -b(ikind+nb13)*b(kjind+nb32)
     &          -b(ikind+nb14)*b(kjind+nb42)
     &          -b(ikind+nb15)*b(kjind+nb52)
     &          -b(ikind+nb16)*b(kjind+nb62)
              b(ijind+nb13)=b(ijind+nb13)
     &          -b(ikind+nb11)*b(kjind+nb13)
     &          -b(ikind+nb12)*b(kjind+nb23)
     &          -b(ikind+nb13)*b(kjind+nb33)
     &          -b(ikind+nb14)*b(kjind+nb43)
     &          -b(ikind+nb15)*b(kjind+nb53)
     &          -b(ikind+nb16)*b(kjind+nb63)
              b(ijind+nb14)=b(ijind+nb14)
     &          -b(ikind+nb11)*b(kjind+nb14)
     &          -b(ikind+nb12)*b(kjind+nb24)
     &          -b(ikind+nb13)*b(kjind+nb34)
     &          -b(ikind+nb14)*b(kjind+nb44)
     &          -b(ikind+nb15)*b(kjind+nb54)
     &          -b(ikind+nb16)*b(kjind+nb64)
              b(ijind+nb15)=b(ijind+nb15)
     &          -b(ikind+nb11)*b(kjind+nb15)
     &          -b(ikind+nb12)*b(kjind+nb25)
     &          -b(ikind+nb13)*b(kjind+nb35)
     &          -b(ikind+nb14)*b(kjind+nb45)
     &          -b(ikind+nb15)*b(kjind+nb55)
     &          -b(ikind+nb16)*b(kjind+nb65)
              b(ijind+nb16)=b(ijind+nb16)
     &          -b(ikind+nb11)*b(kjind+nb16)
     &          -b(ikind+nb12)*b(kjind+nb26)
     &          -b(ikind+nb13)*b(kjind+nb36)
     &          -b(ikind+nb14)*b(kjind+nb46)
     &          -b(ikind+nb15)*b(kjind+nb56)
     &          -b(ikind+nb16)*b(kjind+nb66)


              b(ijind+nb21)=b(ijind+nb21)
     &          -b(ikind+nb21)*b(kjind+nb11)
     &          -b(ikind+nb22)*b(kjind+nb21)
     &          -b(ikind+nb23)*b(kjind+nb31)
     &          -b(ikind+nb24)*b(kjind+nb41)
     &          -b(ikind+nb25)*b(kjind+nb51)
     &          -b(ikind+nb26)*b(kjind+nb61)
              b(ijind+nb22)=b(ijind+nb22)
     &          -b(ikind+nb21)*b(kjind+nb12)
     &          -b(ikind+nb22)*b(kjind+nb22)
     &          -b(ikind+nb23)*b(kjind+nb32)
     &          -b(ikind+nb24)*b(kjind+nb42)
     &          -b(ikind+nb25)*b(kjind+nb52)
     &          -b(ikind+nb26)*b(kjind+nb62)
              b(ijind+nb23)=b(ijind+nb23)
     &          -b(ikind+nb21)*b(kjind+nb13)
     &          -b(ikind+nb22)*b(kjind+nb23)
     &          -b(ikind+nb23)*b(kjind+nb33)
     &          -b(ikind+nb24)*b(kjind+nb43)
     &          -b(ikind+nb25)*b(kjind+nb53)
     &          -b(ikind+nb26)*b(kjind+nb63)
              b(ijind+nb24)=b(ijind+nb24)
     &          -b(ikind+nb21)*b(kjind+nb14)
     &          -b(ikind+nb22)*b(kjind+nb24)
     &          -b(ikind+nb23)*b(kjind+nb34)
     &          -b(ikind+nb24)*b(kjind+nb44)
     &          -b(ikind+nb25)*b(kjind+nb54)
     &          -b(ikind+nb26)*b(kjind+nb64)
              b(ijind+nb25)=b(ijind+nb25)
     &          -b(ikind+nb21)*b(kjind+nb15)
     &          -b(ikind+nb22)*b(kjind+nb25)
     &          -b(ikind+nb23)*b(kjind+nb35)
     &          -b(ikind+nb24)*b(kjind+nb45)
     &          -b(ikind+nb25)*b(kjind+nb55)
     &          -b(ikind+nb26)*b(kjind+nb65)
              b(ijind+nb26)=b(ijind+nb26)
     &          -b(ikind+nb21)*b(kjind+nb16)
     &          -b(ikind+nb22)*b(kjind+nb26)
     &          -b(ikind+nb23)*b(kjind+nb36)
     &          -b(ikind+nb24)*b(kjind+nb46)
     &          -b(ikind+nb25)*b(kjind+nb56)
     &          -b(ikind+nb26)*b(kjind+nb66)



              b(ijind+nb31)=b(ijind+nb31)
     &          -b(ikind+nb31)*b(kjind+nb11)
     &          -b(ikind+nb32)*b(kjind+nb21)
     &          -b(ikind+nb33)*b(kjind+nb31)
     &          -b(ikind+nb34)*b(kjind+nb41)
     &          -b(ikind+nb35)*b(kjind+nb51)
     &          -b(ikind+nb36)*b(kjind+nb61)
              b(ijind+nb32)=b(ijind+nb32)
     &          -b(ikind+nb31)*b(kjind+nb12)
     &          -b(ikind+nb32)*b(kjind+nb22)
     &          -b(ikind+nb33)*b(kjind+nb32)
     &          -b(ikind+nb34)*b(kjind+nb42)
     &          -b(ikind+nb35)*b(kjind+nb52)
     &          -b(ikind+nb36)*b(kjind+nb62)
              b(ijind+nb33)=b(ijind+nb33)
     &          -b(ikind+nb31)*b(kjind+nb13)
     &          -b(ikind+nb32)*b(kjind+nb23)
     &          -b(ikind+nb33)*b(kjind+nb33)
     &          -b(ikind+nb34)*b(kjind+nb43)
     &          -b(ikind+nb35)*b(kjind+nb53)
     &          -b(ikind+nb36)*b(kjind+nb63)
              b(ijind+nb34)=b(ijind+nb34)
     &          -b(ikind+nb31)*b(kjind+nb14)
     &          -b(ikind+nb32)*b(kjind+nb24)
     &          -b(ikind+nb33)*b(kjind+nb34)
     &          -b(ikind+nb34)*b(kjind+nb44)
     &          -b(ikind+nb35)*b(kjind+nb54)
     &          -b(ikind+nb36)*b(kjind+nb64)
              b(ijind+nb35)=b(ijind+nb35)
     &          -b(ikind+nb31)*b(kjind+nb15)
     &          -b(ikind+nb32)*b(kjind+nb25)
     &          -b(ikind+nb33)*b(kjind+nb35)
     &          -b(ikind+nb34)*b(kjind+nb45)
     &          -b(ikind+nb35)*b(kjind+nb55)
     &          -b(ikind+nb36)*b(kjind+nb65)
              b(ijind+nb36)=b(ijind+nb36)
     &          -b(ikind+nb31)*b(kjind+nb16)
     &          -b(ikind+nb32)*b(kjind+nb26)
     &          -b(ikind+nb33)*b(kjind+nb36)
     &          -b(ikind+nb34)*b(kjind+nb46)
     &          -b(ikind+nb35)*b(kjind+nb56)
     &          -b(ikind+nb36)*b(kjind+nb66)




              b(ijind+nb41)=b(ijind+nb41)
     &          -b(ikind+nb41)*b(kjind+nb11)
     &          -b(ikind+nb42)*b(kjind+nb21)
     &          -b(ikind+nb43)*b(kjind+nb31)
     &          -b(ikind+nb44)*b(kjind+nb41)
     &          -b(ikind+nb45)*b(kjind+nb51)
     &          -b(ikind+nb46)*b(kjind+nb61)
              b(ijind+nb42)=b(ijind+nb42)
     &          -b(ikind+nb41)*b(kjind+nb12)
     &          -b(ikind+nb42)*b(kjind+nb22)
     &          -b(ikind+nb43)*b(kjind+nb32)
     &          -b(ikind+nb44)*b(kjind+nb42)
     &          -b(ikind+nb45)*b(kjind+nb52)
     &          -b(ikind+nb46)*b(kjind+nb62)
              b(ijind+nb43)=b(ijind+nb43)
     &          -b(ikind+nb41)*b(kjind+nb13)
     &          -b(ikind+nb42)*b(kjind+nb23)
     &          -b(ikind+nb43)*b(kjind+nb33)
     &          -b(ikind+nb44)*b(kjind+nb43)
     &          -b(ikind+nb45)*b(kjind+nb53)
     &          -b(ikind+nb46)*b(kjind+nb63)
              b(ijind+nb44)=b(ijind+nb44)
     &          -b(ikind+nb41)*b(kjind+nb14)
     &          -b(ikind+nb42)*b(kjind+nb24)
     &          -b(ikind+nb43)*b(kjind+nb34)
     &          -b(ikind+nb44)*b(kjind+nb44)
     &          -b(ikind+nb45)*b(kjind+nb54)
     &          -b(ikind+nb46)*b(kjind+nb64)
              b(ijind+nb45)=b(ijind+nb45)
     &          -b(ikind+nb41)*b(kjind+nb15)
     &          -b(ikind+nb42)*b(kjind+nb25)
     &          -b(ikind+nb43)*b(kjind+nb35)
     &          -b(ikind+nb44)*b(kjind+nb45)
     &          -b(ikind+nb45)*b(kjind+nb55)
     &          -b(ikind+nb46)*b(kjind+nb65)
              b(ijind+nb46)=b(ijind+nb46)
     &          -b(ikind+nb41)*b(kjind+nb16)
     &          -b(ikind+nb42)*b(kjind+nb26)
     &          -b(ikind+nb43)*b(kjind+nb36)
     &          -b(ikind+nb44)*b(kjind+nb46)
     &          -b(ikind+nb45)*b(kjind+nb56)
     &          -b(ikind+nb46)*b(kjind+nb66)


              b(ijind+nb51)=b(ijind+nb51)
     &          -b(ikind+nb51)*b(kjind+nb11)
     &          -b(ikind+nb52)*b(kjind+nb21)
     &          -b(ikind+nb53)*b(kjind+nb31)
     &          -b(ikind+nb54)*b(kjind+nb41)
     &          -b(ikind+nb55)*b(kjind+nb51)
     &          -b(ikind+nb56)*b(kjind+nb61)
              b(ijind+nb52)=b(ijind+nb52)
     &          -b(ikind+nb51)*b(kjind+nb12)
     &          -b(ikind+nb52)*b(kjind+nb22)
     &          -b(ikind+nb53)*b(kjind+nb32)
     &          -b(ikind+nb54)*b(kjind+nb42)
     &          -b(ikind+nb55)*b(kjind+nb52)
     &          -b(ikind+nb56)*b(kjind+nb62)
              b(ijind+nb53)=b(ijind+nb53)
     &          -b(ikind+nb51)*b(kjind+nb13)
     &          -b(ikind+nb52)*b(kjind+nb23)
     &          -b(ikind+nb53)*b(kjind+nb33)
     &          -b(ikind+nb54)*b(kjind+nb43)
     &          -b(ikind+nb55)*b(kjind+nb53)
     &          -b(ikind+nb56)*b(kjind+nb63)
              b(ijind+nb54)=b(ijind+nb54)
     &          -b(ikind+nb51)*b(kjind+nb14)
     &          -b(ikind+nb52)*b(kjind+nb24)
     &          -b(ikind+nb53)*b(kjind+nb34)
     &          -b(ikind+nb54)*b(kjind+nb44)
     &          -b(ikind+nb55)*b(kjind+nb54)
     &          -b(ikind+nb56)*b(kjind+nb64)
              b(ijind+nb55)=b(ijind+nb55)
     &          -b(ikind+nb51)*b(kjind+nb15)
     &          -b(ikind+nb52)*b(kjind+nb25)
     &          -b(ikind+nb53)*b(kjind+nb35)
     &          -b(ikind+nb54)*b(kjind+nb45)
     &          -b(ikind+nb55)*b(kjind+nb55)
     &          -b(ikind+nb56)*b(kjind+nb65)
              b(ijind+nb56)=b(ijind+nb56)
     &          -b(ikind+nb51)*b(kjind+nb16)
     &          -b(ikind+nb52)*b(kjind+nb26)
     &          -b(ikind+nb53)*b(kjind+nb36)
     &          -b(ikind+nb54)*b(kjind+nb46)
     &          -b(ikind+nb55)*b(kjind+nb56)
     &          -b(ikind+nb56)*b(kjind+nb66)


              b(ijind+nb61)=b(ijind+nb61)
     &          -b(ikind+nb61)*b(kjind+nb11)
     &          -b(ikind+nb62)*b(kjind+nb21)
     &          -b(ikind+nb63)*b(kjind+nb31)
     &          -b(ikind+nb64)*b(kjind+nb41)
     &          -b(ikind+nb65)*b(kjind+nb51)
     &          -b(ikind+nb66)*b(kjind+nb61)
              b(ijind+nb62)=b(ijind+nb62)
     &          -b(ikind+nb61)*b(kjind+nb12)
     &          -b(ikind+nb62)*b(kjind+nb22)
     &          -b(ikind+nb63)*b(kjind+nb32)
     &          -b(ikind+nb64)*b(kjind+nb42)
     &          -b(ikind+nb65)*b(kjind+nb52)
     &          -b(ikind+nb66)*b(kjind+nb62)
              b(ijind+nb63)=b(ijind+nb63)
     &          -b(ikind+nb61)*b(kjind+nb13)
     &          -b(ikind+nb62)*b(kjind+nb23)
     &          -b(ikind+nb63)*b(kjind+nb33)
     &          -b(ikind+nb64)*b(kjind+nb43)
     &          -b(ikind+nb65)*b(kjind+nb53)
     &          -b(ikind+nb66)*b(kjind+nb63)
              b(ijind+nb64)=b(ijind+nb64)
     &          -b(ikind+nb61)*b(kjind+nb14)
     &          -b(ikind+nb62)*b(kjind+nb24)
     &          -b(ikind+nb63)*b(kjind+nb34)
     &          -b(ikind+nb64)*b(kjind+nb44)
     &          -b(ikind+nb65)*b(kjind+nb54)
     &          -b(ikind+nb66)*b(kjind+nb64)
              b(ijind+nb65)=b(ijind+nb65)
     &          -b(ikind+nb61)*b(kjind+nb15)
     &          -b(ikind+nb62)*b(kjind+nb25)
     &          -b(ikind+nb63)*b(kjind+nb35)
     &          -b(ikind+nb64)*b(kjind+nb45)
     &          -b(ikind+nb65)*b(kjind+nb55)
     &          -b(ikind+nb66)*b(kjind+nb65)
              b(ijind+nb66)=b(ijind+nb66)
     &          -b(ikind+nb61)*b(kjind+nb16)
     &          -b(ikind+nb62)*b(kjind+nb26)
     &          -b(ikind+nb63)*b(kjind+nb36)
     &          -b(ikind+nb64)*b(kjind+nb46)
     &          -b(ikind+nb65)*b(kjind+nb56)
     &          -b(ikind+nb66)*b(kjind+nb66)




            endif
  340     continue
c         u(i,j)=u(i,j)/l(i,i)
          bksb1=b(ijind+nb11)
          bksb2=b(ijind+nb21)
          bksb3=b(ijind+nb31)
          bksb4=b(ijind+nb41)
          bksb5=b(ijind+nb51)
          bksb6=b(ijind+nb61)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11



          b(ijind+nb11)=bksb1
          b(ijind+nb21)=bksb2
          b(ijind+nb31)=bksb3
          b(ijind+nb41)=bksb4
          b(ijind+nb51)=bksb5
          b(ijind+nb61)=bksb6
          bksb1=b(ijind+nb12)
          bksb2=b(ijind+nb22)
          bksb3=b(ijind+nb32)
          bksb4=b(ijind+nb42)
          bksb5=b(ijind+nb52)
          bksb6=b(ijind+nb62)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11



          b(ijind+nb12)=bksb1
          b(ijind+nb22)=bksb2
          b(ijind+nb32)=bksb3
          b(ijind+nb42)=bksb4
          b(ijind+nb52)=bksb5
          b(ijind+nb62)=bksb6
          bksb1=b(ijind+nb13)
          bksb2=b(ijind+nb23)
          bksb3=b(ijind+nb33)
          bksb4=b(ijind+nb43)
          bksb5=b(ijind+nb53)
          bksb6=b(ijind+nb63)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11




          b(ijind+nb13)=bksb1
          b(ijind+nb23)=bksb2
          b(ijind+nb33)=bksb3
          b(ijind+nb43)=bksb4
          b(ijind+nb53)=bksb5
          b(ijind+nb63)=bksb6
          bksb1=b(ijind+nb14)
          bksb2=b(ijind+nb24)
          bksb3=b(ijind+nb34)
          bksb4=b(ijind+nb44)
          bksb5=b(ijind+nb54)
          bksb6=b(ijind+nb64)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11




          b(ijind+nb14)=bksb1
          b(ijind+nb24)=bksb2
          b(ijind+nb34)=bksb3
          b(ijind+nb44)=bksb4
          b(ijind+nb54)=bksb5
          b(ijind+nb64)=bksb6
          bksb1=b(ijind+nb15)
          bksb2=b(ijind+nb25)
          bksb3=b(ijind+nb35)
          bksb4=b(ijind+nb45)
          bksb5=b(ijind+nb55)
          bksb6=b(ijind+nb65)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11



          b(ijind+nb15)=bksb1
          b(ijind+nb25)=bksb2
          b(ijind+nb35)=bksb3
          b(ijind+nb45)=bksb4
          b(ijind+nb55)=bksb5
          b(ijind+nb65)=bksb6
          bksb1=b(ijind+nb16)
          bksb2=b(ijind+nb26)
          bksb3=b(ijind+nb36)
          bksb4=b(ijind+nb46)
          bksb5=b(ijind+nb56)
          bksb6=b(ijind+nb66)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb5=bksb5
     &      -dd51*bksb1
     &      -dd52*bksb2
     &      -dd53*bksb3
     &      -dd54*bksb4
          bksb6=bksb6
     &      -dd61*bksb1
     &      -dd62*bksb2
     &      -dd63*bksb3
     &      -dd64*bksb4
     &      -dd65*bksb5
          bksb6=bksb6/dd66


          bksb5=bksb5
     &      -dd56*bksb6
          bksb5=bksb5/dd55
          bksb4=bksb4
     &      -dd45*bksb5
     &      -dd46*bksb6
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
     &      -dd35*bksb5
     &      -dd36*bksb6
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
     &      -dd25*bksb5
     &      -dd26*bksb6
          bksb2=bksb2/dd22

          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
     &      -dd15*bksb5
     &      -dd16*bksb6
          bksb1=bksb1/dd11








          b(ijind+nb16)=bksb1
          b(ijind+nb26)=bksb2
          b(ijind+nb36)=bksb3
          b(ijind+nb46)=bksb4
          b(ijind+nb56)=bksb5
          b(ijind+nb66)=bksb6


  350   continue
        im1nq2=im1nq2+36
  360 continue
      return
      end
