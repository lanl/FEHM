      subroutine  sub_ilu4(nel,a,b,na,nb,ncon,nop
     *     ,irb,iirb,npvt,dum1,dum2,dum3,dum4,piv)
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
CD1 To perform incomplete lu factorization for 4 degree of freedom
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilu4.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Thu Sep 12 08:26:46 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.5   Fri Feb 02 13:02:34 1996   hend
CD2 Fixed Error in Comments
CD2 
CD2    Rev 1.4   Fri Feb 02 12:43:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   08/18/95 10:33:14   llt
CD2 ik already defined, removed for cray
CD2 
CD2    Rev 1.2   05/11/94 16:19:24   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:46   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:47:00   pvcs
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
CPS BEGIN sub_ilu4
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
CPS END sub_ilu4
CPS
C**********************************************************************
c
c four degree of freedom solver
c
c     4 degree of freedom ilu preconditioner
     
      implicit none

      integer nel
      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),piv(*)
      real*8 dum1(nel,4),dum2(nel,4),dum3(nel,4),dum4(nel,4)
      integer na(*),nb(*)                    
      integer nelp1
      integer i1, i2, i3, i4
      integer na11, na12, na13, na14
      integer na21, na22, na23, na24
      integer na31, na32, na33, na34
      integer na41, na42, na43, na44
      integer nb11, nb12, nb13, nb14
      integer nb21, nb22, nb23, nb24
      integer nb31, nb32, nb33, nb34
      integer nb41, nb42, nb43, nb44
      integer i, ifin, ij, ijind, ikb, ka, kb
      integer im1nq2, ipvt, ista, ipvtp1, ipvtm1, j, ik, k, kj, l
      integer j2, kpvt, ikind, kjind
      integer ipvind, ip
      real*8 a1, a2, b1, b2, a3, a4
      real*8 dd11, dd12, dd13, dd14
      real*8 dd21, dd22, dd23, dd24
      real*8 dd31, dd32, dd33, dd34
      real*8 dd41, dd42, dd43, dd44
      real*8 bksb1, bksb2, bksb3, bksb4
c gaz 1-18-2002
c function definitions not used in this routine
c     real*8 alm, ali
c
c***     linear algebra     ***
c
c     alm(a1,a2,b1,b2)=a1*b1+a2*b2
c     ali(a1,a2,a3,a4,b1)=b1/(a1*a4-a2*a3)
c       
      nelp1=nel+1
      na11=na( 1)
      na12=na( 2)
      na13=na( 3)
      na14=na( 4)
      na21=na( 5)
      na22=na( 6)
      na23=na( 7)
      na24=na( 8)
      na31=na( 9)
      na32=na(10)
      na33=na(11)
      na34=na(12)
      na41=na(13)
      na42=na(14)
      na43=na(15)
      na44=na(16)
      nb11=nb( 1)
      nb12=nb( 2)
      nb13=nb( 3)
      nb14=nb( 4)
      nb21=nb( 5)
      nb22=nb( 6)
      nb23=nb( 7)
      nb24=nb( 8)
      nb31=nb( 9)
      nb32=nb(10)
      nb33=nb(11)
      nb34=nb(12)
      nb41=nb(13)
      nb42=nb(14)
      nb43=nb(15)
      nb44=nb(16)
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
            dum2(ikb,1)=0.0d00
            dum2(ikb,2)=0.0d00
            dum2(ikb,3)=0.0d00
            dum2(ikb,4)=0.0d00
            dum3(ikb,1)=0.0d00
            dum3(ikb,2)=0.0d00
            dum3(ikb,3)=0.0d00
            dum3(ikb,4)=0.0d00
            dum4(ikb,1)=0.0d00
            dum4(ikb,2)=0.0d00
            dum4(ikb,3)=0.0d00
            dum4(ikb,4)=0.0d00
          enddo
          do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            dum1(kb,1)=a(k-nelp1+na11)
            dum1(kb,2)=a(k-nelp1+na12)
            dum1(kb,3)=a(k-nelp1+na13)
            dum1(kb,4)=a(k-nelp1+na14)
            dum2(kb,1)=a(k-nelp1+na21)
            dum2(kb,2)=a(k-nelp1+na22)
            dum2(kb,3)=a(k-nelp1+na23)
            dum2(kb,4)=a(k-nelp1+na24)
            dum3(kb,1)=a(k-nelp1+na31)
            dum3(kb,2)=a(k-nelp1+na32)
            dum3(kb,3)=a(k-nelp1+na33)
            dum3(kb,4)=a(k-nelp1+na34)
            dum4(kb,1)=a(k-nelp1+na41)
            dum4(kb,2)=a(k-nelp1+na42)
            dum4(kb,3)=a(k-nelp1+na43)
            dum4(kb,4)=a(k-nelp1+na44)
          enddo
          do l=i1,i2
            kb=nop(l)
            b(l-nelp1+nb11)=dum1(kb,1)
            b(l-nelp1+nb12)=dum1(kb,2)
            b(l-nelp1+nb13)=dum1(kb,3)
            b(l-nelp1+nb14)=dum1(kb,4)
            b(l-nelp1+nb21)=dum2(kb,1)
            b(l-nelp1+nb22)=dum2(kb,2)
            b(l-nelp1+nb23)=dum2(kb,3)
            b(l-nelp1+nb24)=dum2(kb,4)
            b(l-nelp1+nb31)=dum3(kb,1)
            b(l-nelp1+nb32)=dum3(kb,2)
            b(l-nelp1+nb33)=dum3(kb,3)
            b(l-nelp1+nb34)=dum3(kb,4)
            b(l-nelp1+nb41)=dum4(kb,1)
            b(l-nelp1+nb42)=dum4(kb,2)
            b(l-nelp1+nb43)=dum4(kb,3)
            b(l-nelp1+nb44)=dum4(kb,4)
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
              b(ijind+nb12)=b(ijind+nb12)
     &          -b(ikind+nb11)*b(kjind+nb12)
     &          -b(ikind+nb12)*b(kjind+nb22)
     &          -b(ikind+nb13)*b(kjind+nb32)
     &          -b(ikind+nb14)*b(kjind+nb42)
              b(ijind+nb13)=b(ijind+nb13)
     &          -b(ikind+nb11)*b(kjind+nb13)
     &          -b(ikind+nb12)*b(kjind+nb23)
     &          -b(ikind+nb13)*b(kjind+nb33)
     &          -b(ikind+nb14)*b(kjind+nb43)
              b(ijind+nb14)=b(ijind+nb14)
     &          -b(ikind+nb11)*b(kjind+nb14)
     &          -b(ikind+nb12)*b(kjind+nb24)
     &          -b(ikind+nb13)*b(kjind+nb34)
     &          -b(ikind+nb14)*b(kjind+nb44)
              b(ijind+nb21)=b(ijind+nb21)
     &          -b(ikind+nb21)*b(kjind+nb11)
     &          -b(ikind+nb22)*b(kjind+nb21)
     &          -b(ikind+nb23)*b(kjind+nb31)
     &          -b(ikind+nb24)*b(kjind+nb41)
              b(ijind+nb22)=b(ijind+nb22)
     &          -b(ikind+nb21)*b(kjind+nb12)
     &          -b(ikind+nb22)*b(kjind+nb22)
     &          -b(ikind+nb23)*b(kjind+nb32)
     &          -b(ikind+nb24)*b(kjind+nb42)
              b(ijind+nb23)=b(ijind+nb23)
     &          -b(ikind+nb21)*b(kjind+nb13)
     &          -b(ikind+nb22)*b(kjind+nb23)
     &          -b(ikind+nb23)*b(kjind+nb33)
     &          -b(ikind+nb24)*b(kjind+nb43)
              b(ijind+nb24)=b(ijind+nb24)
     &          -b(ikind+nb21)*b(kjind+nb14)
     &          -b(ikind+nb22)*b(kjind+nb24)
     &          -b(ikind+nb23)*b(kjind+nb34)
     &          -b(ikind+nb24)*b(kjind+nb44)
              b(ijind+nb31)=b(ijind+nb31)
     &          -b(ikind+nb31)*b(kjind+nb11)
     &          -b(ikind+nb32)*b(kjind+nb21)
     &          -b(ikind+nb33)*b(kjind+nb31)
     &          -b(ikind+nb34)*b(kjind+nb41)
              b(ijind+nb32)=b(ijind+nb32)
     &          -b(ikind+nb31)*b(kjind+nb12)
     &          -b(ikind+nb32)*b(kjind+nb22)
     &          -b(ikind+nb33)*b(kjind+nb32)
     &          -b(ikind+nb34)*b(kjind+nb42)
              b(ijind+nb33)=b(ijind+nb33)
     &          -b(ikind+nb31)*b(kjind+nb13)
     &          -b(ikind+nb32)*b(kjind+nb23)
     &          -b(ikind+nb33)*b(kjind+nb33)
     &          -b(ikind+nb34)*b(kjind+nb43)
              b(ijind+nb34)=b(ijind+nb34)
     &          -b(ikind+nb31)*b(kjind+nb14)
     &          -b(ikind+nb32)*b(kjind+nb24)
     &          -b(ikind+nb33)*b(kjind+nb34)
     &          -b(ikind+nb34)*b(kjind+nb44)
              b(ijind+nb41)=b(ijind+nb41)
     &          -b(ikind+nb41)*b(kjind+nb11)
     &          -b(ikind+nb42)*b(kjind+nb21)
     &          -b(ikind+nb43)*b(kjind+nb31)
     &          -b(ikind+nb44)*b(kjind+nb41)
              b(ijind+nb42)=b(ijind+nb42)
     &          -b(ikind+nb41)*b(kjind+nb12)
     &          -b(ikind+nb42)*b(kjind+nb22)
     &          -b(ikind+nb43)*b(kjind+nb32)
     &          -b(ikind+nb44)*b(kjind+nb42)
              b(ijind+nb43)=b(ijind+nb43)
     &          -b(ikind+nb41)*b(kjind+nb13)
     &          -b(ikind+nb42)*b(kjind+nb23)
     &          -b(ikind+nb43)*b(kjind+nb33)
     &          -b(ikind+nb44)*b(kjind+nb43)
              b(ijind+nb44)=b(ijind+nb44)
     &          -b(ikind+nb41)*b(kjind+nb14)
     &          -b(ikind+nb42)*b(kjind+nb24)
     &          -b(ikind+nb43)*b(kjind+nb34)
     &          -b(ikind+nb44)*b(kjind+nb44)
            endif
  260     continue
  270   continue
c       calculate the inverse of l(i,i)
        ipvind=ipvt-nelp1
        dd11=b(ipvind+nb11)
        dd12=b(ipvind+nb12)
        dd13=b(ipvind+nb13)
        dd14=b(ipvind+nb14)
        dd21=b(ipvind+nb21)
        dd22=b(ipvind+nb22)
        dd23=b(ipvind+nb23)
        dd24=b(ipvind+nb24)
        dd31=b(ipvind+nb31)
        dd32=b(ipvind+nb32)
        dd33=b(ipvind+nb33)
        dd34=b(ipvind+nb34)
        dd41=b(ipvind+nb41)
        dd42=b(ipvind+nb42)
        dd43=b(ipvind+nb43)
        dd44=b(ipvind+nb44)
        dd21=dd21/dd11
        dd31=dd31/dd11
        dd41=dd41/dd11
        dd22=dd22
     &    -dd21*dd12
        dd32=dd32
     &    -dd31*dd12
        dd42=dd42
     &    -dd41*dd12
        dd32=dd32/dd22
        dd42=dd42/dd22
        dd23=dd23
     &    -dd21*dd13
        dd33=dd33
     &    -dd31*dd13
     &    -dd32*dd23
        dd43=dd43
     &    -dd41*dd13
     &    -dd42*dd23
        dd43=dd43/dd33
        dd24=dd24
     &    -dd21*dd14
        dd34=dd34
     &    -dd31*dd14
     &    -dd32*dd24
        dd44=dd44
     &    -dd41*dd14
     &    -dd42*dd24
     &    -dd43*dd34
        piv(im1nq2+ 1)=dd11
        piv(im1nq2+ 2)=dd12
        piv(im1nq2+ 3)=dd13
        piv(im1nq2+ 4)=dd14
        piv(im1nq2+ 5)=dd21
        piv(im1nq2+ 6)=dd22
        piv(im1nq2+ 7)=dd23
        piv(im1nq2+ 8)=dd24
        piv(im1nq2+ 9)=dd31
        piv(im1nq2+10)=dd32
        piv(im1nq2+11)=dd33
        piv(im1nq2+12)=dd34
        piv(im1nq2+13)=dd41
        piv(im1nq2+14)=dd42
        piv(im1nq2+15)=dd43
        piv(im1nq2+16)=dd44
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
              b(ijind+nb12)=b(ijind+nb12)
     &          -b(ikind+nb11)*b(kjind+nb12)
     &          -b(ikind+nb12)*b(kjind+nb22)
     &          -b(ikind+nb13)*b(kjind+nb32)
     &          -b(ikind+nb14)*b(kjind+nb42)
              b(ijind+nb13)=b(ijind+nb13)
     &          -b(ikind+nb11)*b(kjind+nb13)
     &          -b(ikind+nb12)*b(kjind+nb23)
     &          -b(ikind+nb13)*b(kjind+nb33)
     &          -b(ikind+nb14)*b(kjind+nb43)
              b(ijind+nb14)=b(ijind+nb14)
     &          -b(ikind+nb11)*b(kjind+nb14)
     &          -b(ikind+nb12)*b(kjind+nb24)
     &          -b(ikind+nb13)*b(kjind+nb34)
     &          -b(ikind+nb14)*b(kjind+nb44)
              b(ijind+nb21)=b(ijind+nb21)
     &          -b(ikind+nb21)*b(kjind+nb11)
     &          -b(ikind+nb22)*b(kjind+nb21)
     &          -b(ikind+nb23)*b(kjind+nb31)
     &          -b(ikind+nb24)*b(kjind+nb41)
              b(ijind+nb22)=b(ijind+nb22)
     &          -b(ikind+nb21)*b(kjind+nb12)
     &          -b(ikind+nb22)*b(kjind+nb22)
     &          -b(ikind+nb23)*b(kjind+nb32)
     &          -b(ikind+nb24)*b(kjind+nb42)
              b(ijind+nb23)=b(ijind+nb23)
     &          -b(ikind+nb21)*b(kjind+nb13)
     &          -b(ikind+nb22)*b(kjind+nb23)
     &          -b(ikind+nb23)*b(kjind+nb33)
     &          -b(ikind+nb24)*b(kjind+nb43)
              b(ijind+nb24)=b(ijind+nb24)
     &          -b(ikind+nb21)*b(kjind+nb14)
     &          -b(ikind+nb22)*b(kjind+nb24)
     &          -b(ikind+nb23)*b(kjind+nb34)
     &          -b(ikind+nb24)*b(kjind+nb44)
              b(ijind+nb31)=b(ijind+nb31)
     &          -b(ikind+nb31)*b(kjind+nb11)
     &          -b(ikind+nb32)*b(kjind+nb21)
     &          -b(ikind+nb33)*b(kjind+nb31)
     &          -b(ikind+nb34)*b(kjind+nb41)
              b(ijind+nb32)=b(ijind+nb32)
     &          -b(ikind+nb31)*b(kjind+nb12)
     &          -b(ikind+nb32)*b(kjind+nb22)
     &          -b(ikind+nb33)*b(kjind+nb32)
     &          -b(ikind+nb34)*b(kjind+nb42)
              b(ijind+nb33)=b(ijind+nb33)
     &          -b(ikind+nb31)*b(kjind+nb13)
     &          -b(ikind+nb32)*b(kjind+nb23)
     &          -b(ikind+nb33)*b(kjind+nb33)
     &          -b(ikind+nb34)*b(kjind+nb43)
              b(ijind+nb34)=b(ijind+nb34)
     &          -b(ikind+nb31)*b(kjind+nb14)
     &          -b(ikind+nb32)*b(kjind+nb24)
     &          -b(ikind+nb33)*b(kjind+nb34)
     &          -b(ikind+nb34)*b(kjind+nb44)
              b(ijind+nb41)=b(ijind+nb41)
     &          -b(ikind+nb41)*b(kjind+nb11)
     &          -b(ikind+nb42)*b(kjind+nb21)
     &          -b(ikind+nb43)*b(kjind+nb31)
     &          -b(ikind+nb44)*b(kjind+nb41)
              b(ijind+nb42)=b(ijind+nb42)
     &          -b(ikind+nb41)*b(kjind+nb12)
     &          -b(ikind+nb42)*b(kjind+nb22)
     &          -b(ikind+nb43)*b(kjind+nb32)
     &          -b(ikind+nb44)*b(kjind+nb42)
              b(ijind+nb43)=b(ijind+nb43)
     &          -b(ikind+nb41)*b(kjind+nb13)
     &          -b(ikind+nb42)*b(kjind+nb23)
     &          -b(ikind+nb43)*b(kjind+nb33)
     &          -b(ikind+nb44)*b(kjind+nb43)
              b(ijind+nb44)=b(ijind+nb44)
     &          -b(ikind+nb41)*b(kjind+nb14)
     &          -b(ikind+nb42)*b(kjind+nb24)
     &          -b(ikind+nb43)*b(kjind+nb34)
     &          -b(ikind+nb44)*b(kjind+nb44)
            endif
  340     continue
c         u(i,j)=u(i,j)/l(i,i)
          bksb1=b(ijind+nb11)
          bksb2=b(ijind+nb21)
          bksb3=b(ijind+nb31)
          bksb4=b(ijind+nb41)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
          bksb2=bksb2/dd22
          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
          bksb1=bksb1/dd11
          b(ijind+nb11)=bksb1
          b(ijind+nb21)=bksb2
          b(ijind+nb31)=bksb3
          b(ijind+nb41)=bksb4
          bksb1=b(ijind+nb12)
          bksb2=b(ijind+nb22)
          bksb3=b(ijind+nb32)
          bksb4=b(ijind+nb42)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
          bksb2=bksb2/dd22
          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
          bksb1=bksb1/dd11
          b(ijind+nb12)=bksb1
          b(ijind+nb22)=bksb2
          b(ijind+nb32)=bksb3
          b(ijind+nb42)=bksb4
          bksb1=b(ijind+nb13)
          bksb2=b(ijind+nb23)
          bksb3=b(ijind+nb33)
          bksb4=b(ijind+nb43)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
          bksb2=bksb2/dd22
          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
          bksb1=bksb1/dd11
          b(ijind+nb13)=bksb1
          b(ijind+nb23)=bksb2
          b(ijind+nb33)=bksb3
          b(ijind+nb43)=bksb4
          bksb1=b(ijind+nb14)
          bksb2=b(ijind+nb24)
          bksb3=b(ijind+nb34)
          bksb4=b(ijind+nb44)
          bksb2=bksb2
     &      -dd21*bksb1
          bksb3=bksb3
     &      -dd31*bksb1
     &      -dd32*bksb2
          bksb4=bksb4
     &      -dd41*bksb1
     &      -dd42*bksb2
     &      -dd43*bksb3
          bksb4=bksb4/dd44
          bksb3=bksb3
     &      -dd34*bksb4
          bksb3=bksb3/dd33
          bksb2=bksb2
     &      -dd23*bksb3
     &      -dd24*bksb4
          bksb2=bksb2/dd22
          bksb1=bksb1
     &      -dd12*bksb2
     &      -dd13*bksb3
     &      -dd14*bksb4
          bksb1=bksb1/dd11
          b(ijind+nb14)=bksb1
          b(ijind+nb24)=bksb2
          b(ijind+nb34)=bksb3
          b(ijind+nb44)=bksb4
  350   continue
        im1nq2=im1nq2+16
  360 continue
      return
      end
