      subroutine sub_bksub4(nel,b,r,nb,nrhs,nop
     *     ,npvt,piv)         
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
CD1 To perform backsubstitution for a 4 degree of freedom (dof) problem.
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_bksub4.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:10 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:24:46 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:18:56   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:42   pvcs
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
CD4 None
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
CD5 nelm1        integer     nel-1
CD5 nelp1        integer     nel+1
CD5 i            integer     Do loop index
CD5 nr1          integer     Pointer for first unknown in solution
CD5                              vector
CD5 nr2          integer     Pointer for second unknown in solution
CD5                              vector
CD5 nr3          integer     Pointer for third unknown in solution
CD5                              vector
CD5 nr4          integer     Pointer for fourth unknown in solution
CD5                              vector
CD5 nb11         integer     Position in b array
CD5 nb12         integer     Position in b array
CD5 nb13         integer     Position in b array
CD5 nb14         integer     Position in b array
CD5 nb21         integer     Position in b array
CD5 nb22         integer     Position in b array
CD5 nb23         integer     Position in b array
CD5 nb24         integer     Position in b array
CD5 nb31         integer     Position in b array
CD5 nb32         integer     Position in b array
CD5 nb33         integer     Position in b array
CD5 nb34         integer     Position in b array
CD5 nb41         integer     Position in b array
CD5 nb42         integer     Position in b array
CD5 nb43         integer     Position in b array
CD5 nb44         integer     Position in b array
CD5 im1neq       integer     Integer index
CD5 im1nq2       integer     Integer index
CD5 im2neq       integer     Integer index
CD5 im2nq2       integer     Integer index
CD5 im1          integer     Do loop index
CD5 i            integer     Integer index
CD5 ifin         integer     Do loop limit parameter
CD5 ijind        integer     Integer index
CD5 ipvtm1       integer     Do loop limit parameter
CD5 ipvtp1       integer     Do loop limit parameter
CD5 ij           integer     Do loop index
CD5 ista         integer     Do loop limit parameter
CD5 j            integer     Integer index
CD5 bksb1        real*8      Term used in calculation
CD5 bksb2        real*8      Term used in calculation
CD5 bksb3        real*8      Term used in calculation
CD5 bksb4        real*8      Term used in calculation
CD5 dd11         real*8      Pivot value
CD5 dd12         real*8      Pivot value
CD5 dd13         real*8      Pivot value
CD5 dd14         real*8      Pivot value
CD5 dd21         real*8      Pivot value
CD5 dd22         real*8      Pivot value
CD5 dd23         real*8      Pivot value
CD5 dd24         real*8      Pivot value
CD5 dd31         real*8      Pivot value
CD5 dd32         real*8      Pivot value
CD5 dd33         real*8      Pivot value
CD5 dd34         real*8      Pivot value
CD5 dd41         real*8      Pivot value
CD5 dd42         real*8      Pivot value
CD5 dd43         real*8      Pivot value
CD5 dd44         real*8      Pivot value
CD5 sum1         real*8      Term in calculation of factorization
CD5 sum2         real*8      Term in calculation of factorization
CD5 sum3         real*8      Term in calculation of factorization
CD5 sum4         real*8      Term in calculation of factorization
CD5 r1           real*8      Current residual value, first unknown
CD5 r2           real*8      Current residual value, second unknown
CD5 r3           real*8      Current residual value, third unknown
CD5 r4           real*8      Current residual value, fourth unknown
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
CD8 3.2 Solve Linear Equation Set
CD8    3.2.3 Form Approximate Solution
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
CPS BEGIN sub_bksub4
CPS 
CPS Define integer parameters
CPS 
CPS FOR each equation
CPS 
CPS   Identify diagonal terms, set values of pivot parameters
CPS   Multiply right hand side by pivot matrix
CPS   
CPS   FOR each nonzero element of lu factorization matrix up to one...
CPS   ... less than the pivot
CPS     Add to running sum of intermediate value of solution
CPS   ENDFOR
CPS 
CPS ENDFOR each equation
CPS 
CPS Calculate solution for last equation
CPS 
CPS FOR each equation from last -1 to first
CPS   
CPS   FOR each nonzero element of lu factorization matrix from pivot...
CPS   ... plus 1 to last nonzero value
CPS     Add to running sum of intermediate value of solution
CPS   ENDFOR
CPS 
CPS ENDFOR
CPS 
CPS END sub_bksub4
CPS
C**********************************************************************
c
c
c     forward and back substitution for 4 dof problems

      implicit none

      real*8 b(*),r(*),piv(*)
      integer nop(*),npvt(*)
      integer nb(*),nrhs(*)
      integer nel
      integer nelm1, nr1, nr2, nr3, nr4, nb11, nb12, nb13, nb14
      integer nb21, nb22, nb23, nb24, nb31, nb32, nb33, nb34
      integer nb41, nb42, nb43, nb44, im2neq, im2nq2, im1, nelp1, i
      integer im1neq, im1nq2, ifin, ijind, ipvtp1
      integer ij, ista, ipvtm1, j
      real*8 bksb1, bksb2, bksb3, bksb4, sum1, sum2, sum3, sum4
      real*8 dd11, dd12, dd13, dd14, dd21, dd22, dd23, dd24 
      real*8 dd31, dd32, dd33, dd34, dd41, dd42, dd43, dd44 
      real*8 r1, r2, r3, r4
c
c     forward substitution
c need definition for nr1 etc see ilu_4
      nelm1=nel-1
      nelp1=nel+1
      nr1=nrhs(1)
      nr2=nrhs(2)
      nr3=nrhs(3)
      nr4=nrhs(4)
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
      im2neq=0
      im2nq2=0
      do 460 im1=1,nelm1
        dd11=piv(im2nq2+ 1)
        dd12=piv(im2nq2+ 2)
        dd13=piv(im2nq2+ 3)
        dd14=piv(im2nq2+ 4)
        dd21=piv(im2nq2+ 5)
        dd22=piv(im2nq2+ 6)
        dd23=piv(im2nq2+ 7)
        dd24=piv(im2nq2+ 8)
        dd31=piv(im2nq2+ 9)
        dd32=piv(im2nq2+10)
        dd33=piv(im2nq2+11)
        dd34=piv(im2nq2+12)
        dd41=piv(im2nq2+13)
        dd42=piv(im2nq2+14)
        dd43=piv(im2nq2+15)
        dd44=piv(im2nq2+16)
        bksb1=r(im1+nr1)
        bksb2=r(im1+nr2)
        bksb3=r(im1+nr3)
        bksb4=r(im1+nr4)
        bksb2=bksb2
     &    -dd21*bksb1
        bksb3=bksb3
     &    -dd31*bksb1
     &    -dd32*bksb2
        bksb4=bksb4
     &    -dd41*bksb1
     &    -dd42*bksb2
     &    -dd43*bksb3
        bksb4=bksb4/dd44
        bksb3=bksb3
     &    -dd34*bksb4
        bksb3=bksb3/dd33
        bksb2=bksb2
     &    -dd23*bksb3
     &    -dd24*bksb4
        bksb2=bksb2/dd22
        bksb1=bksb1
     &    -dd12*bksb2
     &    -dd13*bksb3
     &    -dd14*bksb4
        bksb1=bksb1/dd11
        r(im1+nr1)=bksb1
        r(im1+nr2)=bksb2
        r(im1+nr3)=bksb3
        r(im1+nr4)=bksb4
        im2neq=im2neq+4
        im2nq2=im2nq2+16
        i=im1+1
        ista=nop(i)+1
        ipvtm1=npvt(i)-1
c       factor row i
        sum1=0.0d00
        sum2=0.0d00
        sum3=0.0d00
        sum4=0.0d00
        do 440 ij=ista,ipvtm1
          j=nop(ij)
          ijind=ij-nelp1
          r1=r(j+nr1)
          r2=r(j+nr2)
          r3=r(j+nr3)
          r4=r(j+nr4)
          sum1=sum1+b(ijind+nb11)*r1
     &      +b(ijind+nb12)*r2
     &      +b(ijind+nb13)*r3
     &      +b(ijind+nb14)*r4
          sum2=sum2+b(ijind+nb21)*r1
     &      +b(ijind+nb22)*r2
     &      +b(ijind+nb23)*r3
     &      +b(ijind+nb24)*r4
          sum3=sum3+b(ijind+nb31)*r1
     &      +b(ijind+nb32)*r2
     &      +b(ijind+nb33)*r3
     &      +b(ijind+nb34)*r4
          sum4=sum4+b(ijind+nb41)*r1
     &      +b(ijind+nb42)*r2
     &      +b(ijind+nb43)*r3
     &      +b(ijind+nb44)*r4
  440   continue
        r(i+nr1)=r(i+nr1)-sum1
        r(i+nr2)=r(i+nr2)-sum2
        r(i+nr3)=r(i+nr3)-sum3
        r(i+nr4)=r(i+nr4)-sum4
  460 continue
      im1neq=im2neq
      im1nq2=im2nq2
      dd11=piv(im1nq2+ 1)
      dd12=piv(im1nq2+ 2)
      dd13=piv(im1nq2+ 3)
      dd14=piv(im1nq2+ 4)
      dd21=piv(im1nq2+ 5)
      dd22=piv(im1nq2+ 6)
      dd23=piv(im1nq2+ 7)
      dd24=piv(im1nq2+ 8)
      dd31=piv(im1nq2+ 9)
      dd32=piv(im1nq2+10)
      dd33=piv(im1nq2+11)
      dd34=piv(im1nq2+12)
      dd41=piv(im1nq2+13)
      dd42=piv(im1nq2+14)
      dd43=piv(im1nq2+15)
      dd44=piv(im1nq2+16)
      bksb1=r(nel+nr1)
      bksb2=r(nel+nr2)
      bksb3=r(nel+nr3)
      bksb4=r(nel+nr4)
      bksb2=bksb2
     &  -dd21*bksb1
      bksb3=bksb3
     &  -dd31*bksb1
     &  -dd32*bksb2
      bksb4=bksb4
     &  -dd41*bksb1
     &  -dd42*bksb2
     &  -dd43*bksb3
      bksb4=bksb4/dd44
      bksb3=bksb3
     &  -dd34*bksb4
      bksb3=bksb3/dd33
      bksb2=bksb2
     &  -dd23*bksb3
     &  -dd24*bksb4
      bksb2=bksb2/dd22
      bksb1=bksb1
     &  -dd12*bksb2
     &  -dd13*bksb3
     &  -dd14*bksb4
      bksb1=bksb1/dd11
      r(nel+nr1)=bksb1
      r(nel+nr2)=bksb2
      r(nel+nr3)=bksb3
      r(nel+nr4)=bksb4
c     back substitution
      do 560 i=nelm1,1,-1
        ipvtp1=npvt(i)+1
        ifin=nop(i+1)
        sum1=0.0d00
        sum2=0.0d00
        sum3=0.0d00
        sum4=0.0d00
        do 540 ij=ipvtp1,ifin
          j=nop(ij)
          ijind=ij-nelp1
          r1=r(j+nr1)
          r2=r(j+nr2)
          r3=r(j+nr3)
          r4=r(j+nr4)
          sum1=sum1+b(ijind+nb11)*r1
     &      +b(ijind+nb12)*r2
     &      +b(ijind+nb13)*r3
     &      +b(ijind+nb14)*r4
          sum2=sum2+b(ijind+nb21)*r1
     &      +b(ijind+nb22)*r2
     &      +b(ijind+nb23)*r3
     &      +b(ijind+nb24)*r4
          sum3=sum3+b(ijind+nb31)*r1
     &      +b(ijind+nb32)*r2
     &      +b(ijind+nb33)*r3
     &      +b(ijind+nb34)*r4
          sum4=sum4+b(ijind+nb41)*r1
     &      +b(ijind+nb42)*r2
     &      +b(ijind+nb43)*r3
     &      +b(ijind+nb44)*r4
  540   continue
        r(i+nr1)=r(i+nr1)-sum1
        r(i+nr2)=r(i+nr2)-sum2
        r(i+nr3)=r(i+nr3)-sum3
        r(i+nr4)=r(i+nr4)-sum4
  560 continue
      return
      end                      
