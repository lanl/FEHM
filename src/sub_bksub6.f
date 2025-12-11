      subroutine sub_bksub6(nel,b,r,nb,nrhs,nop
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
CD1 To perform backsubstitution for a 6 degree of freedom (dof) problem.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 4-6-94       G. Zyvoloski   97      Initial implementation
CD2 6-24-94      B. Robinson    97      Made 6 dof routine from 4 dof 
CD2                                        version
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/sub_bksub6.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 19 12:21:38 1996   llt
CD2 added log history
CD2
CD2    Rev 1.3   09/16/96 11:48:02   robinson
CD2 Prolog change
CD2
CD2    Rev 1.2   09/12/96  08:26:14  robinson
CD2 Prolog Changes
CD2
CD2    Rev 1.1   02/02/96  12:25:58  robinson
CD2 Updated Requirements Traceability
CD2
CD2    Rev 1.0   01/16/96  14:48:26  llt
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
CPS BEGIN sub_bksub6
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
CPS END sub_bksub6
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
      integer nelm1, nr1, nr2, nr3, nr4, nr5, nr6, nb11, nb12
      integer nb13, nb14
      integer nb21, nb22, nb23, nb24, nb31, nb32, nb33, nb34
      integer nb41, nb42, nb43, nb44, im2neq, im2nq2, im1, nelp1, i
      integer nb51, nb52, nb53, nb54, nb55, nb56, nb35, nb36, nb45, nb46
      integer nb61, nb62, nb63, nb64, nb65, nb66, nb15, nb16, nb25, nb26
      integer im1neq, im1nq2, ifin, ijind, ipvtp1
      integer ij, ista, ipvtm1, j
      real*8 bksb1, bksb2, bksb3, bksb4, sum1, sum2, sum3, sum4
      real*8 dd11, dd12, dd13, dd14, dd21, dd22, dd23, dd24 
      real*8 dd31, dd32, dd33, dd34, dd41, dd42, dd43, dd44
      real*8 dd51,dd52,dd53,dd54,dd55,dd56 
      real*8 dd61,dd62,dd63,dd64,dd65,dd66 
      real*8 dd15, dd16, dd25, dd26, dd35, dd36, dd45, dd46
      real*8 r1, r2, r3, r4, r5, r6, bksb5, bksb6, sum5, sum6
c
c     forward substitution
c need definition for nr1 etc see ilu_4
      nelm1=nel-1
      nelp1=nel+1
      nr1=nrhs(1)
      nr2=nrhs(2)
      nr3=nrhs(3)
      nr4=nrhs(4)
      nr5=nrhs(5)
      nr6=nrhs(6)
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
      im2neq=0
      im2nq2=0
      do 460 im1=1,nelm1
        dd11=piv(im2nq2+ 1)
        dd12=piv(im2nq2+ 2)
        dd13=piv(im2nq2+ 3)
        dd14=piv(im2nq2+ 4)
        dd15=piv(im2nq2+ 5)
        dd16=piv(im2nq2+ 6)
        dd21=piv(im2nq2+ 7)
        dd22=piv(im2nq2+ 8)
        dd23=piv(im2nq2+ 9)
        dd24=piv(im2nq2+10)
        dd25=piv(im2nq2+11)
        dd26=piv(im2nq2+12)
        dd31=piv(im2nq2+13)
        dd32=piv(im2nq2+14)
        dd33=piv(im2nq2+15)
        dd34=piv(im2nq2+16)
        dd35=piv(im2nq2+17)
        dd36=piv(im2nq2+18)
        dd41=piv(im2nq2+19)
        dd42=piv(im2nq2+20)
        dd43=piv(im2nq2+21)
        dd44=piv(im2nq2+22)
        dd45=piv(im2nq2+23)
        dd46=piv(im2nq2+24)
        dd51=piv(im2nq2+25)
        dd52=piv(im2nq2+26)
        dd53=piv(im2nq2+27)
        dd54=piv(im2nq2+28)
        dd55=piv(im2nq2+29)
        dd56=piv(im2nq2+30)
        dd61=piv(im2nq2+31)
        dd62=piv(im2nq2+32)
        dd63=piv(im2nq2+33)
        dd64=piv(im2nq2+34)
        dd65=piv(im2nq2+35)
        dd66=piv(im2nq2+36)
        bksb1=r(im1+nr1)
        bksb2=r(im1+nr2)
        bksb3=r(im1+nr3)
        bksb4=r(im1+nr4)
        bksb5=r(im1+nr5)
        bksb6=r(im1+nr6)
        bksb2=bksb2
     &    -dd21*bksb1
        bksb3=bksb3
     &    -dd31*bksb1
     &    -dd32*bksb2
        bksb4=bksb4
     &    -dd41*bksb1
     &    -dd42*bksb2
     &    -dd43*bksb3
        bksb5=bksb5
     &    -dd51*bksb1
     &    -dd52*bksb2
     &    -dd53*bksb3
     &    -dd54*bksb4
        bksb6=bksb6
     &    -dd61*bksb1
     &    -dd62*bksb2
     &    -dd63*bksb3
     &    -dd64*bksb4
     &    -dd65*bksb5
        bksb6=bksb6/dd66
        bksb5=bksb5
     &    -dd56*bksb6
        bksb5=bksb5/dd55
        bksb4=bksb4
     &    -dd45*bksb5
     &    -dd46*bksb6
        bksb4=bksb4/dd44
        bksb3=bksb3
     &    -dd34*bksb4
     &    -dd35*bksb5
     &    -dd36*bksb6
        bksb3=bksb3/dd33
        bksb2=bksb2
     &    -dd23*bksb3
     &    -dd24*bksb4
     &    -dd25*bksb5
     &    -dd26*bksb6
        bksb2=bksb2/dd22
        bksb1=bksb1
     &    -dd12*bksb2
     &    -dd13*bksb3
     &    -dd14*bksb4
     &    -dd15*bksb5
     &    -dd16*bksb6
        bksb1=bksb1/dd11
        r(im1+nr1)=bksb1
        r(im1+nr2)=bksb2
        r(im1+nr3)=bksb3
        r(im1+nr4)=bksb4
        r(im1+nr5)=bksb5
        r(im1+nr6)=bksb6
        im2neq=im2neq+6
        im2nq2=im2nq2+36
        i=im1+1
        ista=nop(i)+1
        ipvtm1=npvt(i)-1
c       factor row i
        sum1=0.0d00
        sum2=0.0d00
        sum3=0.0d00
        sum4=0.0d00
        sum5=0.0d00
        sum6=0.0d00
        do 440 ij=ista,ipvtm1
          j=nop(ij)
          ijind=ij-nelp1
          r1=r(j+nr1)
          r2=r(j+nr2)
          r3=r(j+nr3)
          r4=r(j+nr4)
          r5=r(j+nr5)
          r6=r(j+nr6)
          sum1=sum1+b(ijind+nb11)*r1
     &      +b(ijind+nb12)*r2
     &      +b(ijind+nb13)*r3
     &      +b(ijind+nb14)*r4
     &      +b(ijind+nb15)*r5
     &      +b(ijind+nb16)*r6
          sum2=sum2+b(ijind+nb21)*r1
     &      +b(ijind+nb22)*r2
     &      +b(ijind+nb23)*r3
     &      +b(ijind+nb24)*r4
     &      +b(ijind+nb25)*r5
     &      +b(ijind+nb26)*r6
          sum3=sum3+b(ijind+nb31)*r1
     &      +b(ijind+nb32)*r2
     &      +b(ijind+nb33)*r3
     &      +b(ijind+nb34)*r4
     &      +b(ijind+nb35)*r5
     &      +b(ijind+nb36)*r6
          sum4=sum4+b(ijind+nb41)*r1
     &      +b(ijind+nb42)*r2
     &      +b(ijind+nb43)*r3
     &      +b(ijind+nb44)*r4
     &      +b(ijind+nb45)*r5
     &      +b(ijind+nb46)*r6
          sum5=sum5+b(ijind+nb51)*r1
     &      +b(ijind+nb52)*r2
     &      +b(ijind+nb53)*r3
     &      +b(ijind+nb54)*r4
     &      +b(ijind+nb55)*r5
     &      +b(ijind+nb56)*r6
          sum6=sum6+b(ijind+nb61)*r1
     &      +b(ijind+nb62)*r2
     &      +b(ijind+nb63)*r3
     &      +b(ijind+nb64)*r4
     &      +b(ijind+nb65)*r5
     &      +b(ijind+nb66)*r6
  440   continue
        r(i+nr1)=r(i+nr1)-sum1
        r(i+nr2)=r(i+nr2)-sum2
        r(i+nr3)=r(i+nr3)-sum3
        r(i+nr4)=r(i+nr4)-sum4
        r(i+nr5)=r(i+nr5)-sum5
        r(i+nr6)=r(i+nr6)-sum6
  460 continue
      im1neq=im2neq
      im1nq2=im2nq2
      dd11=piv(im1nq2+ 1)
      dd12=piv(im1nq2+ 2)
      dd13=piv(im1nq2+ 3)
      dd14=piv(im1nq2+ 4)
      dd15=piv(im1nq2+ 5)
      dd16=piv(im1nq2+ 6)
      dd21=piv(im1nq2+ 7)
      dd22=piv(im1nq2+ 8)
      dd23=piv(im1nq2+ 9)
      dd24=piv(im1nq2+10)
      dd25=piv(im1nq2+11)
      dd26=piv(im1nq2+12)
      dd31=piv(im1nq2+13)
      dd32=piv(im1nq2+14)
      dd33=piv(im1nq2+15)
      dd34=piv(im1nq2+16)
      dd35=piv(im1nq2+17)
      dd36=piv(im1nq2+18)
      dd41=piv(im1nq2+19)
      dd42=piv(im1nq2+20)
      dd43=piv(im1nq2+21)
      dd44=piv(im1nq2+22)
      dd45=piv(im1nq2+23)
      dd46=piv(im1nq2+24)
      dd51=piv(im1nq2+25)
      dd52=piv(im1nq2+26)
      dd53=piv(im1nq2+27)
      dd54=piv(im1nq2+28)
      dd55=piv(im1nq2+29)
      dd56=piv(im1nq2+30)
      dd61=piv(im1nq2+31)
      dd62=piv(im1nq2+32)
      dd63=piv(im1nq2+33)
      dd64=piv(im1nq2+34)
      dd65=piv(im1nq2+35)
      dd66=piv(im1nq2+36)
      bksb1=r(nel+nr1)
      bksb2=r(nel+nr2)
      bksb3=r(nel+nr3)
      bksb4=r(nel+nr4)
      bksb5=r(nel+nr5)
      bksb6=r(nel+nr6)
      bksb2=bksb2
     &  -dd21*bksb1
      bksb3=bksb3
     &  -dd31*bksb1
     &  -dd32*bksb2
      bksb4=bksb4
     &  -dd41*bksb1
     &  -dd42*bksb2
     &  -dd43*bksb3
      bksb5=bksb5
     &  -dd51*bksb1
     &  -dd52*bksb2
     &  -dd53*bksb3
     &  -dd54*bksb4
      bksb6=bksb6
     &  -dd61*bksb1
     &  -dd62*bksb2
     &  -dd63*bksb3
     &  -dd64*bksb4
     &  -dd65*bksb5
      bksb6=bksb6/dd66
      bksb5=bksb5
     &  -dd56*bksb6
      bksb5=bksb5/dd55
      bksb4=bksb4
     &  -dd45*bksb5
     &  -dd46*bksb6
      bksb4=bksb4/dd44
      bksb3=bksb3
     &  -dd34*bksb4
     &  -dd35*bksb5
     &  -dd36*bksb6
      bksb3=bksb3/dd33
      bksb2=bksb2
     &  -dd23*bksb3
     &  -dd24*bksb4
     &  -dd25*bksb5
     &  -dd26*bksb6
      bksb2=bksb2/dd22
      bksb1=bksb1
     &  -dd12*bksb2
     &  -dd13*bksb3
     &  -dd14*bksb4
     &  -dd15*bksb5
     &  -dd16*bksb6
      bksb1=bksb1/dd11
      r(nel+nr1)=bksb1
      r(nel+nr2)=bksb2
      r(nel+nr3)=bksb3
      r(nel+nr4)=bksb4
      r(nel+nr5)=bksb5
      r(nel+nr6)=bksb6
c     back substitution
      do 560 i=nelm1,1,-1
        ipvtp1=npvt(i)+1
        ifin=nop(i+1)
        sum1=0.0d00
        sum2=0.0d00
        sum3=0.0d00
        sum4=0.0d00
        sum5=0.0d00
        sum6=0.0d00
        do 540 ij=ipvtp1,ifin
          j=nop(ij)
          ijind=ij-nelp1
          r1=r(j+nr1)
          r2=r(j+nr2)
          r3=r(j+nr3)
          r4=r(j+nr4)
          r5=r(j+nr5)
          r6=r(j+nr6)
          sum1=sum1+b(ijind+nb11)*r1
     &      +b(ijind+nb12)*r2
     &      +b(ijind+nb13)*r3
     &      +b(ijind+nb14)*r4
     &      +b(ijind+nb15)*r5
     &      +b(ijind+nb16)*r6
          sum2=sum2+b(ijind+nb21)*r1
     &      +b(ijind+nb22)*r2
     &      +b(ijind+nb23)*r3
     &      +b(ijind+nb24)*r4
     &      +b(ijind+nb25)*r5
     &      +b(ijind+nb26)*r6
          sum3=sum3+b(ijind+nb31)*r1
     &      +b(ijind+nb32)*r2
     &      +b(ijind+nb33)*r3
     &      +b(ijind+nb34)*r4
     &      +b(ijind+nb35)*r5
     &      +b(ijind+nb36)*r6
          sum4=sum4+b(ijind+nb41)*r1
     &      +b(ijind+nb42)*r2
     &      +b(ijind+nb43)*r3
     &      +b(ijind+nb44)*r4
     &      +b(ijind+nb45)*r5
     &      +b(ijind+nb46)*r6
          sum5=sum5+b(ijind+nb51)*r1
     &      +b(ijind+nb52)*r2
     &      +b(ijind+nb53)*r3
     &      +b(ijind+nb54)*r4
     &      +b(ijind+nb55)*r5
     &      +b(ijind+nb56)*r6
          sum6=sum6+b(ijind+nb61)*r1
     &      +b(ijind+nb62)*r2
     &      +b(ijind+nb63)*r3
     &      +b(ijind+nb64)*r4
     &      +b(ijind+nb65)*r5
     &      +b(ijind+nb66)*r6
  540   continue
        r(i+nr1)=r(i+nr1)-sum1
        r(i+nr2)=r(i+nr2)-sum2
        r(i+nr3)=r(i+nr3)-sum3
        r(i+nr4)=r(i+nr4)-sum4
        r(i+nr5)=r(i+nr5)-sum5
        r(i+nr6)=r(i+nr6)-sum6
  560 continue
      return
      end                      
c    This one is done
