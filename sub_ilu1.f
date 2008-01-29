      subroutine  sub_ilu1(neq,a,b,ncon,nop
     *     ,irb,iirb,npvt,sorthm,dum,piv) 
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
CD1 To perform incomplete lu factorization for 1 degree of freedom
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilu1.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:26:38 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri Feb 02 12:40:28 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:19:16   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:40   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:54   pvcs
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
CD3 ncon         integer  I      Connectivity matrix for solution matrix
CD3 nop          integer  I      Connectivity matrix for factorization
CD3                                 matrix
CD3 irb          integer  I      Renumber array - inew=irb(iold)
CD3 iirb         integer  I      Renumber array - iold=iirb(inew)
CD3 sorthm       real*8   I      Scratch storage array
CD3 dum          real*8   I      Scratch storage array
CD3 npvt         integer  I      Pivot positions in nop
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
CD5 i            integer     Do loop index
CD5 i1           integer     Do loop index
CD5 i2           integer     Do loop index
CD5 i3           integer     Do loop index
CD5 i4           integer     Do loop index
CD5 i5           integer     Do loop index
CD5 i6           integer     Do loop index
CD5 ip           integer     Pointer to reorder array
CD5 ik           integer     Do loop index
CD5 k            integer     Do loop index
CD5 ka           integer     Pointer in connectivity array
CD5 kb           integer     Pointer to reorder array
CD5 l            integer     Do loop index
CD5 npiv         integer     Pivot position
CD5 k1           integer     k+1
CD5 npivk1       integer     Index in pivot array
CD5 nz           integer     Pointer to multiplier of previous row
CD5 jj           integer     Do loop index
CD5 nj           integer     Position in lu pointer array
CD5 npivnj       integer     Index in pivot array
CD5 ii           integer     Do loop index
CD5 nk           integer     Position in factorization connectivity
CD5                             array
CD5 ijk          integer     Do loop index
CD5 pivot        real*8      Pivot value
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
CD8    3.2.1 Perform Incomplete Factorization
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
CPS BEGIN sub_ilu1
CPS 
CPS Define integer parameters
CPS 
CPS FOR each equation
CPS 
CPS   Determine index parameters
CPS   
CPS   For each position in dummy storage array
CPS     Initialize to zero
CPS   ENDFOR
CPS   
CPS   FOR each nonzero position in solution array
CPS     Copy solution matrix value to dummy storage array
CPS   ENDFOR
CPS   
CPS   FOR each nonzero position in lu factorization  array
CPS     Copy dummy storage array value into lu factorization array
CPS   ENDFOR
CPS   
CPS ENDFOR each equation
CPS 
CPS FOR each equation except the last
CPS 
CPS   Identify pivot position
CPS   Invert diagonal element to form pivot
CPS   
CPS   FOR each nonzero element in current equation
CPS     Multiply by pivot
CPS   ENDFOR
CPS   
CPS   Save pivot value in array
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
CPS END sub_ilu1
CPS
C**********************************************************************

c
c one degree of freedom solver
c
c     call ilu factorization at preconditioner stage
 
      implicit none

      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),dum(*),sorthm(*),piv(*)
      integer neq
      integer neqm1, neqp1, i, i1, i2, ip, i3, i4, ik, k, ka, kb, l
      integer npiv, k1, npivk1, nz, i5, i6, jj, nj, npivnj, ii, nk, ijk
      real*8 pivot, c1
c
c     define some parameters
c
      neqm1=neq-1
      neqp1=neq+1
c     
c     factor matrix 
c
c     first load b matrix from a matrix
c
      do  i=1,neq
         i1=nop(i)+1
         i2=nop(i+1)
         ip=irb(i)
         i3=ncon(ip)+1
         i4=ncon(ip+1)
c     zero out space in row for lu factors
         do  ik=i1,i2
            dum(nop(ik))=0.0d00
         enddo
c     load solution matrix values into row srorage 
         do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            dum(kb)=a(k-neqp1)
         enddo
c     load solution matrix into space for lu factorization(b)
         do l=i1,i2
            b(l-neqp1)=dum(nop(l))
         enddo
      enddo
c     
c     start factorization
c     
      do k=1,neqm1
c     identify pivot
         npiv=npvt(k)
         pivot=1./b(npiv-neqp1)
c     multiply row k by pivot
         i2=nop(k+1)
         do  i=npiv,i2
            b(i-neqp1)=b(i-neqp1)*pivot
         enddo
c     save pivot
         piv(k)=pivot
         k1=k+1
c     factor row k1
         npivk1=npvt(k1)-1
         nz=0
         i5=nop(k1)+1
         i6=nop(k1+1)
c     store copy of row k1
         do  i=i5,i6
            dum(nop(i))=b(i-neqp1)
         enddo
         do  jj=i5,npivk1
            nj=nop(jj)
            c1=dum(nj)
            nz=nz+1
            sorthm(nz)=c1
c     add linear combo of row nj
            i4=nop(nj+1)
            npivnj=npvt(nj)
            do ii=npivnj,i4
               nk=nop(ii)
               dum(nk)=dum(nk)-c1*b(ii-neqp1)
            enddo
         enddo
         do  i=npivk1+1,i6
            b(i-neqp1)=dum(nop(i))
         enddo
         nz=0
         do ijk=i5,npivk1
            nz=nz+1
            b(ijk-neqp1)=sorthm(nz)
         enddo
      enddo
      piv(neq)=1./b(nop(neqp1)-neqp1)
c
      return
      end
