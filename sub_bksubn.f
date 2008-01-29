      subroutine sub_bksubn(nel,idof,b,r,nb,nrhs,nop
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
CD1 To perform backsubstitution for a n degree of freedom (dof) problem.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sub_bksubn.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:08   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Thu Sep 12 08:26:16 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.1   Fri Feb 02 14:42:16 1996   hend
CD2 Updated Prolog and Log
CD2
CD2 4-6-94       G. Zyvoloski   97      Initial implementation
CD2 6-24-94      B. Robinson    97      Made 6 dof routine from 4 dof 
CD2                                        version
CD2                                     
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
c
c
c     forward and back substitution for 4 dof problems

      implicit none

      real*8 b(*),r(*),piv(*)
      integer nop(*),npvt(*)
      integer nb(*),nrhs(*)
      integer nel, idof, indexpiv, jdof, jdof1, jdofm1, jdofp1, ksub1
      integer ksub
      integer nelm1, im2neq, im2nq2, im1, nelp1, i
      integer im1neq, im1nq2, ifin, ijind, ipvtp1
      integer ij, ista, ipvtm1, j
      integer idofmax
      parameter(idofmax=50)
      real*8 dd(idofmax,idofmax), bksub(idofmax), rsub(idofmax)
      real*8 sum(idofmax)
c
c     forward substitution
c need definition for nr1 etc see ilu_4
      nelm1=nel-1
      nelp1=nel+1
      im2neq=0
      im2nq2=0
      do 460 im1=1,nelm1
         indexpiv = im2nq2
         do jdof = 1, idof
            do jdof1 = 1, idof
               indexpiv = indexpiv + 1
               dd(jdof,jdof1) = piv(indexpiv)
            end do
         end do
         do jdof = 1, idof
            bksub(jdof) = r(im1 + nrhs(jdof))
         end do


         do jdof = 2, idof
            jdofm1 = jdof - 1
            do jdof1 = 1, jdofm1
               bksub(jdof) = bksub(jdof) - dd(jdof,jdof1) * bksub(jdof1)
            end do
         end do
         bksub(idof) = bksub(idof) / dd(idof,idof)
         
         do jdof = idof-1,1,-1
            jdofp1 = jdof + 1
            do jdof1 = jdofp1,idof
               bksub(jdof) = bksub(jdof) - dd(jdof,jdof1) * bksub(jdof1)
            end do
            bksub(jdof) = bksub(jdof) / dd(jdof,jdof)
         end do
         do jdof = 1, idof
            r(im1+nrhs(jdof)) = bksub(jdof)
         end do
        im2neq=im2neq+idof
        im2nq2=im2nq2+idof*idof
        i=im1+1
        ista=nop(i)+1
        ipvtm1=npvt(i)-1
c       factor row i
         do jdof = 1, idof
            sum(jdof) = 0.0d00
         end do
        do 440 ij=ista,ipvtm1
          j=nop(ij)
          ijind=ij-nelp1
          do jdof = 1, idof
             rsub(jdof) = r(j + nrhs(jdof))
          end do
          ksub1 = -idof
          do jdof = 1, idof
             ksub1 = ksub1 + idof
             ksub = ksub1
             do jdof1 = 1, idof
                ksub = ksub + 1
                sum(jdof) = sum(jdof) + 
     2               b(ijind+nb(ksub)) * rsub(jdof1)
             end do
          end do
  440   continue
        do jdof = 1, idof
           r(i+nrhs(jdof)) = r(i+nrhs(jdof)) - sum(jdof)
        end do
  460 continue
      im1neq=im2neq
      im1nq2=im2nq2
      indexpiv = im1nq2
      do jdof = 1, idof
         do jdof1 = 1, idof
            indexpiv = indexpiv + 1
            dd(jdof,jdof1) = piv(indexpiv)
         end do
      end do
      do jdof = 1, idof
         bksub(jdof) = r(nel + nrhs(jdof))
      end do
      do jdof = 2, idof
         jdofm1 = jdof - 1
         do jdof1 = 1, jdofm1
            bksub(jdof) = bksub(jdof) - dd(jdof,jdof1) * bksub(jdof1)
         end do
      end do
      bksub(idof) = bksub(idof) / dd(idof,idof)
      
      do jdof = idof-1,1,-1
         jdofp1 = jdof + 1
         do jdof1 = jdofp1,idof
            bksub(jdof) = bksub(jdof) - dd(jdof,jdof1) * bksub(jdof1)
         end do
         bksub(jdof) = bksub(jdof) / dd(jdof,jdof)
      end do



      do jdof = 1, idof
         r(nel+nrhs(jdof)) = bksub(jdof)
      end do
c     back substitution
      do 560 i=nelm1,1,-1
        ipvtp1=npvt(i)+1
        ifin=nop(i+1)
        do jdof = 1, idof
           sum(jdof) = 0.0d00
        end do
        do 540 ij=ipvtp1,ifin
          j=nop(ij)
          ijind=ij-nelp1
          do jdof = 1, idof
             rsub(jdof) = r(j+nrhs(jdof))
          end do
          ksub1 = -idof
          do jdof = 1, idof
             ksub1 = ksub1 + idof
             ksub = ksub1
             do jdof1 = 1, idof
                ksub = ksub + 1
                sum(jdof) = sum(jdof) + 
     2               b(ijind+nb(ksub)) * rsub(jdof1)
             end do
          end do
  540   continue
        do jdof = 1, idof
           r(i+nrhs(jdof)) = r(i+nrhs(jdof)) - sum(jdof)
        end do
  560 continue
      return
      end                      
