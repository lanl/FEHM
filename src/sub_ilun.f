      subroutine  sub_ilun(nel,idof,a,b,na,nb,ncon,nop
     *     ,irb,iirb,npvt,dum,piv)
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
CD1 To perform incomplete lu factorization for n degree of freedom
CD1 problem.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sub_ilun.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Thu Sep 12 08:26:54 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.1   Fri Feb 02 14:39:36 1996   hend
CD2 Updated Prolog and Log
CD2
CD2 4-6-94       G. Zyvoloski   97      Initial implementation
CD2 6-24-94      B. Robinson    97      Made 4 dof routine into 6 dof
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
CPS BEGIN sub_ilun
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
CPS END sub_ilun
CPS
C**********************************************************************
c
c n degree of freedom solver
c
c     n degree of freedom ilu preconditioner
     
      implicit none

      integer nel, idof
      integer ncon(*),nop(*),irb(*),iirb(*),npvt(*)
      real*8 a(*),b(*),piv(*)
      real*8 dum(nel,idof,idof)
      integer na(*),nb(*)                    
      integer nelp1
      integer i1, i2, i3, i4
      integer i, ifin, ij, ijind, ikb, ka, kb
      integer im1nq2, ipvt, ista, ipvtp1, ipvtm1, j, ik, k, kj, l
      integer j2, kpvt, ikind, kjind
      integer ipvind, ip
      integer jdof, jdof1, jdof2, jdof3, jdof4, ksub, ksub1, index1
      integer index2, index3
      integer idofmax
      parameter(idofmax=50)
      real*8 dd(idofmax,idofmax), bksub(idofmax)
c       
      nelp1=nel+1
c
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
            do jdof = 1, idof
               do jdof1 = 1, idof
                  dum(ikb, jdof1, jdof) = 0.0d00
               end do
            end do
          enddo
          do k=i3,i4
            ka=ncon(k)
            kb=iirb(ka)
            index1 = k-nelp1
            ksub1 = -idof
            do jdof = 1, idof
               ksub1 = ksub1 + idof
               ksub = ksub1
               do jdof1 = 1, idof
                  ksub = ksub + 1
                  index2 = index1 + na(ksub)
                  dum(kb, jdof1, jdof) = a(index2)
               end do
            end do
          enddo
          do l=i1,i2
            kb=nop(l)
            index1 = l-nelp1
            ksub1 = -idof
            do jdof = 1, idof
               ksub1 = ksub1 + idof
               ksub = ksub1
               do jdof1 = 1, idof
                  ksub = ksub + 1
                  index2 = index1 + nb(ksub)
                  b(index2) = dum(kb, jdof1, jdof)
               end do
            end do
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
              do jdof = 1, idof
                 index1 = (jdof-1)*idof
                 do jdof1 = 1, idof
                    index1 = index1 + 1
                    index2 = (jdof-1)*idof
                    index3 = jdof1 - idof
                    do jdof2 = 1, idof
                       index2 = index2 + 1
                       index3 = index3 + idof
                       b(ijind+nb(index1)) = b(ijind+nb(index1)) 
     2                      - b(ikind+nb(index2))*b(kjind+nb(index3))
                    end do
                 end do
              end do

            endif
  260     continue
  270   continue
c       calculate the inverse of l(i,i)
        ipvind=ipvt-nelp1
        index1 = 0
        do jdof = 1, idof
           do jdof1 = 1, idof
              index1 = index1 + 1
              dd(jdof,jdof1) = b(ipvind+nb(index1))
           end do
        end do
        do jdof = 1, idof - 1
           do jdof1 = jdof+1, idof
              dd(jdof1,jdof) = dd(jdof1,jdof)/dd(jdof,jdof)
           end do
           do jdof1 = 2, idof
              do jdof2 = 1, min(jdof1,jdof+1) - 1
                 dd(jdof1,jdof+1) = dd(jdof1,jdof+1)
     2                - dd(jdof1,jdof2)*dd(jdof2,jdof+1)
              end do
           end do
        end do


        index1 = im1nq2
        do jdof = 1, idof
           do jdof1 = 1, idof
              index1 = index1 + 1
              piv(index1) = dd(jdof, jdof1)
           end do
        end do
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

              do jdof = 1, idof
                 index1 = (jdof-1)*idof
                 do jdof1 = 1, idof
                    index1 = index1 + 1
                    index2 = (jdof-1)*idof
                    index3 = jdof1 - idof
                    do jdof2 = 1, idof
                       index2 = index2 + 1
                       index3 = index3 + idof
                       b(ijind+nb(index1)) = b(ijind+nb(index1)) 
     2                      -b(ikind+nb(index2))*b(kjind+nb(index3))
                    end do
                 end do
              end do



            endif
  340     continue
c         u(i,j)=u(i,j)/l(i,i)

          do jdof3 = 1, idof



             index1 = jdof3 - idof
             do jdof4 = 1, idof
                index1 = index1 + idof
                bksub(jdof4)=b(ijind+nb(index1))
             end do


             do jdof = 2, idof
                do jdof1 = 1, jdof-1
                   bksub(jdof) = bksub(jdof) -
     2                  dd(jdof,jdof1)*bksub(jdof1)
                end do
             end do
             bksub(idof) = bksub(idof) / dd(idof,idof)
             do jdof = idof-1, 1, -1
                do jdof1 = jdof+1, idof
                   bksub(jdof) = bksub(jdof) -
     2                  dd(jdof,jdof1)*bksub(jdof1)
                end do
                bksub(jdof) = bksub(jdof) / dd(jdof,jdof)
             end do

             index1 = jdof3 - idof
             do jdof4 = 1, idof
                index1 = index1 + idof
                b(ijind+nb(index1))=bksub(jdof4)
             end do

          end do




  350   continue
c        im1nq2=im1nq2+36
        im1nq2=im1nq2+idof*idof
  360 continue
      return
      end
