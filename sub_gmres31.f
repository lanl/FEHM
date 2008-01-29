      subroutine sub_gmres31(neq,a,b,r,na,nb,nrhs,ncon,nop
     *     ,north,sorthm,epn
     *     ,irb,iirb,npvt,rw,xtemp,dum,xx,piv
     *     ,h,c,s,g,y,iter,idof,iptty,maxor,icoupl,tollr,overf)
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
CD1 To compute solution of linear set of equations by GMRES
CD1 acceleration (three degrees of freedom with rdof).
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_gmres31.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:10   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Tue Sep 17 09:32:58 1996   robinson
CD2 Changed spelling of orthogalization
CD2 
CD2    Rev 1.10   Mon Sep 16 12:11:40 1996   robinson
CD2 Oh, those nasty prologs
CD2 
CD2    Rev 1.9   Thu Sep 12 08:26:24 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.8   Fri May 31 15:36:54 1996   gaz
CD2 correction for more general dof reordering
CD2 
CD2    Rev 1.7   Tue May 14 14:32:42 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.6   Fri Mar 01 15:11:06 1996   gaz
CD2 added inorth = j-1 near done=.true.
CD2 
CD2    Rev 1.5   Fri Feb 02 12:29:22 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   05/12/94 10:06:08   llt
CD2 corrected pvcs log info
CD2 
CD2    Rev 1.3   05/11/94 16:19:06   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   04/04/94 14:25:52   robinson
CD2 Declared variable done as logical.
CD2
CD2    Rev 1.1   03/18/94 16:05:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:50   pvcs
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
CD3 r            real*8   I/O    Input - right hand side, Output -
CD3                                 solution array
CD3 na           integer  I      Pointer array in a matrix
CD3 nb           integer  I      Pointer array in b matrix
CD3 nrhs         integer  I      Pointer array for solution vector
CD3 ncon         integer  I      Connectivity matrix for solution matrix
CD3 nop          integer  I      Connectivity matrix for factorization
CD3                                 matrix
CD3 north        integer  I      Number of orthogonalizations
CD3 sorthm       real*8   I      GMRES storage array
CD3 epn          real*8   I      Tolerance for solution
CD3 irb          integer  I      Renumber array - inew=irb(iold)
CD3 iirb         integer  I      Renumber array - iold=iirb(inew)
CD3 npvt         integer  I      Pivot positions in nop
CD3 rw           real*8   I      Scratch storage array
CD3 dum          real*8   I      Scratch storage array
CD3 xx           real*8   I      Scratch storage array
CD3 xtemp        real*8   I      Scratch storage array
CD3 piv          real*8   O      Array of pivots
CD3 h            real*8   I      Scratch storage array
CD3 c            real*8   I      Scratch storage array
CD3 s            real*8   I      Scratch storage array
CD3 g            real*8   I      Scratch storage array
CD3 y            real*8   I      Scratch storage array
CD3 iter         integer  I/O    Input - maximum number of iterations,
CD3                                 Output - number of iterations
CD3                                 performed
CD3 idof         integer  I      Number of degrees of freedom
CD3 iptty        integer  I      Unit number for warning message
CD3 maxor        integer  I      Maximum number of orthogonalizations
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
CD4 Identifier      Type      Description
CD4 
CD4 equal_array     N/A       Sets values in an array equal to those
CD4                               in a second array
CD4 renumber_array  N/A       Performs renumbering of array
CD4 constant_value  N/A       Sets values in an array to a given value
CD4 sub_bksub31     N/A       Performs back-substitution, 3 dof
CD4 residual        N/A       Compute l2 norm of residual array
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 Identifier   Type        Description
CD5 
CD5 tols         real*8      Minimum value for norm of solution
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 rnorm        real*8      Square root of sum of the squares of the
CD5                              residuals
CD5 sum          real*8      Intermediate term in calculation
CD5 sum1         real*8      Intermediate term in calculation
CD5 sum2         real*8      Intermediate term in calculation
CD5 sum3         real*8      Intermediate term in calculation
CD5 dum1         real*8      Intermediate term in calculation
CD5 dum2         real*8      Intermediate term in calculation
CD5 dum3         real*8      Intermediate term in calculation
CD5 ad1          real*8      Current term of solution matrix
CD5 ad2          real*8      Current term of solution matrix
CD5 ad3          real*8      Current term of solution matrix
CD5 ad4          real*8      Current term of solution matrix
CD5 ad5          real*8      Current term of solution matrix
CD5 ad6          real*8      Current term of solution matrix
CD5 ad7          real*8      Current term of solution matrix
CD5 ad8          real*8      Current term of solution matrix
CD5 ad9          real*8      Current term of solution matrix
CD5 vhatnm       real*8      Current estimate of rnorm
CD5 tmp1         real*8      Intermediate term in q-r factorization
CD5 tmp2         real*8      Intermediate term in q-r factorization
CD5 sqroot       real*8      square root of vhatnm
CD5 neqm1        integer     neq-1
CD5 neqp1        integer     neq+1
CD5 maxit        integer     Maximum number of iterations allowed
CD5 inorth       integer     Current orthogonalization number
CD5 kk           integer     Do loop index
CD5 i            integer     Do loop index
CD5 j            integer     Do loop index
CD5 j1           integer     Do loop index
CD5 j2           integer     Do loop index
CD5 jj           integer     Do loop index
CD5 kb           integer     Position in solution matrix
CD5 k            integer     Do loop index
CD5 jneq         integer     Position in GMRES storage array
CD5 jm1          integer     j-1
CD5 nr1          integer     Pointer for solution vector
CD5 nr2          integer     Pointer for solution vector
CD5 nr3          integer     Pointer for solution vector
CD5 lu2          integer     Integer index
CD5 lu3          integer     Integer index
CD5 ke1          integer     Integer index
CD5 ke2          integer     Integer index
CD5 ke3          integer     Integer index
CD5 k1           integer     Integer index
CD5 k2           integer     Integer index
CD5 k3           integer     Integer index
CD5 done         logical     Flag denoting if we are done
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
CD8    3.2.2 Perform Orthogonalization Calculation
CD8    3.2.4 Check for Convergence
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
CPS BEGIN sub_gmres3
CPS 
CPS Define integer parameters
CPS 
CPS equal_array - set scratch storage arrays to current estimate of...
CPS ... solution
CPS 
CPS constant_value - initialize scratch storage to zero
CPS renumber_array - renumber the solution vector
CPS sub_bksub3 - obtain new pre-conditioned solution
CPS 
CPS LOOP to obtain solution to equation set
CPS 
CPS   renumber_array - renumber solution array
CPS   
CPS   FOR each orthogonalization
CPS     Initialize scratch storage to 0
CPS   ENDFOR
CPS   
CPS   residual - compute current l2 norm
CPS   
CPS   Set first value in storage to l2 norm
CPS   
CPS   FOR each orthogonalization
CPS     Initialize value in scratch storage to 0
CPS   ENDFOR
CPS   
CPS   FOR each equation
CPS     Store current residual value in GMRES  storage array (sorthm)
CPS   ENDFOR
CPS   
CPS   FOR each orthogonalization
CPS   
CPS     FOR each equation
CPS       FOR each nonzero element in solution matrix
CPS         Form linear combination of previous solutions
CPS       ENDFOR
CPS     ENDFOR
CPS     
CPS     renumber_array - renumber solution array
CPS     sub_bksub3 - compute new pre-conditioned solution
CPS     renumber_array - renumber solution array
CPS     
CPS     FOR each orthogonalization performed so far
CPS       FOR each equation
CPS         Compute inner product
CPS       ENDFOR
CPS     ENDFOR
CPS     
CPS     FOR each orthogonalization performed so far
CPS       FOR each equation
CPS         Modify current estimate of solution
CPS       ENDFOR
CPS     ENDFOR
CPS     
CPS     residual - compute l2 norm of modified current estimate of...
CPS     ... solution
CPS     
CPS     FOR each equation
CPS       Store current residual value in GMRES  storage array (sorthm)
CPS     ENDFOR
CPS     
CPS     FOR orthogonalization performed so far (minus one)
CPS       Calculate terms of q-r factorization matrix
CPS     ENDFOR
CPS     
CPS     Perform intermediate calculations in new estimate of l2 norm
CPS     
CPS     IF projected estimate of l2 norm is less than tolerance
CPS   EXITIF projected estimate of l2 norm is less than tolerance
CPS     ENDIF
CPS     
CPS     IF the maximum number of iterations is reached
CPS       IF we are writing a warning message
CPS         Write warning message
CPS       ENDIF
CPS   EXITIF the maximum number of iterations is reached
CPS     ENDIF
CPS   
CPS   ENDFOR
CPS   
CPS   FOR each orthogonalization from current to first
CPS     For each orthogonalization current to last
CPS       Calculate coefficent for GMRES solution
CPS     ENDFOR
CPS     Calculate coefficent for GMRES solution
CPS   ENDFOR
CPS   
CPS   FOR each equation
CPS     FOR each orthogonalization
CPS       Form linear combination to form correction to solution
CPS     ENDFOR
CPS   ENDFOR
CPS   
CPS   FOR each equation
CPS     Add correction to previous solution estimate
CPS   ENDFOR
CPS   
CPS EXITIF convergence achieved or maximum iterations were taken
CPS 
CPS   FOR each equation
CPS     FOR each nonzero term in solution matrix
CPS       Form new estimate of solution
CPS     ENDFOR
CPS   ENDFOR
CPS   
CPS   FOR each equation
CPS     Compute new residual
CPS   ENDFOR
CPS   
CPS   sub_bksb3 - obtain new pre-conditioned solution
CPS 
CPS ENDLOOP
CPS 
CPS FOR each equation
CPS   Overwrite right hand side array with solution array
CPS ENDFOR
CPS 
CPS END sub_gmres3
CPS
C**********************************************************************
c
c
c
c     acceleration routine
    
      implicit none

      integer neq, maxor
      real*8 xx(neq,3),a(*),b(*),r(*),sorthm(*),dum(*),rw(*)
      integer ncon(*),nop(*)
      real*8 xtemp(*),piv(neq,*)
      integer irb(*),iirb(*),npvt(*), na(*), nb(*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*), epn
      integer nrhs(*), north, iter, idof, iptty
      real*8 tols, rnorm, sum, ad1, ad2, ad3, ad4, ad5, ad6, ad7
      real*8 ad8, ad9, vhatnm, tmp1
      real*8 tmp2, sqroot, sum1, sum2, sum3, dum1, dum2, dum3
      integer neqm1, neqp1, maxit, inorth, kk, i, j, j1, j2, jj, kb
      integer k, jneq, jm1, nr1, nr2, nr3, lu2, lu3, ke1, k1, k2
      integer ke2, ke3, k3
      integer nrhs_dum(3),ib
      integer icoupl
      real*8 tollr,overf
      logical done
      real*8,allocatable :: sto(:)
      parameter (tols=1.d-12)
c     
c     define some parameters
c     
      neqm1=neq-1
      neqp1=neq+1
      lu2=(north+1)*neq*2
      lu3=(north+1)*neq
      nr1=nrhs(1)
      nr2=nrhs(2)
      nr3=nrhs(3)
      nrhs_dum(1)=0
      nrhs_dum(2)=neq
      nrhs_dum(3)=neq+neq
      maxit=iter
      iter=1
c     
c     set some arrays equal
c     
      call equal_array(rw,r,neq,idof,nrhs_dum,nrhs)
      call equal_array(xtemp,r,neq,idof,nrhs_dum,nrhs)
c     
c     zero out some storage arrays
c     
      call constant_value(xx,0.0d00,neq,idof,nrhs_dum)
c     
c     renumber the solution vector
c     
      call renumber_array(r,rw,iirb,neq,idof,nrhs,nrhs_dum)
c     
c     forward and back sub for new solution
c     
      call sub_bksub31(neq,a,b,r,na,nrhs,ncon,nop
     *     ,npvt,piv,dum,icoupl,tollr,overf)
c     
*------------------------------------------------------
*     gmres - three degrees of freedom.
*     
*     written by donn r. hines     11-18-87
*------------------------------------------------------
 1650 done = .false.
      inorth = north
      call renumber_array(dum,r,irb,neq,idof,nrhs_dum,nrhs)
c     zero out y
      do kk=1,north
         y(kk)=0.0d00
      enddo
*     ----------------------
*     find norm of residual
*     ----------------------
      call residual(neq,r,1,rnorm,idof,nrhs)
*     ----------
*     compute g
*     ----------
      g(1) = rnorm
      do 1720 i=2,north+1
         g(i)=0.0d00
 1720 continue
      
*     -------------------------------
*     compute v1 in natural order
*     -------------------------------
      do 1800 i=1,neq
         sorthm(i) = dum(i+nrhs_dum(1))/rnorm
         sorthm(i+lu3) = dum(i+nrhs_dum(2))/rnorm
         sorthm(i+lu2) = dum(i+nrhs_dum(3))/rnorm
 1800 continue
      
      do 3500 j = 1,north
         ke1 = (j-1)*neq
         ke2 = ke1+lu3
         ke3 = ke1+lu2
         
*     -----------
*     compute av
*     -----------
         do 2100 i = 1,neq
            sum1 = 0.0d00
            sum2 = 0.0d00
            sum3 = 0.0d00
            j1 = ncon(i) + 1
            j2 = ncon(i+1)
            do 2000 jj = j1,j2
               kb = ncon(jj)
               ad1 = a(jj - neqp1 + na(1))
               ad2 = a(jj - neqp1 + na(2))
               ad3 = a(jj - neqp1 + na(3))
               ad4 = a(jj - neqp1 + na(4))
               ad5 = a(jj - neqp1 + na(5))
               ad6 = a(jj - neqp1 + na(6))
               ad7 = a(jj - neqp1 + na(7))
               ad8 = a(jj - neqp1 + na(8))
               ad9 = a(jj - neqp1 + na(9))
               sum1 = sum1 + ad1*sorthm(ke1+kb) + ad2*sorthm(ke2+kb)
     *              + ad3*sorthm(ke3+kb)
               sum2 = sum2 + ad4*sorthm(ke1+kb) + ad5*sorthm(ke2+kb)
     *              + ad6*sorthm(ke3+kb)
               sum3 = sum3 + ad7*sorthm(ke1+kb) + ad8*sorthm(ke2+kb)
     *              + ad9*sorthm(ke3+kb)
 2000       continue
            dum(i+nrhs_dum(1)) = sum1
            dum(i+nrhs_dum(2)) = sum2
            dum(i+nrhs_dum(3)) = sum3
 2100    continue
         
*     -----------
*     reorder av
*     -----------
         call renumber_array(rw,dum,iirb,neq,idof,nrhs_dum,nrhs_dum)
         
*     -----------------
*     forward-back sub.
*     -----------------
      call sub_bksub31(neq,a,b,rw,na,nrhs_dum,ncon,nop
     *     ,npvt,piv,dum,icoupl,tollr,overf)
         
*     ----------------------------------
*     put (m-inverse)av in natural order
*     ----------------------------------
         call renumber_array(dum,rw,irb,neq,idof,nrhs_dum,nrhs_dum)
         
*     -----------------------
*     compute inner products
*     -----------------------
         do 2900 i = 1,j
            k1 = (i-1)*neq
            k2 = k1 + lu3
            k3 = k1 + lu2
            h(i,j) = 0.0d00
            do 2800 k = 1,neq
               dum1=dum(k+nrhs_dum(1))
               dum2=dum(k+nrhs_dum(2))
               dum3=dum(k+nrhs_dum(3))
               h(i,j)=h(i,j)+dum1*sorthm(k1+k)+dum2*sorthm(k2+k)
     *              +dum3*sorthm(k3+k)
 2800       continue
 2900    continue
         
*     --------------
*     compute v-hat
*     --------------
         do 3100 i = 1,j
            k1 = (i-1)*neq
            k2 = k1 + lu3
            k3 = k1 + lu2
            do 3000 k = 1,neq
             dum(k+nrhs_dum(1)) = dum(k+nrhs_dum(1)) - h(i,j)
     &                           *sorthm(k1 + k)
             dum(k+nrhs_dum(2)) = dum(k+nrhs_dum(2)) - h(i,j)
     &                           *sorthm(k2 + k)
             dum(k+nrhs_dum(3)) = dum(k+nrhs_dum(3)) - h(i,j)
     &                           *sorthm(k3 + k)
 3000       continue
 3100    continue
         
*     -----------------
*     compute ||v-hat||
*     -----------------
         call residual(neq,dum,1,vhatnm,idof,nrhs_dum)
         h(j+1,j) = vhatnm
         vhatnm=max(vhatnm,1.d-12)
         
*     ---------------
*     compute next v
*     ---------------
         jneq = j*neq
         do 3300 i = 1,neq
            sorthm(jneq + i) = dum(i+nrhs_dum(1))/vhatnm
            sorthm(jneq + lu3 + i) = dum(i+nrhs_dum(2))/vhatnm
            sorthm(jneq + lu2 + i) = dum(i+nrhs_dum(3))/vhatnm
 3300    continue
         
*     -----------------
*     q-r factorization
*     -----------------
         jm1 = j - 1
         do 3400 i = 1,jm1
            tmp1 = c(i)*h(i,j) - s(i)*h(i+1,j)
            tmp2 = s(i)*h(i,j) + c(i)*h(i+1,j)
            h(i,j) = tmp1
            h(i+1,j) = tmp2
 3400    continue
         sqroot = sqrt(h(j,j)*h(j,j)+h(j+1,j)*h(j+1,j))
         sqroot =max(sqroot,1.d-12)
         c(j) = h(j,j)/sqroot
         s(j) = -h(j+1,j)/sqroot
         h(j,j) = sqroot
         h(j+1,j) = 0.0d00
         tmp1 = c(j)*g(j) - s(j)*g(j+1)
         tmp2 = s(j)*g(j) + c(j)*g(j+1)
         g(j) = tmp1
         g(j+1) = tmp2
         
         iter=iter+1
         if(iter.gt.maxit) then
            if(iptty .ne. 0 ) then
               write(iptty,34)
 34            format(/,1x,'Warning issued by 3-1 degree of freedom '
     &              ,'solver:')
               write(iptty,35) maxit
 35            format(2x,'Maximum number of iterations (',i4,') '
     &              ,'exceeded.')
               write(iptty,36) tmp2,epn
 36            format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &              e14.6,/)
            end if
            done=.true.
            inorth = j-1
            go to 3550
         endif
         if (abs(g(j+1)) .le. epn) then
            done = .true.
            inorth = j
            go to 3550
         endif
         
 3500 continue
      
*     -------------------------
*     form approximate solution
*     -------------------------
 3550 do 3700 i = inorth,1,-1
         sum = 0.0d00
         do 3600 j = i+1,inorth
            sum = sum + y(j)*h(i,j)
 3600    continue
         y(i) = (g(i) - sum)/h(i,i)
 3700 continue
      
      do 3900 i = 1,neq
         sum1 = 0.0d00
         sum2 = 0.0d00
         sum3 = 0.0d00
         do 3800 j = 1,inorth
            k1=(j-1)*neq
            k2=k1+lu3
            k3=k1+lu2
            sum1 = sum1 + sorthm(k1+i)*y(j)
            sum2 = sum2 + sorthm(k2+i)*y(j)
            sum3 = sum3 + sorthm(k3+i)*y(j)
 3800    continue
         dum(i+nrhs_dum(1)) = sum1
         dum(i+nrhs_dum(2)) = sum2
         dum(i+nrhs_dum(3)) = sum3
 3900 continue
      
      do 4000 i = 1,neq
         xx(i,1) = xx(i,1) + dum(i+nrhs_dum(1))
         xx(i,2) = xx(i,2) + dum(i+nrhs_dum(2))
         xx(i,3) = xx(i,3) + dum(i+nrhs_dum(3))
 4000 continue
      
      if (done) goto 5000
      
**************************************************
*     ---------------------------------------
*     updated solution (xx) is now available.
*     form  axx.
*     ---------------------------------------
      
      do 4200 i = 1,neq
         sum1 = 0.0d00
         sum2 = 0.0d00
         sum3 = 0.0d00
         j1 = ncon(i) + 1
         j2 = ncon(i+1)
         do 4100 jj = j1,j2
            kb = ncon(jj)
            ad1 = a(jj - neqp1 + na(1))
            ad2 = a(jj - neqp1 + na(2))
            ad3 = a(jj - neqp1 + na(3))
            ad4 = a(jj - neqp1 + na(4))
            ad5 = a(jj - neqp1 + na(5))
            ad6 = a(jj - neqp1 + na(6))
            ad7 = a(jj - neqp1 + na(7))
            ad8 = a(jj - neqp1 + na(8))
            ad9 = a(jj - neqp1 + na(9))
            sum1 = sum1 + ad1*xx(kb,1) + ad2*xx(kb,2) + ad3*xx(kb,3)
            sum2 = sum2 + ad4*xx(kb,1) + ad5*xx(kb,2) + ad6*xx(kb,3)
            sum3 = sum3 + ad7*xx(kb,1) + ad8*xx(kb,2) + ad9*xx(kb,3)
 4100    continue
         dum(i+nrhs_dum(1)) = sum1
         dum(i+nrhs_dum(2)) = sum2
         dum(i+nrhs_dum(3)) = sum3
 4200 continue
      
*     -----------------------------------------
*     form r = xtemp - axx (in reordered form)
*     ------------------------------------------
      do 4300 ib = 1,neq
         i=iirb(ib)
         r(i+nrhs(1))=xtemp(ib+nrhs_dum(1))-dum(ib+nrhs_dum(1))
         r(i+nrhs(2))=xtemp(ib+nrhs_dum(2))-dum(ib+nrhs_dum(2))
         r(i+nrhs(3))=xtemp(ib+nrhs_dum(3))-dum(ib+nrhs_dum(3))
 4300 continue
      
*     --------------------------
*     forward-back substitution
*     ---------------------------
      call sub_bksub31(neq,a,b,r,na,nrhs,ncon,nop
     *     ,npvt,piv,dum,icoupl,tollr,overf)
*     
      go to 1650
 5000 do 5050 i=1,neq
         r(i+nrhs(1))=xx(i,1)
         r(i+nrhs(2))=xx(i,2)
         r(i+nrhs(3))=xx(i,3)
 5050 continue
      return
      end
