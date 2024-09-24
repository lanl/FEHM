      subroutine sub_gmres1(neq,a,b,r,ncon,nop,north,sorthm,epn
     *     ,irb,iirb,npvt,rw,dum,xx,xtemp,piv
     *     ,h,c,s,g,y,iter,idof,iptty,maxor)
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
CD1 acceleration (one degree of freedom).
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_gmres1.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:08   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Tue Sep 17 09:32:48 1996   robinson
CD2 Changed spelling of orthogalization
CD2 
CD2    Rev 1.7   Mon Sep 16 12:11:24 1996   robinson
CD2 Oh, those nasty prologs
CD2 
CD2    Rev 1.6   Thu Sep 12 08:26:18 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.5   Fri May 31 15:31:00 1996   gaz
CD2 correction for more general dof reorgering
CD2 
CD2    Rev 1.4   Tue May 14 14:32:34 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.3   Fri Feb 02 12:27:08 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   05/11/94 16:19:00   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:44   pvcs
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
CD3 temp         real*8   I      Scratch storage array
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
CD4 sub_bksub1      N/A       Performs back-substitution, 1 dof
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
CD5 ad           real*8      Current term of solution matrix
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
CD5 jm1neq       integer     (j-1)*neq
CD5 im1neq       integer     (i-1)*neq
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
CPS BEGIN sub_gmres1
CPS 
CPS Define integer parameters
CPS 
CPS equal_array - set scratch storage arrays to current estimate of...
CPS ... solution
CPS 
CPS constant_value - initialize scratch storage to zero
CPS renumber_array - renumber the solution vector
CPS sub_bksub1 - obtain new pre-conditioned solution
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
CPS     sub_bksub1 - compute new pre-conditioned solution
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
CPS   sub_bksb1 - obtain new pre-conditioned solution
CPS 
CPS ENDLOOP
CPS 
CPS FOR each equation
CPS   Overwrite right hand side array with solution array
CPS ENDFOR
CPS 
CPS END sub_gmres1
CPS
C**********************************************************************
c
c
c acceleration routine

      implicit none

      integer maxor
      real*8 xx(*),a(*),b(*),r(*),sorthm(*),dum(*),rw(*)
      integer ncon(*),nop(*)
      real*8 xtemp(*),piv(*)
      integer irb(*),iirb(*),npvt(*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*), epn
      integer nrhs(1), neq, north, iter, idof, iptty
      real*8 tols, rnorm, sum, ad, vhatnm, tmp1, tmp2, sqroot
      integer neqm1, neqp1, maxit, inorth, kk, i, j, j1, j2, jj, kb
      integer k, jneq, jm1, jm1neq, im1neq
      integer nrhs_dum(1)
      logical done
      parameter (tols=1.d-15)
c     
c     define some parameters
c     
      neqm1=neq-1
      neqp1=neq+1
      nrhs_dum(1)=0
      maxit=iter
      iter=1
      nrhs(1)=0
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
      call sub_bksub1(neq,b,r,nop
     *     ,npvt,piv )        
c     
c     
*--------------------------------------------------------
*     gmres - one degree of freedom
*     
*     written by donn r. hines  11-1-87
*--------------------------------------------------------
 1650 done = .false.
      call renumber_array(dum,r,irb,neq,idof,nrhs_dum,nrhs)
      inorth = north
c     zero out y
      do kk=1,north
         y(kk)=0.0d00
      enddo
      
*     ----------------------
*     find norm of residual
*     ----------------------
      call residual(neq,r,1,rnorm,idof,nrhs)
c     
      g(1) = rnorm
      do 1720 i=2,north+1
         g(i)=0.0d00
 1720 continue
      
*     -------------------------------
*     compute v1 in natural order
*     -------------------------------
      do 1800 i=1,neq
         sorthm(i) = dum(i)/max(rnorm,1.d-20)
 1800 continue
      
      do 3500 j = 1,north
         jm1neq = (j-1)*neq
         
*     -----------
*     compute av
*     -----------
         do 2100 i = 1,neq
            sum = 0.0d00
c     
            j1 = ncon(i) + 1
            j2 = ncon(i+1)
            do 2000 jj = j1,j2
               kb = ncon(jj)
               ad = a(jj - neqp1)
               sum = sum + ad * sorthm(jm1neq + kb)
 2000       continue
            dum(i) = sum
 2100    continue
         
*     -----------
*     reorder av
*     -----------
         call renumber_array(rw,dum,iirb,neq,idof,nrhs_dum,nrhs_dum)
c     
c     forward and back sub for new solution
c     
         call sub_bksub1(neq,b,rw,nop
     *        ,npvt,piv)         
c     
*     ----------------------------------
*     put (m-inverse)av in natural order
*     ----------------------------------
         call renumber_array(dum,rw,irb,neq,idof,nrhs_dum,nrhs_dum)
         
*     -----------------------
*     compute inner products
*     -----------------------
         do 2900 i = 1,j
            im1neq = (i-1)*neq
            h(i,j) = 0.0d00
            do 2800 k = 1,neq
c     
               h(i,j) = h(i,j) + dum(k)*sorthm(im1neq + k)
 2800       continue
 2900    continue
         
*     --------------
*     compute v-hat
*     --------------
         do 3100 i = 1,j
            im1neq = (i-1)*neq
            do 3000 k = 1,neq
c     
               dum(k) = dum(k) - h(i,j)*sorthm(im1neq + k)
 3000       continue
 3100    continue
         
*     -----------------
*     compute ||v-hat||
*     -----------------
         call residual(neq,dum,1,vhatnm,idof,nrhs_dum)
         h(j+1,j) = vhatnm
         vhatnm=max(vhatnm,tols)
         
*     ---------------
*     compute next v
*     ---------------
         jneq = j*neq
         do 3300 i = 1,neq
c     
            sorthm(jneq + i) = dum(i)/vhatnm
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
         sqroot =max(sqroot,tols)
         c(j) = h(j,j)/sqroot
         s(j) = -h(j+1,j)/sqroot
         h(j,j) = sqroot
         h(j+1,j) = 0
         tmp1 = c(j)*g(j) - s(j)*g(j+1)
         tmp2 = s(j)*g(j) + c(j)*g(j+1)
         g(j) = tmp1
         g(j+1) = tmp2
         iter=iter+1
         if (abs(g(j+1)) .le. epn) then
            done = .true.
            inorth = j
            go to 3550
         endif
         
         if(iter.gt.maxit) then
            if(iptty .ne. 0 ) then
               write(iptty,34)
 34            format(/,1x,'Warning issued by 1 degree of freedom '
     &              ,'solver:')
               write(iptty,35) maxit
 35            format(2x,'Maximum number of iterations (',i4,') '
     &              ,'exceeded.')
               write(iptty,36) tmp2,epn
 36            format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &              e14.6,/)
            end if
            done=.true.
            go to 3550
         endif
         
 3500 continue
      
*     -------------------------
*     form approximate solution
*     -------------------------
 3550 do 3700 i = inorth,1,-1
         sum = 0.0d00
         do 3272 j = i+1,inorth
            sum = sum + y(j)*h(i,j)
 3272    continue
         y(i) = (g(i) - sum)/h(i,i)
 3700 continue
      
      do 3900 i = 1,neq
         sum = 0.0d00
c     
         do 3800 j = 1,inorth
            jm1neq = (j-1)*neq
            sum = sum + sorthm(jm1neq+i)*y(j)
 3800    continue
         dum(i) = sum
 3900 continue
      
      do 4000 i = 1,neq
c     
         xx(i) = xx(i) + dum(i)
 4000 continue
      
      if (done) goto 5000
      
**************************************************
*     ---------------------------------------
*     updated solution (xx) is now available.
*     form  axx.
*     ---------------------------------------
      
      do 4200 i = 1,neq
         sum = 0.0d00
c     
         j1 = ncon(i) + 1
         j2 = ncon(i+1)
         do 4100 jj = j1,j2
            kb = ncon(jj)
            ad = a(jj - neqp1)
            sum = sum + ad*xx(kb)
 4100    continue
         dum(i) = sum
 4200 continue
      
*     -----------------------------------------
*     form r = xtemp - axx (in reordered form)
*     ------------------------------------------
      do 4300 i = 1,neq
         r(iirb(i)) = xtemp(i) - dum(i)
 4300 continue
      
*     --------------------------
*     forward-back substitution
*     --------------------------
      call sub_bksub1(neq,b,r,nop
     *     ,npvt,piv)         
c     
      go to 1650
 5000 do 5050 i=1,neq
         r(i)=xx(i)
 5050 continue
      return
      end
         
