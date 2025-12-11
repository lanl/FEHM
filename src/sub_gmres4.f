      subroutine sub_gmres4(neq,a,b,r,na,nb,nrhs
     *     ,ncon,nop,north,sorthm,epn
     *     ,irb,iirb,npvt,rw,xx,xtemp,dum,piv
     *     ,h,c,s,g,y,iter,idof,iptty)
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
CD1 acceleration (four degrees of freedom).
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_gmres4.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:10   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:19:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:46   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.10   Tue Sep 17 09:33:02 1996   robinson
CD2 Changed spelling of orthogalization
CD2 
CD2    Rev 1.9   Mon Sep 16 12:11:46 1996   robinson
CD2 Oh, those nasty prologs
CD2 
CD2    Rev 1.8   Mon Sep 16 11:51:26 1996   robinson
CD2 Prolog change
CD2 
CD2    Rev 1.7   Thu Sep 12 08:26:26 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.6   Fri May 31 15:37:22 1996   gaz
CD2 correction for more general dof reordering
CD2 
CD2    Rev 1.5   Tue May 14 14:32:44 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.4   Fri Feb 02 12:30:24 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   11/28/94 14:19:06   llt
CD2 Corrected max function, so would run on ibm.
CD2 
CD2    Rev 1.2   05/11/94 16:19:12   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 16:05:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:52   pvcs
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
CD3 xx           real*8   I      Scratch storage array
CD3 xtemp        real*8   I      Scratch storage array
CD3 dum          real*8   I      Scratch storage array
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
CD4 sub_bksub4      N/A       Performs back-substitution, 4 dof
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
CD5 vhatnm       real*8      Current estimate of rnorm
CD5 tmp1         real*8      Intermediate term in q-r factorization
CD5 s1           real*8      Intermediate term in calculation
CD5 s2           real*8      Intermediate term in calculation
CD5 s3           real*8      Intermediate term in calculation
CD5 s4           real*8      Intermediate term in calculation
CD5 prdinr       real*8      Intermediate term in calculation
CD5 tmp2         real*8      Intermediate term in q-r factorization
CD5 sqroot       real*8      square root of vhatnm
CD5 sum1         real*8      Intermediate term in calculation
CD5 sum2         real*8      Intermediate term in calculation
CD5 sum3         real*8      Intermediate term in calculation
CD5 sum4         real*8      Intermediate term in calculation
CD5 hh           real*8      Intermediate term in calculation
CD5 summ         real*8      Intermediate term in calculation
CD5 neqm1        integer     neq-1
CD5 neqp1        integer     neq+1
CD5 maxit        integer     Maximum number of iterations allowed
CD5 inorth       integer     Current orthogonalization number
CD5 kk           integer     Do loop index
CD5 i            integer     Do loop index
CD5 nr1          integer     Pointer for solution vector
CD5 nr2          integer     Pointer for solution vector
CD5 nr3          integer     Pointer for solution vector
CD5 nr4          integer     Pointer for solution vector
CD5 lu3          integer     Integer index
CD5 k1           integer     Integer index
CD5 k2           integer     Integer index
CD5 k3           integer     Integer index
CD5 k4           integer     Integer index
CD5 jnorm1       integer     Do lopo limit parameter
CD5 nel          integer     Number of equations
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
CD5 jnel1        integer     Integer index
CD5 jnel2        integer     Integer index
CD5 jnel3        integer     Integer index
CD5 jnel4        integer     Integer index
CD5 inor         integer     Do loop index
CD5 jnorth       integer     Number of orthogonalizations
CD5 istaa        integer     Do loop limit parameter
CD5 ifina        integer     Do loop limit parameter
CD5 ija          integer     Do loop index
CD5 ja           integer     Integer index
CD5 nelm1        integer     Number of equations - 1
CD5 knorth       integer     Integer index
CD5 ijaind       integer     Integer index
CD5 nelp1        integer     Number of equations + 1
CD5 jnor         integer     Do loop index
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
CD7 Note that although this routine uses essentially the same algorithm
CD7 as sub_gmres1, sub_gmres2, and sub_gmres3, the actual
CD7 implementation is somewhat different for the sake of computational
CD7 efficiency.  Thus the code structure is somewhat different than
CD7 these other routines.
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
CPS BEGIN sub_gmres4
CPS 
CPS Define integer parameters
CPS 
CPS equal_array - set scratch storage arrays to current estimate of...
CPS ... solution
CPS 
CPS constant_value - initialize scratch storage to zero
CPS renumber_array - renumber the solution vector
CPS sub_bksub4 - obtain new pre-conditioned solution
CPS 
CPS REPEAT to obtain solution
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
CPS   LOOP for computing next orthogonalization
CPS   EXITIF the maximum number of orthogonalizations is reached or...
CPS   ... convergence is achieved
CPS   
CPS     FOR each equation
CPS       FOR each nonzero element in solution matrix
CPS         Form linear combination of previous solutions
CPS       ENDFOR
CPS     ENDFOR
CPS     
CPS     sub_bksub4 - compute new pre-conditioned solution
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
CPS       Set flag to exit after this iteration
CPS     ELSE
CPS     
CPS       IF the maximum number of iterations is reached
CPS         IF we are writing a warning message
CPS           Write warning message
CPS         ENDIF
CPS         Set flag to exit after this iteration
CPS       ELSE
CPS         Increase orthogonalization counter by 1
CPS       ENDIF
CPS     ENDIF
CPS   
CPS   LOOP for computing next orthogonalization
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
CPS   IF we are not done
CPS   
CPS     FOR each equation
CPS       FOR each nonzero term in solution matrix
CPS         Form new estimate of solution
CPS       ENDFOR
CPS     ENDFOR
CPS   
CPS     sub_bksb4 - obtain new pre-conditioned solution
CPS   
CPS   ENDIF
CPS 
CPS UNTIL we are done
CPS 
CPS FOR each equation
CPS   Overwrite right hand side array with solution array
CPS ENDFOR
CPS 
CPS END sub_gmres4
CPS
C**********************************************************************
c
c
c     acceleration routine
c     note dimension of h is 1
     
      implicit none

      integer neq
      real*8 xx(neq,4),a(*),b(*),r(*),sorthm(*),rw(*)
      integer ncon(*),nop(*), na(*), nb(*)
      real*8 xtemp(*),dum(*),piv(*)
      integer irb(*),iirb(*),npvt(*)
      real*8 h(*),c(*),s(*),g(*),y(*), epn
      integer nrhs(*), north, iter, idof, iptty
      real*8 tols, rnorm
      real*8 vhatnm, tmp1, s1, s2, s3, s4, prdinr
      real*8 tmp2, sqroot, sum1, sum2, sum3, sum4, hh, summ
      integer neqm1, neqp1, maxit, inorth, kk, i
      integer nr1, nr2, nr3, lu3, k1, k2, jnorm1
      integer k3, nel
      integer na11, na12, na13, na14
      integer na21, na22, na23, na24
      integer na31, na32, na33, na34
      integer na41, na42, na43, na44
      integer jnel1, jnel2, jnel3, inor, jnel4, jnorth, istaa, ifina
      integer ija, ja, nelm1, nr4, knorth, ijaind, nelp1, k4, jnor
      integer nrhs_dum(4)
      logical done
      parameter (tols=1.d-12)
c     
c     define some parameters
c     
      neqm1=neq-1
      neqp1=neq+1
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
      nel=neq
      nelp1=nel+1
      nelm1=nel-1
      nr1=nrhs(1)
      nr2=nrhs(2)
      nr3=nrhs(3)
      nr4=nrhs(4)
      nrhs_dum(1)=0
      nrhs_dum(2)=neq
      nrhs_dum(3)=neq+neq
      nrhs_dum(4)=neq+neq+neq
      lu3=(north+1)*neq
      jnel1=0
      jnel2=jnel1+lu3
      jnel3=jnel2+lu3
      jnel4=jnel3+lu3
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
      call sub_bksub4(neq,b,r,nb,nrhs,nop
     *     ,npvt,piv)
c     
c     
*------------------------------------------------------
*     gmres - four degrees of freedom.
*     
*     written by donn r. hines     11-16-87
*------------------------------------------------------
c     30 september 1988.  started modifying to neq degrees of freedom.
      done=.false.
 1000 continue
      inorth=north
      call renumber_array(dum,r,irb,neq,idof,nrhs_dum,nrhs)
c     zero out y
      do kk=1,north
         y(kk)=0.0d00
      enddo
c     compute g
      call residual(neq,r,1,rnorm,idof,nrhs)
c     
      g(1)=rnorm
      do 1030 inor=2,north+1
         g(inor)=0.0d00
 1030 continue
c     compute v1 in natural order
      do 1050 i=1,nel
         sorthm(i+jnel1)=dum(i+nrhs_dum(1))/rnorm
         sorthm(i+jnel2)=dum(i+nrhs_dum(2))/rnorm
         sorthm(i+jnel3)=dum(i+nrhs_dum(3))/rnorm
         sorthm(i+jnel4)=dum(i+nrhs_dum(4))/rnorm
 1050 continue
      jnorth=1
 1100 continue
      if ((jnorth.gt.north).or.done) goto 2000
c     compute av
      knorth=(jnorth-1)*neq
      k1=knorth+jnel1
      k2=knorth+jnel2
      k3=knorth+jnel3
      k4=knorth+jnel4
      do 1160 i=1,nel
         sum1=0.0d00
         sum2=0.0d00
         sum3=0.0d00
         sum4=0.0d00
         istaa=ncon(i)+1
         ifina=ncon(i+1)
         do 1140 ija=istaa,ifina
            ja=ncon(ija)
            ijaind=ija-nelp1
            s1=sorthm(ja+k1)
            s2=sorthm(ja+k2)
            s3=sorthm(ja+k3)
            s4=sorthm(ja+k4)
            sum1=sum1
     &           +a(ijaind+na11)*s1
     &           +a(ijaind+na12)*s2
     &           +a(ijaind+na13)*s3
     &           +a(ijaind+na14)*s4
            sum2=sum2
     &           +a(ijaind+na21)*s1
     &           +a(ijaind+na22)*s2
     &           +a(ijaind+na23)*s3
     &           +a(ijaind+na24)*s4
            sum3=sum3
     &           +a(ijaind+na31)*s1
     &           +a(ijaind+na32)*s2
     &           +a(ijaind+na33)*s3
     &           +a(ijaind+na34)*s4
            sum4=sum4
     &           +a(ijaind+na41)*s1
     &           +a(ijaind+na42)*s2
     &           +a(ijaind+na43)*s3
     &           +a(ijaind+na44)*s4
 1140    continue
         dum(i+nrhs_dum(1))=sum1
         dum(i+nrhs_dum(2))=sum2
         dum(i+nrhs_dum(3))=sum3
         dum(i+nrhs_dum(4))=sum4
 1160 continue
*     -----------
*     reorder av
*     -----------
         call renumber_array(rw,dum,iirb,neq,idof,nrhs_dum,nrhs_dum)
         
c     compute m(-1)*av
      call sub_bksub4(neq,b,rw,nb,nrhs_dum,nop
     *     ,npvt,piv)
c
*     ----------------------------------
*     put (m-inverse)av in natural order
*     ----------------------------------
         call renumber_array(dum,rw,irb,neq,idof,nrhs_dum,nrhs_dum)
         
c     
c     compute inner products
c     
      do 1430 inor=1,jnorth
         knorth=(inor-1)*neq
         k1=knorth+jnel1
         k2=knorth+jnel2
         k3=knorth+jnel3
         k4=knorth+jnel4
         prdinr=0.0d00
         do 1420 i=1,nel
            prdinr=prdinr+dum(i+nrhs_dum(1))*sorthm(k1+i)
     &           +dum(i+nrhs_dum(2))*sorthm(k2+i)
     &           +dum(i+nrhs_dum(3))*sorthm(k3+i)
     &           +dum(i+nrhs_dum(4))*sorthm(k4+i)
 1420    continue
         h((inor-1)*north+jnorth)=prdinr
 1430 continue
c     compute v-hat
c     
c     identify new orthogonal term
c     
      do 1630 inor=1,jnorth
         knorth=(inor-1)*neq
         k1=knorth+jnel1
         k2=knorth+jnel2
         k3=knorth+jnel3
         k4=knorth+jnel4
         hh=h((inor-1)*north+jnorth)
         do 1620 i=1,nel
            dum(i+nrhs_dum(1))=dum(i+nrhs_dum(1))-hh*sorthm(i+k1)
            dum(i+nrhs_dum(2))=dum(i+nrhs_dum(2))-hh*sorthm(i+k2)
            dum(i+nrhs_dum(3))=dum(i+nrhs_dum(3))-hh*sorthm(i+k3)
            dum(i+nrhs_dum(4))=dum(i+nrhs_dum(4))-hh*sorthm(i+k4)
 1620    continue
 1630 continue
c     compute ||v-hat||
      call residual(neq,dum,1,vhatnm,idof,nrhs_dum)
c     
      h(jnorth*north+jnorth)=vhatnm
      vhatnm=max(vhatnm,tols)
c     compute next v
      knorth=jnorth*neq
      k1=knorth+jnel1
      k2=knorth+jnel2
      k3=knorth+jnel3
      k4=knorth+jnel4
      do 1820 i=1,nel
         sorthm(k1+i)=dum(i+nrhs_dum(1))/vhatnm
         sorthm(k2+i)=dum(i+nrhs_dum(2))/vhatnm
         sorthm(k3+i)=dum(i+nrhs_dum(3))/vhatnm
         sorthm(k4+i)=dum(i+nrhs_dum(4))/vhatnm
 1820 continue
c     q-r factorization
c     1 october 1988.  the qr factorization is independent of neq.
      jnorm1=jnorth-1
      do 1910 inor=1,jnorm1
         tmp1=h((inor-1)*north+jnorth)
         tmp2=h(inor*north+jnorth)
         h((inor-1)*north+jnorth)=c(inor)*tmp1-s(inor)*tmp2
         h(inor*north+jnorth)=s(inor)*tmp1+c(inor)*tmp2
 1910 continue
      tmp1=h(jnorm1*north+jnorth)
      tmp2=h(jnorth*north+jnorth)
      sqroot=sqrt(tmp1*tmp1+tmp2*tmp2)
      sqroot =max(sqroot,1.0d-12)
      c(jnorth)=tmp1/sqroot
      s(jnorth)=-tmp2/sqroot
      h(jnorm1*north+jnorth)=sqroot
      h(jnorth*north+jnorth)=0.0d00
      tmp1=c(jnorth)*g(jnorth)-s(jnorth)*g(jnorth+1)
      tmp2=s(jnorth)*g(jnorth)+c(jnorth)*g(jnorth+1)
      g(jnorth)=tmp1
      g(jnorth+1)=tmp2
      if (abs(g(jnorth+1)).le.epn) then
         done=.true.
         inorth=jnorth
      else
         iter=iter+1
         if(iter.gt.maxit) then
            if( iptty .ne. 0 ) then
               write(iptty,34)
 34            format(/,1x,'Warning issued by 4 degree of freedom '
     &              ,'solver:')
               write(iptty,35) maxit
 35            format(2x,'Maximum number of iterations (',i4,') '
     &              ,'exceeded.')
               write(iptty,36) tmp2,epn
 36            format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &              e14.6,/)
            end if
            done=.true.
         else
            jnorth=jnorth+1
         endif
      endif
      goto 1100
 2000 continue
c     form approximate solution
      do 2220 inor=inorth,1,-1
         summ=0.0d00
         do 2210 jnor=inor+1,inorth
            summ=summ+y(jnor)*h((inor-1)*north+jnor)
 2210    continue
         y(inor)=(g(inor)-summ)/h((inor-1)*north+inor)
 2220 continue
      do 2350 i=1,nel
         sum1=0.0d00
         sum2=0.0d00
         sum3=0.0d00
         sum4=0.0d00
         do 2330 inor=1,inorth
            knorth=(inor-1)*neq
            k1=knorth+jnel1
            k2=knorth+jnel2
            k3=knorth+jnel3
            k4=knorth+jnel4
            sum1=sum1+sorthm(i+k1)*y(inor)
            sum2=sum2+sorthm(i+k2)*y(inor)
            sum3=sum3+sorthm(i+k3)*y(inor)
            sum4=sum4+sorthm(i+k4)*y(inor)
 2330    continue
         xx(i,1)=xx(i,1)+sum1
         xx(i,2)=xx(i,2)+sum2
         xx(i,3)=xx(i,3)+sum3
         xx(i,4)=xx(i,4)+sum4
 2350 continue
      if (.not.done) then
c     updated solution (xx) is now available form  axx.
c     form r = xtemp - axx (in reordered form)
         do 2460 i=1,nel
            sum1=0.0d00
            sum2=0.0d00
            sum3=0.0d00
            sum4=0.0d00
            istaa=ncon(i)+1
            ifina=ncon(i+1)
            do 2440 ija=istaa,ifina
               ja=ncon(ija)
               ijaind=ija-nelp1
               sum1=sum1
     &              +a(ijaind+na11)*xx(ja,1)
     &              +a(ijaind+na12)*xx(ja,2)
     &              +a(ijaind+na13)*xx(ja,3)
     &              +a(ijaind+na14)*xx(ja,4)
               sum2=sum2
     &              +a(ijaind+na21)*xx(ja,1)
     &              +a(ijaind+na22)*xx(ja,2)
     &              +a(ijaind+na23)*xx(ja,3)
     &              +a(ijaind+na24)*xx(ja,4)
               sum3=sum3
     &              +a(ijaind+na31)*xx(ja,1)
     &              +a(ijaind+na32)*xx(ja,2)
     &              +a(ijaind+na33)*xx(ja,3)
     &              +a(ijaind+na34)*xx(ja,4)
               sum4=sum4
     &              +a(ijaind+na41)*xx(ja,1)
     &              +a(ijaind+na42)*xx(ja,2)
     &              +a(ijaind+na43)*xx(ja,3)
     &              +a(ijaind+na44)*xx(ja,4)
 2440       continue
*     -----------------------------------------
*     form r = xtemp - axx (in reordered form)
*     ------------------------------------------
            r(iirb(i)+nrhs(1))=xtemp(i+nrhs_dum(1))-sum1
            r(iirb(i)+nrhs(2))=xtemp(i+nrhs_dum(2))-sum2
            r(iirb(i)+nrhs(3))=xtemp(i+nrhs_dum(3))-sum3
            r(iirb(i)+nrhs(4))=xtemp(i+nrhs_dum(4))-sum4
 2460    continue
c     compute m(-1)*av
         call sub_bksub4(neq,b,r,nb,nrhs,nop
     *        ,npvt,piv)
c     
      endif
      if (.not.done) goto 1000
      do 3020 i=1,nel
         r(i+nrhs(1))=xx(i,1)           
         r(i+nrhs(2))=xx(i,2)           
         r(i+nrhs(3))=xx(i,3)           
         r(i+nrhs(4))=xx(i,4)           
 3020 continue
c     
      return
      end
