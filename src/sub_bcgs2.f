      subroutine sub_bcgs2(neq,a,b,r,na,nb,nrhs
     *     ,ncon,nop,north,sorthm,epn
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
CD1 To compute solution of linear set of equation by BCGS
CD1 acceleration (two degrees of freedom).
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
CD2 $Log:   /pvcs.config/fehm90/src/sub_bcgs2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:18:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Sep 12 08:25:50 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.3   Fri May 31 15:40:10 1996   gaz
CD2 correction for more general dof reordering
CD2 
CD2    Rev 1.2   Tue May 14 14:32:28 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.1   Fri Feb 02 12:18:26 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.0   Tue Jan 09 14:14:38 1996   llt
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
CD3 neq          integer  I      Number of entries in the array for
CD3                                  each degree of freedom
CD3 a            real*8   I      Solution matrix
CD3 b            real*8   O      lu factorization matrix
CD3 r            real*8   I/O    Input - right hand side, Output -
CD3                                 solution array
CD3 na           integer  I      Pointer array in a matrix
CD3 nb           integer  I      Pointer array in b matrix
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
CD4                               an a second array
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
CD5 rnorm        real*8      Square root of sum of the sqares of the
CD5                              residuals
CD5 sum          real*8      Intermediate term in calculation
CD5 ad           real*8      Current term of solution matrix
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
c bcgs acceleration routine

      implicit none

      integer maxor
      real*8 xx(*),a(*),b(*),r(*),sorthm(*),dum(*),rw(*)
      integer ncon(*),nop(*)
      real*8 xtemp(*),piv(*)
      integer irb(*),iirb(*),npvt(*), na(*), nb(*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*), epn
      integer nrhs(*), neq, north, iter, idof, iptty
      real*8 tols, sum1, sum2, beta_n          
      real*8 ad1,ad2,ad3,ad4
      integer neqm1, neqp1, maxit, i, j1, j2, jj, kb, neq2, neq3
      integer lu1,lu2,lu3
      integer nrhs_dum(2)
      parameter (tols=1.d-12)
c     
c     define some parameters
c     
      neqm1=neq-1
      neqp1=neq+1
      neq2=neq+neq
      neq3=neq2+neq
      nrhs_dum(1)=0
      nrhs_dum(2)=neq
      lu1=neq+neq
      lu2=lu1+lu1
      lu3=lu2+lu1
      maxit=iter
      iter=1
c     
c     set some arrays equal
c 
c  residual(r)=rw
c  original residual=dum
c  u=sorthm
c  q=sorthm(+neq)    
c  p=sorthm(+neq2)    
c  p-1=sorthm(+neq3)    
c  v=xtemp
c  x=xx
c  iter=n
c  dum=r in natural order
c
c add lu1,lu2,lu3 to sorthm to get 2nd dof space
c add neq to rw,dum,xx,xtemp to get 2nd dof space
c        
c     set some arrays equal
c        
      call equal_array(rw,r,neq,idof,nrhs_dum,nrhs)
c        
c     zero out some storage arrays
c    
      call constant_value(xx,0.0d00,neq,idof,nrhs_dum)
      call constant_value(sorthm(lu2+1),0.0d00,neq,idof,nrhs_dum)
      call constant_value(sorthm(lu3+1),0.0d00,neq,idof,nrhs_dum)
c    
      g(1)=1.0d00
c
c     renumber the solution vector
c
      call renumber_array(r,rw,iirb,neq,idof,nrhs,nrhs_dum)
c
c     forward and back sub for new solution
c
      call sub_bksub2(neq,b,r,nb,nrhs,nop
     *     ,npvt,piv)
      call renumber_array(rw,r,irb,neq,idof,nrhs_dum,nrhs)
c
c save copy of original residual(dum)
c
      call equal_array(dum,rw,neq,idof,nrhs_dum,nrhs_dum)
      call equal_array(sorthm,rw,neq,idof,nrhs_dum,nrhs_dum)
      call equal_array(sorthm(lu2+1),rw,neq,idof,nrhs_dum,nrhs_dum)
      call equal_array(sorthm(lu3+1),rw,neq,idof,nrhs_dum,nrhs_dum)
c
*--------------------------------------------------------
*     bcgs - two degree of freedom
*     
*     written by george zyvoloski 9/25/95
*--------------------------------------------------------
      do iter=1,maxit                                         
c
c loop for iterations
c
      g(2)=0.0d00
      do i=1,lu1
       g(2)=g(2)+dum(i)*rw(i)
      enddo
      beta_n=g(2)/g(1)
c set old norm
      g(1)=g(2)
*     -------------------------------
*     compute u(sorthm) and p(sorthm(+neq2)) in natural order
*     -------------------------------
c <<<<<<<<< skip for iter=1 >>>>>>>>>>>>>>>>>>>>>
      if(iter.gt.1) then
        do i =1,lu1
         sorthm(i)=rw(i)+beta_n*sorthm(i+lu1)
         sorthm(i+lu2)=sorthm(i)+beta_n
     $     *(sorthm(i+lu1)+beta_n*sorthm(i+lu3))
        enddo
*     -------------------------------
*     compute p(n-1)(sorthm(i+neq3))                              
*     -------------------------------
       do i=1,lu1
        sorthm(i+lu3)=sorthm(i+lu2)
       enddo
      endif
c <<<<<<<<< skip for iter=1 >>>>>>>>>>>>>>>>>>>>>
*     ----------------------------------
*     calculate A*p(overwrite v(xtemp) with product)     
*     ----------------------------------

         do i=1,neq
          sum1 = 0.0d00     
          sum2 = 0.0d00     
            j1 = ncon(i) + 1
            j2 = ncon(i+1)
            do jj = j1,j2
             kb = ncon(jj)
             ad1 = a(jj - neqp1 + na(1))
             ad2 = a(jj - neqp1 + na(2))
             ad3 = a(jj - neqp1 + na(3))
             ad4 = a(jj - neqp1 + na(4))
             sum1 = sum1 + ad1*sorthm(lu2+kb) + ad2*sorthm(lu2+neq+kb)
             sum2 = sum2 + ad3*sorthm(lu2+kb) + ad4*sorthm(lu2+neq+kb)
            enddo
            xtemp(i) = sum1
            xtemp(i+neq) = sum2
         enddo         
*     -------------------------------
*     compute v(xtemp)         
*     -------------------------------
      call renumber_array(r,xtemp,iirb,neq,idof,nrhs,nrhs_dum)
      call sub_bksub2(neq,b,r,nb,nrhs,nop
     *     ,npvt,piv)  
      call renumber_array(xtemp,r,irb,neq,idof,nrhs_dum,nrhs)
c
c calculate omega=g(3)
c        
         g(3)=0.0d00
         do i=1,lu1
           g(3)=g(3)+dum(i)*xtemp(i)               
         enddo  
c
c calculate alpha=g(4)               
c
         g(4)=g(2)/g(3)
c
c update q 
c
         do i=1,lu1
          sorthm(i+lu1)=sorthm(i)-g(4)*xtemp(i)                 
          xtemp(i)=  g(4)*(sorthm(i)+sorthm(i+lu1))
         enddo 
c      
         
*     ----------------------------------
*     calculate A*v(overwrite u(sorthm) with product)     
*     ----------------------------------

         do i=1,neq
          sum1 = 0.0d00     
          sum2 = 0.0d00     
            j1 = ncon(i) + 1
            j2 = ncon(i+1)
            do jj = j1,j2
             kb = ncon(jj)
             ad1 = a(jj - neqp1 + na(1))
             ad2 = a(jj - neqp1 + na(2))
             ad3 = a(jj - neqp1 + na(3))
             ad4 = a(jj - neqp1 + na(4))
             sum1 = sum1 + ad1*xtemp(kb) + ad2*xtemp(neq+kb)
             sum2 = sum2 + ad3*xtemp(kb) + ad4*xtemp(neq+kb)
            enddo
            sorthm(i) = sum1
            sorthm(i+neq) = sum2
         enddo         
*     ----------------------------------
*     calculate Pl*A*v                       
*     ----------------------------------
c
c      back sub for new solution
c     
      call renumber_array(r,sorthm,iirb,neq,idof,nrhs,nrhs_dum)
c           
      call sub_bksub2(neq,b,r,nb,nrhs,nop
     *     ,npvt,piv)  
c
         call renumber_array(sorthm,r,irb,neq,idof,nrhs_dum,nrhs)
       do i=1,lu1
        rw(i)=rw(i)-sorthm(i)
       enddo
*     ----------------------
*     update solution xx             
*     ----------------------
       do i=1,lu1
        xx(i)=xx(i)+xtemp(i)
       enddo
*     ----------------------
*     find sum squared of residual
*     ----------------------
c
      call residual(neq,rw,1,g(2),idof,nrhs_dum)
c
c check for convergence
c
      if(g(2).le.epn) go to 5000

      enddo
      if (iter .gt. maxit) then
         if(iptty .ne. 0 ) then
            write(iptty,34)
 34         format(/,1x,'Warning issued by 2 degree of freedom '
     &           ,'solver:')
            write(iptty,35) maxit
 35         format(2x,'Maximum number of iterations (',i4,') '
     &           ,'exceeded.')
            write(iptty,36) g(2),epn
 36         format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &           e14.6,/)
         endif
      end if
 5000 do i=1,neq
         r(iirb(i)+nrhs(1))=xx(i+nrhs_dum(1))
         r(iirb(i)+nrhs(2))=xx(i+nrhs_dum(2))
      enddo
      return
      end
