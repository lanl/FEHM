      subroutine sub_bcgstab1(neq,a,b,r,ncon,nop,north,sorthm,epn
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
CD1 To compute solution of linear set of equation by BCGSTAB
CD1 acceleration (one degree of freedom).
CD1
C**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: ?, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/sub_bcgstab1.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:44:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:18:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:32:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:11:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
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
CD5 su2          real*8      Intermediate term in calculation
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
CD8 
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
c
c bcgstab acceleration routine

      implicit none

      integer maxor
      real*8 xx(*),a(*),b(*),r(*),sorthm(*),dum(*),rw(*)
      integer ncon(*),nop(*)
      real*8 xtemp(*),piv(*)
      integer irb(*),iirb(*),npvt(*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*),epn,epnc
      integer nrhs(1), neq, north, iter, idof, iptty
      real*8 tols, sum, su2, ad, beta          
      real*8 rho,rho_old,omega,alpha
      integer neqm1, neqp1, maxit, i, j1, j2, jj, kb
      integer nrhs_dum(1)
      parameter (tols=1.d-12)
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
c  residual(r)=rw
c  original residual=dum
c  p=sorthm
c  v=sorthm(+neq)    
c  x=xx
c  iter=n
c     
c     set some arrays equal
c     
      call equal_array(rw,r,neq,idof,nrhs_dum,nrhs)
c     
c     zero out some storage arrays
c     
      call constant_value(xx,0.0d00,neq,idof,nrhs_dum)
c     
      rho=1.0d00
c
c     forward and back sub for new solution
c
      call renumber_array(r,rw,iirb,neq,idof,nrhs,nrhs_dum)
      call sub_bksub1(neq,b,r,nop,npvt,piv)
      call renumber_array(rw,r,irb,neq,idof,nrhs_dum,nrhs)
c
c save copy of original residual(dum)
c
      call equal_array(dum,rw,neq,idof,nrhs_dum,nrhs_dum)    
      call equal_array(sorthm,rw,neq,idof,nrhs_dum,nrhs_dum)    
c
c loop for iterations
c
      do iter=1,maxit                                         
c
c set old norm
c
      rho_old=rho
      rho=0.0d00
      do i=1,neq
       rho=rho+dum(i)*rw(i)
      enddo
*     -------------------------------
*     compute p(sorthm(+neq)) in natural order
*     -------------------------------
c <<<<<<<<< skip for iter=1 >>>>>>>>>>>>>>>>>>>>>
      if(iter.gt.1) then
        beta=(rho*alpha)/(rho_old*omega)
        do i=1,neq
          sorthm(i)=rw(i)+beta
     $    *(sorthm(i)-omega*sorthm(i+neq))
        enddo
      endif
c <<<<<<<<< skip for iter=1 >>>>>>>>>>>>>>>>>>>>>
*     -------------------------------
*     compute v
*     -------------------------------
      do i=1,neq
        sum = 0.0d00     
        j1 = ncon(i) + 1
        j2 = ncon(i+1)
        do jj = j1,j2
          kb = ncon(jj)
          ad = a(jj - neqp1)
          sum = sum + ad * sorthm(kb)                        
        enddo
        sorthm(i+neq) = sum
      enddo         
      call renumber_array(r,sorthm(neq+1),iirb,neq,idof,nrhs,nrhs_dum)
      call sub_bksub1(neq,b,r,nop,npvt,piv)
      call renumber_array(sorthm(neq+1),r,irb,neq,idof,nrhs_dum,nrhs)
c
c calculate alpha               
c        
      sum=0.0d00
      do i=1,neq
        sum=sum+dum(i)*sorthm(i+neq)              
      enddo  
      alpha=rho/sum
*     ----------------------
*     update solution xx
*     ----------------------
      do i=1,neq
        xx(i)=xx(i)+alpha*sorthm(i)
      enddo
c
c update residuals
c
      do i=1,neq
        rw(i)=rw(i)-alpha*sorthm(i+neq)               
      enddo 
*     ----------------------
*     find sum squared of residual
*     ----------------------
      call residual(neq,rw,1,epnc,idof,nrhs_dum)
c
c check for convergence
c
      if(epnc.le.epn) then
        go to 5000
      endif
c      
*     ----------------------------------
*     compute t
*     ----------------------------------
      do i=1,neq
        sum = 0.0d00    
        j1 = ncon(i) + 1
        j2 = ncon(i+1)
        do jj = j1,j2
          kb = ncon(jj)
          ad = a(jj - neqp1)
          sum = sum + ad * rw(kb)                        
        enddo
        xtemp(i) = sum
      enddo         
      call renumber_array(r,xtemp,iirb,neq,idof,nrhs,nrhs_dum)
      call sub_bksub1(neq,b,r,nop,npvt,piv)
      call renumber_array(xtemp,r,irb,neq,idof,nrhs_dum,nrhs)
c
c calculate omega
c
      sum=0.0d00
      su2=0.0d00
      do i=1,neq
        sum=sum+xtemp(i)*xtemp(i)
        su2=su2+xtemp(i)*rw(i)
      enddo
      omega=su2/sum
*     ----------------------
*     update solution xx
*     ----------------------
      do i=1,neq
        xx(i)=xx(i)+omega*rw(i)
      enddo
c
c update residuals
c
      do i=1,neq
        rw(i)=rw(i)-omega*xtemp(i)
      enddo
*     ----------------------
*     find sum squared of residual
*     ----------------------
      call residual(neq,rw,1,epnc,idof,nrhs_dum)
c
c check for convergence
c
      if(epnc.le.epn) then
        go to 5000
      endif
      if(omega.lt.tols) then
	go to 5000
      endif

      enddo

      if(iptty .ne. 0 ) then
         write(iptty,34)
 34      format(/,1x,'Warning issued by 1 degree of freedom '
     &        ,'solver:')
         if (iter .gt. maxit) write(iptty,35) maxit
 35      format(2x,'Maximum number of iterations (',i4,') '
     &        ,'exceeded.')
         write(iptty,36) epnc,epn
 36      format(2x,'Final l2 norm = ',e14.6,', Tolerance = ',
     &        e14.6,/)
      endif
 5000 call equal_array(r,xx,neq,idof,nrhs_dum,nrhs)
      return
      end
