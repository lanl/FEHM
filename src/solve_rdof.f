      subroutine solve_rdof(neq,a,b,r,na,nb,nrhs,ncon,nop,north,epn
     * ,irb,iirb,npvt,stor1,dum,piv
     * ,h,c,s,g,y,iter,iback,idof,iptty,maxor,icoupl,tollr,overf,accm)
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
CD1 To perform pre-conditioned conjugate gradient solution of a set of
CD1 linear, algebraic equations with reduced degree of freedom.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 3-2-94       G. Zyvoloski   97      Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/solve_rdof.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:16:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:46 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Thu Sep 12 08:25:34 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.8   Fri Apr 26 16:01:58 1996   gaz
CD2 slightly stricter stopping criteria
CD2 
CD2    Rev 1.7   Fri Feb 02 12:03:28 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   Tue Jan 16 14:32:18 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.5   Tue Jan 09 14:12:24 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.4   05/12/94 10:04:42   llt
CD2 corrected pvcs log info
CD2 
CD2    Rev 1.3   05/11/94 16:18:46   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:43:58   robinson
CD2 Added prologs
CD2
CD2    Rev 1.1   03/18/94 16:05:22   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:32   pvcs
CD2 original version
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 neq          integer  I      Number of nodes in the matrix (total
CD3                                  number of unknowns divided by the
CD3                                  number of degrees of freedom)
CD3 a            real*8   I      Array containing solution matrix
CD3 b            real*8   I      Array used for incomplete factorization
CD3 r            real*8  I/O     On entry, right hand side of matrix,
CD3                                  on exit, the solution vector
CD3 na           integer  I      Pointer integer array for the matrix
CD3                                  a in a multiple degree of
CD3                                  freedom solution
CD3 nb           integer  I      Pointer integer array for the matrix
CD3                                  b in a multiple degree of
CD3                                  freedom solution
CD3 nrhs         integer  I      Pointer integer array for the array
CD3                                  r in a multiple degree of
CD3                                  freedom solution
CD3 ncon         integer  I      Node connectivity array for the
CD3                                  matrix a
CD3 nop          integer  I      Connectivity matrix for factorization
CD3                                  matrix
CD3 north        integer  I      Number of orthogonalizations for the
CD3                                  GMRES solution
CD3 epn          real*8   I      Tolerance for equations
CD3 irb          integer  I      Renumber array used when renumbering is
CD3                                  chosen
CD3 iirb         integer  I      Renumber array used when renumbering is
CD3                                  chosen
CD3 npvt         integer  I      Array of pivot positions in nop matrix
CD3 stor1        real*8   I      Array for GMRES storage
CD3 stor2        real*8   I      Storage array
CD3 stor3        real*8   I      Storage array
CD3 stor4        real*8   I      Storage array
CD3 stor5        real*8   I      Storage array
CD3 piv          integer  I      Array of pivot positions
CD3 h            real*8   I      Array used in factorization
CD3 c            real*8   I      Array used in factorization
CD3 s            real*8   I      Array used in factorization
CD3 g            real*8   I      Array used in factorization
CD3 y            real*8   I      Array used in factorization
CD3 iter        integer   O      Number of iterations needed for the
CD3                                  solution
CD3 iback       integer   I      Flag controlling whether
CD3                                  factorization is to be performed
CD3 idof        integer   I      Number of degrees of freedom in the
CD3                                  equation set
CD3 iptty       integer   I      Flag to set the destination of
CD3                                  warning message output
CD3 maxor       integer   I      Maximum number of orthogonalizations
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 iptty                  O    Destination of warning messages
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
CD4 Identifier     Type     Description
CD4 
CD4 residual       N/A      Computes the sum of the squared residuals
CD4                            of an array
CD4 constant_value N/A      Sets all value in an array to a given value
CD4 sub_ilu1       N/A      Performs ilu factorization - one degree-of-
CD4                            freedom
CD4 sub_ilu2       N/A      Performs ilu factorization - two degrees-of-
CD4                            freedom
CD4 sub_ilu3       N/A      Performs ilu factorization - three degrees-
CD4                            of-freedom
CD4 sub_ilu4       N/A      Performs ilu factorization - four degrees-
CD4                            of-freedom
CD4 sub_gmres1     N/A      Performs GMRES accelerated solution of
CD4                            equations - one degree-of-freedom
CD4 sub_gmres2     N/A      Performs GMRES accelerated solution of
CD4                            equations - two degrees-of-freedom
CD4 sub_gmres3     N/A      Performs GMRES accelerated solution of
CD4                            equations - three degrees-of-freedom
CD4 sub_gmres4     N/A      Performs GMRES accelerated solution of
CD4                            equations - four degrees-of-freedom
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
CD5
CD5 Local variables
CD5
CD5 NONE
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 It is assumed that the developer using this reuse component
CD6 allocates appropriate storage space for all arrays in the calling
CD6 program - no checks are provided.
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7
CD7 This is the main calling routine for the solver, so the
CD7 requirements listed in the next section are only the upper level
CD7 requirements.  Detailed requirements are listed in the routines
CD7 called from solve_new.
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8
CD8 3.2         Solve Linear Equation Set
CD8 3.3         Provide Multiple Degree-of-Freedom Option
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************

      implicit none

      integer idof, neq, maxor
      real*8 a(*),b(*),r(*)
      real*8 stor1(*),piv(neq,*)
      real*8 h(maxor,*),c(*),s(*),g(*),y(*), epn, anorm
      real*8 dum(*)
      integer ncon(*),nop(*)
      integer irb(*),iirb(*),npvt(*)
      integer na(*),nb(*),nrhs(*)
      integer north, iter, iback, iptty
      integer isolve
      integer ibcgs
      integer icoupl
      real*8 tollr,overf
      integer na_save(2)
      character*4 accm

c
c determine acceleration method
c
      if(accm.eq.'bcgs') then
       ibcgs = 1
      elseif (accm.eq.'gmre') then
       ibcgs = 0
      elseif (accm.eq.'    ') then
c leaves the old input format valid
      endif
      if( idof .gt. 6 ) then
         isolve = 1
      else
         isolve = 0
      end if
      if( idof .gt. 6 ) then
         isolve = 1
      else
         isolve = 0
      end if
c
c     check initial residual
c
      call residual(neq,r,1,anorm,abs(idof),nrhs)
c
c     if initial residual less then tolerance(epn) then
c     quit and set correction=0.0
c
      if(anorm.le.epn) then
c GAZ feb 9 96
c         call constant_value(r,0.0d00,neq,idof,nrhs)
c         goto 9000
         epn=anorm*0.1d00    
      endif
c
c     call ilu factorization at preconditioner stage
c   
      if(iback.eq.0) then
         if(isolve .eq. 1 ) then
            call sub_ilun(neq,idof,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,dum,piv)  
         elseif(idof.eq.1) then 
            call sub_ilu1(neq,a,b,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),piv)
         else if(idof.eq.2) then
            call sub_ilu2(neq,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),dum(2*neq+1),piv)
         else if(idof.eq.3) then
            call sub_ilu3(neq,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),dum(3*neq+1),dum(6*neq+1)
     *           ,piv)
         else if(idof.eq.4) then
            call sub_ilu4(neq,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),dum(4*neq+1),dum(8*neq+1)
     *           ,piv)  
         else if(idof.eq.5) then
            call sub_ilu5(neq,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),dum(5*neq+1),dum(10*neq+1)
     *           ,dum(15*neq+1),piv)  
         else if(idof.eq.6) then
               call sub_ilu6(neq,a,b,na,nb,ncon,nop
     *              ,irb,iirb,npvt,stor1,dum(1),dum(6*neq+1)
     *              ,dum(12*neq+1),dum(18*neq+1),dum(24*neq+1),piv)  
         elseif(idof.eq.-3) then 
            call sub_ilu1(neq,a,b,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),piv)
         else if(idof.eq.-6) then
          na_save(1)=na(3)
          na_save(2)=na(4)
          na(3)=na(7)
          na(4)=na(8)
            call sub_ilu2(neq,a,b,na,nb,ncon,nop
     *           ,irb,iirb,npvt,stor1,dum(1),dum(2*neq+1),piv)
          na(3)=na_save(1)
          na(4)=na_save(2)
         endif
      endif
c
c     call acceleration routine
c
      if(isolve .eq. 1 ) then
         if (ibcgs.eq.0) then
         call sub_gmresn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(idof*neq+1),dum(2*idof*neq+1)
     *        ,dum(3*idof*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty)
         else
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *           ,irb,iirb,npvt,dum(1),dum(idof*neq+1),dum(2*idof*neq+1)
     *           ,dum(3*idof*neq+1),piv
     *          ,h,c,s,g,y,iter,idof,iptty,maxor)
         end if
      elseif(idof.eq.1) then
         if (ibcgs.eq.0) then
            call sub_gmres1(neq,a,b,r,ncon,nop,north,stor1,epn
     *           ,irb,iirb,npvt,dum(1),dum(neq+1),dum(2*neq+1)
     *           ,dum(3*neq+1),piv
     *           ,h,c,s,g,y,iter,idof,iptty,maxor)
         else
	  if (north.eq.1) then
c           call sub_bcgstabn1(neq,a,b,r,ncon,nop,north,stor1,epn
c    *           ,irb,iirb,npvt,dum(1),dum(neq+1),dum(2*neq+1)
c    *           ,dum(3*neq+1),piv
c    *           ,h,c,s,g,y,iter,idof,iptty,maxor)
            call sub_bcgstab1(neq,a,b,r,ncon,nop,north,stor1,epn
     *           ,irb,iirb,npvt,dum(1),dum(neq+1),dum(2*neq+1)
     *           ,dum(3*neq+1),piv
     *           ,h,c,s,g,y,iter,idof,iptty,maxor)
          else
             call sub_bcgstabn1(neq,a,b,r,ncon,nop,north,stor1,epn
     *            ,irb,iirb,npvt,dum(1),dum(neq+1),dum(2*neq+1)
     *            ,dum(3*neq+1),piv
     *            ,h,c,s,g,y,iter,idof,iptty,maxor)
!       call sub_bcgstabn1(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
!     *           ,irb,iirb,npvt,dum(1),dum(2*neq+1),dum(4*neq+1)
!     *           ,dum(6*neq+1),piv
!     *          ,h,c,s,g,y,iter,idof,iptty,maxor)
          endif
         endif
      else if(idof.eq.2) then
         if(ibcgs.eq.0) then
           call sub_gmres2(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *           ,irb,iirb,npvt,dum(1),dum(2*neq+1),dum(4*neq+1)
     *           ,dum(6*neq+1),piv
     *           ,h,c,s,g,y,iter,idof,iptty,maxor)
        else
c          call sub_bcgs2(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *           ,irb,iirb,npvt,dum(1),dum(2*neq+1),dum(4*neq+1)
     *           ,dum(6*neq+1),piv
     *          ,h,c,s,g,y,iter,idof,iptty,maxor)
        endif
      else if(idof.eq.3) then
         if(ibcgs.eq.0) then
         call sub_gmres3(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(3*neq+1),dum(6*neq+1)
     *        ,dum(9*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty,maxor)
      else
c        call sub_bcgs3(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(3*neq+1),dum(6*neq+1)
     *        ,dum(9*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty,maxor)
      endif
      else if(idof.eq.4) then
         if(ibcgs.eq.0) then
         call sub_gmres4(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(4*neq+1),dum(8*neq+1)
     *        ,dum(12*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty)
         else
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(4*neq+1),dum(8*neq+1)
     *        ,dum(12*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty,maxor)
         endif
      else if(idof.eq.5) then
         if(ibcgs.eq.0) then
         call sub_gmres5(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(5*neq+1),dum(10*neq+1)
     *        ,dum(15*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty)
         else
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(5*neq+1),dum(10*neq+1)
     *        ,dum(15*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty,maxor)
         endif
      else if(idof.eq.6) then
         if(ibcgs.eq.0) then
         call sub_gmres6(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(6*neq+1),dum(12*neq+1)
     *        ,dum(18*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty)
         else
       call sub_bcgstabnn(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(6*neq+1),dum(12*neq+1)
     *        ,dum(18*neq+1),piv
     *        ,h,c,s,g,y,iter,idof,iptty,maxor)
         endif
      else if(idof.eq.-3) then
         call sub_gmres31(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(3*neq+1),dum(6*neq+1)
     *        ,dum(9*neq+1),piv
     *        ,h,c,s,g,y,iter,3   ,iptty,maxor,icoupl,tollr,overf)
      else if(idof.eq.-6) then
         call sub_gmres62(neq,a,b,r,na,nb,nrhs,ncon,nop,north,stor1,epn
     *        ,irb,iirb,npvt,dum(1),dum(6*neq+1),dum(12*neq+1)
     *        ,dum(18*neq+1),piv
     *        ,h,c,s,g,y,iter,6   ,iptty,icoupl,tollr,overf)
      endif
 9000 continue
      return
      end
