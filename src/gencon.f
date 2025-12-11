      subroutine  gencon(igrp,sia_iter,ncsp)
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
CD1 To solve for the concentrations for the current solute at the
CD1 current Newton iteration.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2              G. Zyvoloski   N/A     Initial implementation -
CD2                                     however, previous non-YMP
CD2                                     versions of FEHM exist, and the
CD2                                     current version may differ from
CD2                                     these
CD2 04-07-93      B. Robinson   N/A     Revised to implement the
CD2                                     chemical reaction models
CD2                                     between solutes
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gencon.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.12   Thu Feb 15 10:34:28 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.11   Wed Jan 17 08:55:56 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.10   Tue Jan 16 14:31:00 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.9   Fri Jan 12 17:48:14 1996   llt
CD2 changed mmgetblk agruments
CD2
CD2    Rev 1.8   12/13/95 10:28:20   robinson
CD2 Minor idof_sol fix
CD2 
CD2    Rev 1.7   08/18/95 10:35:42   llt
CD2 iter already defined, removed for cray
CD2
CD2    Rev 1.6   04/03/95 08:44:34   robinson
CD2 Changed to force at least one call to solver
CD2 
CD2    Rev 1.5   01/31/95 17:11:50   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.4   01/28/95 14:20:18   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.3   05/11/94 16:27:08   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:47:48   robinson
CD2  Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 16:15:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:04   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type     Use   Description
CD3 
CD3 tol_value    int      I/O   Flag indicating the status of the
CD3                                solution for the current solute
CD3 iter_counter int      I     Current iteration number
CD3 full_solution
CD3              int      I     Flag denoting whether full solution is
CD3 igrp         int      I     The current coupled grouping number
CD3                                to be performed
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4 
CD4 ipb, ipdeef, nbd, n0, ipdepf, ipdmpf, ipdmef, ipdq, ipdqh, ipdeqh,
CD4 ipdtpae, ipdtpae, ipdevef, kgmres, neq, bp, a, idualp, nelm, nelmdg,
CD4 npt, stored_residual, stored_derivative, fdum, igrp, iadd, tmch,
CD4 north, lenreal, irdof, b, nop, irb, iirb, nopt, npvt, piv, nbnd,
CD4 iptty, maxor, h, c, ss, g, y, mink, icoupl, overf, bn, iback
CD4 
CD4 Global Constants
CD4
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 solve           N/A      Solve linear equation set to determine
CD4                             new concentrations
CD4 rd1dof          N/A      Solve linear equation set to determine
CD4                             new concentrations using reduced degree
CD4                             of freedom solver
CD4 mmgetblk        N/A      Allocate memory for a array
CD4 mmrelblk        N/A      Release  memory for a array
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 neqp1        int         Number of equations plus 1
CD5 fdum2        real*8      RMS error at currrent iteration
CD5 id           int         Do loop index
CD5 i1           int         Do loop index
CD5 i2           int         Do loop index
CD5 idg          int         Parameter indicating position in nelmdg
CD5                             array
CD5 mi           int         Current place in concentration unknown
CD5                             array
CD5 pive         real*8      Normalization parameter for each equation
CD5 ikd          int         Do loop index
CD5 tolls        real*8      Tolerance parameter
CD5 iter         int         Maximum number of iterations of solver
CD5 gmres        real*8      Storage array for gmres solver
CD5 tollr        real*8      Tolerance parameter
CD5 sto1         real*8      Storage space used in solver
CD5 sto2         real*8      Storage space used in solver
CD5 sto3         real*8      Storage space used in solver
CD5 sto4         real*8      Storage space used in solver
CD5 sto5         real*8      Storage space used in solver
CD5 sto6         real*8      Storage space used in solver
CD5 sto7         real*8      Storage space used in solver
CD5 sto8         real*8      Storage space used in solver
CD5 sto9         real*8      Storage space used in solver
CD5 anz          real*8      Storage space used in solver
c changed to incorporate coupling 8/4/94
CD5 spec_num     int         Species number for variable
CD5                          we are computing rate contribution
CD5                          for (need better description)
CD5 deriv_spec_num int       Species number for variable we are
CD5                          taking the concentration with respect
CD5                          to (need better description)
CD5 igrp       int         Do loop index indicating the current
CD5                          group number
CD5 matnum       int         Index denoting the current submatrix 
CD5                          number
CD5 diag_mat     int         The number of the submatrix along the
CD5                          main diagonal
CD5 kgmres_sol   int         The storage required for gmres
CD5 idof_sol     int         The degrees of freedom for the 
CD5                          current grouping
CD5 
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following is a description of the functions are performed by
CD6 this routine:
CD6 
CD6   Sets residuals at each node to 0 (in do-loop).
CD6   
CD6   Sets all A-matrix values to 0 (in do-loop) after computing the
CD6   total number.
CD6   
CD6   Loops over each node, calling coneq1 to compute the advection
CD6   and dispersion terms for that node.
CD6   
CD6   Calls dual to compute dual porosity terms if a dual porosity
CD6   simulation is being performed.
CD6   
CD6   Loops through each node, computing the values of the residuals
CD6   and the elements of the A-matrix.
CD6   
CD6   Computes the RMS error.
CD6   
CD6   If solutions had converged at the previous iteration, check to
CD6   see if they are still converged.  If they are, set flag to
CD6   indicate convergence, but if they are no longer converged, reset
CD6   flag to perform another iteration.  If solution had not
CD6   converged at previous iteration, prepare for another solution at
CD6   the current iteration.
CD6   
CD6   If a solution is needed, set parameter values, and call either
CD6   solve or rd1dof (for reduced DOF solver) for solution at current
CD6   iteration.
CD6   
CD6   If this is a dual porosity solution, call dual to obtain the
CD6   solution at matrix nodes by back substitution.
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.3.4 Solute-transport equations
CD9 2.5.2 Solve nonlinear equation set at each time step
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN gencon
CPS 
CPS   Set sum of squares residual to 0
CPS   Set sub matrix counter to 0
c changed to incorporate coupling: take coneq1 and dual out of 
c gencon and put it directly into csolve 
CPS   FOR each unknown
CPS     Compute integer indexes
CPS     FOR each coupled species
CPS       Compute submatrix number along the diagonal
CPS       Compute normalization factor
CPS       Compute normalized residual
CPS       For each derivative species
CPS         FOR each unknown connected to this unknown
CPS           Compute normalized A matrix term
CPS         ENDFOR each connected unknown
CPS       END FOR each derivative species
CPS     ENDFOR each coupled species 
CPS   Compute updated sum of squares residual
CPS   ENDFOR each unknown
CPS 
CPS Compute RMS residual
CPS 
CPS IF a convergence check is to be performed using old convergence...
CPS ... criterion
CPS   IF convergence is no longer met
CPS     Set flag to indicate that concentrations must be recalculated
CPS   ELSE convergence is still met
CPS     Set flag to indicate there is still convergence
CPS   ENDIF
CPS ELSE we set the new convergence criterion
CPS   
CPS   IF this is the first outer iteration
CPS     Compute convergence criterion
CPS     IF this is a full solution
CPS       Set flag to indicate that concentrations need to be calculated
CPS     ELSE
CPS       Set flag to indicate only one solution is needed
CPS     ENDIF
CPS   ELSE
CPS     Set convergence criterion to previously determined value
CPS   ENDIF
CPS     
CPS   IF convergence has been achieved
CPS     Set flag to indicate such
CPS   ENDIF
CPS   
CPS ENDIF
CPS 
CPS IF a full solution is called for
CPS 
CPS   IF concentrations are being recalculated
CPS 
CPS     Compute new tolerance parameters
CPS     IF the residual is out of bounds
CPS       Set the convergence criterion to a appropriate value
CPS     ENDIF
CPS     
CPS     mmgetblk - allocate space for gmres storage
CPS     
CPS     IF reduced degree of freedom solver is not to be used
CPS       solve - solve for new concentrations
CPS     ELSE reduced degree of freedom solver is to be used
CPS       rd1dof - solve for new concentrations
CPS     ENDIF
CPS     
CPS     mmrelblk - free space for gmres storage
CPS   
CPS   IF this is a dual porosity simulation
CPS     dual - compute concentrations at matrix nodes
CPS   ENDIF
CPS   
CPS ENDIF
CPS   
CPS END gencon
CPS 
C***********************************************************************

      use comcouple
      use comrxnb
      use comrxni
      use davidi
      use comji
      use comhi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai

      implicit none

      real*8, allocatable :: anz(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: sto5(:)
      real*8, allocatable :: sto6(:)
      real*8, allocatable :: sto7(:)
      real*8, allocatable :: sto8(:)
      real*8, allocatable :: sto9(:)
      real*8, allocatable :: sto10(:)
      real*8, allocatable :: diag(:,:)
      real*8, allocatable :: diagi(:,:)
      real*8, allocatable :: rdum(:)
      integer iflg,itype,mcount, ioption

      integer iter_counter,full_solution,tol_value,neqp1,id,i1,i2,idg
      integer ikd,icode,spec_num,deriv_spec_num,igrp,matnum,diag_mat
      integer kgmres_sol,idof_sol,sia_iter,ncsp,i
      real*8 fdum2,pive,tolls,tollr
c gaz 070922 added toldil 
      real*8 toldil
      parameter(toldil = 1.d-20)
      neqp1=neq+1
      idof_sol = ncsp

c     zero out arrays
c     normalize equations
      allocate(diag(idof_sol,idof_sol))
      allocate(diagi(idof_sol,idof_sol))
      allocate(rdum(idof_sol))

c     Try putting in uncoupled normalization
      ioption = 0
      call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg,
     &     nmatb,idof_sol,diag,diagi,rdum,ioption,fdum2)
      deallocate(diag)
      deallocate(diagi)
      deallocate(rdum)
      if(nmatb(1).lt.0 ) then
         if (iout .ne. 0) then
            write(iout,*) '   '
            write(iout, 100) 
            write(iout,*) '   '
            write(iout,'(i8,3g12.3)')
     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         end if
         if (iptty.gt.0) then
            write(iptty,*) '   '
            write(iptty, 100) 
            write(iptty,*) '   '
            write(iptty,'(i8,3g12.3)') 
     &           abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         endif
         iad=maxit
         return
      endif
 100  format ('* singular matrix found during normalization *')
      fdum=sqrt(fdum2)
c gaz 070922 added residual to call done     
c      if(fdum.le.toldil) go to 101
      if(sia_iter.eq.0)fzero(igrp)=fdum
      tolls=max(fzero(igrp),min(g1*fdum,g2*fdum2))
      tollr=g3*tolls      
c gaz 070922 changed to maxsolve      
c      iter = 3 * north
      iter = maxsolve
      call storage_derivatives(1,1)
      if( irdof .ne. 9 ) then
               
         if (idof_sol.lt.6) then
            allocate(dum(neq*idof_sol*4))
         else if (idof_sol.eq.6) then
            allocate(dum(neq*idof_sol*5))
         else
            allocate(dum(neq*idof_sol*idof_sol))
         endif
         if(gdpm_flag.eq.0) then            
            if (igauss .gt. 1) then
               call solve_new(neq,a,b,bp,nmat_sol,nb_sol,nrhs_sol,
     *              nelm,nop,north,tollr,irb,iirb,npvt,gmres,
     *              dum,piv,h,c,ss,g,y,iter,iback
     *              ,idof_sol,iptty,maxor,accm)
            else
               call solve_new(neq,a,b,bp,nmat_sol,nmat_sol,nrhs_sol,
     *              nelm,nelm,north,tollr,irb,iirb,npvt,gmres,
     *              dum,piv,h,c,ss,g,y,iter,iback
     *              ,idof_sol,iptty,maxor,accm)
            end if
         else
            if (igauss .gt. 1) then
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &              ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &              ,npvt,gmres,dum,piv
     &              ,h,c,ss,g,y,iter,iback,idof_sol,iptty,maxor
     &              ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm,
     &              mdof_sol)
            else
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat,nrhs
     &              ,nelm,nelm_primary,nelm_primary,north,tollr,irb
     &              ,iirb,npvt,gmres,dum,piv
     &              ,h,c,ss,g,y,iter,iback,idof_sol,iptty,maxor
     &              ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm,
     &              mdof_sol)
            end if
         endif
         deallocate(dum)
      else
         write(ierr,*) 'Gencon not yet set for rd1dof.'
         write(ierr,*) 'Stop in gencon'
         stop
         if (igauss .gt. 1) then
            call rdof_new(neq,a,b,bp,nmat_sol,nb_sol,nrhs_sol,
     *           nelm,nop,north,tollr,irb,iirb,npvt,gmres,
     *           dum,piv,h,c,ss,g,y,iter,iback,itype,iptty,maxor
     *           ,overf,irdof,icoupl,mcount,sto5,accm)  
         else   
            call rdof_new(neq,a,b,bp,nmat_sol,nmat_sol,nrhs_sol,
     *           nelm,nelm,north,tollr,irb,iirb,npvt,gmres,
     *           dum,piv,h,c,ss,g,y,iter,iback,itype,iptty,maxor
     *           ,overf,irdof,icoupl,mcount,sto5,accm)
         end if
      endif
      call storage_derivatives(0,1)
c gaz 070922 jump to 101 if solution attained before solver (commented out 080622)
c101   continue
c      do id = 1, idof_sol
c       do i = 1, neq    
c        if(abs(bp(i+nmat_sol(id))).le.toldil) bp(i+nmat_sol(id))= 0.0d0
c       enddo
c      enddo
      if ( idualp .ne. 0) then
         call dual(13)
      endif

      return
      end
