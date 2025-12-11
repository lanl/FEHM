      subroutine gensl_stress_coupled_3D
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
CD1 To solve stress eqs with full jacobian.
CD1 only one way coupling to fluid or heat 
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                     However, previous non-YMP
CD2                                     versions of FEHM exist, and
CD2                                     the current version differs
CD2                                     from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gensl2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:02:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Correct call to accum_derivative
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.25   Mon Jan 05 21:10:48 1998   gaz
CD2 change ieos description for ps=0 from 4 to -1
CD2 
CD2    Rev 1.24   Fri Nov 21 09:57:24 1997   gaz
CD2 changes to help SS problems
CD2 
CD2    Rev 1.22   Mon Apr 14 12:42:38 1997   gaz
CD2 improved rdof solution
CD2 
CD2    Rev 1.20   Fri Apr 26 15:23:02 1996   gaz
CD2 changes for mdnodes
CD2 
CD2    Rev 1.19   Thu Feb 15 10:48:26 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.18   Wed Feb 07 11:19:10 1996   gaz
CD2 added call mdnodes
CD2 
CD2    Rev 1.17   Mon Jan 29 16:24:56 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.16   Wed Jan 17 09:01:38 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.15   Tue Jan 16 14:31:10 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2
CD2    Rev 1.14   Fri Jan 12 17:49:00 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.13   11/15/95 10:35:34   gaz
CD2 changes so rdof works well
CD2 
CD2    Rev 1.12   07/10/95 15:49:32   llt
CD2 moved 999 continue at bottom up before mmrelblks
CD2 
CD2    Rev 1.11   05/02/95 16:18:36   pvcs
CD2 corrected error when using restart (gaz)
CD2 
CD2    Rev 1.10   04/25/95 09:40:46   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.9   03/30/95 09:11:42   gaz
CD2 gaz added more abs(bp
CD2 
CD2    Rev 1.8   03/29/95 15:29:58   gaz
CD2 gaz: changed more bp to abs(bp
CD2 
CD2    Rev 1.7   03/23/95 18:52:42   gaz
CD2 gaz put abs(bp instead of bp in one ifblock
CD2 
CD2    Rev 1.6   03/23/95 18:32:52   gaz
CD2 gaz changed call to rdof_new, new option for tolerance
CD2 
CD2    Rev 1.5   03/10/95 10:47:26   llt
CD2 modified to allow for Richard's Equation - gaz
CD2
CD2    Rev 1.4   01/28/95 13:54:58   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.3   05/11/94 16:03:08   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:38:20   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:41:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:24   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 None
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
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 n0, ipdgle, ipdeef, ipdepf, ipdmpf, ipdemf, ipdq, ipdeqh, ipdtpae,
CD4 ipdtpa, kgmres, lenreal, neq, bp, a, fdum, ieos, nelmdg, sx1, 
CD4 islord, nmat, nmatb, nrhs, nrhsb, iad, tmch, epe, g1, g2, g3,
CD4 iter, north, irdof, nb, nelm, nop, irb, iirb, nopt, npvt, piv,
CD4 nbnd, iback, iptty, maxor, h, c, ss, g, y, itert, mink, minkt,
CD4 overf, icoupl, iout
CD4 
CD4 Global Subprograms
CD4
CD4 Name      Type     Description
CD4 
CD4 mmgetblk  N/A      Allocates space for an array
CD4 mmrelblk  N/A      Deallocates space for an array
CD4 geneq2    N/A      Generate heat and mass equation for a node
CD4 geneq3    N/A      Generate heat equation for a node
CD4 dual      N/A      Computes dual porosity contributions to equations
CD4 switch    N/A      Reorders A matrix values
CD4 switchb   N/A      Reorders right hand side array values
CD4 normal    N/A      Normalizes matrix equations, computes RMS error
CD4 solve2    N/A      Solves two degree of freedom equation set
CD4 rdof_new  N/A      Solves two degree of freedom equation set using
CD4                        reduced degree of freedom method
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
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 icode        int         Returned from call to mmgetblk
CD5 neqp1        int         Number of equations plus 1
CD5 i            int         Do loop index over all nodes
CD5 id           int         Do loop index over all nodes
CD5 nsizea       int         Number of elements in A array
CD5 nsizea1      int         Total number of elements in 2neq x 2neq A
CD5                             array
CD5 nmat_save    int         Temporary storage of nmat(1) position
CD5 fdum2        real*8      RMS error of equations
CD5 facr         real*8      Factor used in setting convergence
CD5                             criterion
CD5 tolls        real*8      Tolerence factor
CD5 tollr        real*8      Tolerence factor
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
CD6 
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
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
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations
CD9 2.5.2 Solve nonlinear equation set at each time step
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN gensl2
CPS 
CPS mmgetblk - allocate space in storage arrays
CPS 
CPS FOR all nodes
CPS   Initialize right hand sides of equations to 0
CPS ENDFOR
CPS 
CPS Compute size of A matrix array
CPS 
CPS FOR all elements in A array
CPS   Initialize value to zero
CPS ENDFOR
CPS 
CPS FOR all nodes
CPS 
CPS   IF heat and mass solution is being performed
CPS     geneq2 - generate equation for this node
CPS   ELSE heat only solution is being performed
CPS     geneq3 - generate equation for this node
CPS   ENDIF
CPS 
CPS ENDFOR all nodes
CPS 
CPS IF this is an air only problem
CPS 
CPS   FOR each node
CPS     Determine factor to normalize with
CPS     FOR each node connected to this node
CPS       Renormalize the A matrix value
CPS     ENDFOR
CPS     Renormalize the residual value
CPS     Add to running total of RMS error
CPS   ENDFOR
CPS   Compute RMS error
CPS   IF this is the first iteration
CPS     Compute equation tolerance
CPS   ENDIF
CPS   
CPS   IF another iteration is needed
CPS     Compute equation tolerance needed for convergence
CPS     Compute the maxmimum number of inner iterations
CPS     mmgetblk - allocate space for gmres array
CPS     solve - obtain new solution to equation set
CPS     mmrelblk - free gmres memory
CPS     Increase number of iterations by one
CPS   ENDIF
CPS   
CPS ELSEIF this is a water only problem
CPS 
CPS   FOR each node
CPS     Determine factor to normalize with
CPS     FOR each node connected to this node
CPS       Renormalize the A matrix value
CPS     ENDFOR
CPS     Renormalize the residual value
CPS     Add to running total of RMS error
CPS   ENDFOR
CPS   Compute RMS error
CPS   IF this is the first iteration
CPS     Compute equation tolerance
CPS   ENDIF
CPS   
CPS   IF another iteration is needed
CPS     Compute equation tolerance needed for convergence
CPS     Compute the maxmimum number of inner iterations
CPS     mmgetblk - allocate space for gmres array
CPS     solve - obtain new solution to equation set
CPS     mmrelblk - free gmres memory
CPS     Increase number of iterations by one
CPS   ENDIF
CPS   
CPS ELSE
CPS   
CPS   dual - factor in dual porosity contribution
CPS 
CPS   IF reordering of the variables is specified
CPS     switch - reorder A matrix values
CPS     switchb - reorder b array values
CPS   ENDIF
CPS 
CPS   normal - normalize matrix equations, compute RMS error
CPS 
CPS   IF this is the first iteration
CPS     Compute convergence criterion
CPS   ENDIF
CPS 
CPS   IF another iteration is needed
CPS 
CPS     Compute equation tolerance needed for convergence
CPS     Compute the maxmimum number of inner iterations
CPS     mmgetblk - allocate space for gmres array
CPS   
CPS     IF a reduced degree of freedom model is not used
CPS       solve2 - solve linear equation set
CPS       Increase number of iterations by 1
CPS     ELSE a reduced degree of freedom model is to be used
CPS       rdof_new - solve linear equation set
CPS       Increase number of iterations by 1
CPS     ENDIF
CPS   
CPS     mmrelblk - release allocated memory of gmres array
CPS 
CPS     IF reordering of the variables is specified
CPS       switch - reorder A matrix values to original ordering
CPS       switchb - reorder b array values to original ordering
CPS     ENDIF
CPS   
CPS     dual - backsubstitute to obtain solution at matrix nodes
CPS     
CPS   ENDIF
CPS 
CPS ENDIF another iteration is needed
CPS 
CPS mmrelblk - release allocated memory for storage arrays
CPS 
CPS END gensl2
CPS 
C**********************************************************************
c
c solve isothermal air-water eqs with full jacobian(unsymmetric,2n by 2
c
      use davidi
      use comhi
      use comcouple
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsi
      use comwt
      use comfem
      implicit none

      real*8, allocatable :: sto5(:,:)
      integer icode,nsat
      integer ndex(2)
      integer neqp1
      integer i,j
      integer id,kb,inorm
      integer nsizea
      integer nsizea1
      integer nmat_save,dry_nodes,jm
      integer np, i1, i2, ij, ijj, ibp_max
      real*8 aij, apiv
      real*8 fdum2, fdumstress, machine
      real*8 facr
      real*8 tolls
      real*8 tollr
      real*8 res(2)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: a_save(:)
      real*8 a11, a22, adiag_tol, bp_max, bp_tot
	real*8 sum1,sum2,sum3,fdums
      parameter(adiag_tol=1.d-10, machine = 1.d-20)
      integer jj
      integer idbg

      integer neq3, indx(100)
	parameter (inorm = 1)
      neqp1=neq+1

c     zero out arrays
      do 10 i=1,neq
	  
        bp(i+nrhs(1))=0.0
        bp(i+nrhs(2))=0.0
        bp(i+nrhs(3))=0.0 
        bp(i+nrhs(4))=0.0 
        bp(i+nrhs(5))=0.0 
	 
 10   continue

      nsizea1=nelm(neqp1)-neqp1
      nsizea=nsizea1

      a=0.0

      call stressctr(4,0)
c
c      calculate volume strain
c
      call stressctr(6,0)

c s kelkar 11feb2011 ihms=15, 16, 17 are a special option
      if(ihms.eq.15.or.ihms.eq.17) then
c ps calculated below is a factor so that 
c pore volume = ps()*initial bulk volume (=sx1d)
c this factor ps() is used in therm w to calculate deni()=ps()*roho
c and the sx1d is multiplied in geneq1_.._..
         call porosity_wrt_displacements
      else
         call stressctr(-7,0)
      endif

c
c     get phase state and properties
c
       call varchk(0,0)
c
c fill NR and RHS in 5 by 5 format
c
      call stress_combine(9,idof_stress)
c
c generate stress equations
c  
      if(ifem.eq.1) then
         call geneq_stress_fem_3D()
        if(iPlastic.eq.1) then
          conv_pstrain = plastic_strain
        endif
      else
         do  id=1,neq
            if(icnl.eq.0) then
               call geneq_stress_uncoupled_3D(id)   
            else
               write(*,*) 'Stopping in gensl_stress_3D (icnl ne 0)'
               stop
            endif 
         enddo
      endif
      
c     
c     add derivatives wrt to p and t
c     
      call stressctr(15,0)
      
c     
c     apply boundary conditions
c     
      call stress_boun3(1,0)
c     
c     
c     get derivatives with respect to fluid and temperature variables
c     
      call stressctr(5,0)
c     
c     test for convergence of stress equations
c     
      
      fdum2=0.0d0
      
      if(tol_stress.ne.0.0d0) then
         bp_stress  = -1.0d0
         ibp_stress = 0
c     fdum1 = 1.0
         do id=1,neq
            bp_tot = abs(bp(id+nrhs(1)))+
     &           abs(bp(id+nrhs(2)))+abs(bp(id+nrhs(3)))
            if(bp_tot.gt.bp_stress) then
               ibp_stress = id
               bp_stress = bp_tot
            endif                  
         enddo
         
         fdum2 = bp_stress*bp_stress
         
         if(iad.eq.0) then
            bp_update = 1.0
            if(tol_stress.gt.0.0d0) then
               tol_stress1 = bp_stress*tol_stress  
            else
               tol_stress1 = abs(tol_stress)  
            endif
         endif
         if (bp_stress.le.tol_stress1.and.iad.ge.1) then
            fdumstress = -998.0
            fdums = bp_stress
            bpx = bp(ibp_stress+nrhs(1))
            bpy = bp(ibp_stress+nrhs(2))
            bpz = bp(ibp_stress+nrhs(3))
c     bp = 0.0d0	     
            go to 1001
            return
         else if(bp_update.lt.machine) then
            fdumstress = -997.0
            fdums = bp_stress
            bpx = bp(ibp_stress+nrhs(1))
            bpy = bp(ibp_stress+nrhs(2))
            bpz = bp(ibp_stress+nrhs(3))	    
c     bp = 0.0d0	     
            go to 1001
            return          
         endif
      endif
      
c     fill NR array and RHS  in 5 by 5 format
c     
 1001 call stress_combine(10,idof_stress)
      
c     generate fluid and energy balance equations equations
      do id=1, neq
         if(ihms.eq.15.or.ihms.eq.17) then  
      !! include derivatives of only pore volume wrt deformation
      !! in the accumulation term, but not those of permeability 
      !! in the flux terms
            call geneq1_stress_coupl_porosi(id)
c         elseif(ihms.eq.16) then
c      !! include derivatives of only  permeability 
c      !! in the flux terms but not of pore volume wrt deformation
c      !! in the accumulation term
c            call geneq1_stress_coupl_femperm(id)
         else
       !!  derives wrt p and T only
            call geneq1(id)
         endif
      enddo

c 
c get derivatives with respect to displacements (flux terms)
c 
c
      if(ihms.eq.16.or.ihms.eq.17) then
         if(ifem.eq.1) then
            call geneq_flowflux_perm()
         else
            do id=1,neq
               call geneq1_stress_coupl(id)
            enddo
         endif
      endif

      do id=1,neq
         if(ps(id).le.0.0) then
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0
         endif
      enddo

c
c     release memory          
c 
      call stress_perm(-2,0)      
c get derivatives with respect to displacements (porosity and cell volume)
c 
c	 call stressctr(7,0)
c
c
c  now reset order of equations in 5 by 5 format
c
      call stress_combine(12,idof_stress)
c
c for debugging..................................................
c      do idbg=1,12400
c         write(89,8989)a(idbg)
c      enddo
c 8989 format(g15.8)
c................................................................

      if(inorm.ne.0) then
c     block normalization
         allocate(dumn(100))
         ndex(1) = 1
         
         call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &        ,ndex,idof_stress,dumn(1),dumn(37),dumn(73),0,fdum2)
         
         deallocate(dumn)
         
         if(ndex(1).lt.0 ) then
            if (iout .ne. 0) then
               write(iout,*) '   '
               write(iout,*) '* singular matrix found during ',
     &              'normalization *'
               write(iout,*) '   '
               write(iout,'(i8,3g12.3)')
     &              abs(ndex(1)),(cord(abs(ndex(1)),i),i=1,3)
            end if
            if(iptty.gt.0)then
               write(iptty,*) '   '
               write(iptty,*) '* singular matrix found during ',
     &              'normalization *'
               write(iptty,*) '   '
               write(iptty,'(i8,3g12.3)') 
     &              abs(ndex(1)),(cord(abs(ndex(1)),i),i=1,3)
            endif
            iad=maxit
            return
         endif
      else
c     no normalization
         fdum2 = 0.0
         do i =1,neq
            do j = 1,5
               fdum2 = fdum2 + bp(i + nrhs(j))*bp(i + nrhs(j))
            enddo
	 enddo
      endif
c     
         fdum=sqrt(fdum2)

c 
          if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
	   else if(iad.ne.0.and.fdum.eq.0.0) then
	    go to 999
         endif    
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
             tollr = min(tol_stress1,0.1*fdum)
             tollr = max(tollr,abs_tol_stress)
             iter=maxsolve      
c allocate scratch storage dof*4*neq
        call storage_derivatives(1,1)
        allocate(dum(20*neq)) 

c    
            if(irdof.le.0) then 

               if(gdpm_flag.eq.0) then
                  if (igauss .gt. 1) then
                     call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop
     &                    ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,5,iptty,maxor,accm)
                  else
                     call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     &                    ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,5,iptty,maxor,accm)
                  end if
               else
                  if (igauss .gt. 1) then
                     call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &                    ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &                    ,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,5,iptty,maxor
     &                    ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                    ,mdof_sol)
                  else
                     call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat
     &                    ,nrhs,nelm,nelm_primary,nelm_primary,north
     &                    ,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,5,iptty,maxor
     &                    ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                    ,mdof_sol)
                  end if
               endif
               itert=itert+iter
               itotals=itotals+iter
            else
               allocate(sto5(neq,3))
               if (igauss .gt. 1) then
                  call  rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *                 ,tollr,irb,iirb,npvt,gmres,dum,piv
     *                 ,h,c,ss,g,y,iter,iback,3,iptty,maxor
     *                 ,overf,irdof,icoupl,0,sto5,accm)
               else
                  call  rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     *                 ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     *                 ,h,c,ss,g,y,iter,iback,3,iptty,maxor
     *                 ,overf,irdof,icoupl,0,sto5,accm)
               end if
               itert=itert+iter
               itotals=itotals+iter
               deallocate(sto5)
            endif
         deallocate(dum)   
         endif
           call storage_derivatives(0,1)
 998     continue
      
 999  continue
      

        
      return
      end


      
