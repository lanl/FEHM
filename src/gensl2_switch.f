      subroutine gensl2_switch
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
CD1 To solve isothermal air-water eqs with full jacobian with variable switch.
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
      use comwt
      implicit none

      real*8, allocatable :: sto5(:,:)
      integer icode,nsat
      integer ndex(2)
      integer neqp1
      integer i,j
      integer id,kb
      integer nsizea
      integer nsizea1
      integer nmat_save,dry_nodes,jm
      integer np, i1, i2, ij, ijj, ibp_max
      real*8 aij, apiv
      real*8 fdum2
      real*8 facr
      real*8 tolls
      real*8 tollr
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: a_save(:)
      real*8 a11, a22, adiag_tol,fdum2_tol, bp_max, bp_maxc, bp_tol
      real*8 tol_temp
c gaz debug 112019  
      real*8 dumdb1, dumdb2
      real*8, allocatable :: dumdb3(:)
      real*8, allocatable :: dumdb4(:)
c gaz debug 110319    fdum2_tol= 1.d-15  to  fdum2_tol= 1.d-30  
      parameter(adiag_tol=1.d-9, fdum2_tol= 1.d-30)
      integer jj
c      
      neqp1=neq+1
c     zero out arrays
      do 10 i=1,neq
         bp(i+nrhs(1))=0.0
         bp(i+nrhs(2))=0.0
 10   continue
      nsizea1=nelm(neqp1)-neqp1
      if(irdof.ne.13) then
         nsizea=4*nsizea1
         allocate(dum(8*neq))
      else
         nsizea=nsizea1
         allocate(dum(4*neq))
      endif
c
      a=0.0d0
c
      do  id=1,neq
c     
c     decide on equation type
c     
         if(ps(id).gt.0.0) then
               call geneq2_switch(id)
               call add_accumulation(id)          
         else if(ps(id).le.0.0) then        
               a(nelmdg(id)-neqp1+nmat(2))=sx1(id)   
c               a(nelmdg(id)-neqp1+nmat(3))=0.0d00    
               a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
c               a(nelmdg(id)-neqp1+nmat(4))=sx1(id)
               bp(id+nrhs(1))=0.0d00
               bp(id+nrhs(2))=0.0d00            
         endif 
      enddo

c gaz debug 112019
c      dumdb1 = 0.0
c      dumdb2 = 0.0
c      if(.not.allocated(dumdb3)) allocate(dumdb3(neq))
c      if(.not.allocated(dumdb4)) allocate(dumdb4(neq))
c      do i = 1, neq
c        dumdb1 = dumdb1 +abs(bp(i))
c        dumdb2 = dumdb2 +abs(bp(i+neq))
c        dumdb3(i) = bp(i)
c        dumdb4(i) = bp(i)
c      enddo
c
c arrange arrays for 1 dof solve (ie variable switching)
c
       call air_combine(1)
c     
c     coding for 1dof  only
c          
         fdum2=fdum2_tol
         ibp_max = 0
         bp_max = 0.0
         do id=1,neq
c gaz_test 080819             
            np=nelmdg(id)
            i1=nelm(id)+1
            i2=nelm(id+1)
            apiv=a(np-neqp1+nmat(1))+adiag_tol
c             if(apiv.lt.1.e-6) pause
            do ijj=i1,i2
               ij=nelm(ijj)
               aij=a(ijj-neqp1+nmat(1))/apiv
               a(ijj-neqp1)=aij
            enddo 

            bp_maxc = abs(bp(id+nrhs(1)))
	      if(bp_maxc.gt.bp_max) then
               bp_max = bp_maxc
               ibp_max = id
            endif
            bp(id+nrhs(1))=bp(id+nrhs(1))/apiv
            fdum2=fdum2+bp(id+nrhs(1))*bp(id+nrhs(1))            
c          endif  
         enddo
         bp_max = bp_max + fdum2_tol      
         fdum=sqrt(fdum2)
c check if fluxes are small enough to quit
         mink=neq_active
         bp_tol = 1.d-30
         if(g1.lt.0.) then
            tol_temp = tmch
            if(bp_max.le.tol_temp.and.iad.gt.iad_min-1) then
               fdum = -999.0
               fdum1 = bp_max
               go to 999
            else
c          bp_tol = bp_max/a(nelmdg(ibp_max)+nmat(1)-neqp1)
c               bp_tol = abs(g1)/bp_max*fdum
               f0=-1.0
c gaz 110519               
c               tmch = 1.e-20
               go to 996
            endif
         else
            bp_tol = 1.d-30
         endif

         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            if(ifree.eq.0) then
c               if(bp_max.gt.tmch) go to 995
               do i=1,neq
                  if(abs(bp(i+nrhs(1))).gt.tmch) go to 995
               enddo
            else
               if(abs(head_ck).gt.head_tol) go to 995
               do i=1,neq
                  if(izone_free_nodes(i).eq.0) then
                     if(abs(bp(i+nrhs(1))).gt.tmch) go to 995
                  endif
               enddo
            endif
            fdum=-1.0
c            minkt=minkt+mink           
c gaz 102919 make sure minimum iteration are done           
            if(iad.gt.iad_min-1) go to 999
 995        continue
            f0=-1.0
         endif
c gaz 101919          
c 996    if(fdum.gt.f0.or.iad.le.iad_min) then
996     if(fdum.gt.f0.or.iad.le.iad_min-1) then         
            facr=1.0
c            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
c            tollr=g3*max(tolls,tmch)
            tollr=g3*max(f0,tmch,bp_tol)
            iter=maxsolve      
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(1,1)
            if(gdpm_flag.eq.0) then
	       if(iback.eq.0) then
                  if (igauss .gt. 1) then
                     call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,
     &                    north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
                  else
                     call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,
     &                    north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
                  end if
	       else
                  if(.not.allocated(bsave)) then
c      form LU factors once			
                     allocate(bsave(nbd))
                     allocate(pivsave(n0))
                     if (igauss .gt. 1) then
                        call solve_new(neq,a,bsave,bp,nmat,nb,nrhs,nelm,
     &                       nop,north,tollr,irb,iirb,npvt,gmres,dum,
     &                       pivsave,h,c,ss,g,y,iter,0,1,iptty,maxor,
     &                       accm)
                     else
                        call solve_new(neq,a,bsave,bp,nmat,nmat,nrhs,
     &                       nelm,nelm,north,tollr,irb,iirb,npvt,gmres,
     &                       dum,pivsave,h,c,ss,g,y,iter,0,1,iptty,
     &                       maxor,accm)
                     endif
                  else
c      use already-formed LU factors 
                     if (igauss .gt. 1) then
                        call solve_new(neq,a,bsave,bp,nmat,nb,nrhs,nelm,
     &                       nop,north,tollr,irb,iirb,npvt,gmres,dum,
     &                       pivsave,h,c,ss,g,y,iter,1,1,iptty,maxor,
     &                       accm)
                     else
                        call solve_new(neq,a,bsave,bp,nmat,nmat,nrhs,
     &                       nelm,nelm,north,tollr,irb,iirb,npvt,gmres,
     &                       dum,pivsave,h,c,ss,g,y,iter,1,1,iptty,
     &                       maxor,accm)
                     endif
                  endif
               endif
            else
               if (igauss .gt. 1) then
                  call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &                 ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &                 ,npvt,gmres,dum,piv
     &                 ,h,c,ss,g,y,iter,iback,1,iptty,maxor
     &                 ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                 ,mdof_sol)
               else
                  call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat,nrhs
     &                 ,nelm,nelm_primary,nelm_primary,north,tollr,irb
     &                 ,iirb,npvt,gmres,dum,piv
     &                 ,h,c,ss,g,y,iter,iback,1,iptty,maxor
     &                 ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                 ,mdof_sol)
               end if
            endif
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(0,1)
            itert=itert+iter
            itotals=itotals+iter
            minkt=minkt+mink
c   sort out variable updates
        call air_combine(2)

       end if
      
 999  continue
      
      deallocate(dum)
      
      return
      end


      
