      subroutine gensl2
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
CD1 To solve isothermal air-water eqs with full jacobian.
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
CD2
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
      use comwt
      use comriv
      use comwellphys
c gaz 040724
      use com_prop_data, only : ihenryiso, ctest, xnl_ngas, xnl_max 
      use com_nondarcy
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
      integer iwellp_chk, iparchek
      real*8 aij, apiv
      real*8 fdum2
      real*8 facr
      real*8 tolls
      real*8 tollr
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: a_save(:)
      real*8 a11, a22, adiag_tol, bp_max, bp_maxc, bp_tol
c gaz debug 112019  
      real*8 dumdb1, dumdb2
      real*8, allocatable :: dumdb3(:)
      real*8, allocatable :: dumdb4(:)   
      integer, allocatable :: ieosdum(:)  
c gaz 112419      
      integer i_uzsz_re, num_unit
      logical op,ex
      parameter(i_uzsz_re=1)
      parameter(adiag_tol=1.d-10)
      parameter(iparchek = 0)
      integer jj
      neqp1=neq+1
c test of coefficients
      jj = phi(1)+s(1)
       if(iparchek.eq.1) then
       write(ierr,*) 'fe coef. upper diag'
       write(ierr,*)
       do i = 1,neq
        i1 = nelmdg(i)+1
        i2 = nelm(i+1)
        do j = i1,i2
         kb = nelm(j)
         iw = istrw(j-neqp1)
         a11 = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
         write(ierr,355)i,kb,iw,a11
        enddo
       enddo   
355   format(i6,1x,i6,1x,i6,1x,f12.3)        
       stop
       endif
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
      fdum2=0.
      bp_max = 0.0
      if(ifree.ne.0) then
         dry_zone=1
         dry_nodes=0;nsat=0

         do id=1,neq
            i1=nelm(id)+1
            i2=nelm(id+1)
            if(rlxyf(id).lt.dry_tol) then
               dry_zone(id)=0
               dry_nodes=dry_nodes+1
            endif
            do j=i1,i2
               kb=nelm(j)
               if(rlxyf(kb).gt.dry_tol) dry_zone(id)=1
            end do
            if(dry_zone(id).eq.0) nsat=nsat+1
         end do
c		write(*,*) dry_nodes, ' nodes are dry'
c	if(ifree.ne.0) write(*,*) nsat,' dry nodes are frozen'
      endif
c GAZ 11/01/08 get drift flux info (must be call after rel perms)
      call wellphysicsctr(1,0) 
c gaz 120323 added coding for boundary flow from gas diffusion
         if(idiff_iso.ne.0)  then
          q_gas = 0.0d0
         endif
c gaz  111824 debug
c gaz  120324 debug
c      nd_test= -1   
c      if(l.eq.1.and.nd_test.ne.0) then   
c       if(fdum.eq.-1.0) then 
c        nd_test = 100
c        iad_min = iad_min_sv 
c        iad = iad_min_sv
c        nd_test = iad_min_sv
c        return
c       endif
c      else
c       nd_test = -1
c       iad_min = iad_min_sv
c      endif
c
      if(nd_test_write.gt.0) then
c gaz debug 151224          
       inquire(file='nd_comp1-21.chk', opened = op)
       inquire  (file='nd_comp1-21.chk', exist = ex)
       if(ex.and.op.and.l.eq.1.and.iad.eq.0) then
        inquire(file='nd_comp1-21.chk', number = num_unit)
        close(num_unit)
        open(unit=55,file='nd_comp1-21.chk',status='unknown')
         write(55,180) 
       else if(.not.op.and.l.eq.1.and.iad.eq.0) then
        open(unit=55,file='nd_comp1-21.chk',status='unknown')
         write(55,180) 
       endif
c
       inquire(file='nd_comp360-380.chk', opened = op)
       inquire  (file='nd_comp360-380.chk', exist = ex)
       if(ex.and.op.and.l.eq.1.and.iad.eq.0) then
        inquire(file='nd_comp360-380.chk', number = num_unit)
        close(num_unit)
        open(unit=56,file='nd_comp360-380.chk',status='unknown')
         write(56,180) 

       else if(.not.op.and.l.eq.1.and.iad.eq.0) then
        open(unit=56,file='nd_comp360-380.chk',status='unknown')
         write(56,180)
       endif
c
       inquire(file='nd_comp365-385.chk', opened = op)
       inquire  (file='nd_comp365-385.chk', exist = ex)
       if(ex.and.op.and.l.eq.1.and.iad.eq.0) then
        inquire(file='nd_comp365-385', number = num_unit)
        close(num_unit)
        open(unit=57,file='nd_comp365-385.chk',status='unknown')
         write(57,180) 
       else if(.not.op.and.l.eq.1.and.iad.eq.0) then
        open(unit=57,file='nd_comp365-385.chk',status='unknown')
         write(57,180) 
       endif
180   format(t5,'ts',t14,'i',t21,'kb',t34,'vxy',t51,'dvapi',t67,'dvapkb'
     &    ,t85,'dvaei',t97,'dvaekb')
      endif
c end nd_test_write
      do 101 id=1,neq
c     
c     decide on equation type
c     
         if (nwellphy .ne. 0) then
            iwellp_chk = iwellp(id)
         else
            iwellp_chk = 0
         end if

         if(ps(id).gt.0.0) then
            if(iwt_uz.ne.0) then
               call geneq2_uz_wt(id)
               call add_accumulation(id)
            else if(iwellp_chk.ne.0) then
c no accumulation term in the wellbore            
               call geneq2_wellphysics(id)	
            else if(ianpe.ne.0.and.irdof.eq.13) then
               call geneq2_ani(id)
               call add_accumulation(id)   
            else if(ianpe.ne.0.and.irdof.ne.13) then
               call geneq2_ani_2p(id)
               call add_accumulation(id)                            
            else if(irich.ne.0) then
               call geneq2_rich(id)
               call add_accumulation(id)	
            else if(ifree.ne.0) then
              if(iriver.eq.2.and.id.gt.neq_primary) then
               call geneq2_wtsi_well(id)
              else
               call geneq2_wtsi(id)
              endif
               call add_accumulation(id)
c set updates to zero if s=0 and all neighbors are s=0
               if(dry_zone(id).eq.0.and.sk(id).eq.0.0) then	
c do we need to make sure this is a wtsi node?						 
                  a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
                  bp(id+nrhs(1))=0.0d00
               endif
            else if(ihenryiso.ne.0) then
               call geneq2_sol(id)
               call add_accumulation(id)
            else if(nd_flow) then
               call geneq2_nondarcy(id)
               call add_accumulation(id)   
            else
               call geneq2(id)
               call add_accumulation(id)
            endif
         else if(ps(id).le.0.0) then
            if(irdof.ne.13) then
               a(nelmdg(id)-neqp1+nmat(2))=0.0d00     
               a(nelmdg(id)-neqp1+nmat(3))=0.0d00    
               a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
               a(nelmdg(id)-neqp1+nmat(4))=sx1(id)
               bp(id+nrhs(1))=0.0d00
               bp(id+nrhs(2))=0.0d00
            else
               a(nelmdg(id)-neqp1)=sx1(id)
               bp(id)=0.0d00
            endif
         endif
 101  continue

           
c     
c
c add correction for saturations over 1.0
c
      if(iflux_ts.ne.0) call cascade_sat(2)
c
      if(irdof.eq.14) then
c
c     first call air_rdof to rearrange equations
c     
         call air_rdof(0,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
c     
         allocate(a_save(nsizea1))

         fdum2=0.0
         do i=1,neq
            a11=max(a(nelmdg(i)-neqp1+nmat(1)),a(nelmdg(i)
     &           -neqp1+nmat(2)))
            a22=max(a(nelmdg(i)-neqp1+nmat(3)),a(nelmdg(i)
     &           -neqp1+nmat(4)))
            i1=nelm(i)+1
            i2=nelm(i+1)
            if(a11.ne.0.0) then
               do jj=i1,i2
                  a(jj-neqp1+nmat(1))=a(jj-neqp1+nmat(1))/a11
                  a(jj-neqp1+nmat(2))=a(jj-neqp1+nmat(2))/a11
               enddo
               fdum2=fdum2+bp(i+nrhs(1))**2
               bp(i+nrhs(1))=bp(i+nrhs(1))/a11
            endif
            if(a22.ne.0.0) then
               do jj=i1,i2
                  a(jj-neqp1+nmat(3))=a(jj-neqp1+nmat(3))/a22
                  a(jj-neqp1+nmat(4))=a(jj-neqp1+nmat(4))/a22
               enddo
               bp(i+nrhs(2))=bp(i+nrhs(2))/a22
               fdum2=fdum2+bp(i+nrhs(2))**2
            endif
         enddo
         fdum=sqrt(fdum2)
         if(fdum.eq.0.0) then
            deallocate(a_save)
            go to 999
         endif
c     
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.tmch) go to 991
               if(abs(bp(i+nrhs(2))).gt.tmch) go to 991
            enddo
            fdum=-1.0
c gaz 042219 mink not calculated correctly (hopefully this works            
c            minkt=minkt+mink
            go to 999
 991        continue
            f0=-1.0
         endif
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
c            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tolls=facr*f0
            tollr=tolls*g3
            if(tmch.gt.tollr) tollr=tmch
            iter=maxsolve      
            call storage_derivatives(1,1)
            allocate(sto5(neq,2))
            if (igauss .gt. 1) then
               call  combine_a(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *              ,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,2,iptty,maxor,overf
     *              ,irdof,icoupl,0,sto5,a_save,nelmdg,accm,mdof_sol)
            else
               call  combine_a(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *              ,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,2,iptty,maxor,overf
     *              ,irdof,icoupl,0,sto5,a_save,nelmdg,accm,mdof_sol)
            end if 
            itert=itert+iter
            itotals=itotals+iter
            deallocate(sto5)
            call storage_derivatives(0,1)
            itert=itert+iter
            minkt=minkt+mink
         end if
         deallocate(a_save)
c     
c     extract proper solution
c     
         call air_rdof(1,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
c     
      else if(irdof.eq.11) then
c     
c     coding for 1dof airflow only
c     
c     
         fdum2=0.0
         do id=1,neq
            np=nelmdg(id)
            i1=nelm(id)+1
            i2=nelm(id+1)
            apiv=a(np-neqp1+nmat(3))
            do ijj=i1,i2
               ij=nelm(ijj)
               aij=a(ijj-neqp1+nmat(3))/apiv
               a(ijj-neqp1)=aij
            enddo
            bp(id+nrhs(1))=bp(id+nrhs(2))/apiv
            fdum2=fdum2+bp(id+nrhs(1))*bp(id+nrhs(1))
         enddo
         fdum=sqrt(fdum2)
         if(fdum.eq.0.0) go to 999
	   mink = n
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.tmch) go to 992
            enddo
            fdum=-1.0
            minkt=minkt+mink
            go to 999
 992        continue
            f0=-1.0
         endif
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
c            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tolls=facr*f0
            tollr=tolls*g3
            if(g3*tmch.gt.tollr) tollr=g3*tmch
            iter=maxsolve      
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(1,1)
            nmat_save=nmat(1)
            nmat(1)=nmat(3)
            if (igauss .gt. 1) then
               call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop
     *              ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
            else
               call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     *              ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
            end if
c***************************Change 3/2/94 gaz       .
            nmat(1)=nmat_save
            call storage_derivatives(0,1)
            itert=itert+iter
            itotals=itotals+iter
            minkt=minkt+mink
c     
c     zero out saturation change
c     
            do i=1,n
               bp(i+nrhs(2))=0.0
            enddo
c     
         end if
      else if(irdof.eq.13) then
c     
c     coding for 1dof waterflow only
c          
         fdum2=0.0
         ibp_max = 0
         bp_max = 0.0
         do id=1,neq
            np=nelmdg(id)
            i1=nelm(id)+1
            i2=nelm(id+1)
            apiv=a(np-neqp1+nmat(1))+adiag_tol     
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
         enddo
       
           
         if(fdum2.eq.0.0) go to 999


         fdum=sqrt(fdum2)
c check if fluxes are small enough to quit
         mink=n
         if(g1.lt.0.) then
c gaz 110519 added coding to insure that at least one iteration is performed           
            if(bp_max.le.abs(g1).and.iad.gt.0) then
               fdum = -999.0
               fdum1 = bp_max
               go to 999
            else
c          bp_tol = bp_max/a(nelmdg(ibp_max)+nmat(1)-neqp1)
               bp_tol = abs(g1)/bp_max*fdum
               f0=-1.0
               tmch = 1.e-20
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


 995        continue
            f0=-1.0
         endif
 996     if(fdum.gt.f0.or.iad.lt.iad_min) then
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
        
c     
c     zero out saturation change
c     
            do i=1,n
               bp(i+nrhs(2))=0.0
            enddo
c     
         end if
      else
         
         
c
c     first call air_rdof to rearrange equations
c  
c gaz debug 112019
c      dumdb1 = 0.0
c      dumdb2 = 0.0
c      if(.not.allocated(dumdb3)) allocate(dumdb3(neq))
c      if(.not.allocated(dumdb4)) allocate(dumdb4(neq))
c      if(.not.allocated(ieosdum)) allocate(ieosdum(neq))
c     do i = 1, neq
c        dumdb1 = dumdb1 +abs(bp(i))
c        dumdb2 = dumdb2 +abs(bp(i+neq))
c        dumdb3(i) = bp(i)
c        dumdb4(i) = bp(i)
c        ieosdum(i) = ieos(i)
c      enddo

         call dual(10)
         call air_rdof(0,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
         if(i_uzsz_re.eq.1) call reorder_uzsz(1)
         if(islord.ne.0) then
            call switch(nmat,nmatb,islord,2,1,nsizea1)
            call switchb(bp,nrhs,nrhsb,islord,2,1,neq)
         end if
c     
         allocate(dumn(100))
         
         call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &        ,ndex,2,dumn(1),dumn(37),dumn(73),0,fdum2)

         deallocate(dumn)
c
         fdum=sqrt(fdum2)
         if(fdum.eq.0.0) go to 999
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.tmch) go to 993
               if(irdof.ne.-11.and.abs(bp(i+nrhs(2))).gt.tmch) go to 993
            enddo
            fdum=-1.0
c gaz 111319 make sure iad_min iterations occur        
c            pause
            if(iad.ge.iad_min) go to 998
 993        continue
            f0=-1.0
         endif
c     
         call explicit(mink)
         minkt=minkt+mink
c 
         if(fdum.gt.f0.or.iad.le.iad_min-1) then
            facr=1.0
c            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
c            tollr=g3*max(tolls,tmch)
            tolls=facr*f0
            tollr=g3*max(tolls,tmch)
c     time=second(ghv)
c     
c     times=second(0.0)
            iter=maxsolve      
            call storage_derivatives(1,1)
c    
            if(irdof.le.0) then 
               if(gdpm_flag.eq.0) then
                  if (igauss .gt. 1) then
                     call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop
     &                    ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
                  else
                     call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     &                    ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
                  end if
               else
                  if (igauss .gt. 1) then
                     call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &                    ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &                    ,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     &                    ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                    ,mdof_sol)
                  else
                     call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat
     &                    ,nrhs,nelm,nelm_primary,nelm_primary,north
     &                    ,tollr,irb,iirb,npvt,gmres,dum,piv
     &                    ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     &                    ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                    ,mdof_sol)
                  end if
               endif
               itert=itert+iter
               itotals=itotals+iter
            else
               allocate(sto5(neq,2))
               if (igauss .gt. 1) then
                  call  rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *                 ,tollr,irb,iirb,npvt,gmres,dum,piv
     *                 ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *                 ,overf,irdof,icoupl,0,sto5,accm)
               else
                  call  rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     *                 ,north,tollr,irb,iirb,npvt,gmres,dum,piv
     *                 ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *                 ,overf,irdof,icoupl,0,sto5,accm)
               end if
               itert=itert+iter
               itotals=itotals+iter
               deallocate(sto5)
            endif
            call storage_derivatives(0,1)
         endif
 998     continue
         call air_rdof(1,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
         if(i_uzsz_re.eq.1) call reorder_uzsz(2)
         if(islord.ne.0) then
            call switch(nmat,nmatb,islord,2,2,nsizea1)
            call switchb(bp,nrhs,nrhsb,islord,2,2,neq)
         end if
         call dual(11)
      end if
      
 999  continue
   
      deallocate(dum)

      
      return
      end

      subroutine add_accumulation(i)        

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsplitts

      implicit none

      integer i, icd, ii1, iz 
      integer jmi, jmia, neqp1
      real*8  sx1d

c first check if saturation change(to 0r from 2-phase) occurs

      neqp1 = neq + 1

      if(i.gt.neq.and.idualp.eq.0) then
         icd=neq
         iz=i - neq
      else
         icd=0
         iz=i
      endif
      ii1=nelm(i-icd)+1
      jmi=nelmdg(i-icd)
      jmia=jmi-neqp1
      sx1d = sx1(i)
      
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)+sk(i)
      a(jmia+nmat(1))=a(jmia+nmat(1))+sx1d*dmpf(i)+dq(i)
	if(jswitch.ne.0) then
         a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*+dmef(i)+dqh(i)
         bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      else if(irdof.ne.13) then
         a(jmia+nmat(2))=a(jmia+nmat(2))+sx1d*dmef(i)+dqh(i)
         a(jmia+nmat(3))=a(jmia+nmat(3))+sx1d*depf(i)+drc(i)
         a(jmia+nmat(4))=a(jmia+nmat(4))+sx1d*deef(i)+deqh(i)
         bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)+qh(i)
      endif
      
      return
      end
      
      subroutine gensl2_explicit        

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comsplitts
      
      implicit none
      
      integer i, ja, mid
      real*8  dels, delp, detdf, dfdum11, dfdum11i, dfdum12, dfdum12i
      real*8  dfdum21, dfdum21i, dfdum22, dfdum22i
      real*8  fdum01, fdum02, sx1d
      real*8  tol_ts, stol, s_up
      parameter(tol_ts=1.d-16,stol=1.0d00,s_up=0.996d00)
c first check if saturation change(to 0r from 2-phase) occurs
c
c perform timestep with explicit flux terms
c

      call airctr(3,0)
      bp = 0.0d00
      do i=1,neq
         call flux_net(i)
         bp(i+nrhs(1)) = bp(i+nrhs(1)) + sk(i)
         bp(i+nrhs(2)) = bp(i+nrhs(2)) + qh(i)             
      enddo

c     loop on grid blocks   
      do mid=1,neq       
c         loop on N-R iteration 
         do ja=1,maxit
c     call EOS
            call accum_derivative(mid,0)
c     form equations
            sx1d = sx1(mid)
            fdum01 = bp(mid+nrhs(1)) + sx1d*deni(mid)
            dfdum11 = sx1d*dmpf(mid)
            if(irdof.ne.13) then
               fdum02 = bp(mid+nrhs(2)) + sx1d*denei(mid)
               dfdum12 = sx1d*dmef(mid)
               dfdum21 = sx1d*depf(mid)
               dfdum22 = sx1d*deef(mid)
            endif
            fdum = sqrt(fdum01*fdum01 + fdum02*fdum02)
c     solve for corrections
            if(irdof.eq.13) then
               delp = fdum01/dfdum11  
               phi(mid) = phi(mid) - delp
               if(fdum.le.abs(tmch)) go to 999
            else
c     need to get this in
               fdum01 = fdum01/dfdum11
               fdum02 = fdum02/dfdum22
               dfdum12 = dfdum12/dfdum11
               dfdum21 = dfdum21/dfdum22
               dfdum11 = 1.d00
               dfdum22 = 1.d00
               detdf = dfdum11*dfdum22 - dfdum12*dfdum21
               dfdum11i = dfdum22/detdf
               dfdum22i = dfdum11/detdf
               dfdum12i = -dfdum12/detdf
               dfdum21i = -dfdum21/detdf
               delp = dfdum11i*fdum01 +dfdum12i*fdum02
               dels = dfdum21i*fdum01 +dfdum22i*fdum02
               phi(mid) = phi(mid) - delp
               s(mid) = s(mid) - dels
               if(fdum.le.abs(tmch)) go to 999
            endif
         enddo
         write(ierr,*) "explicit TS iteration failed : stopping"
         write(ierr,*) "node = ", mid ," fdum01 = ",fdum01        
         stop
 999     continue
      enddo
      isplitts = 2
      return
      end
      subroutine reorder_uzsz(iflg)  
c
c subroutine to reorder uzsz
c on a node by node basis reorders equations      
c  
      use comai
      use comdi
      use combi
      use comgi
      use comei
      use davidi
      implicit none
      integer iflg, idofm, nmatd, neqp1
      integer i, icoupl_iter, id, idl, i1, i2, jj, j, kb
      integer icd,ii1,ii2
c
      real*8 a_line11, a_line21, bpsave
          if(iflg.eq.1) then
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     uzsz  switch columns
c     
            do i = 1,n
               i1=nelm(i)+1
               i2=nelm(i+1)

               do j = i1,i2
                  kb = nelm(j)  
                  if(ieos(kb).ne.1) then 
                     jj = j - neqp1
                     a_line11 = a(jj+nmat(1))
                     a(jj+nmat(1)) = a(jj+nmat(2))
                     a(jj+nmat(2)) = a_line11
                     a_line21 = a(jj+nmat(3))
                     a(jj+nmat(3)) = a(jj+nmat(4))
                     a(jj+nmat(4)) = a_line21                     
                  endif
               enddo 
            enddo
c     
  
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
            do i = 1, n
               if(ieos(i).ne.1) then
                  bpsave = bp(i+nrhs(1))
                  bp(i+nrhs(1)) = bp(i+nrhs(2))
                  bp(i+nrhs(2)) = bpsave    	   
               endif	  
            enddo
        endif    
      return 
      end
      
