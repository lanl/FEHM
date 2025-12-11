      subroutine gentdp(iz,igrp,matnum,spec_num,sia_iter,
     2     ispecies,ncsp)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine sets up for and calls the necessary routines to 
CD1  solve the tracer equations with full jacobian (unsymmetric, 2n
CD1  by 2n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gentdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:16   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.13   Fri May 24 09:55:24 1996   hend
CD2 Updated trac for mdnodes
CD2 
CD2    Rev 1.12   Thu Apr 25 12:17:32 1996   hend
CD2 Added 2 necessary lines after calls to coneq1
CD2 
CD2    Rev 1.11   Thu Mar 21 13:18:24 1996   hend
CD2 Fixed for new coneq1 calling sequence
CD2 
CD2    Rev 1.10   Mon Mar 04 15:04:58 1996   hend
CD2 Fixed for new solver
CD2 
CD2    Rev 1.9   Fri Jan 12 17:49:16 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.8   Wed Jan 10 12:26:18 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.7   08/18/95 10:18:10   llt
CD2 full_solution was already defined, removed for cray
CD2 
CD2    Rev 1.6   08/07/95 13:50:46   awolf
CD2 New call requires matnum and nspecies and iz
CD2 Check is iz = 1 or 0 and devides routine work based on that 
CD2 parameter. dpdpta calls now require and arguement "call dpdpta()
CD2 
CD2    Rev 1.5   03/24/95 00:07:52   gaz
CD2 gaz took out reference to fdum1(not used)
CD2 
CD2    Rev 1.3   05/11/94 16:07:16   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:41:10   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:47:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:30   pvcs
CD2 original version in process of being certified
CD2 
CD2 c 18-mar-94 gaz added mmget..."b"
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.4 Solute-transport equations
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
c Note dpdp is only implemented for one species 

      use comcouple
      use comrxnb
      use comrxni
      use davidi
      use comji
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai

      implicit none

      integer igrp, iz, matnum, ncsp, sia_iter, neqp1, idof_sol
      integer i
      real*8 fdm, fdum2, tolls, tollr
      real*8, allocatable :: sto5(:)
      real*8, allocatable :: dum(:)      
      real*8, allocatable :: diag(:,:)
      real*8, allocatable :: diagi(:,:)
      real*8, allocatable :: rdum(:)
      integer iflg
      parameter(fdm=20.0)
      integer iter_counter,spec_num,full_solution,tol_value
      integer spec_numf
      integer matnumf
      integer spec_numm
      integer matnumm
      integer ispecies
      integer kgmres_sol
      integer totmat
      
      if(iz.eq.0) then 
         call thermc(neq)
c     
c     adjust storage
c     
         
         neqp1=neq+1
         fdum2=0.
         matnumf = matnum
         spec_numf = spec_num
         call coneq1(0,matnumf,spec_numf)
         matnumm = matnum+ncsp*2+1
         spec_numm = spec_num+1
         call coneq1(neq,matnumm,spec_numm)
c        if (imdnode.ne.0) call coneq1mdnode(matnum,spec_num)
c     
c     now get coupling for influence terms
c     
         call dpdpta(ispecies,spec_numf,matnumf,spec_numm,matnumm)
         
      elseif(iz.eq.1) then
c     
c     normalize the equations
c     
         totmat = ncsp*2
         idof_sol = totmat               
         allocate(diag(idof_sol,idof_sol))
         allocate(diagi(idof_sol,idof_sol))
         allocate(rdum(idof_sol))
         call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg,
     &        nmatb,idof_sol,diag,diagi,rdum,0,fdum2)
         deallocate(diag)
         deallocate(diagi)
         deallocate(rdum)
         if(nmatb(1).lt.0 ) then
            if (iout .ne. 0) then
               write(iout,*) '   '
               write(iout, 100)
               write(iout,*) '   '
               write(iout,'(i8,3g12.3)')
     &              abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
            end if
            if(iptty.gt.0)then
               write(iptty,*) '   '
               write(iptty, 100)
               write(iptty,*) '   '
               write(iptty,'(i8,3g12.3)') 
     &              abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
            endif
            iad=maxit
            return
         endif
 100     format ('*singular matrix found during normalization *')

         fdum=sqrt(fdum2)
         if(sia_iter.eq.0)fzero(igrp)=fdum
         tolls = max(fzero(igrp),min(g1*fdum,g2*fdum2))
         tollr=g3*tolls
         iter=maxsolve
         call storage_derivatives(1,1)

         if(irdof.eq.0.or.irdof.ne.0) then
            if (idof_sol.lt.6) then
               allocate(dum(neq*idof_sol*4))
            else if (idof_sol.eq.6) then
               allocate(dum(neq*idof_sol*5))
            else
               allocate(dum(neq*idof_sol*idof_sol))
            endif

            if (igauss .gt.1) then
               call solve_new(neq,a,b,bp,nmat_sol,nb,nrhs_sol,
     1              nelm,nop,north,tollr,irb,iirb,npvt,gmres,
     2              dum,piv,h,c,ss,g,y,iter,iback,idof_sol,
     3              iptty,maxor,accm)
            else
               call solve_new(neq,a,b,bp,nmat_sol,nmat_sol,nrhs_sol,
     1              nelm,nelm,north,tollr,irb,iirb,npvt,gmres,
     2              dum,piv,h,c,ss,g,y,iter,iback,idof_sol,
     3              iptty,maxor,accm)
            end if
            itert=itert+iter
            minkt=minkt+mink
         else
            allocate(sto5(neq*2))
            if (igauss .gt.1) then
               call  rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *              ,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *              ,overf,irdof,icoupl,0,sto5,accm)
            else
               call  rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *              ,tollr,irb,iirb,npvt,gmres,dum,piv
     *              ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *              ,overf,irdof,icoupl,0,sto5,accm)
            end if
            deallocate(sto5)
            itert=itert+iter
            minkt=minkt+mink
         endif
         deallocate(dum)
         call storage_derivatives(0,1)
 999     continue
      endif
      
      return
      end
