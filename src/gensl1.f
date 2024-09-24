      subroutine gensl1
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
CD1  solve the heat and mass transfer equations with full 
CD1  jacobian (unsymmetric, 2n by 2n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gensl1.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:46 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.16   Fri Apr 26 15:20:50 1996   gaz
CD2 coding added for mdnodes
CD2 
CD2    Rev 1.15   Wed Jan 17 08:57:30 1996   hend
CD2 Removed nonconstant parameter statement for IBM
CD2 
CD2    Rev 1.14   Tue Jan 16 14:31:08 1996   hend
CD2 Added capability for 5,6, and n degrees of freedom
CD2 
CD2    Rev 1.13   Fri Jan 12 17:48:52 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.12   Wed Jan 10 12:28:56 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.11   Wed Jan 10 12:01:46 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.10   08/02/95 17:30:40   llt
CD2 moved where sto5 was getting released
CD2 
CD2    Rev 1.9   06/01/95 16:48:06   gaz
CD2 minor corrections f0=-1
CD2 
CD2    Rev 1.8   04/25/95 09:06:34   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.7   03/23/95 18:50:46   gaz
CD2 No change.
CD2 
CD2    Rev 1.6   03/23/95 18:43:56   gaz
CD2 gaz changed bp to abs(bp in one if block
CD2 
CD2    Rev 1.5   03/23/95 18:26:36   gaz
CD2 gaz changes to call to rdof_new and tolerance options
CD2 
CD2    Rev 1.4   03/10/95 10:45:48   llt
CD2 removed unneeded pointers - gaz
CD2
CD2    Rev 1.3   05/11/94 16:14:50   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:43:12   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.1   03/18/94 15:59:32   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:22   pvcs
CD2 original version in process of being certified
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use davidi
      use comflow
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comriv
      use comcouple
      implicit none

      integer ndex(2)
      integer i, id, neqp1, nsizea, nsizea1 
      integer iparchek,kb,jj,i1,i2,i3,i4,j
      real a11     
      real*8 dumn(100),term
      real*8, allocatable :: sto5(:,:)
      real*8, allocatable :: dum(:)
      real*8  facr, fdum2, tollr, tolls, tmch_new, tmch_old, bp_max
      save tmch_old
      parameter(iparchek =0)
      neqp1=neq+1
c test of coefficients
       if(iparchek.eq.1) then
       write(ierr,*) 'fe coef. neighbors'
       write(ierr,*)
       do i = 1,neq
        write(ierr,*)' vol node ', i,' = ', sx1(i)
        i1 = nelm(i)+1
        i2 = nelmdg(i)-1
        do j = i1,i2
         kb = nelm(j)
         i3 = nelmdg(kb)+1
         i4 = nelm(kb+1)
         do jj = i3,i4
         if(nelm(jj).eq.i) then
          iw = istrw(jj-neqp1)
          if(iw.gt.0) then
           a11 = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
          else
           a11 = 0.0
          endif
          write(ierr,355)i,kb,iw,a11
          go to 899
         endif
        enddo
899     continue         
       enddo        
        i1 = nelmdg(i)+1
        i2 = nelm(i+1)        
        do j = i1,i2
         kb = nelm(j)
         iw = istrw(j-neqp1)
          if(iw.gt.0) then
           a11 = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
          else
           a11 = 0.0
          endif
         write(ierr,355)i,kb,iw,a11
        enddo        
       enddo   
355   format(i6,1x,i6,1x,i6,1x,f12.3)        
       stop
       endif    
c     zero out arrays       
c
c  check for possible stopping on variable changes
c
       
      if(iad.ge.1) then
c gaz 113021          
c       if(nr_stop.eq.2) then
        if(nr_stop.gt.2) then        
          fdum=-1.0
        go to 999
       endif
      endif
c           
      do i=1,neq
         bp(i+nrhs(1))=0.0
         bp(i+nrhs(2))=0.0
      enddo
      nsizea1=nelm(neqp1)-neqp1
      nsizea=4*nsizea1
      do i=1,nsizea
         a(i)=0.0
      enddo
      fdum2=0.
      do id=1,neq
        if(iriver.eq.2.and.id.gt.neq_primary) then
         call geneq1_well(id)
        else if(ianpe.ne.0) then 
         call geneq1_ani(id)
        else
         call geneq1(id)
        endif
      enddo
      do id=1,neq
         if(ps(id).le.0.0) then
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0
         endif
      enddo

c     call md_modes to complete equations for multiply defined nodes
c don't need to call geneqmdnode (through md_nodes) gaz
c     if(imdnode.ne.0) call md_nodes(3,0,0)

      call dual(10)
c
c possible correction for rate-limited gdpm model
c
      call gdpm_corr(1)

      if(islord.ne.0) call switch(nmat,nmatb,islord,2,1,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,2,1,neq)
c    
c normalize equations 
c     
      call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &     ,ndex,2,dumn(1),dumn(37),dumn(73),0,fdum2)
      if(ndex(1).lt.0) then
         write(ierr,*) 
     & '>>> cannot normalize, stopping, i=',abs(ndex(1)),'<<<'
         if(iout.ne.0) then
          write(iout,*) 
     & '>>> cannot normalize, stopping, i=',abs(ndex(1)),'<<<'
         endif
         if(iptty.ne.0) then
          write(iptty,*) 
     & '>>> cannot normalize, stopping, i=',abs(ndex(1)),'<<<'
         endif 
c gaz 112521 added more information   
       if(iout.ne.0) then  
        id = abs(ndex(1))
        write(iout,*) 'id ',id,' vol ', sx1(id)
        write(iout,*)
        write(iout, 9010) l, day
        write(iout, 9011) id,cord(id,1),cord(id,2),
     &           cord(id,3)
        write(iout, 9012) phi(id), t(id), s(id)  
        write(iout, 9015) ps(id),pnx(id)*1.d-6,ieos(id),iad,izonef(id)
       endif 
       if(iptty.ne.0) then  
        id = abs(ndex(1))
        write(iptty,*)
        write(iptty, 9010) l, day
        write(iptty, 9011) id,cord(id,1),cord(id,2),
     &           cord(id,3)
        write(iptty, 9012) phi(id), t(id), s(id) 
        write(iptty, 9015) ps(id),pnx(id)*1.d-6,ieos(id),iad,izonef(id)
       endif  
       if(iptty.ne.0) then  
        id = abs(ndex(1))
        write(ierr,*)
        write(ierr, 9010) l, day
        write(ierr, 9011) id,cord(id,1),cord(id,2),
     &           cord(id,3)
        write(ierr, 9012) phi(id), t(id), s(id) 
        write(ierr, 9015) ps(id),pnx(id)*1.d-6,ieos(id),iad,izonef(id)
       endif  
 9010 format (1p,'time step = ', i8, ' time step size = ', g21.14)
 9011 format (1p,'out of bounds : node ', i8,
     &     ' x = ', g12.4, ' y = ', g12.4, ' z = ', g12.4)
 9012 format (1p,' p = ', g16.9, ' t =', g16.9, ' s = ', g16.9)       
 9015 format (1p,' porosity = ',g16.8,' permx ',g16.8,
     & ' phase state ',i3, ' iter ',i3,' zone ', i7)       
         stop
      endif

      fdum=sqrt(fdum2)
      if(fdum.eq.0.0) go to 999
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
	  bp_max = 0.0
         do i=1,neq
            bp_max= max(abs(bp(i+nrhs(1))),abs(bp(i+nrhs(2))),bp_max) 
         enddo
c gaz 113121 modified to work with changes in nr_stop_ctr1         
c	 if(bp_max.lt.tmch.and.nr_stop.eq.0) then         
	 if(bp_max.lt.tmch.and.nr_stop.le.2) then
         fdum=-1.0
         go to 999
       else
         f0=-1.0
	   if(iad.eq.1) tmch_old = bp_max
	   tmch_new = min(bp_max,tmch_old)
	   tmch_old = tmch_new
	 endif
      endif
      if(f0.gt.0.0d00) then
         if(fdum.le.f0.and.iad.ne.0) goto 999
         facr=1.0d00
         tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
         tollr=tolls*g3
         if(tmch*g3.gt.tollr) tollr=tmch*g3 
      else if(f0.le.0.0d00) then
         tollr=tmch_new*g3    
      endif

c     time=second(ghv)
c
c      times=second(0.0)
      iter=maxsolve 
      call storage_derivatives(1,1)

      if(irdof.eq.0) then
c
c       full solution
c

         allocate(dum(neq*8))
         if(gdpm_flag.eq.0) then
            if (igauss .gt. 1) then
               call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,
     &              north,tollr,irb,iirb,npvt,gmres,dum,piv,
     &              h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
            else
               call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,
     &              north,tollr,irb,iirb,npvt,gmres,dum,piv,
     &              h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
            end if
         else
            if (igauss .gt. 1) then
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs,
     &              nelm,nelm_primary,nop,north,tollr,irb,iirb,
     &              npvt,gmres,dum,piv,h,c,ss,g,y,iter,iback,2,
     &              iptty,maxor,igdpm,maxgdpmlayers,ngdpm_layers,
     &              nelmdg,accm,mdof_sol)
            else
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat,nrhs,
     &              nelm,nelm_primary,nelm_primary,north,tollr,irb,iirb,
     &              npvt,gmres,dum,piv,h,c,ss,g,y,iter,iback,2,iptty,
     &              maxor,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,
     &              accm,mdof_sol)
            end if
         endif
	   deallocate(dum)
         itert=itert+iter
         itotals=itotals+iter
         minkt=minkt+mink
      else

         allocate(sto5(neq,2))
         allocate(dum(neq*8))
         if (igauss .gt. 1) then
            call rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm) 
         else
            call rdof_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *           ,overf,irdof,icoupl,0,sto5,accm) 
         end if
         deallocate(sto5)
         deallocate(dum)
         
         itert=itert+iter
         itotals=itotals+iter
         minkt=minkt+mink
      endif

      call storage_derivatives(0,1)
      if(islord.ne.0) call switch(nmat,nmatb,islord,2,2,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,2,2,neq)
      call dual(11)
      
 999  continue
      
      return
      end
c.............................................................

      subroutine extract_adiag
c s kelkar 3/20/12. extract diagonal elements of 'a'
      use comai, only : neq, ierr
      use combi, only : nelmdg, sx1
      use comei, only : a
      use comdi, only : ps, pnx
      use davidi, only : nmat
      implicit none

      integer idisk,j

c      real*8 adiagsk(2*neq)
      write(ierr,99) 
      do idisk=1,neq
c         adiagsk(idisk) = a(nelmdg(idisk)-neq-1+nmat(1))
c         adiagsk(idisk+neq) = a(nelmdg(idisk)-neq-1+nmat(4))
         write(ierr,100) idisk, (a(nelmdg(idisk)-neq-1+nmat(j)),j = 1,4)
     &   ,ps(idisk),sx1(idisk), pnx(idisk)*1.e-6
      enddo
99    format('node',t14,'a11',t29,'a12',t44,'a21',t59,'a22',t74,'ps',
     &        t89,'sx1',t104,'pnx')
100   format(i6,1p,7(1x,g14.6))

      return
      end 
c...........................................................

      subroutine extract_nodal_equations(idisk)
c s kelkar 3/20/12. extract rows of a corrosponding 
c to the node is idisk from the 'a' matrix
      use comai, only : neq, ierr, l
      use combi, only : nelmdg, nelm
      use comei, only : a
      use comgi, only : bp
      use davidi, only : nmat, nrhs
      implicit none

      integer idisk,iblock, jdisk, indsk,ii1,ii2, j,k,iq
      integer neqp1
      integer kbsk(20), inda, kbp

      real*8 aeqsk(20,4)

      neqp1=neq+1

      ii1=nelm(idisk)+1
      ii2=nelm(idisk+1)

      iq=0
      do j=ii1,ii2
         iq=iq+1
         kbsk(iq)=nelm(j)
         do k=1,4
            inda = j-neqp1+nmat(k)
            aeqsk(iq,k)=a(inda)
         enddo
      enddo

      write(ierr,*)' Time step is:', l
      write(ierr,*)' Node I: ', idisk
      write(ierr,*)(kbsk(j),j=1,iq)
 9991 format(30(2x,i10))
      do k=1,4
         kbp=1
         if(k.gt.2)kbp=2
         write(ierr,9992)(aeqsk(j,k),j=1,iq),bp(nrhs(kbp)+idisk)
      enddo
 9992 format(30(2x,g12.5))
         

      return
      end 
c...........................................................
