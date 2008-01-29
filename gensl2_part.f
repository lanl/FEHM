      subroutine gensl2_part(zone)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To solve isothermal air-water eqs with full jacobian.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/gensl2_part.f_a  $
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
!D3 2.5.2 Solve nonlinear equation set at each time step
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************
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
      use com_part
      implicit none
      real*8, allocatable :: sto5(:,:)
      integer icode
      integer index(2)
      integer neqp1
      integer i
      integer id
      integer nsizea
      integer nsizea1
      integer nmat_save
      integer np, i1, i2, ij, ijj
      real*8 aij, apiv
      real*8 fdum2
      real*8 facr
      real*8 tolls
      real*8 tollr
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: a_save(:)
      real*8 a11, a22, adiag_tol
      parameter(adiag_tol=0.d-90)
      integer jj
      integer zone

      neqp1=neq_part(zone)+1
c     zero out arrays
      do 10 i=1,neq_part(zone)
         bp_part(zone,i+nrhs_part(zone,1))=0.0
         bp_part(zone,i+nrhs_part(zone,2))=0.0
 10   continue
      nsizea1=nelm_part(zone,neqp1)-neqp1
      if(irdof.ne.13) then
        nsizea=4*nsizea1
        allocate(dum(8*neq_part(zone)))
      else
        nsizea=nsizea1
        allocate(dum(4*neq_part(zone)))
      endif
c
      a=0.0
c
      fdum2=0.
      do 101 id=1,neq_part(zone)
c     
c     decide on equation type
c     
         if(ps(index_part(zone,id)).gt.0.0) then
           if(ianpe.ne.0) then
            call geneq2_ani(id)
            call add_accumulation_part(id,zone)
           else if(ifree.ne.0) then
            call geneq2_wtsi(id)
            call add_accumulation_part(id,zone)
           else
            call geneq2_part(id,zone)
            call add_accumulation_part(id,zone)
           endif
         else if(ps(index_part(zone,id)).le.0.0) then
            if(irdof.ne.13) then
               a(nelmdg_part(zone,id)-neqp1+nmat_part(zone,2))=0.0d00     
               a(nelmdg_part(zone,id)-neqp1+nmat_part(zone,3))=0.0d00    
               a(nelmdg_part(zone,id)-neqp1+nmat_part(zone,1))=
     &                              sx1(index_part(zone,id))
               a(nelmdg_part(zone,id)-neqp1+nmat_part(zone,4))=
     &                              sx1(index_part(zone,id))
               bp_part(zone,id+nrhs_part(zone,1))=0.0d00
               bp_part(zone,id+nrhs_part(zone,2))=0.0d00
            else
               a(nelmdg_part(zone,id)-neqp1)=
     &                              sx1(index_part(zone,id))
               bp_part(zone,id)=0.0d00
            endif
         endif
c
c check for singularity
c
c              a11 = abs(a(nelmdg(id)-neqp1+nmat(1)))    
c              a22 = abs(a(nelmdg(id)-neqp1+nmat(4)))     
c              if(a11.lt.adiag_tol) then
c               do jj=nelm(id)+1,nelm(id+1)
c                a(jj-neqp1+nmat(1)) =0.0d00
c                a(jj-neqp1+nmat(2)) =0.0d00
c                a(jj-neqp1+nmat(3)) =0.0d00
c               enddo                         
c               a(nelmdg(id)-neqp1+nmat(1)) = 1.0d00
c               bp(id+nrhs(1))=0.0d00
c              else if(a22.lt.adiag_tol) then
c               do jj=nelm(id)+1,nelm(id+1)
c                a(jj-neqp1+nmat(4)) =0.0d00
c                a(jj-neqp1+nmat(2)) =0.0d00
c                a(jj-neqp1+nmat(3)) =0.0d00
c               enddo                         
c               a(nelmdg(id)-neqp1+nmat(4)) = 1.0d00
c               bp(id+nrhs(2))=0.0d00
c              endif
 101  continue
c
c add correction for saturations over 1.0
c
       if(iflux_ts.ne.0) call cascade_sat(2)
c
      if(irdof.eq.14) then
c
c     first call air_rdof to rearrange equations
c     
         call air_rdof(0,irdof,0,nelm,nmat,
     &                 nrhs,a,bp,sx1)
c     
         allocate(a_save(nsizea1))

c
c sort out degenerate cases
c
c        do i=1,neq
c         i1=nelm(i)+1
c         i2=nelm(i+1)
c         a11=a(nelmdg(i)-neqp1+nmat(1))
c         if(a11.eq.0.0) then
c          a(nelmdg(i)-neqp1+nmat(1))=1.0
c          bp(i+nrhs(1))=0.0
c          do jj=i1,i2
c           a(jj-neqp1+nmat(1))=0.0
c          enddo
c         endif
c         a22=a(nelmdg(i)-neqp1+nmat(4))
c         if(a22.eq.0.0) then
c          a(nelmdg(i)-neqp1+nmat(4))=1.0
c          bp(i+nrhs(2))=0.0
c          do jj=i1,i2
c           a(jj-neqp1+nmat(4))=0.0
c          enddo
c         endif
c        enddo
c
         fdum2=0.0
         do i=1,neq_part(zone)
         a11=max(a(nelmdg_part(zone,i)-neqp1+nmat_part(zone,1)),
     &           a(nelmdg_part(zone,i)-neqp1+nmat_part(zone,2)))
         a22=max(a(nelmdg_part(zone,i)-neqp1+nmat_part(zone,3)),
     &           a(nelmdg_part(zone,i)-neqp1+nmat_part(zone,4)))
          i1=nelm_part(zone,i)+1
          i2=nelm_part(zone,(i+1))
           if(a11.ne.0.0) then
            do jj=i1,i2
             a(jj-neqp1+nmat_part(zone,1))=
     &                     a(jj-neqp1+nmat_part(zone,1))/a11
             a(jj-neqp1+nmat_part(zone,2))=
     &                     a(jj-neqp1+nmat_part(zone,2))/a11
            enddo
            fdum2=fdum2+bp_part(zone,i+nrhs_part(zone,1))**2
            bp_part(zone,i+nrhs_part(zone,1))=
     &                bp_part(zone,i+nrhs_part(zone,1))/a11
           endif
           if(a22.ne.0.0) then
            do jj=i1,i2
             a(jj-neqp1+nmat_part(zone,3))=
     &                       a(jj-neqp1+nmat_part(zone,3))/a22
             a(jj-neqp1+nmat_part(zone,4))=
     &                       a(jj-neqp1+nmat_part(zone,4))/a22
            enddo
            bp_part(zone,i+nrhs_part(zone,2))=
     &                 bp_part(zone,i+nrhs_part(zone,2))/a22
           fdum2=fdum2+bp_part(zone,i+nrhs_part(zone,2))**2
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
            do i=1,neq_part(zone)
               if(abs(bp_part(zone,i+nrhs_part(zone,1)))
     &                          .gt.tmch) go to 991
               if(abs(bp_part(zone,i+nrhs_part(zone,2)))
     &                          .gt.tmch) go to 991
            enddo
            fdum=-1.0

            fdum=-1.0
            minkt=minkt+mink
            go to 999
 991        continue
            f0=-1.0
         endif
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tollr=tolls*g3
            if(tmch.gt.tollr) tollr=tmch
            iter=maxsolve      
            call storage_derivatives(1,1)
             allocate(sto5(neq_part(zone),2))
               call  combine_a(neq,a,b,bp,nmat
     *          ,nb,nrhs,nelm,nop,north
     *          ,tollr,irb,iirb,npvt,gmres,dum,piv
     *          ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     *          ,overf,irdof,icoupl,0,sto5,a_save,nelmdg,accm,mdof_sol) 
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
            call air_rdof(1,irdof,0,nelm,
     &            nmat,nrhs,a,bp,sx1)
c     
      else if(irdof.eq.11) then
c     
c     coding for 1dof airflow only
c     
c     
         fdum2=0.0
         do id=1,neq_part(zone)
            np=nelmdg_part(zone,id)
            i1=nelm_part(zone,id)+1
            i2=nelm_part(zone,id+1)
            apiv=a(np-neqp1+nmat_part(zone,3))
            do ijj=i1,i2
               ij=nelm_part(zone,ijj)
               aij=a(ijj-neqp1+nmat_part(zone,3))/apiv
               a(ijj-neqp1)=aij
            enddo
            bp_part(zone,id+nrhs_part(zone,1))=
     &                  bp_part(zone,id+nrhs_part(zone,2))/apiv
            fdum2=fdum2+bp_part(zone,id+nrhs_part(zone,1))*
     &                  bp_part(zone,id+nrhs_part(zone,1))
         enddo
         fdum=sqrt(fdum2)
         if(fdum.eq.0.0) go to 999
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            do i=1,neq_part(zone)
               if(abs(bp_part(zone,i+nrhs_part(zone,1)))
     &                          .gt.tmch) go to 992
            enddo
            fdum=-1.0
            minkt=minkt+mink
            go to 999
 992        continue
            f0=-1.0
         endif
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tollr=tolls*g3
            if(g3*tmch.gt.tollr) tollr=g3*tmch
            iter=maxsolve      
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(1,1)
            nmat_save=nmat_part(zone,1)
            nmat_part(zone,1)=nmat_part(zone,3)
            call solve_new(neq,a,b,bp,nmat
     *           ,nb,nrhs,nelm
     *           ,nelm,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
c***************************Change 3/2/94 gaz       .
            nmat_part(zone,1)=nmat_save
            call storage_derivatives(0,1)
            itert=itert+iter
            itotals=itotals+iter
            minkt=minkt+mink
c     
c     zero out saturation change
c     
            do i=1,n
               bp_part(zone,i+nrhs_part(zone,2))=0.0
            enddo
c     
         end if
      else if(irdof.eq.13) then
c     
c     coding for 1dof waterflow only
c     
c     
         fdum2=0.0
         do id=1,neq_part(zone)
            np=nelmdg_part(zone,id)
            i1=nelm_part(zone,id)+1
            i2=nelm_part(zone,id+1)
            apiv=a(np-neqp1+nmat_part(zone,1))
            do ijj=i1,i2
               ij=nelm_part(zone,ijj)
               aij=a(ijj-neqp1+nmat_part(zone,1))/apiv
               a(ijj-neqp1)=aij
            enddo
            bp_part(zone,id+nrhs_part(zone,1))=
     &                   bp_part(zone,id+nrhs_part(zone,1))/apiv
            fdum2=fdum2+bp_part(zone,id+nrhs_part(zone,1))*
     &               bp_part(zone,id+nrhs_part(zone,1))
         enddo
         fdum=sqrt(fdum2)
         fdum_part(zone)=fdum
         if(fdum.eq.0.0) go to 999
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
         if(fdum1.lt.0.0.and.iad.ne.0) then
            do i=1,neq_part(zone)
               if(abs(bp_part(zone,i+nrhs_part(zone,1)))
     &                     .gt.tmch) go to 995
            enddo
            fdum_part(zone)=-1.0
            go to 999
 995        continue
            f0=-1.0
         endif
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tollr=g3*max(tolls,tmch)
            iter=maxsolve      
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(1,1)
            if(gdpm_flag.eq.0) then 
              call solve_new(neq_part(zone),a,b,bp_part(zone,:)
     &             ,nmat_part(zone,:),nb_part(zone,:)
     &             ,nrhs_part(zone,:),nelm_part(zone,:)
     &             ,nelm_part(zone,:),north
     &             ,tollr,irb,iirb
     &             ,npvt_part(zone,:),gmres,dum,piv
     &             ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)

            else
             call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &            ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &            ,npvt,gmres,dum,piv
     &            ,h,c,ss,g,y,iter,iback,1,iptty,maxor
     &            ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &            ,mdof_sol)
            endif
c***************************Change 3/2/94 gaz       .
            call storage_derivatives(0,1)
c        endif
               itert_part(zone)=itert_part(zone)+iter
               itotals_part(zone)=itotals_part(zone)+iter
               minkt=minkt+mink
c     
c     zero out saturation change
c     
            do i=1,neq_part(zone)
               bp_part(zone,i+nrhs_part(zone,2))=0.0
            enddo
c     
         end if
      else
         
         
c
c     first call air_rdof to rearrange equations
c     
         
         call dual(10)
         call air_rdof(0,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
         if(islord.ne.0) then
            call switch(nmat,nmatb,islord,2,1,nsizea1)
            call switchb(bp,nrhs,nrhsb,islord,2,1,neq)
         end if
c     
       allocate(dumn(100))

       call normal_dof(neq_part(zone),a,bp_part(zone,:)
     &  ,nelm_part(zone,:),nmat_part(zone,:),nrhs_part(zone,:)
     &  ,nelmdg_part(zone,:)
     &  ,index,2,dumn(1),dumn(37),dumn(73),0,fdum2)

       deallocate(dumn)
c
         fdum=sqrt(fdum2)
         fdum_part(zone)=fdum
         if(fdum.eq.0.0) go to 999
         if(iad.eq.0) then
            f0=max(fdum*epe,tmch)
         endif
c        if(fdum1.lt.0.0.and.iad.ne.0) then
         if(iad.ne.0) then
            do i=1,neq_part(zone)
               if(abs(bp_part(zone,i+nrhs_part(zone,1)))
     &                    .gt.tmch) go to 993
               if(irdof.ne.-11.and.
     &              abs(bp_part(zone,i+nrhs_part(zone,2)))
     &                         .gt.tmch) go to 993
            enddo
            fdum_part(zone)=-1.0

            go to 998
 993        continue
            f0=-1.0
         endif
c
         call explicit(mink)
         minkt=minkt+mink
c
         if(fdum.gt.f0.or.iad.eq.0) then
            facr=1.0
            tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
            tollr=g3*max(tolls,tmch)
c     time=second(ghv)
c     
c     times=second(0.0)
            iter=maxsolve      
            call storage_derivatives(1,1)
c    
            if(irdof.le.0) then 
             if(gdpm_flag.eq.0) then 
              call solve_new(neq_part(zone),a,b,bp_part(zone,:)
     &             ,nmat_part(zone,:),nb_part(zone,:)
     &             ,nrhs_part(zone,:),nelm_part(zone,:)
     &             ,nelm_part(zone,:),north
     &             ,tollr,irb,iirb
     &             ,npvt_part(zone,:),gmres,dum,piv
     &             ,h,c,ss,g,y,iter,iback,2,iptty,maxor,accm)
             else
             call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &            ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &            ,npvt,gmres,dum,piv
     &            ,h,c,ss,g,y,iter,iback,2,iptty,maxor
     &            ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &            ,mdof_sol)
             endif
               itert_part(zone)=itert_part(zone)+iter
               itotals_part(zone)=itotals_part(zone)+iter
            endif
            call storage_derivatives(0,1)
         endif
 998     continue
         call air_rdof(1,irdof,0,nelm,nmat,nrhs,a,bp,sx1)
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
      subroutine add_accumulation_part(i,zone)        

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
      use com_part

      implicit none

      integer zone
      integer i, icd, ii1, iz 
      integer jmi, jmia, neqp1
      real*8  sx1d

c first check if saturation change(to 0r from 2-phase) occurs

      neqp1 = neq_part(zone) + 1

      if(i.gt.neq_part(zone).and.idualp.eq.0) then
         icd=neq_part(zone)
         iz=i - neq_part(zone)
      else
         icd=0
         iz=i
      endif
      ii1=nelm_part(zone,i-icd)+1
      jmi=nelmdg_part(zone,i-icd)
      jmia=jmi-neqp1
      sx1d = sx1(index_part(zone,i))

      bp_part(zone,iz+nrhs_part(zone,1))=
     &            bp_part(zone,iz+nrhs_part(zone,1))+sx1d*
     &            deni(index_part(zone,i))+sk(index_part(zone,i))
      a(jmia+nmat_part(zone,1))=a(jmia+
     &              nmat_part(zone,1))+sx1d*
     &              dmpf(index_part(zone,i))+dq(index_part(zone,i))
      if(irdof.ne.13) then
       a(jmia+nmat_part(zone,2))=
     &              a(jmia+nmat_part(zone,2))+sx1d*
     &              dmef(index_part(zone,i))+dqh(index_part(zone,i))
       a(jmia+nmat_part(zone,3))=
     &              a(jmia+nmat_part(zone,3))+sx1d*
     &              depf(index_part(zone,i))+drc(index_part(zone,i))
       a(jmia+nmat_part(zone,4))=
     &              a(jmia+nmat_part(zone,4))+sx1d*
     &              deef(index_part(zone,i))+deqh(index_part(zone,i))
       bp_part(zone,iz+nrhs_part(zone,2))=
     &      bp_part(zone,iz+nrhs_part(zone,2))+sx1d*
     &      denei(index_part(zone,i))+qh(index_part(zone,i))
      endif

      return
      end

   
      

