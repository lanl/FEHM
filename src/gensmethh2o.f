      subroutine gensmethh2o    
C***********************************************************************
CD1
CD1  PURPOSE 
CD1
CD1  This subroutine solves the non-isothermal water-methane equations 
CD1  full jacobian (unsymmetric, 3n by 3n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 
CD2    Rev 1.0   06/20/02 10:24:20   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
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

      use comai
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comhi
      use davidi
      use comdti
      use comcouple
      use commeth
      implicit none

      integer nmatd, index(6)
      integer i, j, icoupl_iter, id, idl, idofm 
      integer jj,ii,idiag
      integer iprint,inorm, j_temp,i_temp
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls, fdum_temp
      parameter (iprint=0,inorm = 0)
      
      neqp1=neq+1
      nmatd=nelm(neqp1)-neqp1
      nsizea=idof*idof*nmatd
      do i=1,nsizea
         a(i)=0.0
      enddo

      do i=1,neq
       do j=1,idof
         bp(i+nrhs(j))=0.0
       enddo
      enddo

c
c    call to icectr to get permeability reduction factor
c
        call icectr(10,0)
c
c    call to ther_meth_h2o to get rel perms for water and methane
c
        call ther_meth_h2o(0,0)
c 
c    hydrate componet
c
c  no rearrangement needed for hydrate since full 6 by 6
c  ie geneq_hyd assumes structure
c     
c     add contribution heat-mass of hydrate
c     C of water equation
c     C of energy equation
c     care not to use derivative terms already used by water
c     
         call ther_meth_h2o(3,0)
c     
      do id = 1,neq
          call geneq_hyd(id)
      enddo
c  
c
c     water component (gaz)
c
c     
c     add contribution heat-mass of water
c     C of water equation
c     C of energy equation
c     
         call ther_meth_h2o(1,0)
c     
      do id = 1,neq
          call geneq_h2o(id)
      enddo
c  
c     add contribution methane            
c     C of methane equation
c     C of energy equation (contribution)
c     new 10-28-2002 gaz
         call ther_meth_h2o(2,0)
c     
c     
c     decide on equation type)
c     
      do id = 1,neq
          call geneq_meth(id)
      enddo
c
c    call to ther_meth_h2o to release memory for rel perms
c
        call ther_meth_h2o(-1,0)
c
c   enforce equilibrium conditions as necessary
c
      call methh2o_combine(1,idofm)

c     print out before manipulating or normalizing arrays
      if(iprint.ne.0) then
       id =1
       write(*,*) 'node = ',id
       write(*,*) 'residuals'
       do jj=1,idof
        write(*,*) jj
        write(*,*) bp(id+nrhs(jj))
        write(*,*) 'diagonal of a'
	  idiag=nelmdg(id)-neqp1
        write(*,200) (a(idiag+nmat(6*(jj-1)+ii)),ii=1,6)
       enddo
       write(*,*) 
       id =neq
       write(*,*) 'node = ',id
       write(*,*) 'residuals'
       do jj=1,idof
        write(*,*) jj
        write(*,*) bp(id+nrhs(jj))
        write(*,*) 'diagonal of a'
	   idiag=nelmdg(id)-neqp1
        write(*,200) (a(idiag+nmat(6*(jj-1)+ii)),ii=1,6)
       enddo
       if(iprint.eq.2) stop
      endif
200    format(1x,1p,6(g14.4))

      if(islord.ne.0) call switch(nmat,nmatb,islord,idofm,1,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,idofm,1,neq)

c
c     normalize the equations
c     
      allocate(dumn(108))
      
c  gaz 2-27-03
       if(inorm.ne.0) then
        call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &      ,index,idofm,dumn(1),dumn(17),dumn(33),0,fdum2)
       else
	 fdum2 = 0.0
	   index(1) = 0
         do i=1,neq
          do j=1,idofm
           fdum2 = fdum2+bp(i+nrhs(j))*bp(i+nrhs(j))
          enddo
         enddo
       endif
c     call mmrelblk ("dumn", "dpdp", ipdumn, icode)	
      deallocate(dumn)
c     
c     return and/or stop if singular
c     
      if(index(1).lt.0 ) then
         if (iout .ne. 0) then
            write(iout,*)
            write(iout,100) 
            write(iout,*)
         end if
c     write(iout,'(i8,3g12.3)')
c     &        abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         if(iptty.gt.0)then
            write(iptty,*)
            write(iptty,100) 
            write(iptty,*)
c     write(iptty,'(i8,3g12.3)') 
c     &        abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         endif
         iad=maxit
         return
      endif
 100  format('* singular matrix found during normalization *')
c     find residual
c     
      fdum=sqrt(fdum2)
      if(fdum.eq.0.0) go to 999
      mink=n
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
	   fdum_temp = 0.0d0
	   j_temp = 0
	   i_temp = 0
         do i=1,neq
          do j=1,idofm
            if(abs(bp(i+nrhs(j))).gt.fdum_temp) then
	       fdum_temp = abs(bp(i+nrhs(j)))
	       i_temp = i
	       j_temp = j
	      endif
          enddo
         enddo
	   if(fdum_temp.lt.tmch) then
          fdum=-1.0
          go to 999
         else
          f0=-1.0  
	   endif  
      endif
      if(f0.gt.0.0) then
         if(fdum.le.f0.and.iad.ne.0) goto 999
         facr=1.0
         tolls=max(facr*f0,min(g1*fdum,g2*fdum2))
         tollr=tolls*g3
         if(g3*tmch.gt.tollr) tollr=g3*tmch
       else if(f0.le.0.0) then
         tollr=tmch *g3  
       endif

c
c     set maximum iterations in solve_new
c     
      iter=maxsolve
      call storage_derivatives(1,1)
      allocate(dumz(neq*3*4))
c
c     set maximum recoupling iterations if necesary
c     
      if(icoupl.gt.0) then
        icoupl_iter = icoupl
      else if(icoupl.lt.0) then
        icoupl_iter =
     &   max(float(abs(icoupl)),
     &   (float(iad)/float(abs(maxit)))*float(abs(2*icoupl)))
      endif
      if(irdof.eq.0) then
c
c     full dof solution
c

             if(gdpm_flag.eq.0) then
              if (igauss .gt. 1) then
               call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop
     *              ,north,tollr,irb,iirb,npvt,gmres,dumz,piv
     *              ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor,accm)
              else
               call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     *              ,north,tollr,irb,iirb,npvt,gmres,dumz,piv
     *              ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor,accm)
              end if 

             else
              if (igauss .gt. 1) then
                  call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &                 ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &                 ,npvt,gmres,dumz,piv
     &                 ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor
     &                 ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                 ,mdof_sol)
               else
                  call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat,nrhs
     &                 ,nelm,nelm_primary,nelm_primary,north,tollr,irb
     &                 ,iirb,npvt,gmres,dumz,piv
     &                 ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor
     &                 ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &                 ,mdof_sol)
               end if
             endif

      else if(irdof.eq.-idofm) then
c     reduced degree of freedom with full grmres
         call solve_rdof(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     &        ,tollr,irb,iirb,npvt,gmres,dumz,piv
     &        ,h,c,ss,g,y,iter,iback,-idofm,iptty
     &        ,maxor,icoupl_iter,tollr,overf,accm)
      else
c     
c     approximate solution
c     
         allocate(sto5(neq*6))
         call rdof_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     &        ,tollr,irb,iirb,npvt,gmres,dumz,piv
     &        ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor
     &        ,overf,irdof,icoupl_iter,0,sto5,accm)
         deallocate(sto5)
      endif
      itert=itert+iter 
      itotals=itotals+iter
      minkt=minkt+mink
      deallocate(dumz)
      call storage_derivatives(0,1)

 999  continue
      
      if(islord.ne.0) call switch(nmat,nmatb,islord,idofm,2,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,idofm,2,neq)

c   enforce thermal equilibrium( set t2 = t1)

      call methh2o_combine(2,idofm)

      return
      end
