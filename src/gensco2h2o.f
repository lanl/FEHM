      subroutine gensco2h2o    
C***********************************************************************
CD1
CD1  PURPOSE 
CD1
CD1  This subroutine solves the non-isothermal air-water-co2 equations 
CD1  full jacobian (unsymmetric, 3n by 3n).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 
CD2    Rev 1.0   02/05/07 First version RJP
CD2 
C************************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  RJP 02/05/07 First version of gensco2h2o subroutine
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
      use comco2
      use comwellphys
      implicit none

      integer nmatd, index(6)
      integer i, j, icoupl_iter, id, idl, idofm 
      integer jj,ii,idiag
      integer iprint,inorm
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls
      real*8 tmch1
      parameter (tmch1=1d-3)
      parameter (iprint=0,inorm = 1)
      
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
c    call to ther_co2_h2o to get rel perms & thermodynamic properties
c    for water-rich phase and co2-rich phase. The thermodynamic properties 
c    are now stored in arrays and are accessed during formation of eqns.
c
        call ther_co2_h2o(0,0)
c
c     water component: add heat-mass contribution of water
c     Conservation of water-mass equation
c     Conservation of water-energy equation
c     
	if(idof_co2.eq.2) goto 10
         call ther_co2_h2o(1,0)
c     
c     RJP 6/29/04 introduced new arrays to carry dil for 'trac'
c     calculations.  It is assumed that tracer/reactants are
c     carried in liquid (water).
c
      diltrac=dil
c RJP 02/05/07
c Generate water specific equation
c follwoing is modified form of original geneq2_h2o subroutine for 
c CO2-water-air problem. 
c
c
c GAZ 11/01/08 get drift flux info (must be call after rel perms)
      call wellphysicsctr(2,0)
c
      do id = 1,neq
          call geneq_h2o_co2(id)
      enddo
	if(idof_co2.lt.2) goto 20
 10	continue
c  
c     CO2 component: add heat-mass contribution of co2            
c     Conservation of CO2-mass equation
c     Conservation of CO2-energy equation
c
         call ther_co2_h2o(2,0)
c     
c     Generate CO2 specific equation
c     
      do id = 1,neq
         if(nwellphy.ne.0) then
          call geneq_co2_wellphysics(id)
         else
          call geneq_co2(id)
         endif
      enddo
c  
c     air component: add heat-mass contribution of air            
c     Conservation of air-mass equation
c     Conservation of air-energy equation
c
	if(idof_co2.le.4) goto 20
c     
c     Generate air specific equation
c     

         call ther_co2_h2o(3,0)
c     
c     Generate air specific equation
c     
      do id = 1,neq
          call geneq_air_co2(id)
      enddo
 
 20	continue
c
c    call to ther_co2_h2o to release memory for rel perms
c
        call ther_co2_h2o(-1,0)
c
c   enforce equilibrium conditions as necessary
c
      call co2h2o_combine(1,idofm)

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

c RJP 04/30/05 added following for zero porosity nodes
      do id=1,neq
         if(ps(id).le.0.0) then
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0
		  if(iprtype.eq.-3) then
            a(nelmdg(id)-neqp1+nmat(4))=sx1(id)
            bp(id+nrhs(2))=0.0
		  else
            a(nelmdg(id)-neqp1+nmat(9))=sx1(id)
            bp(id+nrhs(3))=0.0
		  endif
         endif
      enddo

      if(islord.ne.0) call switch(nmat,nmatb,islord,idofm,1,nsizea1)
      if(islord.ne.0) call switchb(bp,nrhs,nrhsb,islord,idofm,1,neq)

c
c     normalize the equations
c     
      allocate(dumn(108))
c      write (ierr,*) 'l = ',l,' iad = ',iad
c      do id = 1,neq
c       if(kaco2(id).ne.0) then
c       jj = nelmdg(id)-neqp1
c       write (ierr,*) 'node = ',id
c       write (ierr,*) 'phase state ', ices(id)
c       write (ierr,*) 'P ',phi(id),' T ', t(id)
c      write (ierr,*) 'dq(id)',dq(id)
c       write (ierr,*) a(jj+nmat(1)),a(jj+nmat(4)),a(jj+nmat(9))
c       write (ierr,*) bp(id+nrhs(1)),bp(id+nrhs(2)),bp(id+nrhs(3))
c       endif
c       if(abs(bp(id+nrhs(2))).gt.1.e5) then
c        write (ierr,*) 'bad node ',id, bp(id+nrhs(2))
c       endif
c      enddo      
c  gaz 2-27-03
       if(inorm.ne.0) then
        call normal_dof(neq,a,bp,nelm,nmat,nrhs,nelmdg
     &      ,index,idofm,dumn(1),dumn(17),dumn(33),0,fdum2)
       else
	 fdum2 = 0.0
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
         write(iout,*) '   '
         write(iout,*) '* singular matrix found during normalization *'
         write(iout,*) '   '
c         write(iout,'(i8,3g12.3)')
c     &        abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         if(iptty.gt.0)then
         write(iptty,*) '   '
         write(iptty,*) '* singular matrix found during normalization *'
         write(iptty,*) '   '
c         write(iptty,'(i8,3g12.3)') 
c     &        abs(nmatb(1)),(cord(abs(nmatb(1)),i),i=1,3)
         endif
         iad=maxit
         return
      endif
      
c     find residual
c     
      fdum=sqrt(fdum2)
      if(fdum.eq.0.0) go to 999
      mink=n
      if(iad.eq.0) then
         f0=max(fdum*epe,tmch)
      endif
      if(fdum1.lt.0.0.and.iad.ne.0) then
         do i=1,neq
          do j=1,idofm
c            if (idco2(i).eq.1) then
c               if(abs(bp(i+nrhs(j))).gt.tmch1) go to 99
c            else
               if(abs(bp(i+nrhs(j))).gt.tmch) go to 99
c            endif
          enddo
         enddo
         fdum=-1.0
         go to 999
 99      continue
         f0=-1.0      
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
      allocate(dumz(neq*5*4))
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
               call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     &              ,tollr,irb,iirb,npvt,gmres,dumz,piv
     &              ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor,accm)
            else
               call solve_new(neq,a,b,bp,nmat,nmat,nrhs,nelm,nelm
     &              ,north,tollr,irb,iirb,npvt,gmres,dumz,piv
     &              ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
               
            end if
         else
            if (igauss .gt. 1) then
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nb,nrhs
     &              ,nelm,nelm_primary,nop,north,tollr,irb,iirb
     &              ,npvt,gmres,dumz,piv
     &              ,h,c,ss,g,y,iter,iback,idofm,iptty,maxor
     &              ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &              ,mdof_sol)
            else
               call solve_dual(neq_primary,neq,a,b,bp,nmat,nmat,nrhs
     &              ,nelm,nelm_primary,nelm_primary,north,tollr,irb
     &              ,iirb,npvt,gmres,dumz,piv
     &              ,h,c,ss,g,y,iter,iback,1,iptty,maxor
     &              ,igdpm,maxgdpmlayers,ngdpm_layers,nelmdg,accm
     &              ,mdof_sol)
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

      call co2h2o_combine(2,idofm)

      return
      end
