      subroutine varchk_AWH(iflg)
c new approach to phase change
c gaz 022221
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use davidi
c new use file      
      use com_exphase
      use commass_AWH
      implicit none
c      
      integer iflg,i,j,i1,i2,ij,jmia,neqp1
      real*8 pl,tl,pcl
      real*8 pvapor,dpsatt,dpsats,dtsatp,psatl
      real*8 wamass,wmass_l,wmass_v
      
c gaz 040122 removed imass_phase = 1 now set in input
c macro 'mcpc'
c      r123tol = 1.e-6
      if(imass_phase.eq.0) return
c 
      if(iflg.eq.0) then
c read input data
c ivar_mass-identifies mass-based variable
c          
       read(inpt,'(a80)') wdd1
       read(wdd1,*) imass_phase, ivar_mass, r123tol
c gaz 052322 iamss_phase = -1, then change ngas input       
       if(imass_phase.eq.-1) then
       endif
c
      else if(iflg.eq.1) then
c gaz 080622 manage mass conservative phase change
c outline of algorithm
c generate set of new balanace eqs after NR update
c form total mass fraction and energy fraction
       call generate_balance_equations_ngas(0,0,0,0)
       call generate_balance_equations_ngas(1,1,neq,0)
       call generate_balance_equations_ngas(2,1,neq,0)
       call generate_balance_equations_ngas(3,1,neq,0)
       
      else if(iflg.eq.2) then
c       
c establish flow terms without (remove accumulation terms)
c with updated variables
c call equation assembly
      if(.not.allocated(r1_mass)) then
        allocate(r1_mass(n))
        allocate(r2_energy(n))
        allocate(r3_amass(n))
        allocate(wmass_awh(n))
        allocate(energy_awh(n))   
        allocate(amass_awh(n))
        allocate(nphase_chk(n))
      endif  
c get new equation parameter properties  with updated variables   
c quick check on thermal phase 
      do i = 1,n 
       tl = t(i)
       pvapor = psatl(tl,0.,0.,dpsatt,dpsats,0,0.0)
       pl = phi(i)
       pcl = pci(i) 
c need to calculate vapor equivalent of pci for henry's law
       if(pl.le.pvapor.and.ieos(i).eq.1) then
        phi(i) = pvapor
        ieos(i) = 2
       else if(pl-pcl.ge.pvapor.and.ieos(i).eq.3) then
        phi(i) = pvapor+pcl
        ieos(i) = 2
       endif
      enddo    
       call thrmwc(0)      
c form residuals       
       do i = 1, n
        call geneqc(i)
       enddo
c save residual and mass and energy quantities
       do i = 1, n
c save conservation EQs residuals       
        r1_mass(i)    = bp(i+nrhs(1))
        r2_energy(i)  = bp(i+nrhs(2))
        r3_amass(i)   = bp(i+nrhs(3))
c gridblock water mass (kg)   
        wmass_awh(i)  = sx1(i)*(deni(i)*dtot + denh(i))
c gridblock energy (Mj)        
	  energy_awh(i) = sx1(i)*(denei(i)*dtot + deneh(i))
c gridblock air mass(kg)	
        amass_awh(i)  = sx1(i)*(denpci(i)*dtot + denpch(i))
       enddo  
c 
c could check here for P,Psat(T) change. 
c could also incorporate equation residual (kg/sec or Mw) into volume check
c check for possible phase change for equations with high residuals
c
       n_phase_chk = 0
       do i = 1, n
        nphase_chk(i) = 0
        r1min = abs(r1_mass(i))
        r2min = abs(r2_energy(i))
        r3min = abs(r3_amass(i))        
        r123min = max(r1min,r2min,r3min)
c count residuals that are higher than r123tol       
        if(r123min.ge.r123tol) then
         n_phase_chk = n_phase_chk + 1
         if(r1min.eq.r123min) then
          nphase_chk(i) = 1
         else if(r2min.eq.r123min) then
          nphase_chk(i) = 2
         else
          nphase_chk(i) = 3
         endif
        endif
       enddo  

c write out variables before phase check adjustments
        write(ierr,*) '***********************************************'
        write(ierr,*)
        write(ierr,*) 'l =',l,' iad ',iad,
     &  ' grid and variable info after NR update'  
        do i = 1, n         
         if(ieos(i).ne.2) then
          write (ierr,101) ieos(i),nphase_chk(i),phi(i),
     &      t(i),pci(i),r1_mass(i),r2_energy(i),r3_amass(i)
         else if(ieos(i).eq.2) then
          write (ierr,102) i,ieos(i),phi(i),
     &      s(i),t(i),r1_mass(i),r2_energy(i),r3_amass(i)
         endif
        enddo       
c       
      do i = 1,n
       if(nphase_chk(i).ne.0) then
c check for phase change
c wmass_l is the total liquid phase mass (wmass + amass)
c wmass_v is the total vapor phase mass (wmass + amass)
c wamass is total accumulation term water mass (h2o + and ngas component)
       	wmass_l = sx1(i)*ps(i)*rolf(i)
       	wmass_v = sx1(i)*ps(i)*rovf(i)
       	wamass = wmass_awh(i) + amass_awh(i)
       	if(wamass.ge.wmass_l) then
c liquid phase check       	
       	 ieos(i) = 1
       	 strd =1.
       	 s(i) =1.
       	else if(wamass.le.wmass_v) then
c vapor phase check       	
       	 ieos(i) = 3
       	 strd = 1.
       	 s(i) = 0.
       	else
c two phase
          ieos(i) = 2
          strd = 1.
c calculate saturation  
          s(i) = (wamass-wmass_v)/(wmass_l-wmass_v)
         endif
       endif
      enddo
c
c recheck new equation parameter properties 
c
       do i = 1, n
        bp(i+nrhs(1)) = 0.0
        bp(i+nrhs(2)) = 0.0
        bp(i+nrhs(3)) = 0.0 
       enddo
c P,T,PC the same, phase state changes and S new estimate
       call thrmwc(0)      
c form residuals       
       do i = 1, n
        call geneqc(i)
       enddo
       
c save residual and mass and energy quantities
       do i = 1, n
        r1_mass(i)    = bp(i+nrhs(1))
        r2_energy(i)  = bp(i+nrhs(2))
        r3_amass(i)   = bp(i+nrhs(3))
        wmass_awh(i)  = sx1(i)*(deni(i)*dtot + denh(i))
	  energy_awh(i) = sx1(i)*(denei(i)*dtot + deneh(i))
        amass_awh(i)  = sx1(i)*(denpci(i)*dtot + denpch(i))
       enddo  
c  gaz test 12221
       n_phase_chk = 0
       do i = 1, n
        nphase_chk(i) = 0
        r1min = abs(r1_mass(i))
        r2min = abs(r2_energy(i))
        r3min = abs(r3_amass(i))        
        r123min = max(r1min,r2min,r3min)
        if(r123min.ge.r123tol) then
         n_phase_chk = n_phase_chk + 1   
        endif 
        if(r1min.eq.r123min) then
         nphase_chk(i) = 1
        else if(r2min.eq.r123min) then
         nphase_chk(i) = 2
        else
         nphase_chk(i) = 3
        endif       
       enddo

c save residual and mass and energy quantities
       do i = 1, n
        r1_mass(i)    = bp(i+nrhs(1))
        r2_energy(i)  = bp(i+nrhs(2))
        r3_amass(i)   = bp(i+nrhs(3))
        wmass_awh(i)  = sx1(i)*(deni(i)*dtot + denh(i))
	  energy_awh(i) = sx1(i)*(denei(i)*dtot + deneh(i))
        amass_awh(i)  = sx1(i)*(denpci(i)*dtot + denpch(i))
       enddo         
       
       
c write out variables and residuals after phase check adjustments
        write(ierr,*) '***********************************************'
        write(ierr,*)
        write(ierr,*) 'l =',l,' iad ',iad,
     &  ' grid and variable info after NR update'  
        do i = 1, n         
         if(ieos(i).ne.2) then
          write (ierr,101) 1,ieos(i),nphase_chk(i),phi(i),
     &      t(i),pci(i),r1_mass(i),r2_energy(i),r3_amass(i)
         else if(ieos(i).eq.2) then
          write (ierr,102) i,ieos(i),phi(i),
     &      s(i),t(i),r1_mass(i),r2_energy(i),r3_amass(i)
         endif
        enddo  
      endif
101   format(1x,'node ',i6,' ieos ',i3,' phase_chk ',i3, 1p,
     &   ' phi ',g12.4,' t ',g12.4,' pci ',g12.4,
     &   ' r1 ',g12.4, ' r2 ',g12.4, ' r3 ',g12.4)    
102   format(1x,'node ',i6,' ieos ',i3,' phase_chk ',i3, 1p,
     &   ' phi ',g12.4,' s ',g12.4,' t ',g12.4,
     &   ' r1 ',g12.4, ' r2 ',g12.4, ' r3 ',g12.4)        
      return
      end