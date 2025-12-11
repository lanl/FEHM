      subroutine explicit_phase_update(iflg,ij)
c calculated new variables values for potential phase_change gridblocks
c can calculate individual gridblocks with NR scheme.
c gaz 122220
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
      implicit none
      integer iflg,i,j,i1,i2,ij,jmia,neqp1,iter_max
      integer indx(3)
      real*8 dterm,pvapor,dpsatt,dpsats,dtsatp,psatl
      
      if(i_ex_update.le.0) return
c 
      if(iflg.eq.-1) then
c read input data
c i_ex_update-identifies explicit update
c iter_expa - maxiter for explicit update          
       read(inpt,'(a80)') wdd1
       read(wdd1,*) i_ex_update,iter_expa,r123tol
      else if(iflg.eq.0) then
c initialize variables and allocate memory 
c this iflag = 0 called from varchk.f          
       if(.not.allocated(a_ex))then   
        allocate(a_ex(3,3))
        allocate(r_ex(3))
        allocate(r1(neq))
        allocate(r2(neq))
        allocate(r3(neq))
        allocate(nphase_chk(neq))
        a_ex = 0.0d0
        r_ex = 0.0d0
        r1 = 0.0d0
        r2 = 0.0d0
        r3 = 0.0d0
       endif
        nphase_chk = 0
        iphase_chk = 0
      else if(iflg.eq.1) then
c keep track of potential phase change gridblocks 
c this iflag = 1 called from varchk.f
       iphase_chk = iphase_chk+1
       nphase_chk(iphase_chk) = ij 
      else if(iflg.eq.2) then
c call thrmwc.f for fluid properties
c set iteration NR loop
       if(iphase_chk.le.0) return
c establish flow terms without (remove accumulation terms)
c with updated variables
c call equation assembly
       allocate(r1_wo_acc(n))
       allocate(r2_wo_acc(n))
       allocate(r3_wo_acc(n))
       do i = 1, neq
        call geneqc(i)
       enddo
c 
       do i = 1, neq
        r1_wo_acc(i) = bp(i+nrhs(1))-sx1(i)*deni(i)
        r2_wo_acc(i) = bp(i+nrhs(2))-sx1(i)*denei(i)
        r3_wo_acc(i) = bp(i+nrhs(3))-sx1(i)*denpci(i)
       enddo  
c  gaz test 12221
       neq_phase_chk = neq
       do i = 1, neq
        nphase_chk(i) = i
       enddo
c       
       do iter = 1, iter_expa
       r1min = -1.e10
       r2min = -1.e10
       r3min = -1.e10
       do i = 1, neq_phase_chk
        ieq_ex = nphase_chk(i)
c thrmwc needs small mod for loop "mi =1, n"        
        call thrmwc(0)
       enddo
c assemble residual and solve for updates  
       r1 = 0.0
       r2 = 0.0
       r3 = 0.0
       neqp1=neq+1
c possible already done in call to geneqc above       
       do i = 1, neq_phase_chk
         ieq_ex = nphase_chk(i)
c
          r_ex(1) = sx1(ieq_ex)*deni(ieq_ex)+r1_wo_acc(i)
          r_ex(2) = sx1(ieq_ex)*denei(ieq_ex)+r2_wo_acc(i)
          r_ex(3) = sx1(ieq_ex)*denpci(ieq_ex)+r3_wo_acc(i)
          a_ex(1,1) = sx1(ieq_ex)*dmpf(ieq_ex)
	    a_ex(1,2) = sx1(ieq_ex)*dmef(ieq_ex)
	    a_ex(1,3) = sx1(ieq_ex)*dmc(ieq_ex)
	    a_ex(2,1) = sx1(ieq_ex)*depf(ieq_ex)   
	    a_ex(2,2) = sx1(ieq_ex)*deef(ieq_ex)   
          a_ex(2,3) = sx1(ieq_ex)*dec(ieq_ex)
  	    a_ex(3,1) = sx1(ieq_ex)*dcp(ieq_ex)
  	    a_ex(3,2) = sx1(ieq_ex)*dce(ieq_ex)
          a_ex(3,3) = sx1(ieq_ex)*dcc(ieq_ex)         
c
c call linear solver 
c  
c  Ax = b, test - solution must be "1"
c        r_ex(1) = a_ex(1,1)+a_ex(1,2)+a_ex(1,3)
c        r_ex(2) = a_ex(2,1)+a_ex(2,2)+a_ex(2,3)
c        r_ex(3) = a_ex(3,1)+a_ex(3,2)+a_ex(3,3)
c        
          call ludcmp0(a_ex,3,3,indx,dterm)
          call lubksb0(a_ex,3,3,indx,r_ex)
          r1(ieq_ex) = r_ex(1)
          r2(ieq_ex) = r_ex(2)
          r3(ieq_ex) = r_ex(3)
c apply correction
          if(ieos(ieq_ex).ne.2) then
           phi(ieq_ex)=phi(ieq_ex)-r1(ieq_ex)*strd
	     t(ieq_ex)=t(ieq_ex)-r2(ieq_ex)*strd
           pci(ieq_ex)=pci(ieq_ex)-r3(ieq_ex)*strd
          elseif(ieos(ieq_ex).eq.2) then
           phi(ieq_ex)=phi(ieq_ex)-r1(ieq_ex)*strd
           s(ieq_ex)=s(ieq_ex)-r2(ieq_ex)*strd
           t(ieq_ex)=t(ieq_ex)-r3(ieq_ex)*strd
                      pvapor=psatl(t(i),pcp(i),dpcef(i),dpsatt,
     &                dpsats,0,an(i))                        
                     if(pvapor.ge.phi(i)) then
                      t(i)=psatl(phi(i),pcp(i),dpcef(i),
     &                          dtsatp,dpsats,1,an(i)) 
                        pvapor = phi(i)
                     endif     
                     pci(i) = phi(i)-pvapor
          endif
        r1min = max(abs(r1(ieq_ex)),r1min)
        r2min = max(abs(r2(ieq_ex)),r2min)
        r3min = max(abs(r3(ieq_ex)),r3min)
        r123min = max(r1min,r2min,r3min)
       enddo
       if(r123min.le.r123tol) then
        go to 100
       endif
       enddo
100    continue
      deallocate(r1_wo_acc,r2_wo_acc,r3_wo_acc)
      else if(iflg.eq.3) then
c release memory
       deallocate(a_ex,r_ex,r1,r2,r3,nphase_chk)
      else if(iflg.eq.41) then
c gaz test 012221
         if(l.ge.10) iphase_chk = neq
         neq_phase_chk = neq
         do i = 1, neq
          nphase_chk(i) = i
         enddo
        if(iphase_chk.le.0) return
c write out variables before corrections
        write(ierr,*) 'l =',l,' iad ',iad
        write(ierr,*) 'phase change variables before ex corr'
        neq_phase_chk = iphase_chk
        do i = 1, neq_phase_chk
         ieq_ex = nphase_chk(i)
         if(ieos(ieq_ex).ne.2) then
          write (ierr,101) ieq_ex,ieos(ieq_ex),phi(ieq_ex),
     &      t(ieq_ex),pci(ieq_ex)
         else if(ieos(ieq_ex).eq.2) then
          write (ierr,102) ieq_ex,ieos(ieq_ex),phi(ieq_ex),
     &      s(ieq_ex),t(ieq_ex)
         endif
       enddo
      else if(iflg.eq.42) then
c write out variables after  corrections
        if(iphase_chk.le.0) return  
        write(ierr,*) 'l =',l,' iad ',iad
        write(ierr,*) 'phase change variables after ex corr'
        do i = 1, neq_phase_chk
         ieq_ex = nphase_chk(i)
         if(ieos(ieq_ex).ne.2) then
          write (ierr,101) ieq_ex,ieos(ieq_ex),phi(ieq_ex),
     &      t(ieq_ex),pci(ieq_ex)
         else if(ieos(ieq_ex).eq.2) then
          write (ierr,102) ieq_ex,ieos(ieq_ex),phi(ieq_ex),
     &      s(ieq_ex),t(ieq_ex)
         endif
       enddo
      endif
101   format(1x,'node ',i6,' ieos ',i3,1p,
     &   ' phi ',g12.4,' t ',g12.4,' pci ',g12.4)    
102   format(1x,'node ',i6,' ieos ',i3,1p,
     &   ' phi ',g12.4,' s ',g12.4,' t ',g12.4)          
      return
      end