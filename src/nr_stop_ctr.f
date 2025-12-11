       subroutine nr_stop_ctr(iflg)
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
!!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 stop newton-raphson iteration on variable change (not equation tol)
!D1 additional stopping criteria added beginning about 090119 
!D1 gaz 122121 added  subroutine phase_count
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 29-Nov-09, Programmer: George Zyvoloski
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 3.0
!D4 
!**********************************************************************

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai
      use comco2

      implicit none

      integer, allocatable :: kq_dum(:), icount(:)
      integer i,j,ii,ij,jj,kb,i1,i2,neqp1,max_var,iflg
      integer inr_chkp,inr_chks,inr_chks2,inr_chkt,nr_write
      character*8 vard
      logical nr_chk
      real*8 value,pnrmax,snrmax,snr2max,tnrmax,hnrmax,panrmax
      real*8 tol_conserve
      parameter (max_var = 20, tol_conserve = 1.d-6)
      save nr_write
c gaz 080220 debug help next line no effect otherwise
      i = fdum + l + iad
c       
      if(nr_stop.eq.0) return
      if(iflg.eq.0) then
c
c input
c
      nr_write = 0
      backspace inpt
      read(inpt,'(a80)') wdd1
      do i = 5,20
       if(wdd1(i:i).eq.'c'.or.wdd1(i:i).eq.'C') nr_write = 1
      enddo
      nr_chk = .true.
      do i = 1,max_var
       read(inpt,'(a4)') vard
       if(vard.eq.'    ') then
        go to 100
       else
        backspace inpt
        read(inpt,'(a80)') wdd1
        read(wdd1,*) vard
        if(vard(1:2).eq.'pw'.and.vard(3:3).eq.' ') then
         read(wdd1,*) vard(1:2), p_stop
        endif
        if(vard(1:2).eq.'pa'.and.vard(3:3).eq.' ') then
         read(wdd1,*) vard(1:2),pa_stop
        endif        
        if(vard(1:1).eq.'t'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),t_stop
        endif
c gaz 070821 p_tol and t_tol control uniqueness of property calls        
        if(vard(1:1).eq.'p'.and.vard(2:2).eq.'_') then
         read(wdd1,*) vard(1:1),p_tol
        endif 
        if(vard(1:1).eq.'t'.and.vard(2:2).eq.'_') then
         read(wdd1,*) vard(1:1),t_tol
        endif
        if(vard(1:1).eq.'c'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),co2f_stop
        endif 
        if(vard(1:1).eq.'h'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),h_stop
        endif
        if(vard(1:4).eq.'step') then
         read(wdd1,*) vard(1:4), strd_iter
        endif
        if(vard(1:4).eq.'mste') then
         read(wdd1,*) vard(1:4), strd_mass  
        endif
c gaz 111520 added option for strd to be reset to 1 at every iteration (not just at iad eq 0)        
        if(vard(1:4).eq.'intv') then
         read(wdd1,*) vard(1:4), iter_intrvl  
        endif        
        if(vard(1:1).eq.'s'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),s_stop
        else if(vard(1:3).eq.'sat'.and.vard(4:4).eq.' ') then
         read(wdd1,*) vard(1:3),s_stop
        else if(vard(1:1).eq.'s'.and.vard(2:2).eq.'2') then
         read(wdd1,*) vard(1:2),s2_stop  
c gaz 032622 add stopping based on conservation eqs        
        else if(vard(1:4).eq.'eqwm'.and.vard(5:5).eq.' ') then
         read(wdd1,*) vard(1:5),eqwm_stop  
        else if(vard(1:4).eq.'eqen'.and.vard(5:5).eq.' ') then
         read(wdd1,*) vard(1:5),eqen_stop      
        else if(vard(1:4).eq.'eqnm'.and.vard(5:5).eq.' ') then
         read(wdd1,*) vard(1:5),eqnm_stop            
        endif
       endif
      enddo
100   continue  
       var_stop = .true.
       if(p_stop.eq.0.0.and.s_stop.eq.0.0.and.s2_stop.eq.0.0.
     &    and.t_stop.eq.0.0.and.h_stop.eq.0.0.and.pa_stop.eq.0.0) then
        var_stop = .false.
      endif
       con_eq_stop = .true.
       write(ierr,*) 'eqwm_stop ', eqwm_stop,' eqen_stop ',eqen_stop,
     &               ' eqnm_stop ', eqnm_stop   
       if(eqwm_stop.eq.0.0.and.eqen_stop.eq.0.0.and.           
     &  eqnm_stop.eq.0.0) then
        con_eq_stop = .false.
       endif
c check for inconsistency between stopping on variable change or conservation equations     
       if(var_stop.and.con_eq_stop) then
         write(ierr,*) 'NR iterations: cannot stop iterations using',
     &   ' variable chk and conservation eq balance. Stopping.'
         if(iout.ne.0) then
          write(iout,*) 'NR iterations: cannot stop iterations using',
     &    ' variable chk and conservation eq balance. Stopping.'
         endif
         if(iptty.ne.0) then
          write(iptty,*) 'NR iterations: cannot stop iterations using',
     &    ' variable chk and conservation eq balance. Stopping.'
         endif         
        stop 
       endif
c        
c check for sufficient closure terms
c
       if(var_stop) then
       if(icarb.lt.0) then
        if(p_stop.le.0.0) nr_chk = .false.
        if(t_stop.le.0.0) nr_chk = .false.
        if(s_stop.le.0.0) nr_chk = .false.     
        if(s2_stop.le.0.0) nr_chk = .false.                     
       else if(ico2.lt.0) then
        if(p_stop.le.0.0) nr_chk = .false.
       else if(ico2.eq.0) then
        if(p_stop.le.0.0) nr_chk = .false.
        if(t_stop.le.0.0) nr_chk = .false.  
       else if(ico2.gt.0) then
        if(p_stop.le.0.0) nr_chk = .false.
        if(t_stop.le.0.0) nr_chk = .false.
        if(pa_stop.le.0.0) nr_chk = .false.      
        if(s_stop.le.0.0) nr_chk = .false.   
       endif
       else if(con_eq_stop) then
          if(eqwm_stop.le.0) nr_chk = .false.
          if(eqen_stop.le.0) nr_chk = .false.
          if(eqnm_stop.le.0) nr_chk = .false.
       endif
       if(.not.nr_chk) then
        if(iptty.ne.0) then
         write(iptty,*)
         write(iptty,*) '>> insufficient NR closure (macro nrst)  <<'
         write(iptty,*)
        endif
        if(iatty.ne.0) then
         write(iatty,*)
         write(iatty,*) '>> insufficient NR closure (macro nrst)  <<'
         write(iatty,*)
        endif 
          write(ierr,*) '>> insufficient NR closure (macro nrst)  <<' 
            
       endif
            
      else if(iflg.eq.1) then
c 
c  iteration control based on variable change called from varchk.f
c  
       if(.not.var_stop) return  
       nr_stop = 2
       if(icarb.ne.0) then
         pnrmax = 0.0
         snrmax = 0.0
         snr2max = 0.0
         tnrmax = 0.0
         inr_chkp= 0
         inr_chks= 0
         inr_chks2= 0
         inr_chkt= 0      
        do i = 1,n
         if(p_stop.ne.0.0.and.abs(bp(i+nrhs(1))).gt.pnrmax) then
          pnrmax = abs(bp(i+nrhs(1)))
          inr_chkp= i
         endif
         if(ices(i).eq.2) then
           if(s_stop.ne.0.0.and.abs(bp(i+nrhs(2))).gt.snrmax) then
            snrmax = abs(bp(i+nrhs(2)))
            inr_chks= i
           endif
           if(s2_stop.ne.0.0.and.abs(bp(i+nrhs(3))).gt.snr2max) then
            snr2max = abs(bp(i+nrhs(3)))
            inr_chks2= i
          endif         
         else if(ices(i).ne.2) then      
          if(t_stop.ne.0.0.and.abs(bp(i+nrhs(2))).gt.tnrmax) then
           tnrmax = abs(bp(i+nrhs(2)))
           inr_chkt= i
          endif
          if(s_stop.ne.0.0.and.abs(bp(i+nrhs(3))).gt.snrmax) then
            snrmax = abs(bp(i+nrhs(3)))
            inr_chks= i
          endif
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(s_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3 
        if(t_stop.ne.0.0.and.tnrmax.gt.t_stop) nr_stop = 3  
        if(s2_stop.ne.0.0.and.snr2max.gt.s2_stop) nr_stop = 3 
        if(nr_write.ne.0) then
         write(ierr,50) 
     &    l,iad,inr_chkp,pnrmax,inr_chks,snrmax,inr_chkt,tnrmax       
        endif       
       else if(ico2.lt.0) then
         pnrmax = 0.0
         snrmax = 0.0
         hnrmax = 0.0
c gaz 080220 corrected s_stop, snrmax had nrhs(1)         
        do i = 1,n
         pnrmax = max(bp(i+nrhs(1)),pnrmax)
         if(ifree.eq.0.and.irdof.ne.13.and.ieos(i).eq.2) then
          if(s_stop.ne.0) then
           snrmax = max(bp(i+nrhs(2)),snrmax)
          endif
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(s_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3    
        if(h_stop.ne.0.0.and.pnrmax.gt.h_stop) nr_stop = 3      
       else if(ico2.eq.0) then
         pnrmax = 0.0
         snrmax = 0.0
         tnrmax = 0.0
         inr_chkp= 0
         inr_chks= 0
         inr_chkt= 0
        do i = 1,n
         if(p_stop.ne.0.0.and.abs(bp(i+nrhs(1))).gt.pnrmax) then
          pnrmax = abs(bp(i+nrhs(1)))
          inr_chkp= i
         endif
         if(s_stop.ne.0.0.and.abs(bp(i+nrhs(2))).gt.snrmax.and.
     &    ieos(i).eq.2) then
          snrmax = abs(bp(i+nrhs(2)))
          inr_chks= i
         endif
         if(t_stop.ne.0.0.and.abs(bp(i+nrhs(2))).gt.tnrmax.and.
     &    ieos(i).ne.2) then
          tnrmax = abs(bp(i+nrhs(2)))
          inr_chkt= i
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(s_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3 
        if(t_stop.ne.0.0.and.tnrmax.gt.t_stop) nr_stop = 3  
        if(nr_write.ne.0) then
         write(ierr,50) 
     &    l,iad,inr_chkp,pnrmax,inr_chks,snrmax,inr_chkt,tnrmax       
        endif 
50     format(1x,'TS',1x,i5,' NR',1x,i3,' i pres',1x,i7,1x,g12.5,
     &    ' i sat',1x,i7,1x,g12.5,' i t',1x,i7,1x,g12.5)         
       else if(ico2.gt.0) then
         pnrmax = 0.0
         snrmax = 0.0
         tnrmax = 0.0
         panrmax = 0.0
c gaz 032522 added correct 3rd variable for stopping  
c pci uses same stopping criteria as phi         
        do i = 1,n
         pnrmax = max(bp(i+nrhs(1)),pnrmax)
         if(ieos(i).eq.2) then
          snrmax = max(bp(i+nrhs(2)),snrmax)
          tnrmax = max(bp(i+nrhs(3)),tnrmax)
         else if(ieos(i).ne.2) then
          tnrmax = max(bp(i+nrhs(2)),tnrmax)
          panrmax = max(bp(i+nrhs(3)),panrmax)
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(s_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3 
        if(t_stop.ne.0.0.and.tnrmax.gt.t_stop) nr_stop = 3 
        if(pa_stop.ne.0.0.and.panrmax.gt.pa_stop) nr_stop = 3  
       endif             
      else if(iflg.eq.2) then      
c 
c  iteration control based on conservation equations  (descrpency of water mass, energy,etc)
c gaz 032622 first test with AWH, then with other physics
c 

      if(.not.con_eq_stop) return
      nr_stop = 4
      if(ico2.gt.0) then
       water_mass_err_max = 0.0d0
       energy_err_max = 0.0d0
       ngas_mass_err_max = 0.0d0
       do i = 1,neq
c         water_mass_err = bp(i+nrhs(1))/(denh(i)+tol_conserve)  
         water_mass_err = bp(i+nrhs(1))*dtot          
         if(abs(water_mass_err).gt.abs(water_mass_err_max)) then
          water_mass_err_max = water_mass_err
          node_water_err_max = i
         endif
c         energy_err = bp(i+nrhs(2))/(deneh(i)+tol_conserve) 
         energy_err = bp(i+nrhs(2))*dtot        
         if(abs(energy_err).gt.abs(energy_err_max)) then
          energy_err_max = energy_err
          node_energy_err_max = i
         endif
c gaz 040122 use absolute err for ngas         
c         ngas_mass_err = bp(i+nrhs(3))/(denpci(i)+tol_conserve) 
         ngas_mass_err = bp(i+nrhs(3))*dtot
          if(abs(ngas_mass_err).gt.abs(ngas_mass_err_max)) then
          ngas_mass_err_max = ngas_mass_err
          node_ngas_err_max = i
         endif
       enddo
       if(abs(water_mass_err_max).gt.abs(eqwm_stop)) nr_stop = 5
       if(abs(energy_err_max).gt.abs(eqen_stop)) nr_stop = 5
       if(abs(ngas_mass_err_max).gt.abs(eqnm_stop)) nr_stop = 5
      endif
      endif    
      return 
      end
      subroutine phase_count(iflg,is_ch,is_ch_t)
c
c subroutine to manage phase counting (gaz 122121, initial coding)
c  
      
      
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai
      use comco2
      use comwt,only : rlptol

      implicit none
      
      integer iflg, i
      integer is_ch
      integer is_ch_t
      real*8 sl, sl_old
cDEC$ FIXEDFORMLINESIZE:132      
      if(iflg.eq.1) then
               is_ch = 0
               nphase_liq = 0
               nphase_2 = 0
               nphase_gas = 0
               nphase_sc = 0
c gaz 121921 allow phase count for isothermal models             
              if(ico2.ge.0) then 
               do i=1,n
                if(ps(i).gt.0.0d0) then
                  if(ieos(i).eq.1) nphase_liq = nphase_liq + 1
                  if(ieos(i).eq.2) nphase_2 = nphase_2 + 1
                  if(ieos(i).eq.3) nphase_gas = nphase_gas + 1
                  if(ieos(i).eq.4) nphase_sc = nphase_sc + 1
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     if(s(i).lt.1.0.and.so(i).ge.1.0) then
                        is_ch=is_ch +1
                     else if(s(i).gt.0.0.and.so(i).le.0.0) then
                        is_ch=is_ch +1
                     else if(s(i).ge.1.0.and.so(i).lt.1.0) then
                        is_ch=is_ch +1
                     else if(s(i).le.0.0.and.so(i).gt.0.0) then
                        is_ch=is_ch +1
                     endif
                   endif                 
                  endif
               enddo
              else
c gaz 121921 allow phase count for isothermal models  (defined sl)                  
               do i=1,n
                if(ps(i).gt.0.0d0) then
                  if(ifree.ne.0) then
                   sl = rlxyf(i)
                  else if (irdof .ne. 13 ) then
                   sl = s(i) 
                  else
                   sl = 1.0d0
                  end if
                  if(sl.ge.1.0-rlptol) then
                    nphase_liq = nphase_liq + 1
                  else if(sl.lt.1.0-rlptol.and.sl.gt.rlptol) then
                    nphase_2 = nphase_2 + 1
                  else if(sl.le.rlptol) then
                    nphase_gas = nphase_gas + 1
                  endif
                   if(irdof.eq.13.and.ifree.ne.0) then
                    sl_old = so(i)
                   else
                    sl_old = 1.0
                   endif
                     if(sl.lt.1.0-rlptol.and.sl_old.ge.1.0-rlptol) then
                        is_ch=is_ch +1
                     else if(sl.gt.rlptol.and.sl_old.le.rlptol) then
                        is_ch=is_ch +1
                     else if(sl.ge.1.0-rlptol.and.sl_old.lt.1.0-rlptol) then
                        is_ch=is_ch +1
                     else if(sl.le.rlptol.and.sl_old.gt.rlptol) then
                        is_ch=is_ch +1
                     endif
c                   endif                 
                  endif
               enddo 
              endif 
c gaz 082223 calculate total changes differently  
c               is_ch_t = is_ch_t + is_ch
c          
               if(l.ne.0) then
                is_ch = 0
                dnphase_liq = nphase_liq - nphase_liq_0
                dnphase_2 = nphase_2 - nphase_2_0
                dnphase_gas = nphase_gas - nphase_gas_0
                dnphase_sc = nphase_sc - nphase_sc_0                                                                 
                if(dnphase_liq.gt.0) is_ch = is_ch + dnphase_liq
                if(dnphase_2.gt.0) is_ch = is_ch + dnphase_2
                if(dnphase_gas.gt.0) is_ch = is_ch + dnphase_gas
                if(dnphase_sc.gt.0) is_ch = is_ch + dnphase_sc
                is_ch_t = is_ch_t + is_ch
               else
                dnphase_liq = 0
                dnphase_2 = 0
                dnphase_gas = 0
                dnphase_sc = 0
               endif
               nphase_liq_0 = nphase_liq
               nphase_2_0 = nphase_2
               nphase_gas_0 = nphase_gas
               nphase_sc_0 = nphase_sc
      endif
      return
      end
          
