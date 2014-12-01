       subroutine active_nodes_ctr(iflg)
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
!D1 calculate active nodes and connections
!D1
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.20
!D2 
!D2 Initial implementation: Date 01-OCT-13, Programmer: George Zyvoloski
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

      integer, allocatable :: kq_dum(:)
      integer i,j,ii,ij,jj,kb,i1,i2,ic,max_var,iflg,n_activep1
      integer  ic1,neqp1,nmatd,ib,id
      integer open_file
      character*6 vard
      logical av_chk, nr_chk
      real*8 value,pnrmax,snrmax,snr2max,tnrmax,hnrmax,panrmax
      real*8 diffa, time_base, p_base, temp_base
      real*8 resid_i, resid_kb
      
      parameter (max_var = 10)

c
c enabled with keyword 'nact'      
c      
      if(iactive.eq.0) return
      if(iflg.eq.0) then
c
c input
c
      p_act_tol = 0.0d0 
      t_act_tol = 0.0d0
      bp_act_tol = 0.0d0
      ij_change = 0
      resid_tol = 1.e-5

      neqp1 = neq +1 
      av_write = 0
      backspace inpt
      read(inpt,'(a80)') wdd1
      do i = 5,20
       if(wdd1(i:i).eq.'c'.or.wdd1(i:i).eq.'C') av_write = 1
      enddo
      av_chk = .true.
      do i = 1,max_var
       read(inpt,'(a6)') vard
       if(vard.eq.'    ') then
        go to 100
       else
        backspace inpt
        read(inpt,'(a80)') wdd1
        read(wdd1,*) vard(1:6)
        if(vard(1:2).eq.'pw'.and.vard(3:3).eq.' ') then
         read(wdd1,*) vard(1:2), p_act_tol
        endif
        if(vard(1:2).eq.'pa'.and.vard(3:3).eq.' ') then
         read(wdd1,*) vard(1:2), p_act_tol
        endif        
        if(vard(1:1).eq.'t'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1), t_act_tol
        endif
        if(vard(1:1).eq.'c'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1), co2f_stop
        endif 
        if(vard(1:1).eq.'h'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1), h_act_tol
        endif
        if(vard(1:2).eq.'bp') then
         read(wdd1,*) vard(1:2), bp_act_tol
        endif        
        if(vard(1:1).eq.'s'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),s_stop, strd_iter 
        else if(vard(1:3).eq.'sat'.and.vard(4:4).eq.' ') then
         read(wdd1,*) vard(1:3),s_stop, strd_iter
        else if(vard(1:1).eq.'s'.and.vard(2:2).eq.'2') then
         read(wdd1,*) vard(1:2),s2_stop, strd_iter         
        endif
        if(vard(1:5).eq.'resid'.and.vard(6:6).eq.' ') then
         read(wdd1,*) vard(1:5), resid_tol
         ij_change = 1
        endif
       endif
      enddo
100   continue      

c
c find unit number for active node output (contour)
c
       iacfile = open_file('active_contour_file.txt', 'unknown')
     
      else if(iflg.eq.-1) then
c
c   initialize and allocate arrays
c
       allocate (node_active(n))
       allocate (phi_base(n))  
       allocate (time_phi_base(n))  
       allocate (phi_base_grad(n))  
       allocate (t_base(n))  
       allocate (time_t_base(n))  
       allocate (t_base_grad(n)) 
       node_active = 0
       phi_base = 0.0
       time_phi_base = 0.0
       phi_base_grad = 0.0
       t_base = 0.0
       time_t_base = 0.0
       t_base_grad = 0.0
      if(p_act_tol.ne.0.0d0) then
       phi_base = phi
      endif
      if(t_act_tol.ne.0.0d0) then
       t_base = t
      endif
      else if(iflg.eq.1) then
c 
c  check for varaible changes, residuals, and accumulation terms
c  best at iad = 0 in equation generation routine 
c  101713
c   
        iav_chkp = 0 
        iav_chkt = 0 
        iav_chkbp = 0 
c
        node_active = 0
c   
        do i = 1,n
c check main variable pw
         if(p_act_tol.ne.0.0d0) then
          p_base = phi_base(i)
          time_base = time_phi_base(i)-day_tol
          diffa = p_base-phi(i)
          if(abs(diffa).gt.p_act_tol) then
           iav_chkp = iav_chkp + 1
           time_phi_base(i) = days 
           phi_base_grad(i) = diffa/(days-time_base)
           node_active(i) = node_active(i) + 1
          endif
         endif
c check main variable t
         if(t_act_tol.ne.0.0d0) then
          temp_base = t_base(i)
          time_base = time_t_base(i)-day_tol
          diffa = temp_base-t(i)
          if(abs(diffa).gt.t_act_tol) then
           iav_chkt = iav_chkt + 1
           time_t_base(i) = days 
           t_base_grad(i) = diffa/(days-time_base)
           node_active(i) = node_active(i) + 1
          endif
         endif
         if(bp_act_tol.ne.0.0d0) then
c  check residual (should be not normalized- in kg/sec)
          if(abs(bp_base(i)-bp(i)).gt.bp_act_tol) then
           iav_chkbp = iav_chkbp + 1
           node_active(i) = node_active(i) + 1
           bp_base(i)=bp(i)
          endif
         endif
        enddo
      else if(iflg.eq.2) then
c 
c  renumber variables
c  
       ic = 0  
       do i = 1, n
        if(node_active(i).ne.0) then
         ic = ic + 1
         irb(ic) = i
         iirb(i) = ic
        endif
       enddo
       n_active = ic
      else if(iflg.eq.3) then
c 
c  establish time series for inactive nodes
c  this can be used to calculate apparent specified flow nodes (inactive nodes)
c  
       do i = 1, n
        if(node_active(i).eq.0) then
         p_base = phi_base(i)
         phi(i) = p_base + phi_base_grad(i)*(days-time_phi_base(i))
        endif
       enddo
      else if(iflg.eq.4) then
c
c  output 
c best called near contour output
c
       write(iacfile,*) days, iav_chkp, iav_chkt, iav_chkbp
c
       do i = 1,n
        write(89,200) i, cord(i,1), cord(i,2), cord(i,3),
     &       node_active(i), phi_base_grad(i), t_base_grad(i)
       enddo
200    format(1x, i5, 1p,6(1x,g12.3))

      else if(iflg.eq.5) then
      if(ij_change.eq.0) then
c identify in solver routines (igauss = -1)
       igauss = -1
c 
c  generate new nelm and load a and bp based on variable gradients
c  
       neqp1= neq+1
       nmatd=nelm(neqp1)
       if(.not.allocated(nelm_active)) then
        allocate(nelm_active(nmatd))
        allocate(nmat_active(16))
        allocate(nrhs_active(4))
       endif
       ic = n_active+1
       ic1 = 0
       nelm_active(1) = ic
       do i = 1, n
        if(node_active(i).ne.0) then
         ic1 = ic1 +1
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
          kb = nelm(j)
          if(node_active(i).ne.0) then
           ic = ic + 1
           nelm_active(ic) = kb
          endif
         enddo
          nelm_active(ic1+1) = ic
        endif
       enddo
c
c identify new block sizes (just for ico2 = 0)
c
      n_activep1 = n_active+1
      nmat_active(1) = 0
      nmat_active(2) = nelm_active(ic1+1)-n_activep1
      nmat_active(3) = nmat_active(2) + nmat_active(2)
      nmat_active(4) = nmat_active(3) + nmat_active(2)
      nrhs_active(1) = 0
      nrhs_active(2) = n_active
c
c modifiy a array- rhs should be fine
c

      ic = n_active+1
      do i = 1, n
        if(node_active(i).ne.0) then
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
          kb = nelm(j)
          if(node_active(kb).ne.0) then
           ic = ic + 1
           a_active(ic-n_activep1+ nmat_active(1)) = a(j-neqp1+nmat(1))
           a_active(ic-n_activep1+ nmat_active(2)) = a(j-neqp1+nmat(2))
           a_active(ic-n_activep1+ nmat_active(3)) = a(j-neqp1+nmat(3))
           a_active(ic-n_activep1+ nmat_active(4)) = a(j-neqp1+nmat(4))
          endif
         enddo
        endif
       enddo
      else
c
c remove connections based on residuals
c
c 
c  generate new pre-conditioner (nop) based on residuals and their gradients
c  
       neqp1= neq+1
       nmatd=nelm(neqp1)
       if(.not.allocated(nelm_active)) then
        allocate(nelm_active(nmatd))
        allocate(nmat_active(16))
        allocate(nrhs_active(4))
       endif
c
c evaluate residuals on a connection by connection basis
c
       irb = 0
       iirb = 0
       ic = 0
       ic1 = 0
       id = 0
       do i = 1, n
         if(idof.eq.1) then
          resid_i = abs(bp(i+nrhs(1)))
         else if(idof.eq.2) then
          resid_i = max(abs(bp(i+nrhs(1))),abs(bp(i+nrhs(2))))
         else if(idof.eq.3) then
          resid_i = max(abs(bp(i+nrhs(1))),abs(bp(i+nrhs(2))),
     &              abs(bp(i+nrhs(3))))
         endif
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         ib = 0
         do j = i1,i2
          kb = nelm(j)
          if(idof.eq.1) then
           resid_kb = abs(bp(kb+nrhs(1)))
          else if(idof.eq.2) then
           resid_kb = max(abs(bp(kb+nrhs(1))),abs(bp(kb+nrhs(2))))
          else if(idof.eq.3) then
           resid_kb = max(abs(bp(kb+nrhs(1))),abs(bp(kb+nrhs(2))),
     &              abs(bp(kb+nrhs(3))))
          endif
          if(max(resid_i,resid_kb).gt.resid_tol) then
c this is a passive connection (zero out)
           ic = ic + 1
           nelm_active(j) = -kb
           ib = 1
          endif
         enddo
c if ib = 0 then all the nodes are less than resid_tol and node is passive 
c (will be renumbered later)
       if(ib.eq.1) then
        id = id +1
        node_active(id) = i
        irb(id) = i
        iirb(i) = id
       else
       endif 
       enddo
c
c at this point we have identified active nodes and active connections
c total active nodes = id, need to label passive nodes (after the active nodes)
c
       n_active = id
       do i = 1,n
        if(iirb(i).eq.0) then
         id = id + 1
         iirb(i) = id
         irb(id) = i
        endif
       enddo
       ic = neqp1
       nop(1) = neqp1
       do id = 1, n_active
         i = node_active(id)
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
          kb = nelm(j)
          if(nelm_active(j).lt.0) then
           ic = ic + 1
           nop(ic) = irb(abs(nelm_active(j)))
          endif
         enddo
         nop(i+1) = ic
       enddo
c
c identify new block sizes (up to 3 dof)
c ic is block size 
      nb(1) = 0
      nb(2) = nop(ic+1)-neqp1
      nb(3) = nop(2) + nop(2)
      nb(4) = nop(3) + nop(2)
      nb(5) = nop(4) + nop(2)
      nb(6) = nop(5) + nop(2)
      nb(7) = nop(6) + nop(2)
      nb(8) = nop(7) + nop(2)
      nb(9) = nop(8) + nop(2)
c
c modifiy a array- rhs should be fine (with irb etc)
c

      endif
      endif   
      return 
      end
