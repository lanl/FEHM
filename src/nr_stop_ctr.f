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
!D1
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
      character*4 vard
      logical nr_chk
      real*8 value,pnrmax,snrmax,snr2max,tnrmax,hnrmax,panrmax
      parameter (max_var = 10)
      save nr_write
      
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
        if(vard(1:1).eq.'c'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),co2f_stop
        endif 
        if(vard(1:1).eq.'h'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),h_stop
        endif
        if(vard(1:4).eq.'step') then
         read(wdd1,*) vard(1:4), strd_iter 
        endif        
        if(vard(1:1).eq.'s'.and.vard(2:2).eq.' ') then
         read(wdd1,*) vard(1:1),s_stop, strd_iter 
        else if(vard(1:3).eq.'sat'.and.vard(4:4).eq.' ') then
         read(wdd1,*) vard(1:3),s_stop, strd_iter
        else if(vard(1:1).eq.'s'.and.vard(2:2).eq.'2') then
         read(wdd1,*) vard(1:2),s2_stop, strd_iter         
        endif
       endif
      enddo
100   continue      
c
c check for sufficient closure terms
c
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
        stop      
       endif
            
      else if(iflg.eq.1) then
c 
c  iteration control
c   
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
        do i = 1,n
         pnrmax = max(bp(i+nrhs(1)),pnrmax)
         if(ifree.eq.0.and.irdof.ne.13.and.ieos(i).eq.2) then
          if(s_stop.ne.0) then
           snrmax = max(bp(i+nrhs(1)),snrmax)
          endif
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(p_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3    
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
        do i = 1,n
         pnrmax = max(bp(i+nrhs(1)),pnrmax)
         if(s_stop.ne.0.and.ieos(i).eq.2) then
          snrmax = max(bp(i+nrhs(2)),snrmax)
         endif
         if(t_stop.ne.0.and.ieos(i).ne.2) then
          tnrmax = max(bp(i+nrhs(2)),tnrmax)
         endif
        enddo
        if(p_stop.ne.0.0.and.pnrmax.gt.p_stop) nr_stop = 3
        if(s_stop.ne.0.0.and.snrmax.gt.s_stop) nr_stop = 3 
        if(t_stop.ne.0.0.and.tnrmax.gt.t_stop) nr_stop = 3  
       endif             
      endif   
      return 
      end
