       subroutine coef_replace_ctr(iflg)
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
!D1 read in area coefficients, volumes, and relace fehm generated term.
!D1 allows for gdpm and gdkm terms
!D1 will only replace identified connnections
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 07-Mar-11, Programmer: George Zyvoloski
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

      implicit none
c    
      integer, allocatable :: lrp(:) 
      integer, allocatable :: nrp(:)  
      integer i,j,ii,ij,jj,kb,ie,iflg,iieosd
      integer node_add,nei_add,ietype
      integer nlayers,nrepeat
      integer max_iter,imodel
      integer ivol,ia,iitemp
      integer neqp1,i1,i2,kk,kb1,kb2
      character*4 vard
      character*80 dummy_string
      logical null1
      real*8 dis_tol, dummy1, dummy2
      parameter (max_iter = 10, dis_tol = 1.e-8)
c      
      if(icoef_replace.eq.0) return
c      
      if(iflg.eq.-1) then
c count entries for storage allocation      
      ivol = 0
      ia = 0
      do i = 1,max_replace
       read(inpt,'(a11)') dummy_string(1:11)
       if(dummy_string(1:11).eq.'           ') then
        go to 999
       else if(dummy_string(1:7).eq.'volume') then
        backspace inpt
        read(inpt,*) dummy_string(1:7), ii, dummy1, dummy2
        ivol = ivol +1        
       else if(dummy_string(1:4).eq.'area') then
        backspace inpt
        read(inpt,*) dummy_string(1:7), ii, jj, dummy1, dummy2
        ia = ia + 1  
       endif
      enddo
      if(iout.ne.0) then
       write(iout,*) 'not enough memory in coef_replace_ctr ',
     & 'max relplacement connection exceeds ', max_replace
      endif
      if(iptty.ne.0) then
       write(iptty,*) 'not enough memory in coef_replace_ctr ',
     & 'max relplacement connection exceeds ', max_replace
      endif  
       write(ierr,*) 'not enough memory in coef_replace_ctr ',
     & 'max relplacement connection exceeds ', max_replace  
       stop        
999   continue      
c            
      allocate(vol_pri(ivol))
      allocate(vol_sec(ivol))
      allocate(ii_vol(ivol))      
      allocate(area_pri(ia))
      allocate(area_sec(ia))
      allocate(icon_area1(ia))
      allocate(icon_area2(ia))       
c 
      ivol_cnt = ivol
      iarea_cnt = ia 
      else if(iflg.eq.0) then     
      ivol = 0
      ia = 0
      do i = 1,max_replace
       read(inpt,'(a11)') dummy_string(1:11)
       if(dummy_string(1:11).eq.'           ') then
        go to 1000
       else if(dummy_string(1:7).eq.'volume') then
        backspace inpt
        ivol = ivol +1
        read(inpt,*) dummy_string(1:7), ii, vol_pri(ivol),vol_sec(ivol)
        ii_vol(ivol) = ii
       else if(dummy_string(1:4).eq.'area') then
        backspace inpt
        ia = ia + 1
        read(inpt,*) dummy_string(1:4), ii, jj, area_pri(ia),      
     &  area_sec(ia)
        area_pri(ia) = abs(area_pri(ia))
        area_sec(ia) = abs(area_sec(ia))
        if(ii.gt.jj) then
         ii = iitemp
         ii = jj
         jj = iitemp
        endif
        icon_area1(ia) = ii
        icon_area2(ia) = jj
       endif       
      enddo
1000  continue    
c       
      else if(iflg.eq.1) then 
c      
c replace connections and volumes
c called after the coefficients have been generated
c
      do i = 1, ivol_cnt
       ii = ii_vol(i)
       sx1(ii) = vol_pri(ivol)
       ij = nelm(nelm(ii+1))
       if(ij.gt.neq_primary) then
c secondary connection found
        sx1(ij) = vol_sec(ivol) 
       endif
      enddo
      neqp1 = neq + 1
      do ia = 1, iarea_cnt
c set primary-primary connection  
       ii = icon_area1(ia)
       kk = icon_area2(ia)
       if(kk.ne.ii) then
       i1 = nelmdg(ii)+1
       i2 = nelm(ii+1)
       do j = i1,i2
        kb = nelm(j)
        if(kb.eq.kk) then      
         iw = istrw(j-neqp1)
          if(isoy.eq.1) then
           sx(iw,1) = -area_pri(ia)/3.   
          else
           sx(iw,1) = -area_pri(ia) 
           sx(iw,2) = 0.
           sx(iw,3) = 0.
          endif
          go to 1001
        endif
       enddo
c set secondary-secondary connection  
1001    continue
        kb1 = nelm(nelm(ii+1))
        kb2 = nelm(nelm(kk+1))
        if(kb1.gt.neq_primary.and.kb2.gt.neq_primary) then
         i1 = nelmdg(kb1)+1
         i2 = nelm(kb1+1)
         do j = i1,i2
          kb = nelm(j+1)
          if(kb.eq.kb2) then
           iw = istrw(j-neqp1)
          if(isoy.eq.1) then
           sx(iw,1) = -area_sec(ia)/3.   
          else
           sx(iw,1) = -area_sec(ia) 
           sx(iw,2) = 0.
           sx(iw,3) = 0.
          endif          
          endif
         enddo     
        endif
       else 
c kk = ii        
c set possible primary-secondary connection  
        j = nelm(ii+1)
        kb = nelm(jj)
         if(kb.gt.neq_primary) then           
          iw = istrw(j-neqp1)
           if(isoy.eq.1) then
            sx(iw,1) = -area_sec(ia)/3.   
           else
            sx(iw,1) = -area_sec(ia) 
            sx(iw,2) = 0.
            sx(iw,3) = 0.
           endif
          endif 
        endif
       enddo      
      
      else if(iflg.eq.2) then
c      
c deallocate memory      
c
      deallocate(vol_pri,vol_sec,area_pri,area_sec)
      deallocate(ii_vol,icon_area1,icon_area2)
       
      endif   
      return 
      end
