      subroutine ingdpm
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
!D1 To read generalized dual porosity and dual permeability model parameters and check the 
!D1 size of arrays.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2
!D2 Initial implementation: 03-NOV-99, Programmer: B. Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/ingdpm.f_a  $
!D2
!D2    Rev 2.3   14 Nov 2001 13:10:06   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.8 Generalized dual-porosity formulation
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comdi
      use comdti
      use comai
      use combi
      use comki
      implicit none

      integer i, idum_gdpm, imodel
      logical null1
      character*80 dummy_string
      integer n_n_n, j, i_chk_nodes, icount
      real*8 gdpm_left, gdpm_right, vol_2nd
      

c     gdpm_flag: if nonzero, determines the geometry of the matrix. 
c     1: parallel plate fractures
c     2: spherical geometry
c     3: 1-d column geometry (edge)
c     4: 1-d column geometry (edge,based on total area)
c     5: 1-d column geometry (block-centered)
c     6: 1-d column geometry (bc,based on total area)
c    11: parallel plate fractures x dir is orthogonal to fracture
c    12: parallel plate fractures y dir is orthogonal to fracture
c    13: parallel plate fractures y dir is orthogonal to fracture
      i_chk_nodes = 0
c gaz 022617 skip this for now
      i_chk_nodes = 1
      go to 998
999   if(i_chk_nodes.eq.0) then
       call backtrack(1,inpt,icount)
       igdpm = 0
       ngdpm_layers= 1
       go to 901
      else
       call backtrack(2,inpt,icount)
       imodel = 0
       idum_gdpm = 1
       ngdpmnodes = ngdpm_actual
c       neq_primary = n0-ngdpmnodes  # neq_primary already calculated
       go to 1000
      endif
998   continue
      if(.not.gdkm_new) then
       read(inpt,*) idum_gdpm, ngdpmnodes      
       neq_primary = n0-ngdpmnodes
      endif
      imodel = 0      
 1000 continue
      read(inpt,'(a80)') dummy_string
      if(.not. null1(dummy_string)) then
         imodel = imodel + 1
         backspace inpt 
c gaz 091116 
c gdkm_dir = 0 older gdkm model (0,1,11) dpdp equivalent
c gdkm_dir = 1 frac orth x dir
c gdkm_dir = 2 frac orth y dir
c gdkm_dir = 3 frac orth z dir         
       if(gdkm_flag.ne.0) then
        if(gdkm_new) then
c input and directional models (including older model             
         read(inpt,*) gdkm_dir(imodel), vol_2nd,
     2         gdpm_x(imodel,1)
         if(gdkm_dir(imodel).eq.4) gdkm_dir(imodel) = 0
          ngdpm_layers(imodel) = 1
          vfrac_primary(imodel) = 1.0 - vol_2nd
        else
         read(inpt,*) gdkm_dir(imodel), vfrac_primary(imodel),
     2         gdpm_x(imodel,1) 
          gdkm_dir(imodel) = 0
          ngdpm_layers(imodel) = 1
c gaz 031918 print out that old format detected 
          if(iout.ne.0) then
          write(iout,*)'>> old gdkm input, non-directional(dpdp-like)',
     &      ' model'
          endif
          if(iptty.ne.0) then
          write(iptty,*)'>> old gdkm input, non-directional(dpdp-like)',
     &      ' model'              
          endif
        endif
       else
         read(inpt,*) ngdpm_layers(imodel), vfrac_primary(imodel),
     2        (gdpm_x(imodel,i),i=1,ngdpm_layers(imodel))           
       endif
         goto 1000
      end if
c     Change from edge coordinates to block centers
c      if(gdpm_flag.eq.4.or.gdpm_flag.eq.3.or.gdpm_flag.eq.1) then
      if(gdpm_flag.eq.4.or.gdpm_flag.eq.3) then      
         do i = 1, imodel
            gdpm_left = gdpm_x(i,1)
            gdpm_x(i,1) = gdpm_x(i,1)/2.0
            do j = 2, ngdpm_layers(i)
               gdpm_right = gdpm_x(i,j)
               gdpm_x(i,j) = (gdpm_right+gdpm_left)/2.0
               gdpm_left = gdpm_right        
            enddo
         enddo      
      endif
c     special fracture model
c     vfrac_primary is the fracture width
c     the other are normalized lengths
      if(gdpm_flag.eq.11) then      
         do i = 1, imodel
            wgt_length(i) = gdpm_x(i,1)
            gdpm_vol(i,1) = gdpm_x(i,1)
            gdpm_left = gdpm_x(i,1)
            gdpm_x(i,1) = gdpm_x(i,1)/2.0
            do j = 2, ngdpm_layers(i)
               gdpm_vol(i,j) = gdpm_x(i,j)
               wgt_length(i) = wgt_length(i) + gdpm_x(i,j)
               gdpm_right = gdpm_x(i,j)
               gdpm_x(i,j) = gdpm_right/2. + gdpm_left
               gdpm_left = gdpm_right + gdpm_left     
            enddo
            wgt_length(i) =  wgt_length(i)
     &       /(dxrg(i)-vfrac_primary(imodel))
            do j = 2, ngdpm_layers(i)
             gdpm_x(i,j) = gdpm_x(i,j)/wgt_length(i)
            enddo          
         enddo      
      endif
c     Set flag to identify which nodes have each gdpm model
901   continue       
      narrays = 1
      itype(1) = 4
      default(1) = 0
      igroup = 2
      macro = 'gdpm'

      call initdata2 (inpt, ischk, neq_primary, narrays, itype, 
     *     default, macroread(10), macro, igroup, ireturn,
     2     i4_1=igdpm(1:neq_primary)) 


c     Determine how many gdpm nodes there actually are
      ngdpm_actual = 0
      do i = 1, neq_primary
         imodel = igdpm(i)
         ngdpm_actual = ngdpm_actual + ngdpm_layers(imodel)
      end do
c gaz 022617 and 042918
c write out actual number of gdkm  or gdpm nodes          
      if(gdkm_flag.ne.0) then
       if(iptty.ne.0) then
       write(iptty,*)
     &  'number of gdkm nodes ', ngdpm_actual  
       endif
       if(iout.ne.0) then
       write(iout,*)
     &  'number of gdkm nodes ', ngdpm_actual             
       endif  
      endif 
      neq = neq_primary+ngdpm_actual
c      n0 = neq
      if(i_chk_nodes.eq.0) then
       i_chk_nodes = 1
       go to 999
      endif
c     Check to see if ngdpmnodes is set large enough, if
c     not, fatal error. If value is set larger than needed,
c     write warning message
      if(gdkm_new) then
c gaz 040517 do nothing here        
      else if(ngdpm_actual .lt. ngdpmnodes) then

         if(iout .ne. 0) then
            write(iout,*) 'In gdpm macro, ngdpmnodes must be reduced'
            write(iout,*) 'to reduce storage requirements'
            write(iout,*) 'A value of ',ngdpm_actual,' is required'
            write(iout,*) 'The current value set is ',ngdpmnodes
         end if


         if(iptty .ne. 0) then
            write(iptty,*) 'In gdpm macro, ngdpmnodes must be reduced'
            write(iptty,*) 'to reduce storage requirements'
            write(iptty,*) 'A value of ',ngdpm_actual,' is required'
            write(iptty,*) 'The current value set is ',ngdpmnodes
         end if

         stop

      elseif(ngdpm_actual .gt. ngdpmnodes) then

         if(iout .ne. 0) then
            write(iout,*) 'Fatal error in gdpm macro'
            write(iout,*) 'A value of ',ngdpm_actual,' is required'
            write(iout,*) 'The current value set is ',ngdpmnodes
            write(iout,*) 'Increase ngdpmnodes to ',ngdpm_actual
            write(iout,*) 'and restart'
         end if


         if(iptty .ne. 0) then
            write(iptty,*) 'Fatal error in gdpm macro'
            write(iptty,*) 'A value of ',ngdpm_actual,' is required'
            write(iptty,*) 'The current value set is ',ngdpmnodes
            write(iptty,*) 'Increase ngdpmnodes to ',ngdpm_actual
            write(iptty,*) 'and restart'
         end if

         stop

      end if
      if(gdkm_flag.ne.0) then
      n_n_n = neq_primary
      do i = 1, neq_primary
       imodel = igdpm(i)
       if(imodel.gt.0) then
c gaz 032218 a secondary node zone always get the primary node value +100           
        n_n_n = n_n_n +1
        izonef(n_n_n) = izonef(i) + 100
       endif
      enddo
      else if(gdpm_flag.ne.0) then

c     Set zones for GDPM nodes for the case in which zone has
c     already been called
c     Convention: all GDPM nodes are assigned a zone that
c     is 100 + the zone number of the primary node

      n_n_n = neq_primary
      do i = 1, neq_primary

c     Loop over all GDPM nodes for primary node i
c     ngdpm_layers(imodel) = 0 for imodel = 0 (i.e. no GDPM nodes)
         imodel = igdpm(i)
         do j = 1, ngdpm_layers(imodel)
            n_n_n = n_n_n + 1
c     Only assign the zone number this way if
c     it hasn't already been assigned a non-zero value
c     for example, in a zone with the nnum option
            if(izonef(n_n_n).eq.0.and.
     &           j.eq.ngdpm_layers(imodel).and.
     &                ngdpm_layers(imodel).gt.1) then
               izonef(n_n_n) = izonef(i) + 200
            else if(izonef(n_n_n).eq.0) then
               izonef(n_n_n) = izonef(i) + 100
            end if
         end do

      end do
      endif
      return
      end
      subroutine backtrack(iflg,iunit,icount)
c
c counts lines and backspaces lines in file
c
      implicit none
      character*80 dummy_string
      logical null1
      integer iflg,iunit,icount,i, max_lines
      parameter (max_lines = 10000)
      if(iflg.eq.0) then
      else if(iflg.eq.1) then
c  read and count lines until a blank is found  
       icount = 0
       do i = 1, max_lines
        read(iunit,'(a80)') dummy_string
        if(.not. null1(dummy_string)) then
          icount = icount +1 
        else
          return
        endif
       enddo
      else if(iflg.eq.2) then
c back space lines
100     continue
        backspace iunit
        read(iunit,'(a80)') dummy_string(1:80)
        if(dummy_string(1:4).ne.'gdkm') then
        backspace iunit
         go to 100
        endif
c        read(iunit,'(a80)') dummy_string(1:80)
        
        return
      else if(iflg.eq.3) then
c call initdata2 here
      endif
      return
      end

