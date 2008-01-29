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
!D1 To read generalized dual porosity model parameters and check the 
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

      integer i, imodel
      logical null1
      character*80 dummy_string
      integer n_n_n, j
      real*8 gdpm_left, gdpm_right
      

c     gdpm_flag: if nonzero, determines the geometry of the matrix. 
c     1: parallel plate fractures,
c     2: spherical geometry
c     3: 1-d column geometry (edge)
c     4: 1-d column geometry (edge,based on total area)
c     5: 1-d column geometry (block-centered)
c     6: 1-d column geometry (bc,based on total area)
      read(inpt,*) gdpm_flag, ngdpmnodes

      imodel = 0
      neq_primary = n0-ngdpmnodes

 1000 continue
      read(inpt,'(a80)') dummy_string
      if(.not. null1(dummy_string)) then
         imodel = imodel + 1
         backspace inpt
         read(inpt,*) ngdpm_layers(imodel), vfrac_primary(imodel),
     2        (gdpm_x(imodel,i),i=1,ngdpm_layers(imodel))
         goto 1000
      end if
c     Change fron edge coordinates to block centers
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

c     Set flag to identify which nodes have each gdpm model
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

c     Check to see if ngdpmnodes is set large enough, if
c     not, fatal error. If value is set larger than needed,
c     write warning message

      if(ngdpm_actual .lt. ngdpmnodes) then

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
     &           j.eq.ngdpm_layers(imodel)) then
               izonef(n_n_n) = izonef(i) + 200
            else if(izonef(n_n_n).eq.0) then
               izonef(n_n_n) = izonef(i) + 100
            end if
         end do

      end do

      return
      end

