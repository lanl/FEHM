      subroutine setconnarray
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
!D1 Determine and store node connections at which interface reduction
!D1 factor model is applied.
!D1 
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 20-AUG-99, Programmer: B. Robinson  
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/setconnarray.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.4.12 Mass transport at interfaces 
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use combi
      use comdti
      use comai
      use comxi
      implicit none

      integer i, ii1, ii2, idg, neqp1, jmi, jmia, jitfc, kitfc
      integer jm, kb, n_levels, ilevel, add_node, add_istrw
      integer kr, idem_redfac
      integer open_file
      logical exists, match

      if(idpdp.eq.0) then
         n_levels = 1
      else
         n_levels = 2
      end if


c     Loop for either 1 or 2 continua
      do ilevel = 1, n_levels

         if(ilevel .eq. 1) then
            add_node = 0
            add_istrw = 0
         else
            add_node = neq
            add_istrw = nelm(neq+1)-neq-1
         end if

c     Section for determining the interface reduction
c     array
c gaz save interface coordinates for visualization, save ID
c only onpen files for nitfcpairs
      if(nitf_use) then
       num_red_fac = 0
       if(icnl.eq.0) then
        idem_redfac = 3
       else
        idem_redfac = 2
       endif
       do i = 1, nitfcpairs
        ii1 = zone_pair(i,1)
        ii2 = zone_pair(i,2)
        if(ii1.eq.0.and.ii2.eq.0) then
c write error message 
         if(iptty.ne.0) write(iptty,200) 
         if(iout.ne.0) write(iout,200)
         if(ierr.ne.0) write(ierr,200)
         stop
        else if(ii1.eq.0) then
         zone_pair(i,1) = ii2
         zone_pair(i,2) = ii1
        endif
        enddo
200     format('error macro itfc(1): both zones cannot be 0')       
       inquire(file = 'red_fac_itfc.txt', exist = exists)
       if(exists) then       
        i_redfac = open_file('red_fac_itfc.txt', 'old')
       else 
        i_redfac = open_file('red_fac_itfc.txt', 'unknown')
       endif
        write(i_redfac,*) 'itfc information file'
        write(i_redfac,*) 'reduction factor number, node 1, node 2,',
     &   ' factor, x, y, z'
101   format(t20,i7,t30,i8,t40,i8,t50,1p,g14.6,t64,3(2x,g14.6),2(1x,i7))   
      endif
      do i = 1, neq
         
         ii1=nelm(i)+1
         ii2=nelm(i+1)
         idg=nelmdg(i)-ii1+1
         neqp1=neq+1
         jmi=nelmdg(i)
         jmia=jmi-neqp1
         
c         do jm=jmi+1,ii2
         do jm=ii1,ii2
            kb=nelm(jm)+add_node
c     loop over all possible reduction models to search for
c     reduction factor to be applied
         if(nitf_use) then         
          reduction: do jitfc = 1, nitfcpairs
            
            if(izonef_itfc(i+add_node).eq.zone_pair(jitfc,1)) then
c gaz 071219 
                match = .false.
                if(zone_pair(jitfc,2).eq.0.and.
     &           izonef_itfc(kb).ne.izonef_itfc(i+add_node)) then
                 match = .true.
                elseif(izonef_itfc(kb).eq.zone_pair(jitfc,2)) then
                 match = .true.
                endif
                 if(match) then
                  istrw_itfc(jm-neqp1+add_istrw) = jitfc
                  num_red_fac = num_red_fac + 1 
                  write(i_redfac,101) jitfc,i,kb,red_factor(jitfc),
     &            ((cord(i,kr)+cord(kb,kr))/2., kr= 1, idem_redfac),
     &             izonef_itfc(i+add_node), izonef_itfc(kb)             
                 endif
c                 exit reduction
            elseif(izonef_itfc(i+add_node).eq.zone_pair(jitfc,2).or.
     &            zone_pair(jitfc,2).eq.0) then
c gaz 071219
                match = .false.
                if (zone_pair(jitfc,2).eq.0.and.
     &           izonef_itfc(kb).ne.zone_pair(jitfc,1).and.
     &           izonef_itfc(kb).ne.izonef_itfc(i+add_node)) then 
                 match = .true.
                elseif(izonef_itfc(kb).eq.zone_pair(jitfc,2).and.
     &           izonef_itfc(kb).ne.izonef_itfc(i+add_node)) then
                 match = .true.
                endif
                  if(match) then
                   istrw_itfc(jm-neqp1+add_istrw) = jitfc
                  endif
c                  num_red_fac = num_red_fac + 1 
c                  write(i_redfac,101) jitfc,i,kb,red_factor(jitfc),
c    &            ((cord(i,kr)+cord(kb,kr))/2., kr= 1, idem_redfac)      
c                  exit reduction
            end if
          end do reduction
         endif   
         if(ncol_use) then 
          colloid: do kitfc = 1, ncoldpairs
            
            if(izonef_itfc(i+add_node).eq.zonec_pair(kitfc,1)) then
               if(izonef_itfc(kb).eq.zonec_pair(kitfc,2)) then
                  istrw_cold(jm-neqp1+add_istrw) = kitfc
                  exit colloid
               end if
            elseif(izonef_itfc(i+add_node).eq.zonec_pair(kitfc,2)) then
               if(izonef_itfc(kb).eq.zonec_pair(kitfc,1)) then
                  istrw_cold(jm-neqp1+add_istrw) = kitfc
                  exit colloid
               end if
            end if
          end do colloid  
         endif    
      end do
      
      end do

      end do
      if(nitfcpairs.gt.0.and.nitf_use) then
c write output in vtk format          
       call vtk_points(1,i_redfac,num_red_fac,idem_redfac)
       close(i_redfac)
       if(iout.ne.0) write(iout,*) 
     &  '>>> number of modified itfc connections ',num_red_fac,' <<<'
       if(iptty.ne.0) write(iptty,*) 
     &  '>>> number of modified itfc connections ',num_red_fac,' <<<'  
      endif
      end
      subroutine vtk_points(iflg,iunit_num,npoints,idem)
c
c     write out set of points in vtk legacy format
c
      implicit none
      
      integer iflg, iunit_num, npoints, i_redvtk, i
      integer j,kk, idem
      integer open_file
      real*8, allocatable ::   points(:,:) 
      real*8, allocatable ::   ftn_num(:,:) 
      logical exists

      if(iflg.eq.1) then
       inquire(file = 'red_fac_itfc.vtk', exist = exists)
       if(exists) then       
        i_redvtk = open_file('red_fac_itfc.vtk', 'old')
       else 
        i_redvtk = open_file('red_fac_itfc.vtk', 'unknown')
       endif
       allocate(points(npoints,idem))
       allocate (ftn_num(npoints,2))
c read from text file       
       rewind iunit_num
       read(iunit_num,*)
       read(iunit_num,*)
       do i = 1, npoints
        read(iunit_num,*) ftn_num(i,1),kk,kk, 
     &      ftn_num(i,2),(points(i,j),j=1,idem)  
       enddo
c write to vtk file       
       write(i_redvtk,100)
       write(i_redvtk,101) npoints
       do i = 1, npoints
       write(i_redvtk,102) (points(i,j),j=1,idem)
       enddo
       write(i_redvtk,103) npoints, 2*npoints
       do i = 1, npoints
       write(i_redvtk,104) 1, i
       enddo 
       write(i_redvtk,105) npoints
       do i = 1, npoints
       write(i_redvtk,106) 1
       enddo 
       write(i_redvtk,107) npoints
       write(i_redvtk,108) (ftn_num(i,1), i=1,npoints)
       write(i_redvtk,109)
       write(i_redvtk,108) (ftn_num(i,2), i=1,npoints)
       close (i_redvtk)
       deallocate(points,ftn_num) 
      else
      endif
100   format('# vtk DataFile Version 2.0',/,
     &       'FEHM itfc interface points',/,
     &       'ASCII',/,
     &       'DATASET UNSTRUCTURED_GRID') 
101   format('POINTS',3x,i7,1x,'double') 
102   format(1p,g14.6,1x,g14.6,1x,g14.6)  
103   format('CELLS',1x,i7,1x,i7)  
104   format(1x,i7,1x,i7)   
105   format('CELL_TYPES',1x,i7)   
106   format(1x,i7)  
107   format('point_data',1x,i7,/,
     &       'scalars number double 1',/,
     &       'lookup_table default')    
108   format(8(1x,f10.4))  
109   format('scalars itfc_value double 1',/,
     &       'lookup_table default')        
      end
