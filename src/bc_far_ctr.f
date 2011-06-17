       subroutine bc_far_ctr(iflg)
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
!D1 manage farfirld boundary conditions, mass balance,NR iterations
!D1
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 29-Nov-10, Programmer: George Zyvoloski
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
c NOTES gaz 112510

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comco2
      use comci
      use combi
      use comdti
      use comki  
      use comai

      implicit none
c    
      integer, allocatable :: lrp(:) 
      integer, allocatable :: nrp(:)  
      integer, allocatable :: md(:)
      integer i,j,ii,ij,jj,kb,ie,iflg,iieosd,i1,i2
      integer node_add,nei_add,ietype,addnode,ic,iv
      integer iconn,indexa_axy,flxz_flag,idof_real
      integer max_iter,imodel,izone_bc
      character*4 vard
      character*80 dummy_string
      logical null1
      real*8 dis_tol,vold,cfrac,wfrac
      integer inode,inneq,izone,idummy
      real*8 sumfout,sumsink,sumsource,sumboun
      logical matrix_node
      parameter (max_iter = 10, dis_tol = 1.e-8)
c      
      if(ibcfar.eq.0) return
c      
      if(iflg.eq.0) then
c count models in scanin, for now allocate blindly
      allocate(ibc_far(100))      
      imodels_far = 0
 1000 continue
      read(inpt,'(a80)') dummy_string
      if(.not. null1(dummy_string)) then
         imodels_far = imodels_far + 1
         backspace inpt
         read(inpt,*) ibc_far(imodels_far)
         goto 1000
      end if      
      allocate(ibc_far_zone(neq_primary))
      ibc_far_zone = 0
c
c input similiar to gdpm
c
c     Set flag to identify which nodes with enrichment model
      narrays = 1
      itype(1) = 4
      default(1) = 0
      igroup = 2
      macro = 'bcfar'

      call initdata2 (inpt, ischk, neq_primary, narrays, itype, 
     *     default, macroread(10), macro, igroup, ireturn,
     *     i4_1=ibc_far_zone(1:neq_primary)) 
c   
      allocate (sumfar(imodels_far,5))    
      else if(iflg.eq.1) then       
c
c  remove nodes in far field zones from mass balance calculations
     

      do i=1,n 
       if(ibc_far_zone(i).ne.0) then
       endif     
      enddo
c      
      else if(iflg.eq.2) then  
c      
c zero out NR correction for farfirld zone members 
c zero out residuals for farfirld zone members 
c    
      if(icarb.ne.0) then  
       idof_real = 3
       do i = 1,n
        if(ibc_far_zone(i).ne.0) then
         do j = 1, idof_real
          bp(i+nrhs(j)) = 0.0d0
         enddo
        endif
       enddo 
      else
      endif   
      else if(iflg.eq.3) then  
c      
c calculate mass and energy flowing into far field zones
c  
c  only for co2 now
      if(icarb.ne.0) then
       do j = 1, 4
        acorr_far(j)= 0.0
       enddo
       do i = 1,n
        if(ibc_far_zone(i).ne.0) then
         vold = sx1(i)
         cfrac = yc(i)
         wfrac = 1.0-cfrac        
         acorr_far(1) = acorr_far(1) + denh(i)*vold*wfrac
         acorr_far(3) = acorr_far(3) + denh(i)*vold*cfrac
     &    + denco2h(i)*vold
         acorr_far(2) = acorr_far(2) + deneh(i)*vold
         acorr_far(4) = acorr_far(4) + deneco2h(i)*vold
        endif
       enddo 
      else  
      endif
      else if(iflg.eq.4) then  
c      
c  calculate fluxes crossing the far field boundaries boundary  
c  valid only for co2 problems now
c      
      ic = 3
      iv = 4
      do izone_bc = 1, imodels_far
c     Loop over all nodes
         izone = ibc_far(izone_bc)
         do inode = 1, n0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
            if(inode.gt.neq) then
               inneq = inode-neq
               matrix_node = .true.
               addnode = nelm(neq+1)-neq-1
               idummy = neq
            else
               matrix_node = .false.
               inneq = inode
               addnode = 0
               idummy = 0
            end if
c     Determine if node is part of the zone being summed
            if(ibc_far_zone(inode).eq.izone) then
c
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
               i1 = nelm(inneq)+1
               i2 = nelm(inneq+1)
c     loop over every connecting node
               do iconn = i1, i2
                  indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
c     RJP 07/05/07 changed for CO2 flux
c     ZVD 10/28/10 add CO2 in vapor
c    gaz          if (flxz_flag.eq.3) then
                     if(c_axy(indexa_axy).ne.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(ibc_far_zone(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfar(izone,ic) = sumfar(izone,ic) + 
     &                          c_axy(indexa_axy)
                        end if
                     end if
                     if(c_vxy(indexa_axy).ne.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(ibc_far_zone(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfar(izone,iv) = sumfar(izone,iv) + 
     &                          c_vxy(indexa_axy)
                        end if
                      end if
c  gaz            else
                     if(a_axy(indexa_axy).ne.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(ibc_far_zone(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfar(izone,ii) = sumfar(izone,ii) + 
     &                          a_axy(indexa_axy)
                        end if
                     end if
c                    if(vflux_flag) then
                        if (a_vxy(indexa_axy).ne.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                           if(ibc_far_zone(idummy+nelm(iconn))
     2                          .ne.izone.or.nelm(iconn)
     3                          .eq.inneq) then
                              sumfar(izone,iv) = sumfar(izone,iv) + 
     &                             a_vxy(indexa_axy)
                           end if
                       end if
c                    end if
c gaz              end if               
               end do
            end if
         end do
      end do     
      
      
      else if(iflg.eq.-1) then
c      
c print out      
c
               if(ii .gt. max_iter) then
                  if(iout.ne.0) 
     *                 write(iout,*) 'failed in enrich_ctr'
                  if(iptty.ne.0)
     *                 write(iptty,*) 'failed in enrich_ctr'
               end if
       
      endif   
      return 
      end
