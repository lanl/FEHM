      subroutine  flxz(flxz_flag,ptime)
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read control information for writing history files.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 1-JUN-04,  Programmer: Z. Dash
!D2
!D2 $Log:   /pvcs.config/fehm90/src/flxz.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai
      use combi
      use comdi
      use comdti
      use comflow
      use davidi

      implicit none
      integer addnode, iconn, idummy, i1, i2, ipr_vapor
      integer flxz_flag, indexa_axy, inneq, inode, izone, md
      real*8 ptime, sumfout, sumsink, sumsource, sumboun, sum_vap
      character*90, allocatable :: flux_string(:)
      logical matrix_node

      if(irdof.ne.13 .or. ifree.ne.0) ipr_vapor = 1
      if (flxz_flag .eq. 1) then
! Fluxes are written to tty and output file
         if(ipr_vapor.eq.0)then
            if(iatty.ne.0) then
               write(iatty,*)
               write(iatty,*)
     2              'Total Flux (kg/s) Leaving Zone (flxz macro option)'
               write(iatty,*)
               write(iatty,*)'Zone(# nodes)         Source         ',
     &              'Sink          Net     Boundary'
               write(iatty,*)
            end if
            if(ntty.eq.2 .and. iout .ne. 0) then
               write(iout,*)
               write(iout,*)
     2              '      Total Flux Leaving Zone (flxz macro option)'
               write(iout,*)
               write(iout,*)'Zone(# nodes)         Source         Sink',
     &              '          Net     Boundary'
               write(iout,*)
            end if
         else
            if(iatty.ne.0) then
               write(iatty,*)
               write(iatty,*)
     2              'Total Flux (kg/s) Leaving Zone (flxz macro option)'
               write(iatty,*)
               write(iatty,*)'Zone(# nodes)         Source         ',
     &              'Sink          Net     Boundary       Vapor'
               write(iatty,*)
            end if
            if(ntty.eq.2 .and. iout .ne. 0) then
               write(iout,*)
               write(iout,*)
     2              '      Total Flux Leaving Zone (flxz macro option)'
               write(iout,*)
               write(iout,*)'Zone(# nodes)         Source         Sink',
     &              '          Net     Boundary       Vapor'
               write(iout,*)
            end if
         end if
      else
! Fluxes are written to flux history file
         if (.not. allocated (flux_string)) 
     &        allocate (flux_string(nflxz))
         flux_string = ''
      end if
c     Compute fluxes out of zone (>0), leaving out flues
c     into other parts of the zone

      do izone = 1, nflxz
         sumfout = 0.
         sumsink = 0.
         sumsource = 0.
         sumboun = 0.0
         sum_vap = 0.0
         md=0
c     Loop over all nodes
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
            if(izoneflxz(inode).eq.izone) then
               md = md+1
c     Add boundary condition sources
               sumboun=sumboun + sk(inode)
c     calculate vapor out (assume zero if in at this time)
               if(irdof.ne.13.or.ifree.ne.0) then
                  if(sk(inode).gt.0.0d0) then
                     sum_vap = sum_vap + (1.0-s(inode))*sk(inode)
                  endif
               endif
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
               i1 = nelm(inneq)+1
               i2 = nelm(inneq+1)
c     loop over every connecting node
               do iconn = i1, i2
                  indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
                  if(a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                     if(izoneflxz(idummy+nelm(iconn))
     2                    .ne.izone.or.nelm(iconn)
     3                    .eq.inneq) then
                        sumfout = sumfout + a_axy(indexa_axy)
                        if(nelm(iconn).eq.inneq) then
                           sumsink = sumsink +
     2                          a_axy(indexa_axy)
                        end if
                     end if
                  elseif(a_axy(indexa_axy).lt.0.) then
c     add to source sum
                     if(nelm(iconn).eq.inneq) then
                        sumsource = sumsource +
     2                       a_axy(indexa_axy)
                     end if
                  end if
               end do
            end if
         end do
c     Write results
         if (flxz_flag .eq. 1) then
! Fluxes are written to tty and output file
            if(ipr_vapor.eq.0) then
               if(iatty.ne.0) then
                  write(iatty,1045) iflxz(izone), md,
     2                 sumsource, sumsink, sumfout, sumboun
               end if
               if(ntty.eq.2 .and. iout .ne. 0) then
                  write(iout,1045) iflxz(izone), md,
     2                 sumsource, sumsink, sumfout, sumboun
               end if
            else 
               if(iatty.ne.0) then
                  write(iatty,1046) iflxz(izone), md,
     2                 sumsource, sumsink, sumfout, sumboun, sum_vap
               end if
               if(ntty.eq.2 .and. iout .ne. 0) then
                  write(iout,1046) iflxz(izone), md,
     2                 sumsource, sumsink, sumfout, sumboun, sum_vap
               end if    
            end if
         else
! Fluxes are written to flux history file
            if (form_flag .le. 1) then
               if (ipr_vapor .eq. 0) then
                  write (flux_string(izone), 1050) sumsource, sumsink, 
     &                 sumfout, sumboun
               else
                  write (flux_string(izone), 1050) sumsource, sumsink, 
     &                 sumfout, sumboun, sum_vap
               end if
            else
               if (ipr_vapor .eq. 0) then
                  write (flux_string(izone), 1055) sumsource, sumsink, 
     &                 sumfout, sumboun
               else
                  write (flux_string(izone), 1056) sumsource, sumsink, 
     &                 sumfout, sumboun, sum_vap
               end if
            end if
         end if
 1045    format(1x,i4,' (',i6,')',2x,1p,4(1x,e12.5))
 1046    format(1x,i4,' (',i6,')',2x,1p,5(1x,e12.5))
 1050    format(4(g16.9, 1x), g16.9)
 1051    format(g16.9, 1x, i4, 1x, a)
 1055    format(3(g16.9, ", "), g16.9)
 1056    format(4(g16.9, ", "), g16.9)
 1057    format(g16.9, ", ", i4, ", ", a)

      end do

      if (flxz_flag .eq. 2) then
! Write to flux history file
         do i1 = 1, nflxz
            if (form_flag .le. 1) then
               write (ishisfz, 1051) ptime, iflxz(i1), 
     &              trim(flux_string(i1))
            else
               write (ishisfz, 1057) ptime, iflxz(i1), 
     &              trim(flux_string(i1))
            end if
         end do
      end if

      end
