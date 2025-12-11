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
      use comco2, only : c_axy, c_vxy, skco2, icarb
      use comdi
      use comci, only : rolf
      use comdti
      use comfi, only : qc
      use comflow
      use davidi

      implicit none
      integer addnode, iconn, idummy, i1, i2, ishisfzz, ishiscfzz
      integer flxz_flag, indexa_axy, inneq, inode, izone
      integer ii, iv, ic
      integer, allocatable :: md(:)
      real*8, allocatable :: sumfout(:,:), sumfin(:,:), sumsink(:,:)
      real*8, allocatable :: sumsource(:,:), sumboun(:,:)
c gaz 051619      
      real*8, allocatable :: summass(:,:), avgsat(:,:)
      real*8 ptime, out_tol, s_dum
      character*18 value_string
      character*24 form_string
      character*90, allocatable :: flux_string(:), fluxv_string(:)
      logical matrix_node

      parameter (out_tol = 1.e-15)
      
      if (.not. allocated(sumfout)) then
         allocate (md(nflxz))
         if (wflux_flag .and. vflux_flag .or. icarb .ne. 0) then
             allocate (sumfout(nflxz,2), sumsink(nflxz,2))
             allocate (sumsource(nflxz,2), sumboun(nflxz,2))
		   allocate (sumfin(nflxz,2), summass(nflxz,2),avgsat(nflxz,2))
          else
             allocate (sumfout(nflxz,1), sumsink(nflxz,1))
             allocate (sumsource(nflxz,1), sumboun(nflxz,1))
		   allocate (sumfin(nflxz,1), summass(nflxz,1),avgsat(nflxz,2))
         end if
      end if
             
      if (flxz_flag .eq. 1) then
c     Fluxes are written to tty and output file
         if(iatty.ne.0) then
            write(iatty,*)
            write(iatty,*)
     2           'Total Flux (kg/s) Leaving Zone (flxz macro option)'
            write(iatty,*)
            write(iatty, 1000) 
         end if
         if(ntty.eq.2 .and. iout .ne. 0) then
            write(iout,*)
            write(iout,*)
     2           '      Total Flux Leaving Zone (flxz macro option)'
            write(iout,*)
            write(iout, 1000)
         end if
      else
c     Fluxes are written to flux history file
         if (.not. allocated (flux_string)) 
     &        allocate (flux_string(nflxz), fluxv_string(nflxz))
      end if
c     Compute fluxes out of zone (>0), leaving out flues
c     into other parts of the zone

      if (flxz_flag .eq. 3) then
         ic = 1
         iv = 2
      else
         if (wflux_flag .and. vflux_flag) then
            ii = 1
            iv = 2
         else if (wflux_flag) then
            ii = 1
         else if (vflux_flag) then
            iv = 1
         end if
      end if

	sumfin = 0.
      sumfout = 0.
      sumsink = 0.
      sumsource = 0.
      sumboun = 0.0
      summass = 0.0
      avgsat = 0.0
      md = 0
      do izone = 1, nflxz
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
               md(izone) = md(izone) + 1
c gaz 051519         
                     if(irdof.eq.13) then
                         s_dum = 1.
                     else
                         s_dum = s(inneq)
                     endif
c     Add boundary condition sources
               if (flxz_flag.eq.3) then
                  sumboun(izone,ic) = sumboun(izone,ic) + skco2(inode)
               else
                  if (wflux_flag) 
     &                 sumboun(izone,ii) = sumboun(izone,ii) + sk(inode)
c gaz debug hack 041321
                  i1 = l +izonef(inode) + qc(inode)                 
                  if (vflux_flag) then
                   if(ico2.le.0) then
                    sumboun(izone,iv) = 
     &                 sumboun(izone,iv) + (1.0-s_dum)*sk(inode)
                   else
                      sumboun(izone,iv) =                        
c  gaz 041321correction for vapor boundary term for ngas                 
     &                 sumboun(izone,iv) + qc(inode)  
                   endif
                  endif 
               end if
                     
                     summass(izone,1) = summass(izone,1) +
     2               volume(inneq)*s_dum*rolf(inneq)*ps(inneq) 
                     avgsat(izone,1) = avgsat(izone,1) + 
     2               volume(inneq)*s_dum*ps(inneq) 
                     avgsat(izone,2) = avgsat(izone,2) + 
     2               volume(inneq)*ps(inneq)                   
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
                  if (flxz_flag.eq.3) then
                     if(c_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfout(izone,ic) = sumfout(izone,ic) + 
     &                          c_axy(indexa_axy)
c     sum_vap = sum_vap + c_vxy(indexa_axy)
                           if(nelm(iconn).eq.inneq) then
                              sumsink(izone,ic) = sumsink(izone,ic) + 
     &                             c_axy(indexa_axy)
                           end if
                        end if
                     else if(c_axy(indexa_axy).lt.0.) then
c     add to source sum
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
	                     sumfin(izone,ic) = sumfin(izone,ic) + 
     &                          c_axy(indexa_axy)
	                     if(nelm(iconn).eq.inneq) then
			               sumsource(izone,ic) = sumsource(izone,ic) + 
     &                          c_axy(indexa_axy)
				         end if
					  end if
                     end if
                     if(c_vxy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfout(izone,iv) = sumfout(izone,iv) + 
     &                          c_vxy(indexa_axy)
c     sum_vap = sum_vap + c_vxy(indexa_axy)
                           if(nelm(iconn).eq.inneq) then
                              sumsink(izone,iv) = sumsink(izone,iv) + 
     &                             c_vxy(indexa_axy)
                           end if
                        end if
                     else if(c_vxy(indexa_axy).lt.0.) then
c     add to source sum
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
	                     sumfin(izone,iv) = sumfin(izone,iv) + 
     &                          c_vxy(indexa_axy)
                        if(nelm(iconn).eq.inneq) then
                           sumsource(izone,iv) = sumsource(izone,iv) + 
     &                          c_vxy(indexa_axy)
					  end if
					  endif
                     end if
                  else
                     if(wflux_flag .and. a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone.or.nelm(iconn)
     3                       .eq.inneq) then
                           sumfout(izone,ii) = sumfout(izone,ii) + 
     &                          a_axy(indexa_axy)
                           if(nelm(iconn).eq.inneq) then
                              sumsink(izone,ii) = sumsink(izone,ii) +
     2                             a_axy(indexa_axy)
                           end if
                        end if
                     elseif(wflux_flag.and.a_axy(indexa_axy).lt.0.) then
c     add to source sum
                        if(nelm(iconn).eq.inneq) then
                           sumsource(izone,ii) = sumsource(izone,ii) +
     2                          a_axy(indexa_axy)
                        end if
                     end if

                     if(vflux_flag) then
                        if (a_vxy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                           if(izoneflxz(idummy+nelm(iconn))
     2                          .ne.izone.or.nelm(iconn)
     3                          .eq.inneq) then
                              sumfout(izone,iv) = sumfout(izone,iv) + 
     &                             a_vxy(indexa_axy)
                              if(nelm(iconn).eq.inneq) then
                                 sumsink(izone,iv) = sumsink(izone,iv) +
     2                                a_vxy(indexa_axy)
                              end if
                           end if
                        else if (a_vxy(indexa_axy).lt.0.) then
c     add to source sum
                           if(nelm(iconn).eq.inneq) then
                              sumsource(izone,iv) = sumsource(izone,iv)
     2                             + a_vxy(indexa_axy)
                           end if
                       end if
                     end if
                  end if               
               end do
            end if
         end do
      end do
c     Write results
      if (flxz_flag .eq. 1) then
c     Fluxes are written to tty and output file
         if(wflux_flag) then
            do izone = 1, nflxz
               if(iatty.ne.0) then
                  if (izone .eq. 1) write(iatty, *) 'Water'
                  write(iatty,1046) iflxz(izone), md(izone),
     2                 sumsource(izone,ii), sumsink(izone,ii), 
     3                 sumfout(izone,ii), sumboun(izone,ii), 
     4                 summass(izone,ii),
     5                 avgsat(izone,1)/(avgsat(izone,2) + out_tol)
               end if
               if(ntty.eq.2 .and. iout .ne. 0) then
                  if (izone .eq. 1) write(iout, *) 'Water'
                  write(iout,1046) iflxz(izone), md(izone),
     2                 sumsource(izone,ii), sumsink(izone,ii), 
     3                 sumfout(izone,ii), sumboun(izone,ii),
     4                 summass(izone,ii),
     5                 avgsat(izone,1)/(avgsat(izone,2) + out_tol)           
               end if
            end do
         end if
         if (vflux_flag) then
            do izone = 1, nflxz
               if(iatty.ne.0) then
                  if (izone .eq. 1) write(iatty, *) 'Vapor'
                  write(iatty,1045) iflxz(izone), md(izone),
     2                 sumsource(izone,iv), sumsink(izone,iv), 
     3                 sumfout(izone,iv), sumboun(izone,iv)

               end if
               if(ntty.eq.2 .and. iout .ne. 0) then
                  if (izone .eq. 1) write(iout, *) 'Vapor'
                  write(iout,1045) iflxz(izone), md(izone),
     2                 sumsource(izone,iv), sumsink(izone,iv), 
     3                 sumfout(izone,iv), sumboun(izone,iv)
               end if
            end do
         end if
      else if (flxz_flag .eq. 2) then
c     Fluxes are written to flux history file
         form_string = ''
         if (form_flag .le. 1) then
            form_string = '(1x, g16.9)'
         else
            form_string = '(", ", g16.9)'
         end if
         if (wflux_flag) then
            do izone = 1, nflxz
               flux_string = ''
               if (prnt_flxzvar(1)) then
                  write(value_string, form_string) sumsource(izone,ii)
                  flux_string(izone) = trim(flux_string(izone)) // 
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(2)) then
                  write(value_string, form_string) sumsink(izone,ii)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(3)) then
                  write(value_string, form_string) sumfout(izone,ii)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(4)) then
                  write(value_string, form_string) sumboun(izone,ii)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               ishisfzz = ishisfz + izone
               if (form_flag .le. 1) then
                  write (ishisfzz, 1051) ptime, trim(flux_string(izone))
               else
                  write (ishisfzz, 1057) ptime, trim(flux_string(izone))
               end if
               call flush(ishisfzz)
            end do
         end if
         if (vflux_flag) then
            do izone = 1, nflxz
               flux_string = ''
               if (prnt_flxzvar(1)) then
                  write(value_string, form_string) sumsource(izone,iv)
                  flux_string(izone) = trim(flux_string(izone)) // 
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(2)) then
                  write(value_string, form_string) sumsink(izone,iv)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(3)) then
                  write(value_string, form_string) sumfout(izone,iv)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(4)) then
                  write(value_string, form_string) sumboun(izone,iv)
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               ishisfzz = ishisfz + izone + 400
               if (form_flag .le. 1) then
                  write (ishisfzz, 1051) ptime, trim(flux_string(izone))
               else
                  write (ishisfzz, 1057) ptime, trim(flux_string(izone))
               end if
               call flush(ishisfzz)
            end do
         end if
      elseif (flxz_flag.eq.3) then
c     Fluxes are written to flux history file
         flux_string = ''
         do izone = 1, nflxz
            if (abs(sumsource(izone,ic)) .lt. out_tol) 
     &           sumsource(izone,ic) = 0.d0
            if (abs(sumsink(izone,ic)) .lt. out_tol)  
     &           sumsink(izone,ic) = 0.d0
            if (abs(sumfin(izone,ic)) .lt. out_tol)  
     &           sumfin(izone,ic) = 0.d0
            if (abs(sumfout(izone,ic)) .lt. out_tol)  
     &           sumfout(izone,ic) = 0.d0
            if (abs(sumboun(izone,ic)) .lt. out_tol)  
     &           sumboun(izone,ic) = 0.d0
             if (abs(sumsource(izone,iv)) .lt. out_tol) 
     &           sumsource(izone,iv) = 0.d0
            if (abs(sumsink(izone,iv)) .lt. out_tol)  
     &           sumsink(izone,iv) = 0.d0
            if (abs(sumfin(izone,iv)) .lt. out_tol)  
     &           sumfin(izone,iv) = 0.d0
            if (abs(sumfout(izone,iv)) .lt. out_tol)  
     &           sumfout(izone,iv) = 0.d0
c            if (abs(sumboun(izone,iv)) .lt. out_tol)  
c     &           sumboun(izone,iv) = 0.d0
           if (form_flag .le. 1) then
               write (flux_string(izone), 1050) sumsource(izone,ic), 
     &              sumsink(izone,ic), sumfin(izone,ic),  
     &              sumfout(izone,ic), sumboun(izone,ic)
               write (fluxv_string(izone), 1050) sumsource(izone,iv), 
     &              sumsink(izone,iv), sumfin(izone,iv), 
     &			  sumfout(izone,iv)
            else
               write (flux_string(izone), 1055) sumsource(izone,ic), 
     &              sumsink(izone,ic), sumfin(izone,ic),  
     &              sumfout(izone,ic), sumboun(izone,ic)
               write (fluxv_string(izone), 1055) sumsource(izone,iv), 
     &              sumsink(izone,iv), sumfin(izone,iv),
     &			  sumfout(izone,iv)
            end if
            ishiscfzz = ishiscfz + izone
            if (form_flag .le. 1) then
               write (ishiscfzz, 1052) ptime, trim(flux_string(izone)),
     &              trim(fluxv_string(izone))
            else
               write (ishiscfzz, 1057) ptime, trim(flux_string(izone)),
     &              trim(fluxv_string(izone))
            end if
            call flush(ishiscfzz)
         end do
      end if

 1000 format ('Zone   (#  nodes)       Source       Sink',
     &           '         Net       Boundary     Water Mass',
     &       '    Avg Saturation' )
 1045 format(1x,i7,' (',i8,')',1x,1p,4(1x,e12.5))
 1046 format(1x,i7,' (',i8,')',1x,1p,5(1x,e12.5),4x,e12.5)
 1047 format('(g16.9, ', i3, '(a))')
 1050 format(4(g16.9, 1x), g16.9)
 1051 format(g16.9, 1x, a)
 1052 format(g16.9, 2(1x, a))
 1055 format(4(", ", g16.9))
 1056 format(5(", ", g16.9))
 1057 format(g16.9, a)
 1058 format(g16.9, a, ', ', a)

      end
