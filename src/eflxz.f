      subroutine  eflxz(flxz_flag,ptime)
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
      use comci
      use comco2, only : c_axy, c_vxy, skco2
      use comdi
      use comji
      use comdti
      use comflow
      use davidi

      implicit none
      integer addnode, iconn, idummy, i1, i2, ishisfzz, ishiscfzz
      integer flxz_flag, indexa_axy, inneq, inode, izone
      integer ii, iv, ic, jj, i3, i4,iq,i11,kb
      integer, allocatable :: md(:)
      real*8, allocatable :: sumfout(:,:), sumsink(:,:)
      real*8, allocatable :: sumsource(:,:), sumboun(:,:)
      real*8, allocatable :: qheat_tot(:)
      real*8 ptime, out_tol,distol,dis2,sx3c
      real*8 sx2c,thxkb,thykb,thzkb,delx2,dely2,delz2
      real*8 sx2t,sx3t,sxzt,thxi,thyi,thzi
      real*8 qheat
      character*18 value_string
      character*24 form_string
      character*90, allocatable :: flux_string(:)
      logical matrix_node

      parameter (out_tol = 1.e-30,distol = 1.e-15)
      

      if(.not.eflux_flag) return
      if (.not. allocated(sumfout)) then
         allocate (md(nflxz))
         if (eflux_flag) then
             allocate (sumfout(nflxz,2), sumsink(nflxz,2))
             allocate (sumsource(nflxz,2), sumboun(nflxz,2))
             allocate (qheat_tot(nflxz))
         end if
      end if
      ii = 1
      iv = 2            
      if (flxz_flag .eq. 1) then
c     Fluxes are written to tty and output file
         if(iatty.ne.0) then
            write(iatty,*)
            write(iatty,*)
     2           'Total Energy (MW) Leaving Zone (eflxz macro option)'
            write(iatty,*)
            write(iatty,*)'Zone(# nodes)       Source      ',
     &       '  Sink        Net     Boundary     Conduction'
            write(iatty,*)
         end if
         if(ntty.eq.2 .and. iout .ne. 0) then
            write(iout,*)
            write(iout,*)
     2            'Total Energy (MW) Leaving Zone (eflxz macro option)'
            write(iout,*)
            write(iout,*)'Zone(# nodes)       Source      ',
     &       '  Sink         Net     Boundary     Conduction'
            write(iout,*)
         end if
      else
c     Fluxes are written to flux history file
         if (.not. allocated (flux_string)) 
     &        allocate (flux_string(nflxz))
      end if
c     Compute fluxes out of zone (>0), leaving out fluxes
c     into other parts of the zone
      ii = 1

      sumfout = 0.0
      sumsink = 0.0
      sumsource = 0.0
      sumboun = 0.0
      qheat_tot = 0.0
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
               thxi=thx(inode)
               thyi=thy(inode)
               thzi=thz(inode)      
c     Figure out coefficients for heat conduction term
               i1 = nelm(inneq)+1
               i11 = nelmdg(inneq)
               i2 = nelm(inneq+1)
               iq = 0
               do iconn = i1,i11-1
                kb = nelm(iconn)
                i3 = nelmdg(kb)+1
                i4 = nelm(kb+1)
                do jj = i3,i4
                 if(nelm(jj).eq.inneq) then
                  indexa_axy = jj-neq-1
                  iq = iq +1
                  it10(iq)=istrw(indexa_axy)
                  go to 300
                 endif
                enddo
300            continue
               enddo 
               do iconn = i11+1,i2
                indexa_axy = iconn-neq-1
                iq = iq +1
                it10(iq)=istrw(indexa_axy)
               enddo     
                       
               md(izone) = md(izone) + 1
c     Add boundary condition sources
                  sumboun(izone,1) = sumboun(izone,1) + qh(inode)
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node

c     loop over every connecting node
               iq = 0
               do iconn = i1, i2
                  kb = nelm(iconn)
                  indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
c     RJP 07/05/07 changed for CO2 flux
                  if(kb.eq.inode) then
                   if(qh(inode).gt.0.0) sumsink(izone,ii) =
     &             sumsink(izone,ii) + qh(inode)              
                   if(qh(inode).lt.0.0) sumsource(izone,ii) =
     &             sumsource(izone,ii) + qh(inode)              
                  endif                  
c     heat conduction term 
                  if(kb.ne.inode) then
                   iq = iq +1
                   iw = it10(iq)
                   sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                   thxkb=thx(kb)
                   thykb=thy(kb)
                   thzkb=thz(kb)
                   sx2t=2.*thxi*thxkb/(thxi+thxkb)
                   sx3t=2.*thyi*thykb/(thyi+thykb)
                   sxzt=2.*thzi*thzkb/(thzi+thzkb)
                    delx2=max((cord(kb,1)-cord(inode,1))**2,distol)
                    dely2=max((cord(kb,2)-cord(inode,2))**2,distol)
                    delz2=max((cord(kb,3)-cord(inode,3))**2,distol)
                    dis2=delx2+dely2+delz2
                    sx3c=sx2c*dis2/(delx2/sx2t+dely2/sx3t+
     &              delz2/sxzt)
                    qheat = sx3c*(t(kb)-t(inode))
                    qheat_tot(izone) = qheat_tot(izone) + qheat
c                    write(ierr,*) l,izone,inode,kb,qheat,qheat_tot
                  endif                  
                  if (flxz_flag.eq.3) then
                  else
                     if(wflux_flag .and. a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone) then
                           sumfout(izone,ii) = sumfout(izone,ii) + 
     &                          a_axy(indexa_axy)*enlf(inode)
                        end if
                     elseif(wflux_flag.and.a_axy(indexa_axy).lt.0.) then
c     add to source sum
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone) then  
                           sumfout(izone,ii) = sumfout(izone,ii) +
     2                          a_axy(indexa_axy)*enlf(kb)
                        end if
                     end if
                     if(vflux_flag) then
                        if (a_vxy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone) then
                              sumfout(izone,iv) = sumfout(izone,iv) + 
     &                             a_vxy(indexa_axy)*envf(inode)
                           end if
                        else if (a_vxy(indexa_axy).lt.0.) then
c     add to source sum
                        if(izoneflxz(idummy+nelm(iconn))
     2                       .ne.izone) then
                              sumfout(izone,iv) = sumfout(izone,iv)
     2                             + a_vxy(indexa_axy)*envf(kb)
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
         if(eflux_flag) then
            do izone = 1, nflxz
               if(iatty.ne.0) then
                  write(iatty,1045) iflxz(izone), md(izone),
     2                 sumsource(izone,ii)+sumsource(izone,iv), 
     &                 sumsink(izone,ii)+sumsink(izone,iv), 
     3                 sumfout(izone,ii)+sumfout(izone,iv), 
     &                 sumboun(izone,ii)+sumboun(izone,iv),
     &                 qheat_tot(izone)
               end if
               if(ntty.eq.2 .and. iout .ne. 0) then
                  write(iout,1045) iflxz(izone), md(izone),
     2                 sumsource(izone,ii)+sumsource(izone,iv), 
     &                 sumsink(izone,ii)+sumsink(izone,iv), 
     3                 sumfout(izone,ii)+sumfout(izone,iv), 
     &                 sumboun(izone,ii)+sumboun(izone,iv),
     &                 qheat_tot(izone)
               end if
c   gaz 042521  (removed  if(l.eq.1.and.izone.eq.1) then
c                 write(ierr,*) 'days,zone, and conduction heat loss(1)'
c                endif
c               if(izone.eq.1) then
c                write(ierr,*) ' days = ',days
c                do i3 = 1,nflxz 
c                 write(ierr,'(1x,i6,1x,g14.5,1x,i6,g14.5)') 
c     &              iflxz(i3),qheat_tot(i3)  
c                enddo
c               endif
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
         if (eflux_flag) then
            do izone = 1, nflxz
               flux_string = ''
               if (prnt_flxzvar(1)) then
                  write(value_string, form_string) (sumsource(izone,ii)+
     &               sumsource(izone,iv))                 
                  flux_string(izone) = trim(flux_string(izone)) // 
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(2)) then
                  write(value_string, form_string) (sumsink(izone,ii)+
     &               sumsink(izone,iv))                   
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(3)) then
                  write(value_string, form_string) (sumfout(izone,ii)+
     &               sumfout(izone,iv))                    
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(4)) then
                  write(value_string, form_string) (sumboun(izone,ii)+
     &               sumboun(izone,iv))                    
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if
               if (prnt_flxzvar(5)) then
                  write(value_string, form_string) qheat_tot(izone)
               
                  flux_string(izone) = trim(flux_string(izone)) //
     &                 trim(value_string)
               end if              
               ishisfzz = ishisfz + izone + 800
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
            if (abs(sumfout(izone,ic)) .lt. out_tol)  
     &           sumfout(izone,ic) = 0.d0
            if (abs(sumboun(izone,ic)) .lt. out_tol)  
     &           sumboun(izone,ic) = 0.d0
            if (form_flag .le. 1) then
               write (flux_string(izone), 1050) sumsource(izone,ic), 
     &              sumsink(izone,ic), sumfout(izone,ic), 
     &              sumboun(izone,ic)
            else
               write (flux_string(izone), 1055) sumsource(izone,ic), 
     &              sumsink(izone,ic), sumfout(izone,ic), 
     &              sumboun(izone,ic)
            end if
            ishiscfzz = ishiscfz + izone
            if (form_flag .le. 1) then
               write (ishiscfzz, 1051) ptime, trim(flux_string(izone))
            else
               write (ishiscfzz, 1057) ptime, trim(flux_string(izone))
            end if
            call flush(ishiscfzz)
         end do
      end if

 1045 format(1x,i4,' (',i6,')',2x,1p,5(1x,e12.5))
 1046 format(1x,i4,' (',i6,')',2x,1p,5(1x,e12.5))
 1047 format('(g16.9, ', i3, '(a))')
 1050 format(4(g16.9, 1x), g16.9)
 1051 format(g16.9, 1x, a)
 1055 format(4(", ", g16.9))
 1056 format(5(", ", g16.9))
 1057 format(g16.9, a)

      end
