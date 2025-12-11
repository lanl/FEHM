      subroutine area_vol_weightctr(iflg)
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
!D1  PURPOSE
!D1
!D1      Calulate the area on a boundary for use with 
!D1      boundary condition routines
!D1
!***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
CD3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comrxni
      use comii
      use davidi
      implicit none

      integer iflg,i,j,i1,i2,izone,inode,ii,ij,kb,iafile
      real*8 x1,x2,y1,y2,area_total,perm
      integer maxiarea, iareap, jj
      parameter (maxiarea=100)
C
      real*8 a1, b1, c1, d1, e1, f1
      real*8 crosx,crosy,crosz,dotpr,vecmag,line
      integer open_file
      character*100 wfilename
C
C#######################################################################
c
c     Statement functions for cross product, dot product and vector magnitude
c
      crosx(a1,b1,c1,d1,e1,f1)=b1*f1-c1*e1
      crosy(a1,b1,c1,d1,e1,f1)=c1*d1-a1*f1
      crosz(a1,b1,c1,d1,e1,f1)=a1*e1-b1*d1
      dotpr(a1,b1,c1,d1,e1,f1)=a1*d1 + b1*e1 + c1*f1
      vecmag(a1,b1,c1) = 0.5d0*sqrt(a1**2 + b1**2 + c1**2)
      line(a1,b1) = 0.5d0*sqrt(a1**2 + b1**2)
c
C
C#######################################################################
C
C
      if(iflg.eq.0) then
c
c read input 
c izone_area: zone on which to calculate areas etc for this macro
c iareaf=1 : calculate area of each node in zone
c iareaf=2 : calculate area of based on total area given (area01)
c            and portioned by volume size  
c iareaf=3 : calculate area of based on total area given (area01)
c            and portioned by approximate area size        
c iareaf=4 : calculate weighting based on area (used in boun)
c iareaf=5 : calculate weighting based on area*perm (used in boun)
c iareaf=6 : calculate weighting based on vol*perm (used in boun)
c iareaf=7 : read a list of weights for every node in zone
c iareaf=8 : read a file with weights for every node in zone
c area01= input area
c area02= base impedance parameter         
c
c wgt_area= area or area weight
c
         if(iarea.eq.0) then
            allocate(isarea(maxiarea))
            allocate(ifarea(maxiarea))
            allocate(wgt_area(n0))
            read(inpt,*) narea
            iarea = 1
            nall = n0
            isarea(iarea)=1
            ifarea(iarea)=narea
            if(.not.allocated(izone_area)) then
               allocate(izone_area(max(1,narea)))
               allocate(iareaf(max(1,narea)))
               allocate(area01(max(1,narea)))
               allocate(area02(max(1,narea)))
               allocate(izone_area_nodes(nall))   
               izone_area_nodes = 0
            end if
            do i = 1, narea
               read(inpt,*) 
     &              izone_area(i),iareaf(i),area01(i),area02(i)
               if(iareaf(i).eq.7) then
c     read a list of node numbers and weights
                  read(inpt,*) ii, (jj,wgt_area(jj),ij=1,ii)
               elseif(iareaf(i).eq.8) then
                  read(inpt,'(a100)') wfilename
                  iafile = open_file(wfilename,'old')
                  read(iafile,*)  ii, (jj,wgt_area(jj),ij=1,ii)
                  close(iafile)
               endif
            enddo

c     Loop over each zone for determining izone_area array

            do izone = 1, narea
               do inode = 1, n0
                  if(izonef(inode).eq.izone_area(izone)) then
                     izone_area_nodes(inode) = izone_area(izone)
                  end if
               end do
            end do
         else 
            narea0 = narea
            read(inpt,*) narea
            iarea = iarea+1
            nall = nall +n0
            narea = narea + narea0
            isarea(iarea)=narea0 + 1
            ifarea(iarea)=narea

            allocate(izone_area_temp(max(1,narea)))
            allocate(iareaf_temp(max(1,narea)))
            allocate(area01_temp(max(1,narea)))
            allocate(area02_temp(max(1,narea)))
            allocate(izone_area_nodes_temp(nall))   

            izone_area_temp(1:narea0) = izone_area(1:narea0)
            iareaf_temp(1:narea0) = iareaf(1:narea0)
            area01_temp(1:narea0) = area01(1:narea0)
            area02_temp(1:narea0) = area02(1:narea0)
            izone_area_nodes_temp(1:nall-n0) = 
     &           izone_area_nodes(1:nall-n0)

            deallocate(izone_area,iareaf,area01)
            deallocate(area02,izone_area_nodes)               

            allocate(izone_area(max(1,narea)))
            allocate(iareaf(max(1,narea)))
            allocate(area01(max(1,narea)))
            allocate(area02(max(1,narea)))
            allocate(izone_area_nodes(nall))   
            izone_area_nodes = 0

            izone_area = izone_area_temp
            iareaf = iareaf_temp
            area01 = area01_temp
            area02 = area02_temp
            izone_area_nodes = izone_area_nodes_temp

            deallocate(izone_area_temp,iareaf_temp)
            deallocate(area01_temp)
            deallocate(area02_temp)
            deallocate(izone_area_nodes_temp)               


            do i = narea0+1,narea
               read(inpt,*) 
     &              izone_area(i),iareaf(i),area01(i),area02(i)
            enddo

c     Loop over each zone for determining izone_area array

            do izone = narea0+1, narea
               do inode = 1, n0
                  if(izonef(inode).eq.izone_area(izone)) then
                     izone_area_nodes(inode+n0*(iarea-1)) = 
     &                    izone_area(izone)
                  end if
               end do
            end do
         endif
         continue
      else if(iflg.eq.1) then
C     
C     calculate areas or 
C     
         do ii = 1,iarea

            do izone=isarea(ii),ifarea(ii)
               iareap = iareaf(ii)
               if(iareap.eq.1) then
c     calculate areas based on line segments
                  do inode=1,n0
                     if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                    izone_area(izone)) then
                          i1 = nelm(inode)+1
                          i2 = nelm(inode+1)
                          x1 = cord(inode,1)
                          y1 = cord(inode,2)
                          wgt_area(inode) = 0.
                          do ij = i1,i2
                             kb = nelm(ij)
                             if(izone_area_nodes(kb+n0*(ii-1)).eq.
     &                            izone_area(izone)) then
                                x2 = cord(kb,1)
                                y2 = cord(kb,2)
                                wgt_area(inode) = wgt_area(inode) + 
     &                               line(x1-x2,y1-y2) 
                             endif
                          enddo
                       endif
                    enddo
                 else if(iareap.eq.2) then
c     calculate areas based on total input area and weight by node 
c     volume
c     first calculate volume
                    area02(izone) = 0.0
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          area02(izone) = area02(izone) + sx1(inode) 
                       endif
                    enddo
c     now apply volume weights
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then

                          wgt_area(inode) = area01(izone)*sx1(inode)/
     &                         area02(izone)
                       endif
                    enddo
                 else if(iareap.eq.3) then
c     calculate areas based on total input area and weight by node 
c     volume
c     first calculate volume
                    area02(izone) = 0.0
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          perm = sqrt(pnx(inode)**2+pny(inode)**2
     &                         +pnz(inode)**2)
                          area02(izone) = area02(izone) + sx1(inode)
     &                         *perm 
                       endif
                    enddo
c     now apply volume weights
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          perm = sqrt(pnx(inode)**2+pny(inode)**2
     &                         +pnz(inode)**2)
                          wgt_area(inode) = 
     &                         area01(izone)*perm*sx1(inode)/
     &                         area02(izone)

                       endif
                    enddo
                    
                 else if(iareap.eq.6) then
c     calculate weights based on perm*volume weighting
c     first calculate total perm*volume
                    area02(izone) = 0.0
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          perm = sqrt(pnx(inode)**2+pny(inode)**2
     &                         +pnz(inode)**2)
                          area02(izone) = area02(izone) + sx1(inode)*
     &                         perm 
                       endif
                    enddo
c     now apply perm*volume weights
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          perm = sqrt(pnx(inode)**2+pny(inode)**2
     &                         +pnz(inode)**2)
                          wgt_area(inode) = perm*sx1(inode)/
     &                         area02(izone)			   
                       endif
                    enddo
                 endif
                 
              enddo
           enddo

	else if(iflg.eq.2) then
C     
C     apply areas or weights to boundary conditions
C     
           do ii = 1,iarea
              iareap = iareaf(ii)
              if(iareap.eq.1) then
                 do izone=isarea(ii),ifarea(ii)
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          if(ka(inode).eq.-2) then
                             
                             pflow(inode) = wgt_area(inode)
                             
                          endif
                       endif
                    enddo
                 enddo
              else if(iareap.eq.2.or.iareap.eq.3) then
                 do izone=isarea(ii),ifarea(ii)
                    do inode=1,n0
                       if(izone_area_nodes(inode+n0*(ii-1)).eq.
     &                      izone_area(izone)) then
                          if(ka(inode).eq.-2) then
                             
                             pflow(inode) = wgt_area(inode)
                             
                          endif
                       endif
                    enddo
                 enddo
              endif
           enddo
C     *******************************************************
        endif

        return
        end
