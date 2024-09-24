      subroutine structured (iflg)
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
!D1 Control reading of input coordinate data for structured (i,j,k) grid
!D1 Generate connectivity and coefficient file.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 09-FEB-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/structured.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      use comxi
      implicit none

      logical null1, opnd
      integer  mb, mc, ic, inelm, nunit, open_file, froot
      integer ii, i, j, k, nx, ny, nz, nxny
      integer iflg, iwt, iwtlx, iwtly, iwtlz, neq_primary_p1
      integer, allocatable :: nelm_temp(:)
      integer, allocatable :: istrw_temp(:)
c     modflow only input variables
      integer nper, itmuni, lenuni, inmod, nstp
      integer ibcfcb, iwdflg, iwetit, ihdwet, iclayer
      integer, allocatable :: laycbd(:)
      real*8 delxm,delym,delzm,perlen,tsmult,tssize
      real*8 denom,timeconv,alenconv,conv
      real*8 hdry, wetfct, trpy
      real*8 storitiv,transk,vertk
      character(7) :: fdm_status = 'unknown'
      character*120 fdm_name
      character*132 modchar      
c     
      integer mbabs
      real*8 xtemp

      real*8 delx, dely, delz
      real*8 afaczp, afaczpl
      real*8 afacyp, afacypl
      real*8 afacxp, afacxpl
      real*8 afac_tol  
      real*8 x0,y0,z0,xl,yl,zl,z1,z2
      real*8, allocatable :: x(:), y(:), z(:)
      real*8, allocatable :: dx(:), dy(:), dz(:)
      real*8, allocatable :: sx_temp(:,:)
      real*8, allocatable :: ztop(:)
      character*4 macro, keyword, response
      parameter (afac_tol = 1.d-8)
      save x, y, z, dx, dy, dz, ztop, trpy
      save timeconv,alenconv, x0, y0, z0
      save nper, nx, ny, nz, nxny
      save keyword

      if(iflg.eq.0) then
         macro = 'fdm '
c**** node coordinate data ****
         if (iout .ne. 0) write(iout, 6010) macro, 'inpt', incoor
         if (iptty .gt. 0) write(iptty, 6010) 
     &        macro, 'structured grid', incoor
 6010    format(1x, '**** input title : ', a4, ' **** ', a4, ' = ', 
     *        i3, ' ****')

c     ivf = -1 will indicate that the strucured grid (block-centered)
c     will be used
         ivf = -1
         mc = 1
         read (incoor, '(a80)') wdd1
         read(wdd1,*) keyword
         if(keyword.eq.'poin') then 
            read (incoor, *) nx, ny, nz
            allocate(dx(0:nx+1))
            allocate(dy(0:ny+1))
            allocate(dz(0:nz+1))
            dx = 0.
            dy = 0.
            dz = 0.

         else if(keyword.eq.'modf') then
c     read in modflow DIS file finished
            inmod = 99
            read(incoor,*) modchar
            open(inmod,file=modchar,form = 'formatted',status = 'old')
c     read in general information
            read(inmod,*) nz,ny,nx,nper,itmuni,lenuni
            allocate(dx(0:nx+1))
            allocate(dy(0:ny+1))
            allocate(dz(0:nz+1))
            dx = 0.0
            dy = 0.0
            dz = 0.0
c     sort out conversion factors
c     use itmuni and lenuni to change units
            if(itmuni.eq.1) then
               timeconv = 1.0/86400.
            else if(itmuni.eq.2) then
               timeconv = 1.0/1440.
            else if(itmuni.eq.3) then
               timeconv = 1.0/24.  
            else if(itmuni.eq.4) then
               timeconv = 1.0  
            else if(itmuni.eq.5) then
               timeconv = 365.25
            endif
            if(lenuni.eq.1) then
               alenconv = 0.3048 
            else if(lenuni.eq.2) then
               alenconv = 1.0000 
            else if(lenuni.eq.3) then
               alenconv = 1.d-02
            endif
            allocate(laycbd(nz))
            read(inmod,*) (laycbd(i),i=1,nz)
            deallocate(laycbd)
c     read row information
            read(inmod,'(a132)') modchar
            if(modchar(1:8).eq.'CONSTANT') then
               backspace inmod
               read(inmod,*) modchar(1:8),delym
               do i = 1,ny
                  dy(i) = delym*alenconv
               enddo
            else
               backspace inmod
               read(inmod,*) (dy(i),i=1,ny)   
               do i = 1,ny
                  dy(i) = dy(i)*alenconv
               enddo
            endif
c     read column information
            read(inmod,'(a132)') modchar
            if(modchar(1:8).eq.'CONSTANT') then
               backspace inmod
               read(inmod,*) modchar(1:8),delxm
               do i = 1,nx
                  dx(i) = delxm*alenconv
               enddo
            else
               backspace inmod
               read(inmod,*) (dx(i),i=1,nx)   
               do i = 1,nx
                  dx(i) = dx(i)*alenconv
               enddo
            endif
c     read layer information
c     first read top of system
            read(inmod,'(a132)') modchar
            if(modchar(1:8).eq.'CONSTANT') then
               backspace inmod
               read(inmod,*) modchar(1:8), z1     
               z1 = z1*alenconv
               z0 = z1
            else
c     other input forms not supported
               stop
            endif
            do k=1,nz
               read(inmod,'(a132)') modchar
               if(modchar(1:8).eq.'CONSTANT') then
                  backspace inmod
                  read(inmod,*) modchar(1:8), z2        
                  z2 = z2*alenconv
                  dz(k) = z1-z2
                  z1 = z2
               else
c     other input forms not supported
                  stop
               endif
            enddo
c     set coordinates x0 and y0
            x0 = 0.0
            y0 = 0.0
c     read in modflow time information
c     convert to FEHM timestep changes
            tims = 0.0
            icgts = nper
            nstep = 100000
            read(inmod,*) 
     &           perlen,nstp,tsmult,modchar(1:2)
            perlen = perlen*timeconv
            tims = tims + perlen
            iprtout = nstp
            if(tsmult.eq.1.0) then
               tssize = perlen/nstp
            else 
               denom = tsmult**nstp -1.0
               tssize = perlen*(tsmult-1)/(denom)
            endif
            day = tssize
            if(.not.allocated(dit)) allocate(dit(4*nper))
            if(.not.allocated(itc)) allocate(itc(nper))
            do i = 2,nper
               read(inmod,*) 
     &              perlen,nstp,tsmult,modchar(1:2)
               perlen = perlen*timeconv
               tims = tims + perlen
               if(tsmult.eq.1.0) then
                  tssize = perlen/nstp
               else 
                  denom = tsmult**nstp -1.0
                  tssize = perlen*(tsmult-1)/(denom)
               endif
               dit(i) = perlen
               dit(i+nper) = tssize
               dit(i+nper+nper) = 1.0      
               itc(i) = nstp
            enddo
c     close modflow files
            close(inmod)
         else
            read (incoor, *) nx, ny, nz
            read (incoor, *) x0, y0, z0
            allocate(dx(0:nx+1))
            allocate(dy(0:ny+1))
            allocate(dz(0:nz+1))
            dx = 0.
            dy = 0.
            dz = 0.
         endif
         nxny = nx*ny
         neq_primary = nxny*nz

         allocate(x(0:nx+1))
         allocate(y(0:ny+1))
         allocate(z(0:nz+1))
c     read in other data if not a modflow file
         if(keyword.ne.'modf') then
            do i = 1,nx+1
               read (incoor, '(a80)') wdd1
               if (null1(wdd1)) go to 25 
               backspace incoor
               if(keyword.eq.'poin') then
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  x(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(x,mbabs, mc)
               else
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  dx(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(dx,mbabs, mc)
               endif
               mc = mbabs
            enddo        
 25         continue
            do i = 1,ny+1
               read (incoor, '(a80)') wdd1
               if (null1(wdd1)) go to 26 
               backspace incoor
               if(keyword.eq.'poin') then
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  y(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(y,mbabs, mc)
               else
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  dy(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(dy,mbabs, mc)
               endif
               mc = mbabs
            enddo        
 26         continue
            if(keyword.eq.'vart') then
               allocate(thic(0:neq))
               do i = 1,neq  
                  read (incoor, '(a80)') wdd1
                  if (null1(wdd1)) go to 28 
                  backspace incoor
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  thic(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(thic,mbabs, mc)
                  mc = mbabs
               enddo
 28            continue
            else if(keyword.eq.'varz') then
               allocate(thic(0:neq))
               do i = 1,neq  
                  read (incoor, '(a80)') wdd1
                  if (null1(wdd1)) go to 30 
                  backspace incoor
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  thic(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(thic,mbabs, mc)
                  mc = mbabs
               enddo
 30            continue
               allocate(ztop(0:nxny))
               do i = 1,nxny 
                  read (incoor, '(a80)') wdd1
                  if (null1(wdd1)) go to 29 
                  backspace incoor
                  read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  ztop(mbabs)=xtemp
                  if (mb .lt. 0) 
     &                 call interpolate_structured(ztop,mbabs, mc)
                  mc = mbabs
               enddo
 29            continue
            else
               do i = 1,nz+1
                  read (incoor, '(a80)') wdd1
                  if (null1(wdd1)) go to 27 
                  backspace incoor
                  if(keyword.eq.'poin') then
                     read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  z(mbabs)=xtemp
                     if (mb .lt. 0) 
     &                    call interpolate_structured(z,mbabs, mc)
                  else
                     read (incoor, *) mb, xtemp
                  mbabs=iabs(mb)
                  dz(mbabs)=xtemp
                     if (mb .lt. 0) 
     &                    call interpolate_structured(dz,mbabs, mc)
                  endif
                  mc = mbabs
               enddo        
 27            continue
            endif
c     add cells for use later
         endif
         if(keyword.eq.'poin') then
c     by first setting x(0) = 0.0, the coordinates will 
c     be correct for nx,ny,or nz = 1
c     assumes in that case origin = (0,0,0)
c     z assumes z(nz+1)is lowest elevation
            if(nx.eq.1) then
               x(0) = -x(1)
               x(nx+1) = 2.0*x(nx) - x(nx-1)
            else
               x(0) = 2.0*x(1) - x(2)
               x(nx+1) = 2.0*x(nx) - x(nx-1)
            endif
            if(ny.eq.1) then
               y(0) = -y(1)
               y(ny+1) = 2.0*y(ny) - y(ny-1)
            else
               y(0) = 2.0*y(1) - y(2)
               y(ny+1) = 2.0*y(ny) - y(ny-1)
            endif
            if(nz.eq.1) then
               z(nz+1) = -z(1)
               z(0) = 2.0*z(1) - z(2)
            else
               z(0) = 2.0*z(1) - z(2)
               z(nz+1) = 2.0*z(nz) - z(nz-1)
            endif
            x(nx+1) = 2.0*x(nx) - x(nx-1)
            y(ny+1) = 2.0*y(ny) - y(ny-1)
            do i = 1,nx
               dx(i) = 0.5*(x(i+1)-x(i-1))
            enddo
            do i = 1,ny
               dy(i) = 0.5*(y(i+1)-y(i-1))
            enddo
         else
            xl = x0
            do i = 1,nx
               x(i) = xl + 0.5*dx(i)
               xl = xl +dx(i)
            enddo
            x(0) = x0 - 0.5*dx(1)
            x(nx+1) = xl + 0.5*dx(nx)
            yl = y0
            do j = 1,ny
               y(j) = yl + 0.5*dy(j)
               yl = yl +dy(j)
            enddo
            y(0) = y0 - 0.5*dy(1)
            y(ny+1) = yl + 0.5*dy(ny)
            zl = z0
            do k = 1,nz
               z(k) = zl - 0.5*dz(k)
               zl = zl - dz(k)
            enddo
            z(0) = z0 + 0.5*dz(1)
            z(nz+1) = zl - 0.5*dz(nz)
         endif
         if(keyword.eq.'varz') then 
            nxny = nx*ny
            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx
                     ic = i+(j-1)*nx+(k-1)*nxny
                     cord(ic,1) = x(i)
                     cord(ic,2) = y(j)
                     if(k.eq.1) then
                        cord(ic,3) = ztop(ic)-thic(ic)/2.
                     else
                        icm = i+(j-1)*nx+(k-2)*nxny
                        cord(ic,3) = cord(icm,3)-(thic(ic)+thic(icm))/2.
                     endif
                  enddo
               enddo
            enddo
         else if(keyword.eq.'vart') then
            nxny = nx*ny
            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx
                     ic = i+(j-1)*nx+(k-1)*nxny
                     cord(ic,1) = x(i)
                     cord(ic,2) = y(j)
                     if(k.eq.1) then
                        cord(ic,3) = z0-thic(ic)/2.
                     else
                        icm = i+(j-1)*nx+(k-2)*nxny
                        cord(ic,3) = cord(icm,3)-(thic(ic)+thic(icm))/2.
                     endif
                  enddo
               enddo
            enddo
         else  
            ic = 0      
            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx
                     ic = ic+1
                     cord(ic,1) = x(i)
                     cord(ic,2) = y(j)
                     cord(ic,3) = z(k)
                  enddo
               enddo
            enddo
         endif
c set grid origins         
       x_orig = x0
       y_orig = y0
       z_orig = z0   
      else if (iflg.eq.1) then
c     generate connectivity and coefficient file
c     assume block-centered grid
c     assume isotropic stor file
c     assume isotropic stor file
c     set ivf=0 for NO finite volume calcs
c     set ivf = -1 5-27-2002 gav (insure code knows FDM)
         ivf = -1
         isox = 1
         isoy = 1
         isoz = 1
         afacxpl = 0.0
         afacypl = 0.0
         afaczpl = 0.0
         neq_primary_p1 = neq_primary +1
c     estimate approximate number of connections
         allocate(sx_temp(neq_primary*3,1))
         sx_temp = 0.0
         allocate(nelm_temp(neq_primary*7+neq_primary_p1))
         allocate(istrw_temp(neq_primary*7))
c     nelmdg already allocated in allocatemem
c     allocate(nelmdg(neq_primary))
c     start counters
         if(keyword.eq.'varz') keyword = 'vart'
         inelm = neq_primary +1
         nelm_temp(1) = inelm
         iwt = 0
         ic = 0      
         do k = 1, nz
            if(keyword.ne.'vart') delz = -0.5*(z(k+1)-z(k-1))
            do j = 1, ny
               dely = dy(j)                    
               do i = 1, nx
                  delx = dx(i)                  
                  ic = ic+1
                  if(keyword.eq.'vart') delz = thic(ic)
                  sx1(ic) = delx*dely*delz
                  if(k.gt.1) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic-nxny
                     istrw_temp(inelm-neq_primary_p1) = 
     &                    istrw_temp(nelmdg(ic-nxny)+3-neq_primary_p1)
                  endif
                  if(j.gt.1) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic-nx
                     istrw_temp(inelm-neq_primary_p1) = 
     &                    istrw_temp(nelmdg(ic-nx)+2-neq_primary_p1)
                  endif
                  if(i.gt.1) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic-1
                     istrw_temp(inelm-neq_primary_p1) = 
     &                    istrw_temp(nelmdg(ic-1)+1-neq_primary_p1)
                  endif
                  inelm = inelm + 1
                  nelm_temp(inelm) = ic
                  nelmdg(ic) = inelm
                  istrw_temp(inelm-neq_primary_p1) = 0  
                  if(i.lt.nx) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic+1
                     afacxp = -dely*delz/(x(i+1)-x(i))
                     if(abs(afacxpl-afacxp).le.afac_tol) then
                        istrw_temp(inelm-neq_primary_p1) = iwtlx
                     else
                        iwt = iwt +1
                        iwtlx = iwt
                        afacxpl = afacxp
                        sx_temp(iwt,1) = afacxp 
                        istrw_temp(inelm-neq_primary_p1) = iwt
                     endif
                  endif
                  if(j.lt.ny) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic+nx
                     afacyp = -delx*delz/(y(j+1)-y(j))
                     if(abs(afacypl-afacyp).le.afac_tol) then
                        istrw_temp(inelm-neq_primary_p1) = iwtly
                     else
                        iwt = iwt +1
                        iwtly = iwt
                        afacypl = afacyp
                        sx_temp(iwt,1) = afacyp 
                        istrw_temp(inelm-neq_primary_p1) = iwt
                     endif
                  endif
                  if(k.lt.nz) then
                     inelm = inelm + 1
                     nelm_temp(inelm) = ic+nxny
c     afaczp = delx*dely/(z(k+1)-z(k))
c     afaczp = -delx*dely/thic(ic)         
                     afaczp = -delx*dely/delz
                     if(abs(afaczpl-afaczp).le.afac_tol) then
                        istrw_temp(inelm-neq_primary_p1) = iwtlz
                     else
                        iwt = iwt +1
                        iwtlz = iwt
                        afaczpl = afaczp
                        sx_temp(iwt,1) = afaczp 
                        istrw_temp(inelm-neq_primary_p1) = iwt
                     endif
                  endif
                  nelm_temp(ic+1) = inelm
               enddo
            enddo
         enddo

         nr = iwt +1
         allocate(sx(nr,1))
         sx = 0.0
c     mutiply by 1/3(3-D) or 1/2(2-D) to account for isotropy
         if(icnl.eq.0) then
            sx(1:nr,1) = 0.3333333333333d00*sx_temp(1:nr,1)
         else
            sx(1:nr,1) = 0.5d00*sx_temp(1:nr,1)
         endif
         deallocate(sx_temp)
         if(allocated(nelm))then
            deallocate(nelm)
         endif
         allocate(nelm(inelm))
         nelm = nelm_temp(1:inelm)
         deallocate(nelm_temp)
         if(allocated(istrw))then
            deallocate(istrw)
         endif
         allocate(istrw(inelm-neq_primary_p1))
         istrw = istrw_temp(1:inelm-neq_primary_p1)
         deallocate(istrw_temp)

         
c     gaz 11-9-2001 allocation now done in startup
c     if(idpdp.eq.0) then
c     allocate(istrw_itfc(inelm-neq_primary_p1),
c     &             istrw_cold(inelm-neq_primary_p1))
c     else
c     allocate(istrw_itfc(2*inelm-neq_primary_p1),
c     &             istrw_cold(2*inelm-neq_primary_p1))
c     end if
c     istrw_itfc = 0
c     istrw_cold = 0

c     deallocate(x,y,z)
c     if(allocated(dx)) deallocate(dx,dy,dz)
         if(allocated(thic)) deallocate(thic)
         if(allocated(ztop)) deallocate(ztop)
      else if(iflg.eq.2) then
c     read in modflow info from 'BC6' package
         inmod = 99
         read(inpt,*) modchar
         open(inmod,file=modchar,form = 'formatted',status = 'old')
         conv = alenconv*timeconv*9.89d-13/9.66d-6
         read(inmod,*) ibcfcb,hdry,iwdflg,wetfct,iwetit,ihdwet
c     read in some layer information
         allocate(laycbd(nz))
         read(inmod,*) (laycbd(i),i=1,nz)
c     read in anisotropy factor (1 for each layer)
         read(inmod,'(a132)') modchar
         if(modchar(1:8).eq.'CONSTANT') then
            backspace inmod
            read(inmod,*) modchar(1:8), trpy   
         else
c     other input forms not supported
            stop
         endif
         iclayer = 0
         do k=1,nz
c     read in storitivity or specific yield
c     covert to porosity and compressibilities
            if(nper.gt.0) then
               iporos = -1
               ic = iclayer
               read(inmod,'(a132)') modchar
               if(modchar(1:8).eq.'CONSTANT') then
                  backspace inmod
                  read(inmod,*) modchar(1:8), storitiv  
                  do j= 1,ny
                     do i= 1,nx
                        ic = ic +1
                        amgang(ic) = storitiv/alenconv
                     enddo
                  enddo
               else
c     other input forms not supported
                  stop
               endif
            endif
c     read in layer transmissibility or hydraulic connectivity
c     covert to permeability 
            ic = iclayer
            read(inmod,'(a132)') modchar
            if(modchar(1:8).eq.'CONSTANT') then
               backspace inmod
               read(inmod,*) modchar(1:8), transk    
               if(laycbd(k).eq.0.or.laycbd(k).eq.2) then
                  do j= 1,ny
                     do i= 1,nx
                        ic = ic +1
                        pnx(ic) = transk/dz(k)*conv
                        pny(ic) = pnx(ic)*trpy
                     enddo
                  enddo
               else
                  do j= 1,ny
                     do i= 1,nx
                        ic = ic +1
                        pnx(ic) = transk*conv
                        pny(ic) = pnx(ic)*trpy
                     enddo
                  enddo
               endif
            else
c     other input forms not supported
               stop
            endif
c     read in vert. hydraulic connectivity/(layer thickness)
c     covert to porosity and compressibilities
            ic = iclayer
            if(k.lt.nz) then
               read(inmod,'(a132)') modchar
               if(modchar(1:8).eq.'CONSTANT') then
                  backspace inmod
                  read(inmod,*) modchar(1:8), vertk     
                  vertk = vertk*dz(k)
                  do j= 1,ny
                     do i= 1,nx
                        ic = ic +1
                        pnz(ic) = vertk*conv
                     enddo
                  enddo
               else
c     other input forms not supported
                  stop
               endif
               iclayer = ic
            else
               do j= 1,ny
                  do i= 1,nx
                     ic = ic +1
                     pnz(ic) = vertk*conv
                  enddo
               enddo
            endif
         enddo

         deallocate(laycbd)
         close(inmod)
      else if(iflg.eq.3) then
c     write coordinates to a file
         nunit = 0
         fdm_name = ''
! Use grid file name root as default
         if (nmfil(3) .ne. ' ') then
            call file_prefix(nmfil(3), froot)
            if (froot .gt. 100) froot = 100
            fdm_name(1:froot) = nmfil(3)(1:froot)
            fdm_name(froot+1:froot+12) = ".coordinates"
            nunit = open_file(fdm_name,'unknown')
            inquire (file = fdm_name, write = response)
!In case can't write to directory with grid file
            if (response .eq. 'NO') then
               close (nunit)
               if (iout .ne. 0) then 
                  write (iout, *) "Can't write to file:", fdm_name
               end if
               if (iptty .ne. 0) then 
                  write (iptty, *) "Can't write to file:", fdm_name
               end if
               fdm_name = '' 
            end if
         end if
         if (nunit .ne. 0) then
            inquire (nunit, OPENED=opnd)
         else
            opnd = .false.
         end if
         if (.not. opnd) then
            if (null1(root_name)) then
               if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
                  call file_prefix(nmfil(5), froot)
                  if (froot .gt. 100) froot = 100
                  fdm_name(1:froot) = nmfil(5)(1:froot)
               else
                  fdm_name = 'fdm'
                  froot = 3
               end if
            else 
               froot = min(100, len_trim (root_name))
               fdm_name(1:froot) = root_name(1:froot)
            end if
            fdm_name(froot+1:froot+12) = ".coordinates"
            nunit = open_file(fdm_name,'unknown')
         end if
         write (nunit,'(a)') 
     &  'Coordinates & volumes for primary finite difference gridblocks'
         write(nunit,*) neq_primary, nx, ny, nz
         do i=1,neq_primary
            write(nunit,'(i9,1x,1p,3(g15.6,1x),g12.4)') 
     &           i,(cord(i,j),j=1,3),sx1(i)
         enddo
      else if (iflg.eq.4) then
c generate a write element info to file 'fdm_coor_elem.macro'      
         if(nx.gt.1.and.ny.gt.1.and.nz.eq.1) then     
           call generate_elements(nx,ny,nz)  
         else
          if(iout.ne.0) write(iout,411) 

          if(iptty.ne.0) write(iptty,411) 
411       format('>>> generated element connectivity not available for', 
     &          ' non-3D FDM grid, request ignored <<<')
         endif  
      end if

      end
      subroutine interpolate_structured(a,mb,mc)
c     
c     interpolate coodinates(or cell size) for structured grids
c     
      implicit none
      real*8 a(0:*)
      real*8 dela
      integer mb,mc,mdif,mdifm,io

      mdif   =  mb-mc
      mdifm  =  mdif-1
c     
c**** interpolate on coordinates ****
c     
      dela   =  ( a(mb)-a(mc) )/mdif
      do io = 1, mdifm
         a(mc+io) =  a(mc)+io*dela
      end do
      end
      subroutine generate_elements(nx,ny,nz)
c     
c     write out generated element connectivity 
c     for equivalent 3D hex elements
c 
      use comai, only: ns_in, nei_in
      use combi, only: elem_geo, nact_elem_geo
      implicit none
      integer nx,ny,nz,nxny,kk
      integer i,ii,j,k,ne,ns,npoint,il
      integer, allocatable :: nelm_temp(:)
      integer open_file, unit_fdm_elem
      character*14  fdm_elem_file
      character*80 dum_elem
      logical opnd, exists
      fdm_elem_file = 'fdm_elem.macro'
      if(nx.le.1) then
      else if(ny.le.1) then
      else if(nz.le.1) then  
      else
c 3D hex elements
      inquire (file = fdm_elem_file, exist=exists)
      if(.not.exists) then
       il = open_file(fdm_elem_file,'new')
      else
       il = open_file(fdm_elem_file,'old')
       go to 100
      endif
      ne = (nz-1)*(ny-1)*(nx-1)
      ns = 8
      nxny = nx*ny
      allocate(nelm_temp(ne*8))
      kk = 0
      do k = 1, nz-1  
        do j = 1, ny-1 
          do i = 1, nx-1 
           kk = kk + 1
           nelm_temp((kk-1)*ns+1) = kk + nxny
           nelm_temp((kk-1)*ns+2) = kk + nxny + 1
           nelm_temp((kk-1)*ns+3) = kk + nxny + nx + 1
           nelm_temp((kk-1)*ns+4) = kk + nxny + nx
           nelm_temp((kk-1)*ns+5) = kk 
           nelm_temp((kk-1)*ns+6) = kk + 1
           nelm_temp((kk-1)*ns+7) = kk + nx + 1
           nelm_temp((kk-1)*ns+8) = kk + nx
          enddo
        enddo
      enddo   
      write(il,'(a9)') 'elem trad'
      write(il,*) 8, kk
      do i = 1,kk
       write(il,'(9(1x,i8))') i,(nelm_temp((i-1)*ns +j),j=1,8)
      enddo
c gaz 061022 populate elem_geo
       nei_in = kk
       ns_in = 8
       if(.not.allocated(elem_geo)) then
        allocate(elem_geo(nei_in*ns_in))
       else
        deallocate(elem_geo)
        allocate(elem_geo(nei_in*ns_in))
       endif
       nei_in = kk
       ns_in = 8
       elem_geo(1:nei_in*ns_in) = nelm_temp(1:nei_in*ns_in)
      write(il,*)
      write(il,'(a4)') 'stop'
      close(il)
      nact_elem_geo = 1
      endif        
      if (allocated(nelm_temp)) deallocate(nelm_temp)
      return
100   continue
      inquire(file = fdm_elem_file, number = unit_fdm_elem)
      read(unit_fdm_elem,*) dum_elem
      read(unit_fdm_elem,*) ns_in, nei_in
c check ns_in for consistency      
       if(.not.allocated(elem_geo)) then
        allocate(elem_geo(ns_in*nei_in))
       else
        deallocate(elem_geo)
        allocate(elem_geo(ns_in*nei_in))
       endif
       do ii = 1, nei_in
        read (unit_fdm_elem,*) i, (elem_geo((i-1)*ns_in +j),j=1, ns_in)
       enddo
       nact_elem_geo = 1
      return      
      end
