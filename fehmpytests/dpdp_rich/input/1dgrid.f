      program onedgrid
      implicit none
      real(8), allocatable :: sx1(:), a_axy(:)
      real(8), allocatable :: flowin(:), flowout(:), fmflow(:)
      integer, allocatable :: iss(:)
      real(8), allocatable :: ssflow(:)
      real(8), allocatable :: xend(:)
      real(8), allocatable :: x(:)
      real(8), allocatable :: delx(:)
      real(8), allocatable :: sx(:,:)
      integer, allocatable :: ndiv(:)
      integer, allocatable :: nelm(:), istrw(:), nelmdg(:)
      integer neq, ncont, iwtotl, i, ipos, j, nss, k, nskip, idum
      integer nzones, xcount, ndivtot, izone,kb,i1,i2,ic
      integer icoord
      real(8) delxzone, xstart, cross_sectional_area

      open(8,file='1dgrid.stor')
      open(7,file='1dgrid.grid')
      open(5,file='1dgrid_input.dat')

      read(5,*) nzones

      allocate(xend(0:nzones),ndiv(nzones),delx(nzones))

      read(5,*) xstart
      xend(0) = xstart
      do izone = 1, nzones
         read(5,*) ndiv(izone), xend(izone)
         delxzone = xend(izone)-xend(izone-1)
         delx(izone) = delxzone/ndiv(izone)
         ndivtot = ndivtot + ndiv(izone)
      end do
      read(5,*,end=2000) cross_sectional_area, icoord
      goto 2100
 2000 continue
      cross_sectional_area = 1.
      icoord = 1
 2100 continue

      neq = ndivtot+1
      allocate(x(0:neq), sx1(neq))
      x(1) = xstart
      xcount = 1
      write(7,'(a4)') 'coor'
      write(7,1000) neq
      if(icoord.eq.1) then
         write(7,1000) 1, x(1),0., 0.
      elseif(icoord.eq.2) then
         write(7,1000) 1, 0.,x(1),0.
      else
         write(7,1000) 1, 0.,0.,x(1)
      end if

      do izone = 1, nzones
         do j = 1, ndiv(izone)
            xcount = xcount + 1
            x(xcount) = x(xcount-1) + delx(izone)
            if(icoord.eq.1) then
               write(7,1000) xcount, x(xcount),0., 0.
            elseif(icoord.eq.2) then
               write(7,1000) xcount, 0.,x(xcount),0.
            else
               write(7,1000) xcount, 0.,0.,x(xcount)
            end if
         end do
      end do
 1000 format(i7,1x,e13.7,1x,e13.7,1x,e13.7)
      write(7,*)
      write(7,'(a4)') 'elem'
      write(7,1001) 2, 1
      write(7,1001) 1, 2, 1
      write(7,*)
      write(7,'(a4)') 'stop'
      close(7)
 1001 format(i7,1x,i7,1x,i7)
      sx1(1) = abs(0.5*(x(2)-x(1))*cross_sectional_area)
      sx1(neq) = abs(0.5*(x(neq)-x(neq-1))*cross_sectional_area)
      do i = 2, neq-1
         sx1(i) = abs((0.5*(x(i)-x(i-1))+0.5*(x(i+1)-x(i)))
     2        *cross_sectional_area)
      end do



      write(8,*) 'Stor file written by auxiliary code stor1d'
      write(8,*) 'One dimensional stor file w/ correct areas'


      ncont = neq+1+3*(neq-2)+4
      iwtotl = ncont-(neq+1)

      allocate(nelm(ncont),istrw(ncont),nelmdg(neq))
      allocate(sx(iwtotl,3))
      nelm = 0
      istrw = 0
      nelmdg = 0
      sx = 0.
      write(8,'(5i10)') iwtotl, neq, ncont, 1, 3


      ipos = neq+1
      nelm(1) = ipos
      ipos = ipos + 1
      nelm(ipos) = 1
      nelmdg(1) = ipos
      ipos = ipos + 1
      nelm(ipos) = 2



      do i = 2, neq-1
         nelm(i) = ipos
         ipos = ipos + 1
         nelm(ipos) = i-1
         ipos = ipos + 1
         nelm(ipos) = i
         nelmdg(i) = ipos
         ipos = ipos + 1
         nelm(ipos) = i+1
      end do

      nelm(neq) = ipos
      ipos = ipos + 1
      nelm(ipos) = neq-1
      ipos = ipos + 1
      nelm(ipos) = neq
      nelmdg(neq) = ipos
      nelm(neq+1) = ipos


      ic = 0
      do i = 1, neq
         i1 = nelm(i) + 1
         i2 = nelm(i+1)
         do j = i1,i2
            ic = ic + 1
            kb = nelm(j)
            if(i.ne.kb) then
               sx(ic,1) = -cross_sectional_area/abs((x(i)-x(kb)))
            else
               sx(ic,1) = 0.
            end if
            sx(ic,2) = 0.
            sx(ic,3) = 0.
            istrw(j-neq-1) = ic
         end do
      end do

      write(8, '(5(1pe20.10))') (sx1(i), i = 1, neq)
      write(8, '(5i10)'  ) (nelm(i), i = 1, ncont)
      write(8, '(5i10)'  ) (istrw(i), i = 1, ncont)
      write(8, '(5i10)'  ) (nelmdg(i), i = 1, neq)

c      do i=1,3
      i=1
         write(8,'(5(1pe20.10))') (sx(j,i), j = 1, iwtotl)
c      end do

      end
